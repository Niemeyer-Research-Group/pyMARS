r"""Steady perfectly stirred reactor (PSR) extinction-curve solver.

Traces the burning (upper) branch of the steady PSR temperature-vs-residence-time
response curve through the extinction turning point, using pseudo-arclength
continuation so the fold (saddle-node bifurcation) can be passed without the
singular Jacobian or basin-of-attraction problems that defeat naive
time-integration and residence-time-parametrized Newton solves.

The reactor is modeled as adiabatic and constant pressure, so the steady
governing equations are

.. math::

    (Y_k - Y_{k,\mathrm{in}}) - \tau\, \frac{\dot\omega_k W_k}{\rho} = 0
    \qquad (k = 1, \dots, K)

.. math::

    h(T, Y) - h_{\mathrm{in}} = 0

for the species mass fractions :math:`Y_k` and temperature :math:`T`, with the
residence time :math:`\tau` as the continuation parameter. Here
:math:`\dot\omega_k` are the net molar production rates, :math:`W_k` the molecular
weights, :math:`\rho` the mass density, and :math:`h` the mass-specific mixture
enthalpy. The corrector is :func:`scipy.optimize.root` (``hybr``) supplied with an
analytic Jacobian assembled from Cantera's kinetics derivatives, falling back to a
finite-difference damped Newton when needed.

.. moduleauthor:: Kyle Niemeyer
"""

import numpy as np
import cantera as ct
from scipy.optimize import root

#: Temperature scale (K) used to nondimensionalize the unknown vector.
TREF = 1000.0
#: Residual returned for an unphysical state so the corrector backtracks.
PENALTY = 1.0e3
#: Maximum pseudo-arclength step (scaled) on the unstable middle branch. A coarser
#: step there lets the corrector jump off onto a spurious solution; this cap is
#: applied only after the extinction fold, during a full ``stop_at_extinction=False``
#: march toward the ignition turning point.
DS_MAX_MIDDLE = 0.05


def _damped_newton(residual, x0, tol=1.0e-9, maxit=50):
    """Damped Newton with a finite-difference Jacobian and line search.

    Used as a fallback corrector. The backtracking line search retreats from
    unphysical trial states (where ``residual`` returns a large penalty), so the
    iteration stays feasible. Returns ``(x, converged)``.
    """
    u = x0.copy()
    for _ in range(maxit):
        f = residual(u)
        nrm = np.linalg.norm(f)
        if nrm < tol:
            return u, True
        n = u.size
        jac = np.empty((f.size, n))
        for j in range(n):
            du = 1.0e-7 * max(1.0, abs(u[j]))
            up = u.copy()
            up[j] += du
            jac[:, j] = (residual(up) - f) / du
        try:
            step = np.linalg.solve(jac, -f)
        except np.linalg.LinAlgError:
            return u, False
        lam = 1.0
        for _ in range(20):
            if np.linalg.norm(residual(u + lam * step)) < nrm:
                break
            lam *= 0.5
        else:
            return u, False
        u = u + lam * step
    return u, np.linalg.norm(residual(u)) < tol * 100


def _core_residual(gas, Y, T, tau, P, Y_in, h_in, h_scale, Wk):
    """The ``gas.n_species+1`` steady-PSR residuals (species, then energy), well scaled.

    Returns a large finite penalty for unphysical states (non-positive
    temperature, non-finite mass fractions, or a Cantera evaluation error) so the
    corrector backtracks instead of crashing.
    """
    if T <= 0.0 or not np.all(np.isfinite(Y)):
        return np.full(gas.n_species + 1, PENALTY)
    try:
        gas.set_unnormalized_mass_fractions(Y)
        gas.TP = T, P
        wdot = gas.net_production_rates
        rho = gas.density
    except ct.CanteraError:
        return np.full(gas.n_species + 1, PENALTY)
    return np.concatenate(
        [(Y - Y_in) - tau * wdot * Wk / rho, [(gas.enthalpy_mass - h_in) / h_scale]]
    )


def _core_jacobian(gas, Y, temp, tau, pressure, h_scale, Wk):
    """Analytic Jacobian of :func:`_core_residual` with respect to ``(Y, temp, tau)``.

    Returns an array of shape ``(gas.n_species+1, gas.n_species+2)``
    (columns ``Y_1 .. Y_gas.n_species``, ``temp``, ``tau``), assembled from Cantera's
    analytic kinetics derivatives (``net_production_rates_ddCi`` and ``net_production_rates_ddT``)
    combined with the density and enthalpy derivatives via the chain rule.
    """
    gas.set_unnormalized_mass_fractions(Y)
    gas.TP = temp, pressure
    wdot = gas.net_production_rates
    rho = gas.density
    mw = gas.mean_molecular_weight
    ddCi = gas.net_production_rates_ddCi
    ddT = gas.net_production_rates_ddT
    yow = Y / Wk
    drho_dY = -rho * mw / Wk
    dCdY = np.outer(yow, drho_dY)
    dCdY[np.diag_indices(gas.n_species)] += rho / Wk
    ds_dY = Wk[:, None] * ((ddCi @ dCdY) / rho - np.outer(wdot, drho_dY) / rho**2)
    drho_dT = -rho / temp
    ds_dT = Wk * ((ddT + ddCi @ (drho_dT * yow)) / rho - wdot * drho_dT / rho**2)
    jac = np.zeros((gas.n_species + 1, gas.n_species + 2))
    jac[: gas.n_species, : gas.n_species] = np.eye(gas.n_species) - tau * ds_dY
    jac[: gas.n_species, gas.n_species] = -tau * ds_dT
    jac[: gas.n_species, gas.n_species + 1] = -(wdot * Wk / rho)
    jac[gas.n_species, : gas.n_species] = gas.partial_molar_enthalpies / Wk / h_scale
    jac[gas.n_species, gas.n_species] = gas.cp_mass / h_scale
    return jac


def trace_extinction_curve(
    gas,
    stop_at_extinction=True,
    ds0=0.10,
    ds_min=1.0e-3,
    ds_max=0.4,
    tau_start=1.0,
    tau_second=0.5,
    max_steps=6000,
):
    """Trace the PSR response curve through its turning point(s).

    Parameters
    ----------
    gas : cantera.Solution
        Gas object already set to the inlet state (inlet temperature, pressure,
        and composition). It is used as scratch space for kinetics/thermo
        evaluations and its state is modified during the solve.
    stop_at_extinction : bool, optional
        If ``True`` (the default, used by model reduction), stop as soon as the
        extinction turning point is found, returning only the upper (burning)
        branch. If ``False``, continue past the extinction fold along the unstable
        middle branch and, for a sufficiently high inlet temperature,
        around the lower (ignition) turning point onto the weakly-reacting branch,
        returning the complete S-curve. For a cold inlet the ignition turning point
        lies at an impractically large residence time (any reactor below the ignition
        delay sits near the inlet temperature), so the full march instead terminates
        at ``max_steps``.
    ds0, ds_min, ds_max : float, optional
        Initial, minimum, and maximum pseudo-arclength step size (in the scaled
        unknown space).
    tau_start, tau_second : float, optional
        Residence times (s) of the two seed points used to start continuation.
    max_steps : int, optional
        Maximum number of continuation steps.

    Returns
    -------
    dict
        A mapping with two keys:

        - ``branch`` -- a :class:`numpy.ndarray` of shape ``(N, 2)`` whose columns
          are the residence time and temperature ``(tau, T)`` along the traced
          curve. With ``stop_at_extinction=True`` it ends at the extinction
          turning point; with ``stop_at_extinction=False`` it continues onto the
          middle branch (and, for a high enough inlet temperature, around the
          ignition turning point onto the weakly-reacting branch).
        - ``points`` -- a :class:`dict` with keys ``"extinction"``,
          ``"near_0.1s"``, and ``"log_mid"``, each mapping to a tuple
          ``(tau, T, Y)`` of the residence time, temperature, and mass-fraction
          vector at that sample point (all on the upper branch).

    Raises
    ------
    RuntimeError
        If the burning branch cannot be seeded, or the continuation never reaches
        the extinction turning point.

    """
    pressure = gas.P
    Y_in = gas.Y.copy()
    h_in = gas.enthalpy_mass
    Wk = gas.molecular_weights
    h_scale = gas.cp_mass * TREF

    # ---- residual and analytic Jacobian (physical variables Y, T, tau) ----
    def core_res_phys(Y, T, tau):
        return _core_residual(gas, Y, T, tau, pressure, Y_in, h_in, h_scale, Wk)

    def core_jac_phys(Y, T, tau):
        return _core_jacobian(gas, Y, T, tau, pressure, h_scale, Wk)

    # ---- scaled unknown vector u = [Y, theta=T/TREF, sigma=ln tau] ----
    def unpack(u):
        return u[: gas.n_species], u[gas.n_species] * TREF, np.exp(u[gas.n_species + 1])

    def core_res(u):
        Y, T, tau = unpack(u)
        return core_res_phys(Y, T, tau)

    def core_jac(u):
        Y, T, tau = unpack(u)
        jac = core_jac_phys(Y, T, tau).copy()
        jac[:, gas.n_species] *= TREF
        jac[:, gas.n_species + 1] *= tau
        return jac

    def solve_fixed_tau(tau, guess):
        def res(x):
            return core_res_phys(
                x[: gas.n_species], x[gas.n_species], tau
            )  # noqa: E731

        def jac(x):
            return core_jac_phys(x[: gas.n_species], x[gas.n_species], tau)[
                :, : gas.n_species + 1
            ]  # noqa: E731

        sol = root(res, guess, jac=jac, method="hybr", tol=1e-10)
        if sol.success and np.linalg.norm(sol.fun) < 1e-7:
            return sol.x, True
        return _damped_newton(res, guess)

    def correct(res, jac, x0):
        sol = root(res, x0, jac=jac, method="hybr", tol=1e-10)
        if sol.success and np.linalg.norm(sol.fun) < 1e-7:
            return sol.x, True
        return _damped_newton(res, x0)  # fall back to FD damped Newton

    # ---- seed the burning branch from adiabatic equilibrium ----
    gas.equilibrate("HP")
    guess = np.concatenate([gas.Y, [gas.T]])
    x0, ok0 = solve_fixed_tau(tau_start, guess)
    x1, ok1 = solve_fixed_tau(tau_second, x0)
    if not (ok0 and ok1):
        raise RuntimeError("PSR: failed to seed the burning branch")

    u_pp = np.concatenate(
        [x0[: gas.n_species], [x0[gas.n_species] / TREF, np.log(tau_start)]]
    )
    u_prev = np.concatenate(
        [x1[: gas.n_species], [x1[gas.n_species] / TREF, np.log(tau_second)]]
    )

    taus = [tau_start, tau_second]
    temps = [x0[gas.n_species], x1[gas.n_species]]
    states = [x0[: gas.n_species].copy(), x1[: gas.n_species].copy()]
    ds = ds0
    ds_cap = ds_max  # max step, tightened past the extinction fold (see below)

    # Turning points are detected by sign changes of the tangent's sigma = ln(tau)
    # component. The march starts heading down the upper branch (tau decreasing, so
    # sigma_sign < 0); the first flip (tau starts increasing) is the extinction
    # fold, and the second flip (tau starts decreasing again) is the ignition fold,
    # the maximum-tau turning point that closes the S-curve. The second fold is only
    # reachable for a high enough inlet temperature; otherwise the full march runs
    # out at ``max_steps`` along the middle branch.
    sigma = gas.n_species + 1
    sigma_sign = -1.0
    n_folds = 0
    i_ext = None

    for _ in range(max_steps):
        tangent = (u_prev - u_pp) / np.linalg.norm(u_prev - u_pp)

        new_sign = np.sign(tangent[sigma])
        if new_sign != 0.0 and new_sign != sigma_sign:
            sigma_sign = new_sign
            n_folds += 1
            if n_folds == 1:
                # extinction turning point = smallest tau reached so far
                i_ext = int(np.argmin(taus))
                if stop_at_extinction:
                    break
                # the middle branch is steep; tighten the step so the corrector
                # stays on the physical branch en route to the ignition fold
                ds_cap = min(ds_max, DS_MAX_MIDDLE)
                ds = min(ds, ds_cap)
            elif n_folds == 2:
                break  # full march: reached the ignition turning point

        def augmented_res(u, _t=tangent, _up=u_prev, _ds=ds):
            return np.concatenate([core_res(u), [np.dot(_t, u - _up) - _ds]])

        def augmented_jac(u, _t=tangent):
            return np.vstack([core_jac(u), _t])

        u_new, ok = correct(augmented_res, augmented_jac, u_prev + ds * tangent)
        if not ok:
            ds *= 0.5
            if ds < ds_min:
                break
            continue

        Y_new, T_new, tau_new = unpack(u_new)
        taus.append(tau_new)
        temps.append(T_new)
        states.append(Y_new.copy())

        u_pp, u_prev = u_prev, u_new
        ds = min(ds * 1.2, ds_cap)

    if i_ext is None:
        raise RuntimeError(
            "PSR: continuation did not reach the extinction turning point"
        )

    taus = np.asarray(taus)
    temps = np.asarray(temps)
    upper_taus = taus[: i_ext + 1]  # burning branch down to the fold

    # point nearest tau = 0.1 s on the burning branch
    i_p1sec = int(np.argmin(np.abs(upper_taus - 0.1)))
    # log-midpoint between the extinction and 0.1 s points
    tau_mid = np.sqrt(taus[i_ext] * taus[i_p1sec])
    i_mid = int(np.argmin(np.abs(upper_taus - tau_mid)))

    def point(i):
        return (taus[i], temps[i], states[i])

    return {
        "branch": np.column_stack([taus, temps]),
        "points": {
            "extinction": point(i_ext),
            "near_0.1s": point(i_p1sec),
            "log_mid": point(i_mid),
        },
    }
