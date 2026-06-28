"""Tests for the psr_solver module (steady PSR extinction-curve solver)."""

import numpy as np
import cantera as ct
import pytest

from pymars import psr_solver
from pymars.psr_solver import (
    trace_extinction_curve,
    _damped_newton,
    _core_residual,
    _core_jacobian,
)


def _ch4_air_gas():
    """Gas object at the stoichiometric CH4/air inlet state (1 atm, 300 K)."""
    gas = ct.Solution("gri30.yaml")
    gas.set_equivalence_ratio(1.0, "CH4", {"O2": 1.0, "N2": 3.76})
    gas.TP = 300.0, ct.one_atm
    return gas


class TestDampedNewton:
    """The finite-difference damped-Newton fallback corrector."""

    def test_scalar_root(self):
        # x**2 - 2 = 0  ->  sqrt(2)
        x, ok = _damped_newton(lambda u: u**2 - 2.0, np.array([1.0]))
        assert ok
        assert x[0] == pytest.approx(np.sqrt(2.0))

    def test_linear_system(self):
        # [u0 + u1 - 3, u0 - u1 - 1] = 0  ->  (2, 1)
        def fun(u):
            return np.array([u[0] + u[1] - 3.0, u[0] - u[1] - 1.0])

        x, ok = _damped_newton(fun, np.array([0.0, 0.0]))
        assert ok
        assert np.allclose(x, [2.0, 1.0])

    def test_nonlinear_system(self):
        # [u0**2 + u1**2 - 1, u0 - u1] = 0  ->  (1/sqrt2, 1/sqrt2) from this guess
        def fun(u):
            return np.array([u[0] ** 2 + u[1] ** 2 - 1.0, u[0] - u[1]])

        x, ok = _damped_newton(fun, np.array([0.8, 0.6]))
        assert ok
        assert np.allclose(x, [1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0)])

    def test_no_root_returns_false(self):
        # constant non-zero residual: the line search can never reduce the norm
        # (and the FD Jacobian is singular), so convergence must fail gracefully
        _x, ok = _damped_newton(lambda u: np.array([1.0]), np.array([0.0]))
        assert not ok


class TestCoreJacobian:
    """The analytic Jacobian must match a finite-difference Jacobian."""

    def test_matches_finite_difference(self):
        gas = _ch4_air_gas()
        Y_in = gas.Y.copy()
        h_in = gas.enthalpy_mass
        Wk = gas.molecular_weights
        h_scale = gas.cp_mass * psr_solver.TREF

        # a representative burning state on the branch
        gas.equilibrate("HP")
        Y = gas.Y.copy()
        temp, tau = 1800.0, 1.0e-3

        jac_analytic = _core_jacobian(gas, Y, temp, tau, gas.P, h_scale, Wk)
        assert jac_analytic.shape == (gas.n_species + 1, gas.n_species + 2)

        # finite-difference Jacobian with respect to (Y_1...Y_K, temp, tau)
        x = np.concatenate([Y, [temp, tau]])

        def res_x(z):
            return _core_residual(
                gas,
                z[: gas.n_species],
                z[gas.n_species],
                z[gas.n_species + 1],
                gas.P,
                Y_in,
                h_in,
                h_scale,
                Wk,
            )

        f0 = res_x(x)
        jac_fd = np.empty((gas.n_species + 1, gas.n_species + 2))
        for idx in range(gas.n_species + 2):
            dx = 1.0e-7 * max(1.0, abs(x[idx]))
            xp = x.copy()
            xp[idx] += dx
            jac_fd[:, idx] = (res_x(xp) - f0) / dx

        # the analytic and FD Jacobians should agree to well within FD accuracy
        assert np.abs(jac_analytic - jac_fd).max() < 1.0e-5 * np.abs(jac_fd).max()

    def test_residual_penalty_for_unphysical_state(self):
        gas = _ch4_air_gas()
        Y_in = gas.Y.copy()
        h_in = gas.enthalpy_mass
        Wk = gas.molecular_weights
        h_scale = gas.cp_mass * psr_solver.TREF

        # non-positive temperature -> large finite penalty (no crash)
        res = _core_residual(gas, Y_in, -1.0, 1.0e-3, gas.P, Y_in, h_in, h_scale, Wk)
        assert res.shape == (gas.n_species + 1,)
        assert np.all(res == psr_solver.PENALTY)


class TestTraceExtinctionCurve:
    """The full pseudo-arclength continuation of the S-curve."""

    def test_structure_and_regression(self):
        gas = _ch4_air_gas()
        result = trace_extinction_curve(gas)

        branch = result["branch"]
        assert branch.ndim == 2 and branch.shape[1] == 2
        assert branch.shape[0] > 20  # a reasonable number of traced points

        tau, _temp = branch[:, 0], branch[:, 1]
        i_ext = int(np.argmin(tau))
        # S-curve shape: the fold is interior, and tau folds back up past it
        assert 0 < i_ext < len(tau) - 1
        assert tau[0] > tau[i_ext]
        assert tau[-1] > tau[i_ext]

        points = result["points"]
        assert set(points) == {"extinction", "near_0.1s", "log_mid"}

        tau_ext, temp_ext, _Y_ext = points["extinction"]
        # extinction is the minimum-tau point of the branch
        assert tau_ext == pytest.approx(tau[i_ext])
        # deterministic-solver regression anchors (generous margins)
        assert tau_ext == pytest.approx(7.9e-5, rel=0.15)
        assert temp_ext == pytest.approx(1706.0, rel=0.05)

        # the 0.1 s sample sits near 0.1 s
        assert points["near_0.1s"][0] == pytest.approx(0.1, rel=0.3)
        # the log-midpoint is ~ the geometric mean of the other two residence times
        tau_mid_expected = np.sqrt(tau_ext * points["near_0.1s"][0])
        assert points["log_mid"][0] == pytest.approx(tau_mid_expected, rel=0.3)

        # each sample carries a valid mass-fraction vector
        for key in points:
            mass_fractions = points[key][2]
            assert mass_fractions.shape == (gas.n_species,)
            assert mass_fractions.sum() == pytest.approx(1.0, abs=1.0e-4)
            assert mass_fractions.min() > -1.0e-4

    def test_deterministic(self):
        ext1 = trace_extinction_curve(_ch4_air_gas())["points"]["extinction"]
        ext2 = trace_extinction_curve(_ch4_air_gas())["points"]["extinction"]
        assert ext1[0] == pytest.approx(ext2[0])
        assert ext1[1] == pytest.approx(ext2[1])

    def test_inert_mixture_raises(self):
        # a non-reacting inlet (pure N2) has no extinction turning point, so the
        # continuation cannot reach a fold and must raise (max_steps kept small)
        gas = ct.Solution("gri30.yaml")
        gas.TPX = 300.0, ct.one_atm, "N2:1.0"
        with pytest.raises(RuntimeError):
            trace_extinction_curve(gas, max_steps=50)
