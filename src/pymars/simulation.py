"""Simulation classes used by pyMARS.

Provides an abstract :class:`BaseSimulation` along with concrete simulation
types: :class:`IgnitionSimulation` (autoignition) and :class:`FlameSimulation`
(one-dimensional laminar flame). Additional types can be added by subclassing
:class:`BaseSimulation`.

.. moduleauthor:: Kyle Niemeyer, Cailin Moore (laminar flame functionality)
"""

import os
import logging
from abc import ABC, abstractmethod

import numpy as np
import h5py
import cantera as ct


class BaseSimulation(ABC):
    """Common interface and shared behavior for a single simulation case.

    Parameters
    ----------
    idx : int
        Identifier index for case
    properties : InputIgnition or InputLaminarFlame
        Object with initial conditions for the simulation
    model : str
        Filename for Cantera-format model to be used
    phase_name : str, optional
        Optional name for phase to load from YAML file (e.g., 'gas').
    path : str, optional
        Path for location of output files

    """

    #: Number of points sampled along each simulation's thermochemical profile.
    num_sample_points = 20

    def __init__(self, idx, properties, model, phase_name="", path=""):
        self.idx = idx
        self.properties = properties
        self.model = model
        self.phase_name = phase_name
        self.path = path

        #: Path to any intermediate data file written during the run. Simulation
        #: types that write one (e.g., autoignition) set this in ``setup_case``;
        #: those that don't (e.g., laminar flame) leave it ``None``.
        self.save_file = None

    def _setup_gas(self):
        """Create the gas object and set its initial temperature, pressure, and composition.

        Returns
        -------
        cantera.Solution
            The initialized gas object (also stored as ``self.gas``).

        """
        self.gas = ct.Solution(self.model, self.phase_name)

        self.gas.TP = (
            self.properties.temperature,
            self.properties.pressure * ct.one_atm,
        )

        # set initial composition using either equivalence ratio or general reactant composition
        if self.properties.equivalence_ratio:
            self.gas.set_equivalence_ratio(
                self.properties.equivalence_ratio,
                self.properties.fuel,
                self.properties.oxidizer,
            )
        else:
            if self.properties.composition_type == "mole":
                self.gas.TPX = (
                    self.properties.temperature,
                    self.properties.pressure * ct.one_atm,
                    self.properties.reactants,
                )
            else:
                self.gas.TPY = (
                    self.properties.temperature,
                    self.properties.pressure * ct.one_atm,
                    self.properties.reactants,
                )

        return self.gas

    @classmethod
    def _sample_profile(cls, temperatures, pressures, mass_fractions):
        """Sample state at evenly-spaced fractions of the total temperature rise.

        Samples ``num_sample_points`` rows of ``[temperature, pressure, *mass_fractions]``
        at the points where the temperature first crosses each evenly-spaced fraction
        of the total rise from the initial to the maximum (final) temperature.

        Parameters
        ----------
        temperatures : numpy.ndarray
            Temperature at each point along the profile (time or space).
        pressures : numpy.ndarray
            Pressure at each profile point.
        mass_fractions : numpy.ndarray
            Mass fractions indexed as ``[profile_point, species]``.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(num_sample_points, 2 + n_species)``.

        """
        delta = 1.0 / cls.num_sample_points
        deltas = np.arange(delta, 1 + delta, delta)

        n_species = mass_fractions.shape[1]
        temperature_initial = temperatures[0]
        temperature_diff = temperatures[-1] - temperature_initial

        sampled_data = np.zeros((len(deltas), 2 + n_species))

        idx = 0
        for point in range(len(temperatures)):
            if idx >= len(deltas):
                break
            if temperatures[point] >= temperature_initial + (
                deltas[idx] * temperature_diff
            ):
                sampled_data[idx, 0:2] = [temperatures[point], pressures[point]]
                sampled_data[idx, 2:] = mass_fractions[point]
                idx += 1

        return sampled_data

    @abstractmethod
    def setup_case(self):
        """Initialize the simulation case."""

    @abstractmethod
    def run_case(self, restart=False):
        """Run the simulation and return its global metric."""

    @abstractmethod
    def calculate(self):
        """Run the case and return only its global metric, without sampling data."""

    @abstractmethod
    def process_results(self, skip_data=False):
        """Process results, returning the metric and (optionally) sampled data."""

    def clean(self):
        """Remove the intermediate data file written during the run, if any.

        Simulation types that write no data file leave ``save_file`` as ``None``,
        so this is a no-op for them and the same call works regardless of type.
        """
        if self.save_file is not None:
            try:
                os.remove(self.save_file)
            except OSError:
                pass


class IgnitionSimulation(BaseSimulation):
    """Class for ignition delay simulations

    Parameters
    ----------
    idx : int
        Identifer index for case
    properties : InputIgnition
        Object with initial conditions for simulation
    model : str
        Filename for Cantera-format model to be used
    phase_name : str, optional
        Optional name for phase to load from YAML file (e.g., 'gas').
    path : str, optional
        Path for location of output files

    """

    def setup_case(self):
        """Initialize simulation case."""
        self._setup_gas()

        # Default maximum number of steps
        self.max_steps = 10000
        if self.properties.max_steps:
            self.max_steps = self.properties.max_steps

        # By default, simulations will run to steady state, with the maximum number of steps
        # given by ``self.max_steps``. Alternatively, an end time (in seconds) can be
        # given in cases where something specific is needed (e.g., longer than normal)
        self.time_end = 0.0
        if self.properties.end_time:
            self.time_end = self.properties.end_time

        if self.properties.kind == "constant pressure":
            self.reac = ct.IdealGasConstPressureMoleReactor(self.gas, clone=False)
        else:
            self.reac = ct.IdealGasMoleReactor(self.gas, clone=False)

        # Create ``ReactorNet`` network. Use an adaptive preconditioner with the
        # mole-based reactors so the integrator runs with a sparse preconditioned
        # GMRES solver, which accelerates integration of large kinetic models.
        self.sim = ct.ReactorNet([self.reac])
        self.sim.preconditioner = ct.AdaptivePreconditioner()

        # Set file for later data file
        self.save_file = os.path.join(self.path, str(self.idx) + ".h5")
        self.sample_points = []

        self.ignition_delay = 0.0

    def _step(self):
        """Advance the reactor network a single time step.

        Thin wrapper around :meth:`cantera.ReactorNet.step`, factored out so the
        integrator failure path can be exercised in tests (the underlying
        Cantera method is read-only and cannot be patched directly).

        Returns
        -------
        float
            The new simulation time, in seconds.

        Raises
        ------
        cantera.CanteraError
            If the integrator fails (e.g. CVODES error from non-finite
            derivatives, as can happen for a degenerate reduced model).

        """
        return self.sim.step()

    def run_case(self, stop_at_ignition=False, restart=False):
        """Run simulation case set up ``setup_case``.

        If no end time is specified for the integration, the function integrates
        to steady state (or a maximum of 10,000 steps, by default). This is done
        by checking whether the system state changes below a certain threshold,
        with the residual computed using feature checking. This is blatantly stolen
        from Cantera's :meth:`cantera.ReactorNet.advance_to_steady_state` method.

        Parameters
        ----------
        stop_at_ignition : bool
            If ``True``, stop integration at ignition point, don't save data.
        restart : bool
            If ``True``, skip if results file exists.

        Returns
        -------
        self.ignition_delay : float
            Computed ignition delay in seconds

        """

        if restart and os.path.isfile(self.save_file):
            print("Skipped existing case ", self.idx)
            return

        # Collect simulation results in lists, then write in bulk.
        times = []
        temperatures = []
        pressures = []
        mass_fractions = []

        def record():
            times.append(self.sim.time)
            temperatures.append(self.reac.T)
            pressures.append(self.reac.phase.P)
            mass_fractions.append(self.reac.Y.copy())

        # Save initial conditions
        record()

        ignition_flag = False

        # Main time integration loop. A CanteraError here means the integrator
        # failed (e.g. CVODES non-finite derivatives). This path runs the
        # original/baseline model during a sampling run, where a failure to
        # simulate should halt the reduction, so it is re-raised as a
        # RuntimeError.
        try:
            if self.time_end:
                # if end time specified, continue integration until reaching that time
                while self.sim.time < self.time_end:
                    self._step()
                    record()

                    if (
                        self.reac.T >= self.properties.temperature + 400.0
                        and not ignition_flag
                    ):
                        self.ignition_delay = self.sim.time
                        ignition_flag = True

                        if stop_at_ignition:
                            break

            else:
                # otherwise, integrate until steady state, or maximum number of steps reached
                self.sim.reinitialize()
                max_state_values = self.sim.get_state()
                residual_threshold = 10.0 * self.sim.rtol
                absolute_tolerance = self.sim.atol

                for step in range(self.max_steps):
                    previous_state = self.sim.get_state()

                    self._step()
                    record()

                    if (
                        self.reac.T >= self.properties.temperature + 400.0
                        and not ignition_flag
                    ):
                        self.ignition_delay = self.sim.time
                        ignition_flag = True

                        if stop_at_ignition:
                            break

                    state = self.sim.get_state()
                    max_state_values = np.maximum(max_state_values, state)
                    residual = np.linalg.norm(
                        (state - previous_state)
                        / (max_state_values + absolute_tolerance)
                    ) / np.sqrt(self.sim.n_vars)

                    if residual < residual_threshold:
                        break

                if step == self.max_steps - 1:
                    logging.error(
                        "Maximum number of steps reached before "
                        f"convergence for ignition case {self.idx}"
                    )
                    raise RuntimeError(
                        "Maximum number of steps reached before "
                        f"convergence for ignition case {self.idx}"
                    )
        except ct.CanteraError as error:
            logging.error(f"Integration failed for ignition case {self.idx}: {error}")
            raise RuntimeError(
                f"Integration failed for ignition case {self.idx}"
            ) from error

        # Write collected data to HDF5 file
        with h5py.File(self.save_file, "w") as h5file:
            grp = h5file.create_group("simulation")
            grp.create_dataset("time", data=np.array(times))
            grp.create_dataset("temperature", data=np.array(temperatures))
            grp.create_dataset("pressure", data=np.array(pressures))
            grp.create_dataset("mass_fractions", data=np.array(mass_fractions))

        if not ignition_flag:
            logging.error(f"No ignition detected for ignition case {self.idx}")
            raise RuntimeError(f"No ignition detected for ignition case {self.idx}")

        return self.ignition_delay

    def calculate(self):
        """Run simulation case, just for ignition delay.

        Returns an ignition delay of ``0.0`` (rather than raising) when the
        model does not ignite or the integrator fails (e.g., a CVODES error),
        so a reduced candidate that can no longer be integrated is rejected
        through the error metric instead of aborting the reduction. This mirrors
        ``FlameSimulation.calculate`` and the no-flame handling (see issue #69).

        Returns
        -------
        float
            Computed ignition delay in seconds, or ``0.0`` if the model does not
            ignite or the integration fails.

        """
        # Main time integration loop. A CanteraError here means the integrator
        # failed; for a candidate reduced model this must not halt the reduction,
        # so it is treated as a non-igniting (zero ignition delay) case.
        try:
            if self.time_end:
                # if end time specified, continue integration until reaching that time
                while self.sim.time < self.time_end:
                    self._step()
                    if self.reac.T >= self.properties.temperature + 400.0:
                        self.ignition_delay = self.sim.time
                        break
                if not self.ignition_delay:
                    logging.warning(
                        f"No ignition detected before end time for ignition case {self.idx}"
                    )
            else:
                # otherwise, integrate until steady state, or maximum number of steps reached
                for step in range(self.max_steps):
                    self._step()
                    if self.reac.T >= self.properties.temperature + 400.0:
                        self.ignition_delay = self.sim.time
                        break
                if step == self.max_steps - 1:
                    logging.warning(
                        "Maximum number of steps reached before "
                        f"convergence for ignition case {self.idx}"
                    )
        except ct.CanteraError as error:
            logging.warning(
                f"Integration failed for ignition case {self.idx}; "
                f"treating the ignition delay as zero ({error})"
            )
            self.ignition_delay = 0.0

        return self.ignition_delay

    def process_results(self, skip_data=False):
        """Process integration results to sample data

        Parameters
        ----------
        skip_data : bool
            Flag to skip sampling thermochemical data

        Returns
        -------
        tuple of float, numpy.ndarray or float
            Ignition delay, or ignition delay and sampled data

        """
        # Load saved integration results
        self.save_file = os.path.join(self.path, str(self.idx) + ".h5")
        with h5py.File(self.save_file, "r") as h5file:
            grp = h5file["simulation"]
            times = grp["time"][:]
            temperatures = grp["temperature"][:]
            pressures = grp["pressure"][:]
            mass_fractions = grp["mass_fractions"][:]

        temperature_initial = temperatures[0]

        # ignition delay: first time the temperature rises 400 K above its initial value
        self.ignition_delay = 0.0
        for time, temp in zip(times, temperatures):
            if temp >= temperature_initial + 400.0:
                self.ignition_delay = time
                break

        if skip_data:
            return self.ignition_delay

        sampled_data = self._sample_profile(temperatures, pressures, mass_fractions)
        self.sampled_data = sampled_data
        return self.ignition_delay, sampled_data


class FlameSimulation(BaseSimulation):
    """Class for one-dimensional freely-propagating laminar flame simulations.

    .. moduleauthor:: Cailin Moore

    Parameters
    ----------
    idx : int
        Identifier index for case
    properties : InputLaminarFlame
        Object with initial conditions for simulation
    model : str
        Filename for Cantera-format model to be used
    phase_name : str, optional
        Optional name for phase to load from YAML file (e.g., 'gas').
    path : str, optional
        Path for location of output files

    """

    #: Default minimum physically-meaningful laminar flame speed, in m/s. A solved
    #: flame speed at or below the active floor (negative or near-zero) is treated
    #: as "no flame detected" -- a degenerate, non-physical result -- and handled
    #: the same as a solver failure (a non-flammable mixture tends to "solve" to
    #: such a speed rather than raising). The floor can be raised or lowered per
    #: run via the ``min_flame_speed`` constructor argument (e.g., for a fuel with
    #: genuinely low flame speeds).
    min_flame_speed = 0.05

    def __init__(
        self, idx, properties, model, phase_name="", path="", min_flame_speed=None
    ):
        super().__init__(idx, properties, model, phase_name=phase_name, path=path)
        # ``None`` keeps the class-level default; otherwise override the floor.
        if min_flame_speed is not None:
            self.min_flame_speed = min_flame_speed

    def setup_case(self):
        """Initialize simulation case."""
        self._setup_gas()

        self.flame_speed = 0.0

        # Create the freely-propagating flame object
        self.flame = ct.FreeFlame(self.gas, width=self.properties.width)
        self.flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)

    def _solve_flame(self):
        """Solve the freely-propagating flame and store the unburned flame speed.

        Returns
        -------
        float
            Computed laminar flame speed in m/s.

        Raises
        ------
        cantera.CanteraError
            If the flame fails to solve (i.e., no flame is detected).

        """
        self.flame.solve(loglevel=0, refine_grid=True, auto=True)
        self.flame_speed = self.flame.velocity[0]
        return self.flame_speed

    def _flame_detected(self):
        """Solve the flame and return its speed, or ``None`` if no physical flame.

        "No flame" means either the solver raised a ``CanteraError`` or the
        resulting flame speed is non-physical -- negative or at/below
        ``min_flame_speed`` (a degenerate solution, as a non-flammable mixture
        tends to converge to a near-zero speed rather than raising).

        Returns
        -------
        float or None
            The laminar flame speed in m/s, or ``None`` if no flame is detected.

        """
        try:
            speed = self._solve_flame()
        except ct.CanteraError:
            return None
        if speed <= self.min_flame_speed:
            return None
        return speed

    def run_case(self, restart=False):
        """Solve the laminar flame and return the unburned flame speed.

        Raises ``RuntimeError`` if no flame is detected (solver failure or a
        degenerate, non-physical flame speed). This is the path used for the
        original (baseline) model, where an undetectable flame should halt the
        reduction rather than be silently ignored.

        Parameters
        ----------
        restart : bool
            Unused; retained for interface symmetry with ``IgnitionSimulation``.

        Returns
        -------
        float
            Computed laminar flame speed in m/s

        """
        speed = self._flame_detected()
        if speed is None:
            logging.error(f"No flame detected for laminar flame case {self.idx}")
            raise RuntimeError(f"No flame detected for laminar flame case {self.idx}")
        self.flame_speed = speed
        return speed

    def calculate(self):
        """Solve the flame and return only the flame speed.

        Returns a flame speed of ``0.0`` (rather than raising) when no flame is
        detected -- a solver failure or a degenerate, non-physical flame speed --
        so that a reduced model that can no longer sustain a flame is rejected
        through the error metric instead of aborting the reduction, mirroring how
        ``IgnitionSimulation.calculate`` treats a non-igniting model.

        Returns
        -------
        float
            Computed laminar flame speed in m/s, or ``0.0`` if no flame is detected.

        """
        speed = self._flame_detected()
        if speed is None:
            logging.warning(
                f"No flame detected for laminar flame case {self.idx}; "
                "treating the flame speed as zero"
            )
            speed = 0.0
        self.flame_speed = speed
        return speed

    def process_results(self, skip_data=False):
        """Solve the flame and sample data along the flame profile.

        Uses ``run_case``, so an undetectable flame raises ``RuntimeError``; this
        path samples data for the original model, where a missing flame should
        halt the reduction.

        Parameters
        ----------
        skip_data : bool
            Flag to skip sampling thermochemical data

        Returns
        -------
        float, or tuple of float and numpy.ndarray
            Flame speed, or flame speed and sampled data

        """
        self.flame_speed = self.run_case()
        if skip_data:
            return self.flame_speed

        temperatures = self.flame.T
        pressures = np.full(len(temperatures), self.flame.P)
        # flame.Y is indexed [species, grid point]; the sampler expects [point, species]
        mass_fractions = self.flame.Y.T

        sampled_data = self._sample_profile(temperatures, pressures, mass_fractions)
        self.sampled_data = sampled_data
        return self.flame_speed, sampled_data
