"""Tests the simulation classes used by pyMARS."""

import os
import pathlib
from tempfile import TemporaryDirectory

import pytest
import numpy as np
import h5py
import cantera as ct

from pymars.sampling import InputIgnition, InputPSR, InputLaminarFlame
from pymars.simulation import (
    BaseSimulation,
    IgnitionSimulation,
    FlameSimulation,
    PSRSimulation,
)


def relative_location(file):
    return str(pathlib.Path(__file__).parent / file)


class TestIgnitionSimulation:
    def test_setup_case_equivalence_ratio(self):
        """Test setting up case that specifies equivalence ratio."""
        case = InputIgnition(
            kind="constant volume",
            pressure=1.0,
            temperature=1000.0,
            equivalence_ratio=1.0,
            fuel={"CH4": 1.0},
            oxidizer={"O2": 1.0, "N2": 3.76},
        )
        sim = IgnitionSimulation(0, case, "gri30.yaml")
        sim.setup_case()

        assert isinstance(sim.reac, ct.IdealGasMoleReactor)
        assert np.allclose(sim.gas.T, 1000.0)
        assert np.allclose(sim.gas.P, ct.one_atm)

        assert np.allclose(
            sim.gas.X[sim.gas.species_index("CH4")], 1.0 / (1.0 + 2.0 + 7.52)
        )
        assert np.allclose(
            sim.gas.X[sim.gas.species_index("O2")], 2.0 / (1.0 + 2.0 + 7.52)
        )
        assert np.allclose(
            sim.gas.X[sim.gas.species_index("N2")], 7.52 / (1.0 + 2.0 + 7.52)
        )

    def test_setup_case_reactants(self):
        """Test setting up case that specifies reactants."""
        case = InputIgnition(
            kind="constant volume",
            pressure=1.0,
            temperature=1000.0,
            reactants={"CH4": 1.0, "O2": 2.0, "N2": 7.52},
        )
        sim = IgnitionSimulation(0, case, "gri30.yaml")
        sim.setup_case()

        assert isinstance(sim.reac, ct.IdealGasMoleReactor)
        assert np.allclose(sim.gas.T, 1000.0)
        assert np.allclose(sim.gas.P, ct.one_atm)

        assert np.allclose(
            sim.gas.X[sim.gas.species_index("CH4")], 1.0 / (1.0 + 2.0 + 7.52)
        )
        assert np.allclose(
            sim.gas.X[sim.gas.species_index("O2")], 2.0 / (1.0 + 2.0 + 7.52)
        )
        assert np.allclose(
            sim.gas.X[sim.gas.species_index("N2")], 7.52 / (1.0 + 2.0 + 7.52)
        )

    def test_setup_case_reactants_mass(self):
        """Test setting up case that specifies reactants using mass fraction."""
        case = InputIgnition(
            kind="constant volume",
            pressure=1.0,
            temperature=1000.0,
            reactants={"CH4": 0.05518667, "O2": 0.22014124, "N2": 0.7246721},
            composition_type="mass",
        )
        sim = IgnitionSimulation(0, case, "gri30.yaml")
        sim.setup_case()

        assert isinstance(sim.reac, ct.IdealGasMoleReactor)
        assert np.allclose(sim.gas.T, 1000.0)
        assert np.allclose(sim.gas.P, ct.one_atm)

        assert np.allclose(
            sim.gas.X[sim.gas.species_index("CH4")], 1.0 / (1.0 + 2.0 + 7.52)
        )
        assert np.allclose(
            sim.gas.X[sim.gas.species_index("O2")], 2.0 / (1.0 + 2.0 + 7.52)
        )
        assert np.allclose(
            sim.gas.X[sim.gas.species_index("N2")], 7.52 / (1.0 + 2.0 + 7.52)
        )

    def test_run_case_steady_state(self):
        """Test running a case without specifying end time."""
        case = InputIgnition(
            kind="constant volume",
            pressure=1.0,
            temperature=1000.0,
            equivalence_ratio=1.0,
            fuel={"CH4": 1.0},
            oxidizer={"O2": 1.0, "N2": 3.76},
        )
        with TemporaryDirectory() as temp_dir:
            sim = IgnitionSimulation(0, case, "gri30.yaml", path=temp_dir)
            sim.setup_case()
            assert np.allclose(sim.run_case(), 1.066766)

            with h5py.File(sim.save_file, "r") as h5file:
                grp = h5file["simulation"]
                temperatures = grp["temperature"][:]
                pressures = grp["pressure"][:]
                mass_fractions = grp["mass_fractions"][:]

            final_state = np.concatenate(
                (np.array([temperatures[-1], pressures[-1]]), mass_fractions[-1])
            )
            next_to_final_state = np.concatenate(
                (np.array([temperatures[-2], pressures[-2]]), mass_fractions[-2])
            )
            max_state_values = np.maximum(np.zeros(len(final_state)), final_state)
            for row in range(len(temperatures)):
                state = np.concatenate(
                    (np.array([temperatures[row], pressures[row]]), mass_fractions[row])
                )
                max_state_values = np.maximum(max_state_values, state)

            residual = np.linalg.norm(
                (final_state - next_to_final_state) / (max_state_values + 1.0e-15)
            ) / np.sqrt(sim.sim.n_vars - 1)
            assert residual < 1.0e-8

    def test_run_case_end_time(self):
        """Test running a case with a specified end time."""
        case = InputIgnition(
            kind="constant volume",
            pressure=1.0,
            temperature=1000.0,
            equivalence_ratio=1.0,
            fuel={"CH4": 1.0},
            oxidizer={"O2": 1.0, "N2": 3.76},
            end_time=2.0,
        )
        with TemporaryDirectory() as temp_dir:
            sim = IgnitionSimulation(0, case, "gri30.yaml", path=temp_dir)
            sim.setup_case()
            assert np.allclose(sim.run_case(), 1.066766)
            assert sim.sim.time >= case.end_time

    def test_run_case_noignition(self):
        """Test running a case with no ignition."""
        case = InputIgnition(
            kind="constant volume",
            pressure=1.0,
            temperature=1000.0,
            reactants={"N2": 1.0},
        )
        with TemporaryDirectory() as temp_dir:
            sim = IgnitionSimulation(0, case, "gri30.yaml", path=temp_dir)
            sim.setup_case()
            with pytest.raises(RuntimeError) as excinfo:
                sim.run_case()
                assert "No ignition detected for integration case 0" in str(
                    excinfo.value
                )

        case = InputIgnition(
            kind="constant volume",
            pressure=1.0,
            temperature=1000.0,
            reactants={"N2": 1.0},
            end_time=1.0,
        )
        with TemporaryDirectory() as temp_dir:
            sim = IgnitionSimulation(0, case, "gri30.yaml", path=temp_dir)
            sim.setup_case()
            with pytest.raises(RuntimeError) as excinfo:
                sim.run_case()
                assert "No ignition detected for integration case 0" in str(
                    excinfo.value
                )

        case = InputIgnition(
            kind="constant volume",
            pressure=1.0,
            temperature=1000.0,
            reactants={"N2": 1.0},
            max_steps=1,
        )
        with TemporaryDirectory() as temp_dir:
            sim = IgnitionSimulation(0, case, "gri30.yaml", path=temp_dir)
            sim.setup_case()
            with pytest.raises(RuntimeError) as excinfo:
                sim.run_case()
                assert (
                    "Maximum number of steps reached before convergence for integration case 0"
                    in str(excinfo.value)
                )

    def test_process_results(self):
        """Test processing of ignition results using artificial data."""
        with TemporaryDirectory() as temp_dir:
            sim = IgnitionSimulation(0, None, "gri30.yaml", path=temp_dir)

            sim.save_file = os.path.join(sim.path, str(sim.idx) + ".h5")

            time_initial = np.arange(0, 10, 0.02)
            temp_initial = 200 * np.ones(len(time_initial))

            # ignition delay (temp = 600) will be at 10.5 s
            time_ramp = np.arange(10, 11.001, 0.005)
            temp_ramp = 200 + 800 * (time_ramp - 10)

            time_flat = np.arange(11.005, 15, 0.01)
            temp_flat = 1000 * np.ones(len(time_flat))

            times = np.concatenate((time_initial, time_ramp, time_flat))
            temps = np.concatenate((temp_initial, temp_ramp, temp_flat))

            # add a very small number to account for floating-point roundoff error
            idx = len(temp_initial) + int((len(time_ramp) - 1) / 2)
            temps[idx] += 1e-9

            with h5py.File(sim.save_file, "w") as h5file:
                grp = h5file.create_group("simulation")
                grp.create_dataset("time", data=times)
                grp.create_dataset("temperature", data=temps)
                grp.create_dataset("pressure", data=np.ones(len(times)))
                grp.create_dataset("mass_fractions", data=np.ones((len(times), 2)))

            ignition_delay, sampled_data = sim.process_results()

            assert np.allclose(ignition_delay, 10.5)

            initial_temp = 200.0
            delta = 40.0
            for idx in range(20):
                assert np.allclose(sampled_data[idx], [initial_temp + delta, 1, 1, 1])
                delta += 40.0

    def test_process_results_skip_data(self):
        """Test processing of ignition results, skipping data sampling, using artificial data."""
        with TemporaryDirectory() as temp_dir:
            sim = IgnitionSimulation(0, None, "gri30.yaml", path=temp_dir)

            sim.save_file = os.path.join(sim.path, str(sim.idx) + ".h5")

            time_initial = np.arange(0, 10, 0.02)
            temp_initial = 200 * np.ones(len(time_initial))

            # ignition delay (temp = 600) will be at 10.5 s
            time_ramp = np.arange(10, 11.02, 0.02)
            temp_ramp = 200 + 800 * (time_ramp - 10)

            time_flat = np.arange(11.02, 15, 0.02)
            temp_flat = 1000 * np.ones(len(time_flat))

            times = np.concatenate((time_initial, time_ramp, time_flat))
            temps = np.concatenate((temp_initial, temp_ramp, temp_flat))

            # add a very small number to account for floating-point roundoff error
            idx = len(temp_initial) + 25
            temps[idx] += 1e-8

            with h5py.File(sim.save_file, "w") as h5file:
                grp = h5file.create_group("simulation")
                grp.create_dataset("time", data=times)
                grp.create_dataset("temperature", data=temps)
                grp.create_dataset("pressure", data=np.ones(len(times)))
                grp.create_dataset("mass_fractions", data=np.ones((len(times), 2)))

            sim.process_results(skip_data=True)

            assert np.allclose(sim.ignition_delay, 10.5)
            assert not hasattr(sim, "sampled_data")

    def test_clean(self):
        """Test successful cleaning up of data."""
        with TemporaryDirectory() as temp_dir:
            sim = IgnitionSimulation(0, None, "gri30.yaml", path=temp_dir)
            sim.save_file = os.path.join(sim.path, str(sim.idx) + ".h5")

            with h5py.File(sim.save_file, "w") as h5file:
                grp = h5file.create_group("simulation")
                grp.create_dataset("time", data=np.array([1.0]))
                grp.create_dataset("temperature", data=np.array([1.0]))
                grp.create_dataset("pressure", data=np.array([1.0]))
                grp.create_dataset("mass_fractions", data=np.ones((1, 2)))

            sim.clean()
            assert not os.path.isfile(sim.save_file)


class TestIgnitionFailure:
    """An ignition integration that fails should be handled gracefully rather
    than crashing the reduction (see issue #69).

    Two paths are covered: the metric-only ``calculate`` path (used for candidate
    reduced models) must degrade to a 0.0 ignition delay so the candidate is
    rejected through the error metric, while the ``run_case`` path (used for the
    original/baseline model during a sampling run) should raise so a broken
    baseline is caught. The integrator failure is forced deterministically by
    monkeypatching ``_step`` to raise ``ct.CanteraError`` (as CVODES does for
    non-finite derivatives); a real non-igniting (inert) mixture covers the
    no-ignition case.
    """

    def _case(self):
        return InputIgnition(
            kind="constant volume",
            pressure=1.0,
            temperature=1000.0,
            equivalence_ratio=1.0,
            fuel={"H2": 1.0},
            oxidizer={"O2": 1.0, "N2": 3.76},
        )

    def _inert_case(self):
        # pure N2 never ignites; small max_steps keeps the steady-state loop short
        return InputIgnition(
            kind="constant volume",
            pressure=1.0,
            temperature=1000.0,
            reactants={"N2": 1.0},
            max_steps=10,
        )

    @staticmethod
    def _force_integration_failure(_self):
        raise ct.CanteraError("forced failure: CVODES error encountered")

    def test_calculate_returns_zero_on_integration_failure(self, monkeypatch):
        """The metric-only path degrades to 0.0 (so the candidate is rejected)."""
        monkeypatch.setattr(
            IgnitionSimulation, "_step", self._force_integration_failure
        )
        with TemporaryDirectory() as temp_dir:
            sim = IgnitionSimulation(0, self._case(), "h2o2.yaml", path=temp_dir)
            sim.setup_case()
            # must not raise; returns 0.0
            assert sim.calculate() == 0.0
            assert sim.ignition_delay == 0.0

    def test_run_case_raises_on_integration_failure(self, monkeypatch):
        """The original-model sampling path raises, so a broken baseline is caught."""
        monkeypatch.setattr(
            IgnitionSimulation, "_step", self._force_integration_failure
        )
        with TemporaryDirectory() as temp_dir:
            sim = IgnitionSimulation(0, self._case(), "h2o2.yaml", path=temp_dir)
            sim.setup_case()
            with pytest.raises(RuntimeError, match="Integration failed"):
                sim.run_case()

    def test_calculate_returns_zero_on_no_ignition(self):
        """A real non-igniting (inert) mixture yields 0.0 from the metric path."""
        with TemporaryDirectory() as temp_dir:
            sim = IgnitionSimulation(0, self._inert_case(), "h2o2.yaml", path=temp_dir)
            sim.setup_case()
            # must not raise; returns 0.0
            assert sim.calculate() == 0.0
            assert sim.ignition_delay == 0.0


class TestFlameSetupCase:
    """The gas state and composition should be initialized correctly."""

    def test_setup_case_equivalence_ratio(self):
        case = InputLaminarFlame(
            pressure=1.0,
            temperature=300.0,
            equivalence_ratio=1.0,
            fuel={"CH4": 1.0},
            oxidizer={"O2": 1.0, "N2": 3.76},
            width=0.03,
        )
        sim = FlameSimulation(0, case, "gri30.yaml")
        sim.setup_case()

        assert isinstance(sim.flame, ct.FreeFlame)
        assert np.allclose(sim.gas.T, 300.0)
        assert np.allclose(sim.gas.P, ct.one_atm)

        assert np.allclose(
            sim.gas.X[sim.gas.species_index("CH4")], 1.0 / (1.0 + 2.0 + 7.52)
        )
        assert np.allclose(
            sim.gas.X[sim.gas.species_index("O2")], 2.0 / (1.0 + 2.0 + 7.52)
        )
        assert np.allclose(
            sim.gas.X[sim.gas.species_index("N2")], 7.52 / (1.0 + 2.0 + 7.52)
        )

    def test_setup_case_reactants(self):
        case = InputLaminarFlame(
            pressure=1.0,
            temperature=300.0,
            reactants={"CH4": 1.0, "O2": 2.0, "N2": 7.52},
            width=0.03,
        )
        sim = FlameSimulation(0, case, "gri30.yaml")
        sim.setup_case()

        assert isinstance(sim.flame, ct.FreeFlame)
        assert np.allclose(sim.gas.T, 300.0)
        assert np.allclose(sim.gas.P, ct.one_atm)
        assert np.allclose(
            sim.gas.X[sim.gas.species_index("CH4")], 1.0 / (1.0 + 2.0 + 7.52)
        )

    def test_setup_case_reactants_mass(self):
        case = InputLaminarFlame(
            pressure=1.0,
            temperature=300.0,
            reactants={"CH4": 0.05518667, "O2": 0.22014124, "N2": 0.7246721},
            composition_type="mass",
            width=0.03,
        )
        sim = FlameSimulation(0, case, "gri30.yaml")
        sim.setup_case()

        assert isinstance(sim.flame, ct.FreeFlame)
        assert np.allclose(
            sim.gas.X[sim.gas.species_index("CH4")], 1.0 / (1.0 + 2.0 + 7.52)
        )


class TestFlameRun:
    """A small hydrogen flame should solve and produce sensible results."""

    def _hydrogen_case(self):
        return InputLaminarFlame(
            pressure=1.0,
            temperature=300.0,
            equivalence_ratio=1.0,
            fuel={"H2": 1.0},
            oxidizer={"O2": 1.0, "N2": 3.76},
            width=0.03,
        )

    def test_run_case_flame_speed(self):
        sim = FlameSimulation(0, self._hydrogen_case(), "h2o2.yaml")
        sim.setup_case()
        flame_speed = sim.run_case()

        # Stoichiometric H2/air laminar flame speed is roughly 2-3 m/s
        assert 0.5 < flame_speed < 5.0

    def test_process_results_shape(self):
        sim = FlameSimulation(0, self._hydrogen_case(), "h2o2.yaml")
        sim.setup_case()
        flame_speed, sampled_data = sim.process_results()

        gas = ct.Solution("h2o2.yaml")
        assert 0.5 < flame_speed < 5.0
        assert sampled_data.shape == (20, 2 + gas.n_species)
        # temperatures should be monotonically non-decreasing across the profile
        assert np.all(np.diff(sampled_data[:, 0]) >= 0)

    def test_calculate_only(self):
        sim = FlameSimulation(0, self._hydrogen_case(), "h2o2.yaml")
        sim.setup_case()
        assert 0.5 < sim.calculate() < 5.0

    def test_run_case_mass_composition(self):
        """A flame specified by mass fractions should also solve."""
        case = InputLaminarFlame(
            pressure=1.0,
            temperature=300.0,
            # stoichiometric H2/air by mass fraction
            reactants={"H2": 0.0285, "O2": 0.2265, "N2": 0.7450},
            composition_type="mass",
            width=0.03,
        )
        sim = FlameSimulation(0, case, "h2o2.yaml")
        sim.setup_case()
        assert 0.5 < sim.run_case() < 5.0


class TestFlameFailure:
    """A flame that fails to solve, or solves to a degenerate (non-physical)
    result, should be handled gracefully rather than crashing the reduction.

    Two failure modes are covered: a solver error (forced deterministically by
    monkeypatching the solve to raise ``ct.CanteraError``), and a real solve of a
    non-flammable mixture, which converges to a near-zero "flame speed" below
    ``min_flame_speed`` rather than raising.
    """

    def _case(self):
        return InputLaminarFlame(
            pressure=1.0,
            temperature=300.0,
            equivalence_ratio=1.0,
            fuel={"H2": 1.0},
            oxidizer={"O2": 1.0, "N2": 3.76},
            width=0.03,
        )

    def _inert_case(self):
        # pure N2: a real flame solve converges to a degenerate ~0.024 m/s speed,
        # which is below the min_flame_speed floor and so counts as "no flame".
        return InputLaminarFlame(
            pressure=1.0,
            temperature=300.0,
            reactants={"N2": 1.0},
            width=0.03,
        )

    @staticmethod
    def _force_no_flame(_self):
        raise ct.CanteraError("forced failure: no flame detected")

    def test_calculate_returns_zero_on_no_flame(self, monkeypatch):
        """The metric-only path degrades to 0.0 (so the model is later rejected)."""
        monkeypatch.setattr(FlameSimulation, "_solve_flame", self._force_no_flame)
        sim = FlameSimulation(0, self._case(), "h2o2.yaml")
        sim.setup_case()
        # must not raise; returns 0.0
        assert sim.calculate() == 0.0
        assert sim.flame_speed == 0.0

    def test_run_case_raises_on_no_flame(self, monkeypatch):
        """The original-model data path raises, so a missing baseline is caught."""
        monkeypatch.setattr(FlameSimulation, "_solve_flame", self._force_no_flame)
        sim = FlameSimulation(0, self._case(), "h2o2.yaml")
        sim.setup_case()
        with pytest.raises(RuntimeError, match="No flame detected"):
            sim.run_case()

    def test_process_results_raises_on_no_flame(self, monkeypatch):
        """process_results follows run_case, so it also raises on a missing flame."""
        monkeypatch.setattr(FlameSimulation, "_solve_flame", self._force_no_flame)
        sim = FlameSimulation(0, self._case(), "h2o2.yaml")
        sim.setup_case()
        with pytest.raises(RuntimeError, match="No flame detected"):
            sim.process_results()

    def test_calculate_returns_zero_on_degenerate_flame(self):
        """A real solve of a non-flammable (pure inert) mixture yields a degenerate
        near-zero speed, treated as no flame -> 0.0 from the metric-only path."""
        sim = FlameSimulation(0, self._inert_case(), "h2o2.yaml")
        sim.setup_case()
        # sanity-check the premise: the raw solve gives a degenerate speed below
        # the floor (so this exercises the degenerate path, not a solver error)
        raw_speed = sim._solve_flame()
        assert 0.0 < raw_speed < sim.min_flame_speed

        assert sim.calculate() == 0.0
        assert sim.flame_speed == 0.0

    def test_run_case_raises_on_degenerate_flame(self):
        """The original-model path raises on the same degenerate (no-flame) result."""
        sim = FlameSimulation(0, self._inert_case(), "h2o2.yaml")
        sim.setup_case()
        with pytest.raises(RuntimeError, match="No flame detected"):
            sim.run_case()

    def test_min_flame_speed_override_accepts_low_speed(self):
        """Lowering the floor lets an otherwise-degenerate speed count as a flame.

        The default floor (0.05) would reject the inert ~0.024 m/s result; a floor
        of 0.0 accepts any positive speed. Fresh simulations are used for each call
        because re-solving an already-solved degenerate flame drifts to <= 0.
        """
        sim = FlameSimulation(0, self._inert_case(), "h2o2.yaml", min_flame_speed=0.0)
        sim.setup_case()
        assert sim.min_flame_speed == 0.0
        assert sim.run_case() > 0.0  # must not raise; positive speed accepted

        sim2 = FlameSimulation(1, self._inert_case(), "h2o2.yaml", min_flame_speed=0.0)
        sim2.setup_case()
        assert sim2.calculate() > 0.0


class TestSampleProfile:
    """Exercises the shared profile sampler used by every simulation type."""

    def test_sampling(self):
        """Each sampled row is the first grid point at or above the corresponding
        fraction-of-rise threshold, taking T, P, and mass fractions from that
        same point."""
        # Fine linear temperature ramp 300 -> 2300 K (1 K per grid step), so each
        # grid point crosses at most one of the 20 evenly-spaced thresholds.
        n = 2001
        temperatures = np.linspace(300.0, 2300.0, n)
        pressures = np.full(n, 2.0)
        # Row i = [i, i] so the originating grid point is recoverable from a row.
        mass_fractions = np.tile(np.arange(n).reshape(-1, 1), (1, 2)).astype(float)

        data = BaseSimulation._sample_profile(temperatures, pressures, mass_fractions)

        num = BaseSimulation.num_sample_points
        n_species = mass_fractions.shape[1]
        assert data.shape == (num, 2 + n_species)

        # Reproduce the thresholds exactly as ``_sample_profile`` computes them,
        # then find the first grid point at or above each. Because this is the
        # same selection the sampler performs, the expected values are exact grid
        # values and the comparison needs no tolerance.
        delta = 1.0 / num
        deltas = np.arange(delta, 1 + delta, delta)
        thresholds = temperatures[0] + deltas * (temperatures[-1] - temperatures[0])
        expected_idx = np.searchsorted(temperatures, thresholds, side="left")

        # Grid is fine enough that each threshold maps to a distinct grid point;
        # this is what makes the sampler's one-point-per-threshold scan agree with
        # an independent ``searchsorted`` lookup.
        assert np.all(np.diff(expected_idx) > 0)

        for k, point in enumerate(expected_idx):
            if point < n:
                assert data[k, 0] == temperatures[point]  # temperature
                assert data[k, 1] == pressures[point]  # pressure
                assert np.array_equal(data[k, 2:], mass_fractions[point])  # Y
            else:
                # a threshold beyond the maximum temperature leaves the row unset
                assert np.array_equal(data[k], np.zeros(2 + n_species))


class TestPSRRun:
    """Exercises the perfectly stirred reactor (PSR) simulation."""

    def _case(self):
        return InputPSR(
            temperature=300.0,
            pressure=1.0,
            equivalence_ratio=1.0,
            fuel={"CH4": 1.0},
            oxidizer={"O2": 1.0, "N2": 3.76},
        )

    def test_run_case_metrics(self):
        """run_case returns the [tau_ext, T_mid, T_near] metric vector."""
        sim = PSRSimulation(0, self._case(), "gri30.yaml")
        sim.setup_case()
        metrics = sim.run_case()

        assert metrics.shape == (3,)
        tau_ext, temp_mid, temp_near = metrics
        # stoichiometric CH4/air, 1 atm, 300 K (GRI-Mech 3.0); the solver is
        # deterministic, so these are tight regression bounds with a small margin
        # for cross-version kinetics/solver differences
        assert tau_ext == pytest.approx(7.90e-5, rel=0.05)
        assert temp_mid == pytest.approx(2062.0, rel=0.03)
        assert temp_near == pytest.approx(2207.0, rel=0.03)

    def test_calculate_matches_run_case(self):
        """The metric-only path matches the full run_case metrics."""
        sim = PSRSimulation(0, self._case(), "gri30.yaml")
        sim.setup_case()
        assert np.allclose(sim.calculate(), sim.run_case())

    def test_process_results_shape(self):
        """process_results samples the three points' states."""
        gas = ct.Solution("gri30.yaml")
        sim = PSRSimulation(0, self._case(), "gri30.yaml")
        sim.setup_case()
        metrics, data = sim.process_results()

        assert metrics.shape == (3,)
        assert data.shape == (PSRSimulation.num_sample_points, 2 + gas.n_species)
        # sampled pressure is the (constant) inlet pressure in Pa
        assert np.allclose(data[:, 1], 1.0 * ct.one_atm)

    def test_regression_gri_ch4(self):
        """Documented regression anchor: stoichiometric CH4/air, 300 K, 1 atm, GRI 3.0."""
        sim = PSRSimulation(0, self._case(), "gri30.yaml")
        sim.setup_case()
        tau_ext = sim.run_case()[0]
        assert tau_ext == pytest.approx(7.9e-5, rel=0.15)


class TestPSRFailure:
    """A PSR whose extinction curve cannot be traced is handled gracefully.

    The integrator failure is forced deterministically by monkeypatching the
    continuation solver to raise; the metric-only ``calculate`` path degrades to a
    zero vector (so the candidate is rejected via the error metric) while the
    baseline ``run_case`` path raises.
    """

    def _case(self):
        return InputPSR(
            temperature=300.0,
            pressure=1.0,
            equivalence_ratio=1.0,
            fuel={"CH4": 1.0},
            oxidizer={"O2": 1.0, "N2": 3.76},
        )

    @staticmethod
    def _force_failure(_gas):
        raise ct.CanteraError("forced failure: no extinction curve")

    def test_calculate_returns_zeros_on_failure(self, monkeypatch):
        """The metric-only path degrades to zeros (so the candidate is rejected)."""
        monkeypatch.setattr(
            "pymars.simulation.trace_extinction_curve", self._force_failure
        )
        sim = PSRSimulation(0, self._case(), "gri30.yaml")
        sim.setup_case()
        metrics = sim.calculate()

        assert metrics.shape == (3,)
        assert np.all(metrics == 0.0)

    def test_run_case_raises_on_failure(self, monkeypatch):
        """The original-model path raises, so a broken baseline is caught."""
        monkeypatch.setattr(
            "pymars.simulation.trace_extinction_curve", self._force_failure
        )
        sim = PSRSimulation(0, self._case(), "gri30.yaml")
        sim.setup_case()
        with pytest.raises(RuntimeError, match="No extinction curve detected"):
            sim.run_case()
