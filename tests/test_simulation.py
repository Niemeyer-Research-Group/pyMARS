"""Tests the simulation module used by pyMARS"""

import os
import pathlib
from tempfile import TemporaryDirectory

import pytest
import numpy as np
import h5py
import cantera as ct

from pymars.sampling import InputIgnition
from pymars.simulation import Simulation


def relative_location(file):
    return str(pathlib.Path(__file__).parent / file)


class TestSimulation:
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
        sim = Simulation(0, case, "gri30.yaml")
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
        sim = Simulation(0, case, "gri30.yaml")
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
        sim = Simulation(0, case, "gri30.yaml")
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
            sim = Simulation(0, case, "gri30.yaml", path=temp_dir)
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
            sim = Simulation(0, case, "gri30.yaml", path=temp_dir)
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
            sim = Simulation(0, case, "gri30.yaml", path=temp_dir)
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
            sim = Simulation(0, case, "gri30.yaml", path=temp_dir)
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
            sim = Simulation(0, case, "gri30.yaml", path=temp_dir)
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
            sim = Simulation(0, None, "gri30.yaml", path=temp_dir)

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
            sim = Simulation(0, None, "gri30.yaml", path=temp_dir)

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
            sim = Simulation(0, None, "gri30.yaml", path=temp_dir)
            sim.save_file = os.path.join(sim.path, str(sim.idx) + ".h5")

            with h5py.File(sim.save_file, "w") as h5file:
                grp = h5file.create_group("simulation")
                grp.create_dataset("time", data=np.array([1.0]))
                grp.create_dataset("temperature", data=np.array([1.0]))
                grp.create_dataset("pressure", data=np.array([1.0]))
                grp.create_dataset("mass_fractions", data=np.ones((1, 2)))

            sim.clean()
            assert not os.path.isfile(sim.save_file)
