""" Tests the simulation module used by pyMARS """

import os
from tkinter import Grid
import pkg_resources
from tempfile import TemporaryDirectory

import pytest
import numpy as np
import cantera as ct

from ..sampling import InputLaminarFlame
from ..simulation_flame import FlameSimulation

def relative_location(file):
    file_path = os.path.join(file)
    return pkg_resources.resource_filename(__name__, file_path)


class TestSimulation:
    def test_setup_case_equivalence_ratio(self):
        """Test setting up case that specifies equivalence ratio.
        """
        case = InputLaminarFlame(
            kind= 'constant volume', pressure=1.0, temperature=300, equivalence_ratio=1.0,
            fuel={'CH4': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}, width= 0.014
            )
        sim = FlameSimulation(0, case, 'gri30.cti')
        sim.setup_case()

        assert type(sim.reac) == ct.IdealGasReactor
        assert np.allclose(sim.gas.T, 300)
        assert np.allclose(sim.gas.P, ct.one_atm)

        assert np.allclose(sim.gas.X[sim.gas.species_index('CH4')], 1.0 / (1.0 + 2.0 + 7.52))
        assert np.allclose(sim.gas.X[sim.gas.species_index('O2')], 2.0 / (1.0 + 2.0 + 7.52))
        assert np.allclose(sim.gas.X[sim.gas.species_index('N2')], 7.52 / (1.0 + 2.0 + 7.52))
    
    def test_setup_case_reactants(self):
        """Test setting up case that specifies reactants.
        """
        case = InputLaminarFlame(kind= 'constant volume',
            pressure=1.0, temperature=1000.0,
            reactants={'CH4': 1.0, 'O2': 2.0, 'N2': 7.52}, width = 0.014
            )
        sim = FlameSimulation(0, case, 'gri30.cti')
        sim.setup_case()

        assert type(sim.reac) == ct.IdealGasReactor
        assert np.allclose(sim.gas.T, 1000.0)
        assert np.allclose(sim.gas.P, ct.one_atm)

        assert np.allclose(sim.gas.X[sim.gas.species_index('CH4')], 1.0 / (1.0 + 2.0 + 7.52))
        assert np.allclose(sim.gas.X[sim.gas.species_index('O2')], 2.0 / (1.0 + 2.0 + 7.52))
        assert np.allclose(sim.gas.X[sim.gas.species_index('N2')], 7.52 / (1.0 + 2.0 + 7.52))
    
    def test_setup_case_reactants_mass(self):
        """Test setting up case that specifies reactants using mass fraction.
        """
        case = InputLaminarFlame(
            kind='constant volume', pressure=1.0, temperature=1000.0,
            reactants={'CH4': 0.05518667, 'O2': 0.22014124, 'N2': 0.7246721}, composition_type='mass', width=0.1
            )
        sim = FlameSimulation(0, case, 'gri30.cti')
        sim.setup_case()

        assert type(sim.reac) == ct.IdealGasReactor
        assert np.allclose(sim.gas.T, 1000.0)
        assert np.allclose(sim.gas.P, ct.one_atm)

        assert np.allclose(sim.gas.X[sim.gas.species_index('CH4')], 1.0 / (1.0 + 2.0 + 7.52))
        assert np.allclose(sim.gas.X[sim.gas.species_index('O2')], 2.0 / (1.0 + 2.0 + 7.52))
        assert np.allclose(sim.gas.X[sim.gas.species_index('N2')], 7.52 / (1.0 + 2.0 + 7.52))

    def test_run_case_steady_state(self):
        """Test running a case
        """
        case = InputLaminarFlame(kind= 'constant volume',
            pressure=1, temperature=1000.0, equivalence_ratio=1.0,
            fuel={'CH4': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}, width = 0.014
            )
        
        sim = FlameSimulation(0, case, 'gri30.cti')
        sim.setup_case()

        assert np.allclose(sim.run_case(), 4.036245971957828)

        temperatures = sim.flame.T
        x_position = sim.flame.grid
        velocities = sim.flame.velocity

        final_state = np.array([temperatures[-1], x_position[-1], 
            velocities[-1]])
            
        next_to_final_state = np.array([temperatures[-2],  x_position[-2], 
            velocities[-2]])

        max_state_values = np.maximum(np.zeros(len(final_state)), final_state)
        for row in range(len(temperatures)):
            state = np.array([temperatures[row], x_position[row], 
            velocities[row]])
            max_state_values = np.maximum(max_state_values, state)
        

        residual = np.linalg.norm(
            (final_state - next_to_final_state) / (max_state_values + 1.e-15)
            ) / np.sqrt(sim.sim.n_vars -1)

        assert residual < 1.e-8
    
    """def test_process_results(self):
        Test processing of flame results using artificial data. Currently inoperable.
        
        case = InputLaminarFlame(kind= 'constant volume',
            pressure=1, temperature=1000.0, equivalence_ratio=1.0,
            fuel={'CH4': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}, width=0.014
            )
        with TemporaryDirectory() as temp_dir:
            sim = FlameSimulation(0, case , 'gri30.cti', path=temp_dir)
            sim.setup_case()
            sim.save_file = os.path.join(sim.path, str(sim.idx) + '.h5')

            step_initial = np.arange(0, 10, 0.02)
            temp_initial = 200 * np.ones(len(step_initial))

            # flame (temp = 600) will be at 10.5 cm
            step_ramp = np.arange(0.10, 0.11001, 0.0005)
            temp_ramp = 200 + 800 * (step_ramp - 10)

            step_flat = np.arange(0.11005, 0.15, 0.01)
            temp_flat = 1000 * np.ones(len(step_flat))
            
            steps = np.concatenate((step_initial, step_ramp, step_flat))
            temps = np.concatenate((temp_initial, temp_ramp, temp_flat))

            # add a very small number to account for floating-point roundoff error
            idx = len(temp_initial) + int((len(step_ramp) - 1) / 2)
            temps[idx] += 1e-9
            
            print(sim)
            #sim.flame.T = temps
            #sim.flame.P = np.ones(len(temps))
            #sim.flame.velocities = np.ones(len(temps))
            #sim.flame.grid = steps
            #sim.flame.Y = np.ones(len(temps))
            
            flame_speed, sampled_data = sim.process_results()

            assert np.allclose(flame_speed, 4.036245971957828)

            initial_temp = 200.
            delta = 40.
            for idx in range(20):
                assert np.allclose(sampled_data[idx], [initial_temp + delta, 1, 1, 1])
                delta += 40.     
        """
    def test_clean(self):
        """Test successful cleaning up of data.
        """
        case = InputLaminarFlame(kind= 'constant volume',
            pressure=101325, temperature=1000.0, equivalence_ratio=1.0,
            fuel={'CH4': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}, width = 0.014
            )
        with TemporaryDirectory() as temp_dir:
            sim = FlameSimulation(0, case, 'gri30.cti', path=temp_dir)
            sim.save_file = os.path.join(sim.path, str(sim.idx) + '.h5')
        
        sim.setup_case()
        sim.clean()
        assert not os.path.isfile(sim.save_file)

