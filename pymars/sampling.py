"""Module for sampling thermochemical data and global metrics"""
import os
import multiprocessing
import logging
from typing import NamedTuple, Dict

import numpy as np
import yaml
import cantera as ct

from .simulation import Simulation

data_files = {
    'data_ignition': 'ignition_data.dat', 'output_ignition': 'ignition_output.txt',
    'data_psr': 'psr_data.dat', 'output_psr': 'psr_output.txt',
    'data_flame': 'laminarflame_data.dat', 'output_flame': 'laminarflame_output.txt'
    }

class InputIgnition(NamedTuple):
    """Holds input parameters for a single autoignition case.
    """
    kind: str
    temperature: float
    pressure: float

    end_time: float = 0.0
    max_steps: int = 10000

    equivalence_ratio: float = 0.0
    fuel: Dict = {}
    oxidizer: Dict = {}
    reactants: Dict = {}
    composition_type: str = 'mole'


class InputPSR(NamedTuple):
    """Holds input parameters for single PSR simulation.
    """


class InputLaminarFlame(NamedTuple):
    """Holds input parameters for single laminar flame simulation.
    """


def simulation_worker(sim_tuple):
    """Worker for multiprocessing of simulation cases.

    Parameters
    ----------
    sim_tuple : tuple
        Contains Simulation object and other parameters needed to setup
        and run case.

    Returns
    -------
    sim : Simulation
        Object with simulation metadata

    """
    sim, stop_at_ignition = sim_tuple

    sim.setup_case()
    sim.run_case(stop_at_ignition)

    sim = Simulation(sim.idx, sim.properties, sim.model, phase_name=sim.phase_name, path=sim.path)
    return sim


def ignition_worker(sim_tuple):
    """Worker for multiprocessing of ignition delay only cases.

    Parameters
    ----------
    sim_tuple : tuple
        Tuple of Simulation object to be run and identifier

    Returns
    -------
    dict
        Case identifier and calculated ignition delay

    """
    sim, idx = sim_tuple
    sim.setup_case()
    ignition_delay = sim.calculate_ignition()
    return {idx: ignition_delay}


def calculate_error(metrics_original, metrics_test):
    """Calculates error of global metrics between test and original model.

    Parameters
    ----------
    metrics_original : numpy.ndarray
        Metrics serving as basis of error calculation
    metrics_test : numpy.ndarray
        Metrics for which error is being calculated with respect to ``metrics_original``

    Returns
    -------
    error : float
        Maximum error over all metrics
    
    """
    error = 100 * np.max(np.abs(metrics_original - metrics_test) / metrics_original)
    # if any zero ignition delays, set error to 100
    if any(metrics_test == 0.0):
        error = 100.0

    return error


def read_metrics(ignition_conditions, psr_conditions=[], flame_conditions=[]):
    """Reads in stored already-sampled metrics.

    Parameters
    ----------
    ignition_conditions : list of InputIgnition
        List of autoignition initial conditions.
    psr_conditions : list of InputPSR, optional
        List of PSR simulation conditions.
    flame_conditions : list of InputLaminarFlame, optional
        List of laminar flame simulation conditions.

    Returns
    -------
    ignition_delays : numpy.ndarray
        Calculated metrics for model, used for evaluating error

    """
    if ignition_conditions:
        exists_output = os.path.isfile(data_files['output_ignition'])
        if exists_output:
            ignition_delays = np.genfromtxt(data_files['output_ignition'], delimiter=',')
        else:
            raise SystemError('Error, no ignition output file present.')

    if psr_conditions:
        raise NotImplementedError('PSR calculations not currently supported.')
    
    if flame_conditions:
        raise NotImplementedError('Laminar flame calculations not currently supported.')

    return ignition_delays


def sample_metrics(model, ignition_conditions, psr_conditions=[], flame_conditions=[],
                   phase_name='', num_threads=1, path='', reuse_saved=False
                   ):
    """Evaluates metrics used for determining error of reduced model

    Initially, supports autoignition delay only.

    Parameters
    ----------
    model : str
        Filename for Cantera model for performing simulations
    ignition_conditions : list of InputIgnition
        List of autoignition initial conditions.
    psr_conditions : list of InputPSR, optional
        List of PSR simulation conditions.
    flame_conditions : list of InputLaminarFlame, optional
        List of laminar flame simulation conditions.
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    num_threads : int, optional
        Number of CPU threads to use for performing simulations in parallel.
        Optional; default = 1, in which the multiprocessing module is not used.
        If 0, then use the available number of cores minus one. Otherwise,
        use the specified number of threads.
    path : str, optional
        Optional path for writing files
    reuse_saved : bool, optional
        Flag to reuse saved output
    
    Returns
    -------
    ignition_delays : numpy.ndarray
        Calculated metrics for model, used for evaluating error

    """
    # If number of threads specified as 0, use either max number of available
    # cores minus 1, or use 1 if multiple cores not available.
    if not num_threads:
        num_threads = multiprocessing.cpu_count()-1 or 1

    if ignition_conditions:
        ignition_delays = np.zeros(len(ignition_conditions))
        
        exists_output = os.path.isfile(data_files['output_ignition'])
        if reuse_saved and exists_output:
            ignition_delays = np.genfromtxt(data_files['output_ignition'], delimiter=',')
        
        if reuse_saved and len(ignition_delays) == len(ignition_conditions):
            logging.info('Reusing existing autoignition samples for the starting model.')
        else:
            simulations = []
            for idx, case in enumerate(ignition_conditions):
                simulations.append([
                    Simulation(idx, case, model, phase_name=phase_name, path=path), idx
                    ])

            jobs = tuple(simulations)
            if num_threads == 1:
                results = []
                for job in jobs:
                    results.append(ignition_worker(job))
            else:
                pool = multiprocessing.Pool(processes=num_threads)
                results = pool.map(ignition_worker, jobs)
                pool.close()
                pool.join()

            results = {key:val for k in results for key, val in k.items()}
            ignition_delays = np.zeros(len(results))
            for idx, ignition_delay in results.items():
                ignition_delays[idx] = ignition_delay
        
    if psr_conditions:
        raise NotImplementedError('PSR calculations not currently supported.')
    
    if flame_conditions:
        raise NotImplementedError('Laminar flame calculations not currently supported.')

    return ignition_delays


def sample(model, ignition_conditions, psr_conditions=[], flame_conditions=[],
           phase_name='', num_threads=1, path=''
           ):
    """Samples thermochemical data and generates metrics for various phenomena.

    Initially, supports autoignition delay only.

    Parameters
    ----------
    model : str
        Filename for Cantera model for performing simulations
    ignition_conditions : list of InputIgnition
        List of autoignition initial conditions.
    psr_conditions : list of InputPSR, optional
        List of PSR simulation conditions.
    flame_conditions : list of InputLaminarFlame, optional
        List of laminar flame simulation conditions.
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    num_threads : int
        Number of CPU threads to use for performing simulations in parallel.
        Optional; default = 1, in which the multiprocessing module is not used.
        If 0, then use the available number of cores minus one. Otherwise,
        use the specified number of threads.
    path : str, optional
        Optional path for writing files
    
    Returns
    -------
    tuple of numpy.ndarray
        Metrics, and sampled data

    """
    # If number of threads given as 0, use either max number of available
    # cores minus 1, or use 1 if multiple cores not available.
    if not num_threads:
        num_threads = multiprocessing.cpu_count()-1 or 1

    if ignition_conditions:
        ignition_delays = np.zeros(len(ignition_conditions))
        ignition_data = []
        
        # check for presence of data and output files; if present, reuse.
        matches_number = False
        matches_shape = False
        exists_data = os.path.isfile(data_files['data_ignition'])
        exists_output = os.path.isfile(data_files['output_ignition'])
        if exists_data and exists_output:
            ignition_delays = np.genfromtxt(data_files['output_ignition'], delimiter=',')
            ignition_data = np.genfromtxt(data_files['data_ignition'], delimiter=',')
            # need to check that saved data at least matches the number of cases
            matches_number = (
                ignition_delays.size == len(ignition_conditions) and 
                ignition_data.shape[0] / 20 == len(ignition_conditions)
                )
            
            # also check that expected data is right shape (e.g., in case number of species 
            # has changed if running a new model)
            gas = ct.Solution(model, phase_name)
            matches_shape = ignition_data.shape[1] == 2 + gas.n_species
        
        if matches_number and matches_shape:
            logging.info('Reusing existing autoignition samples for the starting model.')
        else:
            logging.info('Running autoignition simulations for starting model.')
            stop_at_ignition = False
            simulations = []
            for idx, case in enumerate(ignition_conditions):
                simulations.append([
                    Simulation(idx, case, model, phase_name=phase_name, path=path), stop_at_ignition
                    ])

            jobs = tuple(simulations)
            if num_threads == 1:
                results = []
                for job in jobs:
                    results.append(simulation_worker(job))
            else:
                pool = multiprocessing.Pool(processes=num_threads)
                results = pool.map(simulation_worker, jobs)
                pool.close()
                pool.join()

            ignition_delays = np.zeros(len(ignition_conditions))
            ignition_data = []     
            for idx, sim in enumerate(results):
                ignition_delays[idx], data = sim.process_results()
                ignition_data += list(data)
                sim.clean()
            ignition_data = np.array(ignition_data)

            np.savetxt(data_files['data_ignition'], ignition_data, delimiter=',')
            np.savetxt(data_files['output_ignition'], ignition_delays, delimiter=',')

    if psr_conditions:
        raise NotImplementedError('PSR calculations not currently supported.')
    
    if flame_conditions:
        raise NotImplementedError('Laminar flame calculations not currently supported.')
    
    return ignition_delays, ignition_data


def parse_ignition_inputs(model, conditions, phase_name=''):
    """Parses input for autoignition simulations, raising an error on any errors.

    Parameters
    ----------
    model : str
        Name of Cantera-format kinetic model
    conditions : dict
        Dictionary with list of autoignition inputs
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    
    Returns
    -------
    list of InputIgnition
        List of validated objects with autoignition input parameters

    """
    gas = ct.Solution(model, phase_name)
    
    inputs = []
    for idx, case in enumerate(conditions):
        pre = f'Ignition input {idx}: '
        
        # check required keys        
        kind = case.get('kind', '')
        temperature = case.get('temperature', 0.0)
        pressure = case.get('pressure', 0.0)

        assert kind in ['constant volume', 'constant pressure'], (
            pre + '"case" needs to be "constant volume" or "constant pressure'
            )
        assert temperature > 0.0, pre + '"temperature" needs to be > 0'
        assert pressure > 0.0, pre + '"pressure" needs to be a number > 0'

        end_time = case.get('end-time', 0)
        max_steps = case.get('max-steps', 10000)
        
        equiv_ratio = case.get('equivalence-ratio', 0.0)
        fuel = case.get('fuel', [])
        oxidizer = case.get('oxidizer', [])

        reactants = case.get('reactants', [])

        assert (bool(equiv_ratio or fuel or oxidizer) + bool(reactants)) == 1, (
            pre + 'should specify either equivalence-ratio/fuel/oxidizer or reactants.'
            )
        
        if equiv_ratio or fuel or oxidizer:
            assert equiv_ratio > 0.0, pre + 'needs non-zero "equivalence-ratio"'

            assert fuel, pre + 'needs "fuel" with at least one entry'
            for entry in fuel:
                assert fuel[entry] > 0, pre + entry + ' value needs to be a number > 0'
                assert entry in gas.species_names, (
                    pre + 'fuel species not in model: ' + entry
                    )
            
            assert oxidizer, pre + 'needs "oxidizer" with at least one entry'
            for entry in oxidizer:
                assert oxidizer[entry] > 0, pre + entry + ' value needs to be a number > 0'
                assert entry in gas.species_names, (
                    pre + 'oxidizer species not in model: ' + entry
                    )
        
        if reactants:
            for entry in reactants:
                assert reactants[entry] > 0, pre + entry + ' value needs to be a number > 0'
                assert entry in gas.species_names, (
                    pre + 'reactant not in model: ' + entry
                    )
        
        composition_type = case.get('composition-type', 'mole')
        assert composition_type in ['mole', 'mass'], pre + 'composition-type must be "mole" or "mass"'
        assert not (composition_type == 'mass' and equiv_ratio), (
            pre + 'composition-type: must be mole when specifying equivalence ratio'
            )
        
        inputs.append(InputIgnition(
            kind, temperature, pressure, end_time, max_steps,
            equiv_ratio, fuel, oxidizer, reactants, composition_type
        ))

    return inputs


def parse_psr_inputs(model, conditions, phase_name=''):
    """Parses input for PSR simulations, raising an error on any errors.

    Parameters
    ----------
    model : str
        Name of Cantera-format kinetic model
    conditions : dict
        Dictionary with list of PSR inputs
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    
    Returns
    -------
    list of InputPSR
        List of validated objects with PSR input parameters

    """
    return None


def parse_flame_inputs(model, conditions, phase_name=''):
    """Parses input for laminar flame simulations, raising an error on any errors.

    Parameters
    ----------
    model : str
        Name of Cantera-format kinetic model
    conditions : dict
        Dictionary with list of laminar flame inputs
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    
    Returns
    -------
    list of InputLaminarFlame
        List of validated objects with laminar flame input parameters

    """
    return None
