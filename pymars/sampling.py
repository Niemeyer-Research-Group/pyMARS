"""Module for sampling thermochemical data and global metrics"""
import os
import multiprocessing
import logging
from typing import NamedTuple, Dict

import numpy as np
import cantera as ct

from .simulation import Simulation

data_files = {
    'data_ignition': 'ignition_data.dat', 'output_ignition': 'ignition_output.txt', 'weights_ignition': 'ignition_weights.txt',
    'data_psr': 'psr_data.dat', 'output_psr': 'psr_output.txt', 'weights_psr': 'psr_weights.txt',
    'data_flame': 'laminarflame_data.dat', 'output_flame': 'laminarflame_output.txt', 'weights_flame': 'laminarflame_weights.txt'
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

    target_weights: list = []


class InputPSR(NamedTuple):
    """Holds input parameters for single PSR simulation.
    """


class InputLaminarFlame(NamedTuple):
    """Holds input parameters for single laminar flame simulation.
    """
    kind: str
    temperature: float
    pressure: float

    transport: str

    width: float
    refine_ratio: float
    refine_slope: float
    refine_curve: float
    refine_prune: float
    
    equivalence_ratio: float = 0.0
    fuel: Dict = {}
    oxidizer: Dict = {}
    reactants: Dict = {}
    composition_type: str = 'mole'

    target_weights: list = []


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
    sim, sim_parameters = sim_tuple

    sim.setup_case()
    sim.run_case(sim_parameters)

    sim = Simulation(sim.sim_type, sim.idx, sim.conditions, sim.model, phase_name=sim.phase_name, path=sim.path)
    return sim


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
    # if any zero ignition delays, set error to 100
    if any(metrics_test == 0.0):
        logging.warning(f"Detected zero metric, metrics_test = {metrics_test}")
        return 100.0
    
    error = 100 * np.max(np.abs(metrics_original - metrics_test) / metrics_original)

    return error


def read_metrics(conditions, sim_type, psr_conditions=[], flame_conditions=[]):
    """Reads in stored already-sampled metrics.

    Parameters
    ----------
    conditions : list of InputIgnition or InputLaminarFlame
        List of autoignition initial conditions.
    sim_type : str
        'ignition' or 'flame'
    psr_conditions : list of InputPSR, optional
        List of PSR simulation conditions.
    flame_conditions : list of InputLaminarFlame, optional
        List of laminar flame simulation conditions.

    Returns
    -------
    ignition_delays : numpy.ndarray
        Calculated metrics for model, used for evaluating error

    """
    if conditions:
        exists_output = os.path.isfile(data_files['output_' + sim_type])
        if exists_output:
            ignition_delays = np.genfromtxt(data_files['output_' + sim_type], delimiter=',')
        else:
            raise SystemError('Error, no ignition output file present.')

    if psr_conditions:
        raise NotImplementedError('PSR calculations not currently supported.')
    
    if flame_conditions:
        raise NotImplementedError('Laminar flame calculations not currently supported.')

    return ignition_delays


def run_simulations(model, sim_type, conditions, phase_name='', num_threads=1,
                    path='', save_data=False):
    """Runs the batch of simulations

    Parameters
    ----------
    model : str
        Filename for Cantera model for performing simulations
    sim_type : str
        'ignition' or 'flame'
    conditions : list of InputIgnition or InputLaminarFlame
        List of autoignition initial conditions.
    phase_name : str, optional
        Optional name for phase (e.g., 'gas'). 
    num_threads : int
        Number of CPU threads to use for performing simulations in parallel.
        Optional; default = 1, in which the multiprocessing module is not used.
        If 0, then use the available number of cores minus one. Otherwise,
        use the specified number of threads.
    path : str, optional
        Optional path for writing files
    save_data : Boolean
        Choose whether to save the 'data_*' files or not
    
    Returns
    -------
    tuple of numpy.ndarray
        Metrics, and sampled data

    """
    sim_metrics = np.zeros(len(conditions))
    sim_data = []
    sim_weights = []
    
    # check for presence of data and output files; if present, reuse.
    matches_number = False
    matches_shape = False
    exists_data = os.path.isfile(os.path.join(path,data_files['data_' + sim_type]))
    exists_output = os.path.isfile(os.path.join(path,data_files['output_' + sim_type]))
    exists_weights = os.path.isfile(os.path.join(path,data_files['weights_' + sim_type]))
    if exists_data and exists_output and exists_weights:
        sim_metrics = np.genfromtxt(os.path.join(path,data_files['output_' + sim_type]), delimiter=',')
        sim_data = np.genfromtxt(os.path.join(path,data_files['data_' + sim_type]), delimiter=',')
        sim_weights = np.genfromtxt(os.path.join(path,data_files['weights_' + sim_type]), delimiter=',')
        # need to check that saved data at least matches the number of cases
        matches_number = (
            sim_metrics.size == len(conditions) and 
            sim_data.shape[0] / 20 == len(conditions)
            )
        
        # also check that expected data is right shape (e.g., in case number of species 
        # has changed if running a new model)
        gas = ct.Solution(os.path.join(path,model), phase_name)
        matches_shape = sim_data.shape[1] == 2 + gas.n_species
    
    if matches_number and matches_shape:
        logging.info(f"Reusing existing {sim_type} samples for the starting model {model}.")
    else:
        sim_metrics = np.zeros(len(conditions))
        sim_data = []
        sim_weights = []

        if sim_type == 'ignition':
            stop_at_ignition = False
            simulations = []
            for idx, case in enumerate(conditions):
                simulations.append([
                    Simulation(sim_type, idx, case, model, phase_name=phase_name, path=path), stop_at_ignition
                    ])
        elif sim_type == 'flame':
            simulations = []
            for idx, case in enumerate(conditions):
                simulations.append([
                    Simulation(sim_type, idx, case, model, phase_name=phase_name, path=path), None
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

        for idx, sim in enumerate(results):
            sim_metrics[idx], data, sim_weight = sim.process_results()
            sim_data += list(data)
            sim_weights += list(sim_weight)
            sim.clean()
        sim_data = np.array(sim_data)
        sim_weights = np.array(sim_weights)

        if save_data:
            np.savetxt(data_files['data_' + sim_type], sim_data, delimiter=',')
            np.savetxt(data_files['output_' + sim_type], sim_metrics, delimiter=',')
            np.savetxt(data_files['weights_' +  sim_type], sim_weights, delimiter=',')

    return sim_metrics, sim_data, sim_weights
    

def sample_metrics(model, ignition_conditions=[], psr_conditions=[], flame_conditions=[],
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

    metrics = []

    if ignition_conditions:
        [ignition_delays, _, _] = run_simulations(model,'ignition',ignition_conditions,phase_name=phase_name,num_threads=num_threads,path=path)
        metrics.append(ignition_delays)
        
    if psr_conditions:
        raise NotImplementedError('PSR calculations not currently supported.')
    
    if flame_conditions:
        [flame_speeds, _, _] = run_simulations(model,'flame',flame_conditions,phase_name=phase_name,num_threads=num_threads,path=path)
        metrics.append(flame_speeds)

    metrics = np.concatenate(metrics)

    return metrics


def sample(model, ignition_conditions=[], psr_conditions=[], flame_conditions=[],
           phase_name='', num_threads=1, path=''
           ):
    """Samples thermochemical data and generates metrics for various phenomena.

    Currently supports autoignition and laminar flame speed only.

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
    metrics: numpy.ndarray
        Metrics to be used to assess the model
    sim_data: numpy.ndarray
        Sampled data to use in reduction algorithm
    target_weights: list
        Species weights for reduction algorithm

    """
    # If number of threads given as 0, use either max number of available
    # cores minus 1, or use 1 if multiple cores not available.
    if not num_threads:
        num_threads = multiprocessing.cpu_count()-1 or 1

    metrics = []
    sim_data = []
    sim_weights = []

    if ignition_conditions:
        [ignition_delays, ignition_data, ignition_weights] = run_simulations(model,'ignition',ignition_conditions,phase_name,num_threads,path,save_data=True)
        metrics.append(ignition_delays)
        sim_data.append(ignition_data)
        sim_weights.append(ignition_weights)

    if psr_conditions:
        raise NotImplementedError('PSR calculations not currently supported.')
    
    if flame_conditions:
        [flame_speeds, flame_data, flame_weights] = run_simulations(model,'flame',flame_conditions,phase_name,num_threads,path,save_data=True)
        metrics.append(flame_speeds)
        sim_data.append(flame_data)
        sim_weights.append(flame_weights)
        
    metrics = np.concatenate(metrics)
    sim_data = np.concatenate(sim_data, axis=0)
    sim_weights = np.concatenate(sim_weights, axis=0)

    return metrics, sim_data, sim_weights


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

        target_weights = case.get('target-weights', [])

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
            equiv_ratio, fuel, oxidizer, reactants, composition_type, target_weights
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
    gas = ct.Solution(model, phase_name)
    
    inputs = []
    for idx, case in enumerate(conditions):
        pre = f'Ignition input {idx}: '
        
        # check required keys        
        kind = case.get('kind', '')
        temperature = case.get('temperature', 0.0)
        pressure = case.get('pressure', 0.0)

        assert kind in ['premixed'], (
            pre + '"case" needs to be "premixed", other types of 1D flames not supported yet'
            )
        assert temperature > 0.0, pre + '"temperature" needs to be > 0'
        assert pressure > 0.0, pre + '"pressure" needs to be a number > 0'

        transport = case.get('transport', 'mixture-averaged')
        
        width = case.get('width', 0.1)
        refine_ratio = case.get('refine_ratio', 10)
        refine_slope = case.get('refine_slope', 0.15)
        refine_curve = case.get('refine_curve', 0.3)
        refine_prune = case.get('refine_prune', 0)
        
        equiv_ratio = case.get('equivalence-ratio', None)
        assert equiv_ratio > 0.0, pre + '"equivalence-ratio" needs to be > 0'

        fuel = case.get('fuel', [])
        oxidizer = case.get('oxidizer', [])

        reactants = case.get('reactants', [])

        target_weights = case.get('target-weights', [])

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
        
        inputs.append(InputLaminarFlame(
            kind, temperature, pressure, transport, width, 
            refine_ratio, refine_slope, refine_curve, refine_prune,
            equiv_ratio, fuel, oxidizer, reactants, composition_type, target_weights
        ))
    return inputs
