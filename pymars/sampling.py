"""Module for sampling thermochemical data and global metrics"""
import os
import multiprocessing
import logging
from typing import NamedTuple

import numpy as np
import yaml

from .simulation import Simulation

class SamplingInputs(NamedTuple):
    """Collects input and output files for sampling (autoignition, PSR, laminar flame)
    """
    input_ignition: str = ''
    data_ignition: str = 'ignition_data.dat'
    output_ignition: str = 'ignition_output.txt'

    input_psr: str = ''
    data_psr: str = 'psr_data.dat'
    output_psr: str = 'psr_output.txt'

    input_laminar_flame: str = ''
    data_laminar_flame: str = 'laminarflame_data.dat'
    output_laminar_flame: str = 'laminarflame_output.txt'


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

    sim = Simulation(sim.idx, sim.properties, sim.model, path=sim.path)
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


def read_metrics(inputs):
    """Reads in stored already-sampled metrics.

    Parameters
    ----------
    inputs : SamplingInputs
        Inputs necessary for sampling

    Returns
    -------
    ignition_delays : numpy.ndarray
        Calculated metrics for model, used for evaluating error

    """
    if inputs.input_ignition:
        exists_output = os.path.isfile(inputs.output_ignition)
        if exists_output:
            ignition_delays = np.genfromtxt(inputs.output_ignition, delimiter=',')
        else:
            raise SystemError('Error, no ignition output file present.')

    elif inputs.input_psr:
        raise NotImplementedError('PSR calculations not currently supported.')
    elif inputs.input_laminar_flame:
        raise NotImplementedError('Laminar flame calculations not currently supported.')
    else:
        raise KeyError('No input files specified for sampling.')

    return ignition_delays


def sample_metrics(inputs, model, num_threads=1, path='', reuse_saved=False):
    """Evaluates metrics used for determining error of reduced model

    Initially, supports autoignition delay only.

    Parameters
    ----------
    inputs : SamplingInputs
        Inputs necessary for sampling
    model : str
        Filename for Cantera model for performing simulations
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

    if inputs.input_ignition:
        with open(inputs.input_ignition, 'r') as the_file:
            conditions = yaml.safe_load(the_file)

        ignition_delays = np.zeros(len(conditions))
        
        exists_output = os.path.isfile(inputs.output_ignition)
        if reuse_saved and exists_output:
            ignition_delays = np.genfromtxt(inputs.output_ignition, delimiter=',')
        
        if reuse_saved and len(ignition_delays) == len(conditions):
            logging.info('Reusing existing autoignition samples for the starting model.')
        else:
            simulations = []
            for idx, properties in enumerate(conditions):
                simulations.append([Simulation(idx, properties, model, path), idx])

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
        
    if inputs.input_psr:
        raise NotImplementedError('PSR calculations not currently supported.')
    if inputs.input_laminar_flame:
        raise NotImplementedError('Laminar flame calculations not currently supported.')

    return ignition_delays


def sample(inputs, model, num_threads=1, path=''):
    """Samples thermochemical data and generates metrics for various phenomena.

    Initially, supports autoignition delay only.

    Parameters
    ----------
    inputs : SamplingInputs
        Inputs necessary for sampling
    model : str
        Filename for Cantera model for performing simulations
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

    if inputs.input_ignition:
        with open(inputs.input_ignition, 'r') as the_file:
            conditions = yaml.safe_load(the_file)
        # TODO: validate input file for correctness.

        ignition_delays = np.zeros(len(conditions))
        ignition_data = []
        
        # check for presence of data and output files; if present, reuse.
        exists_data = os.path.isfile(inputs.data_ignition)
        exists_output = os.path.isfile(inputs.output_ignition)
        if exists_data and exists_output:
            ignition_delays = np.genfromtxt(inputs.output_ignition, delimiter=',')
            ignition_data = np.genfromtxt(inputs.data_ignition, delimiter=',')
            # need to check that saved data at least matches the number of 
        
        if len(ignition_delays) == len(conditions) and len(ignition_data)/20 == len(conditions):
            logging.info('Reusing existing autoignition samples for the starting model.')
        else:
            logging.info('Running autoignition simulations for starting model.')
            stop_at_ignition = False
            simulations = []
            for idx, properties in enumerate(conditions):
                simulations.append([Simulation(idx, properties, model, path), stop_at_ignition])

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

            ignition_delays = np.zeros(len(conditions))
            ignition_data = []     
            for idx, sim in enumerate(results):
                ignition_delays[idx], data = sim.process_results()
                ignition_data += list(data)
                sim.clean()
            ignition_data = np.array(ignition_data)

            np.savetxt(inputs.data_ignition, ignition_data, delimiter=',')
            np.savetxt(inputs.output_ignition, ignition_delays, delimiter=',')

    if inputs.input_psr:
        raise NotImplementedError('PSR calculations not currently supported.')
    
    if inputs.input_laminar_flame:
        raise NotImplementedError('Laminar flame calculations not currently supported.')
    
    return ignition_delays, ignition_data


def check_inputs(inputs):
    """Validates input files for simulations, raising an error if any issues detected.

    Parameters
    ----------
    inputs : SamplingInputs
        Inputs necessary for sampling
    
    Returns
    -------
    bool
        Returns ``True`` if no issues detected, otherwise raises an error.

    """
    if inputs.input_ignition:
        with open(inputs.input_ignition, 'r') as the_file:
            conditions = yaml.safe_load(the_file)

        for idx, case in enumerate(conditions):
            pre = f'Ignition input {idx} '
            # check required keys
            required_keys = ['kind', 'temperature', 'pressure']
            for key in required_keys:
                if key not in case:
                    raise KeyError(pre + 'missing required key ' + key)
            
            if case['kind'] not in ['constant volume', 'constant pressure']:
                raise ValueError(pre + '"case" needs to be "constant volume" or "constant pressure')
            
            if type(case['temperature']) not in [float, int] or case['temperature'] <= 0:
                raise ValueError(pre + 'temperature needs to be a number > 0')
            
            if type(case['pressure']) not in [float, int] or case['pressure'] <= 0:
                raise ValueError(pre + 'pressure needs to be a number > 0')
            
            if 'end-time' in case and case.get('end-time', 0) <= 0:
                raise ValueError(pre + '"end-time" needs to be a number > 0')
            
            equiv_ratio = 'fuel' in case or 'oxidizer' in case or 'equivalence-ratio' in case
            reactants = 'reactants' in case
            if equiv_ratio and reactants:
                raise KeyError(pre + 'should specify either fuel/oxidizer/equivalence ratio or reactants')
            
            if not equiv_ratio and not reactants:
                raise KeyError(pre + 'should specify either fuel/oxidizer/equivalence ratio or reactants')
            
            if equiv_ratio:
                if 'fuel' not in case:
                    raise KeyError(pre + 'needs "fuel"')
                if len(case['fuel']) < 1:
                    raise KeyError(pre + '"fuel" needs at least one entry')
                for entry in case['fuel']:
                    if type(case['fuel'][entry]) not in [float, int] or case['fuel'][entry] <= 0:
                        raise ValueError(pre + entry + ' value needs to be a number > 0')
                
                if 'oxidizer' not in case:
                    raise KeyError(pre + 'needs "oxidizer"')
                if len(case['oxidizer']) < 1:
                    raise KeyError(pre + '"oxidizer" needs at least one entry')
                for entry in case['oxidizer']:
                    if type(case['oxidizer'][entry]) not in [float, int] or case['oxidizer'][entry] <= 0:
                        raise ValueError(pre + entry + ' value needs to be a number > 0')
                
                if 'equivalence-ratio' not in case:
                    raise KeyError(pre + 'needs "equivalence-ratio"')
                if (type(case['equivalence-ratio']) not in [float, int] or 
                    case['equivalence-ratio'] <= 0
                    ):
                    raise ValueError(pre + '"equivalence-ratio" needs to be a number > 0')
            
            if reactants:
                if len(case['reactants']) < 1:
                    raise KeyError(pre + '"reactants" needs at least one entry')
                
                for entry in case['reactants']:
                    if type(case['reactants'][entry]) not in [float, int] or case['reactants'][entry] <= 0:
                        raise ValueError(pre + entry + ' value needs to be a number > 0')

    if inputs.input_psr:
        raise NotImplementedError('PSR calculations not currently supported.')
    
    if inputs.input_laminar_flame:
        raise NotImplementedError('Laminar flame calculations not currently supported.')

    return True
