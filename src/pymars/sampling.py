"""Module for sampling thermochemical data and global metrics"""

import os
import multiprocessing
import logging
from typing import NamedTuple, Dict

import numpy as np
import cantera as ct

from .simulation import IgnitionSimulation, PSRSimulation, FlameSimulation

data_files = {
    "data_ignition": "ignition_data.dat",
    "output_ignition": "ignition_output.txt",
    "data_psr": "psr_data.dat",
    "output_psr": "psr_output.txt",
    "data_flame": "laminarflame_data.dat",
    "output_flame": "laminarflame_output.txt",
}


class InputIgnition(NamedTuple):
    """Holds input parameters for a single autoignition case."""

    kind: str
    temperature: float
    pressure: float

    end_time: float = 0.0
    max_steps: int = 10000

    equivalence_ratio: float = 0.0
    fuel: Dict = {}
    oxidizer: Dict = {}
    reactants: Dict = {}
    composition_type: str = "mole"


class InputPSR(NamedTuple):
    """Holds input parameters for a single perfectly stirred reactor (PSR) case.

    PSR cases are modeled as adiabatic and constant-pressure.
    ``temperature`` is the inlet temperature.
    """

    temperature: float
    pressure: float

    equivalence_ratio: float = 0.0
    fuel: Dict = {}
    oxidizer: Dict = {}
    reactants: Dict = {}
    composition_type: str = "mole"


class InputLaminarFlame(NamedTuple):
    """Holds input parameters for single laminar flame simulation.

    Freely-propagating laminar flames are inherently constant pressure.
    """

    temperature: float
    pressure: float

    width: float = 0.1

    equivalence_ratio: float = 0.0
    fuel: Dict = {}
    oxidizer: Dict = {}
    reactants: Dict = {}
    composition_type: str = "mole"


def ignition_sample_worker(sim_tuple):
    """Worker for multiprocessing of autoignition cases with data sampling.

    Runs, processes, and cleans up the case entirely within the worker, returning
    only the (picklable) metric and sampled data -- mirroring
    ``flame_sample_worker``. The intermediate HDF5 file written by ``run_case`` is
    read by ``process_results`` and removed by ``clean`` here, so the live
    (unpicklable) reactor never needs to cross the process boundary.

    Parameters
    ----------
    sim_tuple : tuple
        Tuple of IgnitionSimulation object to be run and identifier

    Returns
    -------
    dict
        Case identifier mapped to a tuple of the ignition delay and sampled data

    """
    sim, idx = sim_tuple

    sim.setup_case()
    sim.run_case()
    ignition_delay, data = sim.process_results()
    sim.clean()

    return {idx: (ignition_delay, data)}


def ignition_worker(sim_tuple):
    """Worker for multiprocessing of ignition delay only cases.

    Parameters
    ----------
    sim_tuple : tuple
        Tuple of IgnitionSimulation object to be run and identifier

    Returns
    -------
    dict
        Case identifier and calculated ignition delay

    """
    sim, idx = sim_tuple
    sim.setup_case()
    ignition_delay = sim.calculate()
    return {idx: ignition_delay}


def flame_sample_worker(flamesim_tuple):
    """Worker for multiprocessing of laminar flame cases with data sampling.

    Parameters
    ----------
    flamesim_tuple : tuple
        Tuple of FlameSimulation object to be run and identifier

    Returns
    -------
    dict
        Case identifier mapped to a tuple of the flame speed and sampled data

    """
    sim, idx = flamesim_tuple
    sim.setup_case()
    flame_speed, data = sim.process_results()
    return {idx: (flame_speed, data)}


def flame_worker(flamesim_tuple):
    """Worker for multiprocessing of laminar flame speed only cases.

    Parameters
    ----------
    flamesim_tuple : tuple
        Tuple of FlameSimulation object to be run and identifier

    Returns
    -------
    dict
        Case identifier and calculated laminar flame speed

    """
    sim, idx = flamesim_tuple
    sim.setup_case()
    flame_speed = sim.calculate()
    return {idx: flame_speed}


def psr_sample_worker(sim_tuple):
    """Worker for multiprocessing of PSR cases with data sampling.

    Runs, processes, and samples the case entirely within the worker (the PSR
    solver writes no intermediate file), returning the (picklable) metric vector
    and sampled data.

    Parameters
    ----------
    sim_tuple : tuple
        Tuple of PSRSimulation object to be run and identifier

    Returns
    -------
    dict
        Case identifier mapped to a tuple of the PSR metric vector and sampled data

    """
    sim, idx = sim_tuple
    sim.setup_case()
    metrics, data = sim.process_results()
    return {idx: (metrics, data)}


def psr_worker(sim_tuple):
    """Worker for multiprocessing of PSR metric-only cases.

    Parameters
    ----------
    sim_tuple : tuple
        Tuple of PSRSimulation object to be run and identifier

    Returns
    -------
    dict
        Case identifier mapped to the PSR metric vector

    """
    sim, idx = sim_tuple
    sim.setup_case()
    metrics = sim.calculate()
    return {idx: metrics}


def _run_workers(simulations, worker, num_threads):
    """Run ``worker`` over a list of job tuples and merge the per-case results.

    Parameters
    ----------
    simulations : list
        List of ``[sim, idx]`` job tuples to pass to ``worker``.
    worker : callable
        Worker returning a ``{idx: result}`` dict for one job.
    num_threads : int
        Number of processes to use; 1 runs serially.

    Returns
    -------
    dict
        Merged ``{idx: result}`` mapping over all cases.

    """
    jobs = tuple(simulations)
    if num_threads == 1:
        results = [worker(job) for job in jobs]
    else:
        pool = multiprocessing.Pool(processes=num_threads)
        results = pool.map(worker, jobs)
        pool.close()
        pool.join()
    return {key: val for k in results for key, val in k.items()}


def _run_sampling_jobs(simulations, worker, num_threads):
    """Run data-sampling workers and collect their metrics and sampled data.

    For ``ignition_sample_worker`` / ``flame_sample_worker`` / ``psr_sample_worker``,
    whose results are ``{idx: (metric, sampled_data)}``. The per-case metric may be
    a scalar (ignition delay, flame speed) or a vector (the PSR metric triple); it
    is flattened so each case contributes one or more entries, in case order.

    Returns
    -------
    metrics : numpy.ndarray
        1-D array of the per-case global metrics, concatenated in case order.
    data : numpy.ndarray
        Stacked sampled-state rows from all cases.

    """
    results = _run_workers(simulations, worker, num_threads)
    metrics = []
    data = []
    for idx in range(len(results)):
        metric, case_data = results[idx]
        metrics.append(np.atleast_1d(metric))
        data += list(case_data)
    metrics = np.concatenate(metrics) if metrics else np.array([])
    return metrics, np.array(data)


def _run_metric_jobs(simulations, worker, num_threads):
    """Run metric-only workers and collect their per-case metrics.

    For ``ignition_worker`` / ``flame_worker`` / ``psr_worker``, whose results are
    ``{idx: metric}``. The per-case metric may be a scalar (ignition delay, flame
    speed) or a vector (the PSR metric triple); it is flattened so each case
    contributes one or more entries, in case order.

    Returns
    -------
    numpy.ndarray
        1-D array of the per-case global metrics, concatenated in case order.

    """
    results = _run_workers(simulations, worker, num_threads)
    metrics = [np.atleast_1d(results[idx]) for idx in range(len(results))]
    return np.concatenate(metrics) if metrics else np.array([])


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
    # if any zero ignition delays or flame speeds, set error to 100
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
    numpy.ndarray
        Combined metrics for the model (ignition delays followed by flame
        speeds), used for evaluating error

    """
    metrics = []

    if ignition_conditions:
        if os.path.isfile(data_files["output_ignition"]):
            metrics.append(
                np.atleast_1d(
                    np.genfromtxt(data_files["output_ignition"], delimiter=",")
                )
            )
        else:
            raise SystemError("Error, no ignition output file present.")

    if psr_conditions:
        if os.path.isfile(data_files["output_psr"]):
            metrics.append(
                np.atleast_1d(np.genfromtxt(data_files["output_psr"], delimiter=","))
            )
        else:
            raise SystemError("Error, no PSR output file present.")

    if flame_conditions:
        if os.path.isfile(data_files["output_flame"]):
            metrics.append(
                np.atleast_1d(np.genfromtxt(data_files["output_flame"], delimiter=","))
            )
        else:
            raise SystemError("Error, no laminar flame output file present.")

    return np.concatenate(metrics) if metrics else np.array([])


def sample_metrics(
    model,
    ignition_conditions,
    psr_conditions=[],
    flame_conditions=[],
    phase_name="",
    num_threads=1,
    path="",
    reuse_saved=False,
    min_flame_speed=None,
):
    """Evaluates metrics used for determining error of reduced model

    Supports autoignition delay, laminar flame speed, and PSR temperature response curve metrics.

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
        Optional name for phase to load from YAML file (e.g., 'gas').
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
    numpy.ndarray
        Combined metrics for the model (ignition delays, PSR temperature response
        curve, flame speeds), used for evaluating error

    """
    # If number of threads specified as 0, use either max number of available
    # cores minus 1, or use 1 if multiple cores not available.
    if not num_threads:
        num_threads = multiprocessing.cpu_count() - 1 or 1

    ignition_delays = np.array([])
    if ignition_conditions:
        ignition_delays = np.zeros(len(ignition_conditions))

        exists_output = os.path.isfile(data_files["output_ignition"])
        if reuse_saved and exists_output:
            ignition_delays = np.atleast_1d(
                np.genfromtxt(data_files["output_ignition"], delimiter=",")
            )

        if reuse_saved and len(ignition_delays) == len(ignition_conditions):
            logging.info(
                "Reusing existing autoignition samples for the starting model."
            )
        else:
            simulations = []
            for idx, case in enumerate(ignition_conditions):
                simulations.append(
                    [
                        IgnitionSimulation(
                            idx, case, model, phase_name=phase_name, path=path
                        ),
                        idx,
                    ]
                )

            ignition_delays = _run_metric_jobs(
                simulations, ignition_worker, num_threads
            )

    # PSR cases contribute three metrics each (extinction residence time and the
    # two response temperatures), concatenated in case order.
    psr_metrics = np.array([])
    if psr_conditions:
        # each PSR case yields three metrics: extinction tau and two temperatures
        n_psr_metrics = 3 * len(psr_conditions)

        exists_output = os.path.isfile(data_files["output_psr"])
        if reuse_saved and exists_output:
            psr_metrics = np.atleast_1d(
                np.genfromtxt(data_files["output_psr"], delimiter=",")
            )

        if reuse_saved and len(psr_metrics) == n_psr_metrics:
            logging.info("Reusing existing PSR samples for the starting model.")
        else:
            simulations = []
            for idx, case in enumerate(psr_conditions):
                simulations.append(
                    [
                        PSRSimulation(
                            idx, case, model, phase_name=phase_name, path=path
                        ),
                        idx,
                    ]
                )

            psr_metrics = _run_metric_jobs(simulations, psr_worker, num_threads)

    flame_speeds = np.array([])
    if flame_conditions:
        flame_speeds = np.zeros(len(flame_conditions))

        exists_output = os.path.isfile(data_files["output_flame"])
        if reuse_saved and exists_output:
            flame_speeds = np.atleast_1d(
                np.genfromtxt(data_files["output_flame"], delimiter=",")
            )

        if reuse_saved and len(flame_speeds) == len(flame_conditions):
            logging.info(
                "Reusing existing laminar flame samples for the starting model."
            )
        else:
            simulations = []
            for idx, case in enumerate(flame_conditions):
                simulations.append(
                    [
                        FlameSimulation(
                            idx,
                            case,
                            model,
                            phase_name=phase_name,
                            path=path,
                            min_flame_speed=min_flame_speed,
                        ),
                        idx,
                    ]
                )

            flame_speeds = _run_metric_jobs(simulations, flame_worker, num_threads)

    metric_arrays = [
        np.atleast_1d(m) for m in (ignition_delays, psr_metrics, flame_speeds) if m.size
    ]
    return np.concatenate(metric_arrays) if metric_arrays else np.array([])


def sample(
    model,
    ignition_conditions,
    psr_conditions=[],
    flame_conditions=[],
    phase_name="",
    num_threads=1,
    path="",
    min_flame_speed=None,
):
    """Samples thermochemical data and generates metrics for various phenomena.

    Supports autoignition delay, PSR, and laminar flame speed metrics.

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
        Optional name for phase to load from YAML file (e.g., 'gas').
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
        num_threads = multiprocessing.cpu_count() - 1 or 1

    ignition_delays = np.array([])
    ignition_data = []
    if ignition_conditions:
        ignition_delays = np.zeros(len(ignition_conditions))
        ignition_data = []

        # check for presence of data and output files; if present, reuse.
        matches_number = False
        matches_shape = False
        exists_data = os.path.isfile(data_files["data_ignition"])
        exists_output = os.path.isfile(data_files["output_ignition"])
        if exists_data and exists_output:
            ignition_delays = np.genfromtxt(
                data_files["output_ignition"], delimiter=","
            )
            ignition_data = np.genfromtxt(data_files["data_ignition"], delimiter=",")
            # need to check that saved data at least matches the number of cases
            matches_number = ignition_delays.size == len(
                ignition_conditions
            ) and ignition_data.shape[0] / 20 == len(ignition_conditions)

            # also check that expected data is right shape (e.g., in case number of species
            # has changed if running a new model)
            gas = ct.Solution(model, phase_name)
            matches_shape = ignition_data.shape[1] == 2 + gas.n_species

        if matches_number and matches_shape:
            logging.info(
                "Reusing existing autoignition samples for the starting model."
            )
        else:
            logging.info("Running autoignition simulations for starting model.")
            simulations = []
            for idx, case in enumerate(ignition_conditions):
                simulations.append(
                    [
                        IgnitionSimulation(
                            idx, case, model, phase_name=phase_name, path=path
                        ),
                        idx,
                    ]
                )

            ignition_delays, ignition_data = _run_sampling_jobs(
                simulations, ignition_sample_worker, num_threads
            )

            np.savetxt(data_files["data_ignition"], ignition_data, delimiter=",")
            np.savetxt(data_files["output_ignition"], ignition_delays, delimiter=",")

    # PSR cases contribute three metrics each (extinction residence time and the
    # two response temperatures) and ``num_sample_points`` sampled-state rows.
    psr_metrics = np.array([])
    psr_data = []
    if psr_conditions:
        n_rows = PSRSimulation.num_sample_points
        n_psr_metrics = 3 * len(psr_conditions)

        # check for presence of data and output files; if present, reuse.
        matches_number = False
        matches_shape = False
        exists_data = os.path.isfile(data_files["data_psr"])
        exists_output = os.path.isfile(data_files["output_psr"])
        if exists_data and exists_output:
            psr_metrics = np.atleast_1d(
                np.genfromtxt(data_files["output_psr"], delimiter=",")
            )
            psr_data = np.genfromtxt(data_files["data_psr"], delimiter=",")
            matches_number = psr_metrics.size == n_psr_metrics and psr_data.shape[
                0
            ] / n_rows == len(psr_conditions)

            gas = ct.Solution(model, phase_name)
            matches_shape = psr_data.shape[1] == 2 + gas.n_species

        if matches_number and matches_shape:
            logging.info("Reusing existing PSR samples for the starting model.")
        else:
            logging.info("Running PSR simulations for starting model.")
            simulations = []
            for idx, case in enumerate(psr_conditions):
                simulations.append(
                    [
                        PSRSimulation(
                            idx, case, model, phase_name=phase_name, path=path
                        ),
                        idx,
                    ]
                )

            psr_metrics, psr_data = _run_sampling_jobs(
                simulations, psr_sample_worker, num_threads
            )

            np.savetxt(data_files["data_psr"], psr_data, delimiter=",")
            np.savetxt(data_files["output_psr"], psr_metrics, delimiter=",")

    flame_speeds = np.array([])
    flame_data = []
    if flame_conditions:
        flame_speeds = np.zeros(len(flame_conditions))
        flame_data = []

        # check for presence of data and output files; if present, reuse.
        matches_number = False
        matches_shape = False
        exists_data = os.path.isfile(data_files["data_flame"])
        exists_output = os.path.isfile(data_files["output_flame"])
        if exists_data and exists_output:
            flame_speeds = np.atleast_1d(
                np.genfromtxt(data_files["output_flame"], delimiter=",")
            )
            flame_data = np.genfromtxt(data_files["data_flame"], delimiter=",")
            # need to check that saved data at least matches the number of cases
            matches_number = flame_speeds.size == len(
                flame_conditions
            ) and flame_data.shape[0] / 20 == len(flame_conditions)

            # also check that expected data is right shape (e.g., in case number of species
            # has changed if running a new model)
            gas = ct.Solution(model, phase_name)
            matches_shape = flame_data.shape[1] == 2 + gas.n_species

        if matches_number and matches_shape:
            logging.info(
                "Reusing existing laminar flame samples for the starting model."
            )
        else:
            logging.info("Running laminar flame simulations for starting model.")
            simulations = []
            for idx, case in enumerate(flame_conditions):
                simulations.append(
                    [
                        FlameSimulation(
                            idx,
                            case,
                            model,
                            phase_name=phase_name,
                            path=path,
                            min_flame_speed=min_flame_speed,
                        ),
                        idx,
                    ]
                )

            flame_speeds, flame_data = _run_sampling_jobs(
                simulations, flame_sample_worker, num_threads
            )

            np.savetxt(data_files["data_flame"], flame_data, delimiter=",")
            np.savetxt(data_files["output_flame"], flame_speeds, delimiter=",")

    # combine metrics and sampled data from all phenomena (ignition, PSR, flame)
    metric_arrays = [
        np.atleast_1d(m) for m in (ignition_delays, psr_metrics, flame_speeds) if m.size
    ]
    data_arrays = [
        np.asarray(d) for d in (ignition_data, psr_data, flame_data) if len(d)
    ]
    sampled_metrics = np.concatenate(metric_arrays) if metric_arrays else np.array([])
    sampled_data = np.vstack(data_arrays) if data_arrays else np.array([])
    return sampled_metrics, sampled_data


def parse_ignition_inputs(model, conditions, phase_name=""):
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
        pre = f"Ignition input {idx}: "

        # check required keys
        kind = case.get("kind", "")
        temperature = case.get("temperature", 0.0)
        pressure = case.get("pressure", 0.0)

        assert kind in ["constant volume", "constant pressure"], (
            pre + '"case" needs to be "constant volume" or "constant pressure'
        )
        assert temperature > 0.0, pre + '"temperature" needs to be > 0'
        assert pressure > 0.0, pre + '"pressure" needs to be a number > 0'

        end_time = case.get("end-time", 0)
        max_steps = case.get("max-steps", 10000)

        equiv_ratio = case.get("equivalence-ratio", 0.0)
        fuel = case.get("fuel", [])
        oxidizer = case.get("oxidizer", [])

        reactants = case.get("reactants", [])

        assert (bool(equiv_ratio or fuel or oxidizer) + bool(reactants)) == 1, (
            pre + "should specify either equivalence-ratio/fuel/oxidizer or reactants."
        )

        if equiv_ratio or fuel or oxidizer:
            assert equiv_ratio > 0.0, pre + 'needs non-zero "equivalence-ratio"'

            assert fuel, pre + 'needs "fuel" with at least one entry'
            for entry in fuel:
                assert fuel[entry] > 0, pre + entry + " value needs to be a number > 0"
                assert entry in gas.species_names, (
                    pre + "fuel species not in model: " + entry
                )

            assert oxidizer, pre + 'needs "oxidizer" with at least one entry'
            for entry in oxidizer:
                assert oxidizer[entry] > 0, (
                    pre + entry + " value needs to be a number > 0"
                )
                assert entry in gas.species_names, (
                    pre + "oxidizer species not in model: " + entry
                )

        if reactants:
            for entry in reactants:
                assert reactants[entry] > 0, (
                    pre + entry + " value needs to be a number > 0"
                )
                assert entry in gas.species_names, (
                    pre + "reactant not in model: " + entry
                )

        composition_type = case.get("composition-type", "mole")
        assert composition_type in ["mole", "mass"], (
            pre + 'composition-type must be "mole" or "mass"'
        )
        assert not (composition_type == "mass" and equiv_ratio), (
            pre + "composition-type: must be mole when specifying equivalence ratio"
        )

        inputs.append(
            InputIgnition(
                kind,
                temperature,
                pressure,
                end_time,
                max_steps,
                equiv_ratio,
                fuel,
                oxidizer,
                reactants,
                composition_type,
            )
        )

    return inputs


def parse_psr_inputs(model, conditions, phase_name=""):
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
    gas = ct.Solution(model, phase_name)

    inputs = []
    for idx, case in enumerate(conditions):
        pre = f"PSR input {idx}: "

        # check required keys
        temperature = case.get("temperature", 0.0)
        pressure = case.get("pressure", 0.0)

        assert temperature > 0.0, pre + '"temperature" needs to be > 0'
        assert pressure > 0.0, pre + '"pressure" needs to be a number > 0'

        equiv_ratio = case.get("equivalence-ratio", 0.0)
        fuel = case.get("fuel", [])
        oxidizer = case.get("oxidizer", [])

        reactants = case.get("reactants", [])

        assert (bool(equiv_ratio or fuel or oxidizer) + bool(reactants)) == 1, (
            pre + "should specify either equivalence-ratio/fuel/oxidizer or reactants."
        )

        if equiv_ratio or fuel or oxidizer:
            assert equiv_ratio > 0.0, pre + 'needs non-zero "equivalence-ratio"'

            assert fuel, pre + 'needs "fuel" with at least one entry'
            for entry in fuel:
                assert fuel[entry] > 0, pre + entry + " value needs to be a number > 0"
                assert entry in gas.species_names, (
                    pre + "fuel species not in model: " + entry
                )

            assert oxidizer, pre + 'needs "oxidizer" with at least one entry'
            for entry in oxidizer:
                assert oxidizer[entry] > 0, (
                    pre + entry + " value needs to be a number > 0"
                )
                assert entry in gas.species_names, (
                    pre + "oxidizer species not in model: " + entry
                )

        if reactants:
            for entry in reactants:
                assert reactants[entry] > 0, (
                    pre + entry + " value needs to be a number > 0"
                )
                assert entry in gas.species_names, (
                    pre + "reactant not in model: " + entry
                )

        composition_type = case.get("composition-type", "mole")
        assert composition_type in ["mole", "mass"], (
            pre + 'composition-type must be "mole" or "mass"'
        )
        assert not (composition_type == "mass" and equiv_ratio), (
            pre + "composition-type: must be mole when specifying equivalence ratio"
        )

        inputs.append(
            InputPSR(
                temperature,
                pressure,
                equiv_ratio,
                fuel,
                oxidizer,
                reactants,
                composition_type,
            )
        )

    return inputs


def parse_flame_inputs(model, conditions, phase_name=""):
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
        pre = f"Flame input {idx}: "

        # check required keys
        temperature = case.get("temperature", 0.0)
        pressure = case.get("pressure", 0.0)

        assert temperature > 0.0, pre + '"temperature" needs to be > 0'
        assert pressure > 0.0, pre + '"pressure" needs to be a number > 0'

        width = case.get("width", 0.1)
        assert width > 0.0, pre + '"width" needs to be a number > 0'

        equiv_ratio = case.get("equivalence-ratio", 0.0)
        fuel = case.get("fuel", [])
        oxidizer = case.get("oxidizer", [])

        reactants = case.get("reactants", [])

        assert (bool(equiv_ratio or fuel or oxidizer) + bool(reactants)) == 1, (
            pre + "should specify either equivalence-ratio/fuel/oxidizer or reactants."
        )

        if equiv_ratio or fuel or oxidizer:
            assert equiv_ratio > 0.0, pre + 'needs non-zero "equivalence-ratio"'

            assert fuel, pre + 'needs "fuel" with at least one entry'
            for entry in fuel:
                assert fuel[entry] > 0, pre + entry + " value needs to be a number > 0"
                assert entry in gas.species_names, (
                    pre + "fuel species not in model: " + entry
                )

            assert oxidizer, pre + 'needs "oxidizer" with at least one entry'
            for entry in oxidizer:
                assert oxidizer[entry] > 0, (
                    pre + entry + " value needs to be a number > 0"
                )
                assert entry in gas.species_names, (
                    pre + "oxidizer species not in model: " + entry
                )

        if reactants:
            for entry in reactants:
                assert reactants[entry] > 0, (
                    pre + entry + " value needs to be a number > 0"
                )
                assert entry in gas.species_names, (
                    pre + "reactant not in model: " + entry
                )

        composition_type = case.get("composition-type", "mole")
        assert composition_type in ["mole", "mass"], (
            pre + 'composition-type must be "mole" or "mass"'
        )
        assert not (composition_type == "mass" and equiv_ratio), (
            pre + "composition-type: must be mole when specifying equivalence ratio"
        )

        inputs.append(
            InputLaminarFlame(
                temperature,
                pressure,
                width,
                equiv_ratio,
                fuel,
                oxidizer,
                reactants,
                composition_type,
            )
        )

    return inputs
