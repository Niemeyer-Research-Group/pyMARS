import sys
import os
import pkg_resources

import pytest
import numpy as np
import networkx as nx
import cantera as ct

from ..sampling import data_files, InputIgnition
from ..sensitivity_analysis import run_sa

# Taken from http://stackoverflow.com/a/22726782/1569494
try:
    from tempfile import TemporaryDirectory
except ImportError:
    from contextlib import contextmanager
    import shutil
    import tempfile
    import errno

    @contextmanager
    def TemporaryDirectory():
        name = tempfile.mkdtemp()
        try:
            yield name
        finally:
            try:
                shutil.rmtree(name)
            except OSError as e:
                # Reraise unless ENOENT: No such file or directory
                # (ok if directory has already been deleted)
                if e.errno != errno.ENOENT:
                    raise


def relative_location(file):
    file_path = os.path.join(file)
    return pkg_resources.resource_filename(__name__, file_path)


def check_equal(list1, list2):
    """Check whether two lists have the same contents (regardless of order).

    Taken from https://stackoverflow.com/a/12813909

    Parameters
    ----------
    list1 : list
        First list, containing all of a particular type
    list2: list
        Second list, containing all of a particular type

    Returns
    -------
    bool
        ``True`` if lists are equal

    """
    return len(list1) == len(list2) and sorted(list1) == sorted(list2)


class TestRunSA:
    def test_drgepsa(self):
        """Test SA using stored DRGEP result with upper_threshold = 0.5
        """
        starting_model = relative_location(os.path.join('assets', 'drgep_gri30.cti'))
        conditions = [
            InputIgnition(
                kind='constant volume', pressure=1.0, temperature=1000.0, equivalence_ratio=1.0,
                fuel={'CH4': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}
                ),
            InputIgnition(
                kind='constant volume', pressure=1.0, temperature=1200.0, equivalence_ratio=1.0,
                fuel={'CH4': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}
                ),
        ]
        data_files['output_ignition'] = relative_location(
            os.path.join('assets', 'example_ignition_output.txt')
            )
        data_files['data_ignition'] = relative_location(
            os.path.join('assets', 'example_ignition_data.dat')
            )
        
        limbo_species = ['H2', 'H2O2', 'CH2(S)', 'C2H4', 'C2H5', 'C2H6', 'HCCO', 'CH2CO']

        # Get expected model	
        expected_model = ct.Solution(relative_location(os.path.join('assets', 'drgepsa_gri30.cti')))

        # try using initial SA method
        with TemporaryDirectory() as temp_dir:
            reduced_model = run_sa(
                starting_model, 3.22, conditions, [], [], 5.0, ['N2'], 
                algorithm_type='initial', species_limbo=limbo_species[:], num_threads=1, 
                path=temp_dir
                )

        # Make sure models are the same	
        assert check_equal(reduced_model.model.species_names, expected_model.species_names)
        assert reduced_model.model.n_reactions == expected_model.n_reactions
        assert round(reduced_model.error, 2) == 3.20

        # try using greedy SA method
        with TemporaryDirectory() as temp_dir:
            reduced_model = run_sa(
                starting_model, 3.22, conditions, [], [], 5.0, ['N2'], 
                algorithm_type='greedy', species_limbo=limbo_species[:], num_threads=1, 
                path=temp_dir
                )
        
        # Make sure models are the same	
        assert check_equal(reduced_model.model.species_names, expected_model.species_names)
        assert reduced_model.model.n_reactions == expected_model.n_reactions
        assert round(reduced_model.error, 2) == 3.20
