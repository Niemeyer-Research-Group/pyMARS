""" Tests the create output file unit used by pyMARS """

import os
import pkg_resources

import pytest
import numpy as np
import cantera as ct

from ..tools import compare_models
from ..reduce_model import trim
from ..soln2cti import write

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


class TestWrite:
    def test_GRI_write(self):
        """Test writing unmodified GRI Mech 3.0 returns same model.
        """
        solution = ct.Solution('gri30.cti')
        with TemporaryDirectory() as temp_dir:
            output = write(solution, 'pym_gri30.cti', path=temp_dir)
            new_solution = ct.Solution(output)

        assert compare_models(solution, new_solution)

    def test_artificial_write(self):
        """Test writing unmodified artificial model.
        """
        solution = ct.Solution(
            relative_location(os.path.join('assets', 'artificial-mechanism.cti'))
            )

        with TemporaryDirectory() as temp_dir:
            output = write(solution, 'pym_gas.cti', path=temp_dir)
            new_solution = ct.Solution(output)

        assert compare_models(solution, new_solution)
    
    def test_gri_write_red(self):
        """Test writing slightly reduced GRI Mech 3.0
        """
        solution = trim('gri30.cti', ['CH4', 'O2', 'N2'], 'reduced_gri30.cti')
        
        with TemporaryDirectory() as temp_dir:
            output = write(solution, 'reduced_gri30.cti', path=temp_dir)
            new_solution = ct.Solution(output)

        expected_solution = ct.Solution(
            relative_location(os.path.join('assets', 'eout_gri30.cti'))
            )

        assert compare_models(new_solution, expected_solution)

    def test_artificial_write_red(self):
        """Test writing slightly reduced version of artificial model.
        """
        solution = trim(
            relative_location(os.path.join('assets', 'artificial-mechanism.cti')), 
            ['H', 'O2'], 'gas.cti'
            )

        with TemporaryDirectory() as temp_dir:        
            output = write(solution, 'gas.cti', path=temp_dir)
            new_solution = ct.Solution(output)

        expected_solution = ct.Solution(
            relative_location(os.path.join('assets', 'eout_artificial-mechanism.cti'))
            )

        assert compare_models(new_solution, expected_solution)

    def test_gri_chem_plog_cheb(self):
        """Test writing GRI Mech 3.0 with chemically activated, Plog, and Chebyshev reactions.
        """
        old_solution = ct.Solution('gri30.cti')

        R1 = ct.ChemicallyActivatedReaction.fromCti(
        "units(length='cm', quantity='mol')\n" + 
        '''chemically_activated_reaction('CH3 + OH (+ M) <=> CH2O + H2 (+ M)',
            kLow=[2.823201e+02, 1.46878, (-3270.56495, 'cal/mol')],
            kHigh=[5.880000e-14, 6.721, (-3022.227, 'cal/mol')],
            falloff=Troe(A=1.671, T3=434.782, T1=2934.21, T2=3919.0))'''
        )

        R2 = ct.ChemicallyActivatedReaction.fromCti(
        "units(length='cm', quantity='mol')\n" +
        '''chemically_activated_reaction('CH3 (+ M) <=> CH2 + H (+ M)',
            kLow=[2.823201e+02, 1.46878, (-3270.56495, 'cal/mol')],
            kHigh=[5.880000e-14, 6.721, (-3022.227, 'cal/mol')],
            falloff=Troe(A=1.671, T3=434.782, T1=2934.21, T2=3919.0))'''
        )
        
        R3 = ct.PlogReaction.fromCti(
        "units(length='cm', quantity='mol')\n" +
        '''pdep_arrhenius('CH3 + 2 OH <=> CH2O + H2 + OH',
            [(0.001315789, 'atm'), 2.440000e+10, 1.04, (3980.0, 'cal/mol')],
            [(0.039473684, 'atm'), 3.890000e+10, 0.989, (4114.0, 'cal/mol')],
            [(1.0, 'atm'), 3.460000e+12, 0.442, (5463.0, 'cal/mol')],
            [(10.0, 'atm'), 1.720000e+14, -0.01, (7134.0, 'cal/mol')],
            [(100.0, 'atm'), 1.900000e+15, -0.29, (8306.0, 'cal/mol')])'''
        )

        R4 = ct.ChebyshevReaction.fromCti(
        "units(length='cm', quantity='mol')\n" +
        '''chebyshev_reaction('CH3 + CH3 <=> 2 CH2 + H2',
            Tmin=290.0, Tmax=3000.0,
            Pmin=(0.001, 'atm'), Pmax=(100.0, 'atm'),
            coeffs=[[-1.44280e+01,  2.59970e-01, -2.24320e-02, -2.78700e-03],
                    [ 2.20630e+01,  4.88090e-01, -3.96430e-02, -5.48110e-03],
                    [-2.32940e-01,  4.01900e-01, -2.60730e-02, -5.04860e-03],
                    [-2.93660e-01,  2.85680e-01, -9.33730e-03, -4.01020e-03],
                    [-2.26210e-01,  1.69190e-01,  4.85810e-03, -2.38030e-03],
                    [-1.43220e-01,  7.71110e-02,  1.27080e-02, -6.41540e-04]])'''
            )
        solution = ct.Solution(
            species=old_solution.species(), 
            reactions=(old_solution.reactions() + [R1, R2, R3, R4]),
            thermo='IdealGas', kinetics='GasKinetics'
            )

        with TemporaryDirectory() as temp_dir:
            output = write(solution, 'mod_gri30.cti', path=temp_dir)
            new_solution = ct.Solution(output)

        assert compare_models(solution, new_solution)
