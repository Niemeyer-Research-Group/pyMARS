""" Tests the create output file unit used by pyMARS """

import os
import pkg_resources

import pytest
import numpy as np
import cantera as ct

from ..tools import compare_models
from ..reduce_model import trim
from ..soln2cti import write, build_efficiencies

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

class TestBuildEfficiencies:
    def test_regular_efficiencies(self):
        """Test returning regular efficiency string.
        """
        efficiencies = {
            'co': 2.8, 'co2': 1.6, 'h2': 3.7, 'h2o': 0.0, 'h2o2': 7.7, 
            'he': 0.65, 'n2': 1.5, 'o2': 1.2
            }
        species = ['co', 'co2', 'h2', 'h2o', 'h2o2', 'he', 'n2', 'o2']
        efficiency_str = build_efficiencies(efficiencies, species)
        assert efficiency_str == 'co:2.8  co2:1.6  h2:3.7  h2o:0.0  h2o2:7.7  he:0.65  n2:1.5  o2:1.2'
    
    def test_regular_efficiencies(self):
        """Test returning regular efficiency string with removed species.
        """
        efficiencies = {
            'co': 2.8, 'co2': 1.6, 'h2': 3.7, 'h2o': 0.0, 'h2o2': 7.7, 
            'he': 0.65, 'n2': 1.5, 'o2': 1.2
            }
        species = ['co', 'co2', 'h2', 'h2o', 'h2o2', 'n2', 'o2']
        efficiency_str = build_efficiencies(efficiencies, species)
        assert efficiency_str == 'co:2.8  co2:1.6  h2:3.7  h2o:0.0  h2o2:7.7  n2:1.5  o2:1.2'
    
    def test_explicit_third_body(self):
        """Test appropriate handling of reaction with explicit third body.
        """
        rxn = ct.FalloffReaction.fromCti(
        "units(length='cm', quantity='mol')\n" +
        '''falloff_reaction('h2o2 (+ h2o) <=> oh + oh (+ h2o)',
            kf=[2.000000e+12, 0.9, 48749.0],
            kf0=[1.865000e+25, -2.3, 48749.0],
            falloff=Troe(A=0.51, T3=1e-30, T1=1e+30))'''
        )
        efficiency_str = build_efficiencies(
            rxn.efficiencies, ['h2o2', 'h2o', 'oh'], rxn.default_efficiency
            )
        assert not efficiency_str


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

    def test_multiple_pressure_dependent_reactions(self):
        """Test appropriate writing of multiple pressure dependent reactions.

        We can encounter multiple instances of the same reaction, but with
        different explicit third-body species.
        """
        species_h = ct.Species.fromCti(
        '''species(name='h', atoms='H:1',
        thermo=(NASA([200.00, 1000.00],
                     [ 2.50000000E+00,  0.00000000E+00,  0.00000000E+00,
                       0.00000000E+00,  0.00000000E+00,  2.54736600E+04,
                      -4.46682850E-01]),
                NASA([1000.00, 6000.00],
                     [ 2.50000000E+00,  0.00000000E+00,  0.00000000E+00,
                       0.00000000E+00,  0.00000000E+00,  2.54736600E+04,
                      -4.46682850E-01])))'''
        )
        species_o2 = ct.Species.fromCti(
        '''species(name='o2', atoms='O:2',
        thermo=(NASA([200.00, 1000.00],
                     [ 3.78245636E+00, -2.99673416E-03,  9.84730201E-06,
                      -9.68129509E-09,  3.24372837E-12, -1.06394356E+03,
                       3.65767573E+00]),
                NASA([1000.00, 6000.00],
                     [ 3.66096065E+00,  6.56365811E-04, -1.41149627E-07,
                       2.05797935E-11, -1.29913436E-15, -1.21597718E+03,
                       3.41536279E+00])))'''
        )
        species_ho2 = ct.Species.fromCti(
            '''species(name='ho2', atoms='H:1 O:2',
        thermo=(NASA([200.00, 1000.00],
                     [ 4.30179807E+00, -4.74912097E-03,  2.11582905E-05,
                      -2.42763914E-08,  9.29225225E-12,  2.64018485E+02,
                       3.71666220E+00]),
                NASA([1000.00, 5000.00],
                     [ 4.17228741E+00,  1.88117627E-03, -3.46277286E-07,
                       1.94657549E-11,  1.76256905E-16,  3.10206839E+01,
                       2.95767672E+00])))'''
        )
        species_ar = ct.Species.fromCti(
            '''species(name='ar', atoms='Ar:1',
        thermo=(NASA([200.00, 1000.00],
                     [ 2.50000000E+00,  0.00000000E+00,  0.00000000E+00,
                       0.00000000E+00,  0.00000000E+00, -7.45375000E+02,
                       4.37967491E+00]),
                NASA([1000.00, 6000.00],
                     [ 2.50000000E+00,  0.00000000E+00,  0.00000000E+00,
                       0.00000000E+00,  0.00000000E+00, -7.45375000E+02,
                       4.37967491E+00])))'''
        )
        species_he = ct.Species.fromCti(
            '''species(name='he', atoms='He:1',
        thermo=(NASA([200.00, 1000.00],
                     [ 2.50000000E+00,  0.00000000E+00,  0.00000000E+00,
                       0.00000000E+00,  0.00000000E+00, -7.45375000E+02,
                       9.28723974E-01]),
                NASA([1000.00, 6000.00],
                     [ 2.50000000E+00,  0.00000000E+00,  0.00000000E+00,
                       0.00000000E+00,  0.00000000E+00, -7.45375000E+02,
                       9.28723974E-01])))'''
        )
        species_h2o = ct.Species.fromCti(
            '''species(name='h2o', atoms='H:2 O:1',
        thermo=(NASA([200.00, 1000.00],
                     [ 4.19863520E+00, -2.03640170E-03,  6.52034160E-06,
                      -5.48792690E-09,  1.77196800E-12, -3.02937260E+04,
                      -8.49009010E-01]),
                NASA([1000.00, 6000.00],
                     [ 2.67703890E+00,  2.97318160E-03, -7.73768890E-07,
                       9.44335140E-11, -4.26899910E-15, -2.98858940E+04,
                       6.88255000E+00])))'''
        )
        full_species = [species_h, species_o2, species_ho2, species_ar, species_he, species_h2o]
        
        R1 = ct.FalloffReaction.fromCti(
        "units(length='cm', quantity='mol')\n" +
        '''falloff_reaction('h + o2 (+ M) <=> ho2 (+ M)',
            kf=[4.650000e+12, 0.44, 0.0],
            kf0=[1.737000e+19, -1.23, 0.0],
            efficiencies='he:0.0 h2o:10.0 ar:0.0',
            falloff=Troe(A=0.67, T3=1e-30, T1=1e+30, T2=1e+30))'''
        )
        R2 = ct.FalloffReaction.fromCti(
        "units(length='cm', quantity='mol')\n" +
        '''falloff_reaction('h + o2 (+ ar) <=> ho2 (+ ar)',
            kf=[4.650000e+12, 0.44, 0.0],
            kf0=[6.810000e+18, -1.2, 0.0],
            falloff=Troe(A=0.7, T3=1e-30, T1=1e+30, T2=1e+30))'''
        )
        R3 = ct.FalloffReaction.fromCti(
        "units(length='cm', quantity='mol')\n" +
        '''falloff_reaction('h + o2 (+ he) <=> ho2 (+ he)',
            kf=[4.650000e+12, 0.44, 0.0],
            kf0=[9.192000e+18, -1.2, 0.0],
            falloff=Troe(A=0.59, T3=1e-30, T1=1e+30, T2=1e+30))'''
        )

        solution = ct.Solution(
            species=full_species, 
            reactions=[R1, R2, R3],
            thermo='IdealGas', kinetics='GasKinetics'
            )
        
        with TemporaryDirectory() as temp_dir:
            output = write(solution, 'test.cti', path=temp_dir)
            new_solution = ct.Solution(output)

        assert compare_models(solution, new_solution)
