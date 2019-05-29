""" Tests the create output file unit used by pyMARS """

import os
import pkg_resources
from operator import attrgetter

import pytest
import numpy as np
import cantera as ct

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


def compare_models(model1, model2):
    """Checks whether two Cantera models are equivalent.

    Parameters
    ----------
    model1 : ct.Solution
        First Cantera model object
    model2 : ct.Solution
        Second Cantera model object

    Returns
    -------
    bool
        ``True`` if models match, ``False`` if models disagree
    
    """
    # at minimum, numbers of species and reactions should be the same
    if model1.n_species != model2.n_species:
        return False
    if model1.n_reactions != model2.n_reactions:
        return False

    for sp1, sp2 in zip(
        sorted(model1.species(), key=attrgetter('name')), 
        sorted(model2.species(), key=attrgetter('name'))
        ):
        if sp1.name != sp2.name:
            return False

        if sp1.composition != sp2.composition:
            return False
        
        if sp1.thermo.n_coeffs == sp2.thermo.n_coeffs:
            if any(sp1.thermo.coeffs != sp2.thermo.coeffs):
                return False
        else:
            return False
        
        if hasattr(sp1, 'transport') or hasattr(sp2, 'transport'):
            if hasattr(sp1, 'transport') and hasattr(sp2, 'transport'):
                # iterate over transport parameters
                params = [a for a in dir(sp1.transport) 
                          if not a.startswith('__') and 
                          not callable(getattr(sp1.transport, a))
                          ]
                for attr in params:
                    if getattr(sp1.transport, attr) != getattr(sp2.transport, attr, 0.0):
                        return False
            else:
                return False

    for rxn1, rxn2 in zip(
        sorted(model1.reactions(), key=attrgetter('equation')), 
        sorted(model2.reactions(), key=attrgetter('equation'))
        ):
        if type(rxn1) != type(rxn2):
            return False
        
        if rxn1.reactants != rxn2.reactants:
            return False
        if rxn1.products != rxn2.products:
            return False

        if rxn1.duplicate != rxn2.duplicate:
            return False

        # Check rate parameters for elementary and third-body reactions
        if hasattr(rxn1, 'rate') or hasattr(rxn2, 'rate'):
            if hasattr(rxn1, 'rate') and hasattr(rxn2, 'rate'):
                if type(rxn1.rate) != type(rxn2.rate):
                    return False
                if len(dir(rxn1.rate)) != len(dir(rxn2.rate)):
                    return False
                params = [
                    a for a in dir(rxn1.rate) 
                    if not a.startswith('__') and not callable(getattr(rxn1.rate, a))
                    ]
                for attr in params:
                    if getattr(rxn1.rate, attr) != getattr(rxn2.rate, attr, 0.0):
                        return False
            else:
                return False
        
        # For falloff and chemically activated reactions, check low and high rates
        if hasattr(rxn1, 'low_rate') or hasattr(rxn2, 'low_rate'):
            if hasattr(rxn1, 'low_rate') and hasattr(rxn2, 'low_rate'):
                if type(rxn1.low_rate) != type(rxn2.low_rate):
                    return False
                if len(dir(rxn1.low_rate)) != len(dir(rxn2.low_rate)):
                    return False
                params = [
                    a for a in dir(rxn1.low_rate) 
                    if not a.startswith('__') and not callable(getattr(rxn1.low_rate, a))
                    ]
                for attr in params:
                    if getattr(rxn1.low_rate, attr) != getattr(rxn2.low_rate, attr, 0.0):
                        return False
            else:
                return False

        if hasattr(rxn1, 'high_rate') or hasattr(rxn2, 'high_rate'):
            if hasattr(rxn1, 'high_rate') and hasattr(rxn2, 'high_rate'):
                if type(rxn1.high_rate) != type(rxn2.high_rate):
                    return False
                if len(dir(rxn1.high_rate)) != len(dir(rxn2.high_rate)):
                    return False
                params = [
                    a for a in dir(rxn1.high_rate) 
                    if not a.startswith('__') and not callable(getattr(rxn1.high_rate, a))
                    ]
                for attr in params:
                    if getattr(rxn1.high_rate, attr) != getattr(rxn2.high_rate, attr, 0.0):
                        return False
            else:
                return False
        
        # check Plog rates
        if hasattr(rxn1, 'rates') or hasattr(rxn2, 'rates'):
            if hasattr(rxn1, 'rates') and hasattr(rxn2, 'rates'):
                if len(rxn1.rates) != len(rxn2.rates):
                    return False
                for rate1, rate2 in zip(
                    sorted(rxn1.rates, key=lambda rate: rate[0]),
                    sorted(rxn2.rates, key=lambda rate: rate[0]),
                    ):
                    if not np.allclose(rate1[0], rate2[0]):
                        return False
                    params = ['activation_energy', 'pre_exponential_factor', 'temperature_exponent']
                    for param in params:
                        if getattr(rate1[1], param, 0.0) != getattr(rate2[1], param, 0.0):
                            return False
            
        # check Chebyshev parameters
        if hasattr(rxn1, 'coeffs') or hasattr(rxn2, 'coeffs'):
            if hasattr(rxn1, 'coeffs') and hasattr(rxn2, 'coeffs'):
                if rxn1.nPressure != rxn2.nPressure:
                    return False
                if rxn1.nTemperature != rxn2.nTemperature:
                    return False
                if (rxn1.Pmax != rxn2.Pmax) or (rxn1.Pmin != rxn1.Pmin):
                    return False
                if (rxn1.Tmax != rxn2.Tmax) or (rxn1.Tmin != rxn1.Tmin):
                    return False
                if not np.allclose(rxn1.coeffs, rxn2.coeffs):
                    return False
        
        if hasattr(rxn1, 'efficiencies') or hasattr(rxn2, 'efficiencies'):
            if hasattr(rxn1, 'efficiencies') and hasattr(rxn2, 'efficiencies'):
                if rxn1.efficiencies != rxn2.efficiencies:
                    return False
            else:
                return False
        
        # Check falloff parameters if any
        if hasattr(rxn1, 'falloff') or hasattr(rxn2, 'falloff') :
            if hasattr(rxn1, 'falloff') and hasattr(rxn2, 'falloff'):
                if len(rxn1.falloff.parameters) == len(rxn2.falloff.parameters):
                    if any(rxn1.falloff.parameters != rxn2.falloff.parameters):
                        return False
                else:
                    return False
            else:
                return False

    return True


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
        solution = ct.Solution(relative_location('artificial-mechanism.cti'))

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

        expected_solution = ct.Solution(relative_location('eout_gri30.cti'))

        assert compare_models(new_solution, expected_solution)

    def test_artificial_write_red(self):
        """Test writing slightly reduced version of artificial model.
        """
        solution = trim(relative_location('artificial-mechanism.cti'), ['H', 'O2'], 'gas.cti')

        with TemporaryDirectory() as temp_dir:        
            output = write(solution, 'gas.cti', path=temp_dir)
            new_solution = ct.Solution(output)

        expected_solution = ct.Solution(relative_location('eout_artificial-mechanism.cti'))

        assert compare_models(new_solution, expected_solution)

    def test_gri_chem_plog_cheb(self):
        """Test writing GRI Mech 3.0 with chemically activated, Plog, and Chebyshev reactions.
        """
        old_solution = ct.Solution('gri30.cti')

        R1 = ct.ChemicallyActivatedReaction.fromCti('''
        units(length='cm', quantity='mol')
        chemically_activated_reaction('CH3 + OH (+ M) <=> CH2O + H2 (+ M)',
            kLow=[2.823201e+02, 1.46878, (-3270.56495, 'cal/mol')],
            kHigh=[5.880000e-14, 6.721, (-3022.227, 'cal/mol')],
            falloff=Troe(A=1.671, T3=434.782, T1=2934.21, T2=3919.0))'''
        )

        R2 = ct.ChemicallyActivatedReaction.fromCti('''
        units(length='cm', quantity='mol')
        chemically_activated_reaction('CH3 (+ M) <=> CH2 + H (+ M)',
            kLow=[2.823201e+02, 1.46878, (-3270.56495, 'cal/mol')],
            kHigh=[5.880000e-14, 6.721, (-3022.227, 'cal/mol')],
            falloff=Troe(A=1.671, T3=434.782, T1=2934.21, T2=3919.0))'''
        )
        
        R3 = ct.PlogReaction.fromCti('''
        units(length='cm', quantity='mol')
        pdep_arrhenius('CH3 + 2 OH <=> CH2O + H2 + OH',
            [(0.001315789, 'atm'), 2.440000e+10, 1.04, 3980.0],
            [(0.039473684, 'atm'), 3.890000e+10, 0.989, 4114.0],
            [(1.0, 'atm'), 3.460000e+12, 0.442, 5463.0],
            [(10.0, 'atm'), 1.720000e+14, -0.01, 7134.0],
            [(100.0, 'atm'), 1.900000e+15, -0.29, 8306.0])'''
        )

        R4 = ct.ChebyshevReaction.fromCti('''
        units(length='cm', quantity='mol')
        chebyshev_reaction('CH3 + CH3 <=> 2 CH2 + H2',
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
