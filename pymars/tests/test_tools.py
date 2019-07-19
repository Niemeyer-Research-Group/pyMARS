""" Tests the tools module used by pyMARS """

import os
import pkg_resources
import random
from tempfile import TemporaryDirectory

import pytest
import cantera as ct

from ..tools import compare_models, convert

def relative_location(file):
    file_path = os.path.join(file)
    return pkg_resources.resource_filename(__name__, file_path)


class TestCompareModels:
    def test_same(self):
        """Tests that two identical models are the same.
        """
        gas1 = ct.Solution('gri30.cti')
        gas2 = ct.Solution('gri30.cti')

        assert compare_models(gas1, gas2)
    
    def test_missing_reactions(self):
        """Checks that a random missing reaction is detected.
        """
        gas1 = ct.Solution('gri30.cti')

        for __ in range(3):
            gas2 = ct.Solution('gri30.cti')

            idx = random.randint(0, gas2.n_reactions - 1)
            reactions = gas2.reactions()
            del reactions[idx]

            gas2 = ct.Solution(
                species=gas2.species(), reactions=reactions,
                thermo='IdealGas', kinetics='GasKinetics'
                )

        assert not compare_models(gas1, gas2)

    def test_missing_species(self):
        """Checks that a random missing species is detected.
        """
        gas1 = ct.Solution('gri30.cti')

        for __ in range(3):
            gas2 = ct.Solution('gri30.cti')

            idx = random.randint(0, gas2.n_species - 1)
            species = gas2.species()
            sp = species.pop(idx)
            
            reactions = [rxn for rxn in gas2.reactions() 
                         if sp.name not in {**rxn.reactants, **rxn.products}
                         ]

            gas2 = ct.Solution(
                species=species, reactions=reactions,
                thermo='IdealGas', kinetics='GasKinetics'
                )

        assert not compare_models(gas1, gas2)
    
    def test_species_different_thermo(self):
        """Checks detection of different species thermo parameter.
        """
        spA = ct.Species.fromCti(
            '''species(name = "O2",
               atoms = "O:2",
               thermo = (
                NASA( [ 200.00, 1000.00], [ 3.782456360E+00, -2.996734160E-03,
                    9.847302010E-06, -9.681295090E-09, 3.243728370E-12,
                    -1.063943560E+03, 3.657675730E+00] ),
                NASA( [ 1000.00, 3500.00], [ 3.282537840E+00, 1.483087540E-03,
                    -7.579666690E-07, 2.094705550E-10, -2.167177940E-14,
                    -1.088457720E+03, 5.453231290E+00] ) ) )'''
                )
        spB = ct.Species.fromCti(
            '''species(name = "O2",
               atoms = "O:2",
               thermo = (
                NASA( [ 200.00, 1000.00], [ 0.000000000+00, -2.996734160E-03,
                    9.847302010E-06, -9.681295090E-09, 3.243728370E-12,
                    -1.063943560E+03, 3.657675730E+00] ),
                NASA( [ 1000.00, 3500.00], [ 3.282537840E+00, 1.483087540E-03,
                    -7.579666690E-07, 2.094705550E-10, -2.167177940E-14,
                    -1.088457720E+03, 5.453231290E+00] ) ) )'''
                )
        gas1 = ct.Solution(species=[spA], reactions=[], thermo='IdealGas', kinetics='GasKinetics')
        gas2 = ct.Solution(species=[spB], reactions=[], thermo='IdealGas', kinetics='GasKinetics')
        assert not compare_models(gas1, gas2)
    
    def test_species_different_composition(self):
        """Checks detection of different species composition.
        """
        spA = ct.Species.fromCti(
            '''species(name = "O2",
               atoms = "O:2",
               thermo = (
                NASA( [ 200.00, 1000.00], [ 3.782456360E+00, -2.996734160E-03,
                    9.847302010E-06, -9.681295090E-09, 3.243728370E-12,
                    -1.063943560E+03, 3.657675730E+00] ),
                NASA( [ 1000.00, 3500.00], [ 3.282537840E+00, 1.483087540E-03,
                    -7.579666690E-07, 2.094705550E-10, -2.167177940E-14,
                    -1.088457720E+03, 5.453231290E+00] ) ) )'''
            )
        spB = ct.Species.fromCti(
            '''species(name = "O2",
               atoms = "O:1",
               thermo = (
                NASA( [ 200.00, 1000.00], [ 3.782456360E+00, -2.996734160E-03,
                    9.847302010E-06, -9.681295090E-09, 3.243728370E-12,
                    -1.063943560E+03, 3.657675730E+00] ),
                NASA( [ 1000.00, 3500.00], [ 3.282537840E+00, 1.483087540E-03,
                    -7.579666690E-07, 2.094705550E-10, -2.167177940E-14,
                    -1.088457720E+03, 5.453231290E+00] ) ) )'''
            )
        gas1 = ct.Solution(species=[spA], reactions=[], thermo='IdealGas', kinetics='GasKinetics')
        gas2 = ct.Solution(species=[spB], reactions=[], thermo='IdealGas', kinetics='GasKinetics')
        assert not compare_models(gas1, gas2)
    
    def test_species_different_name(self):
        """Checks detection of different species name.
        """
        spA = ct.Species.fromCti(
            '''species(name = "OO",
               atoms = "O:2",
               thermo = (
                NASA( [ 200.00, 1000.00], [ 3.782456360E+00, -2.996734160E-03,
                    9.847302010E-06, -9.681295090E-09, 3.243728370E-12,
                    -1.063943560E+03, 3.657675730E+00] ),
                NASA( [ 1000.00, 3500.00], [ 3.282537840E+00, 1.483087540E-03,
                    -7.579666690E-07, 2.094705550E-10, -2.167177940E-14,
                    -1.088457720E+03, 5.453231290E+00] ) ) )'''
            )
        spB = ct.Species.fromCti(
            '''species(name = "O2",
               atoms = "O:1",
               thermo = (
                NASA( [ 200.00, 1000.00], [ 3.782456360E+00, -2.996734160E-03,
                    9.847302010E-06, -9.681295090E-09, 3.243728370E-12,
                    -1.063943560E+03, 3.657675730E+00] ),
                NASA( [ 1000.00, 3500.00], [ 3.282537840E+00, 1.483087540E-03,
                    -7.579666690E-07, 2.094705550E-10, -2.167177940E-14,
                    -1.088457720E+03, 5.453231290E+00] ) ) )'''
            )
        gas1 = ct.Solution(species=[spA], reactions=[], thermo='IdealGas', kinetics='GasKinetics')
        gas2 = ct.Solution(species=[spB], reactions=[], thermo='IdealGas', kinetics='GasKinetics')
        assert not compare_models(gas1, gas2)

    def test_reaction_different_arrhenius(self):
        gas = ct.Solution('gri30.cti')

        R1 = ct.ElementaryReaction.fromCti(
            '''reaction( "O + H2 <=> H + OH", [3.87000E+04, 2.7, 6260])'''
            )
        gas1 = ct.Solution(species=gas.species(), reactions=[R1], thermo='IdealGas', kinetics='GasKinetics')

        R2 = ct.ElementaryReaction.fromCti(
            '''reaction( "O + H2 <=> H + OH", [6.87000E+04, 2.7, 6260])'''
            )
        gas2 = ct.Solution(species=gas.species(), reactions=[R2], thermo='IdealGas', kinetics='GasKinetics')
        
        assert not compare_models(gas1, gas2)

    def test_reaction_different_stoich(self):
        gas = ct.Solution('gri30.cti')

        R1 = ct.ElementaryReaction.fromCti(
            '''reaction( "O + H2 <=> H + OH", [3.87000E+04, 2.7, 6260])'''
            )
        gas1 = ct.Solution(species=gas.species(), reactions=[R1], thermo='IdealGas', kinetics='GasKinetics')

        R2 = ct.ElementaryReaction.fromCti(
            '''reaction( "O + H2 <=> H2O", [3.87000E+04, 2.7, 6260])'''
            )
        gas2 = ct.Solution(species=gas.species(), reactions=[R2], thermo='IdealGas', kinetics='GasKinetics')
        
        assert not compare_models(gas1, gas2)

        R1 = ct.ElementaryReaction.fromCti(
            '''reaction( "O + H2 <=> H + OH", [3.87000E+04, 2.7, 6260])'''
            )
        gas1 = ct.Solution(species=gas.species(), reactions=[R1], thermo='IdealGas', kinetics='GasKinetics')

        R2 = ct.ElementaryReaction.fromCti(
            '''reaction( "H2O <=> H + OH", [3.87000E+04, 2.7, 6260])'''
            )
        gas2 = ct.Solution(species=gas.species(), reactions=[R2], thermo='IdealGas', kinetics='GasKinetics')
        
        assert not compare_models(gas1, gas2)
    

class TestConvert:
    def test_convert_gri_to_chemkin(self):
        """Test converting Cantera version of gri30 to Chemkin.
        """
        with TemporaryDirectory() as temp_dir:
            output = convert('gri30.cti', path=temp_dir)
            assert output == [
                os.path.join(temp_dir, 'gri30.inp'), 
                os.path.join(temp_dir, 'gri30_thermo.dat'), 
                os.path.join(temp_dir, 'gri30_transport.dat')
                ]
    
    def test_convert_gri_to_cantera(self):
        """Test converting Chemkin version of gri30 to Cantera.
        """
        with TemporaryDirectory() as temp_dir:
            output = convert(relative_location(os.path.join('assets', 'gri30.inp')), path=temp_dir)
            assert output == os.path.join(temp_dir, 'gri30.cti')
            assert compare_models(ct.Solution(output), ct.Solution('gri30.cti'))
