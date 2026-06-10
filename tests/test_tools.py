"""Tests the tools module used by pyMARS"""

import os
import pathlib
import random
from tempfile import TemporaryDirectory

import cantera as ct

from pymars.tools import compare_models, convert


def relative_location(file):
    return str(pathlib.Path(__file__).parent / file)


# Activation energy conversion: cal/mol → J/kmol
_CAL_TO_J_KMOL = 4184.0


class TestCompareModels:
    def test_same(self):
        """Tests that two identical models are the same."""
        gas1 = ct.Solution("gri30.yaml")
        gas2 = ct.Solution("gri30.yaml")

        assert compare_models(gas1, gas2)

    def test_missing_reactions(self):
        """Checks that a random missing reaction is detected."""
        gas1 = ct.Solution("gri30.yaml")

        for __ in range(3):
            gas2 = ct.Solution("gri30.yaml")

            idx = random.randint(0, gas2.n_reactions - 1)
            reactions = gas2.reactions()
            del reactions[idx]

            gas2 = ct.Solution(
                species=gas2.species(),
                reactions=reactions,
                thermo="ideal-gas",
                kinetics="bulk",
            )

        assert not compare_models(gas1, gas2)

    def test_missing_species(self):
        """Checks that a random missing species is detected."""
        gas1 = ct.Solution("gri30.yaml")

        for __ in range(3):
            gas2 = ct.Solution("gri30.yaml")

            idx = random.randint(0, gas2.n_species - 1)
            species = gas2.species()
            sp = species.pop(idx)

            reactions = [
                rxn
                for rxn in gas2.reactions()
                if sp.name not in {**rxn.reactants, **rxn.products}
            ]

            gas2 = ct.Solution(
                species=species,
                reactions=reactions,
                thermo="ideal-gas",
                kinetics="bulk",
            )

        assert not compare_models(gas1, gas2)

    def test_species_different_thermo(self):
        """Checks detection of different species thermo parameter."""
        spA = ct.Species.from_yaml("""
            name: O2
            composition: {O: 2}
            thermo:
              model: NASA7
              temperature-ranges: [200.0, 1000.0, 3500.0]
              data:
              - [3.782456360e+00, -2.996734160e-03, 9.847302010e-06, -9.681295090e-09,
                 3.243728370e-12, -1.063943560e+03, 3.657675730e+00]
              - [3.282537840e+00, 1.483087540e-03, -7.579666690e-07, 2.094705550e-10,
                 -2.167177940e-14, -1.088457720e+03, 5.453231290e+00]
            """)
        spB = ct.Species.from_yaml("""
            name: O2
            composition: {O: 2}
            thermo:
              model: NASA7
              temperature-ranges: [200.0, 1000.0, 3500.0]
              data:
              - [0.000000000e+00, -2.996734160e-03, 9.847302010e-06, -9.681295090e-09,
                 3.243728370e-12, -1.063943560e+03, 3.657675730e+00]
              - [3.282537840e+00, 1.483087540e-03, -7.579666690e-07, 2.094705550e-10,
                 -2.167177940e-14, -1.088457720e+03, 5.453231290e+00]
            """)
        gas1 = ct.Solution(
            species=[spA], reactions=[], thermo="ideal-gas", kinetics="bulk"
        )
        gas2 = ct.Solution(
            species=[spB], reactions=[], thermo="ideal-gas", kinetics="bulk"
        )
        assert not compare_models(gas1, gas2)

    def test_species_different_composition(self):
        """Checks detection of different species composition."""
        spA = ct.Species.from_yaml("""
            name: O2
            composition: {O: 2}
            thermo:
              model: NASA7
              temperature-ranges: [200.0, 1000.0, 3500.0]
              data:
              - [3.782456360e+00, -2.996734160e-03, 9.847302010e-06, -9.681295090e-09,
                 3.243728370e-12, -1.063943560e+03, 3.657675730e+00]
              - [3.282537840e+00, 1.483087540e-03, -7.579666690e-07, 2.094705550e-10,
                 -2.167177940e-14, -1.088457720e+03, 5.453231290e+00]
            """)
        spB = ct.Species.from_yaml("""
            name: O2
            composition: {O: 1}
            thermo:
              model: NASA7
              temperature-ranges: [200.0, 1000.0, 3500.0]
              data:
              - [3.782456360e+00, -2.996734160e-03, 9.847302010e-06, -9.681295090e-09,
                 3.243728370e-12, -1.063943560e+03, 3.657675730e+00]
              - [3.282537840e+00, 1.483087540e-03, -7.579666690e-07, 2.094705550e-10,
                 -2.167177940e-14, -1.088457720e+03, 5.453231290e+00]
            """)
        gas1 = ct.Solution(
            species=[spA], reactions=[], thermo="ideal-gas", kinetics="bulk"
        )
        gas2 = ct.Solution(
            species=[spB], reactions=[], thermo="ideal-gas", kinetics="bulk"
        )
        assert not compare_models(gas1, gas2)

    def test_species_different_name(self):
        """Checks detection of different species name."""
        spA = ct.Species.from_yaml("""
            name: OO
            composition: {O: 2}
            thermo:
              model: NASA7
              temperature-ranges: [200.0, 1000.0, 3500.0]
              data:
              - [3.782456360e+00, -2.996734160e-03, 9.847302010e-06, -9.681295090e-09,
                 3.243728370e-12, -1.063943560e+03, 3.657675730e+00]
              - [3.282537840e+00, 1.483087540e-03, -7.579666690e-07, 2.094705550e-10,
                 -2.167177940e-14, -1.088457720e+03, 5.453231290e+00]
            """)
        spB = ct.Species.from_yaml("""
            name: O2
            composition: {O: 1}
            thermo:
              model: NASA7
              temperature-ranges: [200.0, 1000.0, 3500.0]
              data:
              - [3.782456360e+00, -2.996734160e-03, 9.847302010e-06, -9.681295090e-09,
                 3.243728370e-12, -1.063943560e+03, 3.657675730e+00]
              - [3.282537840e+00, 1.483087540e-03, -7.579666690e-07, 2.094705550e-10,
                 -2.167177940e-14, -1.088457720e+03, 5.453231290e+00]
            """)
        gas1 = ct.Solution(
            species=[spA], reactions=[], thermo="ideal-gas", kinetics="bulk"
        )
        gas2 = ct.Solution(
            species=[spB], reactions=[], thermo="ideal-gas", kinetics="bulk"
        )
        assert not compare_models(gas1, gas2)

    def test_reaction_different_arrhenius(self):
        gas = ct.Solution("gri30.yaml")

        R1 = ct.Reaction(
            equation="O + H2 <=> H + OH",
            rate=ct.ArrheniusRate(3.87e4, 2.7, 6260 * _CAL_TO_J_KMOL),
        )
        gas1 = ct.Solution(
            species=gas.species(), reactions=[R1], thermo="ideal-gas", kinetics="bulk"
        )

        R2 = ct.Reaction(
            equation="O + H2 <=> H + OH",
            rate=ct.ArrheniusRate(6.87e4, 2.7, 6260 * _CAL_TO_J_KMOL),
        )
        gas2 = ct.Solution(
            species=gas.species(), reactions=[R2], thermo="ideal-gas", kinetics="bulk"
        )

        assert not compare_models(gas1, gas2)

    def test_reaction_different_stoich(self):
        gas = ct.Solution("gri30.yaml")

        R1 = ct.Reaction(
            equation="O + H2 <=> H + OH",
            rate=ct.ArrheniusRate(3.87e4, 2.7, 6260 * _CAL_TO_J_KMOL),
        )
        gas1 = ct.Solution(
            species=gas.species(), reactions=[R1], thermo="ideal-gas", kinetics="bulk"
        )

        R2 = ct.Reaction(
            equation="O + H2 <=> H2O",
            rate=ct.ArrheniusRate(3.87e4, 2.7, 6260 * _CAL_TO_J_KMOL),
        )
        gas2 = ct.Solution(
            species=gas.species(), reactions=[R2], thermo="ideal-gas", kinetics="bulk"
        )

        assert not compare_models(gas1, gas2)

        R1 = ct.Reaction(
            equation="O + H2 <=> H + OH",
            rate=ct.ArrheniusRate(3.87e4, 2.7, 6260 * _CAL_TO_J_KMOL),
        )
        gas1 = ct.Solution(
            species=gas.species(), reactions=[R1], thermo="ideal-gas", kinetics="bulk"
        )

        R2 = ct.Reaction(
            equation="H2O <=> H + OH",
            rate=ct.ArrheniusRate(3.87e4, 2.7, 6260 * _CAL_TO_J_KMOL),
        )
        gas2 = ct.Solution(
            species=gas.species(), reactions=[R2], thermo="ideal-gas", kinetics="bulk"
        )

        assert not compare_models(gas1, gas2)


class TestConvert:
    def test_convert_gri_to_chemkin(self):
        """Test converting Cantera version of gri30 to Chemkin."""
        with TemporaryDirectory() as temp_dir:
            output = convert("gri30.yaml", path=temp_dir)
            assert output == [
                os.path.join(temp_dir, "gri30.inp"),
                os.path.join(temp_dir, "gri30_thermo.dat"),
                os.path.join(temp_dir, "gri30_transport.dat"),
            ]

    def test_convert_gri_to_cantera(self):
        """Test converting Chemkin version of gri30 to Cantera."""
        with TemporaryDirectory() as temp_dir:
            output = convert(
                relative_location(os.path.join("assets", "gri30.inp")), path=temp_dir
            )
            assert output == os.path.join(temp_dir, "gri30.yaml")
            assert compare_models(ct.Solution(output), ct.Solution("gri30.yaml"))
