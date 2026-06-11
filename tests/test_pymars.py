"""Tests for the top-level pymars input parsing."""

from pymars.pymars import parse_inputs


def _base_inputs(conditions, **overrides):
    """Build a minimal valid reduction-input dict for the given conditions."""
    inputs = {
        "model": "gri30.yaml",
        "error": 5.0,
        "method": "DRGEP",
        "targets": ["CH4"],
        "autoignition-conditions": conditions,
    }
    inputs.update(overrides)
    return inputs


class TestAutoRetainInputSpecies:
    """Fuel/oxidizer/reactant species should be retained automatically."""

    def test_fuel_and_oxidizer_retained(self):
        input_dict = _base_inputs(
            [
                {
                    "kind": "constant volume",
                    "pressure": 1.0,
                    "temperature": 1000.0,
                    "equivalence-ratio": 1.0,
                    "fuel": {"CH4": 1.0},
                    "oxidizer": {"O2": 1.0, "N2": 3.76},
                }
            ]
        )
        inputs = parse_inputs(input_dict)

        # O2 and N2 are inputs but neither target nor explicitly retained
        assert "O2" in inputs.safe_species
        assert "N2" in inputs.safe_species
        # CH4 is a target, so it should not be added to the retained list
        assert "CH4" not in inputs.safe_species

    def test_reactants_retained(self):
        input_dict = _base_inputs(
            [
                {
                    "kind": "constant volume",
                    "pressure": 1.0,
                    "temperature": 1000.0,
                    "reactants": {"CH4": 1.0, "O2": 2.0, "N2": 7.52},
                }
            ]
        )
        inputs = parse_inputs(input_dict)

        assert "O2" in inputs.safe_species
        assert "N2" in inputs.safe_species
        assert "CH4" not in inputs.safe_species

    def test_no_duplicates_for_already_retained(self):
        input_dict = _base_inputs(
            [
                {
                    "kind": "constant volume",
                    "pressure": 1.0,
                    "temperature": 1000.0,
                    "equivalence-ratio": 1.0,
                    "fuel": {"CH4": 1.0},
                    "oxidizer": {"O2": 1.0, "N2": 3.76},
                }
            ],
            **{"retained-species": ["N2"]},
        )
        inputs = parse_inputs(input_dict)

        assert inputs.safe_species.count("N2") == 1
        assert "O2" in inputs.safe_species

    def test_species_from_all_conditions_retained(self):
        input_dict = _base_inputs(
            [
                {
                    "kind": "constant volume",
                    "pressure": 1.0,
                    "temperature": 1000.0,
                    "equivalence-ratio": 1.0,
                    "fuel": {"CH4": 1.0},
                    "oxidizer": {"O2": 1.0, "N2": 3.76},
                },
                {
                    "kind": "constant volume",
                    "pressure": 1.0,
                    "temperature": 1200.0,
                    "reactants": {"H2": 1.0, "O2": 1.0, "AR": 3.76},
                },
            ]
        )
        inputs = parse_inputs(input_dict)

        for sp in ("O2", "N2", "H2", "AR"):
            assert sp in inputs.safe_species
