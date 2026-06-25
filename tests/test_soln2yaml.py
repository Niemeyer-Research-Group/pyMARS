"""Tests for soln2yaml module"""

import os
import pathlib
from tempfile import TemporaryDirectory

import cantera as ct

from pymars.tools import compare_models
from pymars.reduce_model import trim
from pymars.soln2yaml import write


def relative_location(file):
    return str(pathlib.Path(__file__).parent / file)


class TestWrite:
    def test_GRI_write(self):
        """Test writing unmodified GRI Mech 3.0 returns same model."""
        solution = ct.Solution("gri30.yaml")
        with TemporaryDirectory() as temp_dir:
            output = write(solution, "pym_gri30.yaml", path=temp_dir)
            new_solution = ct.Solution(output)

        assert compare_models(solution, new_solution)

    def test_artificial_write(self):
        """Test writing unmodified artificial model."""
        solution = ct.Solution(
            relative_location(os.path.join("assets", "artificial-mechanism.yaml"))
        )

        with TemporaryDirectory() as temp_dir:
            output = write(solution, "pym_gas.yaml", path=temp_dir)
            new_solution = ct.Solution(output)

        assert compare_models(solution, new_solution)

    def test_gri_write_red(self):
        """Test writing slightly reduced GRI Mech 3.0."""
        solution = trim("gri30.yaml", ["CH4", "O2", "N2"], "reduced_gri30.yaml")

        with TemporaryDirectory() as temp_dir:
            output = write(solution, "reduced_gri30.yaml", path=temp_dir)
            new_solution = ct.Solution(output)

        assert compare_models(solution, new_solution)

    def test_artificial_write_red(self):
        """Test writing slightly reduced version of artificial model."""
        solution = trim(
            relative_location(os.path.join("assets", "artificial-mechanism.yaml")),
            ["H"],
            "gas.yaml",
        )

        with TemporaryDirectory() as temp_dir:
            output = write(solution, "gas.yaml", path=temp_dir)
            new_solution = ct.Solution(output)

        assert compare_models(solution, new_solution)

    def test_third_bodies_write(self):
        """Test writing model with third-body falloff reactions."""
        solution = ct.Solution(
            relative_location(os.path.join("assets", "model-third-bodies.yaml"))
        )

        with TemporaryDirectory() as temp_dir:
            output = write(solution, "test.yaml", path=temp_dir)
            new_solution = ct.Solution(output)

        assert compare_models(solution, new_solution)
