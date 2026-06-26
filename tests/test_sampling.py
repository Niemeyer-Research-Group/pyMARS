"""Tests the sampling module in pyMARS"""

import pathlib

import pytest
import numpy as np
import cantera as ct

from pymars import sampling
from pymars.sampling import (
    parse_ignition_inputs,
    parse_flame_inputs,
    sample,
    sample_metrics,
    InputIgnition,
    InputLaminarFlame,
)

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
    return str(pathlib.Path(__file__).parent / file)


class TestCheckInputs:
    def test_good_example(self):
        """Tests (good) sample input file."""
        inputs = [
            {
                "kind": "constant volume",
                "pressure": 1.0,
                "temperature": 1000.0,
                "fuel": {"CH4": 1.0},
                "oxidizer": {"O2": 1.0, "N2": 3.76},
                "equivalence-ratio": 1.0,
            },
            {
                "kind": "constant volume",
                "pressure": 1.0,
                "temperature": 1200.0,
                "fuel": {"CH4": 1.0},
                "oxidizer": {"O2": 1.0, "N2": 3.76},
                "equivalence-ratio": 1.0,
            },
        ]

        conditions = parse_ignition_inputs("gri30.yaml", inputs)
        for item in conditions:
            assert isinstance(item, InputIgnition)

    def test_good_example_alternate(self):
        """Tests correct sample input file with alternate values."""
        inputs = [
            {
                "kind": "constant pressure",
                "pressure": 1.0,
                "temperature": 1000.0,
                "end time": 10.0,
                "reactants": {"CH4": 1.0, "O2": 1.0, "N2": 3.76},
            },
        ]

        conditions = parse_ignition_inputs("gri30.yaml", inputs)
        for item in conditions:
            assert isinstance(item, InputIgnition)

    @pytest.mark.parametrize(
        "key",
        ["kind", "pressure", "temperature", "fuel", "oxidizer", "equivalence-ratio"],
    )
    def test_missing_keys(self, key):
        """Tests correct errors for missing required keys."""
        case = [
            {
                "kind": "constant volume",
                "pressure": 1.0,
                "temperature": 1000.0,
                "end-time": 10.0,
                "fuel": {"CH4": 1.0},
                "oxidizer": {"O2": 1.0, "N2": 3.76},
                "equivalence-ratio": 1.0,
            }
        ]
        del case[0][key]

        with pytest.raises(AssertionError):
            parse_ignition_inputs("gri30.yaml", case)

    def test_bad_fuel_oxidizer_value(self):
        """Tests correct errors for improper value."""
        case = [
            {
                "kind": "constant volume",
                "pressure": 1.0,
                "temperature": 1000,
                "end-time": 10.0,
                "fuel": {"CH4": 0.0},
                "oxidizer": {"O2": 1.0, "N2": 3.76},
                "equivalence-ratio": 1.0,
            }
        ]
        with pytest.raises(AssertionError):
            parse_ignition_inputs("gri30.yaml", case)

        case = [
            {
                "kind": "constant volume",
                "pressure": 1.0,
                "temperature": 1000,
                "end-time": 10.0,
                "fuel": {"CH4": 1.0},
                "oxidizer": {"O2": 0.0, "N2": 3.76},
                "equivalence-ratio": 1.0,
            }
        ]
        with pytest.raises(AssertionError):
            parse_ignition_inputs("gri30.yaml", case)

    def test_bad_species(self):
        """Tests raising error for species not in model."""
        case = [
            {
                "kind": "constant volume",
                "pressure": 1.0,
                "temperature": 1000,
                "end-time": 10.0,
                "fuel": {"C4H10": 0.0},
                "oxidizer": {"O2": 1.0, "N2": 3.76},
                "equivalence-ratio": 1.0,
            }
        ]
        with pytest.raises(AssertionError):
            parse_ignition_inputs("gri30.yaml", case)

        case = [
            {
                "kind": "constant volume",
                "pressure": 1.0,
                "temperature": 1000,
                "end-time": 10.0,
                "fuel": {"CH4": 1.0},
                "oxidizer": {"O2": 0.0, "HE": 3.76},
                "equivalence-ratio": 1.0,
            }
        ]
        with pytest.raises(AssertionError):
            parse_ignition_inputs("gri30.yaml", case)

        case = [
            {
                "kind": "constant volume",
                "pressure": 1.0,
                "temperature": 1000,
                "end-time": 10.0,
                "reactants": {"C4H10": 1.0, "O2": 1.0, "N2": 3.76},
            }
        ]
        with pytest.raises(AssertionError):
            parse_ignition_inputs("gri30.yaml", case)


class TestParseFlameInputs:
    def test_good_example_equivalence_ratio(self):
        """Tests a valid flame input using an equivalence ratio."""
        inputs = [
            {
                "pressure": 1.0,
                "temperature": 300.0,
                "width": 0.03,
                "fuel": {"CH4": 1.0},
                "oxidizer": {"O2": 1.0, "N2": 3.76},
                "equivalence-ratio": 1.0,
            },
        ]
        conditions = parse_flame_inputs("gri30.yaml", inputs)
        assert len(conditions) == 1
        for item in conditions:
            assert isinstance(item, InputLaminarFlame)
            assert item.width == 0.03

    def test_good_example_reactants(self):
        """Tests a valid flame input using a reactant list."""
        inputs = [
            {
                "pressure": 1.0,
                "temperature": 300.0,
                "reactants": {"CH4": 1.0, "O2": 1.0, "N2": 3.76},
            },
        ]
        conditions = parse_flame_inputs("gri30.yaml", inputs)
        for item in conditions:
            assert isinstance(item, InputLaminarFlame)
            # width should default when not specified
            assert item.width == 0.1

    @pytest.mark.parametrize(
        "key",
        ["pressure", "temperature", "fuel", "oxidizer", "equivalence-ratio"],
    )
    def test_missing_keys(self, key):
        """Tests correct errors for missing required keys."""
        case = [
            {
                "pressure": 1.0,
                "temperature": 300.0,
                "width": 0.03,
                "fuel": {"CH4": 1.0},
                "oxidizer": {"O2": 1.0, "N2": 3.76},
                "equivalence-ratio": 1.0,
            }
        ]
        del case[0][key]
        with pytest.raises(AssertionError):
            parse_flame_inputs("gri30.yaml", case)

    def test_bad_width(self):
        """Tests correct error for non-positive width."""
        case = [
            {
                "pressure": 1.0,
                "temperature": 300.0,
                "width": -0.01,
                "fuel": {"CH4": 1.0},
                "oxidizer": {"O2": 1.0, "N2": 3.76},
                "equivalence-ratio": 1.0,
            }
        ]
        with pytest.raises(AssertionError):
            parse_flame_inputs("gri30.yaml", case)

    def test_bad_species(self):
        """Tests raising error for species not in model."""
        case = [
            {
                "pressure": 1.0,
                "temperature": 300.0,
                "reactants": {"C4H10": 1.0, "O2": 1.0, "N2": 3.76},
            }
        ]
        with pytest.raises(AssertionError):
            parse_flame_inputs("gri30.yaml", case)


class TestFlameSampling:
    """Exercises the laminar flame branches of sample / sample_metrics."""

    def _hydrogen_flame(self):
        return InputLaminarFlame(
            pressure=1.0,
            temperature=300.0,
            equivalence_ratio=1.0,
            fuel={"H2": 1.0},
            oxidizer={"O2": 1.0, "N2": 3.76},
            width=0.03,
        )

    def test_sample_metrics_flame_only(self):
        """sample_metrics should return the laminar flame speed as the metric."""
        metrics = sample_metrics(
            "h2o2.yaml", [], flame_conditions=[self._hydrogen_flame()], num_threads=1
        )
        assert metrics.shape == (1,)
        assert 0.5 < metrics[0] < 5.0

    def test_sample_flame_data_and_reuse(self, tmp_path, monkeypatch):
        """sample should write flame metric/data files and reuse them on a second call."""
        monkeypatch.setitem(
            sampling.data_files, "data_flame", str(tmp_path / "flame_data.dat")
        )
        monkeypatch.setitem(
            sampling.data_files, "output_flame", str(tmp_path / "flame_out.txt")
        )

        flame_conditions = [self._hydrogen_flame()]
        metrics, data = sample(
            "h2o2.yaml", [], flame_conditions=flame_conditions, num_threads=1
        )

        gas = ct.Solution("h2o2.yaml")
        assert metrics.shape == (1,)
        assert 0.5 < metrics[0] < 5.0
        assert data.shape == (
            sampling.IgnitionSimulation.num_sample_points,
            2 + gas.n_species,
        )
        assert (tmp_path / "flame_data.dat").is_file()
        assert (tmp_path / "flame_out.txt").is_file()

        # second call should reuse the saved samples and give identical results
        metrics_reuse, data_reuse = sample(
            "h2o2.yaml", [], flame_conditions=flame_conditions, num_threads=1
        )
        assert np.allclose(metrics_reuse, metrics)
        assert np.allclose(data_reuse, data)

    def test_sample_combined_ignition_and_flame(self, tmp_path, monkeypatch):
        """Combined sampling concatenates ignition metrics/data before flame ones."""
        for key in (
            "data_ignition",
            "output_ignition",
            "data_flame",
            "output_flame",
        ):
            monkeypatch.setitem(sampling.data_files, key, str(tmp_path / key))

        ignition_conditions = [
            InputIgnition(
                kind="constant volume",
                pressure=1.0,
                temperature=1200.0,
                equivalence_ratio=1.0,
                fuel={"H2": 1.0},
                oxidizer={"O2": 1.0, "N2": 3.76},
            )
        ]
        flame_conditions = [self._hydrogen_flame()]

        metrics, data = sample(
            "h2o2.yaml",
            ignition_conditions,
            flame_conditions=flame_conditions,
            num_threads=1,
        )

        gas = ct.Solution("h2o2.yaml")
        n_points = sampling.IgnitionSimulation.num_sample_points

        # one ignition metric then one flame metric, in that order
        assert metrics.shape == (2,)
        assert 0.0 < metrics[0] < 0.5  # ignition delay (s) in the first slot
        assert 0.5 < metrics[1] < 5.0  # flame speed (m/s) in the second slot

        # sampled data is the ignition states stacked on top of the flame states
        assert data.shape == (2 * n_points, 2 + gas.n_species)

    def test_sample_metrics_no_flame_is_zero(self, monkeypatch):
        """A reduced model whose flame fails to solve yields a 0.0 metric.

        ``sample_metrics`` (the reduced-model path) must not crash on an
        undetectable flame; the 0.0 metric then drives a 100% error in
        ``calculate_error``, so the model is rejected during reduction.
        """
        from pymars.simulation import FlameSimulation

        def force_no_flame(_self):
            raise ct.CanteraError("forced failure: no flame detected")

        # num_threads=1 so the monkeypatch applies (it would not reach subprocesses)
        monkeypatch.setattr(FlameSimulation, "_solve_flame", force_no_flame)
        metrics = sample_metrics(
            "h2o2.yaml", [], flame_conditions=[self._hydrogen_flame()], num_threads=1
        )

        assert metrics.shape == (1,)
        assert metrics[0] == 0.0

    @pytest.mark.slow
    def test_parallel_matches_serial(self):
        """The multiprocessing path (num_threads > 1) matches serial results.

        Exercises the ``multiprocessing.Pool`` dispatch end to end, including
        pickling the simulation objects and results across the process boundary.
        """
        flame_conditions = [self._hydrogen_flame(), self._hydrogen_flame()]
        serial = sample_metrics(
            "h2o2.yaml", [], flame_conditions=flame_conditions, num_threads=1
        )
        parallel = sample_metrics(
            "h2o2.yaml", [], flame_conditions=flame_conditions, num_threads=2
        )
        assert serial.shape == (2,)
        assert parallel.shape == (2,)
        assert np.allclose(serial, parallel)


class TestIgnitionSampling:
    """Exercises the autoignition branch of sample / sample_metrics."""

    def _hydrogen_ignition(self):
        return InputIgnition(
            kind="constant volume",
            pressure=1.0,
            temperature=1200.0,
            equivalence_ratio=1.0,
            fuel={"H2": 1.0},
            oxidizer={"O2": 1.0, "N2": 3.76},
        )

    def test_sample_metrics_integration_failure_is_zero(self, monkeypatch):
        """A reduced model whose ignition integration fails yields a 0.0 metric.

        Reproduces the issue #69 scenario: ``sample_metrics`` (the reduced-model
        path) must not crash when the integrator fails (CVODES error); the 0.0
        metric then drives a 100% error in ``calculate_error``, so the model is
        rejected during reduction.
        """
        from pymars.simulation import IgnitionSimulation

        def force_failure(_self):
            raise ct.CanteraError("forced failure: CVODES error encountered")

        # num_threads=1 so the monkeypatch applies (it would not reach subprocesses)
        monkeypatch.setattr(IgnitionSimulation, "_step", force_failure)
        metrics = sample_metrics(
            "h2o2.yaml", [self._hydrogen_ignition()], num_threads=1
        )

        assert metrics.shape == (1,)
        assert metrics[0] == 0.0

    @pytest.mark.slow
    def test_parallel_matches_serial(self):
        """The ignition multiprocessing path matches serial results.

        Verifies ``IgnitionSimulation`` objects survive the spawn-pool pickle
        round-trip; a live reactor net cannot be pickled and is rebuilt per
        worker via ``setup_case()``.
        """
        ignition_conditions = [self._hydrogen_ignition() for _ in range(2)]
        serial = sample_metrics("h2o2.yaml", ignition_conditions, num_threads=1)
        parallel = sample_metrics("h2o2.yaml", ignition_conditions, num_threads=2)
        assert serial.shape == (2,)
        assert parallel.shape == (2,)
        assert np.allclose(serial, parallel)


def _double_worker(job):
    """Trivial picklable worker for testing the parallel dispatch in isolation."""
    _sim, idx = job
    return {idx: idx * 2}


class TestRunWorkers:
    """The shared dispatch helper merges results identically, serial or parallel."""

    def test_parallel_matches_serial(self):
        simulations = [[None, idx] for idx in range(4)]
        serial = sampling._run_workers(simulations, _double_worker, num_threads=1)
        parallel = sampling._run_workers(simulations, _double_worker, num_threads=2)
        assert serial == parallel == {0: 0, 1: 2, 2: 4, 3: 6}
