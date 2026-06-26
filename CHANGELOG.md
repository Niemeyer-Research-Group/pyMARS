# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Added

- Added a global `min-flame-speed` input option (default 0.05 m/s): a solved laminar flame speed at or below this floor is treated as a degenerate, non-physical result ("no flame"). Lower it when studying fuels with genuinely low flame speeds.
- Adds tests for `num_workers` > 1

### Changed

- Sampling workers now have the same behavior; the ignition sampling worker processes and removes the `.h5` files.

### Fixed

- Laminar flame reductions no longer abort when a candidate reduced model cannot sustain a flame. A failed or degenerate (negative/near-zero) flame solve is now treated as "no flame," so the reduced model is rejected via the error metric (mirroring how a non-igniting model is handled) instead of raising. A flame failure for the original/baseline model still raises so a missing baseline is caught.

## [1.2.0] - 2026-06-24

### Added

- Added laminar flame speed support: freely-propagating `FlameSimulation`, the `laminar-flame-conditions` input block, and laminar flame speed as an error metric, enabling reductions driven by flame data alone or combined with autoignition (based on work by Cailin Moore, #82)
- Added `BaseSimulation` abstract base class with shared gas setup and thermochemical-profile sampling, giving `IgnitionSimulation` and `FlameSimulation` a common interface
- Added automatic retention of fuel, oxidizer, and reactant species from the initial conditions so they are never removed during reduction (fixes #12; based on work by Katherine Bronstein, #80)
- Added `soln2yaml` module for writing Cantera `Solution` objects to YAML files, delegating to Cantera's built-in `write_yaml()` method
- Added YAML equivalents of all test asset models (`artificial-mechanism.yaml`, `model-third-bodies.yaml`, reduced GRI 3.0 models)
- Added conda recipe (`conda.recipe/meta.yaml`) updated for Python 3.10+ and Cantera 3.x
- Added GitHub Actions CI workflow with matrix testing on Python 3.10–3.14 (Linux) and Python 3.14 (macOS, Windows), with Coveralls coverage reporting
- Added GitHub Actions publish workflow for PyPI (OIDC trusted publishing) and Anaconda.org; conda publish triggers automatically on `v*` tags after CI passes
- Added GitHub Actions docs workflow that deploys to GitHub Pages after CI passes on `main`

### Changed

- Autoignition simulations now use Cantera mole reactors (`IdealGasMoleReactor` / `IdealGasConstPressureMoleReactor`) with an adaptive preconditioner, accelerating integration of large models (based on work by Anthony Walker, #84)
- Reduction drivers (`run_drg`, `run_drgep`, `run_pfa`, and sensitivity analysis) now accept laminar flame conditions and combine ignition and flame error metrics and sampled state data uniformly, regardless of simulation type
- Input files now require at least one of `autoignition-conditions` or `laminar-flame-conditions` (previously autoignition was mandatory)
- Renamed the `Simulation` class to `IgnitionSimulation` and consolidated the base, ignition, and flame simulation classes into a single `simulation.py`
- `trim()` now preserves the original model's transport model so reduced models remain usable for laminar flame simulations
- Minimum Python version raised from 3.6 to 3.10
- Updated Cantera dependency to 3.x; updated all API calls for Cantera 3.x (reaction type strings, stoichiometric coefficient property names, reactor constructor, phase/kinetics identifiers)
- Replaced `.cti` format throughout with Cantera's YAML format; `soln2cti.py` removed and replaced by `soln2yaml.py`
- Replaced `pkg_resources` with `pathlib` and `importlib.metadata` for Python 3.14 compatibility
- Switched build system from setuptools to hatchling
- Model file references in CLI, documentation, and examples updated from `.cti` to `.yaml`
- `compare_models` now uses `np.isclose()` for floating-point rate parameter comparisons and evaluates thermodynamic properties at reference temperatures instead of comparing raw polynomial coefficients
- `docs/conf.py` updated to use `importlib.metadata` for version detection; removed Travis CI references; updated Cantera intersphinx URL to 3.x docs

### Fixed

- Fixed bug in soln2ck.py where the reaction high rate was being used instead of the reaction low rate
- Fixed Cantera 3.x deprecation warnings: `reactor.thermo` → `reactor.phase`, `IdealGas`/`GasKinetics` → `ideal-gas`/`bulk`, reactor `clone=False`
- Removed broken xfail and skip-decorated tests that referenced undefined functions

## [1.1.0] - 2019-09-06

### Added

- Added checks to ensure any species in inputs are in the starting model
- Added ability to specify `ideal_gas` object name for reduction, if model in model file

### Changed

- Moved all reduction inputs to a YAML file, and updated associated docs
- Autoignition simulations now run to steady state by default, with optional
  options of `end-time` or `max-steps` to override

### Fixed

- Fixed bug in DRGEP overall interaction coefficient calculation with multiple conditions
- Removed erroneous message from ignition simulations about lack of convergence
- Fixed potenetial bug when printing CTI files with hyphens in species names
- Fixed bug that printed efficiency list for pressure-dependent reactions with explicit third body
- Fixed bug that retained reactions with explicit third bodies that were removed.


## [1.0.0] - 2019-06-21

### Added

- Added path flux analysis, so methods now include DRG, DRGEP, PFA, and sensitivity analysis
- Adds test suite for all current methods
- Adds Travis CI and AppVeyor continuous integration
- Added testing on Windows systems via Azure Pipelines
- Added printing support for chemically activated and Chebyshev reactions, and non-Troe falloff
- Adds Travis-based package deploy to PyPI and Anaconda
- Adds initial documentation site

### Changed

- Major code restructuring.
- Now supports Python 3.6+

### Fixed

- Fixed bugs in setup.py script and path imports
- Fixed model printing of Plog reactions

## [0.0.1] - 2017-04-01

### Added

- Added preliminary project files (README, CODE_OF_CONDUCT, etc.) without any actual source code


[Unreleased]: https://github.com/Niemeyer-Research-Group/pyMARS/compare/v1.2.0...HEAD
[1.2.0]: https://github.com/Niemeyer-Research-Group/pyMARS/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/Niemeyer-Research-Group/pyMARS/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/Niemeyer-Research-Group/pyMARS/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/Niemeyer-Research-Group/pyMARS/compare/v0.1.0...v1.0.0
[0.1.0]: https://github.com/Niemeyer-Research-Group/pyMARS/compare/v0.0.1...v0.1.0
