# pyMARS

[![DOI](https://zenodo.org/badge/51664233.svg)](https://zenodo.org/badge/latestdoi/51664233)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Code of Conduct](https://img.shields.io/badge/code%20of%20conduct-contributor%20covenant-green.svg)](http://contributor-covenant.org/version/1/4/)

Python-based (chemical kinetic) Model Automatic Reduction Software (pyMARS) implements multiple techniques for reducing the size and complexity of detailed chemical kinetic models.

pyMARS requires models in the Cantera format. However, running pyMARS with a CHEMKIN file will convert it into a Cantera file, and provides the `--convert` option to convert the reduced model back into the CHEMKIN format.

pyMARS currently consists of four methods for model reduction:

 1. Directed relation graph (DRG)
 2. Directed relation graph with error propagation (DRGEP)
 3. Sensitivity analysis (SA)
 4. Path flux analysis (PFA)

These methods are documented in the literature. Sensitivity analysis must be performed after the completion of another method. Additional reduction stages are under development and testing.

See the following publications for more detail:

 * KE Niemeyer, CJ Sung, and MP Raju. Skeletal mechanism generation for surrogate fuels using directed relation graph with error propagation and sensitivity analysis. *Combust. Flame*, 157(9):1760--1770, 2010. doi:[10.1016/j.combustflflame.2009.12.022](https://doi.org/10.1016/j.combustflflame.2009.12.022)
 * KE Niemeyer and CJ Sung. On the importance of graph search algorithms for DRGEP-based mechanism reduction methods. *Combust. Flame*, 158(8):1439--1443, 2011. doi:[10.1016/j.combustflflame.2010.12.010](https://doi.org/10.1016/j.combustflflame.2010.12.010).
 * KE Niemeyer and CJ Sung. Mechanism reduction for multicomponent surrogates: A case study using toluene reference fuels. *Combust. Flame*, in press, 2014. doi:[10.1016/j.combustflame.2014.05.001](https://doi.org/10.1016/j.combustflame.2014.05.001)
 * TF Lu and CK Law. *Combustion and Flame*, 154:153--163, 2008. doi:[10.1016/j.combustflame.2007.11.013](https://doi.org/10.1016/j.combustflame.2007.11.013)

## Usage

While running pyMARS, make sure that you are using an up to date version of Python 3 with Cantera installed. pyMARS no longer works with Python 2.7.

To install: `python setup.py install`

pyMARS is called from terminal via `pyMARS.py` which can be found in the pyMARS directory.

example:
```
    python __main__.py -m ../example_files/gri30.cti --conditions ../example_files/example_input_file.txt -e 5 --method DRGEP --targets CH4 O2 --retained_species CH4 O2 N2 CO2 H2O
```

This will run pyMARS with the GRI Mech 3.0 model with the initial conditions listed in the example file. pyMARS will record data from the autoignition simulation for each initial condition. Then, pyMARS will reduce the model using the DRGEP method with the given target species until the error reaches 5%. The species listed under species will not be removed from the model under any circumstance.

## Options

Running pyMARS without any options will show a list of all possible options. The options supported by pyMARS are:

  * `--file`: The file given for this option should be the mechanism for pyMARS to act on.
  * `--conditions`: A file that contains all initial conditions for autoignition. See the example given for formatting.
  * `--convert`: Calling this option will convert the given Cantera file to the CHEMKIN format.
  * `--thremo`: This option holds the thermo data file if your CHEMKIN model has one.
  * `--transport`: This option holds the transport data file if your CHEMKIN model has one.
  * `--species`: Any species included in this comma seperated list will not be removed from the model no matter what.
  * `--error`: This value is the maximum error level that will be allowed for the reduced model. Error percentage is calulated by comparing autoignition delays from the original and reduced models.
  * `--target`: Comma-seperated list of target species for model reduction.
  * `--run_drg`: This option will run the DRG method for model reduction on the given model. It requires a given error and target species through the --error and --target options.
  * `--run_drgep`: This option will run the DRGEP method for model reduction on the given model. It requires a given error and target species through the --error and --target options.
  * `--run_sa`: Run a sensitivity analysis after completing another reduction method.
  * `--ep_star`: A float to be used as the ep star value for the sensitivity analysis.

## Citation

Please refer to the CITATION file for information about citing pyMARS when used in a scholarly work.

If you use this package as part of a scholarly publication, it may be appropriate to cite the following papers in addition to this resource:

 * KE Niemeyer, CJ Sung, and MP Raju. Skeletal mechanism generation for surrogate fuels using directed relation graph with error propagation and sensitivity analysis. *Combust. Flame*, 157(9):1760--1770, 2010. doi:[10.1016/j.combustflflame.2009.12.022](https://doi.org/10.1016/j.combustflflame.2009.12.022)
 * KE Niemeyer and CJ Sung. On the importance of graph search algorithms for DRGEP-based mechanism reduction methods. *Combust. Flame*, 158(8):1439--1443, 2011. doi:[10.1016/j.combustflflame.2010.12.010](https://doi.org/10.1016/j.combustflflame.2010.12.010).
 * KE Niemeyer and CJ Sung. Mechanism reduction for multicomponent surrogates: A case study using toluene reference fuels. *Combust. Flame*, in press, 2014. doi:[10.1016/j.combustflame.2014.05.001](https://doi.org/10.1016/j.combustflame.2014.05.001)

## License

pyMARS is released under the MIT license, see LICENSE for details.

## Code of Conduct

To have a more open and welcoming community, pyMARS adheres to a code of conduct adapted from the [Contributor Covenant](http://contributor-covenant.org) code of conduct.

Please adhere to this code of conduct in any interactions you have in the pyMARS community. It is strictly enforced on all official PyKED repositories, websites, and resources. If you encounter someone violating these terms, please let the project lead (@kyleniemeyer) know via email at <kyle.niemeyer@gmail.com> and we will address it as soon as possible.
