# pyMARS
<img src="logo/pymars-logo.png" align="right" width="200" />

[![DOI](https://zenodo.org/badge/51664233.svg)](https://zenodo.org/badge/latestdoi/51664233)
[![Build Status](https://travis-ci.org/Niemeyer-Research-Group/pyMARS.svg?branch=master)](https://travis-ci.org/Niemeyer-Research-Group/pyMARS)
[![Build Status](https://dev.azure.com/kyleniemeyer/pyMARS/_apis/build/status/Niemeyer-Research-Group.pyMARS?branchName=master)](https://dev.azure.com/kyleniemeyer/pyMARS/_build/latest?definitionId=1&branchName=master)
[![codecov](https://codecov.io/gh/Niemeyer-Research-Group/pyMARS/branch/master/graph/badge.svg)](https://codecov.io/gh/Niemeyer-Research-Group/pyMARS)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Code of Conduct](https://img.shields.io/badge/code%20of%20conduct-contributor%20covenant-green.svg)](http://contributor-covenant.org/version/1/4/)

Python-based (chemical kinetic) Model Automatic Reduction Software (pyMARS) implements multiple techniques for reducing the size and complexity of detailed chemical kinetic models.

An installation guide, usage examples, theory details, and API docs are provided in the online documentation: https://Niemeyer-Research-Group.github.io/pyMARS/

pyMARS currently consists of four methods for model reduction:

 1. Directed relation graph (DRG)
 2. Directed relation graph with error propagation (DRGEP)
 3. Path flux analysis (PFA)
 4. Sensitivity analysis (SA)

Sensitivity analysis may be run following one of the first three methods, or directly on the starting
model; however, its computational expense is high, and applying this method alone is not recommended.

## Installation

pyMARS supports Python 3.6 and 3.7, and can be installed easily using conda:

    conda install -c niemeyer-research-group pymars

## Usage

For detailed usage examples, see the [online documentation](https://Niemeyer-Research-Group.github.io/pyMARS/).
Once installed, the list of options can be found with:

    pymars --help

pyMARS requires models in the Cantera format. However, running pyMARS with a CHEMKIN file will convert it 
into a Cantera file. pyMARS also provides the `--convert` option to convert a given model to/from 
the CHEMKIN format.

## Citation

Please refer to the CITATION file for information about citing pyMARS when used in a scholarly work.

If you use this package as part of a scholarly publication, please consider citing the appropriate 
theory/method papers in addition to the software itself.

## License

pyMARS is released under the MIT license; see LICENSE for details.

## Code of Conduct

To ensure an open and welcoming community, pyMARS adheres to a code of conduct adapted from the [Contributor Covenant](http://contributor-covenant.org) code of conduct.

Please adhere to this code of conduct in any interactions you have in the pyMARS community. It is strictly enforced on all official PyKED repositories, websites, and resources. If you encounter someone violating these terms, please let the project lead (@kyleniemeyer) know via email at <kyle.niemeyer@gmail.com> and we will address it as soon as possible.
