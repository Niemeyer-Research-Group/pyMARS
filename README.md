# pyMARS

[![Code of Conduct](https://img.shields.io/badge/code%20of%20conduct-contributor%20covenant-green.svg)](http://contributor-covenant.org/version/1/4/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Python-based (chemical kinetic) Model Automatic Reduction Software (MARS), which consists of multiple techniques for reducing the size and complexity of detailed chemical kinetic models.  pyMARS requires mechanisms to be stored in the Cantera format to be reduced.  However, running pyMARS with a CHEMKIN file will convert it into a Cantera file that can then be used.  The --convert option can then be used to turn the reduced model back into a CHEMKIN file.  

pyMARS currently consists of two methods for model reduction:

 1. Directed relation graph (DRG)
 2. Directed relation graph with error propagation (DRGEP)

Both of these methods are documented in liturature.  Additional reduction stages, including isomer lumping and CSP-based quasi-steady-state (QSS) species reduction, are currently under development and testing.

See the following publications for more detail:

 * KE Niemeyer, CJ Sung, and MP Raju. Skeletal mechanism generation for surrogate fuels using directed relation graph with error propagation and sensitivity analysis. *Combust. Flame*, 157(9):1760--1770, 2010. doi:[10.1016/j.combustflflame.2009.12.022](http://dx.doi.org/10.1016/j.combustflflame.2009.12.022)
 * KE Niemeyer and CJ Sung. On the importance of graph search algorithms for DRGEP-based mechanism reduction methods. *Combust. Flame*, 158(8):1439--1443, 2011. doi:[10.1016/j.combustflflame.2010.12.010](http://dx.doi.org/10.1016/j.combustflflame.2010.12.010).
 * KE Niemeyer and CJ Sung. Mechanism reduction for multicomponent surrogates: A case study using toluene reference fuels.  *Combust. Flame*, in press, 2014. doi:[10.1016/j.combustflame.2014.05.001](http://dx.doi.org/10.1016/j.combustflame.2014.05.001)
 * TF Lu and CK Law. *Combustion and Flame*, 154:153--163, 2008. doi:[10.1016/j.combustflame.2007.11.013](http://dx.doi.org/10.1016/j.combustflame.2007.11.013)

## Usage

While running pyMARS, make sure that you are using an up to date version of Python 2.7.  Python 3 compatibility is currently being developed.  

To install:
     `python setup.py install`

pyMARS is called from terminal via `__main__.py`
__main__.py can be found in the pyMARS directory.

example:
    `python __main__.py --file ../example_files/gri30.cti --run_drgep --conditions ../example_files/example_input.txt --species N2,CO2,H2O --error 5 --target CH4,O2 --plot --write_ai_times`

This will run pyMARS with the gri30.cti mechanism with the initial conditions listed in the example file.  pyMARS will give the autoignition times as well as a plot of the autignition simulation for each inital condition.  Then, pyMARS will reduce the mechanism using the DRGEP method with the given target species until the error reaches 5 percent.  The species listed under species will not be removed from the model under any circumstance. 

## Options

Running pyMARS without any options will show a list of all possible options.  The options supported by pyMARS are:

  * --file: The file given for this option should be the mechanism for pyMARS to act on.  
  * --conditions: A file that contains all of the initial conditions for autoignition.  See the example given for formatting.  
  * --convert: Calling this option will convert the given Cantera file to the CHEMKIN format.  
  * --thremo: This option holds the thermo data file if your CHEMKIN model has one.  
  * --transport:  This option holds the transport data file if your CHEMKIN model has one.  
  * --species:  Any species included in this comma seperated list will not be removed from the model no matter what.
  * --error:  This value is the maximum error level that will be allowed for the reduced model.  Error percentage is calulated by comparing autoignition delays from the original and reduced models.  
  * --target: Comma seperated list of target species for model reduction.  
  * --run_drg: This option will run the DRG method for model reduction on the given model.  It requires a given error and target species through the --error and --target options.  
  * --run_drgep: This option will run the DRGEP method for model reduction on the given model.  It requires a given error and target species through the --error and --target options. 
  * --plot: Plots the autoignition simulations for all of the initial conditions for the original mechanism.  
  * --points: Prints the range of the sampling points on the screen. 
  * --writecsv: This option will create csv files containing tempuratures and times from the autoignition simulations for each inital condition for the original mechanism.  
  * --writehdf5: Writes hdf5 files for the autoignition simulations for each inital condition. 
  * --write_ai_times: Creates a file containing the autoignition times for the original mechanism.  

## License

pyMARS is released under the MIT license, see LICENSE for details.

If you use this package as part of a scholarly publication, please cite the following papers in addition to this resource:

 * KE Niemeyer, CJ Sung, and MP Raju. Skeletal mechanism generation for surrogate fuels using directed relation graph with error propagation and sensitivity analysis. *Combust. Flame*, 157(9):1760--1770, 2010. doi:[10.1016/j.combustflflame.2009.12.022](http://dx.doi.org/10.1016/j.combustflflame.2009.12.022)
 * KE Niemeyer and CJ Sung. On the importance of graph search algorithms for DRGEP-based mechanism reduction methods. *Combust. Flame*, 158(8):1439--1443, 2011. doi:[10.1016/j.combustflflame.2010.12.010](http://dx.doi.org/10.1016/j.combustflflame.2010.12.010).
 * KE Niemeyer and CJ Sung. Mechanism reduction for multicomponent surrogates: A case study using toluene reference fuels.  *Combust. Flame*, in press, 2014. doi:[10.1016/j.combustflame.2014.05.001](http://dx.doi.org/10.1016/j.combustflame.2014.05.001)

## Code of Conduct

In order to have a more open and welcoming community, pyMARS adheres to a code of conduct adapted from the [Contributor Covenant](http://contributor-covenant.org) code of conduct.

Please adhere to this code of conduct in any interactions you have in the pyMARS community. It is strictly enforced on all official PyKED repositories, websites, and resources. If you encounter someone violating these terms, please let the project lead (@kyleniemeyer) know via email at <kyle.niemeyer@gmail.com> and we will address it as soon as possible.
