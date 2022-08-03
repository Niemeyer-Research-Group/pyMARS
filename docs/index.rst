.. pyMARS documentation master file, created by
   sphinx-quickstart on Wed May 29 10:17:58 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyMARS |release|
================

pyMARS is a Python package for reducing the size of chemical kinetic models;
the name stands for Python-based Model Automatic Reduction. pyMARS is licensed
under the permissive, open-source MIT license. The source code is publicly
available on GitHub at https://github.com/Niemeyer-Research-Group/pyMARS.

pyMARS takes models given in the standard Cantera or Chemkin file formats
(the latter is first converted to Cantera), along with a specified error limit
and optional lists of species to always retain, and automatically removes as
many species and reactions as possible using a specified method.

It uses provided conditions to sample thermochemical data for the reduction
analysis and calculate global combustion metrics for determining the error.
For example, typically most reductions will use autoignition simulations
to sample state data and determine autoignition delay over the condition
range desired for the reduced model.

pyMARS currently supports four methods for skeletal model reduction:

 1. Directed relation graph (DRG)
 2. Directed relation graph with error propagation (DRGEP)
 3. Path flux analysis (PFA)
 4. Sensitivity analysis (SA)

Citation
--------

If you use pyMARS, please cite it in any resulting publications. Thank you!


User's Guide
------------

.. toctree::
   :maxdepth: 2

   installation
   usage
   theory


Code API
--------

.. toctree::
   :maxdepth: 2
   :caption: Modules:

   drg
   drgep
   pfa
   pymars
   reduce_model
   sampling
   sensitivity_analysis
   simulation
   soln2ck
   soln2cti
   tools



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
