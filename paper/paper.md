---
title: 'pyMARS: automatically reducing chemical kinetic models in Python'
tags:
  - chemical kinetics
  - model reduction
  - Python
authors:
  - name: Phillip O. Mestas III
    orcid: 0000-0003-4379-3592
    affiliation: 1
  - name: Parker Clayton
    affiliation: 1
  - name: Kyle E. Niemeyer
    orcid: 0000-0003-4425-7097
    affiliation: 1
affiliations:
  - name: School of Mechanical, Industrial, and Manufacturing Engineering, Oregon State University, Corvallis, OR USA 97331
    index: 1
date: 22 June 2019
bibliography: paper.bib
---

# Summary

Chemically reacting fluid flows occur in a diverse range of scientific and engineering fields, including
combustion and fire, atmopsheric and oceanic fluid flows, electrochemical devices, heterogeneous catalysis,
materials processing, and astrophysical fluid dynamics. Numerical simulations of reacting fluid flows rely
on accurate chemical kinetic models to describe the participating chemical species and elementary reactions
through which they interact. However, as models grow in detail they also grow in size, adding more species
and reactions to capture intermediates and pathways. In the field of combustion, kinetic models for molecules
relevant to transportation fuels (e.g., gasoline, diesel, jet fuel) can contain thousands of species and tens
of thousands of elementary reactions [@lu:2009]. Incorporating such models into multidimensional
computational fluid dynamics simulations is practically impossible, due to the associated computational expense.

pyMARS, which stands for "Python-based Model Automatic Reduction Software", is a software package
that implements and applies literature methods for reducing chemical kinetic models, particularly
targeted at combustion applications. pyMARS currently implements four "skeletal" reduction methods
that identify and remove unimportant species and reactions: directed relation graph (DRG) [@lu:2005; @lu:2006a; @lu:2006b], 
directed relation graph with error propagation (DRGEP) [@pepiotdesjardins:2008; @niemeyer:2011],
path flux analysis (PFA) [@sun:2010], and
sensitivity analysis [@sankaran:2007; @zheng:2007; @niemeyer:2010; @niemeyer:2015].
pyMARS succeeds the earlier Fortran-based MARS package [@niemeyer:2010; @niemeyer:2014; @niemeyer:2015],
which used an in-house modified version of the proprietary Chemkin III library [@Kee1996].

# Background and features

DRG, DRGEP, and PFA represent the kinetic system as a graph, where nodes are species and directed, weighted
edges represent the dependence of one species on another, through their participation in reactions.
All three methods define interaction coefficients that approximate the error that would be induced in the
overall production/consumption of one species if the other was removed from the model.
To eliminate species, a cutoff threshold (e.g., 0.01--0.1) is applied to the system to eliminate
unimportant connections or species. The methods differ in their definition of these coefficients
and whether/how indirect relationships between species play a role. For all three methods, pyMARS
iteratively increases the cutoff threshold until reaching a user-specified error limit.

Sensitivity analysis can be applied directly to a starting model, but generally it should be informed
by and follow a graph-based method such as DRG or DRGEP (i.e., DRGASA and DRGEPSA) due to the high cost
associated with using such a brute-force approach on a large model. pyMARS implements two sensitivity
analysis approaches: "initial" and "greedy" [@niemeyer:2015]. The initial algorithm finds the error
induced by individually eliminating all species under consideration, one-by-one, then removes
species in ascending order of induced error until reaching the limit. In contrast, the greedy algorithm
reevaluates the induced error of remaining species at each step, ensuring that it uses the most-current
information; this increases computational expense significantly, but generates a smaller reduced model.

Required inputs for pyMARS to perform model reduction include the starting chemical kinetic model
in the standard Cantera [@cantera] or Chemkin [@Kee1996] formats (it first converts the latter
to the former), a maximum error limit to constrain the reduced model, and a YAML file with initial
conditions for homogeneous autoignition simulations. These simulations are used both to obtain
ignition delay values for gauging the error of a candidate model, and also to sample thermochemical
state data for the reduction methods (where 20 points are sampled during the ignition temperature rise).
A graph-based reduction method can be specified, the ``--run_sensitivity_analysis`` flag can be
given to perform standalone sensitivity analysis, or both can be given to perform DRG/DRGEP/PFA-informed
sensitivity analysis. Additional inputs include target species for DRG/DRGEP/PFA and optionally
specifying a list of species to always retain.

Additional features include:

- Simulations can be parallelized via the ``multiprocessing`` module by adding the ``--num_threads``
option and specifying a value greater than one. (Including the option alone will lead to pyMARS
using the available number of threads minus one.)
- To save (potentially significant) time when performing multiple reductions for the same model,
pyMARS saves and automatically reuses sampled autoignition data when possible.
- pyMARS provides a conversion functionality between the Cantera and Chemkin model formats,
accessible by passing the ``--convert`` option.

pyMARS relies on the Cantera suite [@cantera] to handle the chemical kinetics
and perform autoignition simulations; it temporarily stores full simulation results
as HDF5 files using PyTables [@pytables]. It uses PyYAML [@pyyaml] to parse simulation
input files, and NumPy [@vanderWalt2011] arrays to store and manipulate data.
Graph construction and searching rely on NetworkX [@networkx_paper; @networkx].

# Future work

Future work includes adding more reduction stages, including unimportant reaction elimination
and incorporating the quasi-steady-state approximation for species [@niemeyer:2015].
In addition, we plan to add additional combustion phenomena for sampling and error
evaluation, including one-dimensional laminar flame and perfectly stirred reactor
simulations.

# Acknowledgements

This material is based upon work supported by the National Science Foundation
under grant OAC-1535065.

# References
