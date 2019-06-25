.. _sec_usage:

=====
Usage
=====

.. toctree::
    :maxdepth: 1

Once installed, pyMARS can be run as an executable. Change to the directory
where your input files are located, then run::

    pymars [options]

The following options are available, and can also be seen by using the
``-h`` or ``--help`` option:

.. code-block:: bash

     -h, --help:
        show this help message and exit
     -m, --model:
        input model filename (e.g., mech.cti)
     -e, --error:
        Maximum error percentage for the reduced model
     --method {DRG,DRGEP,PFA}
        skeletal reduction method to use
     --conditions:
        File with list of autoignition initial conditions
     --targets:
        List of target species (e.g., "CH4 O2")
     --retained_species:
        List of non-target species to always retain (e.g., "N2 Ar")
     --sensitivity_analysis:
        Run sensitivity analysis after completing another method
     --sensitivity_type {initial, greedy}
        Sensitivity analysis method to use
     --upper_threshold:
        Upper threshold value used to determine species for sensitivity analysis
     --path:
        Path to directory for writing files
     --num_threads:
        Number of CPU cores to use for running simulations in parallel
     --convert:
        Convert files between Cantera and Chemkin formats (.cti <=> .inp)
     --thermo:
        thermodynamic data filename (only necessary for Chemkin files)
     --transport:
        transport data filename (only necessary for Chemkin files)

Example reduction
=================

To apply the DRGEP method for reducing the GRI Mech 3.0 model (included with
Cantera), with a maximum error of 10%, using the main reactants as targets
(CH\ :sub:`4` and O\ :sub:`2`\ ) and always retaining N\ :sub:`2`\ ,
execute this command:

.. code-block:: bash

    pymars -m gri30.cti --method DRGEP --error 10 \
      --conditions input.yaml --targets CH4 O2 --retained_species N2

Note that you need to create a file ``input.yaml`` with a list of autoignition
initial conditions for the reduction; download an
:download:`example file <example_input_file.yaml>` or see the
:ref:`autoignition` section for a detailed explanation.

**Parallelization:** If you have a large number of initial conditions, the
reduction may be sped up by parallelizing the associated simulations over
multiple CPU cores. You can do this by adding ``--num_threads N``, where ``N``
is the desired number of cores. If you just specify ``--num_threads``, then
pyMARS will use the available numbera of cores minus one.

(pyMARS does not currently support distributed memory parallelization, meaning
across multiple nodes that do not share the same memory.)

**Sensitivity analysis:** To perform sensitivity analysis following DRGEP,
add the option ``--sensitivity_analysis``. Select the type of sensitivity analysis
with the ``--sensitivity_method`` option. Generally it is a good idea to specify an upper
threshold value of 0.2 to ensure important species are not evaluated; do this
by adding the option ``--upper_threshold 0.2``. (In the :ref:`sec_theory`
guide, this variable is referred to as :math:`\epsilon^*`.)

Then, the full command for performing a reduction using DRGEPSA with
parallelized simalations would be:

.. code-block:: bash

   pymars -m gri30.cti --method DRGEP --error 10 \
      --conditions input.yaml --targets CH4 O2 --retained_species N2 \
      --num_threads --sensitivity_analysis --sensitvity_method initial --upper_threshold 0.2

.. _autoignition:

Autoignition conditions
=======================

pyMARS currently uses autoignition simulations to sample thermochemical data
for the reduction and to calculate ignition delays for measuring error of
candidate reduced models.

Initial conditions need to be provided for performing these simulations,
in a YAML file specified with the ``--conditions`` option.
These files look like:

.. code-block:: yaml

    - kind: constant volume
      pressure: 1.0
      temperature: 1000
      end-time: 10.0
      fuel:
          CH4: 1.0
      oxidizer:
          O2: 1.0
          N2: 3.76
      equivalence-ratio: 1.0

    - kind: constant volume
      pressure: 1.0
      temperature: 1200
      end-time: 10.0
      fuel:
          CH4: 1.0
      oxidizer:
          O2: 1.0
          N2: 3.76
      equivalence-ratio: 0.5

This example specifies two constant-volume ignition cases, both with initial
pressures of 1 atm and with initial temperatures of 1000 K and 1200 K. Both
cases also specify maximum end integration times of 10 seconds.

These examples specify the intitial reactant mixture with ``fuel`` and
``oxidizer`` fields, and the mole fractions of their respective compositions,
followed by the equivalence ratio.

In addition to constant-volume autoignition, constant-pressure simulations
can be specified by changing ``kind``  to ``constant pressure``.
Initial reactant composition can also be specified directly without using
equivalence ratio, as shown in this example:

.. code-block:: yaml

    - kind: constant pressure
      pressure: 10.0
      temperature: 1000
      end-time: 10.0
      reactants:
          CH4: 1.0
          O2: 1.0
          N2: 3.76

For convenience, and to save significant runtime, pyMARS will automatically
reuse saved ignition data from a prior run, if the number of cases matches
that in the input file.

.. _conversion:

Conversion
==========

pyMARS provides a tool for converting between Chemkin and Cantera model
formats. (This is used implicity if a Chemkin model is given when running
pyMARS.) Generally this will be used to convert a Cantera reduced model
generated by pyMARS into a Chemkin-format model.

To convert a Cantera model into a Chemkin model, do

.. code-block:: bash

    pymars --convert -m model.cti

pyMARS also provides conversion from Chemkin to Cantera for convenience:

.. code-block:: bash

    pymars --convert -m model.inp --thermo thermo.dat

The ``_thermo`` option is not required if the thermodynamic data is contained
within the model file (i.e., with the ``THERMO ALL`` keyword). The transport
data can also be included in the resulting Cantera file with the
``--transport`` option.
