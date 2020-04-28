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

.. code-block:: none

     -h, --help:
        show this help message and exit
     -i, --input:
        YAML file with reduction inputs
     --path:
        Path to directory for writing files
     --num_threads:
        Number of CPU cores to use for running simulations in parallel.
        If no number, then use available number of cores minus 1.
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
first create a YAML file (let's call it ``reduction_input.yaml``:

.. code-block:: yaml

    model: gri30.cti
    targets:
      - CH4
      - O2
    retained-species:
      - N2
    method: DRGEP
    error: 10.0
    sensitivity-analysis: False
    autoignition-conditions:
      - kind: constant volume
        pressure: 1.0
        temperature: 1000.0
        fuel:
          CH4: 1.0
        oxidizer:
          O2: 1.0
          N2: 3.76
        equivalence-ratio: 1.0

      - kind: constant volume
        pressure: 1.0
        temperature: 1200.0
        fuel:
          CH4: 1.0
        oxidizer:
          O2: 1.0
          N2: 3.76
        equivalence-ratio: 0.5

Then, execute this command:

.. code-block:: bash

    pymars --input reduction_input.yaml

You can also download an annotated version of the
:download:`example input file <example_input_file.yaml>`; refer to the
:ref:`inputfile` section for a detailed explanation.

**Parallelization:** If you have a large number of initial conditions, the
reduction may be sped up by parallelizing the associated simulations over
multiple CPU cores. You can do this by adding ``--num_threads N``, where ``N``
is the desired number of cores. If you just specify ``--num_threads``, then
pyMARS will use the available number of cores minus one.

(pyMARS does not currently support distributed memory parallelization, meaning
across multiple nodes that do not share the same memory.)

**Sensitivity analysis:** To perform sensitivity analysis following DRGEP,
change the ``sensitivity-analysis`` key to ``True`` in the input file,
and choose the type of sensitivity analysis with the ``sensitivity-type`` field
(either ``initial`` or ``greedy``). Let's use the "initial" method for now, since
it is less computationally expensive.

Generally it is a good idea to specify an upper
threshold value of 0.2 to ensure important species are not evaluated; do this
by adding the line ``upper-threshold: 0.2``. (In the :ref:`sec_theory`
guide, this variable is referred to as :math:`\epsilon^*`.)

Your new input file (called ``drgepsa_input.yaml``) would then look like:

.. code-block:: yaml

    model: gri30.cti
    targets:
      - CH4
      - O2
    retained-species:
      - N2
    method: DRGEP
    error: 10.0
    sensitivity-analysis: True
    sensitivity-type: initial
    upper-threshold: 0.2
    autoignition-conditions:
      - kind: constant volume
        pressure: 1.0
        temperature: 1000.0
        fuel:
          CH4: 1.0
        oxidizer:
          O2: 1.0
          N2: 3.76
        equivalence-ratio: 1.0

      - kind: constant volume
        pressure: 1.0
        temperature: 1200.0
        fuel:
          CH4: 1.0
        oxidizer:
          O2: 1.0
          N2: 3.76
        equivalence-ratio: 0.5

Then, the command for performing a reduction using DRGEPSA with
parallelized simalations would be:

.. code-block:: bash

   pymars --input drgepsa_input.yaml --num_threads


.. _inputfile:

Reduction input file
====================

You control the model reduction process in pyMARS through a YAML input file,
indicated with the ``--input`` or ``-i`` command-line argument. Keys include:

- ``model:`` filename of kinetic model being reduced (Chemkin or Cantera)
- ``phase-name:`` Optional name of phase in Cantera CTI file to be reduced
- ``targets:`` List of one or more target species; required for DRG, DRGEP,
  and PFA methods
- ``retained-species:`` Optional list of one or more species to never remove.
- ``method``: Reduction method; one of ``DRG``, ``DRGEP``, or ``PFA``
- ``error``: Maximum error limit of reduced model, given as percentage
  (e.g., ``10.0`` for 10%).
- ``sensitivity-analysis``: Specify ``True`` to perform sensitivity analysis,
  either alone or following a method given by ``method``
- ``sensitivity-type``: Type of sensitivity analysis, either
  ``initial`` or ``greedy``
- ``upper-threshold``: Upper threshold value for species to be considered for
  sensitivity analysis; only used when following one of the graph-based
  reduction methods
- ``autoignition-conditions``: List of initial conditions for autoignition
  simulations, described in more detail next

Species given in ``targets``, ``retained-species``, or
``fuel``/``oxidizer``/``reactants`` must be present in the model specified in
``model``, spelling must match exactly (including case).

**Autoignition parameters:** pyMARS currently uses autoignition simulations to
sample thermochemical data for the reduction and to calculate ignition delays
for measuring error of candidate reduced models. Initial conditions need to be
provided for performing these simulations, in the ``autoignition-conditions``
field of the input file.

These initial conditions are given as a list, with these required keys:

- ``kind``: Type of homogeneous autoignition simulation; either
  ``constant volume`` or ``constant pressure``
- ``pressure``: initial pressure, given in atm
- ``temperature``: initial temperature, given in K

The initial reactant mixture can be given using either an equivalence
ratio with separate fuel and oxidizer specifications, or as list of
reactants.

To specify the mixture using an equivalence ratio, you must give lists
of species in the fuel and oxidizer, with the mole fraction/number of
the species in each (these will be automatically normalized):

.. code-block:: yaml

    fuel:
      CH4: 1.0
    oxidizer:
      O2: 1.0
      N2: 3.76
    equivalence-ratio: 1.0

To specify the mixture using a list of reactants, just give the number of
moles of each species in the initial mixture:

.. code-block:: yaml

    reactants:
      CH4: 1.0
      O2: 2.0
      N2: 7.52

When giving the composition as a list of reactants, you can also specify the
mass fraction of the mixture using ``composition-type: mass``:

.. code-block:: yaml

    reactants:
      CH4: 0.05518632
      O2: 0.22014867
      N2: 0.724665
    composition-type: mass

**Note:** By default, pyMARS automatically integrates each autoignition case
to steady state, or to a maximum of 10,000 integration steps. This can be
bypassed by specifying either a different number of maximum steps with
``max-steps:`` or a maximum integration end time with ``end-time:``
(in seconds). During initial sampling, pyMARS will raise an error if it
does not detect autoignition, based on reaching the initial temperature +
400 K.

For convenience, and to save significant runtime when reducing the same
model with different parameters, pyMARS will automatically
reuse saved ignition data from a prior run. It semi-intelligently checks
if the number of cases matches that in the input file, but to be safe
output files should be cleaned between applications.


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
