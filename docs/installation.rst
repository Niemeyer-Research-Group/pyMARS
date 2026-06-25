.. _sec_installation:

==================
Installation Guide
==================

pyMARS requires Python 3.10+ and is available on Linux, macOS, and Windows.

Via pip (recommended)
---------------------

Install from PyPI::

    pip install nrg-pymars

.. note::

   On PyPI, pyMARS is distributed under the name ``nrg-pymars``, because the
   ``pymars`` name belongs to an unrelated, active project. The import package
   and command-line tool remain ``pymars`` (i.e., ``pip install nrg-pymars``,
   then ``import pymars`` or run ``pymars``).

To install the latest development version directly from GitHub::

    pip install git+https://github.com/Niemeyer-Research-Group/pyMARS.git

Or clone the repository and install from the local copy::

    git clone https://github.com/Niemeyer-Research-Group/pyMARS/
    cd pyMARS
    pip install .

Via conda
---------

Install pyMARS from the ``niemeyer-research-group`` channel, with
``conda-forge`` providing the remaining dependencies (including Cantera)::

    conda install -c niemeyer-research-group -c conda-forge pymars

Development
-----------

To install an editable (development) copy, clone the repository and use
``pip``::

    git clone https://github.com/Niemeyer-Research-Group/pyMARS/
    cd pyMARS
    pip install -e .

To uninstall::

    pip uninstall pymars
