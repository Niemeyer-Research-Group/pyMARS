.. _sec_installation:

==================
Installation Guide
==================

pyMARS requires Python 3.10+ and is available on Linux, macOS, and Windows.

Via conda (recommended)
-----------------------

Install pyMARS from the ``niemeyer-research-group`` channel, with
``conda-forge`` providing the remaining dependencies (including Cantera)::

    conda install -c niemeyer-research-group -c conda-forge pymars

Via pip
-------

Install directly from GitHub::

    pip install git+https://github.com/Niemeyer-Research-Group/pyMARS.git

Or clone the repository and install from the local copy::

    git clone https://github.com/Niemeyer-Research-Group/pyMARS/
    cd pyMARS
    pip install .

.. note::

   The ``pymars`` name on PyPI is currently held by an unrelated project.
   A name reclaim request is in progress per
   `PEP 541 <https://peps.python.org/pep-0541/>`_.
   We will publish to PyPI once the name is recovered.

Development
-----------

To install an editable (development) copy, clone the repository and use
``pip``::

    git clone https://github.com/Niemeyer-Research-Group/pyMARS/
    cd pyMARS
    pip install -e .

To uninstall::

    pip uninstall pymars
