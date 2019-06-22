.. _sec_installation:

==================
Installation Guide
==================

pyMARS is available for Python 3.6+ on Linux, macOS, and Windows via
``conda`` and ``pip``. pyMARS depends on ``cantera``, so you should
first install that via the ``cantera`` channel, then install
pyMARS via the ``niemeyer-research-group`` channel::

    conda install -c cantera cantera
    conda install -c niemeyer-research-group pymars

Note that you might need to include the ``conda-forge`` channel
by editing your conda configuration::

    conda config --append channels conda-forge
    conda install -c niemeyer-research-group pymars

You can also install using ``pip`` by downloading the source code and changing
into that directory::

    git clone https://github.com/Niemeyer-Research-Group/pyMARS/
    cd pyMARS
    pip install pymars

Unfortunately, it is not currently possible to install using ``pip``
from PyPI, because (1) there is no Cantera package on PyPI and
(2) someone else claimed the ``pymars`` name.

Development
-----------

pyMARS can be installed from source by cloning the git repository
and changing into that directory::

    git clone https://github.com/Niemeyer-Research-Group/pyMARS/
    cd pyMARS

Then run::

    conda develop .

if you're using ``conda`` (you may need to install ``conda-build`` first).
To uninstall, run::

    conda develop . --uninstall

Note that this doesn't install the standalone converter scripts. With
``pip``, installing is done by::

    pip install -e .

To uninstall with ``pip``::

    pip uninstall pymars

``pip`` does install the standalone scripts.
