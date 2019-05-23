from ._version import __version__  # noqa: F401
import logging
from cantera import suppress_thermo_warnings

# Avoid long warnings from Cantera about thermodynamic polynomials
suppress_thermo_warnings() 

handlers = [logging.FileHandler('info.log'), logging.StreamHandler()]
logging.basicConfig(level=logging.INFO, format='%(message)s', handlers=handlers)
