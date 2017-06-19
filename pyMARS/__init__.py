from create_trimmed_model import trim
from pyMARS import readin
from convert_chemkin_file import convert
from write_to_cti import write
from autoignition_module import run_sim
import os
os.environ['Cantera_Data'] =os.getcwd()
