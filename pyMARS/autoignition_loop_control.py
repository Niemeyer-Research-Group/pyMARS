import cantera as ct
from autoignition_module import run_sim

def autoignition_loop_control(solution_object, args):
    """Controls autoignition module

    :param solution_object:
        Cantera solution object
    :param args:
        Arguments from terminal
    """
