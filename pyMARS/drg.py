import networkx as nx
import cantera as ct


def make_graph(data_file):
    """ Use the Direct Relation Graph (DRG) method to choose species to
    eliminate from reaction mechanism.

    Parameters:
    -------------
    data_file: Mechanism file in the form of .cti

    Returns:
    -------------
    exclusion_list: List of species to trim from mechanism
    """

    
