"""Writes a Cantera Solution object to a YAML file.

Delegates to Cantera's built-in write_yaml() method (available in Cantera 3.0+).
"""

import os


def write(solution, output_filename="", path=""):
    """Write a Cantera solution object to a YAML file.

    Parameters
    ----------
    solution : cantera.Solution
        Model to be written
    output_filename : str, optional
        Name of file to be written; if not provided, use ``solution.name``
    path : str, optional
        Path for writing file.

    Returns
    -------
    output_filename : str
        Name of output model file (.yaml)

    Examples
    --------
    >>> gas = cantera.Solution('gri30.yaml')
    >>> soln2yaml.write(gas, 'copy_gri30.yaml')
    copy_gri30.yaml

    """
    if output_filename:
        output_filename = os.path.join(path, output_filename)
    else:
        output_filename = os.path.join(path, f"{solution.name}.yaml")

    if os.path.isfile(output_filename):
        os.remove(output_filename)

    solution.write_yaml(output_filename)
    return output_filename
