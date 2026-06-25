"""Module with miscellanous tools."""

import os
from operator import attrgetter
import logging

import numpy as np
import cantera as ct
from cantera import ck2yaml

from . import soln2ck


def compare_models(model1, model2):
    """Checks whether two Cantera models are equivalent.

    Parameters
    ----------
    model1 : cantera.Solution
        First Cantera model object
    model2 : cantera.Solution
        Second Cantera model object

    Returns
    -------
    bool
        ``True`` if models match, ``False`` if models disagree

    """
    # at minimum, numbers of species and reactions should be the same
    if model1.n_species != model2.n_species:
        return False
    if model1.n_reactions != model2.n_reactions:
        return False

    for sp1, sp2 in zip(
        sorted(model1.species(), key=attrgetter("name")),
        sorted(model2.species(), key=attrgetter("name")),
    ):
        if sp1.name != sp2.name:
            return False

        if sp1.composition != sp2.composition:
            return False

        for T_ref in (500.0, 1000.0, 2000.0):
            if not np.isclose(sp1.thermo.cp(T_ref), sp2.thermo.cp(T_ref), rtol=1e-5):
                return False
            if not np.isclose(sp1.thermo.h(T_ref), sp2.thermo.h(T_ref), rtol=1e-5):
                return False
            if not np.isclose(sp1.thermo.s(T_ref), sp2.thermo.s(T_ref), rtol=1e-5):
                return False

        if hasattr(sp1, "transport") or hasattr(sp2, "transport"):
            if hasattr(sp1, "transport") and hasattr(sp2, "transport"):
                # iterate over transport parameters
                params = [
                    a
                    for a in dir(sp1.transport)
                    if not a.startswith("__")
                    and not callable(getattr(sp1.transport, a))
                ]
                for attr in params:
                    if getattr(sp1.transport, attr) != getattr(
                        sp2.transport, attr, 0.0
                    ):
                        return False
            else:
                return False

    for rxn1, rxn2 in zip(
        sorted(model1.reactions(), key=attrgetter("equation")),
        sorted(model2.reactions(), key=attrgetter("equation")),
    ):
        if rxn1.reaction_type != rxn2.reaction_type:
            return False

        if rxn1.reactants != rxn2.reactants:
            return False
        if rxn1.products != rxn2.products:
            return False

        if rxn1.duplicate != rxn2.duplicate:
            return False

        # Check primary rate object (present on all reactions in Cantera 3)
        if hasattr(rxn1, "rate") and hasattr(rxn2, "rate"):
            rate1, rate2 = rxn1.rate, rxn2.rate
            if type(rate1) is not type(rate2):
                return False

            # Falloff / chemically-activated: compare low-pressure and high-pressure rates
            if hasattr(rate1, "low_rate") and hasattr(rate2, "low_rate"):
                lo1, lo2 = rate1.low_rate, rate2.low_rate
                hi1, hi2 = rate1.high_rate, rate2.high_rate
                for attr in (
                    "pre_exponential_factor",
                    "temperature_exponent",
                    "activation_energy",
                ):
                    if not np.isclose(getattr(lo1, attr, 0.0), getattr(lo2, attr, 0.0)):
                        return False
                    if not np.isclose(getattr(hi1, attr, 0.0), getattr(hi2, attr, 0.0)):
                        return False
                # falloff parameters (Troe, SRI, …)
                if hasattr(rate1, "falloff_coeffs") and hasattr(
                    rate2, "falloff_coeffs"
                ):
                    if not np.allclose(rate1.falloff_coeffs, rate2.falloff_coeffs):
                        return False

            # Plog: compare pressure-dependent rate table
            elif hasattr(rate1, "rates") and hasattr(rate2, "rates"):
                table1 = rate1.rates
                table2 = rate2.rates
                if len(table1) != len(table2):
                    return False
                for (p1, r1), (p2, r2) in zip(
                    sorted(table1, key=lambda x: x[0]),
                    sorted(table2, key=lambda x: x[0]),
                ):
                    if not np.isclose(p1, p2):
                        return False
                    for attr in (
                        "pre_exponential_factor",
                        "temperature_exponent",
                        "activation_energy",
                    ):
                        if not np.isclose(
                            getattr(r1, attr, 0.0), getattr(r2, attr, 0.0)
                        ):
                            return False

            # Chebyshev: compare polynomial coefficients and ranges
            elif hasattr(rate1, "data") and hasattr(rate2, "data"):
                if not np.allclose(rate1.pressure_range, rate2.pressure_range):
                    return False
                if not np.allclose(rate1.temperature_range, rate2.temperature_range):
                    return False
                if not np.allclose(rate1.data, rate2.data):
                    return False

            # Simple Arrhenius (Arrhenius / three-body-Arrhenius)
            else:
                for attr in (
                    "pre_exponential_factor",
                    "temperature_exponent",
                    "activation_energy",
                ):
                    if not np.isclose(
                        getattr(rate1, attr, 0.0), getattr(rate2, attr, 0.0)
                    ):
                        return False

        # Third-body efficiencies (Cantera 3: accessed via reaction.third_body)
        tb1 = getattr(rxn1, "third_body", None)
        tb2 = getattr(rxn2, "third_body", None)
        if (tb1 is None) != (tb2 is None):
            return False
        if tb1 is not None and tb2 is not None:
            if tb1.default_efficiency != tb2.default_efficiency:
                return False
            if tb1.efficiencies != tb2.efficiencies:
                return False

    return True


def convert(model_file, thermo_file=None, transport_file=None, path=""):
    """Function to convert between Cantera and Chemkin model formats.

    Parameters
    ----------
    model_file : str
        Input model file (Cantera .yaml or Chemkin)
    thermo_file : str, optional
        Chemkin thermodynamic properties file
    transport_file : str, optional
        Chemkin transport data file
    path : str, optional
        Path for writing file

    Returns
    -------
    str or list
        Path to converted file, or list of files (for Chemkin)

    Example
    -------
    >>> convert('gri30.inp')
    gri30.yaml

    >>> convert('gri30.yaml')
    [gri30.inp, gri30_thermo.dat, gri30_transport.dat]

    """
    basename = os.path.splitext(os.path.basename(model_file))[0]
    extension = os.path.splitext(os.path.basename(model_file))[1]

    # Chemkin files can have multiple extensions, so easier to check if Cantera
    if extension in (".yaml", ".yml", ".cti"):
        # Convert from Cantera to Chemkin format.
        logging.info("Converter detected Cantera input model: " + model_file)
        logging.info("Converting to Chemkin format.")

        solution = ct.Solution(model_file)
        converted_files = soln2ck.write(solution, basename + ".inp", path=path)
        return converted_files
    else:
        # Convert from Chemkin to Cantera YAML format.
        logging.info("Converter detected Chemkin input model: " + model_file)
        logging.info("Converting to Cantera format.")

        converted_file = os.path.join(path, basename + ".yaml")

        # calls ck2yaml based on given files
        args = [f"--input={model_file}"]
        if thermo_file:
            args.append(f"--thermo={thermo_file}")
        if transport_file:
            args.append(f"--transport={transport_file}")
        args.append(f"--output={converted_file}")

        # generally Chemkin files have issues (redundant species, etc.) that require this argument
        args.append("--permissive")

        ck2yaml.main(args)
        return converted_file
