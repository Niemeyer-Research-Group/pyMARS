"""writes a solution object to a chemkin inp file

currently only works for Elementary, Falloff and ThreeBody Reactions
"""

import os
from textwrap import fill

import cantera as ct

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0

# Conversion from 1 debye to coulomb-meters
DEBEYE_CONVERSION = 3.33564e-30


def build_arrhenius(rate, reaction_order, reaction_type):
    """Builds Arrhenius coefficient string based on reaction type.

    Parameters
    ----------
    rate : cantera.ArrheniusRate
        Arrhenius-form reaction rate coefficient
    reaction_order : int or float
        Order of reaction (sum of reactant stoichiometric coefficients)
    reaction_type : str
        Reaction type string (e.g. 'Arrhenius', 'three-body-Arrhenius',
        'pressure-dependent-Arrhenius')

    Returns
    -------
    str
        String with Arrhenius coefficients

    """
    if reaction_type in ("Arrhenius", "pressure-dependent-Arrhenius"):
        pre_exponential_factor = rate.pre_exponential_factor * 1e3 ** (
            reaction_order - 1
        )

    elif reaction_type == "three-body-Arrhenius":
        pre_exponential_factor = rate.pre_exponential_factor * 1e3**reaction_order

    elif reaction_type.startswith("falloff-") or reaction_type.startswith(
        "chemically-activated"
    ):
        raise ValueError(
            "Function does not support falloff or chemically activated reactions"
        )
    else:
        raise NotImplementedError("Reaction type not supported: ", reaction_type)

    arrhenius = [
        f"{pre_exponential_factor:.4e}",
        f"{rate.temperature_exponent:.3f}",
        f"{(rate.activation_energy / CALORIES_CONSTANT):.2f}",
    ]
    return "  ".join(arrhenius)


def build_falloff_arrhenius(rate, reaction_order, reaction_type, pressure_limit):
    """Builds Arrhenius coefficient strings for falloff and chemically-activated reactions.

    Parameters
    ----------
    rate : cantera.ArrheniusRate
        Arrhenius-form reaction rate coefficient
    reaction_order : int or float
        Order of reaction (sum of reactant stoichiometric coefficients)
    reaction_type : str
        Reaction type string (e.g. 'falloff-Troe', 'falloff-Lindemann', 'chemically-activated')
    pressure_limit : {'high', 'low'}
        string designating pressure limit

    Returns
    -------
    str
        Arrhenius coefficient string

    """
    assert pressure_limit in ["low", "high"], "Pressure range needs to be high or low"

    if reaction_type.startswith("falloff-"):
        if pressure_limit == "low":
            pre_exponential_factor = rate.pre_exponential_factor * 1e3 ** (
                reaction_order
            )
        elif pressure_limit == "high":
            pre_exponential_factor = rate.pre_exponential_factor * 1e3 ** (
                reaction_order - 1
            )

    elif reaction_type.startswith("chemically-activated"):
        if pressure_limit == "low":
            pre_exponential_factor = rate.pre_exponential_factor * 1e3 ** (
                reaction_order - 1
            )
        elif pressure_limit == "high":
            pre_exponential_factor = rate.pre_exponential_factor * 1e3 ** (
                reaction_order - 2
            )
    else:
        raise ValueError("Reaction type not supported: ", reaction_type)

    arrhenius = [
        f"{pre_exponential_factor:.4e}",
        f"{rate.temperature_exponent:.3f}",
        f"{(rate.activation_energy / CALORIES_CONSTANT):.3e}",
    ]
    return "  ".join(arrhenius)


def build_falloff(parameters, falloff_function):
    """Creates falloff reaction Troe parameter string

    Parameters
    ----------
    parameters : numpy.ndarray
        Array of falloff parameters; length varies based on ``falloff_function``
    falloff_function : {'Troe', 'SRI'}
        Type of falloff function

    Returns
    -------
    falloff_string : str
        String of falloff parameters

    """
    if falloff_function == "Troe":
        falloff_string = (
            "TROE / " + f"{parameters[0]}  {parameters[1]}  "
            f"{parameters[2]}  {parameters[3]} /\n"
        )
    elif falloff_function == "SRI":
        falloff_string = (
            "SRI / "
            + f"{parameters[0]}  {parameters[1]}  "
            + f"{parameters[2]}  {parameters[3]}  {parameters[4]} /\n"
        )
    else:
        raise NotImplementedError(f"Falloff function not supported: {falloff_function}")

    return falloff_string


def write_thermo_data(species_list, filename="generated_thermo.dat"):
    """Writes thermodynamic data to Chemkin-format file.

    Parameters
    ----------
    species_list : list of cantera.Species
        List of species objects
    filename : str, optional
        Filename for new Chemkin thermodynamic database file

    """
    with open(filename, "w") as the_file:

        the_file.write("THERMO\n" + "   300.000  1000.000  5000.000\n")

        # write data for each species in the Solution object
        for species in species_list:

            composition_string = "".join(
                [f"{s:2}{int(v):>3}" for s, v in species.composition.items()]
            )

            species_string = (
                f"{species.name:<18}"
                + 6 * " "
                + f"{composition_string:<20}"
                + "G"
                + f"{species.thermo.min_temp:10.3f}"
                + f"{species.thermo.max_temp:10.3f}"
                + f"{species.thermo.coeffs[0]:8.2f}"
                + 6 * " "
                + "1\n"
            )

            species_string += (
                "".join([f"{c:15.8e}" for c in species.thermo.coeffs[1:6]])
                + "    "
                + "2\n"
            )

            species_string += (
                "".join([f"{c:15.8e}" for c in species.thermo.coeffs[6:8]])
                + "".join([f"{c:15.8e}" for c in species.thermo.coeffs[8:11]])
                + "    "
                + "3\n"
            )

            species_string += (
                "".join([f"{c:15.8e}" for c in species.thermo.coeffs[11:15]])
                + 19 * " "
                + "4\n"
            )

            the_file.write(species_string)

        the_file.write("END\n")


def write_transport_data(species_list, filename="generated_transport.dat"):
    """Writes transport data to Chemkin-format file.

    Parameters
    ----------
    species_list : list of cantera.Species
        List of species objects
    filename : str, optional
        Filename for new Chemkin transport database file

    """
    geometry = {"atom": "0", "linear": "1", "nonlinear": "2"}

    with open(filename, "w") as the_file:

        # write data for each species in the Solution object
        for species in species_list:
            species_string = (
                f"{species.name:<16}"
                + f"{geometry[species.transport.geometry]:>4}"
                + f"{(species.transport.well_depth / ct.boltzmann):>10.3f}"
                + f"{(species.transport.diameter * 1e10):>10.3f}"
                + f"{(species.transport.dipole / DEBEYE_CONVERSION):>10.3f}"
                + f"{(species.transport.polarizability * 1e30):>10.3f}"
                + f"{species.transport.rotational_relaxation:>10.3f}"
                + "\n"
            )

            the_file.write(species_string)


def write(
    solution, output_filename="", path="", skip_thermo=False, skip_transport=False
):
    """Writes Cantera solution object to Chemkin-format file.

    Parameters
    ----------
    solution : cantera.Solution
        Model to be written
    output_filename : str, optional
        Name of file to be written; if not provided, use ``solution.name``
    path : str, optional
        Path for writing file.
    skip_thermo : bool, optional
        Flag to skip writing thermo data in separate file
    skip_transport : bool, optional
        Flag to skip writing transport data in separate file

    Returns
    -------
    output_file_name : str
        Name of output model file (.inp)

    Examples
    --------
    >>> gas = cantera.Solution('gri30.yaml')
    >>> soln2ck.write(gas)
    reduced_gri30.inp

    """
    if output_filename:
        output_filename = os.path.join(path, output_filename)
    else:
        output_filename = os.path.join(path, f"{solution.name}.inp")

    if os.path.isfile(output_filename):
        os.remove(output_filename)

    with open(output_filename, "w") as the_file:

        # Write title block to file
        the_file.write(
            f"!Chemkin file converted from solution object: {solution.name}\n\n"
        )

        # write species and element lists to file
        element_names = "  ".join(solution.element_names)
        the_file.write("ELEMENTS\n" + f"{element_names}\n" + "END\n\n")
        species_names = fill(
            "  ".join(solution.species_names),
            width=60,
            break_long_words=False,
            break_on_hyphens=False,
        )
        the_file.write("SPECIES\n" + f"{species_names}\n" "END\n\n")

        # Write reactions to file
        the_file.write("REACTIONS\n")

        for reaction in solution.reactions():

            rxn_type = reaction.reaction_type
            reaction_string = f"{reaction.equation:<51}"

            # The Arrhenius parameters that follow the equation string on the main line
            # depend on the type of reaction.
            if rxn_type in ("Arrhenius", "three-body-Arrhenius"):
                arrhenius = build_arrhenius(
                    reaction.rate, sum(reaction.reactants.values()), rxn_type
                )

            elif rxn_type.startswith("falloff-"):
                # high-pressure limit is on the main reaction line
                arrhenius = build_falloff_arrhenius(
                    reaction.rate.high_rate,
                    sum(reaction.reactants.values()),
                    rxn_type,
                    "high",
                )

            elif rxn_type.startswith("chemically-activated"):
                # low-pressure limit is on the main reaction line
                arrhenius = build_falloff_arrhenius(
                    reaction.rate.low_rate,
                    sum(reaction.reactants.values()),
                    rxn_type,
                    "low",
                )

            elif rxn_type == "Chebyshev":
                arrhenius = "1.0e0  0.0  0.0"

            elif rxn_type == "pressure-dependent-Arrhenius":
                arrhenius = build_arrhenius(
                    reaction.rate.rates[0][1],
                    sum(reaction.reactants.values()),
                    "pressure-dependent-Arrhenius",
                )
            else:
                raise NotImplementedError(f"Unsupported reaction type: {rxn_type}")

            reaction_string += arrhenius + "\n"

            # print third-body efficiencies when a third body is present
            if reaction.third_body is not None:
                tb = reaction.third_body
                reduced_efficiencies = {
                    s: v
                    for s, v in tb.efficiencies.items()
                    if s in solution.species_names
                }
                efficiencies_str = "  ".join(
                    [f"{s}/{v}/" for s, v in reduced_efficiencies.items()]
                )
                if efficiencies_str:
                    reaction_string += efficiencies_str + "\n"

            # write auxiliary information
            if rxn_type.startswith("falloff-"):
                arrhenius = build_falloff_arrhenius(
                    reaction.rate.low_rate,
                    sum(reaction.reactants.values()),
                    rxn_type,
                    "low",
                )
                reaction_string += f"LOW / {arrhenius} /\n"

                if reaction.rate.falloff_coeffs.size > 0:
                    falloff_str = build_falloff(
                        reaction.rate.falloff_coeffs, reaction.rate.sub_type
                    )
                    reaction_string += falloff_str

            elif rxn_type.startswith("chemically-activated"):
                arrhenius = build_falloff_arrhenius(
                    reaction.rate.high_rate,
                    sum(reaction.reactants.values()),
                    rxn_type,
                    "high",
                )
                reaction_string += f"HIGH / {arrhenius} /\n"

                if reaction.rate.falloff_coeffs.size > 0:
                    falloff_str = build_falloff(
                        reaction.rate.falloff_coeffs, reaction.rate.sub_type
                    )
                    reaction_string += falloff_str

            elif rxn_type == "pressure-dependent-Arrhenius":
                for pressure, rate in reaction.rate.rates:
                    pressure_str = f"{pressure / ct.one_atm}"
                    arrhenius = build_arrhenius(
                        rate,
                        sum(reaction.reactants.values()),
                        "pressure-dependent-Arrhenius",
                    )
                    reaction_string += f"PLOG / {pressure_str} {arrhenius} /\n"

            elif rxn_type == "Chebyshev":
                cheb = reaction.rate
                reaction_string += (
                    f"TCHEB / {cheb.temperature_range[0]}  {cheb.temperature_range[1]} /\n"
                    + f"PCHEB / {cheb.pressure_range[0] / ct.one_atm}  "
                    f"{cheb.pressure_range[1] / ct.one_atm} /\n"
                    + f"CHEB / {cheb.n_temperature}  {cheb.n_pressure} /\n"
                )
                for coeffs in cheb.data:
                    coeffs_row = " ".join([f"{c:.6e}" for c in coeffs])
                    reaction_string += f"CHEB / {coeffs_row} /\n"

            if reaction.duplicate:
                reaction_string += "DUPLICATE\n"

            the_file.write(reaction_string)

        the_file.write("END")

    basename = os.path.splitext(output_filename)[0]
    outputs = [output_filename]

    # write thermo data
    if not skip_thermo:
        write_thermo_data(solution.species(), basename + "_thermo.dat")
        outputs.append(basename + "_thermo.dat")

    if not skip_transport and all(sp.transport for sp in solution.species()):
        write_transport_data(solution.species(), basename + "_transport.dat")
        outputs.append(basename + "_transport.dat")

    return outputs
