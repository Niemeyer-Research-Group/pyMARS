import cantera as ct

def test(original_file, new_file):
    """Test written cti file against original cti file.

    Arguments
    -------------
    Original cti file
    Newly written cti file

    Returns
    -----------
    List of file properties that do not match


    """
    original = ct.Solution(original_file)
    new = ct.Solution(new_file)
    comparison_list=[]
    for i, species in enumerate(new.species_names):
        if species in original.species_names:
            comparison_list.append(species)


    def test_species_def():
        num =0
        #make sure comparing same species
        for i, name1 in enumerate(new.species_names):

            #start comparison with same species
            new_species= new.species(i)
            for j, name2 in enumerate(original.species_names):
                if name1.upper().lower() == name2.upper().lower():
                    original_species=original.species(j)
                    num += 1
            if original_species.name.upper().lower() != new_species.name.upper().lower():
                print (j, original_species, i, new_species,)


            assert original_species.composition == new_species.composition

            assert original_species.thermo.coeffs.all() == new_species.thermo.coeffs.all()

            try:
                assert original_species.transport.geometry == new_species.transport.geometry
            except AttributeError:
                pass

            try:
                assert original_species.transport.diameter == new_species.transport.diameter
            except AttributeError:
                pass

            try:
                assert original_species.transport.well_depth == new_species.transport.well_depth
            except AttributeError:
                pass

            try:
                assert original_species.transport.polarizability == new_species.transport.polarizability
            except AttributeError:
                pass
            try:
                assert original_species.transport.rotational_relaxation == new_species.transport.rotational_relaxation
            except AttributeError:
                pass

            try:
                assert original_species.transport.dipole == new_species.transport.dipole
            except AttributeError:
                pass

        print ('done with testing species definition into \n\n\n')
    def test_reactions():
        c=4184.0
        num = 0
        print 'Any errors shown below:\n'
        for k, name1 in enumerate(new.reaction_equations()):
            num += 1
            new_reaction=new.reaction(k)
            new_equation_type = type(new_reaction).__name__
            for l, name2 in enumerate(original.reaction_equations()):
                if original.reaction(l).equation == new_reaction.equation:
                    original_reaction=original.reaction(l)
                    original_equation_type=type(original_reaction).__name__

            assert new_equation_type == original_equation_type
            try:
                if new_reaction.rate.pre_exponential_factor != original_reaction.rate.pre_exponential_factor:
                    #if new_reaction.rate.pre_exponential_factor/ original_reaction.rate.pre_exponential_factor > .004:
                    print (k, (new_reaction.rate.pre_exponential_factor/ original_reaction.rate.pre_exponential_factor), new_reaction.reaction_type, new_reaction.rate.temperature_exponent, (new_reaction.rate.activation_energy/c) , new_reaction )
            except AttributeError:
                pass
            #assert new_reaction.efficiencies == original_reaction.efficiencies
        print ('\ndone with testing equation info ')
    test_species_def()
    test_reactions()


test('gri301.cti', 'pym_gri30.cti')
