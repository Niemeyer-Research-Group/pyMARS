import cantera as ct

def test(original_file, new_file):
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

            assert original_species.transport.geometry == new_species.transport.geometry

            assert original_species.transport.diameter == new_species.transport.diameter

            assert original_species.transport.well_depth == new_species.transport.well_depth

            assert original_species.transport.polarizability == new_species.transport.polarizability

            assert original_species.transport.rotational_relaxation == new_species.transport.rotational_relaxation

            assert original_species.transport.dipole == new_species.transport.dipole

    def test_reactions():
        num = 0
        for k, name1 in enumerate(new.reaction_equations()):
            num += 1
            new_reaction=new.reaction(k)
            new_equation_type = type(new_reaction).__name__
            for l, name2 in enumerate(original.reaction_equations()):
                if original.reaction(l).equation == new_reaction.equation:
                    original_reaction=original.reaction(l)
                    original_equation_type=type(original_reaction).__name__

            assert new_equation_type == original_equation_type

            if new_equation_type == 'ThreeBodyReaction':
                if new_reaction.efficiencies != original_reaction.efficiencies:
                    print (k, new_reaction, l, original_reaction)
            if new_equation_type == 'ElementaryReaction':
                if new_reaction.rate.pre_exponential_factor != original_reaction.rate.pre_exponential_factor:
                        print (k, new_reaction, l, original_reaction)
            #assert new_reaction.efficiencies == original_reaction.efficiencies

    test_species_def()
    test_reactions()

test('gri30.cti', 'trimmed_gri30_mod.cti')
