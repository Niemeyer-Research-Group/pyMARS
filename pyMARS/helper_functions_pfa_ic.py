#####################
# File: helper_functions_pfa_ic.py
# Description:  Functions to be used when calculating interaction coefficients between species in the PFA reducation method.  
# Author: Phillip Mestas
# Date: 1/30/18
# Input: Model being reduced 
# Output: Various dictionaries with calculations along the way to getting the direct interaction coefficients between species.  
######################

import cantera as ct

######################
# Function: get_PA
# Description: Gets the PA (and CA) values of all species in a given solution.  
# Parameters: Must be a valid model stored in the cantera format with corresponding production rates.
# Input: A solution class from the Cantera library and an array of production rates.
# Output: Two dictionaries both keyed by species name.  One for PA, one for CA.  
# Returns: PA and CA dictionaries.  
#######################

def get_PA(new_solution, new_reaction_production_rates):
	PA = {} #Dictionary that will hold the PA values for each species.
	CA = {} #Dictionary that will hold the CA values for each species.
	
	s_names = new_solution.species_names
	for species in s_names: #For species in the solutuion
		PA[species] = 0
		CA[species] = 0
		
		for i, reac in enumerate(new_solution.reactions()): #For all reactions
			reac_prod_rate = float(new_reaction_production_rates[i]) #Set up values
			all_species = reac.reactants 
			all_species.update(reac.products)
		
			if reac_prod_rate != 0:
				if species in all_species: #If the species is in a productive reaction, sum up omega times v in the appropriate dictionary.
					add = float(reac_prod_rate * all_species[species])
					if add > 0:
						PA[species] += abs(add)
					else:
						CA[species] += abs(add)
	
	return PA,CA


######################
# Function: get_PAB
# Description: Gets the PAB (and CAB) values of all species in a given solution.  
# Parameters: Must be a valid model stored in the cantera format with corresponding production rates.
# Input: A solution class from the Cantera library and an array of production rates.
# Output: Two dictionaries both keyed by species name.  One for PAB, one for CAB.  
# Returns: PAB and CAB dictionaries.  
#######################

def get_PAB(new_solution, new_reaction_production_rates):
	PAB = {} #Set up dictionaries
	CAB = {}

	s_names = new_solution.species_names
	for species_a in s_names: #For every pair of species A and B in the solution
		for species_b in s_names:
			if species_a != species_b:
				full_name = species_a + "_" + species_b
				PAB[full_name] = 0
				CAB[full_name] = 0

				for i, reac in enumerate(new_solution.reactions()): #For all reactions
					reac_prod_rate = float(new_reaction_production_rates[i]) #Set up values
					all_species = reac.reactants 
					all_species.update(reac.products)

					if reac_prod_rate != 0: #If both species exsist in the reaction, add the calculated value to the correct dictionary.  
						if species_a in all_species:
							if species_b in all_species: 
								add = float(reac_prod_rate * all_species[species_a])
								if add > 0:
									PAB[full_name] += abs(add)
								else:
									CAB[full_name] += abs(add)

	return PAB, CAB

######################
# Function: get_rAB_1
# Description: Gets the rAB_p1 (and rAB_c1) values of all species in a given solution.  
# Parameters: Must be a valid model stored in the cantera format with corresponding production rates and 4 dictionaries.
# Input: A solution class from the Cantera library and all of the PA,CA,PAB, and CAB values.
# Output: Two dictionaries both keyed by species name.  One for rAB_p1, one for rAB_c1.  
# Returns: rAB_p1 and rAB_c1 dictionaries.  
#######################

def get_rAB_1(new_solution,PA,CA,PAB,CAB):
	rAB_p1 = {} #Set up dictionaries
	rAB_c1 = {}

	s_names = new_solution.species_names
	for species_a in s_names: #For all pairs of species
		for species_b in s_names:
			if species_a != species_b:
				full_name = species_a + "_" + species_b #Set up
				rAB_p1[full_name] = 0
				rAB_c1[full_name] = 0
	
				top_p = PAB[full_name] #Get top
				top_c = CAB[full_name]
	
				if (PA[species_a] > CA[species_a]): #Get bot
					bot = PA[species_a]
				else:
					bot = CA[species_a]
				
				if (bot != 0): #Calculate 
					rAB_p1[full_name] = top_p/bot
					rAB_c1[full_name] = top_c/bot

	return rAB_p1, rAB_c1


######################
# Function: get_rAB_2
# Description: Gets the rAB_p2 (and rAB_c2) values of all species in a given solution.  
# Parameters: Must be a valid model stored in the cantera format with corresponding production rates and 2 dictionaries.
# Input: A solution class from the Cantera library and all of the rAB_p1 and rAB_c1 values.
# Output: Two dictionaries both keyed by species name.  One for rAB_p2, one for rAB_c2.  
# Returns: rAB_p2 and rAB_c2 dictionaries.  
#######################

def get_rAB_2(new_solution,rAB_p1,rAB_c1):
	rAB_p2 = {} #Set up dictionaries
	rAB_c2 = {}

	s_names = new_solution.species_names
	for species_a in s_names: #For all pairs of species
		for species_b in s_names:
			if species_a != species_b:
				full_name = species_a + "_" + species_b #Set up
				rAB_p2[full_name] = 0
				rAB_c2[full_name] = 0
	
				for species_m in s_names: #Look through all possible middle step species
					if (species_m != species_a and species_m != species_b):
						am_name = species_a + "_" + species_m
						mb_name = species_m + "_" + species_b
						
						add_p = rAB_p1[am_name] * rAB_p1[mb_name] #Get what to add for species_m
						add_c = rAB_c1[am_name] * rAB_c1[mb_name]
		
						rAB_p2[full_name] += add_p #Add that value
						rAB_c2[full_name] += add_c

	return rAB_p2,rAB_c2
