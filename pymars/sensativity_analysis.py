from create_trimmed_model import trim
from simulation import Simulation
from drgep import make_dic_drgep
from drgep import get_rates
from readin_initial_conditions import readin_conditions
import numpy as np
import os
import helper

#################
## Function: create_limbo
## Description: Creates a list of speices in limbo for use during a sensativity analysis.
## Input: The reduced version of the original model, the ep star value for sensativity analysis, and the drgep coefficient dictionary from the original model.
## Parameters: A Cantera gas solution, an integer, and a dictionary with species names as keys and overall interaction coeffiecents from the original model as values.
## Post-conditions: A list that may or may not be empty.
## Return: A list of all species in limbo.
#################

def create_limbo(reduced_model, ep_star, drgep_coeffs,safe):
	limbo = []
	reduc_species = []
	species_objex = reduced_model.species()
	for sp in species_objex:
		reduc_species.append(sp.name)
	for sp in reduc_species:
		if sp in drgep_coeffs and drgep_coeffs[sp] < ep_star and (not sp in limbo) and (not sp in safe): #All species that fit the condition of being in limbo are added to a list.
			limbo.append(sp)
	return limbo

#################
## Function: get_limbo_dic
## Description: Creates a dictionary of all of the species in limbo and their errors for sensativity analysis.
## Input: The original model, the reduced version of the original model, the array of speceis in limbo, the final error from the reduced model, and the command line arguements.
## Parameters: Two Cantera gas solutions, a list of strings, an integer, and a user defined args class.
## Post-conditions: A dictionary is keyed by species that stores the error value when that species is removed minus the final error.
## Return: A dictionary with species error to be used for sensativity anaylsis.
#################

def get_limbo_dic(original_model,reduced_model,limbo,final_error,args,id_detailed,conditions_array):
	dic = {}

	og_excl = [] #For information on how this is set up, refer to run_sa function.
	keep = []
	og_sn = []
	new_sn = []

	species_objex = reduced_model.species()
	for sp in species_objex:
		new_sn.append(sp.name)

	species_objex = original_model.species()
	for sp in species_objex:
		og_sn.append(sp.name)

	for sp in og_sn:
		if sp in new_sn:
			keep.append(sp)

	for sp in original_model.species():
		if not (sp.name in keep):
			og_excl.append(sp.name)

	for sp in limbo: #For all species in limbo
		excluded = [sp]
		for p in og_excl:
			excluded.append(p) #Add that species to the list of exclusion.
		new_sol_obs = trim(original_model,excluded,"sa_trim.cti") #Remove species from the model.
		new_sol = new_sol_obs[1]

		#simulated reduced solution
		new_sim = helper.setup_simulations(conditions_array,new_sol) #Create simulation objects for reduced model for all conditions
		id_new = helper.simulate(new_sim) #Run simulations and process results
		error = (abs(id_new - id_detailed)/id_detailed)*100
		error = round(np.max(error), 2)
		print(sp + ": " + str(error))
		error = abs(error - final_error)
		dic[sp] = error #Add adjusted error to dictionary.
	return dic

#################
## Function: dic_lowest
## Description: Gets the key with the lowest value in the dictionary.
## Input: The dictionary to be operated on.
## Parameters: A dictionary that has at least one key value pair where the value is an integer.
## Post-conditions: A string that represents the key with the lowest value in the dictionary.
## Return: The key with the lowest value in the dictionary.
#################

def dic_lowest(dic):
	lowest = 100000000
	s = "error"
	for sp in dic:
		if dic[sp] < lowest:
			lowest = dic[sp]
			s = sp
	return s

#################
## Function: run_sa
## Description: Runs a sensativity analysis on a resulting reduced model.
## Input: The original and reduced models, the ep star value, the final error, and the command like arguments.
## Parameters: Two Cantera gas solutions, two integers, and an args object.
## Post-conditions: A Cantera gas soltion.
## Return: The model after the sensativity analysis has been preformed on it.
#################

def run_sa(original_model,reduced_model,ep_star,final_error,args):
	print(final_error)

	if args.conditions_file:
		conditions_array = readin_conditions(str(args.conditions_file))
	elif not args.conditions_file:
		print("Conditions file not found")
		exit()

	sim_array = helper.setup_simulations(conditions_array,original_model) #Turn conditions array into unran simulation objects for the original solution
	id_detailed = helper.simulate(sim_array) #Run simulations and process results

	rate_edge_data = get_rates(sim_array,original_model) #Get edge weight calculation data.
	drgep_coeffs = make_dic_drgep(original_model, rate_edge_data, args.target) #Make a dictionary of overall interaction coefficients.
	if (id_detailed.all() == 0): #Ensure that ignition occured
		print("Original model did not ignite.  Check initial conditions.")
		exit()
	old = reduced_model

	while True:

		og_sn = [] #Original species names
		new_sn = [] #Species names in current reduced model
		keep = [] #Species retained from removals
		og_excl = [] #Species that will be excluded from the final model (Reduction will be preformed on original model)

		species_objex = old.species()
		for sp in species_objex:
			new_sn.append(sp.name)

		species_objex = original_model.species()
		for sp in species_objex:
			og_sn.append(sp.name)

		for sp in og_sn:
			if sp in new_sn:
				keep.append(sp)

		for sp in original_model.species():
			if not (sp.name in keep):
				og_excl.append(sp.name)

		limbo = create_limbo(old,ep_star,drgep_coeffs,args.keepers) #Find all the species in limbo.

		if len(limbo) == 0:
			return old

		print("In limbo:")
		print(limbo)

		dic = get_limbo_dic(original_model,old,limbo,final_error,args,id_detailed,conditions_array) #Calculate error for removing each limbo species.
		rm = dic_lowest(dic) #Species that should be removed (Lowest error).
		exclude = [rm]

		for sp in og_excl: #Add to list of species that should be excluded from final model.
			exclude.append(sp)

		print()
		print("attempting to remove " + rm)
		new_sol_obs = trim(original_model,exclude,"sa_trim.cti") #Remove exclusion list from original model
		new_sol = new_sol_obs[1]

		#simulated reduced solution
		new_sim = helper.setup_simulations(conditions_array,new_sol) #Create simulation objects for reduced model for all conditions
		id_new = helper.simulate(new_sim) #Run simulations and process results

		error = (abs(id_new - id_detailed)/id_detailed)*100
		error = round(np.max(error), 2)
		print("Error of: " + str(error))
		print()

		if error > args.error: #If error is greater than allowed, previous reduced model was final reduction.
			print("Final Solution:")
			print(str(old.n_species) + " Species")
			return old

		else: #If error is still within allowed limit, loop through again to further reduce.
			old = new_sol
