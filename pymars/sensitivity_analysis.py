import cantera as ct

from create_trimmed_model import trim
from simulation import Simulation
from drgep import get_rates
from readin_initial_conditions import readin_conditions
import numpy as np
import os
import helper


def create_limbo(reduced_model, ep_star, drgep_coeffs, safe):

	"""
	Creates a list of species in limbo for use during a sensitivity analysis.

	Parameters
	----------

	reduced_model: The model reduced by the previous reduction
	ep_star: Epsilon star value for the sensitivity analysis
	drgep_coeffs: The dictionary of direct interaction coefficients
	safe: species that are safe from being removed under any condition

	Returns
	-------

	A list of all species in limbo.
	
	"""
	
	limbo = []
	reduc_species = []
	species_objex = reduced_model.species()
	for sp in species_objex:
		reduc_species.append(sp.name)
	for sp in reduc_species:
		# All species that fit the condition of being in limbo are added to a list.
		if sp in drgep_coeffs and drgep_coeffs[sp] < ep_star and (not sp in limbo) and (not sp in safe):
			limbo.append(sp)
	return limbo

def get_limbo_dic(original_model, reduced_model, limbo, final_error, id_detailed, conditions_array):

	"""
	Creates a dictionary of all of the species in limbo and their errors for sensitivity analysis.

	Parameters
	----------

	original_model: The original version of the model being reduced
	reduced_model: The model produced by the previous reduction
	limbo: A list of the species in limbo
	final_error: Error percentage between the reduced and origanal models
	id_detailed: The ignition delays for each simulation of the original model
	conditions_array: An array holding the initial conditions for simulations
	
	Returns
	-------

	A dictionary with species error to be used for sensitivity anaylsis.

	"""

	dic = {}

	# For information on how this is set up, refer to run_sa function.
	og_excl = []
	keep = []
	og_sn = []
	new_sn = []

	# Append species names
	species_objex = reduced_model.species()
	for sp in species_objex:
		new_sn.append(sp.name)

	# Append species names
	species_objex = original_model.species()
	for sp in species_objex:
		og_sn.append(sp.name)

	# If its in original model, and new model then keep
	for sp in og_sn:
		if sp in new_sn:
			keep.append(sp)
	
	# If its not being kept, exclude (all that were removed in original reduction)
	for sp in original_model.species():
		if not (sp.name in keep):
			og_excl.append(sp.name)

	for sp in limbo: # For all species in limbo
		excluded = [sp]
		for p in og_excl:
			excluded.append(p) # Add that species to the list of exclusion.
		# Remove species from the model.
		new_sol_obs = trim(original_model,excluded,"sa_trim.cti")
		new_sol = new_sol_obs[1]

		# Simulated reduced solution
		new_sim = helper.setup_simulations(conditions_array,new_sol) # Create simulation objects for reduced model for all conditions
	
		try:	
			id_new = helper.simulate(new_sim) # Run simulations and process results
		except ct.CanteraError:
			limbo.remove(sp)
			id_new = 0
	
		error = (abs(id_new - id_detailed)/id_detailed)*100
		error = round(np.max(error), 2)
		print(sp + ": " + str(error))
		error = abs(error - final_error)
		dic[sp] = error # Add adjusted error to dictionary.
	return dic

def dic_lowest(dic):

	"""
	Gets the key with the lowest value in the dictionary.

	Parameters
	----------

	dic: The dictionary to get the lowest value out of

	Returns
	-------

	The key with the lowest value in the dictionary.

	"""

	lowest = 100000000
	s = "error"
	for sp in dic:
		if dic[sp] < lowest:
			lowest = dic[sp]
			s = sp
	return s


def run_sa(original_model, reduced_model, final_error, conditions_file, target, keepers, error_limit, limbo):
	"""Runs a sensitivity analysis on a resulting reduced model.
	
	Parameters
	----------
	original_model: The original version of the model being reduced
	reduced_model: The model produced by the previous reduction
	final_error: Error percentage between the reduced and origanal models
	conditions_file: The file holding the initial conditions for simulations
	target: The target species for the reduction
	keepers: A list of species that should be retained no matter what
	error_limit: The maximum allowed error between the reduced and original models
	limbo: A list of species to be considered for reduction by the sensativity analysis
	
	Returns
	-------
	The model after the sensitivity analysis has been preformed on it.

	"""

	if conditions_file:
		conditions_array = readin_conditions(str(conditions_file))
	elif not conditions_file:
		print("Conditions file not found")
		exit()

	# Turn conditions array into unran simulation objects for the original solution
	sim_array = helper.setup_simulations(conditions_array, original_model)
	id_detailed = helper.simulate(sim_array)  # Run simulations and process results

	rate_edge_data = get_rates(sim_array, original_model) # Get edge weight calculation data.
	
	if (id_detailed.all() == 0): # Ensure that ignition occured
		print("Original model did not ignite.  Check initial conditions.")
		exit()
	old = reduced_model

	while True:

		og_sn = []  # Original species names
		new_sn = []  # Species names in current reduced model
		keep = []  # Species retained from removals
		og_excl = []  # Species that will be excluded from the final model (Reduction will be preformed on original model)

		# Append species names
		species_objex = old.species()
		for sp in species_objex:
			new_sn.append(sp.name)

		# Append species names
		species_objex = original_model.species()
		for sp in species_objex:
			og_sn.append(sp.name)

		# If its in original model, and new model
		for sp in og_sn:
			if sp in new_sn:
				keep.append(sp)

		# If its not being kept, exclude (all that were removed in original reduction)
		for sp in original_model.species():
			if not (sp.name in keep):
				og_excl.append(sp.name)

		if len(limbo) == 0:
			return old

		print("In limbo:")
		print(limbo)

		# Calculate error for removing each limbo species.
		dic = get_limbo_dic(original_model,old,limbo,final_error, id_detailed,conditions_array)
		rm = dic_lowest(dic)  # Species that should be removed (Lowest error).
		exclude = [rm]
		limbo.remove(rm) # Remove species from limbo

		for sp in og_excl:  # Add to list of species that should be excluded from final model.
			exclude.append(sp)

		print()
		print("attempting to remove " + rm)
		
		# Remove exclusion list from original model
		new_sol_obs = trim(original_model,exclude,"sa_trim.cti")
		new_sol = new_sol_obs[1]

		# Simulated reduced solution
		new_sim = helper.setup_simulations(conditions_array,new_sol)  # Create simulation objects for reduced model for all conditions
		id_new = helper.simulate(new_sim)  # Run simulations and process results

		error = (abs(id_new - id_detailed)/id_detailed)*100
		error = round(np.max(error), 2)
		print("Error of: " + str(error))
		print()

		# If error is greater than allowed, previous reduced model was final reduction.
		if error > error_limit:
			print("Final Solution:")
			print(str(old.n_species) + " Species")
			return old

		else:  # If error is still within allowed limit, loop through again to further reduce.
			old = new_sol
