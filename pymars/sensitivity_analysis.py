import numpy as np
import cantera as ct

from . import helper
from .create_trimmed_model import trim
from .simulation import Simulation
from .drgep import get_rates
from .readin_initial_conditions import readin_conditions


def create_limbo(reduced_model, ep_star, drgep_coeffs, safe):
	"""Creates a list of species in limbo for use during sensitivity analysis.

	Parameters
	----------

	reduced_model : Cantera.Solution
		The model reduced by the previous reduction
	ep_star	: float
		Epsilon star (upper threshold) value for sensitivity analysis
	drgep_coeffs : dict
		The dictionary of direct interaction coefficients
	safe : list of Cantera.Species
		Species that are protected from removal under any condition

	Returns
	-------
	limbo : list of Cantera.Species
		List of all species in limbo
	
	"""
	
	limbo = []
	reduc_species = []
	species_objex = reduced_model.species()
	for sp in species_objex:
		reduc_species.append(sp.name)
	for sp in reduc_species:
		# All species that fit the condition of being in limbo are added to a list.
		if (sp in drgep_coeffs and drgep_coeffs[sp] < ep_star and 
			(not sp in limbo) and (not sp in safe)
			):
			limbo.append(sp)
	return limbo


def get_limbo_dic(original_model, reduced_model, limbo, final_error, 
				  id_detailed, conditions_array
				  ):
	"""Creates dictionary of all species in limbo and their errors for sensitivity analysis.

	Parameters
	----------
	original_model : cantera.Solution
		Original model being reduced
	reduced_model : cantera.Solution
		Model produced by the previous reduction stage
	limbo : list of Cantera.Species
		List of species in limbo
	final_error : float
		Error percentage between the reduced and origanal models
	id_detailed : numpy.array
		The ignition delays for each simulation of the original model
	conditions_array : list of Condition
		Array holding the initial conditions for simulations
	
	Returns
	-------
	Dictionary with species error to be used for sensitivity anaylsis.

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
		new_sol = trim(original_model, excluded, "sa_trim.cti")

		# Simulated reduced solution
		new_sim = helper.setup_simulations(conditions_array, new_sol) # Create simulation objects for reduced model for all conditions
	
		try:	
			id_new = helper.simulate(new_sim) # Run simulations and process results
		except ct.CanteraError:
			limbo.remove(sp)
			id_new = 0
	
		error = 100 * abs(id_new - id_detailed) / id_detailed
		error = round(np.max(error), 2)
		print(sp + ": " + str(error))
		error = abs(error - final_error)
		dic[sp] = error # Add adjusted error to dictionary.
	return dic


def dic_lowest(dic):
	"""Gets the key with the lowest value in the dictionary.

	Parameters
	----------
	dic : dict
		Dictionary to get the lowest value out of

	Returns
	-------
	The key with the lowest value in the dictionary

	"""

	lowest = 100000000
	s = "error"
	for sp in dic:
		if dic[sp] < lowest:
			lowest = dic[sp]
			s = sp
	return s


def run_sa(model_file, final_error, conditions_file, target, keepers, error_limit, limbo):
	"""Runs a sensitivity analysis on a resulting reduced model.
	
	Parameters
	----------
	original_model : Cantera.Solution
		The original version of the model being reduced
	reduced_model : Cantera.Solution
		The model produced by the previous reduction
	final_error : float
		Error percentage between the reduced and origanal models
	conditions_file : str
		Filename of file holding the initial conditions for simulations
	target : 
		The target species for the reduction
	keepers : 
		A list of species that should be retained no matter what
	error_limit : float
		The maximum allowed error between the reduced and original models
	limbo : list of Cantera.Species
		List of species to be considered for reduction by the sensitivity analysis
	
	Returns
	-------
	The model after the sensitivity analysis has been performed on it

	"""
	original_model = ct.Solution(model_file)

	if conditions_file:
		conditions_array = readin_conditions(str(conditions_file))
	elif not conditions_file:
		print("Conditions file not found")
		exit()

	# Turn conditions array into unran simulation objects for the original solution
	sim_array = helper.setup_simulations(conditions_array, original_model)
	id_detailed = helper.simulate(sim_array)  # Run simulations and process results
	
	if (id_detailed.all() == 0): # Ensure that ignition occured
		print("Original model did not ignite.  Check initial conditions.")
		exit()
	old = reduced_model

	# TODO: replace this potentially infinite while loop
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
		new_sol = trim(original_model, exclude, "sa_trim.cti")

		# Simulated reduced solution
		# Create simulation objects for reduced model for all conditions
		new_sim = helper.setup_simulations(conditions_array, new_sol)
		# Run simulations and process results
		id_new = helper.simulate(new_sim)

		error = 100 * abs(id_new - id_detailed) / id_detailed
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
