""" Tests the get rates method for the DRG algorithm """

import sys
import os
import math 

import cantera as ct
import pytest

from .. import drg
from .. import helper
from ..readin_initial_conditions import readin_conditions

def relative_location(file):
	file_path = os.path.join(file)
	return pkg_resources.resource_filename(__name__, file_path)

@pytest.mark.xfail
def testArtificial():
	
	# Load model
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	# Set up simulation data from givin initial coniditions
	conditions = relative_location("example_input_artificial.txt") # Locate conditions file
	conditions_array = readin_conditions(conditions) # Load conditions
	sim_array = helper.setup_simulations(conditions_array, solution_object) # Set up simulation
	helper.simulate(sim_array) # Run autoignition simulation

	rate_edge_data = drg.get_rates_drg(sim_array, solution_object) # Run unit

	print(rate_edge_data)

	# Pull out timestep one denomenator and numerator dicts
	ic_one = rate_edge_data[list(rate_edge_data.keys())[0]]
	tstep_one = ic_one[list(ic_one.keys())[0]]
	denoms = tstep_one[0]
	numers = tstep_one[1]

	# Expected values for denomenators
	expected_denoms = {}
	expected_denoms["H2O"] = 1.9573216e-13
	expected_denoms["H2"] = .00025854374
	expected_denoms["O2"] = 9.7866081e-14
	expected_denoms["H"] = .00051708749

	assert math.isclose(expected_denoms["H2O"],denoms["H2O"],abs_tol=1.0e-17)
	assert math.isclose(expected_denoms["H2"],denoms["H2"],abs_tol=1.0e-10)
	assert math.isclose(expected_denoms["O2"],denoms["O2"],abs_tol=1.0e-18)
	assert math.isclose(expected_denoms["H"],denoms["H"],abs_tol=1.0e-10)

	expected_numers = {}
	expected_numers["H2O_H2"] = 1.9573216e-13
	expected_numers["H2O_O2"] = 1.9573216e-13
	expected_numers["H2_O2"] = 1.9573216e-13
	expected_numers["H2_H2O"] = 1.9573216e-13
	expected_numers["O2_H2"] = 9.7866081e-14
	expected_numers["O2_H2O"] = 9.7866081e-14
	expected_numers["H2_H"] = .00025854374
	expected_numers["H_H2"] = .00051708749
	
	assert math.isclose(expected_numers["H2O_H2"],numers["H2O_H2"],abs_tol=1.0e-17)
	assert math.isclose(expected_numers["H2O_O2"],numers["H2O_O2"],abs_tol=1.0e-17)
	assert math.isclose(expected_numers["H2_O2"],numers["H2_O2"],abs_tol=1.0e-17)
	assert math.isclose(expected_numers["H2_H2O"],numers["H2_H2O"],abs_tol=1.0e-17)
	assert math.isclose(expected_numers["O2_H2"],numers["O2_H2"],abs_tol=1.0e-18)
	assert math.isclose(expected_numers["O2_H2O"],numers["O2_H2O"],abs_tol=1.0e-18)
	assert math.isclose(expected_numers["H2_H"],numers["H2_H"],abs_tol=1.0e-18)
	assert math.isclose(expected_numers["H_H2"],numers["H_H2"],abs_tol=1.0e-18)
