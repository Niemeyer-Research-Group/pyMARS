import os, sys, argparse
import cantera as ct
os.environ['Cantera_Data'] =os.getcwd()
import soln2ck
import soln2cti
from get_rate_data_drg import get_rates_drg
from drg_loop_control import drg_loop_control
from autoignition_loop_control import autoignition_loop_control
import numpy as np

def run_drg(args, solution_object):           
	"""
    Function to run the drgep method for model reduction

    Parameters
    ----------
    solution_object: ct._Solutuion object
	An object that represents the solution to be trimmed. 
    args: args object
	An object that contains the following:
    file:
        Input mechanism file (ex. file='gri30.cti')
    species:
        Species to eliminate (ex. species='H, OH')
    thermo:
        Thermo data file if Chemkin format (ex. thermo= 'thermo.dat')
    transport:
        Transport data file if Chemkin format
    plot:
        plot ignition curve (ex. plot='y')
    points:
        print ignition point and sample range (ex. points='y')
    writecsv:
        write data to csv (ex. writecsv='y')
    writehdf5:
        write data to hdf5 (ex. writehdf5='y')
    run_drg:
        Run DRG model reduction
    run_drgep:
	Run drgep model.
    error:
	Maximum ammount of error allowed. 
    keepers: list of strings
	The string names of the species that should be kept in the model no matter what.
    targets: list of strings
	The string names of the species that should be used as target species.

	"""

	if len(args.target) == 0: #If the target species are not specified, puke and die.
		print("Please specify a target species.")
		exit()
	done = [] #Singleton to hold wether or not any more species can be cut from the simulation.  
	done.append(False)
	threshold = .1 #Starting threshold value
	threshold_i = .1
	n = 1
	error = [10.0] #Singleton to hold the error value of the previously ran simulation.
	
	try:
		os.system('rm mass_fractions.hdf5')
	except Exception:
		pass
	
	detailed_result = autoignition_loop_control(solution_object, args) #The simulation needs to be ran to make the mass_fractions file which has the info to calucalte edge weights I think?
	detailed_result.test.close()
	ignition_delay_detailed = np.array(detailed_result.tau_array)
	rate_edge_data = get_rates_drg('mass_fractions.hdf5', solution_object) #Get edge weight calculation data. 
	
	print("Testing for starting threshold value")
	drg_loop_control(solution_object, args, error, threshold, done, rate_edge_data) #Trim the solution at that threshold and find the error.
	while error[0] != 0: #While the error for trimming with that threshold value is greater than allowed.
		threshold = threshold / 10 #Reduce the starting threshold value and try again.
		threshold_i = threshold_i / 10
		n = n + 1
		drg_loop_control(solution_object, args, error, threshold, done, rate_edge_data)
	
	print("Starting with a threshold value of " + str(threshold))
	sol_new = solution_object
	past = 0 #An integer representing the error introduced in the past simulation.  
	done[0] = False

	while not done[0] and error[0] < args.error: #Run the simulation until nothing else can be cut. 
		sol_new = drg_loop_control( solution_object, args, error, threshold, done, rate_edge_data) #Trim at this threshold value and calculate error.
		if args.error > error[0]: #If a new max species cut without exceeding what is allowed is reached, save that threshold.
                        max_t = threshold
		#if (past == error[0]): #If error wasn't increased, increase the threshold at a higher rate. 
		#	threshold = threshold + (threshold_i * 4)
		past = error[0] 
		#if (threshold >= .01):
                #        threshold_i = .01
		threshold = threshold + threshold_i
		threshold = round(threshold, n)

	print("\nGreatest result: ")
	sol_new = drg_loop_control( solution_object, args, error, max_t, done, rate_edge_data)
	drgep_trimmed_file = soln2cti.write(sol_new) #Write the solution object with the greatest error that isn't over the allowed ammount.
	
	try:
		os.system('rm production_rates.hdf5')
	except Exception:
		pass
