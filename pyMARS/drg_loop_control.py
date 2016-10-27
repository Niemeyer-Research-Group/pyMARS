from create_trimmed_model import trim
from convert_chemkin_file import convert
from soln2cti import write
from autoignition_loop_control import autoignition_loop_control
from get_rate_data import get_rates
from drg import make_graph
import os
from get_error import get_error
from graph_search import graph_search
from numpy import genfromtxt
import numpy as np

def drg_loop_control(solution_object, args):
    """ Controls repeated use of drg Function

        :param solution_object:
            Cantera solution object

        :param hdf5_file:
            string of hdf5 file containing production rates

        :param args:
            function arguments object

        """
    #get user input
    target_species = str(raw_input('Enter target starting species: '))

    #run detailed mechanism and retain initial conditions
    args.multiple_conditions = True
    detailed_result = autoignition_loop_control(solution_object, args)
    detailed_result.test.close()
    ignition_delay_detailed = np.array(detailed_result.tau_array)
    get_rates('mass_fractions.hdf5', solution_object)
    os.system('rm mass_fractions.hdf5')

    if args.threshold_values is None:
        try:
            threshold = float(raw_input('Enter threshold value: '))
        except ValueError:
            print 'try again'
            threshold = float(raw_input('Enter threshold value: '))
        drg = make_graph(solution_object, 'production_rates.hdf5', threshold)
        exclusion_list = graph_search(solution_object, drg, target_species)
        new_solution_objects = trim(solution_object, exclusion_list, args.data_file)

        reduced_result = autoignition_loop_control(new_solution_objects[1], args)
        reduced_result.test.close()
        ignition_delay_reduced = np.array(reduced_result.tau_array)
        error = (abs(ignition_delay_reduced-ignition_delay_detailed)/ignition_delay_detailed)*100
        print 'Maximum error: %s%' % max(error)
        #get_error()
        n_species_retained = len(new_solution_objects[1].species())
        print 'Number of species in reduced model: %s' %n_species_retained
        os.system('rm mass_fractions.hdf5')
    else:
        threshold_values = genfromtxt(args.threshold_values, delimiter=',')
        species_retained = []
        printout = ''
        print 'Threshold     Species in Mech      Error'
        for threshold in threshold_values:
            drg = make_graph(solution_object, 'production_rates.hdf5', threshold)
            exclusion_list = graph_search(solution_object, drg, target_species)
            new_solution_objects = trim(solution_object, exclusion_list, args.data_file)
            species_retained.append(len(new_solution_objects[1].species()))

            reduced_result = autoignition_loop_control(new_solution_objects[1], args)
            reduced_result.test.close()
            ignition_delay_reduced = np.array(reduced_result.tau_array)
            error = (abs(ignition_delay_reduced-ignition_delay_detailed)/ignition_delay_detailed)*100
            os.system('rm mass_fractions.hdf5')
            printout += str(threshold) + '  ' + str(len(new_solution_objects[1].species())) + '  '+  str(error) + '\n'

        print printout
    return new_solution_objects
