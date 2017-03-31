from create_trimmed_model import trim
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

        Parameters
        ----------
        solution_object : obj
            Cantera solution object
        args : obj
            function arguments object

        Returns
        -------
        new_solution_objects : obj
            Cantera solution object with skeletal mechanism
        """
    #get user input
    target_species = str(raw_input('Enter target starting species: '))

    #run detailed mechanism and retain initial conditions
    args.multiple_conditions = True
    detailed_result = autoignition_loop_control(solution_object, args)
    detailed_result.test.close()
    ignition_delay_detailed = np.array(detailed_result.tau_array)
    get_rates('mass_fractions.hdf5', solution_object)
    #print 'triggered'
    #print args.threshold_values

    if args.threshold_values is None:
        try:
            threshold = float(raw_input('Enter threshold value: '))
        except ValueError:
            print 'try again'
            threshold = float(raw_input('Enter threshold value: '))

        #run DRG and create new reduced solution
        drg = make_graph(solution_object, 'production_rates.hdf5', threshold)
        exclusion_list = graph_search(solution_object, drg, target_species)
        new_solution_objects = trim(solution_object, exclusion_list, args.data_file)

        #simulate reduced solution
        reduced_result = autoignition_loop_control(new_solution_objects[1], args)
        reduced_result.test.close()
        ignition_delay_reduced = np.array(reduced_result.tau_array)
        error = (abs(ignition_delay_reduced-ignition_delay_detailed)/ignition_delay_detailed)*100
        print 'Error index: %s' %error
        #get_error()
        n_species_retained = len(new_solution_objects[1].species())
        print 'Number of species in reduced model: %s' %n_species_retained
        os.system('rm mass_fractions.hdf5')
    else:
        threshold_values = genfromtxt(args.threshold_values, delimiter=',')
        species_retained = []
        printout = ''
        print 'Threshold     Species in Mech      Error'
        try:
            os.system('rm mass_fractions.hdf5')
        except Exception:
            pass
        if type(threshold_values) is list:
            for threshold in threshold_values:
                try:
                    os.system('rm mass_fractions.hdf5')
                except Exception:
                    pass
                #run DRG and create new reduced solution
                drg = make_graph(solution_object, 'production_rates.hdf5', threshold)
                exclusion_list = graph_search(solution_object, drg, target_species)
                new_solution_objects = trim(solution_object, exclusion_list, args.data_file)
                species_retained.append(len(new_solution_objects[1].species()))

                #simulated reduced solution
                reduced_result = autoignition_loop_control(new_solution_objects[1], args)
                reduced_result.test.close()
                ignition_delay_reduced = np.array(reduced_result.tau_array)
                error = (abs(ignition_delay_reduced-ignition_delay_detailed)/ignition_delay_detailed)*100
                printout += str(threshold) + '  ' + str(len(new_solution_objects[1].species())) + '  '+  str(np.max(error)) + '\n'
        else:
            if os.path.exists('mass_fractions.hdf5'):
                os.system('rm mass_fractions.hdf5')
            #run DRG and create new reduced solution
            drg = make_graph(solution_object, 'production_rates.hdf5', threshold_values)
            exclusion_list = graph_search(solution_object, drg, target_species)
            new_solution_objects = trim(solution_object, exclusion_list, args.data_file)
            species_retained.append(len(new_solution_objects[1].species()))

            #simulated reduced solution
            reduced_result = autoignition_loop_control(new_solution_objects[1], args)
            reduced_result.test.close()
            ignition_delay_reduced = np.array(reduced_result.tau_array)
            error = (abs(ignition_delay_reduced-ignition_delay_detailed)/ignition_delay_detailed)*100
            printout += str(threshold_values) + '  ' + str(len(new_solution_objects[1].species())) + '  '+  str(np.max(error)) + '\n'
        print printout
        print 'Detailed soln ign delay %0.05f'  %ignition_delay_detailed
        print 'Reduced soln ign delay %0.05f'   %ignition_delay_reduced
    return new_solution_objects
