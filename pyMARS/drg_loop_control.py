from create_trimmed_model import trim
from convert_chemkin_file import convert
from soln2cti import write
from autoignition_module import run_sim
from get_rate_data import get_rates
from drg import make_graph

def loop_control(solution_object, args):
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



    #run first sim and retain initial conditions
    args.initiial_sim = True
    sim1_result = run_sim(solution_object, args)
    tau1 = sim1_result.tau
    args.frac = sim1_result.frac
    args.Temp = sim1_result.Temp
    get_rates('mass_fractions.hdf5', solution_object)
    args.initial_sim = False


    if args.iterate is True:
        #set initial 0 cases
        error = 0.0
        threshold = 0.00
        loop_number = 0
        error_limit = float(raw_input('Acceptable Error Limit (%): '))
        while error < error_limit:
            loop_number += 1
            threshold +=.01
            drg_exclusion_list = make_graph(solution_object, 'production_rates.hdf5', threshold, target_species)
            new_solution_objects = trim(solution_object, drg_exclusion_list, args.data_file)
            #run second sim
            sim2_result = run_sim(new_solution_objects[1], args)
            tau2 = sim2_result.tau
            error = float((abs((tau1-tau2)/tau1))*100.0)
        threshold -= .01
        drg_exclusion_list = make_graph(solution_object, 'production_rates.hdf5', threshold, target_species)
        new_solution_objects = trim(solution_object, drg_exclusion_list, args.data_file)
        #run second sim
        sim2_result = run_sim(new_solution_objects[1], args)
        #compare error
        tau2 = sim2_result.tau
        error = float((abs((tau1-tau2)/tau1))*100.0)
        print 'Loop number: ' + str(loop_number) + '|| error: ' + str(error) + ' %'
        print 'Number of loops: %s' %loop_number
        print 'Final max threshold value: %s' %threshold
        print 'Error: %s ' %error
    else:
        threshold = float(raw_input('Enter threshold value: '))
    return new_solution_objects
