from autoignition_loop_control import autoignition_loop_control

def get_error(solution_object, args):
    """Get error from ignition delay




    """
    args.initial_sim = False
    if args.tau_array:
        detailed_tau_array = args.tau_array
        error_array = []
        for index, initial_condition in enumerate(args.initial_temperature_array):
            tau_detailed = detailed_tau_array[index]
            args.Temp = initial_condition
            try:
                solution_result = autoignition_loop_control(solution_object, args)
                tau_skeletal = solution_result.tau
                error_array.append(float((abs((tau_detailed-tau_skeletal)/tau_detailed))*100.0))
            except Exception:
                tau_skeletal = 0.0
                error_array.append(0.0)

        print error_array
        print 'Error: %s%%' %"{0:.2f}"''.format(max(error_array))
    else:
        tau_detailed = args.tau
        solution_result = autoignition_loop_control(solution_object, args)
        tau_skeletal = solution_result.tau
        error = float((abs((tau_detailed-tau_skeletal)/tau_detailed))*100.0)
        print 'Error: %s%%' %"{0:.2f}"''.format(max(error))
