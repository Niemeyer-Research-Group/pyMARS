import cantera as ct
import numpy as np

def check_error(original_mechanism, new_mechanism):
    original_solution = ct.Solution(original_mechanism)
    reduced_solution = ct.Solution(new_mechanism)
    mole_fractions = raw_input('Enter mole fractions (ex.CH4:1, O2:2, N2:7.52 for Gri30 Stoich) :  ') # (ex.CH4:1, O2:2, N2:7.52 for Gri30 Stoich)and 100 for TPY where Y is mass fractions
    initial_temperature = float(raw_input('Enter Solution Temperature in (K):'))
    original_solution.TPX = initial_temperature, ct.one_atm, mole_fractions #1001.0
    reduced_solution.TPX = initial_temperature, ct.one_atm, mole_fractions #1001.0

    def advance(solution_object):
        reactor = ct.Reactor(solution_object)
        sim = ct.ReactorNet([reactor])
        current_time = 0.0
        end_time = 5.0e-3
        time_array = []
        temperature_array = []
        data_index = 0
        while current_time < end_time:
            data_index += 1
            current_time = sim.step(end_time)
            time_array.append(current_time)
            temperature_array.append(reactor.T)
        time_array = np.array(time_array)
        temperature_array = np.array(temperature_array)
        temperature_profile = np.vstack((time_array, temperature_array)).T
        return [temperature_profile, time_array, temperature_array]

    def get_ignition_delay(time_array, temperature_array):
        dt = np.ones(len(time_array)-1)*(time_array[1]-time_array[0])
        dT = np.diff(temperature_array)
        deriv = dT/dt
        tau = float(time_array[deriv.argmax()])
        return tau

    original_simulation = advance(original_solution)
    reduced_simulation = advance(reduced_solution)
    original_delay = get_ignition_delay(original_simulation[1], original_simulation[2])
    reduced_delay = get_ignition_delay(reduced_simulation[1], reduced_simulation[2])

    print ('Original ign delay: %s [s]') %original_delay
    print ('Reduced ign delay: %s [s]') %reduced_delay
    error = ((original_delay-reduced_delay)/original_delay)*100.0
    print ('Error: %s %') %error
