import numpy as np

def get_range(times, temps, sdata, production_data):
    """Get sample range information from autoignition data

    :param times:

    :param temps:

    :param sdata:

    :param production_data:
    """

    times = np.array(times)
    temps = np.array(temps)
    time_temperature = np.vstack((times,temps)).T
    species_data = np.hstack((time_temperature, production_data))

    T = np.array(temps)
    dt = np.ones(len(times)-1)*(times[1]-times[0])
    dT = np.diff(T)
    derivative = dT/dt
    index = derivative.argmax()
    derivative_max = [times[index], T[index], index]
    tau = times[index]

    try:
        initial_point = [times[index-20], T[index-20], index-20]
    except IndexError:
        initial_point = [times[0], T[0], 0]
        error_string = ('not enough timesteps before ignition\n',
                        'timesteps before ignition: ' + str(index) +'\n',
                        'total timesteps: ' + str(len(T)))
    
        print error_string
    try:
        final_point = [times[index+20], T[index+20], index+20]
    except IndexError:
        final_point = [times[len(times)-1], T[len(T)-1], (len(T)-1)]
        error_string = ('not enough timesteps after ignition\n' +
                        'timesteps after ignition' + str(len((T))-index) +'\n' +
                        'total timesteps: ' + str(len(T)))
        print error_string

    class sample_data:
        def __init__(self, tau, index, times, temps, species_data,
                    production_data, derivative_max, initial_point, final_point):
            self.tau = tau
            self.index = index
            self.derivative_max = derivative_max
            self.initial_point = initial_point
            self.final_point = final_point

            self.times_total = times
            self.temps_total = temps
            self.species_data_total = species_data
            self.production_data_total = production_data

            self.times = times[ index-20:index+20 ]
            self.temps = temps[index-20:index+20 ]
            self.species_data = species_data[index-20:index+20, :]
            self.production_data = production_data[index-20:index+20, :]
    return sample_data(tau, index, times, temps, species_data, production_data,
                        derivative_max, initial_point, final_point)
