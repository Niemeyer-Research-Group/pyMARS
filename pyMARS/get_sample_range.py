import numpy as np

def get_range(times, temps, sdata, production_data):
    """Get sample range information from autoignition data

    Parameters
    ----------
    times : numpy matrix
        Times from autoignition
    temps : numpy matrix
        Autoigntion temperatures
    sdata : numpy matrix
        Species mass fraction data at a point in time
    production_data : numpy matrix
        Species total production data  at a point in time

    Returns
    -------
    sample_data : obj
        .tau
        .index
        .derivative_max
        .initial_point
        .final_point

        .times_total
        .temps_total
        .species_data_total
        .production_data_total

        .times
        .temps
        .species_data
        .production_data
    """
    times = np.array(times)
    temps = np.array(temps)
    time_temperature = np.vstack((times,temps)).T
    species_data = np.hstack((time_temperature, sdata)) #was production_data

    T = np.array(temps)
    dt = np.ones(len(times)-1)*(times[1]-times[0])
    dT = np.diff(T)
    derivative = dT/dt
    index = 0
    for i in range (0, len(temps)):
        if index == 0 and temps[i] >= temps[0] + 400:
            index = i
    derivative_max = [times[index], T[index], index]
    tau = times[index]

    try:
        initial_point = [times[index-20], T[index-20], index-20]
    except IndexError:
        initial_point = [times[0], T[0], 0]
        error_string = ('not enough timesteps before ignition\n',
                        'timesteps before ignition: ' + str(index) +'\n',
                        'total timesteps: ' + str(len(T)))

        #print error_string
    try:
        final_point = [times[index+20], T[index+20], index+20]
    except IndexError:
        final_point = [times[len(times)-1], T[len(T)-1], (len(T)-1)]
        error_string = ('not enough timesteps after ignition\n' +
                        'timesteps after ignition: ' + str(len((T))-index) +'\n' +
                        'total timesteps: ' + str(len(T)))
        #print error_string
    delta = temps[len(temps)-1] - temps[0]
    time_trim = []
    temp_trim = []
    spec_trim = []
    prod_trim = []
    ind_trim = []
    time_trim.append(times[0])
    temp_trim.append(temps[0])
    spec_trim.append(species_data[0])
    prod_trim.append(production_data[0])
    ind_trim.append(0)
    i = .05
    while i < 1: #store values every time the tempreture reaches five percent more of its total change for 20 total data points
        j = 0
        while temps[j] < temps[0] + (delta * i):
            j = j + 1
        time_trim.append(times[j])
        temp_trim.append(temps[j])
        spec_trim.append(species_data[j])
        prod_trim.append(production_data[j])
        ind_trim.append(j)
        i = i + .05
    #time_trim.append(times[len(times)-1])
    #temp_trim.append(temps[len(times)-1])
    #spec_trim.append(species_data[len(times)-1])
    #prod_trim.append(production_data[len(times)-1])
    #ind_trim.append(len(times)-1)
    initial_point = [times[0], T[0], 0]
    f = len(times)-1
    final_point = [times[f], T[f], f]
    class sample_data:
        def __init__(self, tau, index, times, temps, species_data,
                    production_data, derivative_max, initial_point, final_point,ti,te,sp,pr,i):
            self.tau = tau
            self.index = index
            self.derivative_max = derivative_max
            self.initial_point = initial_point
            self.final_point = final_point

            self.times_total = times
            self.temps_total = temps
            self.species_data_total = species_data
            self.production_data_total = production_data

            self.times = ti
            self.temps = te
            self.species_data = sp
            self.production_data =pr
            self.index = i
    return sample_data(tau, index, times, temps, species_data, production_data,
                        derivative_max, initial_point, final_point,time_trim,temp_trim,spec_trim,prod_trim,ind_trim)
