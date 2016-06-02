

import sys
import numpy as np
import cantera as ct


#read in solution file, and make constant pressure adiabatic reactor
solution1 = ct.Solution('gri30.cti')  #trimmed_h2_v1b_mech.cti
solution1.TPX = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'
r1 = ct.IdealGasConstPressureReactor(solution1)
sim1 = ct.ReactorNet([r1])

time1 = 0.0


#run sim to find ignition delay
times1 = np.zeros(100)
data1 = np.zeros((100,2)) #first is time, second is temperature

for n in range(100):
    time1 += 1.e-5
    sim1.advance(time1)
    times1[n] = time1 * 1e3  # time in ms
    data1[n,0] = r1.T
    data1[n, 1] = r1.thermo['OH'].Y

T=np.array(data1[:,0])
dt= np.ones(len(times1)-1)*(times1[1]-times1[0])
dT= np.diff(T)
deriv= dT/dt
i=deriv.argmax()
deriv_max=[times1[i], T[i]]
tau= times1[i]

initial_point=[times1[i-2], T[i-2]]


#read in solution file, and make constant pressure adiabatic reactor
solution2 = ct.Solution('gri30.cti')  #trimmed_h2_v1b_mech.cti
solution2.TPX = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'
r2 = ct.IdealGasConstPressureReactor(solution2)
sim2 = ct.ReactorNet([r2])

time2= tau
times2= np.zeros(100)
n_species=len(solution2.species_names)
data2=np.zeros((100, n_species+2))
species_array=np.ones(n_species)


for n in range(100):
    time2 +=1.e-5
    sim2.advance(time2)
    times2[n] = time2 * 1e3  # time in ms
    data2[n, 0] = time2
    data2[n, 1] = r2.T
    for k in species_array:
        species=solution2.species(k).name
        data2[n, k+2] = r2.thermo[species].Y


# Plot the ignition delay
if '--TAU' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    plt.clf()

    #plot combustion point
    plt.plot(deriv_max[0], deriv_max[1], 'ro')
    plt.plot(initial_point[0], initial_point[1], 'rx')
    #plot temp vs time
    plt.plot(times1, data1[:,0])
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.show()

if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    plt.clf()

    #plot temp vs time
    plt.plot(times2, data2[:,1])
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')

    plt.show()
