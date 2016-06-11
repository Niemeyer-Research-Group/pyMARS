

import sys
import numpy as np
import cantera as ct


#read in solution file, and make constant volume adiabatic reactor
solution1 = ct.Solution('trimmed_h2_v1b_mech.cti')
solution1.TPX = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'
r1 = ct.IdealGasReactor(solution1)
sim1 = ct.ReactorNet([r1])

time1 = 0.0


"""----------------------------------------------------------------------------
run sim to find ignition delay from dT/dt max
-----------------------------------------------------------------------------"""


times1 = np.zeros(100)
data1 = np.zeros((100,2)) #first column is time, second is temperature

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


#finds initial sample point
for j, dti in enumerate(dT):
    if dti > 5:     #when dT > 5 degrees kelvin
        initial_point=[times1[j], T[j], j] # initial point, Temp, index
        break


#finds final sample point
for k, dti in enumerate(dT):
    if k > i:
        if dti < 1:
            final_point=[times1[k], T[k], k]
            break

#prints points of interest
if '--points' in sys.argv[1:]:
    print("Time[ms]    Temp[K]    Index        Point")
    print( str(initial_point[0]) +  "       " + str("{0:.2f}".format(initial_point[1]))\
     + "       " + str(initial_point[2]) + "     " + "Initial sample point")
    print(str(tau) + "        " + str("{0:.2f}".format(T[i])) + "       " + str(i)\
            + "     " + "Ignition point")
    print( str(final_point[0]) +  "       " + str("{0:.2f}".format(final_point[1]))\
     + "       " + str(final_point[2]) + "     " + "Final sample point")

"""----------------------------------------------------------------------------
run sim to get data around ignition
-----------------------------------------------------------------------------"""


#read in solution file, and make constant volume adiabatic reactor
solution2 = ct.Solution('trimmed_h2_v1b_mech.cti')
solution2.TPX = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'
r2 = ct.IdealGasReactor(solution2)
sim2 = ct.ReactorNet([r2])

time2= tau*10**-3
times2= np.zeros(100)
n_species=len(solution2.species_names)
data2=np.zeros((100, n_species+2))
species_array=np.ones(n_species)


for n in range(100):
    time2 +=1.e-6
    sim2.advance(time2)
    times2[n] = time2 * 1e3  # time in ms
    data2[n, 0] = time2
    data2[n, 1] = r2.T
    for k in species_array:
        species=solution2.species(k).name
        data2[n, k+2] = r2.thermo[species].Y


"""----------------------------------------------------------------------------
plot temperature vs time, and ignition point
-----------------------------------------------------------------------------"""



# Plot the ignition delay
if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    plt.clf()

    #plot combustion point
    plt.plot(deriv_max[0], deriv_max[1], 'ro')
    #plot points where dT > 5 degrees kelvin
    plt.plot(initial_point[0], initial_point[1], 'rx', ms=5, mew=2)
    plt.plot(final_point[0], final_point[1], 'rx', ms=5, mew=2)
    #plot temp vs time
    plt.plot(times1, data1[:,0])
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.axis([0, 1.2, 900, 2800])
    plt.show()
