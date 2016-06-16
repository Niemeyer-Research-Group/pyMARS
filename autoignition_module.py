import sys
import numpy as np
import cantera as ct
import os


#read in solution file, and make constant volume adiabatic reactor
solution1 = ct.Solution('gri30.cti') #trimmed_h2_v1b_mech.cti
solution1.TPY = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'
r1 = ct.Reactor(solution1)
sim1 = ct.ReactorNet([r1])
time1 = 0.0


"""----------------------------------------------------------------------------
run sim to find ignition delay from dT/dt max
-----------------------------------------------------------------------------"""


times1 = np.zeros(100)
data1 = np.zeros((100,2)) #first column is time, second is temperature

for n in range(100):
    time1 += 1.e-4
    sim1.advance(time1)
    times1[n] = time1 * 1e3  # time in ms
    data1[n,0] = r1.T


#get ignition point from dT/dt
T=np.array(data1[:,0])
dt= np.ones(len(times1)-1)*(times1[1]-times1[0])
dT= np.diff(T)
deriv= dT/dt
i=deriv.argmax()
deriv_max=[times1[i], T[i]]
tau= times1[i]


#find initial sample point
for j, dTi in enumerate(dT):
    if dTi > 15.0:     #when dT > 5 degrees kelvin
        initial_point=[times1[j], T[j], j] # initial point, Temp, index
        break
    else:
        #sys.exit("initial sample point not found") #alternative sys exit option
        initial_point=[times1[i-5], T[i-5], j]
        print("\nInitial Sample Point not found based on dTmax. Alternative method"+\
            " of 5 steps before tau is used.")
        break

#find final sample point
for k, dti in enumerate(dT):
    if k > i:
        if dti < 1:
            final_point=[times1[k], T[k], k]    #final point, Temp, index
            break

#prints points of interest
if '--points' in sys.argv[1:]:
    print("\nTime[ms]    Temp[K]    Index        Point")
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
solution2 = ct.Solution('gri30.cti')  #trimmed_h2_v1b_mech.cti
solution2.TPY = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'
r2 = ct.Reactor(solution2)
sim2 = ct.ReactorNet([r2])

time2= initial_point[0]*10e-4 #initial_point[0]

refined_steps= 200
times2= np.zeros(refined_steps)
n_species=len(solution2.species_names)
data2=np.zeros((refined_steps, n_species+2))
species_array=np.ones(n_species)
for n in range(refined_steps):
    time2 +=1.0e-5
    sim2.advance(time2)
    times2[n] = time2 * 1e3  # time in ms
    data2[n, 0] = time2
    data2[n, 1] = r2.T
    for i, val in enumerate(r2.Y):                     #get species data
        species=solution2.species(i).name
        data2[n, int(i+2)] = val
    if time2 > final_point[0]*10e-4:        #breaks sim, and trims array
        data2 = data2[0:n, :]
        times2=times2[0:n]
        break



"""----------------------------------------------------------------------------
plot temperature vs time, and ignition point
-----------------------------------------------------------------------------"""

# Plot the ignition delay
if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    plt.clf()

    #plot combustion point
    plt.plot(deriv_max[0], deriv_max[1], 'ro', label= 'ignition point')

    #plot initial and final sample points
    plt.plot(initial_point[0], initial_point[1], 'rx', ms=5, mew=2)
    plt.plot(final_point[0], final_point[1], 'rx', ms=5, mew=2)

    #plot temp vs time
    plt.plot(times1, data1[:,0], lw=1.5)
    plt.plot(times2, data2[:,1], ls='--', lw=2.5, label='sample range')
    plt.xlabel('Time (ms)')
    plt.title('Mixture Temperature vs Time')
    plt.legend()
    plt.ylabel('Temperature (K)')
    #plt.axis([0, 1.2, 900, 2800])
    plt.show()


"""----------------------------------------------------------------------------
write data to csv
-----------------------------------------------------------------------------"""
if '--writecsv' in sys.argv[1:]:
    #format matrix for csv
    names=str(solution2.species_names)
    tt=['Time (ms)', 'Temp (K)']
    names = solution2.species_names
    name_array = np.append(tt, names)
    data2=data2.astype('|S10')
    file_data= np.vstack((name_array, data2))
    #open and write to file
    with open('test.csv', 'wb') as f:
        np.savetxt(f, file_data, fmt=('%+12s'),  delimiter=',')


    os.system('atom test.csv')
