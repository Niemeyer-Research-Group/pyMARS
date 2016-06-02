

import sys
import numpy as np
import cantera as ct


#read in solution file, and make constant pressure adiabatic reactor
gri3 = ct.Solution('gri30.cti')  #trimmed_h2_v1b_mech.cti
gri3.TPX = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'
r = ct.IdealGasConstPressureReactor(gri3)
sim = ct.ReactorNet([r])

time = 0.0


#run sim to find ignition delay
times = np.zeros(100)
data = np.zeros((100,2))

for n in range(100):
    time += 1.e-5
    sim.advance(time)
    times[n] = time * 1e3  # time in ms
    data[n,0] = r.T

T=np.array(data[:,0])

dt= np.ones(len(times)-1)*(times[1]-times[0])
dT= np.diff(T)
deriv= dT/dt
i=deriv.argmax()
deriv_max=[times[i], T[i]]
tau= times[i]
print(deriv_max)


# Plot the results when --plot is called from terminal
if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    plt.clf()

    #plot combustion point
    plt.plot(deriv_max[0], deriv_max[1], 'ro')
    #plot temp vs time
    plt.plot(times, data[:,0])
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')

    plt.show()
