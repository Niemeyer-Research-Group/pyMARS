


import sys
import numpy as np
import cantera as ct
import os
import progressbar as pb
import h5py




def run_sim(solution_objects, sys_args):
    """Function to run Cantera reactor simulation

    Parameters
    ----------
    Cantera Solution Object
    Command Line arguments
    -------
        Plot of Temp vs Time
        CSV file
        Hdf5 file
    """


    data_file=sys_args.file
    solution1 = solution_objects[0]
    solution1.TPY = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'

    """-------------------------------------------------------------------------
    setup progressbar
    -------------------------------------------------------------------------"""

    #initialize widgets

    widgets = ['Time for loop of 1471 iterations: ', pb.Percentage(), ' ',
                pb.Bar(marker=pb.RotatingMarker()), ' ', pb.ETA()]
    #initialize timer
    timer = pb.ProgressBar(widgets=widgets, maxval=1471).start() #1471

    """-------------------------------------------------------------------------
    run sim to find ignition delay from dT/dt max
    -------------------------------------------------------------------------"""

    #read in solution file, and make constant volume adiabatic reactor


    r1 = ct.Reactor(solution1)
    sim1 = ct.ReactorNet([r1])
    tnow=0.0
    tfinal=5.0e-3
    index1 = 0
    times1 = []
    temps = [] #first column is time, second is temperature
    sdata=np.zeros([0, len(r1.Y)])
    while tnow < tfinal:

        index1 += 1
        tnow = sim1.step(tfinal)
        times1.append(tnow)
        temps.append(r1.T)

        species_data=np.array(r1.Y)
        species_data=species_data[:,np.newaxis].T #translate from [1,n] to [n, 1]
        sdata= np.vstack((sdata, species_data))
        timer.update(index1)
    timer.finish

    #concatenate time and temperature values
    times1=np.array(times1)
    temps=np.array(temps)
    timetemp=np.vstack((times1, temps)).T
    sdata= np.hstack((timetemp, sdata))

    #get ignition point from dT/dt
    T=np.array(temps)
    dt= np.ones(len(times1)-1)*(times1[1]-times1[0])
    dT= np.diff(T)
    deriv= dT/dt
    i=deriv.argmax()
    deriv_max=[times1[i], T[i], i]
    tau= times1[i]


    """-------------------------------------------------------------------------
    find initial and final sample points
    -------------------------------------------------------------------------"""

    #find initial sample point
    for j, dTi in enumerate(dT):
        if dTi > .2:     #when dT > 5 degrees kelvin
            initial_point=[times1[j], T[j], j] # initial point, Temp, index
            break
    try:
        initial_point
    except NameError:
        #sys.exit("initial sample point not found") #alternative sys exit option
        initial_point=[times1[i-20], T[i-20], j]
        print("\nInitial Sample Point not found based on dTmax. Alternative method"+\
        " of 5 steps before tau is used.")

    #find final sample point
    for k, dti in enumerate(dT):
        if k > i:
            if dti < .1:
                final_point=[times1[k], T[k], k]    #final point, Temp, index
                break


    """-------------------------------------------------------------------------
    plot temperature vs time, and ignition point
    -------------------------------------------------------------------------"""

    # Plot the ignition delay
    if sys_args.plot:
        import matplotlib.pyplot as plt
        plt.clf()

        #plot combustion point
        plt.plot(deriv_max[0], deriv_max[1], 'ro', ms=7, label= 'ignition point')
        #plot initial and final sample points
        plt.plot(initial_point[0], initial_point[1], 'rx', ms=5, mew=2)
        plt.plot(final_point[0], final_point[1], 'rx', ms=5, mew=2)
        #plot temp vs time
        plt.plot(times1, temps)

        plt.xlabel('Time (ms)')
        plt.title('Mixture Temperature vs Time')
        plt.legend()
        plt.ylabel('Temperature (K)')
        #plt.axis([0, 1.2, 900, 2800])
        plt.show()


    """-------------------------------------------------------------------------
    write data to csv
    -------------------------------------------------------------------------"""
    if sys_args.writecsv:
        #format matrix for csv
        names=str(solution1.species_names)
        tt=['Time (ms)', 'Temp (K)']
        names = solution1.species_names
        name_array = np.append(tt, names)
        sdata=sdata.astype('|S10')
        file_data= np.vstack((name_array, sdata))
        #open and write to file
        input_file_name_stripped=os.path.splitext(data_file)[0]
        output_file_name=os.path.abspath('Output_Data_Files/'+ 'species_data_' + input_file_name_stripped + '.csv')
        with open(output_file_name, 'wb') as f:
            np.savetxt(f, file_data, fmt=('%+12s'),  delimiter=',')


        os.system('atom '+ output_file_name)

    """-------------------------------------------------------------------------
    write data to hdf5
    -------------------------------------------------------------------------"""

    if sys_args.writehdf5:
        #format matrix for hdf5
        names=str(solution1.species_names)
        tt=['Time (ms)', 'Temp (K)']
        names = solution1.species_names
        name_array = np.append(tt, names)
        sdata=sdata.astype('|S10')
        file_data= np.vstack((name_array, sdata))

        #open and write to file
        input_file_name_stripped=os.path.splitext(data_file)[0]
        output_file_name=os.path.abspath('Output_Data_Files/'+ 'species_data_' + input_file_name_stripped + '.hdf5')
        with h5py.File(output_file_name, 'w') as f:
            Times = f.create_dataset("Times", data=times1)
            Temps = f.create_dataset("Temps", data=temps)
            sgroup= f.create_group('Species_Data')

            for i, sp in enumerate(solution1.species_names):
                    sgroup.create_dataset(sp, data=sdata[:,i+2])
            #nco=np.array(f.get('Species_Data/H2'))
            #print(nco)


    """-------------------------------------------------------------------------
    prints out points of interest
    -------------------------------------------------------------------------"""

    if sys_args.points:
        print("\nTime[ms]    Temp[K]    Index        Point")
        print( str(initial_point[0]) +  "       " + str("{0:.2f}".format(initial_point[1]))\
         + "       " + str(initial_point[2]) + "     " + "Initial sample point")
        print(str(tau) + "        " + str("{0:.2f}".format(deriv_max[1])) + "       " + str(deriv_max[2])\
                + "     " + "Ignition point")
        print( str(final_point[0]) +  "       " + str("{0:.2f}".format(final_point[1]))\
         + "       " + str(final_point[2]) + "     " + "Final sample point")
