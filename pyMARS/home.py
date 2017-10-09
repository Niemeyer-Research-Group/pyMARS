import os, sys, argparse
import cantera as ct
os.environ['Cantera_Data'] =os.getcwd()
from convert_chemkin_file import convert
import soln2ck
import soln2cti
import numpy as np
from run_drgep import run_drgep
from run_drg import run_drg
from autoignition_loop_control import autoignition_loop_control

def readin(args='none', **argv):
    '''
    Main function for pyMARS

    Parameters
    ----------
    file:
        Input mechanism file (ex. file='gri30.cti')
    species:
        Species to eliminate (ex. species='H, OH')
    thermo:
        Thermo data file if Chemkin format (ex. thermo= 'thermo.dat')
    transport:
        Transport data file if Chemkin format
    plot:
        plot ignition curve (ex. plot='y')
    points:
        print ignition point and sample range (ex. points='y')
    writecsv:
        write data to csv (ex. writecsv='y')
    writehdf5:
        write data to hdf5 (ex. writehdf5='y')
    write_ai_times:
        write autoignition times for each inital condition
    run_drg:
        Run DRG model reduction
    run_drgep:
	Run drgep model.
    error:
	Maximum ammount of error allowed. 
    keepers: list of strings
	The string names of the species that should be kept in the model no matter what.
    targets: list of strings
	The string names of the species that should be used as target species.

    Returns
    -------
        Converted mechanism file
        Trimmed Solution Object
        Trimmed Mechanism file

    Examples
    --------
    readin(file='gri30.cti', plot='y', species='OH, H')
    '''

    class args():

        #package from terminal use case
        if args is not 'none':
            plot = args.plot
            points = args.points
            writecsv = True
            writehdf5 = True
            data_file= args.file
            thermo_file = args.thermo
            transport_file = args.transport
            run_drg = args.run_drg
            conditions_file = args.conditions
            convert = args.convert
            error = args.error
            run_drgep = args.run_drgep
            write_ai_times = args.write_ai_times
            target = 0
            if args.species is None:
                keepers = []
            else:
                keepers = [str(item) for item in args.species.split(',')]
                #strip spaces
                for i, sp in enumerate(keepers):
                    keepers[i]=sp.strip()
            if args.target is None:
                target = []
            else:
                target = [str(item) for item in args.target.split(',')]
                #strip spaces
                for i, sp in enumerate(target):
                    target[i]=sp.strip()

    file_extension= os.path.splitext(args.data_file)[1]

    if file_extension == ".cti" or file_extension == ".xml": #If the file is a Cantera file.
        print("\nThis is a Cantera xml or cti file\n")
        solution_object = ct.Solution(args.data_file)
    
        #runs simulation once with additional features on
        if args.plot is True or args.writecsv is True or args.points is True or args.writehdf5 is True or args.write_ai_times is True:
            if os.path.exists('mass_fractions.hdf5'):
                os.system('rm mass_fractions.hdf5')
            if args.write_ai_times is True:
                if os.path.exists('autoignition_times.txt'):
                    os.system('rm autoignition_times.txt')
            print('running simulation\n')
            sim_result = autoignition_loop_control(solution_object, args, True)
	
        if args.run_drg is True:
        	run_drg(args, solution_object)
	
        if args.convert is True:
            soln2ck.write(solution_object)
	
        if args.run_drgep is True: #If the user wants to run drgep and specifies it as a command line argument.
            run_drgep(args, solution_object)

    elif file_extension == ".inp" or file_extension == ".dat" or file_extension == ".txt":
        print("\n\nThis is a Chemkin file")
        #convert file to cti
        converted_file_name = convert(args.data_file, args.thermo_file, args.transport_file)

    else:
        print("\n\nFile type not supported")

