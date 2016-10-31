import os, sys, argparse
import cantera as ct
os.environ['Cantera_Data'] =os.getcwd()
from create_trimmed_model import trim
from convert_chemkin_file import convert
import soln2ck
import soln2cti
from autoignition_module import run_sim
from get_rate_data import get_rates
from drg import make_graph
from drg_loop_control import drg_loop_control
from autoignition_loop_control import autoignition_loop_control


def readin(args='none', **argv):
    """Main function for pyMARS

    :param file:
        Input mechanism file (ex. file='gri30.cti')
    :param species:
        Species to eliminate (ex. species='H, OH')
    :param thermo:
        Thermo data file if Chemkin format (ex. thermo= 'thermo.dat')
    :param transport:
        Transport data file if Chemkin format
    :param plot:
        plot ignition curve (ex. plot='y')
    :param points:
        print ignition point and sample range (ex. points='y')
    :param writecsv:
        write data to csv (ex. writecsv='y')
    :param writehdf5:
        write data to hdf5 (ex. writehdf5='y')
    :param run_drg:
        Run DRG model reduction

    :returns:
        Converted mechanism file
        Trimmed Solution Object
        Trimmed Mechanism file

    >>>readin(file='gri30.cti', plot='y', species='OH, H')
    """

    class args():

        #package from terminal use case
        if args is not 'none':
            plot = args.plot
            points = args.points
            writecsv = args.writecsv
            writehdf5 = args.writehdf5
            data_file= args.file
            thermo_file = args.thermo
            transport_file = args.transport
            run_drg = args.run_drg
            multiple_conditions = True
            threshold_values = args.thresholds
            conditions_file = args.conditions
            convert = args.convert
            if args.species is None:
                exclusion_list = []
            else:
                exclusion_list = [str(item) for item in args.species.split(',')]
                #strip spaces
                for i, sp in enumerate(exclusion_list):
                    exclusion_list[i]=sp.strip()

    file_extension= os.path.splitext(args.data_file)[1]

    if file_extension == ".cti" or file_extension == ".xml":
        print("\nThis is an Cantera xml or cti file\n")
        solution_object = ct.Solution(args.data_file)
        if args.plot is True or args.writecsv is True or args.points is True or args.writehdf5 is True:
            print 'running sim'
            sim_result = autoignition_loop_control(solution_object, args)
        if args.run_drg is True:
            new_solution_objects = drg_loop_control(solution_object, args)
            os.system('rm production_rates.hdf5')
            os.system('rm mass_fractions.hdf5')
            drg_trimmed_file = soln2cti.write(new_solution_objects[1])
            try:
                os.system('rm production_rates.hdf5')
            except Exception:
                pass
        if args.convert is True:
            soln2ck.write(solution_object)


    elif file_extension == ".inp" or file_extension == ".dat" or file_extension == ".txt":
        print("\n\nThis is a Chemkin file")
        #convert file to cti
        converted_file_name = convert(args.data_file, args.thermo_file, args.transport_file)

    else:
        print("\n\nFile type not supported")
