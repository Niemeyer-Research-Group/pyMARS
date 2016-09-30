import os, sys, argparse
import cantera as ct
os.environ['Cantera_Data'] =os.getcwd()

from create_trimmed_model import trim
from convert_chemkin_file import convert
from soln2cti import write
from autoignition_module import run_sim
from get_rate_data import get_rates
from drg import make_graph
from drg_loop_control import loop_control


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
        #direct use case
            if args is 'none':
                argparse.Namespace()
                plot = False
                points = False
                writecsv = False
                writehdf5 = False
                data_file = argv['file']
                thermo_file = None
                transport_file = None
                run_drg = None
                initial_sim = True
                iterate = False

                if 'thermo' in argv:
                    thermo_file = argv['thermo']
                if 'transport' in argv:
                    transport_file = argv['transport']
                if 'species' in argv:
                    species = argv['species']
                    exclusion_list = [str(item) for item in species.split(',')]
                    #strip spaces
                    for i, sp in enumerate(exclusion_list):
                        exclusion_list[i]=sp.strip()
                if 'species' not in argv:
                    exclusion_list=[]
                if 'plot' in argv:
                    plot = True
                if 'writecsv' in argv:
                    writecsv = True
                if 'writehdf5' in argv:
                    writehdf5 = True
                if 'points' in argv:
                    points = True
                if 'run_drg' in argv:
                    run_drg = True
                if 'iterate' in argv:
                    iterate = True
                x ='arg_none'

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
                initial_sim = True
                iterate = args.iterate
                if args.species is None:
                    exclusion_list = []
                else:
                    exclusion_list = [str(item) for item in args.species.split(',')]
                    #strip spaces
                    for i, sp in enumerate(exclusion_list):
                        exclusion_list[i]=sp.strip()
                x='args_not_none'

    file_extension= os.path.splitext(args.data_file)[1]

    if file_extension == ".cti" or file_extension == ".xml":
        print("\n\nThis is an Cantera xml or cti file\n")
        solution_object = ct.Solution(args.data_file)
        #trims file
        #need case if no trim necessary
        solution_objects = trim(solution_object, args.exclusion_list, args.data_file)
        if args.run_drg is False:
            trimmed_file = write(solution_objects[1])
        if args.plot is True or args.writecsv is True or args.points is True or args.writehdf5 is True:
            print 'running sim'
            sim_result = run_sim(solution_object, args)
        if args.run_drg is True:
            new_solution_objects = loop_control(solution_object, args)
            n_species_eliminated = len(solution_object.species())-len(new_solution_objects[1].species())
            print 'Number of species eliminated: %s' %n_species_eliminated
            drg_trimmed_file = write(new_solution_objects[1])

    elif file_extension == ".inp" or file_extension == ".dat" or file_extension == ".txt":
        print("\n\nThis is a Chemkin file")
        #convert file to cti
        converted_file_name = convert(args.data_file, args.thermo_file, args.transport_file)
        #trims newly converted file
        converted_file_object = ct.Solution(converted_file_name)
        solution_objects = trim(converted_file_object, args.exclusion_list, args.data_file)
        trimmed_file = write(solution_objects[1])

        if "plot" or "points" or "writecsv" or "writehdf5" in args:
            print 'running sim'
            #run_sim(converted_file_name, args)

        if 'run_drg' in args:
            print 'running sim'
            run_sim(converted_file_name, args)
            get_rates('mass_fractions.hdf5', converted_file_name)
            print 'running DRG'
            drg_exclusion_list = make_graph(converted_file_name, 'production_rates.hdf5')
            new_solution_objects = trim(converted_file_name, drg_exclusion_list)
            args.data_file = write(new_solution_objects[1])


    else:
        print("\n\nFile type not supported")
