import os, sys, argparse
import cantera as ct
os.environ['Cantera_Data'] =os.getcwd()

from create_trimmed_model import trim
from convert_chemkin_file import convert
from write_to_cti import write
from autoignition_module import run_sim


def readin(args='none', **argv):
    """Main function for pyMARS

    Arguments
        file: Input mechanism file (ex. file='gri30.cti')
        species: Species to eliminate (ex. species='H, OH')
        thermo: Thermo data file if Chemkin format (ex. thermo= 'thermo.dat')
        transport: Transport data file if Chemkin format
        plot: plot ignition curve (ex. plot='y')
        points: print ignition point and sample range (ex. points='y')
        writecsv: write data to csv (ex. writecsv='y')
        writehdf5: write data to hdf5 (ex. writehdf5='y')

    ----------
    Returns
        Converted mechanism file
        Trimmed Solution Object
        Trimmed Mechanism file
    ----------
    Example
        readin(file='gri30.cti', plot='y', species='OH, H')
    """

    "--------------------------------------------------------------------------"
    "--------------------------------------------------------------------------"

    class args():
            if args is 'none':
                argparse.Namespace()
                plot = False
                points = False
                writecsv = False
                writehdf5 = False
                data_file = argv['file']
                thermo_file = None
                transport_file = None
                if 'thermo' in argv:
                    thermo_file = argv['thermo']
                if 'transport' in argv:
                    transport_file = argv['transport']
                if 'species' in argv:
                    species = argv['species']
                    exclusion_list = [str(item) for item in species.split(',')]
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
                x ='arg_none'
            if args is not 'none':
                plot = args.plot
                points = args.points
                writecsv = args.writecsv
                writehdf5 = args.writehdf5
                data_file= args.file
                thermo_file = args.thermo
                transport_file=args.transport
                if args.species is None:
                    exclusion_list=[]
                else:
                    exclusion_list=[str(item) for item in args.species.split(',')]
                    print(exclusion_list)
                x='args_not_none'
    ext= os.path.splitext(args.data_file)[1]

    if ext == ".cti" or ext == ".xml":
        print("\n\nThis is an Cantera xml or cti file\n")
        #trims file
        solution_objects=trim(args.data_file, args.exclusion_list)
        args.data_file=write(solution_objects[1])
        if args.plot is True:
            print 'running sim'
            run_sim(args.data_file, args)

    elif ext == ".inp" or ext == ".dat" or ext == ".txt":
        print("\n\nThis is a Chemkin file")
        #convert file to cti
        converted_file_name = convert(args.data_file, args.thermo_file, args.transport_file)
        #trims newly converted file
        solution_objects=trim(converted_file_name, args.exclusion_list)
        args.data_file=write(solution_objects[1])

        if "plot" or "points" or "writecsv" or "writehdf5" in args:
            print 'running sim'
            run_sim(args.data_file, args)

        #run_sim(converted_file_name, args)



    else:
        print("\n\nFile type not supported")
