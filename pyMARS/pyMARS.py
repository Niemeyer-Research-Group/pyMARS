
import os, sys, argparse


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

    if args is 'none':
        args=argparse.Namespace()
        args.plot = False
        args.points = False
        args.writecsv = False
        args.writehdf5 = False

        data_file = argv['file']
        if 'thermo' in argv:
            thermo_file = argv['thermo']
        if 'transport' in argv:
            transport_file = argv['transport']
        if 'species' in argv:
            species = argv['species']
            exclusion_list = [str(item) for item in species.split(',')]
            print exclusion_list
        else:
            exclusion_list=[]
        if 'plot' in argv:
            args.plot = True
        if 'writecsv' in argv:
            args.writecsv = True
        if 'writehdf5' in argv:
            args.writehdf5 = True
        if 'points' in argv:
            args.points = True
        print args

    else:
        data_file= args.file
        thermo_file = args.thermo
        transport_file=args.transport
        if args.species == None:
            exclusion_list=[]
        else:
            exclusion_list=[str(item) for item in args.species.split(',')]
            print(exclusion_list)


    print args
    ext= os.path.splitext(data_file)[1]

    if ext == ".cti" or ext == ".xml":

        print("\nThis is an Cantera xml or cti file\n")
        #trims file
        solution_objects=trim(data_file, exclusion_list)
        write(solution_objects[1])
        run_sim(data_file, args)

    elif ext == ".inp" or ext == ".dat" or ext == ".txt":

        print("This is a Chemkin file")
        #convert file to cti
        converted_file_name = convert(data_file, thermo_file, transport_file)

        #trims newly converted file
        solution_objects=trim(converted_file_name, \
                                    exclusion_list)
        print(converted_file_name)
        run_sim(converted_file_name, args)


    else:
        print("File type not supported")
