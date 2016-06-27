
import os, sys, argparse
import cantera as ct


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

    class args_none:
        def __init__(self, argv):
            self=argparse.Namespace()
            self.plot = False
            self.points = False
            self.writecsv = False
            self.writehdf5 = False
            self.data_file = argv['file']
        if 'thermo' in argv:
            self.thermo_file = argv['thermo']
        if 'transport' in argv:
            self.transport_file = argv['transport']
        if 'species' in argv:
            self.species = argv['species']
            self.exclusion_list = [str(item) for item in species.split(',')]
            print self.exclusion_list
        #if 'species' not in argv:
            #self.exclusion_list=[]

        if 'plot' in argv:
            self.plot = True
        if 'writecsv' in argv:
            self.writecsv = True
        if 'writehdf5' in argv:
            self.writehdf5 = True
        if 'points' in argv:
            self.points = True

    class args_not_none:
        def __init__(self, args):
            self.plot = args.plot
            self.points = args.points
            self.writecsv = args.writecsv
            self.writehdf5 = args.writehdf5
            self.data_file= args.file
            self.thermo_file = args.thermo
            self.transport_file=args.transport
            if args.species == None:
                self.exclusion_list=[]
            else:
                self.exclusion_list=[str(item) for item in args.species.split(',')]
                print(self.exclusion_list)

    #if function is used directly
    if args is 'none':
        args = args_none(argv)
    #function used from command line
    else:
        args = args_not_none(args)
    ext= os.path.splitext(args.data_file)[1]

    if ext == ".cti" or ext == ".xml":
        print("\nThis is an Cantera xml or cti file\n")
        #trims file
        solution_objects=trim(args.data_file, args.exclusion_list)
        write(solution_objects[1])
        run_sim(args.data_file, args)

    elif ext == ".inp" or ext == ".dat" or ext == ".txt":
        print("This is a Chemkin file")
        #convert file to cti
        converted_file_name = convert(args.data_file, args.thermo_file, args.transport_file)
        #trims newly converted file
        solution_objects=trim(converted_file_name, exclusion_list)
        print(converted_file_name)
        run_sim(converted_file_name, args)

    else:
        print("File type not supported")
