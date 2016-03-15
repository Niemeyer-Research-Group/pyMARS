#read in data file, determine type and call the appropriate function


import cantera as ct


#used when calling the function locally for testing
#exclusion_list=['O2']

def readin(data_file):
	if data_file.endswith(".xml") or data_file.endswith(".cti"): 
		print("This is an Cantera xml or cti file")
		from readin_func import datareadin				
		datareadin(data_file, exclusion_list)		
	elif data_file.endswith(".inp"):
		print("This is a Chemkin inp file")
	else:
		print("File type not supported")
		

#used when calling the function locally for testing
#readin('gri30.cti')

