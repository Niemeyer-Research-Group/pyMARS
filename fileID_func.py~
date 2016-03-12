# this function reads in the data file , and determines what type of file it is
#it then calls the appropriate function depending on data filetype
#the data file called must be in the same directory


import cantera as ct
import os
import shutil

#used when calling the function locally for testing
exclusion_list=['H2', 'H']

def readin(data_file):
	if data_file.endswith(".xml") or data_file.endswith(".cti"): 
		print("This is an Cantera xml or cti file")
		#from xml_readin_func import xmlreadin				
		#xmlreadin(data_file, exclusion_list)		
	elif data_file.endswith(".inp"):
		print("This is a Chemkin inp file")
	else:
		print("File type not supported")
		

#used when calling the function locally for testing
#readin('gri30.cti')