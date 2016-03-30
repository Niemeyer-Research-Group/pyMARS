
#local import
import cantera as ct

def readin(data_file):
	"""Function to import data file and identify format.

	Parameters
	----------
	data_file:
		Local Chemkin or Cantera data file containing mechanism information

	Returns
	-------
		Calls function corresponding to data file
	"""

	if data_file.endswith(".xml") or data_file.endswith(".cti"):
		print("This is an Cantera xml or cti file")
		from readin_func import datareadin
		datareadin(data_file, exclusion_list)
	elif data_file.endswith(".inp"):
		print("This is a Chemkin inp file")
        ck2cti --input=gri30.inp
	#else:
		#print("File type not supported")

readin('gri30.inp')
