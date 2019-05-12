
def name_trim(string):
	"""Changes string with file path and returns file alone.

	Parameters
	----------
	string : str
		Full path to a file
	
	Returns
	-------
	string : str
		Filename extracted from full path
	"""
	first_index = 0
	for i in range (0,len(string)):
		if (string[i] == '/'):
			first_index = i+1
	if (string[len(string) - 1] == '/'):
		print("Invalid file name.")
		exit()
	string = string[first_index:len(string)]
	
#	last_index = 0
#	for i in range (0,len(string)):
#		if (string[i] == '.'):
#			last_index = i
#	if (string[len(string) - 1] == '.'):
#		print("Invalid file name.")
#		exit()
#	string = string[0:last_index]

	return string
