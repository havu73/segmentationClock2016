import sys
import shared
from collections import Counter

def write_parameters(outF, parameters):
	"""
	outFname: output file object
	parameters: list of strings that indicate the parameters found by sres
	This function write into outFname the parameters we read from sres' results
	"""
	for i, param in enumerate(parameters[:-1]):
		outF.write(param + ",")
	outF.write(parameters[-1])
	
def process_sres_good_folders(folderIn, inFname, run_indices, score_thres, outFname):
	print "processing input files to filter out good parameter sets..."
	outF = open(outFname, 'w')
	for i, index in enumerate(run_indices):
		fIn = open(folderIn + inFname + str(index) + ".txt", 'r')
		for line in fIn: 
			parameters = line.split(",")
			score = float(parameters[0]) # The first number of a good parameter set line was written in sres indicates the score that
			# the parameter sets receive
			if score <= score_thres: # a good parameter set will have score <= the threshold provided by user
				write_parameters(outF, parameters[1:])
	outF.close()
	
def get_param_value(args, index):
	try:
		value = args[index]
	except:
		print "Program not inputted correctly. Exiting...."
		usage()
		exit(1)
	return value

def check_missing_parameter(mendatory_param, cnt):
	"""
	Check whether or not all the required parameters are provided by the user
	"""
	for i, param_name in enumerate(mendatory_param):
		if cnt[param_name] == 0:
			print "Missing parameter: ", param_name
			usage()
			exit(1)
	
def main():
	print 'Reading command-line arguments...'
	args = sys.argv[1:] # Remove the name of the program from the arguments
	num_args = len(args)
	mendatory_param = ['folderIn', 'inFname', 'start_file_index', 'end_file_index', 'out_fname', 'score_threshold'] # List of names of parameter that the users have to provide
	c = Counter([]) # this is used to keep track of what parameters has been inputted by the user
	i = 0
	while i < len(args):
		arg = args[i]			
		if arg == '-f': 
			folderIn = get_param_value(args, i + 1)
			c['folderIn'] = 1
			i += 2
		elif arg == '-ifn':
			inFname = get_param_value(args, i + 1)
			c['inFname'] = 1
			i += 2
		elif arg == '-sfi':
			start_file_index = int(get_param_value(args, i + 1))
			c['start_file_index'] = 1
			i += 2
		elif arg == '-efi':
			end_file_index = int (get_param_value(args, i + 1))
			c['end_file_index'] = 1
			i += 2
		elif arg == '-of':
			out_fname = get_param_value(args, i + 1)
			c['out_fname'] = 1
			i += 2
		elif arg == '-st':
			score_threshold = float(get_param_value(args, i + 1))
			c['score_threshold'] = 1
			i += 2
		elif arg == '-h':
			usage()
			exit(0)
		elif arg == '-l':
			license
			exit(0)
	check_missing_parameter(mendatory_param, c)
	process_sres_good_folders(folderIn, inFname, range(start_file_index, end_file_index), score_threshold, out_fname)
	print "Done"

def usage():
	print "This file is used to get all parameter sets havings scores below a threshold and combined them into one file, with the score eliminated (only parameters left)."
	print "After running parameter search using sres, each process will print out good parameter sets (with scores as the first numbers of each line) into separate files."
	print "This program will combine them all into one file, accept only sets with scores less than an user-specified threshold"
	print "This program output will be used as input for refine_ranges.py"
	print "-f			:	path to folder that contains all output of good set files of sres search, relative to this file."
	print "-ifn			:	initial names of files. Most of the time, files have the same beginning. The only difference is the indices of files"
	print "-sfi			:	the index of the starting file. For example, if my files are wt0.txt, wt1.txt, wt2.txt, then start_file_index shoud be 0."
	print "-efi			:	the index of the ending file. For example, if my files are wt0.txt, wt1.txt, wt2.txt, then end_file_index shoud be 3."
	print "-of			:	the name of the output file."
	print "-st			: 	score threshold. If parameter sets have score less than this threshold, it will be printed into output file."
	print "-h			: 	get guidance about how to run this program."
	print "-l			:	get information about licenses of this program."
	
def license():
	print "This program is created by Ha Vu and Ahmet Ay at Colgate University, 2017."
				
main()
