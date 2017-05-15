
import sys, os
from subprocess import call
sto_sim = "../../stochastic/stochastic/stochastic"

def run_stochastic(param_file, num_sets, num_rounds, out_dir):
	sto_command = [sto_sim, "-i", param_file, "-ns", str(num_sets), "-pc", "-od"]
	for i in range(num_rounds):
		command = sto_command + [out_dir + str(i)]
		if 1 == call(command):
			exit(1)
		
def plot_overlapping_stochastic(in_dir_format, num_round, num_sets, out_dir):
	for i in range(num_sets):
		command = ["/Applications/MATLAB_R2013a.app/bin/matlab", "-nodesktop", "-nosplash", "-nodisplay", \
		"-r"]
		function_call = "plot_overlapping_stochastic(" +str(num_round) +","
		for j in range(num_round):
			file_name = "\'" + in_dir_format + str(j) + "/set_" + str(i) + "/wildtype/mHer1.txt\'"
			function_call += file_name + ","
		out_file = "\'" + out_dir + "/overlapping_set_" + str(i) +".fig\'"
		function_call += out_file + ")"
		command.append(function_call)
		if 1 == call(command):
			exit(1)
		

def main():
	out_dir = "../sto/check_random_round"
	params_file = "../params_100_compare.txt"
	#run_stochastic(params_file, 3, 5, out_dir)
	plot_overlapping_stochastic("../sto/check_random_round", 5, 3, "../sto/check_random")
main()
