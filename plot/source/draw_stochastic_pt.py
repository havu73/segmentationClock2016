import sys, os
from subprocess import call

sto_sim = "../../stochastic/stochastic/stochastic"
good_sets = range(1)
dir_head = "../sto/set_"
dir_tail = "/wildtype"
out_dir = "../peak_trough_fig"

def run_stochastic(param_file, num_sets, out_dir):
	sto_command = [sto_sim, "-i", param_file, "-ns", str(num_sets), "-pc", "-od", out_dir]
	if 1 == call(sto_command):
		exit(1)

def create_param_file(index_set, in_name, out_name):
	fin = open(in_name, "r")
	fout = open(out_name, "w")
	lines = fin.readlines()
	for index in index_set:
		fout.write(lines[index])
	fin.close()
	fout.close()
	

def draw_sto_data(in_dir, out_dir, set_indices):
	for i, set_index in enumerate(set_indices):
		sto_dir = in_dir + "/set_" + str(i) + "/wildtype"
		raw = sto_dir + "/mHer1.txt"
		smooth = sto_dir + "/mher1_smoothed.txt"
		peak = sto_dir + "/mHer1_peaks.txt"
		trough = sto_dir + "/mHer1_troughs.txt"
		outfile = out_dir + "/sto_" + str(set_index) + "_fig.fig"
		command = ["/Applications/MATLAB_R2013a.app/bin/matlab", "-nodesktop", "-nosplash", "-nodisplay", \
		"-r", "plot_smooth_data(\'" + raw + "\',\'" + smooth + "\',1,\'" + peak + "\',\'" + trough\
		+ "\',\'" + outfile + "\')"]
		print command
		if 1 == call(command):
			print "Problems running ", command
			exit(1)

def main():
	param_100 = "../params_100_compare.txt"
	param_file  = "../experimental_sets.txt"
	sto_out_dir = "../sto"
	out_dir = "../sto/check_new_smooth"
	#create_param_file(good_sets, param_100, param_file)
	#run_stochastic(param_file, 1, out_dir)
	draw_sto_data(out_dir, out_dir, good_sets)
	
	#for i in range(len(good_sets)):
		#command = ["rm", "-r", sto_out_dir + "/set_" + str(i)]
		#if 1 == call(command):
			#exit(1)
	#command = ["rm", param_file]
	#if 1 == call(command):
		#exit(1)

#main()
draw_sto_data("../../stochastic/stochastic", "../../stochastic/stochastic/set_0", [0])

