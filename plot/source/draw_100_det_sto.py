
import sys, os
from subprocess import call

sto_sim = "../../stochastic/stochastic/stochastic"
det_sim = "../../deterministic2013Version/deterministic/deterministic1"
sto_dir = "../sto"
det_dir = "../det"
ds_fig_dir = "../ds_figures"

def main():
	sto_command = [sto_sim, "-i", "../params_100_compare.txt", "-ns", str(100), "-pc", "-od", sto_dir]
	det_command = [det_sim, "-i", "../params_100_compare.txt", "-ns", str(100), "-pc", "-od", det_dir]
	if 1 == call(sto_command):
		exit(1)
	#if 1 == call(det_command):
	#	exit(1)
	for i in range(100):
		command = ["/Applications/MATLAB_R2013a.app/bin/matlab", "-nodesktop", "-nosplash", "-nodisplay", "-r", "plot_det_sto_data(\'" + det_dir + "/set_" + str(i) + "/wildtype.txt\',\'" + sto_dir+ "/set_" + str(i+1)+ "/wildtype/mHer1_smoothed.txt\'," + str(0) + ", \'" + ds_fig_dir + "/ds_fig_" + str(i) + ".fig\')"]
		if 1 == call(command):
			exit(1)
main()
