import sys, os
from subprocess import call
import numpy, math
import matplotlib.pyplot as plt
import matplotlib.legend_handler as handler
import matplotlib.patches as patches
import matplotlib.lines as mlines
from matplotlib import rc # text style 
rc('text', usetex=True) # activate latex text rendering

noise_colors = ['#2A6EFF','m','g', '#000000']
amp_colors = ['b','r']
black = '#000000'
sto_sim = "../../stochastic/stochastic/stochastic"

def determineTickInterval(r,l): # determine tick interval given a range (r)
	# r: range
	# l: limit (increase l for more ticks)
	candidates = [1,2,5,10,20,30,50,100]	
	for candidate in candidates:
		if r/candidate<l:
			return candidate
	return 1

def updateTicklabels(ax):
	xlabels = [format(label, r',.0f') for label in ax.get_xticks()]
	ax.set_xticklabels(xlabels)
	ax.tick_params(axis='x', pad=10)
	ylabels = [format(label, r',.1f') for label in ax.get_yticks()]
	ax.set_yticklabels(ylabels)

def updateYTicklabels(ax):
	ylabels = [format(label, r',.0f') for label in ax.get_yticks()] # intergers
	ax.set_yticklabels(ylabels)

def updateLogTicklabels(ax): # notation needed for log-log plot
	xlabels = [r'$\textbf{%.0f}$' % 10**label for label in ax.get_xticks()]
	ax.set_xticklabels(xlabels)
	ax.tick_params(axis='x', pad=10)
	ylabels = [r'$\textbf{%.2f}$' % 10**label for label in ax.get_yticks()] 
	ax.set_yticklabels(ylabels)

def read_features_data(filename):
	f = open(filename, "r")
	lines = f.readlines()
	basic_data = (lines[0]).split(",")
	num_cells = int(basic_data[0])
	num_bins = int(basic_data[1])
	amplitude_h1 = []
	amplitude_h7 = []
	period_h1 = []
	period_h7 = []
	avg_cons = []
	in_noise = []
	ex_noise = []
	# get lines of data
	ah1 = (lines[1]).split(",")
	ah7 = (lines[2]).split(",")
	ph1 = (lines[3]).split(",")
	ph7 = (lines[4]).split(",")
	for i in range(num_cells):
		amplitude_h1.append(float(ah1[i]))
		amplitude_h7.append(float(ah7[i]))
		period_h1.append(float(ph1[i]))
		period_h7.append(float(ph7[i]))
	binned_cons = (lines[5]).split(",")
	binned_in = (lines[6]).split(",")
	binned_ex = (lines[7]).split(",")
	for i in range(num_bins):
		avg_cons.append(float(binned_cons[i]))
		in_noise.append(float(binned_in[i]))
		ex_noise.append(float(binned_ex[i]))
	return num_cells, num_bins, amplitude_h1, amplitude_h7, period_h1, period_h7, avg_cons, in_noise, ex_noise

EXP_AVG_CONS = [31.33330984, 77.41325691, 123.208366]
EXP_TOT_NOISE = [0.468385227, 0.286143297, 0.207439948]
EXP_IN_NOISE = [0.200901981, 0.076300351, 0.036878339]
EXP_EX_NOISE = [0.267483246, 0.209842946, 0.17056161]
EXP_NOISE = [EXP_TOT_NOISE, EXP_IN_NOISE, EXP_EX_NOISE]
EXP_NUM_BINS = 3

def calculate_exp_bounds(noise_list, multiplier):
	result = []
	for num in noise_list:
		result.append(num* multiplier)
	return result

def plot_expected_binned_noise (ax, num_bins):
	assert num_bins == EXP_NUM_BINS, "The number of bins that you provided is not the same as the number of bins in experimental data"
	x_cons_corner = calculate_exp_bounds(EXP_AVG_CONS, 0.9)
	x_cons_width = calculate_exp_bounds(EXP_AVG_CONS, 0.2)
	y_noise_corner = []
	y_noise_width = []
	for i in range(EXP_NUM_BINS):
		y_noise_corner.append(calculate_exp_bounds(EXP_NOISE[i], 0.9))
		y_noise_width.append(calculate_exp_bounds(EXP_NOISE[i], 0.2))
	#draw
	for i in range(len(EXP_NOISE)):
		for j in range(EXP_NUM_BINS):
			p = patches.Rectangle((x_cons_corner[j], y_noise_corner[i][j]), \
			x_cons_width[j], y_noise_width[i][j], color = noise_colors[i], alpha = 0.2)
			ax.add_patch(p)
	
	
def plot_binned_noise(out_fname, num_bins, avg_cons, in_noise, ex_noise):
	fig = plt.figure(figsize = (6,4), dpi=300)
	ax = fig.add_subplot(111)
	total_noise = [i + e for i,e in zip(in_noise, ex_noise)]
	plot_expected_binned_noise(ax, num_bins)
	tot_pt = ax.scatter(avg_cons, total_noise, s = 22, edgecolors = "none", c= noise_colors[0], label = "Total")
	in_pt = ax.scatter(avg_cons, in_noise, s = 22, edgecolors = "none", c=noise_colors[1], label = "Intrinsic")
	ex_pt = ax.scatter(avg_cons, ex_noise, s = 22, edgecolors = "none", c=noise_colors[2], label = "Extrinsic")
	
	plt.legend([tot_pt, in_pt, ex_pt], ["Total", "Intrinsic", "Extrinsic"], ncol = 3, \
	scatterpoints = 1, bbox_to_anchor = (0.5, 1.15), loc= 9, fontsize = 12)
	
	xmin = min([float(avg_cons[0] - 10), EXP_AVG_CONS[0] - 10])
	xmax = max([float(avg_cons[num_bins -1 ] + 10), float(EXP_AVG_CONS[EXP_NUM_BINS - 1] + 10)])
	
	ymin = min([min(in_noise) * 0.9, min(ex_noise) * 0.9, min(EXP_IN_NOISE) * 0.9, min(EXP_EX_NOISE) * 0.9])
	ymax = max([max(total_noise) * 1.1, max(EXP_TOT_NOISE) * 1.1])
	
	ax.tick_params(labelsize = 8)
	ax.xaxis.set_ticks(numpy.arange(0, xmax, determineTickInterval(xmax,10)))
	if xmax>=1: # need more space on the x-axis
		ax.set_xlim(0,xmax+xmin)
	else:
		ax.set_xlim(0,1)
	ax.yaxis.set_ticks(numpy.arange(0, ymax, .2), determineTickInterval(ymax, 10))
	ax.set_ylim(ymin - 0.05, ymax + 0.05)
	ax.set_ylabel("Noise (coefficient of variation squared)")
	ax.set_xlabel(r"Total $\textit{her}$ mRNA", fontsize = 12)
	# save figure
	fig.subplots_adjust(left=0.125, bottom=0.15, right=.95, top=.875, wspace=None, hspace=0.3)
	fig.savefig(out_fname, format = "jpg", dpi=300)
	
def run_stochastic(param_file, num_sets, out_dir):
	sto_command = [sto_sim, "-i", param_file, "-ns", str(num_sets), "-pf", "-od", out_dir]
	if 1 == call(sto_command):
		exit(1)

H1_AMP = [30,50]
H7_AMP = [40,55]

def plot_amplitude(out_fname, amplitude_h1, amplitude_h7):
	fig = plt.figure(figsize = (3,4), dpi=300)
	ax = fig.add_subplot(111)
	width = 0.4
	# calculate mean amplitude
	avg_amp_h1 = numpy.mean(amplitude_h1)
	avg_amp_h7 = numpy.mean(amplitude_h7)
	# draw bars based on simulation data
	ax.bar(0, avg_amp_h1, width, color = amp_colors[0], alpha = 0.4)
	ax.bar(1, avg_amp_h7, width, color = amp_colors[1], alpha = 0.4)
	# draw lines for the expected range 
	ax.plot([0 , 0 + width], [H1_AMP[0], H1_AMP[0]], color = black)
	ax.plot([0 , 0 + width], [H1_AMP[1], H1_AMP[1]], color = black)
	ax.plot([1, 1 + width], [H7_AMP[0], H7_AMP[0]], color = black)
	bounds = ax.plot([1, 1 + width], [H7_AMP[1], H7_AMP[1]], color = black)
	# setting ticks for the x axis
	ax.set_xticks([0+width/2, 1+width/2])
	ax.set_xticklabels((r'\textit{her1}',r'\textit{her7}'))
	ax.tick_params(axis='x', pad=5)
	ax.set_xlim(-0.3,1.6)
	# setting ticks for the y axis
	ymax = max([avg_amp_h1 + 5, avg_amp_h7 + 5, H1_AMP[1] + 5, H7_AMP[1] + 5])	
	ax.set_ylim(0, ymax)
	ax.set_yticks(numpy.arange(0, ymax, determineTickInterval(ymax, 10)))
	ax.set_ylabel(r'Spatial amplitude')
	updateYTicklabels(ax)
	bg_lg = mlines.Line2D([],[], color = '#000000')
	plt.legend([bg_lg], ['Expected bounds'], ncol = 1, \
	bbox_to_anchor = (0.5, 1.15), loc = 9, fontsize = 8)

	fig.subplots_adjust(left=0.2, bottom=0.075, right=0.95, top=.875,  wspace=None, hspace=0.2)	
	fig.savefig(out_fname, format = "jpg", dpi=300)

PER_BOUNDS = [22,38]
def plot_period(out_fname, period_h1, period_h7):
	fig = plt.figure(figsize = (3,4), dpi=300)
	ax = fig.add_subplot(111)
	width = 0.4
	# calculate mean amplitude
	avg_per_h1 = numpy.mean(period_h1)
	avg_per_h7 = numpy.mean(period_h7)
	# draw bars based on simulation data
	ax.bar(0, avg_per_h1, width, color = amp_colors[0], alpha = 0.4)
	ax.bar(1, avg_per_h7, width, color = amp_colors[1], alpha = 0.4)
	# draw lines for the expected range 
	ax.plot([0 , 0 + width], [PER_BOUNDS[0], PER_BOUNDS[0]], color = black)
	ax.plot([0 , 0 + width], [PER_BOUNDS[1], PER_BOUNDS[1]], color = black)
	ax.plot([1, 1 + width], [PER_BOUNDS[0], PER_BOUNDS[0]], color = black)
	bounds = ax.plot([1, 1 + width], [PER_BOUNDS[1], PER_BOUNDS[1]], color = black)
	# setting ticks for the x axis
	ax.set_xticks([0+width/2, 1+width/2])
	ax.set_xticklabels((r'\textit{her1}',r'\textit{her7}'))
	ax.tick_params(axis='x', pad=5)
	ax.set_xlim(-0.3,1.6)
	# setting ticks for the y axis
	ymax = max([avg_per_h1 + 5, avg_per_h7 + 5, PER_BOUNDS[1] + 5])	
	ax.set_ylim(0, ymax)
	ax.set_yticks(numpy.arange(0, ymax, determineTickInterval(ymax, 10)))
	ax.set_ylabel(r'Period')
	updateYTicklabels(ax)
	bg_lg = mlines.Line2D([],[], color = '#000000')
	plt.legend([bg_lg], ['Expected bounds'], ncol = 1, \
	bbox_to_anchor = (0.5, 1.15), loc = 9, fontsize = 8)

	fig.subplots_adjust(left=0.2, bottom=0.075, right=0.95, top=0.875,  wspace=None, hspace=0.2)	
	fig.savefig(out_fname, format = "jpg", dpi=300)

def draw_overlapping_cells_smooth(num_cells, save_fname, smooth_fname):
	command = ["/Applications/MATLAB_R2013a.app/bin/matlab", "-nodesktop", "-nosplash", "-nodisplay", \
		"-r", "plot_overlapping_cells(\'" + smooth_fname + "\', " + str(num_cells) + ", \'" + save_fname  + "\')"]
	if 1 == call(command):
		exit(1)
		
def main():
	"""
	param_100 = "../params_100_compare.txt"
	sto_out_dir = "../sto"
	# run stochastic simulation to get features
	num_set = 15
	#run_stochastic(param_100, num_set, sto_out_dir)
	for i in range(num_set):
		feature_filename = sto_out_dir + "/set_" + str(i) + "/wildtype/features.txt"
		num_cells, num_bins, amplitude_h1, amplitude_h7, period_h1, period_h7, \
		avg_cons, in_noise, ex_noise = read_features_data(feature_filename)
		plot_binned_noise(sto_out_dir + "/set_" + str(i) + "/wildtype/longNoise_her" + str(i) + ".jpg", num_bins, avg_cons, in_noise, ex_noise)
	"""
	feature_filename = "../../stochastic/stochastic/set_0/wildtype/features.txt"
	num_cells, num_bins, amplitude_h1, amplitude_h7, period_h1, period_h7, \
	avg_cons, in_noise, ex_noise = read_features_data(feature_filename)
	plot_binned_noise("../../stochastic/stochastic/set_0/wildtype/logNoise_her" + str(0) + ".jpg", 3, avg_cons, in_noise, ex_noise)
	plot_amplitude("../../stochastic/stochastic/set_0/wildtype/amplitude" + str(0) + ".jpg", amplitude_h1, amplitude_h7)
	plot_period("../../stochastic/stochastic/set_0/wildtype/period" + str(0) + ".jpg", period_h1, period_h7)
	draw_overlapping_cells_smooth(4, "../../stochastic/stochastic/set_0/wildtype/cells_smooth.fig", "../../stochastic/stochastic/set_0/wildtype/mHer1_smoothed.txt")
main()
