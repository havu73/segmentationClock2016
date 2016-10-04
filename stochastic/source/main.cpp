#include "main.hpp" // Function declarations
#include "structs.hpp"
#include "init.hpp"
#include "sim.hpp"

using namespace std;

extern terminal* term; // Declared in init.cpp

/* main is called when the program is run and performs all program functionality
	parameters:
		argc: the number of command-line arguments
		argv: the array of command-line arguments
	returns: 0 on success, a positive integer on failure
	notes:
		Main should only delegate functionality; let the functions it calls handle specific tasks. This keeps the function looking clean and helps maintain the program structure.
	todo:
*/
int main(int argc, char** argv) {
	input_params ip;
	ip.print_states = new int[2];
	(ip.print_states)[0] = PH1;
	(ip.print_states)[1] = PH7;
	ip.num_print_states = 2;
	
	concentrations cons(ip);
	
	for (int i = 0; i < 4; i++){
		(*(cons.cons_data)[MH1])[i] = i;
		(*(cons.cons_data)[MH7])[i] = i;
		(*(cons.cons_data)[MD])[i] = i;
		(*(cons.cons_data)[PH1])[i] = i;
		(*(cons.cons_data)[PH1])[i] = i;
	}
	
	for (int i = 1; i < 3; i++){
		cout << (*(cons.cons_data)[MH1])[i] << endl;
		cout << (*(cons.cons_data)[MH7])[i] << endl;
		cout << (*(cons.cons_data)[MD])[i] << endl;
		cout << (*(cons.cons_data)[PH1])[i] << endl;
		cout << (*(cons.cons_data)[PH1])[i] << endl;
	}
	cons.reset();
	for (int i = 1; i < 3; i++){
		cout << (*(cons.cons_data)[MH1])[i] << endl;
		cout << (*(cons.cons_data)[MH7])[i] << endl;
		cout << (*(cons.cons_data)[MD])[i] << endl;
		cout << (*(cons.cons_data)[PH1])[i] << endl;
		cout << (*(cons.cons_data)[PH1])[i] << endl;
	}
}

/* usage prints the usage information and, optionally, an error message and then exits
	parameters:
		message: an error message to print before the usage information (set message to NULL or "\0" to not print any error)
	returns: nothing
	notes:
		This function exits after printing the usage information.
		Note that accept_input_params in init.cpp handles actual command-line input and that this information should be updated according to that function.
	todo:
		TODO somehow free memory even with the abrupt exit
*/
void usage (const char* message) {
	cout << endl;
	bool error = message != NULL && message[0] != '\0';
	if (error) {
		cout << term->red << message << term->reset << endl << endl;
	}
	cout << "Usage: [-option [value]]. . . [--option [value]]. . ." << endl;
	cout << "-i, --params-file        [filename]   : the relative filename of the parameter sets input file, default=none" << endl;
	cout << "-R, --ranges-file        [filename]   : the relative filename of the parameter ranges input file, default=none" << endl;
	cout << "-u, --perturb-file       [filename]   : the relative filename of the perturbations input file, default=none" << endl;
	cout << "-r, --gradients-file     [filename]   : the relative filename of the gradients input file, default=none" << endl;
	cout << "-o, --print-passed       [filename]   : the relative filename of the passed sets output file, default=none" << endl;
	cout << "-t, --print-cons         [N/A]        : print concentration values to the specified output directory, default=unused" << endl;
	cout << "-B, --binary-cons-output [N/A]        : print concentration values as binary numbers rather than ASCII, default=unused" << endl;
	cout << "-f, --print-osc-features [filename]   : the relative filename of the file summarizing all the oscillation features, default=none" << endl;
	cout << "-D, --directory-path     [directory]  : the relative directory where concentrations or anterior oscillation features files will be printed, default=none" << endl;
	cout << "-A, --anterior-feats     [N/A]        : print in depth oscillation features for the anterior cells over time, default=unused" << endl;
	cout << "-P, --posterior-feats    [N/A]        : print in depth oscillation features for the posterior cells over time, default=unused" << endl;
	cout << "-W, --print-conditions   [filename]   : the relative filename of the passed and failed conditions file, default=none" << endl;
	cout << "-E, --print-scores       [filename]   : the relative filename of the mutant scores file, default=none" << endl;
	cout << "-L, --print-cells        [int]        : the number of columns of cells to print for plotting of single cells on top of each other, min=0, default=0" << endl;
	cout << "-b, --big-granularity    [int]        : the granularity in time steps with which to store data, min=1, default=1" << endl;
	cout << "-g, --small-granularity  [int]        : the granularity in time steps with which to simulate data, min=1, default=1" << endl;
	cout << "-x, --total-width        [int]        : the tissue width in cells, min=3, default=3" << endl;
	cout << "-w, --initial-width      [int]        : the tissue width in cells before anterior growth, min=3, max=total width, default=3" << endl;
	cout << "-y, --height             [int]        : the tissue height in cells, min=1, default=1" << endl;
	cout << "-S, --step-size          [float]      : the size of the timestep to be used for solving the DDEs using Euler's method, default=0.01" << endl;
	cout << "-m, --total-time         [int]        : the number of minutes to simulate before ending, min=1, default=1200" << endl;
	cout << "-T, --split-time         [int]        : the number of minutes it takes for cells to split, min=1, default=6, 0=never" << endl;
	cout << "-G, --time-til-growth    [int]        : the number of minutes to wait before allowing cell growth, min=0, default=600" << endl;
	cout << "-p, --parameter-sets     [int]        : the number of parameters for which to simulate the model, min=1, default=1" << endl;
	cout << "-s, --seed               [int]        : the seed to generate random numbers, min=1, default=generated from the time and process ID" << endl;
	cout << "-X, --reset-seed         [N/A]        : reset the seed after each parameter set so the initial seed is used each time, default=unused" << endl;
	cout << "-d, --parameters-seed    [int]        : the seed to generate random parameter sets, min=1, default=generated from the time and process ID" << endl;
	cout << "-e, --print-seeds        [filename]   : the relative filename of the seed output file, default=none" << endl;
	cout << "-a, --max-con-threshold  [float]      : the concentration threshold at which to fail the simulation, min=1, default=infinity" << endl;
	cout << "-C, --short-circuit      [N/A]        : stop simulating a parameter set after a mutant fails, default=unused" << endl;
	cout << "-M, --mutants            [int]        : the number of mutants to run for each parameter set, min=1, max=" << ", default=" << endl;
	cout << "-I, --pipe-in            [file desc.] : the file descriptor to pipe data from (usually passed by the sampler), default=none" << endl;
	cout << "-O, --pipe-out           [file desc.] : the file descriptor to pipe data into (usually passed by the sampler), default=none" << endl;
	cout << "-c, --no-color           [N/A]        : disable coloring the terminal output, default=unused" << endl;
	cout << "-v, --verbose            [N/A]        : print detailed messages about the program and simulation state, default=unused" << endl;
	cout << "-q, --quiet              [N/A]        : hide the terminal output, default=unused" << endl;
	cout << "-l, --licensing          [N/A]        : view licensing information (no simulations will be run)" << endl;
	cout << "-h, --help               [N/A]        : view usage information (i.e. this)" << endl;
	cout << endl << term->blue << "Example: ./simulation -i parameters.csv --parameters 10 -m 2000 --no-color" << term->reset << endl << endl;
	if (error) {
		exit(EXIT_INPUT_ERROR);
	} else {
		exit(EXIT_SUCCESS);
	}
}

/* licensing prints the program's copyright and licensing information and then exits
	parameters:
	returns: nothing
	notes:
	todo:
*/
void licensing () {
	cout << endl;
	cout << "Simulation for zebrafish segmentation" << endl;
	cout << "Copyright (C) 2016 -----" << endl;
	cout << "This program comes with ABSOLUTELY NO WARRANTY" << endl;
	cout << "This is free software, and you are welcome to redistribute it under certain conditions;" << endl;
	cout << "You can use this code and modify it as you wish under the condition that you refer to the article: \"Short-lived Her proteins drive robust synchronized oscillations in the zebrafish segmentation clock\" (Development 2013 140:3244-3253; doi:10.1242/dev.093278)" << endl;
	cout << endl;
	exit(EXIT_SUCCESS);
}
