#include "main.hpp" // Function declarations
#include "structs.hpp"
#include "init.hpp"
#include "sim.hpp"
#include "debug.hpp"

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
	test_embryo_neighbor_network(ip);
	/*
	init_terminal();
	accept_input_params (argc, argv, ip);
	check_input_params(ip);
	init_verbosity(ip);
	
	// Calculate the maximum conditional score that we can get
	if (ip.num_mutants != NUM_MUTANTS){ // Since the default is to run all mutants available, and ip.max_cond_score is default MAX_SCORE, we only need to calculate score when users specify something different. 
										// Dont worry about cases when users enter the same number of mutants but they are wrongly entered, because that gets checked in accept_input_params
		calculate_max_cond_score(ip);
	}
	
	// declare input_data objects based on users' input about input file names. The buffer of these objects is empty now, and will be filled in right after this declaration section
	input_data params_data(ip.params_file);
	input_data ranges_data(ip.ranges_file);
	*/
	
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
	cout << "-i,  --input-params-file  [filename]   : the relative filename of the parameter sets input file, default=none" << endl;
	cout << "-r,  --ranges-file        [filename]   : the relative filename of the parameter ranges input file, default=none" << endl;
	cout << "-ps, --print-states	   [string]		: the space separated list of state indices that users want to keep track of over time. MH1, MH7, MD are default. Ex: -ps \"4 5 6\"" << endl;
	cout << "-pc, --print-cons         [N/A]        : print concentration values to the specified output directory, default=false" << endl;
	cout << "-od, --output-directory   [directory]  : the relative directory where output of the program is saved, default=none" << endl;
	cout << "-ns, --num-parameter-sets [int]        : the number of parameters for which to simulate the model, min=1, default=1" << endl;
	cout << "-st, --total-time         [int]        : the number of minutes to simulate before ending, min=1, default=600" << endl;
	cout << "-s,  --seed               [int]        : the seed to generate random numbers, min=1, default=generated from the time and process ID" << endl;
	cout << "-rs, --reset-seed         [N/A]        : reset the seed after each parameter set so the initial seed is used each time, default=unused" << endl;
	cout << "-paseed, --parameters-seed[int]        : the seed to generate random parameter sets, min=1, default=generated from the time and process ID" << endl;
	cout << "-prseed, --print-seeds    [filename]   : the relative filename of the seed output file, default=none" << endl;
	cout << "-m,  --mutants            [int]        : the space-separated list of mutants indices to run for each parameter set, default= list of all possible mutants" << endl;
	cout << "-pi, --pipe-in            [file desc.] : the file descriptor to pipe data from (usually passed by the sampler), default=none" << endl;
	cout << "-po, --pipe-out           [file desc.] : the file descriptor to pipe data into (usually passed by the sampler), default=none" << endl;
	cout << "-v,  --verbose            [N/A]        : print detailed messages about the program and simulation state, default=unused" << endl;
	cout << "-q,  --quiet              [N/A]        : hide the terminal output, default=unused" << endl;
	cout << "-nc, --no-color           [N/A]        : disable coloring the terminal output, default=unused" << endl;
	cout << "-l,  --licensing          [N/A]        : view licensing information (no simulations will be run)" << endl;
	cout << "-h,  --help               [N/A]        : view usage information (i.e. this)" << endl;
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
	cout << "Stochastic Simulation for zebrafish segmentation" << endl;
	cout << "Copyright (C) 2016 -----" << endl;
	cout << "This program comes with ABSOLUTELY NO WARRANTY" << endl;
	cout << "This is free software, and you are welcome to redistribute it under certain conditions;" << endl;
	cout << "You can use this code and modify it as you wish under the condition that you refer to the article: \"Short-lived Her proteins drive robust synchronized oscillations in the zebrafish segmentation clock\" (Development 2013 140:3244-3253; doi:10.1242/dev.093278)" << endl;
	cout << endl;
	exit(EXIT_SUCCESS);
}
