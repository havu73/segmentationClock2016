#include <time.h>
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
	//clock_t t;
	//t = clock();
	
	input_params ip;
	
	init_terminal();
	accept_input_params (argc, argv, ip);
	check_input_params(ip);
	
	init_verbosity(ip);
	
	
	// declare input_data objects based on users' input about input file names. The buffer of these objects is empty now, and will be filled in right after this declaration section
	input_data params_data(ip.params_file);
	input_data ranges_data(ip.ranges_file);
	
	// process parameters
	parameters pr(ip.num_sets);
	read_sim_params(ip, pr, params_data, ranges_data);
	
	simulate_all_params_sets(ip, pr);
	
	//t = clock() - t; 
	//cout << "Execution time: " << ((float) t) / CLOCKS_PER_SEC << endl;
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
	cout << "-i,  --input-params-file  	[filename]   : the relative filename of the parameter sets input file, default=none" << endl;
	cout << "-r,  --ranges-file        	[filename]   : the relative filename of the parameter ranges input file, default=none" << endl;
	cout << "-pc, --print-cons         	[N/A]        : print concentration values to the specified output directory, default=false" << endl;
	cout << "-od, --output-directory   	[directory]  : the relative directory where output of the program is saved, default=none" << endl;
	cout << "-ns, --num-parameter-sets 	[int]        : the number of parameters for which to simulate the model, min=1, default=1" << endl;
	cout << "-st, --total-time         	[int]        : the number of minutes to simulate before ending, min=1, default=600" << endl;
	cout << "-s,  --seed               	[int]        : the seed to generate random numbers, min=1, default=generated from the time and process ID" << endl;
	cout << "-rs, --reset-seed         	[N/A]        : reset the seed after each parameter set so the initial seed is used each time, default=unused" << endl;
	cout << "-paseed, --parameters-seed	[int]        : the seed to generate random parameter sets, min=1, default=generated from the time and process ID" << endl;
	cout << "-prseed, --print-seeds    	[filename]   : the relative filename of the seed output file, default=none" << endl;
	cout << "-pi, --pipe-in            	[file desc.] : the file descriptor to pipe data from (usually passed by the sampler), default=none" << endl;
	cout << "-po, --pipe-out           	[file desc.] : the file descriptor to pipe data into (usually passed by the sampler), default=none" << endl;
	cout << "-v,  --verbose            	[N/A]        : print detailed messages about the program and simulation state, default=unused" << endl;
	cout << "-q,  --quiet              	[N/A]        : hide the terminal output, default=unused" << endl;
	cout << "-rg, --record-granularity	[int]		 : number of reactions fired in between two times that we record the concentrations of states we care about. Default = 30" << endl;
	cout << "We require that check-done-granularity and record-granularity are both at least 1, and that check-done-granularity is invisible by record-granularity" << endl;
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
