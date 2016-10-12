#include <unistd.h> // Needed for getpid
#include <stdlib.h> // needed for atoi

#include "init.hpp" // Function declarations
#include "io.hpp"
#include "propensity.hpp"
#include "main.hpp"
#include "memory.hpp"
using namespace std;

terminal* term = NULL; // The global terminal struct

/* init_terminal creates and initializes a new terminal struct
	parameters:
	returns: nothing
	notes:
	todo:
*/
void init_terminal () {
	if (term != NULL) {
		delete term;
	}
	term = new terminal();
}

/* free_terminal frees the terminal from memory and resets the terminal text color to its default value
	parameters:
	returns: nothing
	notes:
	todo:
*/
void free_terminal () {
	cout << term->reset;
	delete term;
}

/* init_verbosity sets the verbose stream to /dev/null if verbose mode is not enabled
	parameters:
		ip: the program's input parameters
	returns: nothing
	notes:
	todo:
*/
void init_verbosity (input_params& ip) {
	if (!ip.verbose) {
		term->set_verbose_streambuf(ip.null_stream->rdbuf());
	}
}

void init_propensities(propensities& prop){
	prop.prop_funs[RPSH1] = &propensityRPSH1;
	prop.prop_funs[RPSH7] = &propensityRPSH7;
	prop.prop_funs[RPSD] = &propensityRPSD;
	
	prop.prop_funs[RPDH1] = &propensityRPDH1;
	prop.prop_funs[RPDH7] = &propensityRPDH7;
	prop.prop_funs[RPDD] = &propensityRPDD;
	prop.prop_funs[RPDH11] = &propensityRPDH11;
	prop.prop_funs[RPDH17] = &propensityRPDH17;
	prop.prop_funs[RPDH77] = &propensityRPDH77;
	
	prop.prop_funs[RDAH11] = &propensityRDAH11;
	prop.prop_funs[RDAH17] = &propensityRDAH17;
	prop.prop_funs[RDAH77] = &propensityRDAH77;
	
	prop.prop_funs[RDDH11] = &propensityRDDH11;
	prop.prop_funs[RDDH17] = &propensityRDDH17;
	prop.prop_funs[RDDH77] = &propensityRDDH77;
	
	prop.prop_funs[RMDH1] = &propensityRMDH1;
	prop.prop_funs[RMDH7] = &propensityRMDH7;
	prop.prop_funs[RMDD] = &propensityRMDD;
	
	prop.prop_funs[RMSH1] = &propensityRMSH1;
	prop.prop_funs[RMSH1N] = &propensityRMSH1N;
	
	prop.prop_funs[RAG1PH11] = &propensityRAG1PH11;
	prop.prop_funs[RDG1PH11] = &propensityRDG1PH11;
	
	prop.prop_funs[RAG1N] = &propensityRAG1N;
	prop.prop_funs[RDG1N] = &propensityRDG1N;
	
	prop.prop_funs[RMSH7] = &propensityRMSH7;
	prop.prop_funs[RMSH7N] = &propensityRMSH7N;
	
	prop.prop_funs[RAG7PH11] = &propensityRAG7PH11;
	prop.prop_funs[RDG7PH11] = &propensityRDG7PH11;
	
	prop.prop_funs[RAG7N] = &propensityRAG7N;
	prop.prop_funs[RDG7N] = &propensityRDG7N;
	
	prop.prop_funs[RMSD] = &propensityRMSD;
	prop.prop_funs[RAGDPH11] = &propensityRAGDPH11;
	prop.prop_funs[RDGDPH11] = &propensityRDGDPH11;
}

/* 
 * Put users input into a structures called input_params& ip that store data related to the
 * simulation
 */
void accept_input_params (int num_args, char** args, input_params& ip) {
	string o;
	string v;

	if (num_args >1){
		for (int i = 1; i < num_args; i += 2){
			
			//Process option and values. Changed compared to segmentation clock code because of some
			// weird bugs
			o = args[i];
			if (i < num_args - 1) {
				v = args[i + 1];
			} else {
				v = "";
			}
			
			char * option = new char [o.length() +1];
			char * value = new char [v.length() +1];
			strcpy(option, o.c_str());
			if (v.length() != 0){
				strcpy(value, v.c_str());
			} else {
				value = NULL; 
			}
			
			// Accept command-line arguments in both short and long form
			if (option_set(option, "-i", "--input-params-file")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.params_file), value);
				ip.read_params = true;
			} else if (option_set(option, "-r", "--ranges-file")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.ranges_file), value);
				ip.read_ranges = true;
			} else if (option_set(option, "-ps", "--print-states")){
				ensure_nonempty(option, value);
				get_print_states_data(ip, value);
				ip.print_cons = true;
			} else if (option_set(option, "-pc", "--print-cons")){ // Users only needs to specify this when they want to print any states other than the default MH1, MH7, MD
				ip.print_cons = true;
				i--;
			} else if (option_set(option, "-od", "--output-directory")){
				ensure_nonempty(option, value);
				store_filename(&(ip.out_dir), value);
			} else if (option_set(option, "-ns", "--num-parameter-sets")) {
				ensure_nonempty(option, value);
				ip.num_sets = atoi(value);
				if (ip.num_sets < 1) {
					usage("The number of parameters to run must be a positive integer. Set -p or --parameter-sets to at least 1.");
				}
			} else if (option_set(option, "-st", "--simulation-time")) {
				ensure_nonempty(option, value);
				ip.time_total = atoi(value);
				if (ip.time_total < 1) {
					usage("The simulation must be run for positive number of minutes. Set -m or --minutes to at least 1.");
				}
			} else if (option_set(option, "-s", "--seed")) {
				ensure_nonempty(option, value);
				ip.seed = atoi(value);
				if (ip.seed <= 0) {
					usage("The seed to generate random numbers must be a positive integer. Set -s or --seed to at least 1.");
				}
				ip.reset_seed = true;
			} else if (option_set(option, "-rs", "--reset-seed")) {
				ip.reset_seed = true;
				i--;
			} else if (option_set(option, "-paseed", "--parameters-seed")) {
				ensure_nonempty(option, value);
				ip.pseed = atoi(value);
				ip.store_pseed = true;
				if (ip.pseed <= 0) {
					usage("The seed to generate random parameters must be a positive integer. Set -d or --parameters-seed to at least 1.");
				}
			} else if (option_set(option, "-prseed", "--print-seeds")) {
				ensure_nonempty(option, value);
				store_filename(&(ip.seed_file), value);
				ip.print_seeds = true;
			} else if (option_set(option, "-m", "--mutants")) {
				ensure_nonempty(option, value);
				get_mutant_data(ip, value);
			} else if (option_set(option, "-pi", "--pipe-in")) {
				ensure_nonempty(option, value);
				ip.piping = true;
				ip.pipe_in = atoi(value);
				if (ip.pipe_in <= 0) {
					usage("The file descriptor to pipe data from must be a positive integer. Set -I or --pipe-in to be at least 1.");
				}
			} else if (option_set(option, "-po", "--pipe-out")) {
				ensure_nonempty(option, value);
				ip.piping = true;
				ip.pipe_out = atoi(value);
				if (ip.pipe_out <= 0) {
					usage("The file descriptor to pipe data into must be a positive integer. Set -O or --pipe-out to be at least 1.");
				}
			} else if (option_set(option, "-v", "--verbose")) {
				if (!ip.verbose) {
					ip.verbose = true;
				}
				i--;
			} else if (option_set(option, "-q", "--quiet")) {
				if (!ip.quiet) {
					ip.quiet = true;
					ip.cout_orig = cout.rdbuf();
					cout.rdbuf(ip.null_stream->rdbuf());
					term->set_verbose_streambuf(ip.null_stream->rdbuf());
				}
				i--;
			} else if (option_set(option, "-nc", "--no-color")) {
				free(term->blue);
				free(term->red);
				free(term->reset);
				strcpy(term->blue, "");
				strcpy(term->red, "");
				strcpy(term->reset, "");
				i--;
			} else if (option_set(option, "-h", "--help")) {
				usage("");
				i--;
			} else if (option_set(option, "-l", "--licensing")) { 
				licensing();
				i--;
			} else { // Do not ignore invalid arguments; exit with an error to let the user know this argument is problematic
				const char* message_0 = "'";
				const char* message_1 = "' is not a valid option! Please check that every argument matches one available in the following usage information.";
				char* message = (char*)mallocate(sizeof(char) * (strlen(message_0) + strlen(option) + strlen(message_1) + 1));
				sprintf(message, "%s%s%s", message_0, option, message_1);
				usage(message);
			}
			delete [] option;
			delete [] value;
		}
	} else {
		usage("Please provide parameters to the program");
	}
}

void ensure_nonempty (const char* option, const char* arg) {
	if (arg == NULL) {
		char* message = (char*)mallocate(strlen("Missing the argument for the '' option.") + strlen(option) + 1);
		sprintf(message, "Missing the argument for the '%s' option.", option);
		usage(message);
	}
}

inline bool option_set (const char* option, const char* short_name, const char* long_name) {
	return strcmp(option, short_name) == 0 || strcmp(option, long_name) == 0;
}

/* Process user's input to store the indices of the states into ip.print_states
 * ip.print_states is a list of indices of users want to keep track of over time
 * To see the indicies of states hard-coded in this program, go to macros.hpp
 * example users' input: -ps "0 1 2" means users want to keep track of MH1, MH7, MD
 * This function makes sure that indices of MH1, MH7, MD are default, and will not be 
 * included into ip.print_states, because we always include tracking MH1, MH7, MD in our 
 * concentrations structs
 */
void get_print_states_data(input_params& ip, char* input){
	char * result;
	int indices [NUM_STATES];
	int index = 0; 
	result = strtok(input, " \t\n"); // Each mutant is separated from each other by space
	while (result != NULL){
		if (!((strcmp(result, "0") == 0) || (strcmp(result, "1") == 0) || (strcmp(result, "2") == 0))){
			int state_index = atoi(result);
			if (state_index < 0 || state_index >= NUM_STATES){
				usage("Print state index is out of bound.");
			}
			indices[index] = state_index;
			index ++; 
		}
		result = strtok(NULL, " \t\n");
	}
	
	ip.num_print_states = index;
	ip.print_states = new int [index];
	for (int i = 0; i < index; i++){
		(ip.print_states)[i] = indices[i];
	}
}

/* Process user's input to store the mutants they want to simulate into ip.mutants
 * ip.mutants is a list of integers corresponding to indices of mutants
 * To see the indicies of mutants hard-coded in this program, go to macros.hpp
 * example users' input: -M "wt p1 p2", or -M "0 1 2 3"
 */
void get_mutant_data(input_params& ip, char* input){
	/*
	char * result;
	int index = 0; 
	result = strtok(input, " \t\n"); // Each mutant is separated from each other by space
	while (result != NULL){
		if ((strcmp(result, "wt") == 0) || (strcmp(result, "0") == 0)){
			ip.mutants[index] = WT;
		}
		
		else {
			usage("Mutant input is errornous");
		}
		
		index++;
		result = strtok(NULL, " \t\n");
	}
	ip.num_mutants = index;
	*/
}


/* Checking some parameters in ip
 * For security reasons for the program
 */
void check_input_params (input_params& ip){
	if (ip.piping && (ip.pipe_in == 0 || ip.pipe_out == 0)) {
		usage("If one end of a pipe is specified, the other must be as well. Set the file descriptors for both the pipe in (-I or --pipe-in) and the pipe out (-O or --pipe-out).");
	}
	if (!(ip.piping || ip.read_params || ip.read_ranges)) {
		usage("Parameter must be piped in via -I or --pipe-in, read from a file via -i or --params-file, or generated from a ranges file and number of sets via -R or --ranges-file and -p or --parameter-sets, respectively.");
	}
	if (ip.reset_seed) {
		init_seeds(ip, 0, false, false);
	}
}

/* Seeds is used in generating random parameter if users passed in ranges file.
 * Params: ip: structure to store users' input params
 * 		   set_num: index of the set currently being simulated (each set has its own seed)
 * 		   append: whether or not to indent in the file (not really important)
 * 		   indent_message: whether or not to indent message to print out in terminal (not important)
 */ 
void init_seeds (input_params& ip, int set_num, bool append, bool indent_message) {
	// If the file is being created then generate the parameter set seed
	if (!append && ip.store_pseed) {
		ip.pseed = generate_seed();
		if (indent_message) {
			term->verbose() << "  ";
		}
		term->verbose() << term->blue << "Using seed " << term->reset << ip.pseed << " for parameter set generation" << endl;
	}
	
	// If the seed was not specified by command-line argument or a new one should be generated, generate one
	if (ip.seed == 0 || append) {
		ip.seed = generate_seed();
		if (indent_message) {
			term->verbose() << "  ";
		}
		term->verbose() << term->blue << "Using seed " << term->reset << ip.seed << " for each run" << endl;
	}
	
	// If the user specified printing seeds to a file then print them to the specified file
	if (ip.print_seeds) {
		ofstream seed_file;
		term->verbose() << "  ";
		open_file(&seed_file, ip.seed_file, append);
		if (!append && ip.store_pseed) {
			seed_file << "pseed: " << ip.pseed << endl;
		}
		seed_file << "seed " << set_num << ": " << ip.seed << endl;
	}
}

/* generate_seed generates a seed based on the current UNIX time and process ID (PID) of the program instance
	parameters:
	returns: the generated seed
	notes:
		Generating a seed based on both the time and PID ensures concurrently executed versions of this program do not generated the same seed.
	todo:
*/
int generate_seed () {
	return abs(((time(0) * 181) * ((getpid() - 83) * 359)) % 805306457);
}

void calculate_max_cond_score (input_params& ip){}
