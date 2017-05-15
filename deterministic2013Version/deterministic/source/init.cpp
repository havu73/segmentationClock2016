#include <unistd.h> // Needed for getpid
#include <stdlib.h> // needed for atoi

#include "init.hpp" // Function declarations
#include "io.hpp"
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
			} else if (option_set(option, "-pc", "--print-cons")){ // Users only needs to specify this when they want to print any states other than the default MH1, MH7, MD
				ip.print_cons = true;
				i--;
			} else if (option_set(option, "-od", "--output-directory")){
				ensure_nonempty(option, value);
				store_filename(&(ip.out_dir), value);
				ip.has_out_dir = true; 
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
			}else if (option_set(option, "-rg", "--record-granularity")){
				ensure_nonempty(option, value);
				ip.record_gran = atoi(value);
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
	if (ip.record_gran <= 0){
		usage("Record_granularity needs to be at least 1.");
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


/* Read in the parameter set from a file/ a range file or through piping (method used by sres)
 * The functions that this function call can be found in io.cpp
 * Params: ip
 * 		   params_data: input_data for users' parameters file, declared in main.cpp
 * 		   pr: declared in main.cpp
 * 		   ranges_date: input_data for users' ranges file, declared in main.cpp 
 */
void read_sim_params(input_params& ip, parameters& pr, input_data& params_data, input_data& ranges_data){
	cout << term->blue;
	// OPTION 1: PIPING
	if (ip.piping) { 
		//cout << "Reading pipe " << term->reset << "(file descriptor " << ip.pipe_in << ") . . . ";
		read_pipe(pr, ip); 
		//term->done();
	// OPTION 2: PARAMETER FILE
	} else if (ip.read_params) { // If the user specified a parameter sets input file
		read_file(&params_data);
		for (int i = 0; i < ip.num_sets; i++) {
			if (params_data.index < params_data.size) { // Parse only as many lines as specified, even if the file is longe
				if(!parse_param_line(pr, i , params_data.buffer, params_data.index)) { // Parse each line until the file is empty or the required number of sets have been found
					ip.num_sets = i;
				}
			} else {
				ip.num_sets = i;
			}
		}
	// OPTION 3: RANGES FILE
	} else if (ip.read_ranges) { // If the user specified a ranges input file to generate random numbers from
		cout << "Generating " << term->reset << ip.num_sets << " random parameter sets according to the ranges in " << ranges_data.filename << " . . ." << endl;
		cout << "  ";
		read_file(&ranges_data);
		pair <double, double> ranges[NUM_PARAMS];
		parse_ranges_file(ranges, ranges_data.buffer); // get the ranges from users' ranges file
		srand(ip.pseed);
		for (int i = 0; i < ip.num_sets; i++) {
			for (int j = 0; j < NUM_PARAMS; j++) {
				pr.data[i][j] = random_double(ranges[j]); // Generate random numbers as paramters within the ranges
			}
		}
		term->done();
	} else {
		usage("Parameter must be piped in via -I or --pipe-in, read from a file via -i or --params-file, or generated from a ranges file and number of sets via -R or --ranges-file and -p or --parameter-sets, respectively.");
	} 
}

/* random_double generates a random double in the range specified by the given pair of doubles
	parameters:
		range: a pair of doubles that specifies the lower and upper bounds in that order
	returns: the random double
	notes:
	todo:
*/
double random_double (pair<double, double> range) {
	return range.first + (range.second - range.first) * rand() / (RAND_MAX + 1.0);
}
