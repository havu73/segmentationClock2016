#ifndef STRUCTS_HPP
#define STRUCTS_HPP
///INCLUDE HEADER FILES
#include <cstdlib> // Needed for cmath
#include <cstring> // Needed for strlen, memset, memcpy
#include <iostream> // Needed for cout
#include <bitset> // Needed for bitset
#include <fstream> // Needed for ofstream
#include <string>
#include <map> // Needed for map
#include <vector>
#include "macros.hpp"
#include "memory.hpp"


//use standard name_space here
using namespace std;


struct terminal {
	// Escape codes
	const char* code_blue;
	const char* code_red;
	const char* code_reset;
	
	// Colors
	char * blue;
	char * red;
	char * reset;
	
	// Verbose stream
	streambuf* verbose_streambuf;
	ostream* verbose_stream;
	
	terminal () {
		this->code_blue = "\x1b[34m";
		this->code_red = "\x1b[31m";
		this->code_reset = "\x1b[0m";
		this->blue = new char [10];
		this->red = new char [10];
		this->reset = new char [10];
		strcpy(this->blue, this->code_blue);
		strcpy(this->red, this->code_red);
		strcpy(this->reset, this->code_reset);
		this->verbose_stream = new ostream(cout.rdbuf());
	}
	
	~terminal () {
		delete []this->blue;
		delete []this->red;
		delete []this->reset;
		delete verbose_stream;
	}
	
	// Indicates a task is done (pass terminal->verbose() into this function to print only with verbose mode on)
	void done (ostream& stream) {
		stream << this->blue << "Done" << this->reset << endl;
	}
	
	// Indicates a task is done
	void done () {
		done(cout);
	}
	
	// Indicates the program is out of memory
	void no_memory () {
		cout << this->red << "Not enough memory!" << this->reset << endl;
	}
	
	// Indicates the program couldn't read from a pipe
	void failed_pipe_read () {
		cout << this->red << "Couldn't read from the pipe!" << this->reset << endl;
	}
	
	// Indicates the program couldn't write to a pipeinput_params
	void failed_pipe_write () {
		cout << this->red << "Couldn't write to the pipe!" << this->reset << endl;
	}
	
	// Returns the verbose stream that prints only when verbose mode is on
	ostream& verbose () {
		return *(this->verbose_stream);
	}
	
	// Sets the stream buffer for verbose mode
	void set_verbose_streambuf (streambuf* sb) {
		this->verbose_stream->rdbuf(sb);
		this->verbose_streambuf = this->verbose_stream->rdbuf();
	}
};

struct input_params {
	// Input and output files' paths and names (either absolute or relative)
	char* out_dir; // the path of the output directory
	bool has_out_dir; // whether user specified an output directory
	char* params_file; // The path and name of the parameter sets file, default=input.params
	bool read_params; // Whether or not the read the parameter sets file, default=false
	char* ranges_file; // The path and name of the parameter ranges file, default=none
	bool read_ranges; // Whether or not to read the ranges file, default=false
	char* dir_path; // The path of the output directory for concentrations or oscillation features, default=none
	bool print_cons; // Whether or not to print concentrations, default=false
	bool binary_cons_output; // Whether or not to print the binary or ASCII value of numbers in the concentrations output files
	
	// Sets
	int num_sets; // The number of parameter sets to simulate, default=1
	

	// Time
	int time_total; // The number of minutes to run each simulation for, default=1200
	double step_size; // The time step in minutes used for Euler's method, default=0.01
	int record_gran; // The granularity in time steps with which to store data, default=1
	
	//Seed
	int seed; // The seed, used for generating random numbers, default=generated from the time and process ID
	bool reset_seed; // Whether or not to reset the seed after each parameter set, default=false
	int pseed; // The seed, used for generating random parameter sets, default=generated from the time and process ID
	bool store_pseed; // Whether or not to store the parameter generation seed, pseed, default=false
	
	//mutant management
	bool short_circuit; // Whether or not to stop simulating a parameter set after a mutant fails
	int num_active_mutants;
	
	// Piping data
	bool piping; // Whether or not input and output should be piped (as opposed to written to disk), default=false
	int pipe_in; // The file descriptor to pipe data from, default=none (0)
	int pipe_out; // The file descriptor to pipe data into, default=none (0)
	
	// Output stream data
	bool verbose; // Whether or not the program is verbose, i.e. prints many messages about program and simulation state, default=false
	bool quiet; // Whether or not the program is quiet, i.e. redirects cout to /dev/null, default=false
	streambuf* cout_orig; // cout's original buffer to be restored at program completion
	ofstream* null_stream; // A stream to /dev/null that cout is redirected to if quiet mode is set
	
	input_params () {
		//input and output files
		this->out_dir = new char [30];
		this->has_out_dir = false;
		this->params_file = new char[30];
		this->read_params = false;
		this->ranges_file = new char [30];
		this->read_ranges = false;
		this->dir_path = new char [30];
		this->print_cons = false;
		this->binary_cons_output = false;
		
		//timing
		this->num_sets = 1;
		this->step_size = 0.01;
		this->time_total = 600;
		this->record_gran = 100;
		
		//seeds
		this->seed = 0;
		this->reset_seed = false;
		this->pseed = 0;
		this->store_pseed = false;
		
		//mutant management
		this->short_circuit = false;
		this->num_active_mutants = 1;
		
		//piping data
		this->piping = false;
		this->pipe_in = 0;
		this->pipe_out = 0;
		
		//output stream data
		this->verbose = false;
		this->quiet = false;
		this->cout_orig = NULL;
		this->null_stream = new ofstream("/dev/null");
	}
	
	~input_params () {
		delete [] this->params_file;
		delete [] this->dir_path;
		delete this->null_stream;
	}
};

struct rates {
	//rates bases and rates for mutants
	double* data;  // Base rates taken from the current parameter set
	
	explicit rates () {
		this->data = new double [NUM_RATES];
		memset(this->data, 0, sizeof(double) * NUM_RATES);
	}
	
	void clear (){
		delete [](this->data);
	}
	
	void reset(){
		for (int i = 0; i < NUM_RATES; i ++){
			this->data[i] = 0;
		}
	}
	
	~rates () {
		delete [](this->data);
	}
};

struct con_levels {
	int num_sim_steps;
	double ** sim_data; //[states][timestep]
	
	con_levels (int num_sim_steps) {
		this->num_sim_steps = num_sim_steps;
		this->sim_data = new double* [NUM_STATES];
		for (int i = 0; i < NUM_STATES; i ++){
			(this->sim_data)[i] = new double [this->num_sim_steps];
		}
		this->reset();
	}
	
	// Sets every value in the struct to 0, except genes should be 2, but does not free any memory
	void reset () {
		for (int i = 0 ; i < NUM_STATES; i++){
			for (int j = 0; j < this->num_sim_steps; j++){
				(this->sim_data)[i][j] = 0;
			}
		}
	}
	
	// Frees the memory used by the struct
	void clear () {
		for (int i = 0; i < NUM_STATES; i++){
			delete [] (this->sim_data)[i];
		}
		delete [] this->sim_data;
	}
	
	~con_levels () {
		this->clear();
	}
};

struct input_data {
	char* filename; // The path and name of the file
	char* buffer; // A buffer to store the file's contents
	int size; // The number of bytes the file's contents take up
	int index; // The current index to access the buffer from
	
	explicit input_data (char* filename) {
		this->filename = filename;
		this->buffer = NULL;
		this->size = 0;
		this->index = 0;
	}
	void clear(){
		free(this->buffer);
	}
	
	~input_data () {
		free(this->buffer);
	}
};


struct parameters {
	//rates bases and rates for mutants
	double** data;  // Base rates taken from the current parameter set
	int num_sets; 
	explicit parameters (int num_sets) {
		this->num_sets = num_sets;
		this->data = new double* [num_sets];
		for (int i = 0; i < num_sets; i ++){
			this->data[i] = new double[NUM_PARAMS];
			memset(this->data[i], 0, sizeof(double) * NUM_PARAMS);
		}
	}
	
	void clear (){
		for (int i = 0; i < this->num_sets; i++) {
				delete[] this->data[i];
		}
		delete[] this->data;
	}
	
	~parameters () {
		this->clear();
	}
};

struct growin_array {
	int* array; // The array of integers
	int size; // The size of the array

	growin_array() {}

	void initialize(int size) {
		array = new int[size];
		this->size = size;
	}
	
	explicit growin_array (int size) {
		this->initialize(size);
	}

	int& operator[] (int index) {
		if (index >= this->size) {
			int* new_array = new int[2 * this->size];
			for (int i = 0; i < this->size; i++) {
				new_array[i] = array[i];
			}
			delete[] array;
			this->size = 2 * this->size;
			array = new_array;
		}
		return array[index];
	}
	
	void reset(int new_size) {
		this->size = new_size;
	}
	
	int get_size() {
		return this->size;
	}
	
	~growin_array () {
		delete[] array;
	}
};

struct peak_trough {
	std::vector<int> peaks;
	std::vector<int> troughs;
	
	void initialize(){
		this->peaks.clear();
		this->troughs.clear();
	}
	
	explicit peak_trough(){
		this->initialize();
	}
	
	void reset(){
		initialize();
	}
};
#endif
