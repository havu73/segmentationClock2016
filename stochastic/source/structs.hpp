#ifndef STRUCTS_HPP
#define STRUCTS_HPP
///INCLUDE HEADER FILES
#include <cstdlib> // Needed for cmath
#include <cstring> // Needed for strlen, memset, memcpy
#include <iostream> // Needed for cout
#include <bitset> // Needed for bitset
#include <fstream> // Needed for ofstream
#include <string>
#include <vector> // Needed for concentration levels
#include <queue> // Needed for priority queue inside complete_delay
#include "macros.hpp"

//use standard name_space here
using namespace std;

/*
 * This struct contains repeated used phrases and colors to print onto terminal in our program
 */
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

struct input_params{
	// Input and output files' paths and names (either absolute or relative)
	char* params_file; // The path and name of the parameter sets file, default=input.params
	bool read_params; // Whether or not the read the parameter sets file, default=false
	char* gradients_file; // The path and name of the gradients file, default=none
	bool read_gradients; // Whether or not to read the gradients file, default=false
	char* ranges_file; // The path and name of the parameter ranges file, default=none
	bool read_ranges; // Whether or not to read the ranges file, default=false
	char* perturb_file; // The path and name of the perturbation file, default=none
    bool read_perturb; // Whether or not to read the perturbation file, default=false
	
	// What states whose concentrations to store
	int* print_states; // list of indices of states that we want to store concentrations over time for plotting and testing
	int num_print_states; // number of states to print the concentrations over time
	bool print_cons; // whther or not program is requested to print concentrations over time. False when the size of this->print_states is 0 
	char* out_dir; // The path of the output directory for concentrations or oscillation features, default=none	
	
	// Sets
	int num_sets; // The number of parameter sets to simulate, default=1
	
	// Timing parameters
	int time_total; // The number of minutes to run each simulation for, default=600
	
	// Parameters about seeds
	unsigned int seed; // The seed, used for generating random numbers, default=generated from the time and process ID
	bool reset_seed; // Whether or not to reset the seed after each parameter set, default=false
	int pseed; // The seed, used for generating random parameter sets, default=generated from the time and process ID
	bool store_pseed; // Whether or not to store the parameter generation seed, pseed, default=false
	char* seed_file; // Default=none
	bool print_seeds; // Whether or not to print the seeds used to the seed file
	
	// Mutant parameters
	int * mutants; // List of index of mutants that users want to simulate and test, always need to have wildtype in it
	
	// Piping data
	bool piping; // Whether or not input and output should be piped (as opposed to written to disk), default=false
	int pipe_in; // The file descriptor to pipe data from, default=none (0)
	int pipe_out; // The file descriptor to pipe data into, default=none (0)
	
	// Output stream data
	bool verbose; // Whether or not the program is verbose, i.e. prints many messages about program and simulation state, default=false
	bool quiet; // Whether or not the program is quiet, i.e. redirects cout to /dev/null, default=false
	streambuf* cout_orig; // cout's original buffer to be restored at program completion
	ofstream* null_stream; // A stream to /dev/null that cout is redirected to if quiet mode is set
	
	input_params(){
		// IO files
		this->params_file = new char[30];
		this->read_params = false;
		this->ranges_file = new char [30];
		this->read_ranges = false;		
		this->gradients_file = new char[30];
		this->read_gradients = false;
		this->perturb_file = new char [30];
		this->read_perturb = false; 
		
		// printing concentrations over time
		this->print_states = 0; // NULL, need to be initialzied using new int[size] later if users want to print
		this->num_print_states = 0;
		this->print_cons = false;
		this->out_dir = new char [30];
		
		//sets
		this->num_sets = 1;
		
		// timing
		this->time_total = 600;
		
		//seeds
		this->seed = 0;
		this->reset_seed = false;
		this->pseed = 0;
		this->store_pseed = false;
		this->seed_file = new char [30];
		this->print_seeds = false;		
		
		//mutants
		this->mutants = 0; // NULL
		
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
	
	~input_params(){
		delete [] this->params_file;
		delete [] this->ranges_file;
		delete [] this->perturb_file;
		delete [] this->gradients_file;
		delete [] this->mutants;
		delete [] this->seed_file;
		delete this->null_stream;
		delete [] this->out_dir;
		if (this->print_states != 0){
			delete [] this->print_states;
		}
	}
};

struct concentrations{
	bool initialized; // Whether or not this struct's data have been initialized
	vector<int>** cons_data; // a list of vectors (dynamically growing lists) containing the concentrations of different states over time. 
							 // Most vectors have length 1, aka we do not care about its progression over time. 
						     // For some states that users require to keep tract of, we will keep record of all cahnges in concentrations  over time
	vector<double> ** time_data; 	// a list of (#print_states) vectors recording at what time is the concentration of each states that we keep record of
	
	concentrations(){
		this->initialized = false;
	}
	
	concentrations(input_params& ip){
		this->initialize(ip);
	}
	
	void initialize(input_params& ip){
		if (this->initialized){
			this->reset();
		}
		else{
			this->cons_data = new vector<int>* [NUM_STATES];
			this->time_data = new vector<double>* [NUM_KEEP_STATES + ip.num_print_states];
			int time_index = 0;
			int keep_index = 0;
			for (int i = 0; i < NUM_STATES; i++){
				if (i == KEEPMH1 || i == KEEPMH7 || i == KEEPMD){
					(this->cons_data)[i] = new vector<int>(40000, 0); // we want to keep tract 40000
					(this->time_data)[time_index] = new vector<double>(40000, 0);
					time_index ++; 
				}
				else if (ip.num_print_states > 0 && keep_index < ip.num_print_states && i == ip.print_states[keep_index]){
					(this->cons_data)[i] = new vector<int>(40000,0);
					(this->time_data)[time_index] = new vector<double>(40000,0);
					time_index ++;
					keep_index ++;
				}
				else{
					(this->cons_data)[i] = new vector<int>(1,0);
				}
			}
			this->initialized = true;
		}
	}
	
	void reset(){
		for (int i = 0; i < NUM_STATES; i++){
			memset(&((*((this->cons_data)[i]))[0]), 0, sizeof(int) * ((this->cons_data)[i])->size());
		}
		
		int time_size = sizeof(this->time_data) / sizeof((this->time_data)[0]);
		for (int i = 0; i < time_size; i++){
			memset(&((*((this->time_data)[i]))[0]), 0, sizeof(int) * ((this->time_data)[i])->size());
		}
	}
	
	~concentrations(){
		for (int i = 0; i < NUM_STATES; i++){
			delete (this->cons_data)[i];
		}
		int time_size = sizeof(this->time_data) / sizeof((this->time_data)[0]);
		for (int i = 0; i < time_size; i++){
			delete (this->time_data)[i];
		}
	}
};

struct delay_reaction{
	int reaction_index;
	double complete_time;
	delay_reaction(int index, double time){
		this->reaction_index = index;
		this->complete_time = time;
	}
};

struct compare_sooner_completion{
	bool operator()(delay_reaction& lhs, delay_reaction& rhs){
		return (lhs.complete_time >= rhs.complete_time);
	}
};

struct complete_delay{
	priority_queue<delay_reaction, vector<delay_reaction>, compare_sooner_completion>* pq;
	complete_delay(){
		this->pq = new priority_queue<delay_reaction, vector<delay_reaction>, compare_sooner_completion>;
	}
	
	~complete_delay(){
		delete this->pq;
	}
	
	bool is_empty(){
		return (this->pq)->empty();
	}
	
	// return the time of the soonest complete delay reaction
	// need to ensure that the queue is not empty before calling this function
	double see_soonest(){ 
		return ((this->pq)->top()).complete_time;
	}
	
	void complete_soonest(){
		(this->pq)->pop();
	}
	
	void initiate_delay(int function_index, double next_complete){
		delay_reaction dr (function_index, next_complete);
		(this->pq)->push(dr);
	}
};

struct rates {
	//rates bases and rates for mutants
	double* rates_base;  // Base rates taken from the current parameter set
	
	explicit rates () {
		this->rates_base = new double [NUM_RATES];
		memset(this->rates_base, 0, sizeof(double) * NUM_RATES);
	}
	
	void clear (){
		delete [](this->rates_base);
	}
	
	void reset(){
		memset(this->rates_base, 0, sizeof(double) * NUM_RATES);
	}
	
	~rates () {
		delete [](this->rates_base);
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
			this->data[i] = new double[NUM_RATES];
			memset(this->data[i], 0, sizeof(double) * NUM_RATES);
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


struct cell{
	concentrations * cons;		// concentrations of different states inside the cell
	double * propen;		// propensities of reactions inside the cell
	double * next_internal;		// next internal time of reactions in the system
	double * current_internal; // current internal time of reactions in the system
	complete_delay * cdelay; // the priority queue of the complettion time of delay reactions going on at a specific moment in the cell
	double absolute_time; // the current absolute time of the cell
	int index; 	//the index of the cell inside the embryo
	
	cell(int index, input_params& ip){
		this->cons = new concentrations(ip);
		this->propen = new double[NUM_REACTIONS];
		this->next_internal = new double [NUM_REACTIONS];
		this->current_internal = new double [NUM_REACTIONS];
		this->cdelay = new complete_delay();
		this->absolute_time = 0;
		this->index = index;
	}
	
	~cell(){
		delete this->cons;
		delete [] this->propen;
		delete [] this->next_internal;
		delete [] this->current_internal;
		delete this->cdelay;
	}
};

struct dependency_graph{
	int ** dependency; // 2D arrays containing the indicies of reactions that if a reactions occur, would need to recalculate the propensity function
	// Notes: there are two reactions: RPSD and RPDD are reactions that we also need to calculate the neighbor cells too
	
	int * dependency_size; // 1D array containing the size of each array of dependent reactions for each reaction
						   // Ex: reaction 1 happens  requires updating 2 other reactions' propensity --> dependency_size[0] = 2
	dependency_graph(){
		this->dependency = new int* [NUM_REACTIONS];
		this->dependency_size = new int [NUM_REACTIONS];
		for (int i = 0; i < NUM_REACTIONS; i++){
			(this->dependency)[i] = new int [8];
			memset((this->dependency)[i], 0, sizeof(int) * 8);
			this->dependency_size[i] = 0;
		}
		this->construct_dependency();
	}
	
	void construct_dependency(){
		(this->dependency)[RPSH1][0] = RPDH1;
		(this->dependency)[RPSH1][1] = RDAH11;
		(this->dependency)[RPSH1][2] = RDAH17;
		(this->dependency_size)[RPSH1] = 3;
		
		(this->dependency)[RPSH7][0] = RPDH7;
		(this->dependency)[RPSH7][1] = RDAH17;
		(this->dependency)[RPSH7][2] = RDAH77;
		(this->dependency_size)[RPSH7] = 3;
		
		(this->dependency)[RPSD][0] = RPDD;
		(this->dependency)[RPSD][1] = RAG1N; // We actually update RAG1N propensities for neighbors, not the cell itself
		(this->dependency)[RPSD][2] = RAG7N; // Same above
		(this->dependency_size)[RPSD] = 3;
		
		(this->dependency)[RPDH1][0] = RPDH1;
		(this->dependency)[RPDH1][1] = RDAH11;
		(this->dependency)[RPDH1][2] = RDAH17;
		(this->dependency_size)[RPDH1] = 3;
		
		(this->dependency)[RPDH7][0] = RPDH7;
		(this->dependency)[RPDH7][1] = RDAH17;
		(this->dependency)[RPDH7][2] = RDAH77;
		(this->dependency_size)[RPDH7] = 3;
		
		(this->dependency)[RPDD][0] = RPDD;
		(this->dependency)[RPDD][1] = RAG1N; // We actually update RAG1N propensities for neighbors, not the cell itself
		(this->dependency)[RPDD][2] = RAG7N; // Same above
		(this->dependency_size)[RPDD] = 3;
		
		(this->dependency)[RPDH11][0] = RPDH11;
		(this->dependency)[RPDH11][1] = RDDH11;
		(this->dependency)[RPDH11][2] = RAG1PH11;
		(this->dependency)[RPDH11][3] = RAG7PH11;
		(this->dependency)[RPDH11][4] = RAGDPH11;
		(this->dependency_size)[RPDH11] = 5;
		
		(this->dependency)[RPDH17][0] = RPDH17;
		(this->dependency)[RPDH17][1] = RDDH17;
		(this->dependency_size)[RPDH17] = 2;
		
		(this->dependency)[RPDH77][0] = RPDH77;
		(this->dependency)[RPDH77][1] = RDDH77;
		(this->dependency_size)[RPDH77] = 2;
		
		(this->dependency)[RDAH11][0] = RPDH1;
		(this->dependency)[RDAH11][1] = RDAH11;
		(this->dependency)[RDAH11][2] = RDAH17;
		(this->dependency)[RDAH11][3] = RPDH11;
		(this->dependency)[RDAH11][4] = RDDH11;
		(this->dependency)[RDAH11][5] = RAG1PH11;
		(this->dependency)[RDAH11][6] = RAG7PH11;
		(this->dependency)[RDAH11][7] = RAGDPH11;
		(this->dependency_size)[RDAH11] = 8;
		
		(this->dependency)[RDAH17][0] = RPDH1;
		(this->dependency)[RDAH17][1] = RDAH11;
		(this->dependency)[RDAH17][2] = RDAH17;
		(this->dependency)[RDAH17][3] = RPDH7;
		(this->dependency)[RDAH17][4] = RDAH77;
		(this->dependency)[RDAH17][5] = RPDH17;
		(this->dependency)[RDAH17][6] = RDDH17;
		(this->dependency_size)[RDAH17] = 7;
		
		(this->dependency)[RDAH77][0] = RPDH7;
		(this->dependency)[RDAH77][1] = RDAH17;
		(this->dependency)[RDAH77][2] = RDAH77;
		(this->dependency)[RDAH77][3] = RPDH77;
		(this->dependency)[RDAH77][4] = RDDH77;
		(this->dependency_size)[RDAH77] = 5;
		
		(this->dependency)[RDDH11][0] = RPDH1;
		(this->dependency)[RDDH11][1] = RDAH11;
		(this->dependency)[RDDH11][2] = RDAH17;
		(this->dependency)[RDDH11][3] = RPDH11;
		(this->dependency)[RDDH11][4] = RDDH11;
		(this->dependency)[RDDH11][5] = RAG1PH11;
		(this->dependency)[RDDH11][6] = RAG7PH11;
		(this->dependency)[RDDH11][7] = RAGDPH11;
		(this->dependency_size)[RDDH11] = 8;
		
		(this->dependency)[RDDH17][0] = RPDH1;
		(this->dependency)[RDDH17][1] = RDAH11;
		(this->dependency)[RDDH17][2] = RDAH17;
		(this->dependency)[RDDH17][3] = RPDH7;
		(this->dependency)[RDDH17][4] = RDAH77;
		(this->dependency)[RDDH17][5] = RPDH17;
		(this->dependency)[RDDH17][6] = RDDH17;
		(this->dependency_size)[RDDH17] = 7;
		
		(this->dependency)[RDDH77][0] = RPDH7;
		(this->dependency)[RDDH77][1] = RDAH17;
		(this->dependency)[RDDH77][2] = RDAH77;
		(this->dependency)[RDDH77][3] = RPDH77;
		(this->dependency)[RDDH77][4] = RDDH77;
		(this->dependency_size)[RDDH77] = 5;
		
		(this->dependency)[RMDH1][0] = RPSH1;
		(this->dependency)[RMDH1][1] = RMDH1;
		(this->dependency_size)[RMDH1] = 2;
		
		(this->dependency)[RMDH7][0] = RPSH7;
		(this->dependency)[RMDH7][1] = RMDH7;
		(this->dependency_size)[RMDH7] = 2;
		
		(this->dependency)[RMDD][0] = RPSD;
		(this->dependency)[RMDD][1] = RMDD;
		(this->dependency_size)[RMDD] = 2;
		
		(this->dependency)[RMSH1][0] = RPSH1;
		(this->dependency)[RMSH1][1] = RMDH1;
		(this->dependency_size)[RMSH1] = 2;
		
		(this->dependency)[RMSH1N][0] = RPSH1;
		(this->dependency)[RMSH1N][1] = RMDH1;
		(this->dependency_size)[RMSH1N] = 2;
		
		(this->dependency)[RAG1PH11][0] = RMSH1;
		(this->dependency)[RAG1PH11][1] = RAG1PH11;
		(this->dependency)[RAG1PH11][2] = RAG1N;
		(this->dependency)[RAG1PH11][3] = RPDH11;
		(this->dependency)[RAG1PH11][4] = RDDH11;
		(this->dependency)[RAG1PH11][5] = RAG7PH11;
		(this->dependency)[RAG1PH11][6] = RAGDPH11;
		(this->dependency)[RAG1PH11][7] = RDG1PH11;
		(this->dependency_size)[RAG1PH11] = 8;
		
		(this->dependency)[RDG1PH11][0] = RMSH1;
		(this->dependency)[RDG1PH11][1] = RAG1PH11;
		(this->dependency)[RDG1PH11][2] = RAG1N;
		(this->dependency)[RDG1PH11][3] = RPDH11;
		(this->dependency)[RDG1PH11][4] = RDDH11;
		(this->dependency)[RDG1PH11][5] = RAG7PH11;
		(this->dependency)[RDG1PH11][6] = RAGDPH11;
		(this->dependency)[RDG1PH11][7] = RDG1PH11;
		(this->dependency_size)[RDG1PH11] = 8;
		
		(this->dependency)[RAG1N][0] = RMSH1;
		(this->dependency)[RAG1N][1] = RAG1PH11;
		(this->dependency)[RAG1N][2] = RAG1N;
		(this->dependency)[RAG1N][3] = RMSH1N;
		(this->dependency)[RAG1N][4] = RDG1N;
		(this->dependency_size)[RAG1N] = 5;
		
		(this->dependency)[RDG1N][0] = RMSH1;
		(this->dependency)[RDG1N][1] = RAG1PH11;
		(this->dependency)[RDG1N][2] = RAG1N;
		(this->dependency)[RDG1N][3] = RMSH1N;
		(this->dependency)[RDG1N][4] = RDG1N;
		(this->dependency_size)[RDG1N] = 5;
		
		(this->dependency)[RMSH7][0] = RPSH7;
		(this->dependency)[RMSH7][1] = RMDH7;
		(this->dependency_size)[RMSH7] = 2;
		
		(this->dependency)[RMSH7N][0] = RPSH7;
		(this->dependency)[RMSH7N][1] = RMDH7;
		(this->dependency_size)[RMSH7N] = 2;
		
		(this->dependency)[RAG7PH11][0] = RMSH7;
		(this->dependency)[RAG7PH11][1] = RAG7PH11;
		(this->dependency)[RAG7PH11][2] = RAG7N;
		(this->dependency)[RAG7PH11][3] = RPDH11;
		(this->dependency)[RAG7PH11][4] = RDDH11;
		(this->dependency)[RAG7PH11][5] = RAG1PH11;
		(this->dependency)[RAG7PH11][6] = RAGDPH11;
		(this->dependency)[RAG7PH11][7] = RDG7PH11;
		(this->dependency_size)[RAG7PH11] = 8;
		
		(this->dependency)[RDG7PH11][0] = RMSH7;
		(this->dependency)[RDG7PH11][1] = RAG7PH11;
		(this->dependency)[RDG7PH11][2] = RAG7N;
		(this->dependency)[RDG7PH11][3] = RPDH11;
		(this->dependency)[RDG7PH11][4] = RDDH11;
		(this->dependency)[RDG7PH11][5] = RAG1PH11;
		(this->dependency)[RDG7PH11][6] = RAGDPH11;
		(this->dependency)[RDG7PH11][7] = RDG7PH11;
		(this->dependency_size)[RDG7PH11] = 8;
		
		(this->dependency)[RAG7N][0] = RMSH7;
		(this->dependency)[RAG7N][1] = RAG7PH11;
		(this->dependency)[RAG7N][2] = RAG7N;
		(this->dependency)[RAG7N][3] = RMSH7N;
		(this->dependency)[RAG7N][4] = RDG7N;
		(this->dependency_size)[RAG7N] = 5;
		
		(this->dependency)[RDG7N][0] = RMSH7;
		(this->dependency)[RDG7N][1] = RAG7PH11;
		(this->dependency)[RDG7N][2] = RAG7N;
		(this->dependency)[RDG7N][3] = RMSH7N;
		(this->dependency)[RDG7N][4] = RDG7N;
		(this->dependency_size)[RDG7N] = 5;
		
		(this->dependency)[RMSD][0] = RPSD;
		(this->dependency)[RMSD][1] = RMDD;
		(this->dependency_size)[RMSD] = 2;
		
		(this->dependency)[RAGDPH11][0] = RMSD;
		(this->dependency)[RAGDPH11][1] = RAGDPH11;
		(this->dependency)[RAGDPH11][2] = RPDH11;
		(this->dependency)[RAGDPH11][3] = RDDH11;
		(this->dependency)[RAGDPH11][4] = RAG1PH11;
		(this->dependency)[RAGDPH11][5] = RAG7PH11;
		(this->dependency)[RAGDPH11][6] = RDGDPH11;
		(this->dependency_size)[RAGDPH11] = 7;
		
		(this->dependency)[RDGDPH11][0] = RMSD;
		(this->dependency)[RDGDPH11][1] = RAGDPH11;
		(this->dependency)[RDGDPH11][2] = RPDH11;
		(this->dependency)[RDGDPH11][3] = RDDH11;
		(this->dependency)[RDGDPH11][4] = RAG1PH11;
		(this->dependency)[RDGDPH11][5] = RAG7PH11;
		(this->dependency)[RDGDPH11][6] = RDGDPH11;
		(this->dependency_size)[RDGDPH11] = 7;
	}
	
	~dependency_graph(){
		delete[] this->dependency_size;
		for (int i = 0; i < NUM_REACTIONS; i++){
			delete[] (this->dependency)[i];
		}
		delete[] this->dependency;
	}
};

struct propensities{
	void (*prop_funs[NUM_REACTIONS]) (cell&, rates&);
	dependency_graph * dgraph;
	
	propensities(dependency_graph* dg){
		memset(this->prop_funs, 0, sizeof(this->prop_funs));
		this->dgraph = dg;
	}
	
	~propensities(){
		
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
#endif
