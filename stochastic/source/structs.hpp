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
#include <cmath> // Need for INFINITY in priority queue
#include "macros.hpp"
#include "memory.hpp"
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
	char* ranges_file; // The path and name of the parameter ranges file, default=none
	bool read_ranges; // Whether or not to read the ranges file, default=false
	
	// What states whose concentrations to store
	int* print_states; // list of indices of states that we want to store concentrations over time for plotting and testing
	int num_print_states; // number of states to print the concentrations over time
	bool print_cons; // whther or not program is requested to print concentrations over time. False when ip->num_print_states = 0 and users do not specify to print the default MH1, MH7, MD; 
	char* out_dir; // The path of the output directory for concentrations or oscillation features, default=none	
	bool has_out_dir; // whether user specify an output directory or not
	
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
	int * mutants; // List of index of mutants that users want to simulate and test. Default : all possible mutants
	int num_mutants; // number of mutants to be simulated. Default = NUM_MUTANTS
	int max_cond_score; // maximum scores when all conditions specified by users are passed.
	
	// Piping data
	bool piping; // Whether or not input and output should be piped (as opposed to written to disk), default=false
	int pipe_in; // The file descriptor to pipe data from, default=none (0)
	int pipe_out; // The file descriptor to pipe data into, default=none (0)
	
	// Output stream data
	bool verbose; // Whether or not the program is verbose, i.e. prints many messages about program and simulation state, default=false
	bool quiet; // Whether or not the program is quiet, i.e. redirects cout to /dev/null, default=false
	streambuf* cout_orig; // cout's original buffer to be restored at program completion
	ofstream* null_stream; // A stream to /dev/null that cout is redirected to if quiet mode is set
	
	// number of cells in the embryo
	int num_cells;
	
	// number of bin to do noise analysis
	int num_bin;
	
	// Parameters related to the speed and data collection of the simulation
	int check_done_granularity; // every 120 reactions fired, we will check each cell and see whether or not they have reached a time point
									  // above the time we need to simulate
	int record_granularity;	// Every 30 reactions, we will record the concentrations and time points of states that users want to keep track of 
	// For efficiency reasons, we require that check_done_granularity is divisible by record_granularity
	input_params(){
		// IO files
		this->params_file = new char[30];
		this->read_params = false;
		this->ranges_file = new char [30];
		this->read_ranges = false;		
		
		// printing concentrations over time
		this->print_states = 0; // NULL, need to be initialzied using new int[size] later if users want to print
		this->num_print_states = 0;
		this->print_cons = false;
		this->out_dir = new char [30];
		this->has_out_dir = false;
		
		//sets
		this->num_sets = 1;
		
		// timing
		this->time_total = 610;
		
		//seeds
		this->seed = 0;
		this->reset_seed = false;
		this->pseed = 0;
		this->store_pseed = false;
		this->seed_file = new char [30];
		this->print_seeds = false;		
		
		//mutants
		this->mutants = new int [NUM_MUTANTS]; // NULL
		for (int i = 0; i< NUM_MUTANTS; i++){
			(this->mutants)[i] = i;
		}
		this->num_mutants = NUM_MUTANTS;
		this->max_cond_score = TOTAL_SC;
		
		//piping data
		this->piping = false;
		this->pipe_in = 0;
		this->pipe_out = 0;
		
		//output stream data
		this->verbose = false;
		this->quiet = false;
		this->cout_orig = NULL;
		this->null_stream = new ofstream("/dev/null");
		
		// number of cells in the embryo
		this->num_cells = 4;
		
		// number of bins to do data analysis
		this->num_bin = DEFAULT_NUM_BIN;
		
		// parameters about data transfer and time control of the simulation
		this->check_done_granularity = 2000;
		this->record_granularity = 1000;
	}
	
	~input_params(){
		delete [] this->params_file;
		delete [] this->ranges_file;
		delete [] this->mutants;
		delete [] this->seed_file;
		//delete [] this->print_states;
		delete this->null_stream;
		delete [] this->out_dir;
		if (this->print_states != 0){
			delete [] this->print_states;
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
	
	int soonest_reaction(){
		return ((this->pq)->top()).reaction_index;
	}
	
	void complete_soonest(){
		(this->pq)->pop();
	}
	
	void initiate_delay(int function_index, double next_complete){
		delay_reaction dr (function_index, next_complete);
		(this->pq)->push(dr);
	}
	
	int size(){
		return (this->pq)->size();
	}
	
	void reset(){
		while(!(this->pq)->empty()){
			(this->pq)->pop();
		}
		// Need to always keep an infinity fake reaction because I have 
		// to find soonest reaction in core_simulation(sim.cpp)  too often, so I dont want to check
		// pq empty more often than necessary
		delay_reaction dr (NUM_REACTIONS, INFINITY);
		(this->pq)->push(dr);
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
		memset(this->data, 0, sizeof(double) * NUM_RATES);
	}
	
	~rates () {
		delete [](this->data);
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
	bool initialized; // Whether or not this struct's data have been initialized
	int * current_cons;
	vector<int>** cons_record; // a list of vectors (dynamically growing lists) containing the concentrations of different states over time. 
							   // Needed to keep records of concentrations of some states over time
									// that the concentration of each states that we keep record of is recorded
	int num_records;
	
	double * propen;		// propensities of reactions inside the cell
	double * next_internal;		// next internal time of reactions in the system
	double * current_internal; // current internal time of reactions in the system
	complete_delay * cdelay; // the priority queue of the complettion time of delay reactions going on at a specific moment in the cell
	int index; 	//the index of the cell inside the embryo
	
	cell(){
		this->initialized = false;
	}
	
	cell(int index, input_params& ip){
		if (this->initialized){
			this->reset();
		}
		else{
			this->num_records = NUM_KEEP_STATES + ip.num_print_states;
			this->current_cons = new int [NUM_STATES];
			this->cons_record = new vector<int> * [NUM_KEEP_STATES + ip.num_print_states];

			for (int i = 0; i < NUM_KEEP_STATES; i++){
				(this->cons_record)[i] = new vector<int>;
				(this->cons_record)[i]->reserve(50000);
			}
			
			for (int i = 0; i < ip.num_print_states; i++){
				(this->cons_record)[i + NUM_KEEP_STATES] = new vector<int>;
				(this->cons_record)[i + NUM_KEEP_STATES]->reserve(50000);
			}
			
			this->propen = new double[NUM_REACTIONS];
			this->next_internal = new double [NUM_REACTIONS];
			this->current_internal = new double [NUM_REACTIONS];
			this->cdelay = new complete_delay();
			this->index = index;
			this->initialized = true;
		}
	}
	
	void reset(){
		memset(this->current_cons, 0, sizeof(int) * NUM_STATES);
		
		for (int i = 0; i < (this->num_records); i++){
			(this->cons_record)[i]->clear();
		}
		
		memset(this->propen, 0, sizeof(double) * NUM_REACTIONS);
		memset(this->next_internal, 0, sizeof(double) * NUM_REACTIONS);
		memset(this->current_internal, 0, sizeof(double) * NUM_REACTIONS);
		(this->cdelay)->reset();
	}
	
	void transfer_record(int real_index, int record_index){
		(this->cons_record[record_index])->push_back(this->current_cons[real_index]);
	}
	
	~cell(){
		delete [] this->current_cons;
		for (int i = 0; i < (this->num_records); i++){
			delete (this->cons_record)[i];
		}
		delete [] this->cons_record;
		delete [] this->propen;
		delete [] this->next_internal;
		delete [] this->current_internal;
		delete this->cdelay;
	}
};

struct embryo{
	cell ** cell_list;
	int ** neighbors;
	vector<double>* time_record;
	double absolute_time;
	int neighbor_per_cell; // number of neighbors that each cell has
	int num_cells;
	
	embryo(input_params& ip){
		this->time_record = new vector<double>;
		(this->time_record)->reserve(50000);
		this->num_cells = ip.num_cells;
		this->cell_list = new cell * [this->num_cells];
		for (int i = 0; i < num_cells; i++){
			this->cell_list [i] = new cell(i, ip);
		}
		this->neighbors = new int* [this->num_cells];
		
		if (this->num_cells == 2){
			for (int i = 0; i < num_cells; i++){
				(this->neighbors)[i] = new int [1];
			}
		}
		else if (this->num_cells == 4){
			for (int i = 0; i < num_cells; i++){
				(this->neighbors)[i] = new int [FOUR_NEIGHBORS];
			}
		}
		else if (this->num_cells == 16){
			for (int i = 0; i < this->num_cells; i++){
				(this->neighbors)[i] = new int [SIXTEEN_NEIGHBORS];
				memset((this->neighbors)[i], 0, sizeof(int) * SIXTEEN_NEIGHBORS);
			}
		}
		this->construct_neighbors();
	}
	
	void construct_neighbors(){
		if (this->num_cells == 2){
			this->neighbor_per_cell = TWO_NEIGHBORS;
			(this->neighbors)[0][0] = 1;
			(this->neighbors)[1][0] = 0;
		}
		else if (this->num_cells == 4){
			this->neighbor_per_cell = FOUR_NEIGHBORS;
			int zero [] = {1,2,3};
			int one [] = {0,2,3};
			int two [] = {0,1,3};
			int three [] = {0,1,2};
			memcpy((this->neighbors)[0], zero, sizeof(int) * FOUR_NEIGHBORS);
			memcpy((this->neighbors)[1], one, sizeof(int) *	FOUR_NEIGHBORS);
			memcpy((this->neighbors)[2], two, sizeof(int) * FOUR_NEIGHBORS);
			memcpy((this->neighbors)[3], three, sizeof(int) * FOUR_NEIGHBORS);
		}
		else if(this->num_cells == 16){
			this->neighbor_per_cell = SIXTEEN_NEIGHBORS;
			int zero [] = {1,4,3,7,12,13};
			int one [] = {0,4,5,2,12,13};
			int two [] = {1,3,5,6,13,14};
			int three [] = {2,6,7,14,15,0};
			int four [] = {0,1,5,9,8,7};
			int five [] = {1,2,6,4,9,10};
			int six [] = {2,3,7,11,10,5};
			int seven [] = {3,6,11,0,4,8};
			int eight [] = {4,9,12,7,11,15};
			int nine [] = {4,5,10,13,12,8};
			int ten [] = {5,6,11,14,13,9};
			int eleven [] = {6,7,8,15,14,10};
			int twelve [] = {8,9,13,1,0,15};
			int thirteen [] = {9,10,14,2,1,12};
			int fourteen [] = {10,11,15,3,2,13};
			int fifteen [] = {11,8,12,3,2,14};
			memcpy((this->neighbors)[0], zero, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[1], one, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[2], two, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[3], three, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[4], four, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[5], five, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[6], six, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[7], seven, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[8], eight, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[9], nine, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[10], ten, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[11], eleven, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[12], twelve, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[13], thirteen, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[14], fourteen, sizeof(int) * SIXTEEN_NEIGHBORS);
			memcpy((this->neighbors)[15], fifteen, sizeof(int) * SIXTEEN_NEIGHBORS);
		}
	}
	
	void transfer_time_record(){
		(this->time_record)->push_back(this->absolute_time);	
	}
	
	void reset(){
		for (int i = 0; i < this->num_cells; i ++){
			((this->cell_list)[i])->reset();
		}
		(this->time_record)->clear();
	}
	
	~embryo(){
		delete [] this->cell_list;
		/*
		for (int i = 0; i < this->num_cells; i++){
			delete [] this->neighbors[i];
		}
		* */
		delete [] this->neighbors;
	}
};




struct sim_data{
	void (*reac_funs[NUM_REACTIONS]) (embryo&, int, sim_data&, rates&, bool);
	void (*prop_funs[NUM_REACTIONS]) (embryo&, int, rates&);
	int ** dependency; // 2D arrays containing the indicies of reactions that if a reactions occur, would need to recalculate the propensity function
	// Notes: there are two reactions: RPSD and RPDD are reactions that we also need to calculate the neighbor cells too
	
	int * dependency_size; // 1D array containing the size of each array of dependent reactions for each reaction
						   // Ex: reaction 1 happens  requires updating 2 other reactions' propensity --> dependency_size[0] = 2
	
	sim_data(){
		// 1. reac_funs
		memset(this->reac_funs, 0, sizeof(this->reac_funs));
		
		// 2. propensity functions
		memset(this->prop_funs, 0, sizeof(this->reac_funs));
		
		// 3. dependency graph
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
		(this->dependency_size)[RPSD] = 1;
		
		(this->dependency)[RPDH1][0] = RPDH1;
		(this->dependency)[RPDH1][1] = RDAH11;
		(this->dependency)[RPDH1][2] = RDAH17;
		(this->dependency_size)[RPDH1] = 3;
		
		(this->dependency)[RPDH7][0] = RPDH7;
		(this->dependency)[RPDH7][1] = RDAH17;
		(this->dependency)[RPDH7][2] = RDAH77;
		(this->dependency_size)[RPDH7] = 3;
		
		(this->dependency)[RPDD][0] = RPDD;
		(this->dependency_size)[RPDD] = 1;
		
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

	~sim_data(){
		delete[] this->dependency_size;
		for (int i = 0; i < NUM_REACTIONS; i++){
			delete[] (this->dependency)[i];
		}
		delete[] this->dependency;
	}
};

struct peak_trough{
	double* smooth_cons;
	vector<int>* peaks;
	vector<int>* troughs;
	int capacity;
	int num_cons;
	
	peak_trough(int smooth_size){
		this->capacity = smooth_size;
		this->num_cons = 0;
		this->smooth_cons = new double [smooth_size];
		this->peaks = new vector<int>;
		this->troughs = new vector<int>;
		(this->peaks)->reserve(200);
		(this->troughs)->reserve(200);
	}
	
	void reset(){
		this->clear(); 
	}
	
	void clear(){
		this->num_cons = 0;
		memset(this->smooth_cons, 0, sizeof(double) * this->capacity);
		(this->peaks)->clear();
		(this->troughs)->clear();
	}
	
	~peak_trough(){
		delete [] this->smooth_cons;
		delete this->peaks;
		delete this->troughs;
	}
};

struct features{
	double** avg_amplitude;
	double** mid_ptt;
	double** last_ptt;
	double* intrinsic;
	double* extrinsic;
	double* avg_cons;
	int num_cells;
	int num_bin;
	
	features(int num_cells, int num_bin){
		this->num_cells = num_cells; 
		this->avg_amplitude = new double* [NUM_KEEP_STATES];
		this->mid_ptt = new double* [NUM_KEEP_STATES];
		this->last_ptt = new double* [NUM_KEEP_STATES];
		for (int i = 0; i < NUM_KEEP_STATES; i ++){
			this->avg_amplitude[i] = new double [num_cells];
			this->mid_ptt[i] = new double [num_cells];
			this->last_ptt[i] = new double [num_cells];
		}
		
		this->intrinsic = new double[num_bin];
		this->extrinsic = new double[num_bin];
		this->avg_cons = new double [num_bin];
		this->num_bin = num_bin;
	}
	
	void reset(){
		for (int i = 0; i < NUM_KEEP_STATES; i ++){
			memset(this->avg_amplitude[i], 0, sizeof(double) * this->num_cells);
			memset(this->mid_ptt[i], 0, sizeof(double) * this->num_cells);
			memset(this->last_ptt[i], 0, sizeof(double) * this->num_cells);
		}
		memset(this->intrinsic, 0, sizeof(double) * this->num_bin);
		memset(this->extrinsic, 0, sizeof(double) * this->num_bin);
		memset(this->avg_cons, 0, sizeof(double) * this->num_bin);
	}
	
	~features(){
		for (int i = 0; i < NUM_KEEP_STATES; i ++){
			delete [] (this->avg_amplitude)[i];
			delete [] (this->mid_ptt)[i];
			delete [] (this->last_ptt)[i];
		}
		delete [] this->avg_amplitude;
		delete [] this->mid_ptt;
		delete [] this->last_ptt;
		
		delete [] this->intrinsic;
		delete [] this->extrinsic;
		delete [] this->avg_cons;
	}
};

struct binned_data{
	vector<double> ** her;
	vector<double> ** in_noise;
	vector<double> ** ex_noise;
	int num_bin;
	
	binned_data(){
		this->initialize(DEFAULT_NUM_BIN);
	}
	
	binned_data(int num_bin){
		this->initialize(num_bin);
	}
	
	void initialize (int num_bin){
		this->her = new vector<double> * [num_bin];
		this->in_noise = new vector<double> * [num_bin];
		this->ex_noise = new vector<double> * [num_bin];
		for (int i = 0; i < num_bin; i ++){
			this->her[i] = new vector<double>;
			(this->her[i])->reserve(40);
			this->in_noise[i] = new vector<double>;
			(this->in_noise[i])->reserve(40);
			this->ex_noise[i] = new vector<double>;
			(this->ex_noise[i])->reserve(40);
		}
		this->num_bin = num_bin;
	}
	~binned_data(){
		for (int i = 0; i < this->num_bin; i ++){
			delete (this->her)[i];
			delete (this->in_noise)[i];
			delete (this->ex_noise)[i];
		}
		delete [] this->her;
		delete [] this->in_noise;
		delete [] this->ex_noise;
	}
};

struct slices{
	double ** her1; // concentrations of her1 mRNA of each cell in the slice
	double ** her7; // concentrations of her7 mRNA of each cell in the slice
	double * avg_h1;
	double * avg_h7;
	double * avg_her;
	double * in_noise; // intrinsic noise in this slice
	double * ex_noise; // extrinsic noise in this slice
	double min_her;
	double max_her;
	int total_cells;
	slices(int num_cells){
		this->her1 = new double* [NUM_SLICES];
		this->her7 = new double* [NUM_SLICES];
		this->avg_h1 = new double [NUM_SLICES];
		this->avg_h7 = new double [NUM_SLICES];
		this->avg_her = new double [NUM_SLICES];
		this->in_noise = new double [NUM_SLICES];
		this->ex_noise = new double [NUM_SLICES];
		this->total_cells = num_cells;
		for (int i = 0; i < NUM_SLICES; i++){
			this->her1[i] = new double[num_cells];
			this->her7[i] = new double[num_cells];
		}
		this->initiate_values();
	}
	
	void initiate_values (){
		for (int i = 0; i < NUM_SLICES; i++){
			for (int j = 0;j < this->total_cells; j ++){
				this->her1[i][j] = -1;
				this->her7[i][j] = -1;
			}
			this->avg_h1[i] = -1;
			this->avg_h7[i] = -1;
			this->avg_her[i] = -1;
			this->in_noise[i] = -1;
			this->ex_noise[i] = -1;
		}
	}
	~slices(){
		for (int i = 0; i < NUM_SLICES; i++){
			delete [] this->her1[i];
			delete [] this->her7[i];
		}
		delete [] this->her1;
		delete [] this->her7;
		delete [] this->avg_h1;
		delete [] this->avg_h7;
		delete [] this->avg_her;
		delete [] this->in_noise;
		delete [] this->ex_noise;
	}
};

struct next_reaction{
	int reaction_index;
	int cell_index;
	double delta;
	bool delay_complete;
	
	next_reaction(){
		this->delta= INFINITY;
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
