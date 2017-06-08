#ifndef IO_HPP
#define IO_HPP
#include "structs.hpp"
using namespace std;
void store_filename (char**,const char*);
void open_file (ofstream*, char*, bool);
void open_file_quiet (ofstream* , char*);
void read_pipe (parameters&, input_params&);
void read_pipe_int (int, int*);
void read_file (input_data*);
bool parse_param_line (parameters&, int, char*, int&);
void parse_ranges_file (pair <double, double>[], char*);
bool not_EOL (char);
void create_set_directory (int, input_params&);
void create_mutant_directory(int, int, input_params&);
void print_concentrations(input_params&, int, int, embryo&);
void create_concentrations_file_name(int, input_params&, char*, char**);
void print_one_state_concentrations(embryo&, int, char**);
void write_pipe(double*, input_params&);
void print_smooth_data (input_params&, embryo&, peak_trough&, int, int);
void print_features_data (input_params&, embryo&, features&, int, int);
void print_perturb_rates(input_params&, rates&, parameters&, int);
#endif
