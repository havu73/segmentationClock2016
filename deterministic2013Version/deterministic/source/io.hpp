#ifndef IO_HPP
#define IO_HPP
#include "structs.hpp"
using namespace std;
void store_filename (char**,const char*);
void open_file (ofstream*, char*, bool);
void read_pipe (parameters&, input_params&);
void read_pipe_int (int, int*);
void read_file (input_data*);
bool parse_param_line (parameters&, int, char*, int&);
void parse_ranges_file (pair <double, double>[], char*);
bool not_EOL (char);
void create_set_directory (int, input_params&);
void print_concentrations(input_params&, int, int, con_levels&);
void write_pipe(double*, input_params&);
#endif
