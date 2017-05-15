#ifndef INIT_HPP
#define INIT_HPP

#include "structs.hpp"
using namespace std;

void init_terminal();
void free_terminal();
void init_verbosity(input_params&);
void accept_input_params (int, char**, input_params&);
void ensure_nonempty (const char*, const char*);
inline bool option_set (const char* , const char* , const char* );
void check_input_params (input_params&);
void init_seeds (input_params&, int, bool, bool);
int generate_seed ();
void read_sim_params(input_params&, parameters&, input_data&, input_data&);
double random_double (pair<double, double>);
#endif
