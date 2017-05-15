#ifndef SIM_HPP
#define SIM_HPP
#include "structs.hpp"
void simulate_all_params_sets(input_params&, parameters&);
void process_rates(rates&, parameters&, input_params&, int);
int simulate_one_param_set (input_params&, rates&, int);
bool core_simulation (input_params&, rates&, con_levels&);
bool check_for_nan (con_levels&, int);
void find_sim_steps (rates&, input_params&);
int simulate_mutant(input_params&, rates&, con_levels&, int, int);
#endif
