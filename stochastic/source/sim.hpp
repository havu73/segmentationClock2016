#ifndef SIM_HPP
#define SIM_HPP
#include "structs.hpp"
void simulate_all_params_sets(input_params&, parameters&);
int simulate_one_param_set(input_params&, sim_data&, rates&, int, embryo&);
int simulate_mutant(input_params&, sim_data&, rates&, embryo&, int, int);
void core_simulation(input_params&, sim_data&, rates&, embryo&);
void update_initial_concentrations (embryo&);
void calculate_initial_propensity(sim_data&, rates&, embryo&);
void calculate_initial_next_internal(embryo&);
void find_next_reaction(next_reaction&, embryo&);
void update_current_internal(embryo&, double);
void find_next_firing(cell*, int);
void transfer_to_record(embryo&, input_params&);
bool check_done(embryo&, input_params&);
#endif

