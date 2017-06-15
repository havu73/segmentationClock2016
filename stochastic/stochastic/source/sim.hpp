#ifndef SIM_HPP
#define SIM_HPP
#include "structs.hpp"
void simulate_all_params_sets(input_params&, parameters&);
double simulate_one_param_set(input_params&, sim_data&, rates&, int, embryo&);
void process_rates(rates&, parameters&, int);
double simulate_mutant(input_params&, sim_data&, rates&, embryo&, features&, int, int);
double simulate_dapt_mutant(input_params&, sim_data&, rates&, embryo&, features&, int, int);
void change_rates_based_on_mutants(rates&, int);
bool core_simulation_dapt(input_params&, sim_data&, rates&, embryo&);
bool core_simulation(input_params&, sim_data&, rates&, embryo&);
void calculate_gradient_delay_rates(int, rates&, double*);
void calculate_gradient_per_step (input_params&, rates&, double*);
bool states_out_of_bound(embryo&);
void update_initial_concentrations (embryo&);
void calculate_initial_propensity(sim_data&, rates&, embryo&);
void calculate_initial_next_internal(input_params&, embryo&);
void find_next_reaction(next_reaction&, embryo&);
void update_current_internal(embryo&, double);
void find_next_firing(input_params&, cell*, int);
void transfer_to_record(embryo&, input_params&);
bool check_done(embryo&, input_params&);
#endif

