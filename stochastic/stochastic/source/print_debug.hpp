#include "structs.hpp"
#include "macros.hpp"

void open_debug_steam(input_params&);
void print_rates(input_params&, rates&);
void print_propensities(input_params&, embryo&);
void print_current_internal(input_params&, embryo&);
void print_next_internal(input_params&, embryo&);
void print_states(input_params&, embryo&);
void print_initial_simulation(input_params&, embryo&);
void print_one_round_simulation(embryo&, input_params&, next_reaction&);
void open_random_steam(input_params&);
void print_random_number(input_params&, double*);
