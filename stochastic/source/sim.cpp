#include <random>
#include <math.h>
#include "sim.hpp"
#include "debug.hpp"
#include "init.hpp"
#include "io.hpp"
#include "randDist.hpp"

extern terminal* term; // Declared in init.cpp
default_random_engine generator;
uniform_real_distribution <double> distribution (0.0, 1.0);

void simulate_all_params_sets(input_params& ip, parameters& pr){
	embryo em(ip);
	rates rs;
	// declare reactions and propensities: lists of function calls to use throughout the simulation
	sim_data sd;
	init_propensities(sd);
	init_reactions(sd);
	
	for (int i = 0; i < ip.num_sets; i++){
		//Create set directories
		if (ip.print_cons){
			create_set_directory(i + 1, ip);
		}
		// copy params from parameters to rates
		memcpy(rs.data, pr.data[i], sizeof(double) * NUM_RATES);
		// simulate each set
		int set_score = simulate_one_param_set(ip, sd, rs, i, em);
		// reset rates after each set
		rs.reset();
	}
}

int simulate_one_param_set(input_params& ip, sim_data& sd, rates& rs, int set_index, embryo& em){
	int score = 0; 
	cout << term->blue << "Simulating set " << term->reset << set_index << " . . ." << endl;
	
	for (int i = 0; i < ip.num_mutants; i++){
		// create mutant directory
		if (ip.print_cons){
			create_mutant_directory(set_index + 1, ip.mutants[i], ip);
		}
		// change rates based on mutant
		// simulate_mutant
		score += simulate_mutant(ip, sd, rs, em, set_index, ip.mutants[i]);
		// revert_rates
	}
	// print score to pipe if necessary
	return score;
}

int simulate_mutant(input_params& ip, sim_data& sd, rates& rs, embryo& em, int set_index, int mutant_index){
	// reset embryo : for each cell:
	// cons: current_cons to all 0; last_change_time to all 0; time_record to all 0; cons_record to all 0
	// propen[NUM_REACTIONS] all to 0
	// next_internal [NUM_REACTIONS]  all to 0
	// current_internal [NUM_REACTIONS] all to 0
	// cdelay empty
	// internal_time = 0
	em.reset();
	// core simulation
	core_simulation(ip, sd, rs, em);
	// test conditions
	// print data if necessary
	return 0;
}

void core_simulation(input_params& ip, sim_data& sd, rates& rs, embryo& em){
	// update initial conditions : all gene has current concentrations to 2
	update_initial_concentrations(em);
	// calculate propensity for each reaction
	calculate_initial_propensity(sd, rs, em);
	calculate_initial_next_internal(em);
	test_next_firing(em);
}


void update_initial_concentrations (embryo& em){
	for (int i = 0; i < em.num_cells; i++){
		(((em.cell_list[i])->cons)->current_cons)[G1] = 2;
		(((em.cell_list[i])->cons)->current_cons)[G7] = 2;
		(((em.cell_list[i])->cons)->current_cons)[GD] = 2;
	}
}

void calculate_initial_propensity(sim_data& sd, rates& rs, embryo& em){
	for (int i = 0; i < em.num_cells; i++){
		for (int j = 0; j < NUM_REACTIONS; j++){
			// propensity functions will update the propensity into cell.propen[NUM]
			sd.prop_funs[j](em, i, rs);
		}
	}
}

void calculate_initial_next_internal(embryo& em){
	double rand_num;
	for (int i = 0; i < em.num_cells; i++){
		for (int j = 0; j < NUM_REACTIONS; j++){
			rand_num = distribution(generator);
			((em.cell_list[i])->propen)[j] = log(1.0/rand_num);
		}
	}
}

double generate_unif_rand_firing_time(){
	double rand_num = distribution (generator);
	return log(1.0/rand_num);
}
