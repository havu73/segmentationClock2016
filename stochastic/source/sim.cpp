#include <cmath>
#include "sim.hpp"
#include "debug.hpp"
#include "init.hpp"
#include "io.hpp"
#include "randDist.hpp"
#include "tests.hpp"

extern terminal* term; // Declared in init.cpp


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
	double sresScore;
	//cout << term->blue << "Simulating set " << term->reset << set_index << " . . ." << endl;
	
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
	// real score to pass into SRES
	sresScore = (double) MAX_SCORE - score;
	// print score to pipe if necessary
	if (ip.piping){
		write_pipe(&sresScore, ip);
	}
	return score;
}

int simulate_mutant(input_params& ip, sim_data& sd, rates& rs, embryo& em, int set_index, int mutant_index){
	int mutant_score = 0;
	// reset embryo : for each cell:
	// cons: current_cons to all 0; cons_record to all 0
	// propen[NUM_REACTIONS] all to 0
	// next_internal [NUM_REACTIONS]  all to 0
	// current_internal [NUM_REACTIONS] all to 0
	// cdelay empty
	// internal_time = 0
	em.reset();
	// core simulation
	core_simulation(ip, sd, rs, em);
	// test conditions
	features wtf(ip.num_cells, ip.num_bin);
	
	if (mutant_index == WT){
		mutant_score = test_wildtype(ip, em, wtf);
	}	
	// print data if necessary
	if (ip.print_cons){
		print_concentrations(ip, set_index, mutant_index, em);
	}
	return mutant_score;
}

void core_simulation(input_params& ip, sim_data& sd, rates& rs, embryo& em){
	// update initial conditions : all gene has current concentrations to 2
	update_initial_concentrations(em);
	// calculate propensity for each reaction
	calculate_initial_propensity(sd, rs, em);
	calculate_initial_next_internal(em);
	
	// time keeping variables
	bool done = false;
	int record_per_done = ip.check_done_granularity / ip.record_granularity; // How many time to checks record per time checking done
	
	// next reaction to denote parameters about what reaction in what cell and whether it is a complete delay reaction
	// that is firing next
	next_reaction nr;

	// enter simulation loop
	while (!done){
		for (int i = 0; i < record_per_done; i++){
			for (int j = 0; j < ip.record_granularity; j++){
				// 1.find the next reaction
				find_next_reaction(nr, em);
			
				// 2. update internal time for each reaction : Tk = Tk + Ak * DELTA
				update_current_internal(em, nr.delta);
					
				// 3. update absolute time of embryo
				em.absolute_time = em.absolute_time + nr.delta;
				
				// 4. reactions happen
				sd.reac_funs[nr.reaction_index](em, nr.cell_index, sd, rs, nr.delay_complete);
				
				// 5. find the next firing time of the reaction that just initiated
				if (! (nr.delay_complete)){
					find_next_firing(em.cell_list[nr.cell_index], nr.reaction_index);
				}
			}
			// transfer data from current time and concentration to record time and concentration
			transfer_to_record(em, ip);
		}
		// check whether done or not
		done = check_done(em, ip);
		
	}
}


void update_initial_concentrations (embryo& em){
	for (int i = 0; i < em.num_cells; i++){
		((em.cell_list[i])->current_cons)[G1] = 2;
		((em.cell_list[i])->current_cons)[G7] = 2;
		((em.cell_list[i])->current_cons)[GD] = 2;
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
	for (int i = 0; i < em.num_cells; i++){
		for (int j = 0; j < NUM_REACTIONS; j++){
			((em.cell_list[i])->next_internal)[j] = pk_dist();
		}
	}
}

void find_next_reaction(next_reaction& nr, embryo& em){
	double min_delta = INFINITY;
	double current_delta;
	int c_index;
	int r_index;
	bool d_complete = false;
	cell* cc;
	for (int i = 0; i < em.num_cells; i++){
		cc = em.cell_list[i];
		for (int j = 0; j < NUM_REACTIONS; j++){
			if ((cc->propen)[j] > 0){
				// deltak = (Pk - Tk) / ak (next internal - current_internal) / propensity
				current_delta = ((cc->next_internal)[j] - (cc->current_internal)[j]) / ((cc->propen)[j]);
				if (min_delta > current_delta){
					min_delta = current_delta;
					c_index = i;
					r_index = j;
				}
			}
		}
	}
	
	for (int i = 0; i < em.num_cells; i++){
		cc = em.cell_list[i];
		if (!((cc->cdelay)->is_empty())){
			current_delta = (cc->cdelay)->see_soonest() - em.absolute_time;
			if (min_delta > current_delta){
				min_delta = current_delta;
				c_index = i;
				r_index = (cc->cdelay)->soonest_reaction();
				d_complete = true;
			}
		}
	}
	
	nr.cell_index = c_index;
	nr.reaction_index = r_index;
	nr.delay_complete = d_complete;
	nr.delta = min_delta;
}



void update_current_internal(embryo& em, double delta){
	cell* cc;
	// For each reaction in each cell: Tk = Tk + Ak * Delta
	for (int i = 0; i < em.num_cells; i++){
		cc = em.cell_list[i];
		for (int j = 0; j < NUM_REACTIONS; j++){
			(cc->current_internal)[j] = (cc->current_internal)[j] + (cc->propen)[j] * delta;
		}
	}
}

void find_next_firing(cell* current_cell, int reaction_index){
	(current_cell->next_internal)[reaction_index] = (current_cell->next_internal)[reaction_index] + pk_dist(); 
}

void transfer_to_record(embryo& em, input_params& ip){
	for (int i = 0; i < ip.num_cells; i++){
		(em.cell_list[i])->transfer_record(MH1, KEEPMH1);
		(em.cell_list[i])->transfer_record(MH7, KEEPMH7);
		for (int j = 0; j < ip.num_print_states; j++){
			(em.cell_list[i])->transfer_record(ip.print_states[j], NUM_KEEP_STATES + j);
		}
	}
	em.transfer_time_record();
}

bool check_done(embryo& em, input_params& ip){
	if (em.absolute_time >= ip.time_total){
		return true;
	}
	return false;
}
