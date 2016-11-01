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
	if (ip.piping){
		write_pipe((double *) &score, ip);
	}
	return score;
}

int simulate_mutant(input_params& ip, sim_data& sd, rates& rs, embryo& em, int set_index, int mutant_index){
	int mutant_score = 0;
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
	cout << em.cell_list[1]->absolute_time << endl;
	cout << ((em.cell_list[1]->cons_record)[KEEPMH1])->size() << endl;
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
	
	// Delta K array: how much time will pass from the current time for each reaction to occur, plus the firing time for the soonest-complete delayed reaction 
	double* deltaKArray = new double[NUM_REACTIONS + 1];
	
	// some parameters used repeatedly during simulation
	int next_fire; // index of reaction that will fire next
	double delta; // time from current time a reaction will occur
	int complete_delay_index;
	
	// enter simulation loop
	while (!done){
		for (int i = 0; i < record_per_done; i++){
			for (int j = 0; j < ip.record_granularity; j++){
				for (int k = 0; k < ip.num_cells; k++){
					// 1.find Delta
					next_fire = update_deltaK_array(&deltaKArray, em.cell_list[k]);
					delta = deltaKArray[next_fire];
					
					// 2. update internal time for each reaction : Tk = Tk + Ak * DELTA
					update_current_internal(em.cell_list[k], delta);
					
					// 3. update absolute time for cell
					(em.cell_list[k])->absolute_time = (em.cell_list[k])->absolute_time + delta;
					
					// 4. take actions based on what reaction is chosen to fire: update concentrations and also update propensity for involved reactions
					if (next_fire != NUM_REACTIONS){ // if either non-delay reaction or delay reaction initiate 
						// Call reaction: update concentrations and also propensity
						sd.reac_funs[next_fire](em, k, sd, rs, false); // complete = false because for delay reactions, this is just initiation
						// find the next firing time of the reaction that just fired
						find_next_firing(em.cell_list[k], next_fire);
					}
					else{ // if a delayed reaction is complete: update concentrations, propen, get rid of the reaction in the priority queue
						complete_delay_index = (em.cell_list[k]->cdelay)->soonest_reaction();
						sd.reac_funs[complete_delay_index](em, k, sd, rs, true);
					}
				}
			}
			// transfer data from current time and concentration to record time and concentration
			transfer_to_record(em, ip);
		}
		// check whether done or not
		done = check_done(em, ip);
		
	}
	
	delete [] deltaKArray;
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

int update_deltaK_array(double** deltaKArray, cell* current_cell){
	double current_min = INFINITY;
	double min_index = NUM_REACTIONS;
	for (int i = 0; i < NUM_REACTIONS; i++){
		if ((current_cell->propen)[i] != 0){
			// deltak = (Pk - Tk) / ak (next internal - current_internal) / propensity;
			(*deltaKArray)[i] = ((current_cell->next_internal)[i] - (current_cell->current_internal)[i]) / (current_cell->propen)[i];
			if ((*deltaKArray)[i] < current_min){
				current_min = (*deltaKArray)[i];
				min_index = i;
			}
		}
		else{
			(*deltaKArray)[i] = INFINITY;
		}
	}
	(*deltaKArray)[NUM_REACTIONS] = ((current_cell)->cdelay)->see_soonest() - (current_cell->absolute_time);
	if ((*deltaKArray)[NUM_REACTIONS] < current_min){
		min_index = NUM_REACTIONS;
	}
	return min_index;
}

void update_current_internal(cell* current_cell, double delta){
	// For each reaction: Tk = Tk + Ak * Delta
	for (int i = 0; i < NUM_REACTIONS; i++){
		(current_cell->current_internal)[i] = (current_cell->current_internal)[i] + (current_cell->propen)[i] * delta;
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
}

bool check_done(embryo& em, input_params& ip){
	for (int i = 0; i < ip.num_cells; i++){
		if (em.cell_list[i]->absolute_time >= ip.time_total){
			cout << "First done cell: " << i+ 1 << endl;
			cout << "Absolute time: " << em.cell_list[i]->absolute_time;
			return true;
		}
	}
	return false;
}
