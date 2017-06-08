#include <cmath>
#include "sim.hpp"
#include "debug.hpp"
#include "init.hpp"
#include "io.hpp"
#include "randDist.hpp"
#include "tests.hpp"
#include "print_debug.hpp"
#include <ctime>
extern terminal* term; // Declared in init.cpp


void simulate_all_params_sets(input_params& ip, parameters& pr){
	embryo em(ip);
	rates rs (ip);
	// declare reactions and propensities: lists of function calls to use throughout the simulation
	sim_data sd;
	init_propensities(sd);
	init_reactions(sd);
	//init_seeds(ip);
	
	// create debug file if necessary
	/*
	if (ip.print_debug){
		open_debug_steam(ip);
	}
	*/
	/*
	if (ip.print_random){
		open_random_steam(ip);
	}
	*/
	for (int i = 0; i < ip.num_sets; i++){
		// process the rates so that we can run the model
		process_rates(rs, pr, i);
		cout << "Done processing rates" << endl;
		//Create set directories
		if (ip.print_cons || ip.print_features || ip.print_rates){
			create_set_directory(i , ip);
		}
		// if users sepcify that they want to print out the rates
		if (ip.print_rates) { 
			print_perturb_rates(ip, rs, pr, i);
		}
		
		// simulate each set
		simulate_one_param_set(ip, sd, rs, i, em);
		// reset rates after each set
		rs.reset();
	}
}


void process_rates(rates& rs, parameters& pr, int set_index){
	double low_perturb = 1 - (pr.data[set_index][PERTURB] / (double) 100);
	double high_perturb = 1 + (pr.data[set_index][PERTURB] / (double) 100);
	for (int i = 0; i < rs.num_cells; i++){ // for each cell
		for (int j = 0; j < NUM_NETWORK_RATES; j ++) { // for each rate related to our system
			double pert_factor = doubleRand(low_perturb, high_perturb);  // geenrate a random number in the allowed perturb range
			rs.perturb_rates[i][j] = pert_factor; 
			rs.data[i][j] = pr.data[set_index][j] * pert_factor; // set the rate for this cell based on the perturb rates we specified
		}
	}
}

double simulate_one_param_set(input_params& ip, sim_data& sd, rates& rs, int set_index, embryo& em){
	double score = 0; 
	double sresScore;
	//cout << term->blue << "Simulating set " << term->reset << set_index << " . . ." << endl;
	for (int i = 0; i < ip.num_mutants; i++){
		// create mutant directory
		if (ip.print_cons || ip.print_features){
			create_mutant_directory(set_index, ip.mutants[i], ip);
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
	if (ip.verbose){
		cout << "Set " << set_index << " received score: " << sresScore << endl;
	}
	return sresScore;
}

double simulate_mutant(input_params& ip, sim_data& sd, rates& rs, embryo& em, int set_index, int mutant_index){
	double mutant_score = 0;
	// reset embryo : for each cell:
	// cons: current_cons to all 0; cons_record to all 0
	// propen[NUM_REACTIONS] all to 0
	// next_internal [NUM_REACTIONS]  all to 0
	// current_internal [NUM_REACTIONS] all to 0
	// cdelay empty
	// internal_time = 0

	em.reset();
	// core simulation
	bool out_of_bound = core_simulation(ip, sd, rs, em); // simulate, and return whether or not the simulations results in moments where the states levels are above our accepted bounds
	
	// test conditions
	features wtf(ip.num_cells, ip.num_bin);
	if (mutant_index == WT && (!out_of_bound)){
		mutant_score = test_wildtype(ip, em, wtf, set_index);
	}	
	// print data if necessary
	if (ip.print_cons){
		print_concentrations(ip, set_index, mutant_index, em);
	}
	return mutant_score;
}

/*
 * Return whether or not the states are found out of bound.
 */
bool core_simulation(input_params& ip, sim_data& sd, rates& rs, embryo& em){
	// update initial conditions : all gene has current concentrations to 2
	update_initial_concentrations(em);
	// calculate propensity for each reaction
	calculate_initial_propensity(sd, rs, em);
	calculate_initial_next_internal(ip, em);
	/*
	if (ip.print_debug){
		print_initial_simulation(ip, em);
	}
	*/
	// time keeping variables
	bool done = false;
		
	// next reaction to denote parameters about what reaction in what cell and whether it is a complete delay reaction
	// that is firing next
	next_reaction nr;
	
	// Keep the record at time 0
	transfer_to_record(em, ip);
	// Initialize the time that the record is taken
	double current_record_time = 0;
	int record_per_check = int((double)MINUTE_PER_SLICE / ip.record_granularity);
	int check_done_time = 0;
	
	// number indicating the time step between two intervals to check for states_out_of_bound
	int bound_check_index = 0;
	// enter simulation loop
	while (!done){
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
			find_next_firing(ip, em.cell_list[nr.cell_index], nr.reaction_index);
		}
		/*
		if (ip.print_debug){
			print_one_round_simulation(em, ip, nr);
		}
		*/
		
		if ((em.absolute_time - current_record_time) >= ip.record_granularity){
			// ip.record_granularity time has passed, now transfer record and update new current_record_time
			// transfer data from current time and concentration to record time and concentration
			transfer_to_record(em, ip);
			current_record_time += ip.record_granularity;
			bound_check_index += 1;
			
			if (bound_check_index == record_per_check){ // MINUTE_PER_SLICE minutes have passed
				// 1. Check that the states are not out of bound every MINUTE_PER_SLICE minutes
				bool out_of_bound = states_out_of_bound(em);
				if (out_of_bound){
					return out_of_bound;
				}
				bound_check_index = 0;
				// 2. Check that time for simulation is up
				check_done_time += (int) MINUTE_PER_SLICE;
				if (check_done_time == (int)ip.time_total){ //time is up
					done = true;
				}
			}
		}
	}
	return false;
}

bool states_out_of_bound(embryo& em){
	for (int i = 0; i < em.num_cells; i ++){
		if (((em.cell_list[i])->current_cons)[MH1] > MRNA_LIM){ return true; }
		if (((em.cell_list[i])->current_cons)[MH7] > MRNA_LIM){ return true; }
		if (((em.cell_list[i])->current_cons)[MD] > MRNA_LIM){ return true; }
		if (((em.cell_list[i])->current_cons)[PH1] > PROTEIN_LIM){ return true; }
		if (((em.cell_list[i])->current_cons)[PH7] > PROTEIN_LIM){ return true; }
		if (((em.cell_list[i])->current_cons)[PD] > PROTEIN_LIM){ return true; }
	}
	return false;
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

void calculate_initial_next_internal(input_params& ip, embryo& em){
	for (int i = 0; i < em.num_cells; i++){
		for (int j = 0; j < NUM_REACTIONS; j++){
			double random_number = pk_dist();
			((em.cell_list[i])->next_internal)[j] = random_number;
			/*
			if (ip.print_random){
				print_random_number(ip, &random_number);
			}
			*/
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

void find_next_firing(input_params& ip, cell* current_cell, int reaction_index){
	//double random_number = pk_dist();
	(current_cell->next_internal)[reaction_index] = (current_cell->next_internal)[reaction_index] + pk_dist(); 
	/*
	if (ip.print_random){
		print_random_number(ip, &random_number);
	}
	*/
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
