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
	// create output directory if the user specified so
	create_output_directory (ip);

	// create debug files if needed here
	
	for (int i = 0; i < ip.num_sets; i++){
		// process the rates so that we can run the model
		process_rates(rs, pr, i);
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
	features wtf(ip.num_cells, ip.num_bin); // structures to store wildtype features, because this will be many times throughout the simulation
	//cout << term->blue << "Simulating set " << term->reset << set_index << " . . ." << endl;
	for (int i = 0; i < ip.num_mutants; i++){
		// create mutant directory
		if (ip.print_cons || ip.print_features){
			create_mutant_directory(set_index, ip.mutants[i], ip);
		}
		// change rates based on mutant
		change_rates_based_on_mutants(rs, i);
		if (ip.mutants[i] != DAPT_MUTANT){
			// simulate_mutant, then test conditions of such mutants, then report back the score
			score += simulate_mutant(ip, sd, rs, em, wtf, set_index, ip.mutants[i]);
		}
		else if (ip.mutants[i] == DAPT_MUTANT){
			score += simulate_dapt_mutant(ip, sd, rs, em, wtf, set_index, ip.mutants[i]);
		}
		// revert_rates
	}
	// real score to pass into SRES
	sresScore = (double) ip.max_cond_score - score;
	// print score to pipe if necessary
	if (ip.piping){
		write_pipe(&sresScore, ip);
	}
	if (ip.verbose){
		cout << "Set " << set_index << " received score: " << sresScore << endl;
	}
	return sresScore;
}

void change_rates_based_on_mutants(rates& rs, int mutant_index){
	if (mutant_index == DELTA_MUTANT) {
		for (int i = 0; i < rs.num_cells; i ++){
			rs.data[i][PSD] = 0;
		}
	}
}

/*
 * 1, Simulate the mutant
 * 2, Test mutant conditions
 * 3, Report mutant score
 */
double simulate_mutant(input_params& ip, sim_data& sd, rates& rs, embryo& em, features& wtf, int set_index, int mutant_index){
	double mutant_score = 0;
	// reset embryo will do the following for each cell:
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
	features fts(ip.num_cells, ip.num_bin);
	if (mutant_index == WT && (!out_of_bound)){
		mutant_score = test_wildtype(ip, em, wtf, set_index);
	}	
	else if (mutant_index == DELTA_MUTANT && (!out_of_bound)){
		mutant_score = test_delta(ip, em, fts, wtf, set_index);
	}
	// print data if necessary
	if (ip.print_cons){
		print_concentrations(ip, set_index, mutant_index, em);
	}
	return mutant_score;
}

double simulate_dapt_mutant(input_params& ip, sim_data& sd, rates& rs, embryo& em, features& wtf, int set_index, int mutant_index){
	double mutant_score = 0;
	// reset embryo will do the following for each cell:
	// cons: current_cons to all 0; cons_record to all 0
	// propen[NUM_REACTIONS] all to 0
	// next_internal [NUM_REACTIONS]  all to 0
	// current_internal [NUM_REACTIONS] all to 0
	// cdelay empty
	// internal_time = 0
	em.reset();
	core_simulation_dapt(ip, sd, rs, em);
	return mutant_score;
}

bool core_simulation_dapt(input_params& ip, sim_data& sd, rates& rs, embryo& em){
	// Core simulation for dapt
	// update initial conditions : all gene has current concentrations to 2
	update_initial_concentrations(em);
	// calculate propensity for each reaction
	calculate_initial_propensity(sd, rs, em);
	calculate_initial_next_internal(ip, em);
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
	
	// an array of psd rates from cells in the embryo, used to store the original value whild delta is knocked down.
	double* revert_rates = new double[ip.num_cells];
	
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
					delete[] revert_rates;
					return out_of_bound;
				}
				bound_check_index = 0;
				// 2. Check that time for simulation is up
				check_done_time += (int) MINUTE_PER_SLICE;
				if (check_done_time == (int)ip.time_total){ //time is up
					done = true;
				}
				// if it is the start time of the dapt mutant, reverse all the values of psd to be 0
				if(check_done_time == (int) INDUCTION_DAPT_TIME){
					// get every psd rates to be 0
					for (int i = 0; i < em.num_cells; i++){
						revert_rates[i] = rs.data[i][PSD];
						rs.data[i][PSD] = 0;
					}
				}
				else if (check_done_time == (int) WITHDRAW_DAPT_TIME){
					// revert the psd rates
					for (int i = 0; i < em.num_cells; i++){
						rs.data[i][PSD] = revert_rates[i];
					}
				}
			}
		}
	}
	delete[] revert_rates;
	return false;
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
	// Find the gradient to increase/decrease protein delay rates after each time interval from posterior to anterior
	double* gradient_per_step = new double [em.num_cells]; 
	calculate_gradient_per_step(ip, rs, gradient_per_step);
	// number indicating the time step between two intervals to check for states_out_of_bound
	int bound_check_index = 0;
	
	// number of gradient degree steps the protein synthesis delay rates should be, given the current time
	// First it has to be 0 until the end of ABSOLUTE_RATE_TIME
	// Then, after every MINUTE_PER_SLICE minutes, it will be incremented by one
	int gradient_degree = 0;
	
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
				// 3. Increase the gradient if necessary
				if (check_done_time >= ABSOLUTE_RATE_TIME && !done){
					gradient_degree += 1;
					calculate_gradient_delay_rates(gradient_degree, rs, gradient_per_step);
				}
			}
		}
	}
	delete [] gradient_per_step;
	return false;
}

void calculate_gradient_delay_rates(int gradient_degree, rates& rs, double* gradient_per_step) {
	for (int i = 0; i < rs.num_cells; i ++){
		rs.data[i][NPH1] = rs.data[i][NPH1] * (1 + gradient_degree * gradient_per_step[i]);
		rs.data[i][NPD] = rs.data[i][NPD] * (1 + gradient_degree * gradient_per_step[i]);
	}
}

void calculate_gradient_per_step (input_params& ip, rates& rs, double* gradient_per_step){
	int num_steps = int ((ip.time_total - ABSOLUTE_RATE_TIME) / MINUTE_PER_SLICE);
	for (int i = 0 ; i  < rs.num_cells; i ++){
		gradient_per_step[i] = (rs.data[i][GRADIENT] - 1) / num_steps;
	}
}

bool states_out_of_bound(embryo& em){
	for (int i = 0; i < em.num_cells; i ++){
		if (((em.cell_list[i])->current_cons)[MH1] > MRNA_LIM){ return true; }
		if (((em.cell_list[i])->current_cons)[MH7] > MRNA_LIM){ return true; }
		if (((em.cell_list[i])->current_cons)[MD] > MRNA_LIM){ return true; }
		if (((em.cell_list[i])->current_cons)[PH1] > PROTEIN_LIM){ return true; }
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
