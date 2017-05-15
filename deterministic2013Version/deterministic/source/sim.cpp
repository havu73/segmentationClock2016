#include <cmath>
#include "sim.hpp"
#include "init.hpp"
#include "io.hpp"
#include "tests.hpp"
#include "calculation.hpp"

extern terminal* term; // Declared in init.cpp


void simulate_all_params_sets(input_params& ip, parameters& pr){
	rates rs;
	for (int i = 0; i < ip.num_sets; i++){
		//Create set directories
		if (ip.print_cons){
			create_set_directory(i , ip);
		}
		// process params from parameters to rates
		process_rates(rs, pr, ip, i);
		// simulate each set
		int set_score = simulate_one_param_set(ip, rs, i);
		// reset rates
		rs.reset();
	}
}

void process_rates(rates& rs, parameters& pr, input_params& ip, int set_index){
	for (int i = 0; i < NUM_CONSTANT_RATES; i ++){
		rs.data[i] = pr.data[set_index][i];
	}
	for (int i = NMH1; i <= NPD; i++){
		rs.data[i] = (double)((int)(pr.data[set_index][i] / ip.step_size));
	}
	rs.data[CRITPDH1] = pr.data[set_index][KDG1PN] / pr.data[set_index][KAG1PN];
	rs.data[CRITPH11H1] = pr.data[set_index][KDG1PH11] / pr.data[set_index][KAG1PH11];
	rs.data[CRITPDH7] = pr.data[set_index][KDG7PN] / pr.data[set_index][KAG7PN];
	rs.data[CRITPH11H7] = pr.data[set_index][KDG7PH11] / pr.data[set_index][KAG7PH11];
	rs.data[CRITPH11D] = pr.data[set_index][KDGDPH11] / pr.data[set_index][KAGDPH11];
}

int simulate_one_param_set (input_params& ip, rates& rs, int set_index){
	double set_score;
	// change the rates of delay 
	// find_sim_steps (rs, ip);
	// find the size of record_steps in con_levels
	int steps = (int) (ip.time_total / ip.step_size);
	cout << "Steps: " << steps << endl;
	// declare con_levels
	con_levels cons (steps);
	for (int i = 0; i < NUM_MUTANTS; i++){
		set_score = (double)MAX_SCORE - (double)simulate_mutant(ip, rs, cons, set_index, i);
		if (ip.pipe_out){
			write_pipe( &set_score, ip);
		}
		cons.reset();
	}
	return set_score;
}

int simulate_mutant(input_params& ip, rates& rs, con_levels& cons, int set_index, int mutant_index){
	int mutant_score = 0;
	bool found_nan = core_simulation(ip, rs, cons);
	if (!found_nan){
		mutant_score += 1;
		// print 
		if (ip.print_cons){
			print_concentrations(ip, set_index, mutant_index, cons);
		}
		// test
		mutant_score += test_wildtype(cons, ip);
	}
	return mutant_score; 
}

bool core_simulation (input_params& ip, rates& rs, con_levels& cons){
	bool found_nan;
	for (int i = 1; i < cons.num_sim_steps; i++){
		// calculation
		calculateMH1(ip, cons, rs, i);
		calculateMH7(ip, cons, rs, i);
		calculateMD(ip, cons, rs, i);
		calculatePH1(ip, cons, rs, i);
		calculatePH7(ip, cons, rs, i);
		calculatePD(ip, cons, rs, i);
		calculatePH11(ip, cons, rs, i);
		calculatePH17(ip, cons, rs, i);
		calculatePH77(ip, cons, rs, i);
		// check for recording and check for nan
		if (i % ip.record_gran == 0){
			found_nan = check_for_nan (cons, i);
			if (found_nan){
				return true;
			}
		}
	}
	return false;
}

bool check_for_nan (con_levels& cons, int record_index){
	for (int i = 0; i < NUM_RECORD_STATES; i++){
		if (isnan(cons.sim_data[i][record_index])){
			cout << term->red << "Found nan..... " << term->reset << endl;
			return  true;
		}
	}
	return false;
}


void find_sim_steps (rates& rs, input_params& ip){
	int steps;
	for (int i = NMH1; i <= NMD; i++){
		steps = (int) (rs.data[i] / ip.step_size);
		rs.data[i] = (double) steps;
	}
}
