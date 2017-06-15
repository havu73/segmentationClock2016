#include <math.h>
#include "tests.hpp"
#include "io.hpp"
#include "macros.hpp"
#include "feats.hpp"

// Levels of covarriance squared in wildtype in 5 bins of covarriance 
double WT_CVS_BINS [5] = {0.225448775, 0.259325483, 0.278990378, 0.325803876, 0.42486995}; 
double DELTA_CVS_BINS [5] = {};
/*
 * Max_score : 5
 */
double test_wildtype(input_params& ip, embryo& em, features& wtf, int set_index){
	double score = 0;

	// process peaks and troughs, find cell-to-cell correlations
	process_wt_smooth_data(ip, em, wtf, set_index);
	
	// first we have to make sure that the cells are in synchronization 
	// sustained oscillation : 1
	score += check_WT_sustained(ip, wtf);
	
	// check period : 1
	score += check_WT_period(ip, wtf);
	
	// H1, H7 spatial amplitude : 1
	score += check_WT_amplitude(ip, wtf);
	
	// check whether or not the simulation passed condition related to the correlation scores. 
	// If not, stop testing and move on. No need to test noise condition anymore.
	bool corr_pass = true;//check_correlations(wtf);
	
	if (!corr_pass){
		if (ip.verbose){
			cout << "Correlation score less than required: " << wtf.avg_correlation << endl;
		}
		return score;
	}
	// This will be executed if it passes the correlations conditions
		
	// calculate noises, and bins
	process_wt_noise_data(em, wtf, ip);
	if (ip.print_features){
		print_features_data(ip, em, wtf, set_index, WT);
	}
	// average mRNA concentrations tests : 1
	score += check_WT_avg_cons(ip, wtf);
	// intrisic noise test
	score += check_WT_intrinsic (ip, wtf);
	// extrinsic noise test
	score += check_WT_extrinsic (ip, wtf);
	// test that the levels of cv squared increases from posterior to anterior : 1
	score += test_wt_cvs(ip, wtf);

	return score;
}

/*
 * Check that for each bin, the cv squared is within the range of 20% from experimental data.
 * This is to check that cv squared increase from posterior to anterior
 * Max score: 0.2 * 5 bins = 1 score
 */
double test_wt_cvs(input_params& ip, features& wtf){
	double score = 0;
	for (int i = 0; i < CVS_BIN; i++){
		double low = WT_CVS_BINS[i] * 0.6;
		double high = WT_CVS_BINS[i] * 1.4;
		if (wtf.cvs_bin[i] <= high && wtf.cvs_bin[i] >= low){
			score += 1.0 / (double) CVS_BIN;
		}
		else {
			if (ip.verbose){
				cout << "CV squared condition does not pass in wildtype, bin " 
				<< i << " ___ " << wtf.cvs_bin[i] << endl;
			}
		}
	}
	return score;
}

int test_delta(input_params& ip, embryo& em, features& fts, features& wtf, int set_index){
	int score = 0;
	// Calculate the average noise level across slices
	process_delta_noise_data(em, fts, ip);
	// Check that the total noise levels across slices in delta mutant is 40% higher than compared to wildtype.
	// score = 1
	double noise_ratio = fts.slices_tot_noise / wtf.slices_tot_noise;
	if ((noise_ratio >= (DELTA_NOISE_WT_RATIO * 0.8)) && (noise_ratio <= (DELTA_NOISE_WT_RATIO * 1.2))) {
		score += 1;
	}
	else{
		if (ip.verbose){
			cout << "Delta / WT ratio: " << noise_ratio << endl;
		}
	}
	return score;
}


/*
 * Max score = 2
 */
double check_WT_sustained (input_params& ip, features& wtf){
	double score = 0;
	for (int i = 0; i < wtf.num_cells; i ++){
		if (wtf.mid_ptt[KEEPMH1][i] >= LOWER_BOUND_PTT){
			score += (double) 1 / ((double) 3 * wtf.num_cells);
		}
		
		else {
			if (ip.verbose){
				cout << "wtf.mid_ptt[KEEPMH1][" << i << "]: " << wtf.mid_ptt[KEEPMH1][i] << endl;
			}
		}
		
		if (wtf.last_ptt[KEEPMH1][i] >= LOWER_BOUND_PTT){
				score += (double) 1 / ((double) 3 * wtf.num_cells);
		}
		
		else {
			if (ip.verbose){
				cout << "wtf.last_ptt[KEEPMH1][" << i << "]: " << wtf.last_ptt[KEEPMH1][i] << endl;
			}
		}
		
		if ((wtf.mid_ptt[KEEPMH1][i] / wtf.last_ptt[KEEPMH1][i])  <= UPPER_BOUND_MTL){
			score += (double) 1 / ((double) 3 * wtf.num_cells);
		}
		
		else {
			if (ip.verbose){
				cout << "wtf.mid_ptt[KEEPMH1][" << i << "] / wtf.last_ptt[KEEPMH1][" << i << "]: " << wtf.mid_ptt[KEEPMH1][i] / wtf.last_ptt[KEEPMH1][i] << endl;
			}
		}
		
	}	
	return score;
}


double check_WT_amplitude (input_params& ip, features& wtf){
	double score = 0;
	double h1_amplitude = 0;
	double h7_amplitude = 0;
	int num_h1 = 0;
	int num_h7 = 0;
	for (int i = 0; i < wtf.num_cells; i ++){
		if (wtf.avg_amplitude[KEEPMH1][i] > 0){
			h1_amplitude += wtf.avg_amplitude[KEEPMH1][i];
			num_h1 += 1;
		}
		if (wtf.avg_amplitude[KEEPMH7][i] > 0){
			h7_amplitude += wtf.avg_amplitude[KEEPMH7][i];
			num_h7 += 1;
		}
	}
	h1_amplitude = h1_amplitude / (double) num_h1;
	h7_amplitude = h7_amplitude / (double) num_h7;
	
	if (h1_amplitude >= WT_H1_AMP_LOW && h1_amplitude <= WT_H1_AMP_HIGH){
		score += 0.5;
	}
	
	else{
		if (ip.verbose){
			cout << "H1 amplitude: " << h1_amplitude << endl;
		}
	}
	
	if (h7_amplitude >= WT_H7_AMP_LOW && h7_amplitude <= WT_H7_AMP_HIGH){
		score += 0.5;
	}
	
	else{
		if (ip.verbose){
			cout << "H7 amplitude: " << h7_amplitude << endl;
		}
	}
	
	return score;
}

/*
 * Check amplitude within range, score right now is 0
 */
double check_wt_alternative_amplitude(input_params& ip, features& wtf){
	double score = 0;
	//check her1
	if (wtf.alternative_amplitude[KEEPMH1] >= WT_H1_ALT_AMP_LOW && wtf.alternative_amplitude[KEEPMH1] <= WT_H1_ALT_AMP_HIGH){
		score += 0;
	}
	else{
		if (ip.verbose){
			cout << "Alternative amplitude her1 failed: " << wtf.alternative_amplitude[KEEPMH1];
		}
	}
	//check her7
	if (wtf.alternative_amplitude[KEEPMH7] >= WT_H7_ALT_AMP_LOW && wtf.alternative_amplitude[KEEPMH7] <= WT_H7_ALT_AMP_HIGH){
		score += 0;
	}
	else{
		if (ip.verbose){
			cout << "Alternative amplitude her7 failed: " << wtf.alternative_amplitude[KEEPMH7];
		}
	}
	return score;
}
/*
 * Max_score : 1
 */
double check_WT_period (input_params& ip, features& wtf){
	double score = 0;
	
	for (int i = 0; i < wtf.num_cells; i++){
		double h1_per = wtf.avg_period[KEEPMH1][i];
		double h7_per = wtf.avg_period[KEEPMH7][i];
		if (h1_per <= WT_HIGH_PER && h1_per >= WT_LOW_PER){
			score += (0.5 / wtf.num_cells);
		}
		
		else{
			if (ip.verbose){
				cout << "wtf.avg_period[KEEPMH1][" << i << "]: " << wtf.avg_period[KEEPMH1][i] << endl;
			}
		}
		
		if (h7_per <= WT_HIGH_PER && h7_per >= WT_LOW_PER){
			score += (0.5 / wtf.num_cells);
		}
		
		else{
			if (ip.verbose){
				cout << "wtf.avg_period[KEEPMH7][" << i << "] :" << wtf.avg_period[KEEPMH7][i] << endl;
			}
		}
		
	}
	
	return score;
}

int check_delta_amplitude (features& fts){
	double h1_amplitude = 0;
	double h7_amplitude = 0;
	int num_h1 = 0;
	int num_h7 = 0; 
	for (int i = 0; i < fts.num_cells; i++){
		if (fts.avg_amplitude[KEEPMH1][i] > 0){
			h1_amplitude += fts.avg_amplitude[KEEPMH1][i];
			num_h1 += 1;
		}
		if (fts.avg_amplitude[KEEPMH7][i] > 0){
			h7_amplitude += fts.avg_amplitude[KEEPMH7][i];
			num_h7 += 1;
		}
	}
	h1_amplitude = h1_amplitude / (double) num_h1;
	h7_amplitude = h7_amplitude / (double) num_h7;
	
	if (h1_amplitude < DELTA_H1_AMP_LOW){
		return 0;
	}
	if (h1_amplitude > DELTA_H1_AMP_HIGH){
		return 0;
	}
	if (h7_amplitude < DELTA_H7_AMP_LOW){
		return 0;
	}
	if (h7_amplitude > DELTA_H7_AMP_HIGH){
		return 0;
	}
	return 1;
}

double check_WT_avg_cons(input_params& ip, features& wtf){
	double score = 0;
	// Bin 1
	if (wtf.avg_cons[0] >=  LOW_WT_HER_BIN1 && wtf.avg_cons[0] <= UP_WT_HER_BIN1){
		score += 1.0 / (double) wtf.num_bin;
	}
	
	else{
		if (ip.verbose){
			cout << "Bin 0 avg_cons: " << wtf.avg_cons[0] << endl;
		}
	}
	
	
	// Bin 2
	if (wtf.avg_cons[1] >= LOW_WT_HER_BIN2 && wtf.avg_cons[1] <= UP_WT_HER_BIN2){
		score += 1.0 / (double) wtf.num_bin;
	}
	
	else{
		if (ip.verbose){
			cout << "Bin 1 avg_cons: " << wtf.avg_cons[1] << endl;  
		}
	}
	
	
	// Bin 3
	if (wtf.avg_cons[2] >= LOW_WT_HER_BIN3 && wtf.avg_cons[2] <= UP_WT_HER_BIN3){
		score += 1.0 / (double) wtf.num_bin;
	}
	
	else {
		if (ip.verbose){
			cout << "Bin 2 avg_cons: " << wtf.avg_cons[2] << endl;
		}
	}
	
	return score;
}

int check_delta_avg_cons(features& fts){
	// Bin 1
	if (fts.avg_cons[0] < LOW_DELTA_HER_BIN1){
		return 0;
	}
	if (fts.avg_cons[0] > UP_DELTA_HER_BIN1){
		return 0;
	}
	// Bin 2
	if (fts.avg_cons[1] < LOW_DELTA_HER_BIN2){
		return 0;
	}
	if (fts.avg_cons[1] > UP_DELTA_HER_BIN2){
		return 0;
	}
	// Bin 3
	if (fts.avg_cons[2] < LOW_DELTA_HER_BIN3){
		return 0;
	}
	if (fts.avg_cons[2] > UP_DELTA_HER_BIN3){
		return 0;
	}
	return 1;
}
double check_WT_intrinsic(input_params& ip, features& wtf){
	double score = 0;
	// Bin 1
	if (wtf.intrinsic[0] >= LOW_WT_IN_NOISE_BIN1 && wtf.intrinsic[0] <= UP_WT_IN_NOISE_BIN1){
		score += 1.0 / (double) wtf.num_bin;
	}
	
	else {
		if (ip.verbose){
			cout << "Intrinsic noise bin 0: " << wtf.intrinsic[0] << endl;
		}
	}
	
	// Bin 2
	if (wtf.intrinsic[1] >= LOW_WT_IN_NOISE_BIN2 && wtf.intrinsic[1] <= UP_WT_IN_NOISE_BIN2){
		score += 1.0 / (double) wtf.num_bin;
	}
	
	else {
		if (ip.verbose){
			cout << "Intrinsic noise bin 1: " << wtf.intrinsic[1] << endl;
		}
	}
	
	// Bin 3
	if (wtf.extrinsic[2] >= LOW_WT_IN_NOISE_BIN3 && wtf.intrinsic[2] <= UP_WT_IN_NOISE_BIN3){
		score += 1.0 / (double) wtf.num_bin;
	}
	
	else{
		if (ip.verbose){
			cout << "Intrinsic noise bin 2: " << wtf.intrinsic[2] << endl;
		}
	}
	
	return score;
}

int check_delta_intrinsic(features& fts){
	// Bin 1
	if (fts.intrinsic[0] < LOW_DELTA_IN_NOISE_BIN1){
		return 0;
	}
	if (fts.intrinsic[0] > UP_DELTA_IN_NOISE_BIN1){
		return 0;
	}
	// Bin 2
	if (fts.intrinsic[1] < LOW_DELTA_IN_NOISE_BIN2){
		return 0;
	}
	if (fts.intrinsic[1] > UP_DELTA_IN_NOISE_BIN2){
		return 0;
	}
	// Bin 3
	if (fts.intrinsic[2] < LOW_DELTA_IN_NOISE_BIN3){
		return 0;
	}
	if (fts.intrinsic[2] > UP_DELTA_IN_NOISE_BIN3){
		return 0;
	}
	return 1;
}

double check_WT_extrinsic (input_params& ip, features& wtf){
	double score = 0;
	// Bin 1
	if (wtf.extrinsic[0] >= LOW_WT_EX_NOISE_BIN1 && wtf.extrinsic[0] <= UP_WT_EX_NOISE_BIN1){
		score += 1.0 / (double) wtf.num_bin;
	}
	
	else{
		if (ip.verbose){
			cout << "Extrinsic noise bin 0: " << wtf.extrinsic[0] << endl;
		}
	}
	
	// Bin 2
	if (wtf.extrinsic[1] >= LOW_WT_EX_NOISE_BIN2 && wtf.extrinsic[1] <= UP_WT_EX_NOISE_BIN2){
		score += 1.0 / (double) wtf.num_bin;
	}
	
	else{
		if (ip.verbose){
			cout << "Extrinsic noise bin 1: " << wtf.extrinsic[1] << endl;
		}
	}
	
	// Bin 3
	if (wtf.extrinsic[2] >= LOW_WT_EX_NOISE_BIN3 && wtf.extrinsic[2] <= UP_WT_EX_NOISE_BIN3){
		score += 1.0 / (double) wtf.num_bin;
	}
	
	else{
		if (ip.verbose){
			cout << "Extrinsic noise bin 2: " << wtf.extrinsic[2] << endl;
		}
	}
	
	return score;
}

int check_delta_extrinsic (features& fts){
	// Bin 1
	if (fts.extrinsic[0] < LOW_DELTA_EX_NOISE_BIN1){
		return 0;
	}
	if (fts.extrinsic[0] > UP_DELTA_EX_NOISE_BIN1){
		return 0;
	}
	// Bin 2
	if (fts.extrinsic[1] < LOW_DELTA_EX_NOISE_BIN2){
		return 0;
	}
	if (fts.extrinsic[1] > UP_DELTA_EX_NOISE_BIN2){
		return 0;
	}
	// Bin 3
	if (fts.extrinsic[2] < LOW_DELTA_EX_NOISE_BIN3){
		return 0;
	}
	if (fts.extrinsic[2] > UP_DELTA_EX_NOISE_BIN3){
		return 0;
	}
	return 1;
}



