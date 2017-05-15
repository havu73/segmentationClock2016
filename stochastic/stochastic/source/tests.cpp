#include <math.h>
#include "tests.hpp"
#include "io.hpp"
#include "macros.hpp"
/*
 * Max_score : 5
 */
double test_wildtype(input_params& ip, embryo& em, features& wtf, int set_index){
	double score = 0;
	/*
	// check mRNA boundaries
	bool mRNA_in_bound = true; //check_mRNA_boundaries(em);
	if (! mRNA_in_bound){
		return 0;	
	}
	* */

	// process peaks and troughs
	process_smooth_data(ip, em, wtf, set_index);

	// sustained oscillation : 1
	score += check_WT_sustained(ip, wtf);
	
	// check period : 1
	score += check_WT_period(ip, wtf);
	
	// H1, H7 spatial amplitude : 1
	score += check_WT_amplitude(ip, wtf);
	
	/*
	// Before we go on to check noise levels, we need to check the correlations among cells first
	if (!check_correlations(wtf)){
		return score;
	}
	*/
	// calculate noises
	process_noise_data(em, wtf, ip);
	if (ip.print_features){
		print_features_data(ip, em, wtf, set_index, WT);
	}
	// average mRNA concentrations tests : 1
	score += check_WT_avg_cons(ip, wtf);
	// intrisic noise test
	score += check_WT_intrinsic (ip, wtf);
	// extrinsic noise test
	score += check_WT_extrinsic (ip, wtf);
	
	return score;
}

int test_delta(input_params& ip, embryo& em, features& fts, int set_index){
	int score = 0;
	// process peaks and troughs
	process_smooth_data(ip, em, fts, set_index);
	// H1, H7 spatial amplitudes are in a specific range
	score += check_delta_amplitude (fts);
	// calculate noises
	process_noise_data(em, fts, ip);
	// average mRNA concentrations tests
	score += check_delta_avg_cons (fts);
	// intrinsic noise test
	score += check_delta_intrinsic (fts);
	// extrinsic noise test
	score += check_delta_extrinsic (fts);
	return score;
}

bool check_mRNA_boundaries (embryo& em){
	int num_steps = em.time_record->size();
	for (int i = 0; i < em.num_cells; i++){
		vector<int> * cch1 = (em.cell_list[i])->cons_record[KEEPMH1];
		vector<int> * cch7 = (em.cell_list[i])->cons_record[KEEPMH7];
		if (cch1->size() != num_steps || cch7->size() != num_steps){
			cout << "Encounter a problem with size of cons and time record: num_steps: " << num_steps << "____ cons: " << cch1->size() << "____" << cch7->size() << endl;
			return false;
		}
		for (int j = 0; j < num_steps; j++){
			if ((cch1->at(j) + cch7->at(j)) > UPPER_BOUND_HER){
				return false;
			}
		}
	}
	return true;
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
	if (wtf.extrinsic[1] >= LOW_WT_EX_NOISE_BIN2 && wtf.extrinsic[0] <= UP_WT_EX_NOISE_BIN2){
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

void process_smooth_data(input_params& ip, embryo& em, features& fts, int set_index){
	// find the good capcaity for peaks_troughs
	int capacity = ((em.cell_list[0])->cons_record[0])->size() - WINDOW_SIZE * 2;
	//int capacity = em.time_record->size() - WINDOW_SIZE * 2;
	peak_trough pt (capacity, ip.num_cells);
	for (int j = 0; j < NUM_KEEP_STATES; j ++){
		smooth_data(em, pt, j);
		for (int i = 0; i < em.num_cells; i ++){
			find_peaks_and_troughs(em, pt,fts, i, j);
		}
		if (ip.print_cons){
			print_smooth_data (ip, em, pt, j, set_index);
		}
		cell_to_cell_correlations(pt, fts, j);
		pt.reset();
	}
}


void smooth_data(embryo& em, peak_trough& pt, int con_index){
	int sum [pt.num_cells];
	memset(sum, 0, sizeof(int) * pt.num_cells); 
	
	int avg_divisor = WINDOW_SIZE * 2 + 1;
	for (int j = 0; j < pt.num_cells; j++){
		vector<int> * current_vector = (em.cell_list[j]->cons_record)[con_index];
		for (int i = 0; i < avg_divisor; i ++){
			sum[j] += current_vector->at(i);
		}
		pt.smooth_cons[j][0] = double(sum[j]) / double(avg_divisor);
	}
	int remove_index = 0;
	int add_index = WINDOW_SIZE * 2 + 1;
	
	for (int i = 1; i < (pt.capacity); i ++){
		for (int j = 0; j < pt.num_cells; j++){
			vector<int> * current_vector = (em.cell_list[j]->cons_record)[con_index];
			sum[j] -= current_vector->at(remove_index);
			sum[j] += current_vector->at(add_index);
			(pt.smooth_cons[j])[i] = double(sum[j]) / double(avg_divisor);
		}
		remove_index += 1;
		add_index += 1;
	}
}

void find_peaks_and_troughs(embryo& em, peak_trough& pt, features& fts, int cell_index, int con_index){
	// amplitudes and ratios
	int max_cons = pt.capacity;
	int i = 1; 
	bool is_peak = false;
	bool is_trough = false;
	double avg_peak = 0;
	double avg_trough = 0;
	double avg_period = 0;
	while (i < (max_cons - 1)){
		// potential peak
		if (pt.smooth_cons[cell_index][i] >= pt.smooth_cons[cell_index][i - 1] && pt.smooth_cons[cell_index][i] >= pt.smooth_cons[cell_index][i + 1]){
			is_peak = check_peak(pt, cell_index, i);
		}
		// potential trough
		if (pt.smooth_cons[cell_index][i] <= pt.smooth_cons[cell_index][i - 1] && pt.smooth_cons[cell_index][i] <= pt.smooth_cons[cell_index][i + 1]){
			is_trough = check_trough(pt, cell_index, i);
		}
		if (is_peak){
			((pt.peaks)[cell_index])->push_back(i);
			avg_peak += pt.smooth_cons[cell_index][i];
			i += WING_CHECK_SIZE_LEVEL_1 + WING_CHECK_SIZE_LEVEL_2 + WING_CHECK_SIZE_LEVEL_3;
			is_peak = false;
		}
		if (is_trough){
			((pt.troughs)[cell_index])->push_back(i);
			avg_trough += pt.smooth_cons[cell_index][i];
			i += WING_CHECK_SIZE_LEVEL_1 + WING_CHECK_SIZE_LEVEL_2 + WING_CHECK_SIZE_LEVEL_3;
			is_trough = false;
		}
		i += 1;
	}
	
	// if we have too few peaks and troughs, set to 0
	// so that later on, I can choose what to take the average of when calculating
	// average of all cells 
	if ((pt.peaks[cell_index])->size() < low_num_pt() || (pt.troughs[cell_index])->size() < low_num_pt()){
		fts.avg_amplitude[con_index][cell_index] = 0;
		fts.avg_period[con_index][cell_index] = 0;
		fts.mid_ptt[con_index][cell_index] = 0;
		fts.last_ptt[con_index][cell_index] = 0;
	}
	else {
		// calculate the amplitude of the concentrations of this cell
		// amplitude 
		avg_peak = avg_peak / ((pt.peaks[cell_index])->size());
		avg_trough = avg_trough / ((pt.troughs[cell_index])->size());
		fts.avg_amplitude[con_index][cell_index] = avg_peak - avg_trough;
		// period
		for (int i = 1; i < pt.peaks[cell_index]->size(); i++){
			int current_index = (pt.peaks[cell_index])->at(i) + WINDOW_SIZE;
			int prev_index = (pt.peaks[cell_index])->at(i - 1) + WINDOW_SIZE;
			avg_period += ((em.time_record)->at(current_index) - (em.time_record)->at(prev_index));
		}
		avg_period = avg_period / ((pt.peaks[cell_index])->size() - 1);
		fts.avg_period[con_index][cell_index] = avg_period;
		//calculate peak_to_trough
		calculate_p2t(pt, fts, con_index, cell_index);
	}
}

bool check_peak(peak_trough& pt, int cell_index, int current_index){
	double current_cons = pt.smooth_cons[cell_index][current_index];
	double* smooth_record = pt.smooth_cons[cell_index];
	vector<int> * peaks = pt.peaks[cell_index];
	vector<int> * troughs = pt.troughs[cell_index];
	int upper_index_1 = current_index + WING_CHECK_SIZE_LEVEL_1;
	int lower_index_1 = current_index - WING_CHECK_SIZE_LEVEL_1;
	int upper_index_2 = upper_index_1 + WING_CHECK_SIZE_LEVEL_2;
	int lower_index_2 = lower_index_1 - WING_CHECK_SIZE_LEVEL_2;
	int upper_index_3 = upper_index_2 + WING_CHECK_SIZE_LEVEL_3;
	int lower_index_3 = lower_index_2 - WING_CHECK_SIZE_LEVEL_3;
	
	if (troughs->size() > 0){
		if (peaks->size() > 0 && peaks->back() > troughs->back()){
			//No two peaks without any trough in between
			//cout << "Check peak: 2 peaks in a row" << endl;
			return false;
		}
		if (smooth_record[troughs->back()] >= current_cons){
			// previous trough cannot be higher than this peak
			//cout << "check peak: previous trough is >= this peak" << endl;
			return false;
		}
	}else{ // no troughs yet
		if (peaks->size() > 0){
			//cout << "Check peak: no troughs yet but this is the second peaks found" << endl;
			// No two peaks without any trough in between
			return false;
		}
	}
	// left wing check_level_1 : if any of the left wing points is greater than the current point: return false immediately
	for (int i = current_index - 1; (i >= lower_index_1 && i >= 0); i --){
		if (smooth_record[i] > current_cons){
			//cout << "Check_peak: Left wing level 1 failed: " << i << endl;
			return false;
		}
	}
	// right wing check_level_1; if any of the right wing points is greater than the current point: return false immediately
	for (int i = current_index + 1; (i <= upper_index_1 && i < pt.capacity); i++){
		if (smooth_record[i] > current_cons){
			//cout << "Check_peak: Right wing level 1 failed: " << i << endl;
			return false;
		}
	}
	// left_wing_check_level_2: if any of the left wing points is greater than the current point: return false immediately
	for (int i = lower_index_1 - 1; (i >= lower_index_2 && i >= 0); i --){
		if (smooth_record[i] > current_cons){
			//cout << "Check_peak: Left wing level 2 failed: " << i << endl;
			return false;
		}
	}
	// right wing check_level_2: if any of the right wing points is greater than the current point: return false immediately
	for (int i = upper_index_1 + 1; (i <= upper_index_2 && i < pt.capacity); i++){
		if (smooth_record[i] > current_cons){
			//cout << "Check_peak: Right wing level 2 failed: " << i << endl;
			return false;
		}
	}
	// left_wing_check_level_3: if any of the left wing points is greater than the current point: return false immediately
	for (int i = lower_index_2 - 1; (i >= lower_index_3 && i >= 0); i--){
		if (smooth_record[i] > current_cons){
			//cout << "Check_peak: Left wing level 3 failed: " << i << endl;
			return false;
		}
	}
	// right wing check_level_3: if any of the right wing points is greater than the current point: return false immediately
	for (int i = upper_index_2 + 1; (i <= upper_index_3 && i < pt.capacity); i++){
		if (smooth_record[i] > current_cons){
			//cout << "Check_peak: Right wing level 3 failed: " << i << endl;
			return false;
		}
	}
	return true;
}

bool check_trough(peak_trough& pt, int cell_index, int current_index){
	double current_cons = pt.smooth_cons[cell_index][current_index];
	double* smooth_record = pt.smooth_cons[cell_index];
	vector<int> * peaks = pt.peaks[cell_index];
	vector<int> * troughs = pt.troughs[cell_index];
	int upper_index_1 = current_index + WING_CHECK_SIZE_LEVEL_1;
	int lower_index_1 = current_index - WING_CHECK_SIZE_LEVEL_1;
	int upper_index_2 = upper_index_1 + WING_CHECK_SIZE_LEVEL_2;
	int lower_index_2 = lower_index_1 - WING_CHECK_SIZE_LEVEL_2;
	int upper_index_3 = upper_index_2 + WING_CHECK_SIZE_LEVEL_3;
	int lower_index_3 = lower_index_2 - WING_CHECK_SIZE_LEVEL_3;
	
	if (peaks->size() > 0){
		if (troughs->size() > 0 && troughs->back() > peaks->back()){
			//cout << "False bc 2 troughs in a row" << endl; 
			// No two troughs without any peak in between
			return false;
		}
		if (smooth_record[peaks->back()] <= current_cons){
			//cout << "False be previous peak is lower than this trough" << endl;
			// previous peak cannot be lower that this trough
			return false;
		}
	}else { // no peaks yet
		if (troughs->size() > 0){
			//cout << "False bc no peaks yet but this is the second trough found" << endl;
			// No two troughs without any peak in between
			return false;
		}
	}
	
	// left wing check_level_1: if any of the left wing points is smaller than the current point: false
	for (int i = current_index - 1; (i >= lower_index_1 && i >= 0); i --){
		if (smooth_record[i] < current_cons){
			//cout << "Check trough: Left wing level 1 failed: " << i << endl;
			return false;
		}
	}	
	// right wing check_level_1: if any of the right wing points is greater than the current point: false
	for (int i = current_index + 1; (i <= upper_index_1 && i < pt.capacity); i ++){
		if (smooth_record[i] < current_cons){
			//cout << "Check trough: Right wing level 1 falied: " << i << endl;
			return false;
		}
	}
	// left_wing_check_level_2: if any of the left wing points is greater than the current point: return false immediately
	for (int i = lower_index_1 - 1; (i >= lower_index_2 && i >= 0); i --){
		if (smooth_record[i] < current_cons){
			//cout << "Check trough: Left wing level 2 failed: " << i << endl;
			return false;
		}
	}
	// right wing check_level_2: if any of the right wing points is greater than the current point: return false immediately
	for (int i = upper_index_1 + 1; (i <= upper_index_2 && i < pt.capacity); i++){
		if (smooth_record[i] < current_cons){
			//cout << "check trough: right wing level 2 failed: " << i << endl;
			return false;
		}
	}
	// left_wing_check_level_3: if any of the left wing points is greater than the current point: return false immediately
	for (int i = lower_index_2 - 1; (i >= lower_index_3 && i >= 0); i--){
		if (smooth_record[i] < current_cons){
			//cout << "Check trough: Left wing level 3 failed: " << i << endl;
			return false;
		}
	}
	// right wing check_level_3: if any of the right wing points is greater than the current point: return false immediately
	for (int i = upper_index_2 + 1; (i <= upper_index_3 && i < pt.capacity); i++){
		if (smooth_record[i] < current_cons){
			//cout << "Check trough: Right wing level 3 failed: " << i << endl;
			return false;
		}
	}
	return true;
}

double calculate_middle_peak_or_trough(double* smooth_record, vector<int>* pt_indices){
	int mid_index = pt_indices->size() / 2;
	double avg_mid_cons = smooth_record[mid_index];
	double num_mid_pt_real = 1.0;
	int num_mid_pt_side = NUM_MID_PT / 2;
	for (int i = 0; i < num_mid_pt_side; i++){
		avg_mid_cons += smooth_record[pt_indices->at(mid_index - i - 1)];
		avg_mid_cons += smooth_record[pt_indices->at(mid_index + i + 1)];
		num_mid_pt_real += 2.0;
	}
	avg_mid_cons = avg_mid_cons / num_mid_pt_real;
	if (avg_mid_cons == 0){
		return 1; // so that when we calculate peak to trough ratio, we don't get error
	}
	return avg_mid_cons;
}

double calculate_last_peak_or_trough(double* smooth_record, vector<int>* pt_indices){
	double avg_last_cons = 0;
	int num_peak_or_trough = pt_indices->size();
	for (int i = 0; i < NUM_LAST_PT; i++){
		avg_last_cons += smooth_record[pt_indices->at(num_peak_or_trough - 1 - i)];
	}
	avg_last_cons = avg_last_cons / NUM_LAST_PT;
	if (avg_last_cons == 0){
		return 1; // so that when we calculate peak to trough ratio, we dont get error
	}
	return avg_last_cons;
}

void calculate_p2t(peak_trough& pt, features& fts, int con_index, int cell_index){
	// Find data about this cell
	double* smooth_record = pt.smooth_cons[cell_index];
	vector<int> * peaks = pt.peaks[cell_index];
	vector<int> * troughs = pt.troughs[cell_index];
	
	// Process the middle ptt	
	// find the average concentrations of middle peaks and middle troughs
	double avg_mid_peak = calculate_middle_peak_or_trough(smooth_record, peaks);
	double avg_mid_trough = calculate_middle_peak_or_trough(smooth_record, troughs);
	
	// report middle peak and trough
	fts.mid_ptt[con_index][cell_index] = (avg_mid_peak / avg_mid_trough);
	
	// Process the last ptt
	double avg_last_peak = calculate_last_peak_or_trough(smooth_record, peaks);
	double avg_last_trough = calculate_last_peak_or_trough(smooth_record, troughs);
	
	// report last peak and trough
	fts.last_ptt[con_index][cell_index] = (avg_last_peak / avg_last_trough);
}

void cell_to_cell_correlations(peak_trough& pt, features& fts, int cons_index){
	for (int i = 0; i < (pt.num_cells - 1); i++){
		fts.correlation[cons_index][i] = calculate_correlation(pt.smooth_cons[0], pt.smooth_cons[i + 1], pt.capacity);
	}
}

double calculate_correlation(double* X, double* Y, int length){
	double sum_XY = 0;
	double sum_X = 0;
	double sum_Y = 0;
	double sum_X_squared = 0;
	double sum_Y_squared = 0;
	for (int i = 0; i < length; i++){
		sum_XY += (X[i] * Y[i]);
		sum_X += X[i];
		sum_Y += Y[i];
		sum_X_squared += (pow(X[i], 2.0));
		sum_Y_squared += (pow(Y[i], 2.0));
	}
	double nominator = sum_XY - ((sum_X * sum_Y) / (double)length); 
	double X_var = sum_X_squared - (pow(sum_X, 2.0) / (double)length);
	double Y_var = sum_Y_squared - (pow(sum_Y, 2.0) / (double)length);
	double denominator = pow(X_var * Y_var, 0.5);
	double r = nominator / denominator;
	return r;
}

bool check_correlations(features& fts){
	double avg_r = 0;
	for (int i = 0; i < NUM_KEEP_STATES; i++){
		for(int j = 0; j < (fts.num_cells - 1); j++){
			avg_r += fts.correlation[i][j];
		}
	}
	avg_r = avg_r / (NUM_KEEP_STATES * (fts.num_cells - 1));
	if (avg_r >= LOW_BOUND_WT_CORR){
		return true;
	}
	return false;
}

void process_noise_data(embryo& em, features& fts, input_params& ip){
	int num_slices = (int)((em.absolute_time - 2 * ip.record_granularity) / MINUTE_PER_SLICE);
	slices sl (em.num_cells, num_slices);
	process_slices (em, sl, ip);
	calculate_slice_noise(em, sl);
	binned_data bd (fts.num_bin);
	double * bounds = new double [bd.num_bin];
	calculate_bounds(sl, bd, bounds);
	put_data_to_bin(sl, bd, bounds); 
	calculate_bins (bd, fts);
}

/*
 * Given the simulationd data, calculate the concentration in each slice of each cell
 * Each slice is defined as the average concentrations of 3 consecutive records, each ip.record_granularity apart 
 * Each slice is MINUTE_PER_SLICE mins apart
 * For example, if simulation time is 600 mins and MINUTE_PER_SLICE is 5, the number of slices would be 120.
 * The first slice is the average of time 4.9, 5 and 5.1. So on...
 * This function fills in slice.her1 and slice.her7 tables
 */
void process_slices (embryo& em, slices& sl, input_params& ip){
	int steps_per_slice = (int)(MINUTE_PER_SLICE / ip.record_granularity);
	int before_offset = steps_per_slice - 1;
	int after_offset = steps_per_slice + 1;
	for (int i = 0; i < sl.num_slices; i++){
		int before_index = i * steps_per_slice + before_offset;
		int current_index = (i + 1) * steps_per_slice;
		int after_index = i * steps_per_slice + after_offset;
		for (int j = 0; j < em.num_cells; j++){
			cell* cc = em.cell_list[j];
			sl.her1[i][j] = (((cc->cons_record)[KEEPMH1])->at(before_index) 
								+ ((cc->cons_record)[KEEPMH1])->at(current_index)
								+ ((cc->cons_record)[KEEPMH1])->at(after_index)) / 3.0;
			sl.her7[i][j] = (((cc->cons_record)[KEEPMH7])->at(before_index) 
								+ ((cc->cons_record)[KEEPMH7])->at(current_index)
								+ ((cc->cons_record)[KEEPMH7])->at(after_index)) / 3.0;
		}
	}
}
	

/*
 * This function calculates the intrinsic and extrinsic noise of each slice throughout simulation.
 * It fills in slice.avg_h1, slice.avg_h7, avg_her, in_noise, ex_noise, min_her, max_her
 */
void calculate_slice_noise(embryo& em, slices& sl){
	double min_her = INFINITY;
	double max_her = -1;
	
	// Calculate the average concentrations in each slice
	for (int i = 0; i < sl.num_slices; i++){
		double h1 = 0; 
		double h7 = 0; 
		for (int j = 0; j < sl.num_cells; j ++){
			h1 += sl.her1[i][j];
			h7 += sl.her7[i][j];
		}
		sl.avg_h1[i] = h1 / sl.num_cells;
		sl.avg_h7[i] = h7 / sl.num_cells;
		sl.avg_her[i] = (h1 + h7) / sl.num_cells;
		if (sl.avg_her[i] > max_her){
			max_her = sl.avg_her[i];
		}
		if (sl.avg_her[i] < min_her){
			min_her = sl.avg_her[i];
		}
	}
	
	// Calculate intrinsic and extrinsic noise in each slice
	for (int i = 0; i < sl.num_slices; i ++){
		double in_slice_noise = 0;
		double ex_slice_noise = 0;
		for (int j = 0; j < sl.num_cells; j++){
			double cell_in_term = (sl.her1[i][j] / sl.avg_h1[i]) - (sl.her7[i][j] / sl.avg_h7[i]);
			in_slice_noise += pow(cell_in_term,2);
			double cell_ex_term = sl.her1[i][j] * sl.her7[i][j];
			ex_slice_noise += cell_ex_term;
		}
		in_slice_noise = (in_slice_noise / sl.num_cells) / 2.0;
		sl.in_noise[i] = in_slice_noise;
		ex_slice_noise = ex_slice_noise / sl.num_cells;
		ex_slice_noise = ex_slice_noise - (sl.avg_h1[i] * sl.avg_h7[i]);
		ex_slice_noise = ex_slice_noise / (sl.avg_h1[i] * sl.avg_h7[i]);
		sl.ex_noise[i] = ex_slice_noise;
	}
	
	sl.min_her = min_her;
	sl.max_her = max_her;
}

void calculate_bounds(slices& sl, binned_data& bd, double* bounds){
	double interval = (sl.max_her - sl.min_her) / bd.num_bin;
	for (int i = 0; i < bd.num_bin; i++){
		bounds[i] = sl.min_her + interval * (i + 1);
	}
}

/*
 * put slices' data into the appropriate bin
 */
void put_data_to_bin(slices& sl, binned_data& bd, double* bounds){
	for (int i = 0; i < sl.num_slices; i++){
		//cout << sl.avg_her[i] << endl;
		int j = 0; 
		while (sl.avg_her[i] >= bounds[j] && j < bd.num_bin){
			j++;
		}
		if (sl.avg_h1[i] > 0 && sl.avg_h7[i] > 0 && j < bd.num_bin){
			(bd.her[j])->push_back(sl.avg_her[i]);
			(bd.in_noise[j])->push_back(sl.in_noise[i]);
			(bd.ex_noise[j])->push_back(sl.ex_noise[i]);
		}
	}
}

/*
 * Calculate the average concetrations, in and ex noise levels of each bin  
 * and fill in fts.avg_cons, intrinsic and extrinsic tables
 */
void calculate_bins(binned_data& bd, features& fts){
	double her = 0;
	double in = 0;
	double ex = 0;
	for (int i = 0; i < bd.num_bin; i++){
		her = 0;
		in = 0; 
		ex = 0;
		for (unsigned int j = 0; j < ((bd.her[i])->size()); j++){
			her += (bd.her[i])->at(j);
			in += (bd.in_noise[i])->at(j);
			ex += (bd.ex_noise[i])->at(j);
		}
		if (((bd.her[i])->size()) > 0){
			fts.avg_cons[i] = her / (bd.her[i])->size();
			fts.intrinsic[i] = in / (bd.her[i])->size();
			fts.extrinsic[i] = ex / (bd.her[i])->size();
		}
	}
}
