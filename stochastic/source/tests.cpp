#include <math.h>
#include "tests.hpp"
/*
 * Max_score : 5
 */
int test_wildtype(input_params& ip, embryo& em, features& wtf){
	int score = 0;
	// check mRNA boundaries
	bool mRNA_in_bound = true; //check_mRNA_boundaries(em);
	if (! mRNA_in_bound){
		return 0;
	}
	// process peaks and troughs
	process_smooth_data(em, wtf);
	// sustained oscillation
	score += check_WT_sustained(wtf);
	// H1, H7 spatial amplitude
	score += check_WT_amplitude(wtf);
	// calculate noises
	process_noise_data(em, wtf);
	// average mRNA concentrations tests
	score += check_WT_avg_cons(wtf);
	// intrisic noise test
	score += check_WT_intrinsic (wtf);
	// extrinsic noise test
	score += check_WT_extrinsic (wtf);
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

int check_WT_sustained (features& wtf){
	for (int i = 0; i < wtf.num_cells; i ++){
		if (wtf.mid_ptt[KEEPMH1][i] < LOWER_BOUND_PTT){
			return 0;
		}
		if (wtf.last_ptt[KEEPMH7][i] < LOWER_BOUND_PTT){
			return 0;
		}
		if ((wtf.mid_ptt[KEEPMH1][i] / wtf.last_ptt[KEEPMH1][i]) > UPPER_BOUND_MTL){
			return 0;
		}
	}	
	return 1;
}


int check_WT_amplitude (features& wtf){
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
	
	if (h1_amplitude < WT_H1_AMP_LOW){
		return 0;
	}
	if (h1_amplitude > WT_H1_AMP_HIGH){
		return 0;
	}
	if (h7_amplitude < WT_H7_AMP_LOW){
		return 0;
	}
	if (h7_amplitude > WT_H7_AMP_LOW){
		return 0;
	}
	return 1;
}

int check_WT_avg_cons(features& wtf){
	// Bin 1
	if (wtf.avg_cons[0] < LOW_WT_HER_BIN1){
		return 0;
	}
	if (wtf.avg_cons[0] > UP_WT_HER_BIN1){
		return 0;
	}
	// Bin 2
	if (wtf.avg_cons[1] < LOW_WT_HER_BIN2){
		return 0;
	}
	if (wtf.avg_cons[1] > UP_WT_HER_BIN2){
		return 0;
	}
	// Bin 3
	if (wtf.avg_cons[2] < LOW_WT_HER_BIN3){
		return 0;
	}
	if (wtf.avg_cons[2] > UP_WT_HER_BIN3){
		return 0;
	}
	return 1;
}

int check_WT_intrinsic(features& wtf){
	// Bin 1
	if (wtf.intrinsic[0] < LOW_WT_IN_NOISE_BIN1){
		return 0;
	}
	if (wtf.intrinsic[0] > UP_WT_IN_NOISE_BIN1){
		return 0;
	}
	// Bin 2
	if (wtf.intrinsic[1] < LOW_WT_IN_NOISE_BIN2){
		return 0;
	}
	if (wtf.intrinsic[1] > UP_WT_IN_NOISE_BIN2){
		return 0;
	}
	// Bin 3
	if (wtf.extrinsic[2] < LOW_WT_IN_NOISE_BIN3){
		return 0;
	}
	if (wtf.extrinsic[2] > UP_WT_IN_NOISE_BIN3){
		return 0;
	}
	return 1;
}

int check_WT_extrinsic (features& wtf){
	// Bin 1
	if (wtf.extrinsic[0] < LOW_WT_EX_NOISE_BIN1){
		return 0;
	}
	if (wtf.extrinsic[0] > UP_WT_EX_NOISE_BIN1){
		return 0;
	}
	// Bin 2
	if (wtf.extrinsic[1] < LOW_WT_EX_NOISE_BIN2){
		return 0;
	}
	if (wtf.extrinsic[1] > UP_WT_EX_NOISE_BIN2){
		return 0;
	}
	// Bin 3
	if (wtf.extrinsic[2] < LOW_WT_EX_NOISE_BIN3){
		return 0;
	}
	if (wtf.extrinsic[2] > UP_WT_EX_NOISE_BIN3){
		return 0;
	}
	return 1;
}

void process_smooth_data(embryo& em, features& fts){
	// find the good capcaity for peaks_troughs
	int capacity = ((em.cell_list[0])->cons_record[0])->size();
	//int capacity = em.time_record->size() - WINDOW_SIZE * 2;
	peak_trough pt (capacity);
	for (int i = 0; i < em.num_cells; i++){
		for (int j = 0; j < NUM_KEEP_STATES; j ++){
			smooth_data(em, pt, i, j);
			find_peaks_and_troughs(pt,fts, i, j);
			pt.reset();
		}
	}
}


void smooth_data(embryo& em, peak_trough& pt, int cell_index, int con_index){
	vector<int> * current_vector = (em.cell_list[cell_index]->cons_record)[con_index];
	pt.num_cons = current_vector->size() - WINDOW_SIZE * 2;
	int sum = 0; 
	int avg_divisor = WINDOW_SIZE * 2 + 1;
	for (int i = 0; i < avg_divisor; i ++){
		sum += current_vector->at(i);
	}
	pt.smooth_cons[0] = double(sum) / double(avg_divisor);
	int remove_index = 0;
	int add_index = WINDOW_SIZE * 2 + 1;
	
	for (int i = 1; i < (pt.num_cons); i ++){
		sum -= current_vector->at(remove_index);
		sum += current_vector->at(add_index);
		pt.smooth_cons[i] = double(sum) / double(avg_divisor);
		remove_index += 1;
		add_index += 1;
	}
	
}

void find_peaks_and_troughs(peak_trough& pt, features& fts, int cell_index, int con_index){
	// amplitudes and ratios
	int max_cons = pt.num_cons;
	int i = 1; 
	bool is_peak = false;
	bool is_trough = false;
	double avg_peak = 0;
	double avg_trough = 0;
	while (i < max_cons){
		// potential peak
		if (pt.smooth_cons[i] >= pt.smooth_cons[i - 1] && pt.smooth_cons[i] >= pt.smooth_cons[i + 1]){
			is_peak = check_peak(pt, i);
		}
		// potential trough
		else if (pt.smooth_cons[i] <= pt.smooth_cons[i - 1] && pt.smooth_cons[i] <= pt.smooth_cons[i + 1]){
			is_trough = check_trough(pt, i);
		}
		if (is_peak){
			pt.peaks->push_back(i);
			avg_peak += pt.smooth_cons[i];
			i += WING_CHECK_SIZE;
			is_peak = false;
		}
		if (is_trough){
			pt.troughs->push_back(i);
			avg_trough += pt.smooth_cons[i];
			i += WING_CHECK_SIZE;
			is_trough = false;
		}
		i += 1;
	}
	
	// if we have too few peaks and troughs, set to 0
	// so that later on, I can choose what to take the average of when calculating
	// average of all cells 
	if (pt.peaks->size() < LOW_NUM_PT || pt.troughs->size() < LOW_NUM_PT){
		fts.avg_amplitude[con_index][cell_index] = 0;
		fts.mid_ptt[con_index][cell_index] = 0;
		fts.last_ptt[con_index][cell_index] = 0;
	}
	else {
		// calculate the amplitude of the concentrations of this cell
		avg_peak = avg_peak / (pt.peaks->size());
		avg_trough = avg_trough / (pt.troughs->size());
		fts.avg_amplitude[con_index][cell_index] = avg_peak - avg_trough;
		//calculate peak_to_trough
		calculate_p2t(pt, fts, con_index, cell_index);
	}
}

bool check_peak(peak_trough& pt, int current_index){
	double current_cons = pt.smooth_cons[current_index];
	int upper_index = current_index + WING_CHECK_SIZE;
	int lower_index = current_index - WING_CHECK_SIZE;
	if (pt.troughs->size() > 0){
		if (pt.peaks->size() > 0 && pt.peaks->back() > pt.troughs->back()){
			//No two peaks without any trough in between
			return false;
		}
		if (pt.smooth_cons[pt.troughs->back()] >= current_cons){
			// previous trough cannot be higher than this peak
			return false;
		}
	}else{ // no troughs yet
		if (pt.peaks->size() > 0){
			// No two peaks without any trough in between
			return false;
		}
	}
	// left wing check : if any of the left wing points is greater than the current point: return false immediately
	for (int i = current_index - 1; (i >= lower_index && i >= 0); i --){
		if (pt.smooth_cons[i] > current_cons){
			return false;
		}
	}
	// right wing check; if any of the right wing points is greater than the current point: return false immediately
	for (int i = current_index + 1; (i <= upper_index && i < pt.num_cons); i++){
		if (pt.smooth_cons[i] > current_cons){
			return false;
		}
	}
	return true;
}

bool check_trough(peak_trough& pt, int current_index){
	double current_cons = pt.smooth_cons[current_index];
	int upper_index = current_index + WING_CHECK_SIZE;
	int lower_index = current_index - WING_CHECK_SIZE;
	if (pt.peaks->size() > 0){
		if (pt.troughs->size() > 0 && pt.troughs->back() > pt.peaks->back()){
			// No two troughs without any peak in between
			return false;
		}
		if (pt.smooth_cons[pt.peaks->back()] <= current_cons){
			// previous peak cannot be lower that this trough
			return false;
		}
	}else { // no peaks yet
		if (pt.troughs->size() > 0){
			// No two troughs without any peak in between
			return false;
		}
	}
	
	// left wing check: if any of the left wing points is smaller than the current point: false
	for (int i = current_index - 1; (i >= lower_index && i >= 0); i --){
		if (pt.smooth_cons[i] < current_cons){
			return false;
		}
	}
	
	// right wing check: if any of the right wing points is greater than the current point: false
	for (int i = current_index + 1; (i <= upper_index && i < pt.num_cons); i ++){
		if (pt.smooth_cons[i] < current_cons){
			return false;
		}
	}
	return true;
}

void calculate_p2t(peak_trough& pt, features& fts, int con_index, int cell_index){
	int mid_peak = pt.peaks->size() / 2;
	int mid_trough = pt.troughs->size() / 2;
	int mp_index1 = pt.peaks->at(mid_peak);
	int mp_index2 = pt.peaks->at(mid_peak - 1);
	int mt_index1 = pt.troughs->at(mid_trough);
	int mt_index2 = pt.troughs->at(mid_trough - 1);
	if (mid_peak == mid_trough && pt.peaks->size() > pt.troughs->size()){
		double peak_con = pt.smooth_cons[mp_index1];
		double trough1 = pt.smooth_cons[mt_index1];
		double trough2 = pt.smooth_cons[mt_index2];
		fts.mid_ptt[con_index][cell_index] = peak_con / ((trough1 + trough2) / 2);
	}
	else if (mid_peak == mid_trough && pt.peaks->size() < pt.troughs->size()){
		double trough = pt.smooth_cons[mt_index1];
		double peak1 = pt.smooth_cons[mp_index1];
		double peak2 = pt.smooth_cons[mp_index2];
		fts.mid_ptt[con_index][cell_index] = ((peak1 + peak2) / 2) / trough;
	}
	else if (mid_peak == mid_trough && pt.peaks->size() == pt.troughs->size()){
		if (pt.peaks->size() % 2 == 0){
			double peak1 = pt.smooth_cons[mp_index1];
			double peak2 = pt.smooth_cons[mp_index2];
			double trough1 = pt.smooth_cons[mt_index1];
			double trough2 = pt.smooth_cons[mt_index2];
			fts.mid_ptt[con_index][cell_index] = ((peak1 + peak2) / 2) / ((trough1 + trough2) / 2);
		}
		else{
			double peak = pt.smooth_cons[mp_index1];
			double trough = pt.smooth_cons[mt_index1];
			fts.mid_ptt[con_index][cell_index] = peak / trough;
		}
	}
	else if (mid_peak > mid_trough){
		double trough = pt.smooth_cons[mt_index1];
		double peak1 = pt.smooth_cons[mp_index1];
		double peak2 = pt.smooth_cons[mp_index2];
		fts.mid_ptt[con_index][cell_index] = ((peak1 + peak2) / 2) / trough;
	}
	else if (mid_peak < mid_trough){
		double trough1 = pt.smooth_cons[mt_index1];
		double trough2 = pt.smooth_cons[mt_index2];
		double peak = pt.smooth_cons[mp_index1];
		fts.mid_ptt[con_index][cell_index] = peak / ((trough1 + trough2) / 2);
	}
	
	int lp_index1 = pt.peaks->back();
	int lp_index2 = pt.peaks->at(pt.peaks->size() - 2);
	int lt_index1 = pt.troughs->back();
	int lt_index2 = pt.troughs->at(pt.troughs->size() - 2);
	if (lp_index1 < lt_index1){ // last peak after last trough
		double peak = pt.smooth_cons[lp_index1];
		double trough1 = pt.smooth_cons[lt_index1];
		double trough2 = pt.smooth_cons[lt_index2];
		fts.last_ptt[con_index][cell_index] = peak / ((trough1 + trough2) / 2);
	}
	else {
		double peak1 = pt.smooth_cons[lp_index1];
		double peak2 = pt.smooth_cons[lp_index2];
		double trough = pt.smooth_cons[lt_index1];
		fts.last_ptt[con_index][cell_index] = ((peak1 + peak2) / 2) / trough;
	}
}

void process_noise_data(embryo& em, features& fts){
	slices sl (em.num_cells);
	process_slices (em, sl);
	calculate_slice_noise(em, sl);
	binned_data bd (fts.num_bin);
	double * bounds = new double [bd.num_bin];
	calculate_bounds(sl, bd, bounds);
	put_data_to_bin(sl, bd, bounds); 
	calculate_bins (bd, fts);
}

void process_slices (embryo& em, slices& sl){
	double low_time;
	double high_time;
	int current_index = 0;
	double current_time = em.time_record->at(current_index);
	double* total_her1 = new double [em.num_cells];
	double* total_her7 = new double [em.num_cells];
	int total_cells;
	for (int i = 0; i < NUM_SLICES; i++){
		low_time = MINUTE_PER_SLICE * (i + 1) - SLICE_EXTENSION;
		high_time = MINUTE_PER_SLICE * (i + 1) + SLICE_EXTENSION;
		memset(total_her1, 0, sizeof(double) * em.num_cells);
		memset(total_her7, 0, sizeof(double) * em.num_cells);
		total_cells = 0;
		while (current_time < high_time){
			if (current_time > low_time){
				for (int j = 0; j < em.num_cells; j++){
					total_her1[j] += ((em.cell_list[j])->cons_record[KEEPMH1])->at(current_index);
					total_her7[j] += ((em.cell_list[j])->cons_record[KEEPMH7])->at(current_index);
				}
				total_cells += 1;
			}
			current_index += 1;
			current_time = em.time_record->at(current_index);
		}
		
		if (total_cells > 0){
			for (int j = 0; j < em.num_cells; j++){
				sl.her1[i][j] = total_her1[j] / total_cells;
				sl.her7[i][j] = total_her7[j] / total_cells;
			}
		}
	}
	

	delete [] total_her1;
	delete [] total_her7;
}

void calculate_slice_noise(embryo& em, slices& sl){
	double min_her = INFINITY;
	double max_her = -1;
	// calculate the average concentrations in each slice
	for (int i = 0; i < NUM_SLICES; i++) {
		if (sl.her1[i][0] > -1){
			sl.avg_h1[i] = 0;
			sl.avg_h7[i] = 0;
			for (int j = 0; j < em.num_cells; j ++){
				sl.avg_h1[i] += sl.her1[i][j];
				sl.avg_h7[i] += sl.her7[i][j];
			}
			sl.avg_h1[i] = sl.avg_h1[i] / em.num_cells;
			sl.avg_h7[i] = sl.avg_h7[i] / em.num_cells;
			sl.avg_her[i] = sl.avg_h1[i] + sl.avg_h7[i];
		}
	}
	
	
	// calculate intrinsic and extrinsic noise for each slice
	for (int i = 0; i < NUM_SLICES; i ++){
		double intrinsic = 0;
		double extrinsic = 0;
		if (sl.avg_h1[i] > 0 && sl.avg_h7[i] > 0){
			double ah1 = sl.avg_h1[i];
			double ah7 = sl.avg_h7[i];
			if (sl.avg_her[i] < min_her){
				min_her = sl.avg_her[i];
			}
			if (sl.avg_her[i] > max_her){
				max_her = sl.avg_her[i];
			}
			for (int j = 0; j < em.num_cells; j ++){
				intrinsic += pow(((sl.her1[i][j] / ah1) - (sl.her7[i][j] / ah7)), 2.0);
				extrinsic += sl.her1[i][j] * sl.her7[i][j]; 
			}
			
			// intrinsic noise of slice
			intrinsic = intrinsic / em.num_cells;
			intrinsic = intrinsic * 0.5;
			sl.in_noise[i] = intrinsic;
			// extrinsic noise of slice
			extrinsic = extrinsic / em.num_cells;
			extrinsic -= ah1 * ah7;
			extrinsic /= (ah1 * ah7);
			sl.ex_noise[i] = extrinsic;
		}
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

void put_data_to_bin(slices& sl, binned_data& bd, double* bounds){
	for (int i = 0; i < NUM_SLICES; i++){
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

void calculate_bins(binned_data& bd, features& fts){
	double her = 0;
	double in = 0;
	double ex = 0;
	for (int i = 0; i < bd.num_bin; i++){
		her = 0;
		in = 0; 
		ex = 0;
		for (int j = 0; j < (bd.her[i])->size(); j++){
			her += (bd.her[i])->at(j);
			in += (bd.in_noise[i])->at(j);
			ex += (bd.ex_noise[i])->at(j);
		}
		if ((bd.her[i])->size() > 0){
			fts.avg_cons[i] = her / (bd.her[i])->size();
			fts.intrinsic[i] = in / (bd.her[i])->size();
			fts.extrinsic[i] = ex / (bd.her[i])->size();
		}
	}
}
