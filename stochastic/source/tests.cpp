#include <math.h>
#include "tests.hpp"

int test_wildtype(input_params& ip, embryo& em, features& wtf){
	// check mRNA boundaries
	// process peaks and troughs
	process_smooth_data(em, wtf);
	// sustained oscillation
	// H1, H7 spatial amplitude
	// binned data used for all types of noises calculation
	binned_data bd(ip.num_bin);
	// calculate noises
	process_binned_data(bd, em, wtf);
	// intrisic noise test
	// extrinsic noise test
	
	return 0;
}

void process_smooth_data(embryo& em, features& fts){
	// find the good capcaity for peaks_troughs
	int capacity = find_peaks_troughs_capacity(em);
	peak_trough pt (capacity);
	for (int i = 0; i < em.num_cells; i++){
		for (int j = 0; j < NUM_KEEP_STATES; j ++){
			smooth_data(em, pt, i, j);
			find_peaks_and_troughs(pt,fts, i, j);
			pt.reset();
		}
	}
}

int find_peaks_troughs_capacity(embryo& em){
	int capacity = 0; 
	for (int i = 0; i < em.num_cells; i++){
		for (int j = 0; j < NUM_KEEP_STATES; j++){
			if (capacity < (em.cell_list[i]->cons_record)[j]->size()){
				capacity = (em.cell_list[i]->cons_record)[j]->size() - WINDOW_SIZE * 2;
			}
		}
	}
	return capacity;
}

void smooth_data(embryo& em, peak_trough& pt, int cell_index, int con_index){
	vector<int> * current_vector = (em.cell_list[cell_index]->cons_record)[con_index];
	pt.num_cons = current_vector->size() - WINDOW_SIZE * 2;
	int sum = 0; 
	int avg_divisor = WINDOW_SIZE * 2;
	for (int i = 0; i < avg_divisor; i ++){
		sum += current_vector->at(i);
	}
	pt.smooth_cons[0] = double(sum) / double(avg_divisor);
	int remove_index = 0;
	int add_index = WINDOW_SIZE * 2 + 1;
	
	for (int i = 1; i < pt.num_cons; i ++){
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
	// calculate the amplitude of the concentrations of this cell
	avg_peak = avg_peak / (pt.peaks->size());
	avg_trough = avg_trough / (pt.troughs->size());
	fts.avg_amplitude[con_index][cell_index] = avg_peak - avg_trough;

	//calculate peak_to_trough
	calculate_p2t(pt, fts, con_index, cell_index);
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
	for (int i = current_index + 1; (i >= upper_index && i < pt.num_cons); i++){
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

void process_binned_data(binned_data& bd, embryo& em, features& fts){
	double * bounds = new double [bd.num_bin];
	find_bin_bounds(bounds, em, bd.num_bin);
	put_data_into_bins(bounds, em, bd);
	calculate_noise_features(bd, fts);
	delete [] bounds;
}

void find_bin_bounds(double* bounds, embryo& em, int num_bin){
	int max = -1 ;
	int min = INFINITY;
	int total_cons = 0;
	vector<int>* her1;
	vector<int>* her7;
	for (int i = 0; i < em.num_cells; i++){
		her1 = ((em.cell_list[i])->cons_record)[KEEPMH1];
		her7 = ((em.cell_list[i])->cons_record)[KEEPMH7];
		for (int j = 0; j < her1->size(); j++){
			total_cons = her1->at(j) + her7->at(j);
			if (total_cons < min){
				min = total_cons;
			}
			if (total_cons > max){
				max = total_cons;
			}
		}
	}
	
	double bin_size = ((double) max - (double) min)	/ (double)num_bin;
	for (int i = 0; i < num_bin; i++){
		bounds[i] = (double) min + bin_size * (i + 1);
	}
}

void put_data_into_bins(double* bounds, embryo& em, binned_data& bd){
	vector<int>* her1;
	vector<int>* her7;
	int current_bound_index = 0; 
	double total_cons;
	for (int i = 0; i < em.num_cells; i++){
		her1 = ((em.cell_list[i])->cons_record)[KEEPMH1];
		her7 = ((em.cell_list[i])->cons_record)[KEEPMH7];
		for (int j = 0; j < her1->size(); j++){
			total_cons = her1->at(j) + her7->at(j);
			while (total_cons > bounds[current_bound_index]){
				current_bound_index += 1;
			}
			(bd.her1_data[current_bound_index])->push_back(her1->at(j));
			(bd.her7_data[current_bound_index])->push_back(her7->at(j));
			bd.total_her1[current_bound_index] += her1->at(j);
			bd.total_her7[current_bound_index] += her7->at(j);
			current_bound_index = 0;
		}
	}
}

void calculate_noise_features(binned_data& bd, features& fts){
	// calculate average concentrations of different bins
	for (int i = 0; i < bd.num_bin; i++){
		bd.avg_her1[i] = (double)bd.total_her1[i] / (double)(bd.her1_data[i])->size();
		bd.avg_her7[i] = (double)bd.total_her7[i] / (double)(bd.her7_data[i])->size();
	}
	double avg1;
	double avg7;
	vector<int> * her1;
	vector<int> * her7;
	double intrinsic_total = 0;
	double extrinsic_total = 0;
	for (int i = 0; i < bd.num_bin; i++){
		avg1 = bd.avg_her1[i];
		avg7 = bd.avg_her7[i];
		her1 = bd.her1_data[i];
		her7 = bd.her7_data[i];
		for (int j = 0; j < her1->size(); j++){
			intrinsic_total += pow((her1->at(j) / avg1) - (her7->at(j) /avg7), 2.0);
			extrinsic_total += (double) her1->at(j) * (double) her7->at(j);
		}
		// intrinsic noise of bin
		fts.intrinsic[i] = 0.5 * intrinsic_total / ((double)her1->size());
		// extrinsic noise of bin
		fts.extrinsic[i] = (extrinsic_total / (her1->size())) - avg1 * avg7;
		fts.extrinsic[i] /= (avg1 * avg7);
		// average mRNAs levels of bin
		fts.avg_cons[i] = avg1 + avg7;
	}
}
