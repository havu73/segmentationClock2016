#include "tests.hpp"
void find_peaks_troughs (con_levels& cons, peak_trough& pt, int index){
	for (int i = 1; i < (cons.num_sim_steps - 1); i ++){
		//Peak
		if ((cons.sim_data[index][i] > cons.sim_data[index][i + 1]) 
				&& (cons.sim_data[index][i] > cons.sim_data[index][i - 1])){
			pt.peaks.push_back(i);
		}
		//Trough
		else if ((cons.sim_data[index][i] < cons.sim_data[index][i + 1]) 
				&& (cons.sim_data[index][i] < cons.sim_data[index][i - 1])){
			pt.troughs.push_back(i);
		}
	}
}

int test_sustained_oscillation(peak_trough& pt, con_levels& cons, int con_index){
	double first_p2t = 0;
	double mid_p2t = 0;
	double last_p2t = 0;
	
	//If we could not find enough peaks and troughs to make the simulation reasonable --> score = 0
	if (pt.peaks.size() < LOW_BOUND_NUM_AMP) {
		return 0;
	}
	
	calculate_p2t(pt, cons, con_index, &first_p2t, &mid_p2t, &last_p2t);
	bool passed = true;
	if (mid_p2t < LOW_BOUND_P2T){
		passed = false;
	}
	if (last_p2t < LOW_BOUND_P2T){
		passed = false;
	}
	if ((mid_p2t / last_p2t) > UP_BOUND_M2L){
		passed = false;
	}
	if (passed){
		return 1;
	}
	return 0;
}

void calculate_p2t(peak_trough& pt, con_levels& cons, int con_index, double* first_p2t, double* mid_p2t, double* last_p2t){
	if (pt.peaks[0] < pt.troughs[0]){
		double peak_con_1 = cons.sim_data[con_index][pt.peaks[0]];
		double peak_con_2 = cons.sim_data[con_index][pt.peaks[1]];
		double trough_con = cons.sim_data[con_index][pt.troughs[0]];
		*first_p2t = ((peak_con_1 + peak_con_2)/2) / trough_con;
	}
	else if (pt.peaks[0] > pt.troughs[0]){
		double tr_con_1 = cons.sim_data[con_index][pt.troughs[0]];
		double tr_con_2 = cons.sim_data[con_index][pt.troughs[1]];
		double peak_con = cons.sim_data[con_index][pt.peaks[0]];
		*first_p2t = peak_con / ((tr_con_1 + tr_con_2) / 2);
	}
	
	
	int mid_peak = pt.peaks.size() / 2;
	int mid_trough = pt.troughs.size() / 2;
	if (mid_peak == mid_trough && pt.peaks.size() > pt.troughs.size()){
		double peak_con = cons.sim_data[con_index][pt.peaks[mid_peak]];
		double tr_con_1 = cons.sim_data[con_index][pt.troughs[mid_trough]];
		double tr_con_2 = cons.sim_data[con_index][pt.troughs[mid_trough - 1]];
		*mid_p2t = peak_con / ((tr_con_1 + tr_con_2) / 2);
	}
	else if (mid_peak == mid_trough && pt.peaks.size() < pt.troughs.size()){
		double tr_con = cons.sim_data[con_index][pt.troughs[mid_trough]];
		double peak_1 = cons.sim_data[con_index][pt.peaks[mid_peak]];
		double peak_2 = cons.sim_data[con_index][pt.peaks[mid_peak - 1]];
		*mid_p2t = ((peak_1 + peak_2) / 2) / tr_con;
	}
	else if (mid_peak == mid_trough && pt.peaks.size() == pt.troughs.size() && (pt.peaks.size() % 2 == 0)){
		double peak_1 = cons.sim_data[con_index][pt.peaks[mid_peak - 1]];
		double peak_2 = cons.sim_data[con_index][pt.peaks[mid_peak]];
		double tr_1 = cons.sim_data[con_index][pt.troughs[mid_trough - 1]];
		double tr_2 = cons.sim_data[con_index][pt.troughs[mid_trough]];
		*mid_p2t = ((peak_1 + peak_2) / 2) / ((tr_1 + tr_2) / 2);
	}
	else if (mid_peak == mid_trough && pt.peaks.size() == pt.troughs.size() && (pt.peaks.size() % 2 == 1)){
		double peak = cons.sim_data[con_index][pt.peaks[mid_peak]];
		double tr = cons.sim_data[con_index][pt.troughs[mid_trough]];
		*mid_p2t = peak / tr;
	}
	else if (mid_peak > mid_trough){
		double tr_con = cons.sim_data[con_index][pt.troughs[mid_trough]];
		double peak_1 = cons.sim_data[con_index][pt.peaks[mid_peak]];
		double peak_2 = cons.sim_data[con_index][pt.peaks[mid_peak - 1]];
		*mid_p2t = ((peak_1 + peak_2) / 2) / tr_con;
	}
	else if (mid_peak < mid_trough){
		double tr_1 = cons.sim_data[con_index][pt.troughs[mid_trough]];
		double tr_2 = cons.sim_data[con_index][pt.troughs[mid_trough - 1]];
		double peak = cons.sim_data[con_index][pt.peaks[mid_peak]];
		*mid_p2t = peak / ((tr_1 + tr_2) / 2);
	}
	
	if (pt.peaks[pt.peaks.size() - 1] < pt.troughs[pt.troughs.size() - 1]){
		double peak = cons.sim_data[con_index][pt.peaks[pt.peaks.size() - 1]];
		double tr_1 = cons.sim_data[con_index][pt.troughs[pt.troughs.size() - 2]];
		double tr_2 = cons.sim_data[con_index][pt.troughs[pt.troughs.size() - 1]];
		*last_p2t = peak / ((tr_1 + tr_2) / 2);
	}else{
		double peak_1 = cons.sim_data[con_index][pt.peaks[pt.peaks.size() - 2]];
		double peak_2 = cons.sim_data[con_index][pt.peaks[pt.peaks.size() - 1]];
		double tr = cons.sim_data[con_index][pt.troughs[pt.troughs.size() - 1]];
		*last_p2t = ((peak_1 + peak_2) / 2) / tr;
	}
}

double calculate_average_period(peak_trough& pt, input_params& ip){
	double per = 0;
	for (int i = 1 ; i < pt.peaks.size(); i ++){
		per += (double)(pt.peaks[i] - pt.peaks[i - 1]);
	}
	
	per *= (ip.step_size);
	per /= (pt.peaks.size() - 1);
	
	return per;
}

int test_wildtype(con_levels& cons, input_params& ip) {
	int score = 0;
	peak_trough pt;
	find_peaks_troughs(cons, pt, MH1);
	if (pt.peaks.size() >= LOW_BOUND_NUM_PEAKS){
		score += test_sustained_oscillation(pt, cons, MH1);
		double mh1Per = calculate_average_period(pt, ip);
		cout << "mh1 period: " << mh1Per << endl;
		if (LOW_BOUND_PERIOD <= mh1Per && mh1Per <= UP_BOUND_PERIOD){
			score += 1;
		}
	}
	
	
	pt.reset();
	find_peaks_troughs(cons, pt, MH7);
	if (pt.peaks.size() >= LOW_BOUND_NUM_PEAKS){
		score += test_sustained_oscillation(pt, cons, MH7);
		double mh7Per = calculate_average_period(pt, ip);
		cout << "mh7 period: " << mh7Per << endl;
		if (LOW_BOUND_PERIOD <= mh7Per && mh7Per <= UP_BOUND_PERIOD){
			score += 1;
		}
	}
	
	return score;
}
