#include "debug.hpp"
#include "macros.hpp"
#include "tests.hpp"
void test_dependency_graph(sim_data& sd){
	for (int i = 0; i < NUM_REACTIONS; i++){
		cout << "Reaction " << i << "  ";
		for (int j = 0; j < (sd.dependency_size)[i]; j++){
			cout << (sd.dependency)[i][j] << " ";
		}
		cout << endl;
	}
}

void test_complete_delay(complete_delay& cdS){
	cdS.initiate_delay(0, 12.5);
	cdS.initiate_delay(1, 34);
	cdS.initiate_delay(2, 1.45);
	cdS.initiate_delay(3, 56.7);
	cout << "Testing complete_delay priority queue... " ; 
	while (! cdS.is_empty()){
		cout << cdS.see_soonest() << "---";
		cout << cdS.soonest_reaction() << "   ";
		cdS.complete_soonest();
	}
	cout << endl;
}

/* Test embryo checks whether or not the embryo object is function properly.
 * Correct output to the following function (When we comment out the part where neighbor network is checked )should be: 
1  2  3  4  5  0  Just printed concentrations ...
1  2  3  4  5  0  Just printed time...
Top of complete delay priority queue: 2
Size of complete delay: 5
1  2  3  4  5  0  Just printed concentrations ...
1  2  3  4  5  0  Just printed time...
Top of complete delay priority queue: 2
Size of complete delay: 5
Just reset .....
0  0  0  0  0  0  Just printed concentrations ...
0  0  0  0  0  0  Just printed time...
Size of complete delay: 0
0  0  0  0  0  0  Just printed concentrations ...
0  0  0  0  0  0  Just printed time...
Size of complete delay: 0
0  0  0  0  0  0  Just printed concentrations ...
0  0  0  0  0  0  Just printed time...
Size of complete delay: 0
 * 
 */

void test_embryo(input_params& ip){
	embryo em(ip);
	int neighbor_per_cell = 0;
	if (ip.num_cells == 16){
		neighbor_per_cell = SIXTEEN_NEIGHBORS;
	}
	else if (ip.num_cells == 4){
		neighbor_per_cell = FOUR_NEIGHBORS;
	}
	else if (ip.num_cells == 2){
		neighbor_per_cell = TWO_NEIGHBORS;
	}
	for (int i = 0; i < ip.num_cells; i++){
		for (int j = 0; j < neighbor_per_cell; j++){
			cout << (em.neighbors)[i][j] << "  ";
		}
		cout << endl;
	}
	
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 5; j++){
			(*(em.cell_list[i])->cons_record[0])[j] = j+ 1;
			(*(em.cell_list[i])->time_record[0])[j] = j+ 1;
			(((em.cell_list)[i])->cdelay)->initiate_delay(j, j+2);
		}
	}
	
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 6; j++){
			cout << (*(em.cell_list[i])->cons_record)[0][j] << "  ";
		}
		cout << "Just printed concentrations ..." << endl;
		for (int j = 0; j < 6 ; j++){
			cout << (*(em.cell_list[i])->time_record)[0][j] << "  ";
		}
		cout << "Just printed time..." << endl;
		cout << "Top of complete delay priority queue: " << (((em.cell_list)[i])->cdelay)->see_soonest() << endl;
		cout << "Size of complete delay: " << (((em.cell_list)[i])->cdelay)->size() << endl;
	}
	em.reset();
	cout << "Just reset ....." << endl;
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 6; j++){
			cout << (*(em.cell_list[i])->cons_record)[0][j] << "  ";
		}
		cout << "Just printed concentrations ..." << endl;
		for (int j = 0; j < 6 ; j++){
			cout << (*(em.cell_list[i])->time_record)[0][j] << "  ";
		}
		cout << "Just printed time..." << endl;
		cout << "Size of complete delay: " << (((em.cell_list)[i])->cdelay)->size() << endl;
	}
}


void test_input_processing(input_params& ip, parameters& pr){
	for (int i = 0; i < ip.num_sets; i++){
		cout << "Set " << i << "---" ;
		for (int j = 0; j < NUM_RATES; j++){
			cout << pr.data[i][j] << "--";
		}
		cout << endl;
	}
}

void test_embryo_concentration(embryo& em){
	cout << "testing initial conditions" << endl;
	for (int i = 0; i < NUM_STATES; i++){
		cout << ((em.cell_list[0])->current_cons)[i] << "---";
	}
	cout << endl;
}

void test_next_firing(embryo& em){
	cout << "Next firing time: " << endl;
	for (int i = 0; i < NUM_REACTIONS; i++){
		cout << ((em.cell_list[0])->next_internal)[i] << "___";
	}
	cout << endl;
}

/*
 * This function should only be called after em.reset()
 */ 
void test_cdelay_inifinity (embryo& em){
	cout << "Should be: " << NUM_REACTIONS << "----" << ((em.cell_list[0])->cdelay)->soonest_reaction() << endl;
	bool infinite = (((em.cell_list[0])->cdelay)->see_soonest() == INFINITY);
	cout << "Should be true: " << infinite << endl;	
}

void test_deltaK_array(double** deltaKArray){
	for (int i = 0; i < (NUM_REACTIONS + 1); i++){
		if ((*deltaKArray)[i] == INFINITY){
			cout << -1 << "___";
		}else{
			cout << (*deltaKArray)[i] << "___";
		}
	}
	cout << endl;
}

void test_transfer_record (embryo& em){
	int steps = ((em.cell_list[0])->cons_record[0])->size();
	/*
	cout << "size: " << steps << endl;
	cout << "Capacity: " << (em.cell_list[0])->cons_record[0]->capacity() << endl;
	cout << "Max size: " << (em.cell_list[0])->cons_record[0]->max_size() << endl;
	*/
	for (int i = 0; i < steps; i++){
		cout << "concentrations: " << (*(em.cell_list[0])->cons_record[MH1])[i] << endl;
		cout << "time: " << (*(em.cell_list[0])->time_record[MH1])[i] << endl;
	}

}

void test_peaks_troughs(){
	peak_trough pt (5);
	for (int i = 0; i < 5; i ++){
		(pt.smooth_cons)[i] = i;
		(pt.peaks)->push_back(i);
		(pt.troughs)->push_back(i);
	}
	cout << "Size of peaks/troughs: " << (pt.peaks)->size() << endl;
	pt.clear();
	cout << "Size of peaks/ troughs: " << (pt.peaks)->size() << endl;
	for (int i = 0; i < 5; i ++){
		(pt.smooth_cons)[i] = i;
		(*(pt.peaks))[i] = (i + 5);
		(*(pt.troughs))[i] = (i + 5);
	}
	cout << "Size of peaks/troughs (5) : " << (pt.peaks)->size() << endl; 
	cout << "Take away: reset a vector using clear(), not by memset all elements to 0" << endl;
	cout << "After reset, use push_back()  to add elemnts to vector, not by (pt.peaks)[0] = 0." << endl;
	cout << "If you do that, the size of the vector is still 0" << endl;
	for (int i = 0; i < 5; i ++){
		cout << (*(pt.peaks))[i] << "   " ;
	}
	cout << endl;
	
}

/*
 * In order to use this function:
 * change NUM_KEEP_STATES to 1
 * change WINDOW_SIZE to 2
 * change WING_CHECK_SIZE to 2
 */
void test_process_smooth_data(){
	input_params ip;
	ip.num_cells = 1;
	embryo em (ip);
	cell* cc = em.cell_list[0];
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	features fts (1, DEFAULT_NUM_BIN);
	process_smooth_data(em, fts);
	for (int i = 0; i < ip.num_cells; i ++){
		cout << "********* Cell number: " << i + 1 << endl;
		for (int j = 0; j < NUM_KEEP_STATES; j ++){
			cout << "Amplitude: " << fts.avg_amplitude[j][i] << endl;
			cout << "Mid_ptt: " << fts.mid_ptt[j][i] << endl;
			cout << "Last_ptt: " << fts.last_ptt[j][i] << endl;
		}
	}
}

/*
 * In order to use this function:
 */
void test_process_binned_data(){
	input_params ip;
	ip.num_cells = 1;
	embryo em (ip);
	cell* cc = em.cell_list[0];
	// her1
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(9);
	(cc->cons_record[0])->push_back(8);
	(cc->cons_record[0])->push_back(7);
	(cc->cons_record[0])->push_back(6);
	
	// her7
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(9);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(6);
	
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(9);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(6);
	
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(9);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(6);
	
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(9);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(6);
	
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(9);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(6);
	
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(9);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(6);
	
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(9);
	(cc->cons_record[1])->push_back(8);
	(cc->cons_record[1])->push_back(7);
	(cc->cons_record[1])->push_back(6);
	features fts (1, DEFAULT_NUM_BIN);
	binned_data bd(DEFAULT_NUM_BIN);
	process_binned_data(bd, em, fts);
	for (int i = 0; i < DEFAULT_NUM_BIN; i ++){
		cout << "Intrinsic noise: " << fts.intrinsic[i] << endl;
		cout << "Extrinsic noise: " << fts.extrinsic[i] << endl;
		cout << "Avg-cons: " << fts.avg_cons[i] << endl;
		cout << endl;
	}
}
