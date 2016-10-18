#include "debug.hpp"
#include "macros.hpp"
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
	
	for (int i = 0; i < ip.num_cells; i++){
		for (int j = 0; j < MAX_NEIGHBORS; j++){
			cout << (em.neighbors)[i][j] << "  ";
		}
		cout << endl;
	}
	
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 5; j++){
			(*((em.cell_list[i])->cons)->cons_record[0])[j] = j+ 1;
			(*((em.cell_list[i])->cons)->time_record[0])[j] = j+ 1;
			(((em.cell_list)[i])->cdelay)->initiate_delay(j, j+2);
		}
	}
	
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 6; j++){
			cout << (*((em.cell_list[i])->cons)->cons_record)[0][j] << "  ";
		}
		cout << "Just printed concentrations ..." << endl;
		for (int j = 0; j < 6 ; j++){
			cout << (*((em.cell_list[i])->cons)->time_record)[0][j] << "  ";
		}
		cout << "Just printed time..." << endl;
		cout << "Top of complete delay priority queue: " << (((em.cell_list)[i])->cdelay)->see_soonest() << endl;
		cout << "Size of complete delay: " << (((em.cell_list)[i])->cdelay)->size() << endl;
	}
	em.reset();
	cout << "Just reset ....." << endl;
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 6; j++){
			cout << (*((em.cell_list[i])->cons)->cons_record)[0][j] << "  ";
		}
		cout << "Just printed concentrations ..." << endl;
		for (int j = 0; j < 6 ; j++){
			cout << (*((em.cell_list[i])->cons)->time_record)[0][j] << "  ";
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
		cout << (((em.cell_list[0])->cons)->current_cons)[i] << "---";
	}
	cout << endl;
}

void test_next_firing(embryo& em){
	cout << "Next firing time: " << endl;
	for (int i = 0; i < NUM_REACTIONS; i++){
		cout << ((em.cell_list[0])->propen)[i] << "___";
	}
	cout << endl;
}
