#include "debug.hpp"
#include "macros.hpp"
void test_dependency_graph(dependency_graph& dg){
	for (int i = 0; i < NUM_REACTIONS; i++){
		cout << "Reaction " << i << "  ";
		for (int j = 0; j < (dg.dependency_size)[i]; j++){
			cout << (dg.dependency)[i][j] << " ";
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
