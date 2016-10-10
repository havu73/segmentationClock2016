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
