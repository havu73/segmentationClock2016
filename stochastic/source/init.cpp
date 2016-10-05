#include <unistd.h> // Needed for getpid

#include "init.hpp" // Function declarations
#include "io.hpp"
#include "propensity.hpp"
using namespace std;

terminal* term = NULL; // The global terminal struct

/* init_terminal creates and initializes a new terminal struct
	parameters:
	returns: nothing
	notes:
	todo:
*/
void init_terminal () {
	if (term != NULL) {
		delete term;
	}
	term = new terminal();
}

/* free_terminal frees the terminal from memory and resets the terminal text color to its default value
	parameters:
	returns: nothing
	notes:
	todo:
*/
void free_terminal () {
	cout << term->reset;
	delete term;
}

/* init_verbosity sets the verbose stream to /dev/null if verbose mode is not enabled
	parameters:
		ip: the program's input parameters
	returns: nothing
	notes:
	todo:
*/
void init_verbosity (input_params& ip) {
	if (!ip.verbose) {
		term->set_verbose_streambuf(ip.null_stream->rdbuf());
	}
}

void init_propensities(propensities& prop){
	prop.prop_funs[RPSH1] = &propensityRPSH1;
	prop.prop_funs[RPSH7] = &propensityRPSH7;
	prop.prop_funs[RPSD] = &propensityRPSD;
	
	prop.prop_funs[RPDH1] = &propensityRPDH1;
	prop.prop_funs[RPDH7] = &propensityRPDH7;
	prop.prop_funs[RPDD] = &propensityRPDD;
	prop.prop_funs[RPDH11] = &propensityRPDH11;
	prop.prop_funs[RPDH17] = &propensityRPDH17;
	prop.prop_funs[RPDH77] = &propensityRPDH77;
	
	prop.prop_funs[RDAH11] = &propensityRDAH11;
	prop.prop_funs[RDAH17] = &propensityRDAH17;
	prop.prop_funs[RDAH77] = &propensityRDAH77;
	
	prop.prop_funs[RDDH11] = &propensityRDDH11;
	prop.prop_funs[RDDH17] = &propensityRDDH17;
	prop.prop_funs[RDDH77] = &propensityRDDH77;
	
	prop.prop_funs[RMDH1] = &propensityRMDH1;
	prop.prop_funs[RMDH7] = &propensityRMDH7;
	prop.prop_funs[RMDD] = &propensityRMDD;
	
	prop.prop_funs[RMSH1] = &propensityRMSH1;
	prop.prop_funs[RMSH1N] = &propensityRMSH1N;
	
	prop.prop_funs[RAG1PH11] = &propensityRAG1PH11;
	prop.prop_funs[RDG1PH11] = &propensityRDG1PH11;
	
	prop.prop_funs[RAG1N] = &propensityRAG1N;
	prop.prop_funs[RDG1N] = &propensityRDG1N;
	
	prop.prop_funs[RMSH7] = &propensityRMSH7;
	prop.prop_funs[RMSH7N] = &propensityRMSH7N;
	
	prop.prop_funs[RAG7PH11] = &propensityRAG7PH11;
	prop.prop_funs[RDG7PH11] = &propensityRDG7PH11;
	
	prop.prop_funs[RAG7N] = &propensityRAG7N;
	prop.prop_funs[RDG7N] = &propensityRDG7N;
	
	prop.prop_funs[RMSD] = &propensityRMSD;
	prop.prop_funs[RAGDPH11] = &propensityRAGDPH11;
	prop.prop_funs[RDGDPH11] = &propensityRDGDPH11;
}
