#include <cerrno> // Needed for errno, EEXIST
#include <cstdio> // Needed for fopen, fclose, fseek, ftell, rewind
#include <sys/stat.h> // Needed for mkdir
#include <unistd.h> // Needed for read, write, close
#include <stdio.h>
#include <string.h>

#include "print_debug.hpp"
#include "io.hpp"
#include "macros.hpp"
char * states_name[] = {"MH1", "MH7", "MD", "PH1", "PH7", "PD"
			, "PH11", "PH17", "PH77", "G1", "G1N", "G1PH11", "G7", "G7N", "G7PH11", "GD" , "GDPH11"};

char * rates_name [] = {"PSH1", "PSH7", "PSD",\
			"PDH1", "PDH7", "PDD",\
			"MSH1", "MSH7", "MSD",\
			"MDH1", "MDH7", "MDD",\
			"DAH11", "DAH17", "DAH77",\
			"DDH11", "DDH17", "DDH77",\
			"PDH11", "PDH17", "PDH77",\
			"NMH1", "NMH7", "NMD",\
			"NPH1", "NPH7", "NPD",\
			"KAG1PN", "KAG1PH11", "KAG7PN", "KAG7PH11", "KAGDPH11",\
			"KDG1PN", "KDG1PH11", "KDG7PN", "KDG7PH11", "KDGDPH11"};
			
char * reactions_name[] = {"RPSH1", "RPSH7", "RPSD",\
				"RPDH1", "RPDH7", "RPDD",\
				"RPDH11", "RPDH17", "RPDH77",\
				"RDAH11", "RDAH17", "RDAH77",\
				"RDDH11", "RDDH17", "RDDH77",\
				"RMDH1", "RMDH7", "RMDD",\
				"RMSH1", "RMSH1N",\
				"RAG1PH11", "RDG1PH11",\
				"RAG1N", "RDG1N",\
				"RMSH7", "RMSH7N",\
				"RAG7PH11", "RDG7PH11",\
				"RAG7N", "RDG7N",\
				"RMSD",\
				"RAGDPH11", "RDGDPH11"};
				
void open_debug_steam(input_params& ip){
	open_file(&ip.debug_stream, ip.debug_file, true);
	ip.debug_stream << "=====" << endl;
}
				
void print_rates(input_params& ip, rates& rs){
	for (int i = 0; i < NUM_RATES; i ++){
		ip.debug_stream << rates_name[i] << " --- " << rs.data[i] << "\n";
	}
}

void print_propensities(input_params& ip, embryo& em){
	ip.debug_stream << "Propensity" << endl;
	cell* cc;
	for (int i = 0; i < em.num_cells; i++){
		ip.debug_stream << "Cell:" << i << endl;
		cc = em.cell_list[i];
		for (int j = 0; j < NUM_REACTIONS; j ++ ){
			ip.debug_stream << reactions_name[j] << ":" << cc->propen[j] << endl;
		}
	}
}

void print_current_internal(input_params& ip, embryo& em){
	ip.debug_stream << "Current Internal" << endl;
	cell* cc;
	for (int i = 0; i < em.num_cells; i++){
		ip.debug_stream << "Cell:" << i << endl;
		cc = em.cell_list[i];
		for (int j = 0; j < NUM_REACTIONS; j ++ ){
			ip.debug_stream << reactions_name[j] << ":" << cc->current_internal[j] << endl;
		}
	}
}	

void print_next_internal(input_params& ip, embryo& em){
	ip.debug_stream << "Next Internal" << endl;
	cell* cc;
	for (int i = 0; i < em.num_cells; i++){
		ip.debug_stream << "Cell:" << i << endl;
		cc = em.cell_list[i];
		for (int j = 0; j < NUM_REACTIONS; j ++ ){
			ip.debug_stream << reactions_name[j] << ":" << cc->next_internal[j] << endl;
		}
	}
}

void print_states(input_params& ip, embryo& em){
	ip.debug_stream << "States: " << endl;
	cell* cc;
	for (int i = 0; i < em.num_cells; i++){
		ip.debug_stream << "Cell:" << i << endl;
		cc = em.cell_list[i];
		for (int j = 0; j < NUM_STATES; j ++ ){
			ip.debug_stream << states_name[j] << ":" << cc->current_cons[j] << endl;
		}
	}
}

void print_initial_simulation(input_params& ip, embryo& em){
	print_states(ip, em);
	print_current_internal(ip, em);
	print_next_internal(ip, em);
	print_propensities(ip, em);
}

void print_one_round_simulation(embryo& em, input_params& ip, next_reaction& nr){
	ip.debug_stream << "Time: " << em.absolute_time << endl;
	ip.debug_stream << "Next cell:" << nr.cell_index << endl;
	ip.debug_stream << "Next reaction:" << reactions_name[nr.reaction_index] << endl;
	ip.debug_stream << "Delay Complete:" << nr.delay_complete << endl; 
	ip.debug_stream << "Delta:" << nr.delta << endl;
	print_states(ip, em);
	print_current_internal(ip,em);
	print_next_internal(ip, em);
	print_propensities(ip, em);
}

void open_random_steam(input_params& ip){
	open_file(&ip.random_stream, ip.random_file, true);
	ip.random_stream << "=====" << endl;
}

void print_random_number(input_params& ip, double* randNum){
	ip.random_stream << *(randNum) << endl;
}
