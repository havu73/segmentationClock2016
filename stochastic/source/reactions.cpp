#include "reactions.hpp"
#include "propensity.hpp"
// General reactions:
// 1. Update cell.current_cons[states_index]
// 2. Update cell.last_absolute_change_time[states_index]
// 3. Recalculate propensities

// Delay reactions
// If complete:
//		1,2,3
//		4. pop the reaction from the cdelay priority queue
// If initiated:
//		put the reaction onto the queue

// Reactions RPSD and RPDD require the calculations of propensities of RAG1N adn RAG7N for neighboring cells

void recalculate_propensities(int reaction_index, embryo& em, int cell_index, sim_data& sd, rates& rs){
	int num_dependants = sd.dependency_size[reaction_index]; 
	int dependant;
	for (int i = 0; i < num_dependants; i++){
		dependant = sd.dependency[reaction_index][i];
		sd.prop_funs[dependant](em, cell_index, rs);
	}
}

/*
 * H1 protein synthesis
 */
void reactionRPSH1 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	if (complete){ // delay reaction completed
		// 1
		(em.cell_list[cell_index]->current_cons)[PH1] += 1;
		// 2
		(em.cell_list[cell_index]->last_absolute_change_time)[PH1] = em.cell_list[cell_index]->absolute_time;
		// 3
		recalculate_propensities(RPSH1, em, cell_index, sd, rs);
		// 4
		(em.cell_list[cell_index]->cdelay)->complete_soonest();
	}	
	else{// delay reaction initiated
		double next_complete_time = em.cell_list[cell_index]->absolute_time + NPH1;
		(em.cell_list[cell_index]->cdelay)->initiate_delay(RPSH1, next_complete_time);
	}
}

/*
 * H7 protein synthesis
 */
void reactionRPSH7 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	if (complete){ // delay reaction completed
		// 1
		(em.cell_list[cell_index]->current_cons)[PH7] += 1;
		// 2
		(em.cell_list[cell_index]->last_absolute_change_time)[PH7] = em.cell_list[cell_index]->absolute_time;
		// 3
		recalculate_propensities(RPSH7, em, cell_index, sd, rs);
		// 4
		(em.cell_list[cell_index]->cdelay)->complete_soonest();
	} else{ // delay reaction initiated
		double next_complete_time = em.cell_list[cell_index]->absolute_time + NPH7;
		(em.cell_list[cell_index]->cdelay)->initiate_delay(RPSH7, next_complete_time);
	}
}
void reactionRPSD (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	if (complete){
		// 1
		(em.cell_list[cell_index]->current_cons)[PD] += 1;
		// 2
		(em.cell_list[cell_index]->last_absolute_change_time)[PD] = em.cell_list[cell_index]->absolute_time;
		// 3
		recalculate_propensities(RPSD, em, cell_index, sd, rs);
		// RPSD change the concentrations of PD, so the propensity of RAG1N and RAG7N of neighbors need to be recalculated
		int neigh_index;
		for (int i = 0; i < em.neighbor_per_cell; i++){
			neigh_index = em.neighbors[cell_index][i];
			propensityRAG1NandRAG7N(em, cell_index, rs);
		}
		// 4
		(em.cell_list[cell_index]->cdelay)->complete_soonest();
	} else{
		double next_complete_time = em.cell_list[cell_index]->absolute_time + NPD;
		(em.cell_list[cell_index]->cdelay)->initiate_delay(RPSD, next_complete_time);
	}
}

void reactionRPDH1 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[PH1] -= 1;
	// 2
	(em.cell_list[cell_index]->last_absolute_change_time)[PH1] = em.cell_list[cell_index]->absolute_time;
	// 3
	recalculate_propensities(RPDH1, em, cell_index, sd, rs);
}

void reactionRPDH7 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[PH7] -= 1;
	// 2
	(em.cell_list[cell_index]->last_absolute_change_time)[PH7] = em.cell_list[cell_index]->absolute_time;
	// 3
	recalculate_propensities(RPDH7, em, cell_index, sd, rs);
}

void reactionRPDD (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[PD] -= 1;
	// 2
	(em.cell_list[cell_index]->last_absolute_change_time)[PD] = em.cell_list[cell_index]->absolute_time;
	// 3
	recalculate_propensities(RPDD, em, cell_index, sd, rs);
	// RPDD change the concentrations of PD, so the propensity of RAG1N and RAG7N for neighbors cells need to be recalculated
	int neigh_index;
	for (int i = 0; i < em.neighbor_per_cell; i++){
		neigh_index = em.neighbors[cell_index][i];
		propensityRAG1NandRAG7N(em, cell_index, rs);
	}
}

void reactionRPDH11 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[PH11] -= 1;
	// 2
	(em.cell_list[cell_index]->last_absolute_change_time)[PH11] = em.cell_list[cell_index]->absolute_time;
	// 3
	recalculate_propensities(RPDH11, em, cell_index, sd, rs);
}

void reactionRPDH17 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[PH17] -= 1;
	// 2
	(em.cell_list[cell_index]->current_cons)[PH17] = em.cell_list[cell_index]->absolute_time;
	// 3
	recalculate_propensities(RPDH17, em, cell_index, sd, rs);
}

void reactionRPDH77 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){	
	// 1
	(em.cell_list[cell_index]->current_cons)[PH77] -= 1;
	// 2
	(em.cell_list[cell_index]->current_cons)[PH77] = em.cell_list[cell_index]->absolute_time;
	// 3
	recalculate_propensities(RPDH77, em, cell_index, sd, rs);
}

void reactionRDAH11 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[PH1] -= 2;
	(em.cell_list[cell_index]->current_cons)[PH11] += 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH1] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH11] = time;
	// 3
	recalculate_propensities(RDAH11, em, cell_index, sd, rs);
}

void reactionRDAH17 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[PH1] -= 1;
	(em.cell_list[cell_index]->current_cons)[PH7] -= 1;
	(em.cell_list[cell_index]->current_cons)[PH17] += 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH1] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH7] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH17] = time;
	// 3
	recalculate_propensities(RDAH17, em, cell_index, sd, rs);
}
void reactionRDAH77 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[PH7] -= 2;
	(em.cell_list[cell_index]->current_cons)[PH77] += 1;
	// 2
	double time = (em.cell_list[cell_index])->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH7] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH77] = time;
	// 3
	recalculate_propensities(RDAH77, em, cell_index, sd, rs);
}

void reactionRDDH11 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[PH1] += 2;
	(em.cell_list[cell_index]->current_cons)[PH11] -= 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH1] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH11] = time;
	// 3
	recalculate_propensities(RDDH11, em, cell_index, sd, rs);
}

void reactionRDDH17 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[PH1] += 1;
	(em.cell_list[cell_index]->current_cons)[PH7] += 1;
	(em.cell_list[cell_index]->current_cons)[PH17] -= 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH1] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH7] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH17] = time;
	// 3
	recalculate_propensities(RDDH17, em, cell_index, sd, rs);
}

void reactionRDDH77 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[PH7] += 2;
	(em.cell_list[cell_index]->current_cons)[PH77] -= 1;
	// 2
	double time = (em.cell_list[cell_index])->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH7] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH77] = time;
	// 3
	recalculate_propensities(RDDH77, em, cell_index, sd, rs);
}

void reactionRMDH1 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[MH1] -= 1;
	// 2
	(em.cell_list[cell_index]->last_absolute_change_time)[MH1] = em.cell_list[cell_index]->absolute_time;
	// 3
	recalculate_propensities(RMDH1, em, cell_index, sd, rs);
}
void reactionRMDH7 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[MH7] -= 1;
	// 2
	(em.cell_list[cell_index]->last_absolute_change_time)[MH7] = em.cell_list[cell_index]->absolute_time;
	// 3
	recalculate_propensities(RMDH7, em, cell_index, sd, rs);
}
void reactionRMDD (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[MD] -= 1;
	// 2
	(em.cell_list[cell_index]->last_absolute_change_time)[MD] = em.cell_list[cell_index]->absolute_time;
	// 3
	recalculate_propensities(RMDD, em, cell_index, sd, rs);
}

void reactionRMSH1 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	if (complete){
		// 1
		(em.cell_list[cell_index]->current_cons)[MH1] += 1;
		// 2 
		(em.cell_list[cell_index]->last_absolute_change_time)[MH1] = em.cell_list[cell_index]->absolute_time;
		// 3
		recalculate_propensities(RMSH1, em, cell_index, sd, rs);
		// 4
		(em.cell_list[cell_index]->cdelay)->complete_soonest();
	} else{
		double next_complete_time = em.cell_list[cell_index]->absolute_time + NMH1;
		(em.cell_list[cell_index]->cdelay)->initiate_delay(RMSH1, next_complete_time);
	}
}

void reactionRMSH1N (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	if (complete){
		// 1
		(em.cell_list[cell_index]->current_cons)[MH1] += 1;
		// 2 
		(em.cell_list[cell_index]->last_absolute_change_time)[MH1] = em.cell_list[cell_index]->absolute_time;
		// 3
		recalculate_propensities(RMSH1N, em, cell_index, sd, rs);
		// 4
		(em.cell_list[cell_index]->cdelay)->complete_soonest();
	} else{
		double next_complete_time = em.cell_list[cell_index]->absolute_time + NMH1;
		(em.cell_list[cell_index]->cdelay)->initiate_delay(RMSH1N, next_complete_time);
	}	
}

void reactionRAG1PH11 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[G1] -= 1;
	(em.cell_list[cell_index]->current_cons)[PH11] -= 1;
	(em.cell_list[cell_index]->current_cons)[G1PH11] += 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G1] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH11] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G1PH11] = time;
	// 3
	recalculate_propensities(RAG1PH11, em, cell_index, sd, rs);
}

void reactionRDG1PH11 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[G1] += 1;
	(em.cell_list[cell_index]->current_cons)[PH11] += 1;
	(em.cell_list[cell_index]->current_cons)[G1PH11] -= 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G1] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH11] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G1PH11] = time;
	// 3
	recalculate_propensities(RDG1PH11, em, cell_index, sd, rs);
}

void reactionRAG1N (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[G1] -= 1;
	(em.cell_list[cell_index]->current_cons)[G1N] += 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G1] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G1N] = time;
	// 3
	recalculate_propensities(RAG1N, em, cell_index, sd, rs);
}
void reactionRDG1N (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[G1] += 1;
	(em.cell_list[cell_index]->current_cons)[G1N] -= 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G1] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G1N] = time;
	// 3
	recalculate_propensities(RDG1N, em, cell_index, sd, rs);
}

void reactionRMSH7 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	if (complete){
		// 1
		(em.cell_list[cell_index]->current_cons)[MH7] += 1;
		// 2 
		(em.cell_list[cell_index]->last_absolute_change_time)[MH7] = em.cell_list[cell_index]->absolute_time;
		// 3
		recalculate_propensities(RMSH7, em, cell_index, sd, rs);
		// 4
		(em.cell_list[cell_index]->cdelay)->complete_soonest();
	} else{
		double next_complete_time = em.cell_list[cell_index]->absolute_time + NMH7;
		(em.cell_list[cell_index]->cdelay)->initiate_delay(RMSH7, next_complete_time);
	}
}

void reactionRMSH7N (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	if (complete){
		// 1
		(em.cell_list[cell_index]->current_cons)[MH7] += 1;
		// 2 
		(em.cell_list[cell_index]->last_absolute_change_time)[MH7] = em.cell_list[cell_index]->absolute_time;
		// 3
		recalculate_propensities(RMSH7N, em, cell_index, sd, rs);
		// 4
		(em.cell_list[cell_index]->cdelay)->complete_soonest();
	} else{
		double next_complete_time = em.cell_list[cell_index]->absolute_time + NMH7;
		(em.cell_list[cell_index]->cdelay)->initiate_delay(RMSH7N, next_complete_time);
	}
}

void reactionRAG7PH11 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[G7] -= 1;
	(em.cell_list[cell_index]->current_cons)[PH11] -= 1;
	(em.cell_list[cell_index]->current_cons)[G7PH11] += 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G7] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH11] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G7PH11] = time;
	// 3
	recalculate_propensities(RAG7PH11, em, cell_index, sd, rs);
}
void reactionRDG7PH11 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[G7] += 1;
	(em.cell_list[cell_index]->current_cons)[PH11] += 1;
	(em.cell_list[cell_index]->current_cons)[G7PH11] -= 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G7] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH11] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G7PH11] = time;
	// 3
	recalculate_propensities(RDG7PH11, em, cell_index, sd, rs);
}

void reactionRAG7N (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[G7] -= 1;
	(em.cell_list[cell_index]->current_cons)[G7N] += 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G7] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G7N] = time;
	// 3
	recalculate_propensities(RAG7N, em, cell_index, sd, rs);
}
void reactionRDG7N (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[G7] += 1;
	(em.cell_list[cell_index]->current_cons)[G7N] -= 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G7] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[G7N] = time;
	// 3
	recalculate_propensities(RDG7N, em, cell_index, sd, rs);
}

void reactionRMSD (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	if (complete){
		// 1
		(em.cell_list[cell_index]->current_cons)[MD] += 1;
		// 2 
		(em.cell_list[cell_index]->last_absolute_change_time)[MD] = em.cell_list[cell_index]->absolute_time;
		// 3
		recalculate_propensities(RMSD, em, cell_index, sd, rs);
		// 4
		(em.cell_list[cell_index]->cdelay)->complete_soonest();
	} else{
		double next_complete_time = em.cell_list[cell_index]->absolute_time + NMD;
		(em.cell_list[cell_index]->cdelay)->initiate_delay(RMSD, next_complete_time);
	}
}

void reactionRAGDPH11 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[GD] -= 1;
	(em.cell_list[cell_index]->current_cons)[PH11] -= 1;
	(em.cell_list[cell_index]->current_cons)[GDPH11] += 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[GD] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH11] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[GDPH11] = time;
	// 3
	recalculate_propensities(RAGDPH11, em, cell_index, sd, rs);
}

void reactionRDGDPH11 (embryo& em, int cell_index, sim_data& sd, rates& rs, bool complete){
	// 1
	(em.cell_list[cell_index]->current_cons)[GD] += 1;
	(em.cell_list[cell_index]->current_cons)[PH11] += 1;
	(em.cell_list[cell_index]->current_cons)[GDPH11] -= 1;
	// 2
	double time = em.cell_list[cell_index]->absolute_time;
	(em.cell_list[cell_index]->last_absolute_change_time)[GD] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[PH11] = time;
	(em.cell_list[cell_index]->last_absolute_change_time)[GDPH11] = time;
	// 3
	recalculate_propensities(RDGDPH11, em, cell_index, sd, rs);	
}
