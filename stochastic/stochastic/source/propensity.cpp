#include "propensity.hpp"
void propensityRPSH1(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPSH1] = rs.data[cell_index][PSH1]
												* (em.cell_list[cell_index]->current_cons)[MH1];
}

void propensityRPSD(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPSD] = rs.data[cell_index][PSD] 
												* (em.cell_list[cell_index]->current_cons)[MD];
}

void propensityRPDH1(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPDH1] = rs.data[cell_index][PDH1] 
												* (em.cell_list[cell_index]->current_cons)[PH1];
}

void propensityRPDD(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPDD] = rs.data[cell_index][PDD] 
												* (em.cell_list[cell_index]->current_cons)[PD];
}

void propensityRPDH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPDH11] = rs.data[cell_index][PDH11] 
												* (em.cell_list[cell_index]->current_cons)[PH11];
}


void propensityRDAH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDAH11] = rs.data[cell_index][DAH11]
												* (em.cell_list[cell_index]->current_cons)[PH1]
												* ((em.cell_list[cell_index]->current_cons)[PH1] - 1) / 2.0;
}


void propensityRDDH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDDH11] = rs.data[cell_index][DDH11]
												* (em.cell_list[cell_index]->current_cons)[PH11];
}


void propensityRMDH1(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMDH1] = rs.data[cell_index][MDH1]
												* (em.cell_list[cell_index]->current_cons)[MH1];
}

void propensityRMDH7(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMDH7] = rs.data[cell_index][MDH7]
												* (em.cell_list[cell_index]->current_cons)[MH7];
}

void propensityRMDD(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMDD] = rs.data[cell_index][MDD]
												* (em.cell_list[cell_index]->current_cons)[MD];
}

void propensityRMSH1(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMSH1] = rs.data[cell_index][MSH1]
												* (em.cell_list[cell_index]->current_cons)[G1];
}
void propensityRMSH1N(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMSH1N] = rs.data[cell_index][MSH1]
												* (em.cell_list[cell_index]->current_cons)[G1N];
}

void propensityRAG1PH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RAG1PH11] = rs.data[cell_index][KAG1PH11]
												* (em.cell_list[cell_index]->current_cons)[G1]
												* (em.cell_list[cell_index]->current_cons)[PH11]
												* ((em.cell_list[cell_index]->current_cons)[PH11] - 1) / 2.0;
}
void propensityRDG1PH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDG1PH11] = rs.data[cell_index][KDG1PH11]
												* (em.cell_list[cell_index]->current_cons)[G1PH11];
}

void propensityRAG1N(embryo& em, int cell_index, rates& rs){
	int neighbor_pd = 0;
	int neighbor_index;
	double num_neighbors = (double) em.neighbor_per_cell;
	for (int i = 0 ; i < em.neighbor_per_cell; i++){
		neighbor_index = em.neighbors[cell_index][i];
		neighbor_pd += (em.cell_list[neighbor_index]->current_cons)[PD];
	}
	(em.cell_list[cell_index]->propen)[RAG1N] = rs.data[cell_index][KAG1PN]
												* (em.cell_list[cell_index]->current_cons)[G1]
												* (1.0/num_neighbors) * neighbor_pd;
}
void propensityRDG1N(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDG1N] = rs.data[cell_index][KDG1PN]
												* (em.cell_list[cell_index]->current_cons)[G1N];
}

void propensityRMSH7(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMSH7] = rs.data[cell_index][MSH7]
												* (em.cell_list[cell_index]->current_cons)[G7];
}
void propensityRMSH7N(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMSH7N] = rs.data[cell_index][MSH7]
												* (em.cell_list[cell_index]->current_cons)[G7N];
}

void propensityRAG7PH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RAG7PH11] = rs.data[cell_index][KAG7PH11]
												* (em.cell_list[cell_index]->current_cons)[G7]
												* (em.cell_list[cell_index]->current_cons)[PH11]
												* ((em.cell_list[cell_index]->current_cons)[PH11] - 1) / 2.0;
}
void propensityRDG7PH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDG7PH11] = rs.data[cell_index][KDG7PH11]
												* (em.cell_list[cell_index]->current_cons)[G7PH11];
}

void propensityRAG7N(embryo& em, int cell_index, rates& rs){
	int neighbor_pd = 0;
	int neighbor_index;
	double num_neighbors = (double) em.neighbor_per_cell;
	for (int i = 0 ; i < em.neighbor_per_cell; i++){
		neighbor_index = em.neighbors[cell_index][i];
		neighbor_pd += (em.cell_list[neighbor_index]->current_cons)[PD];
	}
	(em.cell_list[cell_index]->propen)[RAG7N] = rs.data[cell_index][KAG7PN]
												* (em.cell_list[cell_index]->current_cons)[G7]
												* (1.0/num_neighbors) * neighbor_pd;
}
void propensityRDG7N(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDG7N] = rs.data[cell_index][KDG7PN]
												* (em.cell_list[cell_index]->current_cons)[G7N];
}

void propensityRMSD(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMSD] = rs.data[cell_index][MSD]
												* (em.cell_list[cell_index]->current_cons)[GD];
}

void propensityRAGDPH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RAGDPH11] = rs.data[cell_index][KAGDPH11]
												* (em.cell_list[cell_index]->current_cons)[GD]
												* (em.cell_list[cell_index]->current_cons)[PH11]
												* ((em.cell_list[cell_index]->current_cons)[PH11] - 1) / 2.0;
}
void propensityRDGDPH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDGDPH11] = rs.data[cell_index][KDGDPH11]
												* (em.cell_list[cell_index]->current_cons)[GDPH11];
}

void propensityRAG1NandRAG7N(embryo& em, int cell_index, rates& rs){
	int neighbor_pd = 0;
	int neighbor_index;
	double num_neighbors = (double) em.neighbor_per_cell;
	for (int i = 0 ; i < em.neighbor_per_cell; i++){
		neighbor_index = em.neighbors[cell_index][i];
		neighbor_pd += (em.cell_list[neighbor_index]->current_cons)[PD];
	}
	(em.cell_list[cell_index]->propen)[RAG1N] = rs.data[cell_index][KAG1PN]
												* (em.cell_list[cell_index]->current_cons)[G1]
												* (1.0/num_neighbors) * (double)neighbor_pd;
	(em.cell_list[cell_index]->propen)[RAG7N] = rs.data[cell_index][KAG7PN]
												* (em.cell_list[cell_index]->current_cons)[G7]
												* (1.0/num_neighbors) * (double)neighbor_pd;
}
