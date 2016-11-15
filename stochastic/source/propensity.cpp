#include "propensity.hpp"
void propensityRPSH1(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPSH1] = rs.data[PSH1]
												* (em.cell_list[cell_index]->current_cons)[MH1];
}

void propensityRPSH7(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPSH7] = rs.data[PSH7] 
												* (em.cell_list[cell_index]->current_cons)[MH7];
}

void propensityRPSD(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPSD] = rs.data[PSD] 
												* (em.cell_list[cell_index]->current_cons)[MD];
}

void propensityRPDH1(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPDH1] = rs.data[PDH1] 
												* (em.cell_list[cell_index]->current_cons)[PH1];
}

void propensityRPDH7(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPDH7] = rs.data[PDH7] 
												* (em.cell_list[cell_index]->current_cons)[PH7];
}

void propensityRPDD(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPDD] = rs.data[PDD] 
												* (em.cell_list[cell_index]->current_cons)[PD];
}

void propensityRPDH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPDH11] = rs.data[PDH11] 
												* (em.cell_list[cell_index]->current_cons)[PH11];
}

void propensityRPDH17(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPDH17] = rs.data[PDH17] 
												* (em.cell_list[cell_index]->current_cons)[PH17];
}

void propensityRPDH77(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RPDH77] = rs.data[PDH77] 
												* (em.cell_list[cell_index]->current_cons)[PH77];
}

void propensityRDAH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDAH11] = rs.data[DAH11]
												* (em.cell_list[cell_index]->current_cons)[PH1]
												* ((em.cell_list[cell_index]->current_cons)[PH1] - 1) / 2.0;
}

void propensityRDAH17(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDAH17] = rs.data[DAH17]
												* (em.cell_list[cell_index]->current_cons)[PH1]
												* (em.cell_list[cell_index]->current_cons)[PH7];
}

void propensityRDAH77(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDAH77] = rs.data[DAH77]
												* (em.cell_list[cell_index]->current_cons)[PH7]
												* ((em.cell_list[cell_index]->current_cons)[PH7] - 1) / 2.0;
}

void propensityRDDH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDDH11] = rs.data[DDH11]
												* (em.cell_list[cell_index]->current_cons)[PH11];
}

void propensityRDDH17(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDDH17] = rs.data[DDH17]
												* (em.cell_list[cell_index]->current_cons)[PH17];
}

void propensityRDDH77(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDDH77] = rs.data[DDH77]
												* (em.cell_list[cell_index]->current_cons)[PH77];
}

void propensityRMDH1(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMDH1] = rs.data[MDH1]
												* (em.cell_list[cell_index]->current_cons)[MH1];
}

void propensityRMDH7(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMDH7] = rs.data[MDH7]
												* (em.cell_list[cell_index]->current_cons)[MH7];
}

void propensityRMDD(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMDD] = rs.data[MDD]
												* (em.cell_list[cell_index]->current_cons)[MD];
}

void propensityRMSH1(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMSH1] = rs.data[MSH1]
												* (em.cell_list[cell_index]->current_cons)[G1];
}
void propensityRMSH1N(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMSH1N] = rs.data[MSH1]
												* (em.cell_list[cell_index]->current_cons)[G1N];
}

void propensityRAG1PH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RAG1PH11] = rs.data[KAG1PH11]
												* (em.cell_list[cell_index]->current_cons)[G1]
												* (em.cell_list[cell_index]->current_cons)[PH11];
}
void propensityRDG1PH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDG1PH11] = rs.data[KDG1PH11]
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
	(em.cell_list[cell_index]->propen)[RAG1N] = rs.data[KAG1PN]
												* (em.cell_list[cell_index]->current_cons)[G1]
												* (1.0/num_neighbors) * neighbor_pd;
}
void propensityRDG1N(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDG1N] = rs.data[KDG1PN]
												* (em.cell_list[cell_index]->current_cons)[G1N];
}

void propensityRMSH7(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMSH7] = rs.data[MSH7]
												* (em.cell_list[cell_index]->current_cons)[G7];
}
void propensityRMSH7N(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMSH7N] = rs.data[MSH7]
												* (em.cell_list[cell_index]->current_cons)[G7N];
}

void propensityRAG7PH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RAG7PH11] = rs.data[KAG7PH11]
												* (em.cell_list[cell_index]->current_cons)[G7]
												* (em.cell_list[cell_index]->current_cons)[PH11];
}
void propensityRDG7PH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDG7PH11] = rs.data[KDG7PH11]
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
	(em.cell_list[cell_index]->propen)[RAG7N] = rs.data[KAG7PN]
												* (em.cell_list[cell_index]->current_cons)[G7]
												* (1.0/num_neighbors) * neighbor_pd;
}
void propensityRDG7N(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDG7N] = rs.data[KDG7PN]
												* (em.cell_list[cell_index]->current_cons)[G7N];
}

void propensityRMSD(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RMSD] = rs.data[MSD]
												* (em.cell_list[cell_index]->current_cons)[GD];
}

void propensityRAGDPH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RAGDPH11] = rs.data[KAGDPH11]
												* (em.cell_list[cell_index]->current_cons)[GD]
												* (em.cell_list[cell_index]->current_cons)[PH11];
}
void propensityRDGDPH11(embryo& em, int cell_index, rates& rs){
	(em.cell_list[cell_index]->propen)[RDGDPH11] = rs.data[KDGDPH11]
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
	(em.cell_list[cell_index]->propen)[RAG1N] = rs.data[KAG1PN]
												* (em.cell_list[cell_index]->current_cons)[G1]
												* (1.0/num_neighbors) * neighbor_pd;
	(em.cell_list[cell_index]->propen)[RAG7N] = rs.data[KAG7PN]
												* (em.cell_list[cell_index]->current_cons)[G7]
												* (1.0/num_neighbors) * neighbor_pd;
}
