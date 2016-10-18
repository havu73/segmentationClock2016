#ifndef REACTIONS_HPP
#define REACTIONS_HPP
#include "structs.hpp"

void reactionRPSH1 (embryo&, int, sim_data&, rates&, bool);
void reactionRPSH7 (embryo&, int, sim_data&, rates&, bool);
void reactionRPSD (embryo&, int, sim_data&, rates&, bool);

void reactionRPDH1 (embryo&, int, sim_data&, rates&, bool);
void reactionRPDH7 (embryo&, int, sim_data&, rates&, bool);
void reactionRPDD (embryo&, int, sim_data&, rates&, bool);
void reactionRPDH11(embryo&, int, sim_data&, rates&, bool);
void reactionRPDH17 (embryo&, int, sim_data&, rates&, bool);
void reactionRPDH77 (embryo&, int, sim_data&, rates&, bool);

void reactionRDAH11 (embryo&, int, sim_data&, rates&, bool);
void reactionRDAH17 (embryo&, int, sim_data&, rates&, bool);
void reactionRDAH77 (embryo&, int, sim_data&, rates&, bool);

void reactionRDDH11 (embryo&, int, sim_data&, rates&, bool);
void reactionRDDH17 (embryo&, int, sim_data&, rates&, bool);
void reactionRDDH77 (embryo&, int, sim_data&, rates&, bool);

void reactionRMDH1 (embryo&, int, sim_data&, rates&, bool);
void reactionRMDH7 (embryo&, int, sim_data&, rates&, bool);
void reactionRMDD (embryo&, int, sim_data&, rates&, bool);

void reactionRMSH1 (embryo&, int, sim_data&, rates&, bool);
void reactionRMSH1N (embryo&, int, sim_data&, rates&, bool);

void reactionRAG1PH11 (embryo&, int, sim_data&, rates&, bool);
void reactionRDG1PH11 (embryo&, int, sim_data&, rates&, bool);

void reactionRAG1N (embryo&, int, sim_data&, rates&, bool);
void reactionRDG1N (embryo&, int, sim_data&, rates&, bool);

void reactionRMSH7 (embryo&, int, sim_data&, rates&, bool);
void reactionRMSH7N (embryo&, int, sim_data&, rates&, bool);

void reactionRAG7PH11 (embryo&, int, sim_data&, rates&, bool);
void reactionRDG7PH11 (embryo&, int, sim_data&, rates&, bool);

void reactionRAG7N (embryo&, int, sim_data&, rates&, bool);
void reactionRDG7N (embryo&, int, sim_data&, rates&, bool);

void reactionRMSD (embryo&, int, sim_data&, rates&, bool);

void reactionRAGDPH11 (embryo&, int, sim_data&, rates&, bool);
void reactionRDGDPH11 (embryo&, int, sim_data&, rates&, bool);
#endif
