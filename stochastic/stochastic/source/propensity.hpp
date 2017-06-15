#ifndef PROPENSITY_HPP
#define PROPENSITY_HPP
#include "structs.hpp"
void propensityRPSH1(embryo&, int, rates&);
void propensityRPSD(embryo&, int, rates&);

void propensityRPDH1(embryo&, int, rates&);
void propensityRPDD(embryo&, int, rates&);
void propensityRPDH11(embryo&, int, rates&);

void propensityRDAH11(embryo&, int, rates&);

void propensityRDDH11(embryo&, int, rates&);

void propensityRMDH1(embryo&, int, rates&);
void propensityRMDH7(embryo&, int, rates&);
void propensityRMDD(embryo&, int, rates&);

void propensityRMSH1(embryo&, int, rates&);
void propensityRMSH1N(embryo&, int, rates&);

void propensityRAG1PH11(embryo&, int, rates&);
void propensityRDG1PH11(embryo&, int, rates&);

void propensityRAG1N(embryo&, int, rates&);
void propensityRDG1N(embryo&, int, rates&);

void propensityRMSH7(embryo&, int, rates&);
void propensityRMSH7N(embryo&, int, rates&);

void propensityRAG7PH11(embryo&, int, rates&);
void propensityRDG7PH11(embryo&, int, rates&);

void propensityRAG7N(embryo&, int, rates&);
void propensityRDG7N(embryo&, int, rates&);

void propensityRMSD(embryo&, int, rates&);

void propensityRAGDPH11(embryo&, int, rates&);
void propensityRDGDPH11(embryo&, int, rates&);

void propensityRAG1NandRAG7N(embryo&, int, rates&);
#endif

