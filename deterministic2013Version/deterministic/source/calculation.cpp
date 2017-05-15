#include "calculation.hpp"
#include "macros.hpp"
#include <math.h>
void calculateMH1 (input_params& ip, con_levels& cons, rates& rs, int ci){
	int pi = ci - 1;
	int di = ((int)(ci - rs.data[NMH1]));
	if (di < 0) {
		di = di + cons.num_sim_steps;
		//cout << "G1 di: " << di << "____" << cons.sim_data[G1][di] << "___" << cons.sim_data[PH11][pi] << "___" << cons.sim_data[PD][pi] << endl;
	}
	cons.sim_data[MH1][ci] = cons.sim_data[MH1][pi] 
							+ ip.step_size * (rs.data[MSH1] * ((1 + (cons.sim_data[PD][di] / rs.data[CRITPDH1])) 
							/ (1 + (cons.sim_data[PD][di] / rs.data[CRITPDH1]) + pow(cons.sim_data[PH11][di]/rs.data[CRITPH11H1], 2)))
							- rs.data[MDH1] * cons.sim_data[MH1][pi]);
}

void calculateMH7 (input_params& ip, con_levels& cons, rates& rs, int ci){
	int pi = ci - 1;
	int di = ((int)(ci - rs.data[NMH7]));
	if (di < 0) {
		di = di + cons.num_sim_steps;
	}
	cons.sim_data[MH7][ci] = cons.sim_data[MH7][pi]
							+ ip.step_size * (rs.data[MSH7] * ((1 + (cons.sim_data[PD][di] / rs.data[CRITPDH7]))
							/ (1 + (cons.sim_data[PD][di] / rs.data[CRITPDH7]) + pow(cons.sim_data[PH11][di]/rs.data[CRITPH11H7],2)))
							- rs.data[MDH7] * cons.sim_data[MH7][pi]);
}

void calculateMD (input_params& ip, con_levels& cons, rates& rs, int ci){
	int pi = ci - 1;
	int di = ((int)(ci - rs.data[NMD]));
	if (di < 0) {
		di = di + cons.num_sim_steps;
	}
	cons.sim_data[MD][ci] = cons.sim_data[MD][pi]
							+ ip.step_size * (rs.data[MSD] * (1 / 
							(1 + pow(cons.sim_data[PH11][di] / rs.data[CRITPH11D] , 2)))
							- rs.data[MDD] * cons.sim_data[MD][pi]);
}

void calculatePH1 (input_params& ip, con_levels& cons, rates& rs, int ci){
	int pi = ci - 1;
	int di = ((int)(ci - rs.data[NPH1]));
	if (di < 0) {
		di = di + cons.num_sim_steps;
	}
	cons.sim_data[PH1][ci] = cons.sim_data[PH1][pi]
							+ ip.step_size * (rs.data[PSH1] * cons.sim_data[MH1][di]
							- rs.data[PDH1] * cons.sim_data[PH1][pi]
							- 2 * rs.data[DAH11] * cons.sim_data[PH1][pi] * cons.sim_data[PH1][pi]
							- rs.data[DAH17] * cons.sim_data[PH1][pi] * cons.sim_data[PH7][pi]
							+ 2 * rs.data[DDH11] * cons.sim_data[PH11][pi]
							+ rs.data[DDH17] * cons.sim_data[PH17][pi]);
}

void calculatePH7 (input_params& ip, con_levels& cons, rates& rs, int ci){
	int pi = ci - 1;
	int di = ((int)(ci - rs.data[NPH7]));
	if (di < 0) {
		di = di + cons.num_sim_steps;
	}
	cons.sim_data[PH7][ci] = cons.sim_data[PH7][pi]
							+ ip.step_size * (rs.data[PSH7] * cons.sim_data[MH7][di]
							- rs.data[PDH7] * cons.sim_data[PH7][pi]
							- rs.data[DAH17] * cons.sim_data[PH1][pi] * cons.sim_data[PH7][pi]
							- 2 * rs.data[DAH77] * cons.sim_data[PH7][pi] * cons.sim_data[PH7][pi]
							+ rs.data[DDH17] * cons.sim_data[PH17][pi]
							+ 2 * rs.data[DDH77] * cons.sim_data[PH77][pi]); 
}

void calculatePD (input_params& ip, con_levels& cons, rates& rs, int ci){
	int pi = ci - 1;
	int di = ((int)(ci - rs.data[NPD]));
	if (di < 0) {
		di = di + cons.num_sim_steps;
	}
	cons.sim_data[PD][ci] = cons.sim_data[PD][pi]
							+ ip.step_size * (rs.data[PSD] * cons.sim_data[MD][di]
							- rs.data[PDD] * cons.sim_data[PD][pi]);
}

void calculatePH11 (input_params& ip, con_levels& cons, rates& rs, int ci){
	int pi = ci - 1;
	cons.sim_data[PH11][ci] = cons.sim_data[PH11][pi]
							+ ip.step_size * (rs.data[DAH11] * cons.sim_data[PH1][pi] * cons.sim_data[PH1][pi]
							- rs.data[PDH11] * cons.sim_data[PH11][pi]
							- rs.data[DDH11] * cons.sim_data[PH11][pi]);
}

void calculatePH17 (input_params& ip, con_levels& cons, rates& rs, int ci){
	int pi = ci - 1;
	cons.sim_data[PH17][ci] = cons.sim_data[PH17][pi]
							+ ip.step_size * (rs.data[DAH17] * cons.sim_data[PH1][pi] * cons.sim_data[PH7][pi]
							- rs.data[PDH17] * cons.sim_data[PH17][pi]
							- rs.data[DDH17] * cons.sim_data[PH17][pi]);
}

void calculatePH77 (input_params& ip, con_levels& cons, rates& rs, int ci){
	int pi = ci - 1;
	cons.sim_data[PH77][ci] = cons.sim_data[PH77][pi]
							+ ip.step_size * (rs.data[DAH77] * cons.sim_data[PH7][pi] * cons.sim_data[PH7][pi]
							- rs.data[PDH77] * cons.sim_data[PH77][pi]
							- rs.data[DDH77] * cons.sim_data[PH77][pi]);
}

