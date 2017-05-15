#ifndef TESTS_HPP
#define TESTS_HPP
#include "structs.hpp"
void find_peaks_troughs (con_levels&, peak_trough&, int);
int test_sustained_oscillation(peak_trough&, con_levels&, int);
void calculate_p2t(peak_trough& pt, con_levels&, int, double*, double*, double*);
int test_wildtype(con_levels&, input_params&);
#endif
