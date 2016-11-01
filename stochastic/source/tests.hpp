#ifndef TESTS_HPP
#define TESTS_HPP
#include "structs.hpp"
int test_wildtype(input_params&, embryo&, features&);
void process_smooth_data(embryo&, features&);
int find_peaks_troughs_capacity(embryo&);
void smooth_data(embryo& em, peak_trough&, int, int);
void find_peaks_and_troughs(peak_trough&, features&, int, int);
bool check_peak(peak_trough&, int);
bool check_trough(peak_trough&, int);
void calculate_p2t(peak_trough&, features&, int, int);
void process_binned_data(binned_data&, embryo&, features&);
void find_bin_bounds(double*, embryo&, int);
void put_data_into_bins(double*, embryo&, binned_data&);
void calculate_noise_features(binned_data&, features&);
#endif
