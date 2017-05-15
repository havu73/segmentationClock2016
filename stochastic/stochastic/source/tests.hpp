#ifndef TESTS_HPP
#define TESTS_HPP
#include "structs.hpp"
double test_wildtype(input_params&, embryo&, features&, int);
int test_delta(input_params&, embryo&, features&);
bool check_mRNA_boundaries (embryo&);
double check_WT_sustained (input_params&, features&);
double check_WT_amplitude (input_params&, features&);
double check_WT_period (input_params&, features&);
int check_delta_amplitude (features&);
double check_WT_avg_cons(input_params&, features&);
int check_delta_avg_cons(features&);
double check_WT_intrinsic(input_params&, features&);
int check_delta_intrinsic(features&);
double check_WT_extrinsic (input_params&, features&);
int check_delta_extrinsic (features&);
void process_smooth_data(input_params&, embryo&, features&, int);
void smooth_data(embryo&, peak_trough&, int);
void find_peaks_and_troughs(embryo&, peak_trough&, features&, int, int);
bool check_peak(peak_trough& pt, int, int);
bool check_trough(peak_trough&, int, int);
double calculate_middle_peak_or_trough(double*, vector<int>*);
void calculate_p2t(peak_trough&, features&, int, int);
void cell_to_cell_correlations(peak_trough&, features&, int);
double calculate_correlation(double*, double*, int);
bool check_correlations(features&);
void process_noise_data(embryo&, features&, input_params&);
void process_slices (embryo&, slices&, input_params&);
void calculate_slice_noise(embryo&, slices&);
void calculate_bounds(slices&, binned_data&, double*);
void put_data_to_bin(slices&, binned_data&, double*);
void calculate_bins(binned_data&, features&);
#endif
