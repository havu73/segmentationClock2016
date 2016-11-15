#ifndef TESTS_HPP
#define TESTS_HPP
#include "structs.hpp"
int test_wildtype(input_params&, embryo&, features&);
bool check_mRNA_boundaries (embryo&);
int check_WT_sustained (features&);
int check_WT_amplitude (features&);
int check_WT_avg_cons(features&);
int check_WT_intrinsic(features&);
int check_WT_extrinsic (features&);
void process_smooth_data(embryo&, features&);
void smooth_data(embryo& em, peak_trough&, int, int);
void find_peaks_and_troughs(peak_trough&, features&, int, int);
bool check_peak(peak_trough&, int);
bool check_trough(peak_trough&, int);
void calculate_p2t(peak_trough&, features&, int, int);
void process_noise_data(embryo&, features&);
void process_slices (embryo&, slices&);
void calculate_slice_noise(embryo&, slices&);
void calculate_bounds(slices&, binned_data&, double*);
void put_data_to_bin(slices&, binned_data&, double*);
void calculate_bins(binned_data&, features&);
#endif
