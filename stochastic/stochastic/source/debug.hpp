#ifndef DEBUG_HPP
#define DEBUG_HPP
#include "structs.hpp"

void test_dependency_graph(sim_data&);
void test_complete_delay(complete_delay&);
void test_embryo(input_params&);
void test_input_processing(input_params&, parameters&);
void test_embryo_concentration(embryo&);
void test_next_firing(embryo&);
void test_deltaK_array(double**);
void test_transfer_record (embryo&);
void test_peaks_troughs();
void test_process_smooth_data();
void test_process_noise_data();
void test_slices();
void test_process_slices(embryo&, features&);
#endif

