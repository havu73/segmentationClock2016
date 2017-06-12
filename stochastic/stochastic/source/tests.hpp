#ifndef TESTS_HPP
#define TESTS_HPP
#include "structs.hpp"
double test_wildtype(input_params&, embryo&, features&, int);
int test_delta(input_params&, embryo&, features&, features&, int);
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
double test_wt_cvs(input_params&, features&);
#endif
