#ifndef functions_pi0_H
#define functions_pi0_H
#include "non_linear_fit.hpp"

double lhs_function_pi0_eg(int j, double**** in, int t, struct fit_type fit_info);
double** sub_vev(int j, double**** in, int t, struct fit_type fit_info);
#endif
