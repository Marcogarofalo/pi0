#define functions_pi0_C
#include "functions_pi0.hpp"
#include "tower.hpp"

double lhs_function_pi0_eg(int j, double**** in, int t, struct fit_type fit_info) {
    double r = in[j][0][t][0];
    return r;
}

double** sub_vev(int j, double**** in, int t, struct fit_type fit_info) {
    double** r = malloc_2<double>(fit_info.N, 2);
    // printf("N=%d\n",fit_info.N);
    r[0][0] = in[j][fit_info.corr_id[0]][t][0] - fit_info.ext_P[0][j];
    r[0][1] = 0;
    // printf("%g  %g\n", fit_info.ext_P[0][j] ,r[0][1]);
    return r;
}


double** add_connect(int j, double**** in, int t, struct fit_type fit_info) {
    double** r = malloc_2<double>(fit_info.N, 2);
    r[0][0] = in[j][fit_info.corr_id[0]][t][0] + 2 * in[j][fit_info.corr_id[1]][t][0];
    r[0][1] = 0;
    return r;
}