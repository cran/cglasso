#include <R.h>

void F77_SUB(glasso_trace1)(int *k, double *rho, int *nit) {
    Rprintf("\nglasso model number %d\n\t  rho = %11.6f\n\tsteps = %4d\n", *k, *rho, *nit);
}

void F77_SUB(glasso_trace2_1)(int *k, double *rho) {
    Rprintf("\n*************************************************************************\n");
    Rprintf("Fitting glasso model number %d with rho = %f\n", *k, *rho);
}

void F77_SUB(glasso_trace2_2)(int *i, int *k) {
    Rprintf("\tconnected component number %d/%d\n", *i, *k);
}

void F77_SUB(glasso_trace2_3_1)(void) {
    Rprintf("\t\tbcd step\tlasso steps\t||Sgm_old - Sgm_new||_inf\n");
}
void F77_SUB(glasso_trace2_3_2)(int *nit, int *lnit, double *dSgm12) {
    Rprintf("\t\t%8d\t%11d\t%25.10f\n", *nit, *lnit, *dSgm12);
}
