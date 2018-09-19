#include <R.h>

void F77_SUB(mglasso_trace1)(int *k, double *rho, int *em_nit, int *nit) {
    Rprintf("\nmglasso model number %d\n\t     rho = %11.6f\n\tEM-steps = %4d\n\t   steps = %4d\n", *k, *rho, *em_nit, *nit);
}

void F77_SUB(mglasso_trace2_1)(int *k, double *rho) {
    Rprintf("\n*************************************************************************\n");
    Rprintf("Fitting mglasso model number %d with rho = %f\n", *k, *rho);
}
