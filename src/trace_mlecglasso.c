#include <R.h>

void F77_SUB(mlecglasso_trace1)(int *k, int *em_nit, int *nit) {
    Rprintf("\nFitted cggm model number %d\n\tEM-steps = %4d\n\t   steps = %4d\n", *k, *em_nit, *nit);
}

void F77_SUB(mlecglasso_trace2_1)(int *k) {
    Rprintf("\n*************************************************************************\n");
    Rprintf("Fitting cggm model number %d\n", *k);
}

void F77_SUB(mlecglasso_trace2_4)(void) {
    Rprintf("\tM-step: fitting cggm model\n");
}
