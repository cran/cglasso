#include <R.h>

void F77_SUB(mleglasso_trace1)(int *k, int *nit) {
    Rprintf("Fitted ggm model number %d (steps = %4d)\n", *k, *nit);
}

void F77_SUB(mleglasso_trace2)(int *k) {
    Rprintf("\n*************************************************************************\n");
    Rprintf("Fitting ggm model number %d\n", *k);
}
