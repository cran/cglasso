#include <R.h>

void F77_SUB(cglasso_trace1)(int *k, double *rho, int *em_nit, int *nit) {
    Rprintf("\ncglasso model number %d\n\t     rho = %11.6f\n\tEM-steps = %4d\n\t   steps = %4d\n", *k, *rho, *em_nit, *nit);
}

void F77_SUB(cglasso_trace2_1)(int *k, double *rho) {
    Rprintf("\n*************************************************************************\n");
    Rprintf("Fitting cglasso model number %d with rho = %f\n", *k, *rho);
}

void F77_SUB(cglasso_trace2_2)(int *i) {
    Rprintf("\t\nstarting EM algorithm step = %d\n", *i);
}

void F77_SUB(cglasso_trace2_3)(void) {
    Rprintf("\n\tE-step: completed!\n");
}

void F77_SUB(cglasso_trace2_4)(void) {
    Rprintf("\tM-step: fitting glasso model\n");
}

void F77_SUB(cglasso_trace2_5)(void) {
    Rprintf("\tM-step: completed!\n");
}

void F77_SUB(cglasso_trace2_6)(double *thr_em, double *dmu, double *dTht) {
    Rprintf("\n\tChecking convergence criterion (threshold = %f)\n", *thr_em);
    Rprintf("\t ||mu_old - mu_new||_inf  = %f\n", *dmu);
    Rprintf("\t||Tht_old - Tht_new||_inf = %f\n", *dTht);
    if(*dmu <= *thr_em && *dTht <= *thr_em) Rprintf("\tConvergence criterion is met!\n");
}
