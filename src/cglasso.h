#include <R_ext/RS.h>

void
F77_SUB(setup)(int *n, int *p, double *X, double *lo, double *up,
               int *R, int *startmis);

void
F77_SUB(glasso)(int *p, double *S, double *w, int *pendiag, int *nrho,
                double *rhoratio, double *rho, double *maxR2, int *maxit, double *thr,
                double *Sgm, double *Tht, int *Adj, int *df, double *R2,
                int *ncomp, int *Ck, int *pk, int *nit, int *conv,
                int *trace);

void
F77_SUB(cglasso)(int *n, int *p, double *X, int *R, int *startmis,
                 double *lo, double *up, double *w, double *xm, double *vm,
                 int *nrho, double *rhoratio, double *rho, double *maxR2, int *maxit_em,
                 double *thr_em, int *maxit_bcd, double *thr_bcd, double *Xipt, double *S,
                 double *mu, double *Sgm, double *Tht, int *Adj, int *df,
                 double *R2, int *ncomp, int *Ck, int *pk, int * nit,
                 int *conv, int *subrout, int *trace);

void
F77_SUB(mleglasso)(int *p, double *S, int *nrho, int *maxit, double *thr,
                   double *Sgm, double *Tht, double *R2, int *nit, int *conv,
                   int *trace);

void
F77_SUB(mlecglasso)(int *n, int *p, double * X, int *R, int *startmis,
                    double *lo, double *up, int *nrho, int *maxit_em, double *thr_em,
                    int *maxit_bcd, double *thr_bcd, double *Xipt, double *S, double *mu,
                    double *Sgm, double *Tht, double *R2, int *nit, int *conv,
                    int *subrout, int *trace);
