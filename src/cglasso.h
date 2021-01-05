#include <R_ext/RS.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

void
F77_SUB(setup)(int *n, int *p, double *Y, double *lo, double *up,
               double *Yna, int *R, int *startmis, int *order);

void
F77_SUB(fitmcgm)(int *n, int *p, double *Y, double *lo, double *up,
                 int *R, int *nstp, double *eps, double *ym, double *yv,
                 int *conv);

void
F77_SUB(cglasso_v1)(int *n, int *p, double *Y, int *Id, int *nP,
                    int *InfoP, double *lo, double *up, int *pendiag, double *wTht,
                    double *ym, double *yv, int *nrho, double *rhoratio, double *rho,
                    int *maxit_em, double *thr_em, int *maxit_bcd, double *thr_bcd, double *Yipt,
                    double *B, double *mu, double *R, double *S, double *Sgm,
                    double *Tht, int *Adj_yy, int *dfB, int *dfTht, int *ncomp,
                    int *Ck, int *pk, int *nit, int *conv, int *subrout,
                    int *trace);

void
F77_SUB(cglasso_v2)(int *n, int *q, double *X, int *p, double *Y,
                    int *Id, int *nP, int *InfoP, double *lo, double *up,
                    double *wB, int *pendiag, double *wTht, double *ym, double *yv,
                    int *nlambda, double *lambdaratio, double *lambda, int *nrho, double *rhoratio,
                    double *rho, int *maxit_em, double *thr_em, int *maxit_bcd, double *thr_bcd,
                    double *Yipt, double *B, double *mu, double *R, double *S,
                    double *Sgm, double *Tht, int *Adj_yy, int *dfB, int *dfTht,
                    int *ncomp, int *Ck, int *pk, int *Adj_xy, int *nit,
                    int *conv, int *subrout, int *trace);

void
F77_SUB(predict)(double *newrho, double *newlambda, int *nrho, double *rho, int *nlambda,
                 double *lambda, int *d1, int *d2, double *Min, double *Mout);

void
F77_SUB(impute)(double *newrho, double *newlambda, int *nrho, double *rho, int *nlambda,
                double *lambda, int *n, int *p, double *Yin, int *Id,
                double *Yout);

void
F77_SUB(cggm_v1)(int *n, int *p, double *Y, int *Id, int *nP,
                 int *InfoP, double *lo, double *up, double *ym, double *yv,
                 double *pendiag, double *wTht, int *ntp, double *rho, int *maxit_em,
                 double *thr_em, int *maxit_bcd, double *thr_bcd, double *Yipt, double *B,
                 double *mu, double *R, double *S, double *Sgm, double *Tht,
                 int *nit, int *conv, int *subrout, int *trace);

void
F77_SUB(cggm_v2)(int *n, int *q, double *X, int *p, double *Y,
                 int *Id, int *nP, int *InfoP, double *lo, double *up,
                 double *ym, double *yv, double *wB, int *pendiag, double *wTht,
                 int *ntp, double *lambda, double *rho, int *maxit_em, double *thr_em,
                 int *maxit_bcd, double *thr_bcd, double *Yipt, double *B, double *mu,
                 double *R, double *S, double *Sgm, double *Tht, int *nit,
                 int *conv, int *subrout, int *trace);
