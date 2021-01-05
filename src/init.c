#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "cglasso.h"

static const R_FortranMethodDef FortEntries[] = {
    {"setup", (DL_FUNC) &F77_SUB(setup), 9},
    {"fitmcgm", (DL_FUNC) &F77_SUB(fitmcgm), 11},
    {"cglasso_v1", (DL_FUNC) &F77_SUB(cglasso_v1), 36},
    {"cglasso_v2", (DL_FUNC) &F77_SUB(cglasso_v2), 43},
    {"predict", (DL_FUNC) &F77_SUB(predict), 10},
    {"impute", (DL_FUNC) &F77_SUB(impute), 11},
    {"cggm_v1", (DL_FUNC) &F77_SUB(cggm_v1), 29},
    {"cggm_v2", (DL_FUNC) &F77_SUB(cggm_v2), 33},
    {NULL, NULL, 0}
};

void attribute_visible R_init_cglasso(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
