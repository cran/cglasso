#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "cglasso.h"

static const R_FortranMethodDef FortEntries[] = {
    {"setup", (DL_FUNC) &F77_SUB(setup), 7},
    {"glasso", (DL_FUNC) &F77_SUB(glasso), 21},
    {"cglasso", (DL_FUNC) &F77_SUB(cglasso), 33},
    {"mleglasso", (DL_FUNC) &F77_SUB(mleglasso), 11},
    {"mlecglasso", (DL_FUNC) &F77_SUB(mlecglasso), 22},
    {NULL, NULL, 0}
};

void attribute_visible R_init_cglasso(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
