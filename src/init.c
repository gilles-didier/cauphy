// Generated with : tools::package_native_routine_registration_skeleton(".", "scr/init.c")
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP getLogDensityTipsCauchy(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getPosteriorLogDensityAncestralCauchy(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getPosteriorLogDensityIncrementCauchy(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP printRTree(SEXP);
extern SEXP SimulateTipsCauchy(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"getLogDensityTipsCauchy",               (DL_FUNC) &getLogDensityTipsCauchy,               7},
    {"getPosteriorLogDensityAncestralCauchy", (DL_FUNC) &getPosteriorLogDensityAncestralCauchy, 8},
    {"getPosteriorLogDensityIncrementCauchy", (DL_FUNC) &getPosteriorLogDensityIncrementCauchy, 8},
    {"printRTree",                            (DL_FUNC) &printRTree,                            1},
    {"SimulateTipsCauchy",                    (DL_FUNC) &SimulateTipsCauchy,                    3},
    {NULL, NULL, 0}
};

void R_init_cauphy(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
