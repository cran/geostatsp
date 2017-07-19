#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void matern(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void maternAniso(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void maternArasterBpoints(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void maternLogL(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void maternLogLGivenChol(void *, void *, void *, void *, void *, void *, void *);
extern void maternLogLOpt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void maternRaster(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP gmrfEdge(SEXP, SEXP, SEXP);
extern SEXP gmrfLik(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP maternDistance(SEXP, SEXP, SEXP);
extern SEXP maternPoints(SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"matern",               (DL_FUNC) &matern,                9},
    {"maternAniso",          (DL_FUNC) &maternAniso,          12},
    {"maternArasterBpoints", (DL_FUNC) &maternArasterBpoints, 15},
    {"maternLogL",           (DL_FUNC) &maternLogL,           13},
    {"maternLogLGivenChol",  (DL_FUNC) &maternLogLGivenChol,   7},
    {"maternLogLOpt",        (DL_FUNC) &maternLogLOpt,        13},
    {"maternRaster",         (DL_FUNC) &maternRaster,         13},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"gmrfEdge",       (DL_FUNC) &gmrfEdge,       3},
    {"gmrfLik",        (DL_FUNC) &gmrfLik,        6},
    {"maternDistance", (DL_FUNC) &maternDistance, 3},
    {"maternPoints",   (DL_FUNC) &maternPoints,   3},
    {NULL, NULL, 0}
};

void R_init_geostatsp(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
