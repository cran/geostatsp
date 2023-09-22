
#include"geostatsp.h"
#include<R_ext/Visibility.h>


// .C
// matern
static R_NativePrimitiveArgType matern_t[] = {
    REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP
};

// maternAniso
static R_NativePrimitiveArgType maternAniso_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP
};

// maternArasterBpoints
static R_NativePrimitiveArgType maternArasterBpoints_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP,  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

// maternLogL
static R_NativePrimitiveArgType maternLogL_t[] = {
    REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

// maternLogLopt
static R_NativePrimitiveArgType maternLogLOpt_t[] = {
    REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, STRSXP
};

static R_NativePrimitiveArgType maternRaster_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

static R_NativePrimitiveArgType maternLogLGivenChol_t[] = {
    REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static const R_CMethodDef cMethods[] = {
   {"matern", (DL_FUNC) &matern, 9, matern_t},
   {"maternAniso", (DL_FUNC) &maternAniso, 12, maternAniso_t},
   {"maternRaster", (DL_FUNC) &maternRaster, 13, maternRaster_t},
   {"maternLogLGivenChol", (DL_FUNC) &maternLogLGivenChol, 7, maternLogLGivenChol_t},
   {"maternArasterBpoints", (DL_FUNC) &maternArasterBpoints, 15, maternArasterBpoints_t},
   {"maternLogL", (DL_FUNC) &maternLogL, 13, maternLogL_t},
   {"maternLogLOpt", (DL_FUNC) &maternLogLOpt, 13, maternLogLOpt_t},
   {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[]  = {
  {"maternPoints", (DL_FUNC) &maternPoints, 4},
  {"maternDistance", (DL_FUNC) &maternDistance, 4},
  {"gmrfLik", (DL_FUNC) &gmrfLik, 6},
  {"gmrfEdge", (DL_FUNC) &gmrfEdge, 4},
  {NULL, NULL, 0}
};


void attribute_visible R_init_geostatsp(DllInfo *info) {
	R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
}
