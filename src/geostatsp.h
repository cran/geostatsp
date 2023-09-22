/*
 see https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Fortran-character-strings
 */
#define USE_FC_LEN_T

#include<R.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#include<Rmath.h>
#include<R_ext/Lapack.h>
#include<R_ext/Applic.h>
#include<R_ext/Print.h>
#include<R_ext/Utils.h>
#include<Matrix.h>
#include <Rdefines.h>


/* type = 0, distances is a vector
 type = 1, distances is a symmetric matrix,
 only lower triangle is used
 type = 2, as type=1 but return the cholesky
 type = 3, as type=1, but return inverse
 nugget is only used if type >1
 and added to the diagonals

If the cholkesy is performed (type > 1):
 halfLogDet is the log determinant if the cholesky matrix
 type is info from dpotrf
if the precision is computed type is info from dpotrfi
... and if chol of precision is computed, type is from dtrtri
 */
void maternAniso(
		const double *x,
		const double *y,
		const int *N,
		double *result,
		const double *range,
		const double *shape,
		const double *variance,
		const double *anisoRatio,
		const double *anisoAngleRadians,
		const double *nugget,
		int *type,
		double *halfLogDet
) ;

void matern(
		const double *distance,
		const int *N,
		double *result,
		const double *range,
		const double *shape,
		const double *variance,
		const double *nugget,
		int *type,
		double *halfLogDet) ;


void maternForL(
		const double *xcoord,
		const double *ycoord,
		const int *N,
		double *corMat,
		const double *param,
		// nugget, variance,
		// range, shape,
		// anisoRatio, ansioAngleRadians
		const int *aniso,
		const int *withoutNugget,
		int *type,
		double *halfLogDet
);



SEXP maternDistance(
		SEXP distance,
		SEXP result,
		SEXP params,
		// range, shape,
		// variance, nugget,
		int *type
		//c('variance','cholesky','precision','inverseCholesky')
);


SEXP maternPoints(
		SEXP points,
		SEXP param,
		// range,
		// shape,
		// variance,
		// anisoRatio,
		// anisoAngleRadians,
		// nugget,
		SEXP result,
		SEXP type);

void maternArasterBpoints(
    double *Axmin, double *Axres,
    int *AxN,
    double *Aymax, double *Ayres, int *AyN,
    double *Bx, double *By, int *BN,
    double *result,
    double  *range, double*shape, double *variance,
    double *anisoRatio, double *anisoAngleRadians);

void computeBoxCox(
		double *obsCov,
		// number of observations, number of datasets
		const int *N, // Nobs, Nrep
		double *boxcox, //  Nrep by 3
		// c 1: boxcox par; c 2: sum log(Y), c 3: two log jacobian
		const int boxcoxType
		// 0= do nothing
		// 1= do box-cox
		// 2= put y in 1st column and log(y) in 2nd column of obscov
		// 3= same as 2 but 2nd column of Y and
		//    rows 2,3 of boxcox already computed
		// 4= everything's pre-computed
);

void maternLogLGivenChol(
		double *obsCov,
		const int *N,  // Nobs, Nrep, Ncov,
		const double *cholVariance,
		double *totalSsq, // a 1 by Nrep matrix
		double *betaHat, // an Ncov by Nrep matrix
		double *varBetaHat, // an Ncov by Ncov by Nrep array
		double *determinants, // detVarHalf, detCholCovInvXcrossHalf
		double *LxLy);

void maternLogLGivenVarU(
		double *varMat, // variance matrix
		const double *varDiag, // new entries for diagnoal
		double *obsCov, // Y, X
		const int *N, // Nobs, Nrep, Ncov,
		double *totalSsq, // length Nrep
		double *betaHat, // an Ncov by Nrep matrix
		double *varBetaHat, // an Ncov by Ncov by Nrep array
		double *determinants // detVarHalf, detCholCovInvXcrossHalf
);

void logLfromComponents(
		const int *N,
		const double *boxcox,
		const int boxcoxType,
		double *totalSsq,// matrix Nrep by 2
		const double *determinants,
		const int *Ltype
		// 0=ml, var estimated
		// 1=reml, var estimated
		//2=ml, var fixed
		// 3=reml, var fixed
);

void maternLogL(
		const double *xcoord,
		const double *ycoord,
		const double *param,// nugget, variance,
        // range, shape,
        // anisoRatio, ansioAngleRadians
		const int *aniso,
		double *obsCov,
		const int *N,// Nobs, Nrep, Ncov
		double *boxcox,
		int *boxcoxType,
		double *logL,// length N[1]
		double *totalVarHat,
		double *betaHat,
		double *varBetaHat,
		int *Ltype
		// 0=ml, var estimated
		// 1=reml, var estimated
		//2=ml, var fixed
		// 3=reml, var fixed
		// on exit, info from chol of matern
);

void maternLogLOpt(
		double *fullParam,// nugget, variance,
        // range, shape,
        // anisoRatio, ansioAngleRadians, boxcox
		int *Sparam,
		// vector of 0 and 1, depending on whether corresponding parameter
		// in fullParam is to be optimized
		double *obsCov,
		const double *xcoord,
		const double *ycoord,
		const int *aniso,
		const int *N,// Nobs, Nrep, Ncov
		int *Ltype,
		int *scalarsInt,
		double *scalarsF,
		double *parLim,
		int *limType,
		char **msg
		// 0=ml, var estimated
		// 1=reml, var estimated
		//2=ml, var fixed
		// 3=reml, var fixed
		// on exit, info from chol of matern
);

SEXP gmrfLik(
    SEXP QR,
    SEXP obsCovR,
    SEXP xisqTausq, // first element must be zero
    SEXP reml,
    SEXP YrepAddR,// Jacobian, length used for Nrep
    SEXP optParam  // set first element -999 for verbose
);

SEXP gmrfEdge(
		SEXP LinvQab, // dense rectangular matrix
		SEXP points, // SpatialPoints*
		SEXP params,
		SEXP result
);

void maternRaster(
    double *Axmin, 
    double *Axres,
    int *AxN,
    double *Aymax, 
    double *Ayres, 
    int *AyN,
    double *result,
    double  *range, 
    double *shape, 
    double *variance,
    double *anisoRatio, 
    double *anisoAngleRadians,
    int *type);

// from Matrix_stubs.c
int attribute_hidden M_cholmod_solve2(
		int sys,
		CHM_FR L,
		CHM_DN B, // right
		CHM_DN *X,// solution
		CHM_DN *Yworkspace,
		CHM_DN *Eworkspace,
		cholmod_common *c);
