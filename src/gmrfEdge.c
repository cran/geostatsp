#include"geostatsp.h"

SEXP gmrfEdge(
		SEXP LinvQab, // dense rectangular matrix
		SEXP points, // SpatialPoints*
		SEXP params,
		SEXP result
){

	int Nrow, Ncol, three=3;
	double one = 1.0, halfLogDet;
	int N= Rf_nrows(points);

	Nrow=INTEGER(getAttrib(
			LinvQab,
			R_DimSymbol))[0];
	Ncol=INTEGER(getAttrib(
			LinvQab,
			R_DimSymbol))[1];

//	maternPoints(
//			points,
//			result,
//			params,
//			3);
	
	maternAniso(
	  REAL(points), // x
	  &REAL(points)[N], // y
                &N,
                REAL(GET_SLOT(result, install("x"))),
                &REAL(params)[0],// range,
                &REAL(params)[1],// shape,
                &REAL(params)[2],// variance,
                &REAL(params)[3],// anisoRatio,
                &REAL(params)[4],// anisoAngleRadians,
                &REAL(params)[5],// nugget,
                &three,
                &halfLogDet
	);

	//	result = crossprod(LinvQab) + result
	// blas DSYRK https://www.math.utah.edu/software/lapack/lapack-blas/dsyrk.html
	F77_NAME(dsyrk)(
			"L","T", &Ncol, &Nrow,
			&one, REAL(LinvQab), &Nrow,
			&one, REAL(GET_SLOT(result, install("x"))), 
			&Ncol
			FCONE FCONE);

	return R_NilValue;
}
