#include"geostatsp.h"

SEXP gmrfEdge(
		SEXP LinvQab, // dense rectangular matrix
		SEXP points, // SpatialPoints*
		SEXP params,
		SEXP result
){

	int Nrow, Ncol;
	double one = 1.0;

	Nrow=INTEGER(getAttrib(
			LinvQab,
			R_DimSymbol))[0];
	Ncol=INTEGER(getAttrib(
			LinvQab,
			R_DimSymbol))[1];

	maternPoints(
			points,
			result,
			params,
			ScalarInteger(3));

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
