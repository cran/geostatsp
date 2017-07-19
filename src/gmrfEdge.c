#include"geostatsp.h"

SEXP gmrfEdge(
		SEXP LinvQab, // dense rectangular matrix
		SEXP points, // SpatialPoints*
		SEXP params
){

	SEXP result, typePrecision; // dense symmetric
	int Nrow, Ncol;
	double one = 1.0;

	Nrow=INTEGER(getAttrib(
			LinvQab,
			R_DimSymbol))[0];
	Ncol=INTEGER(getAttrib(
			LinvQab,
			R_DimSymbol))[1];

	PROTECT(typePrecision = NEW_CHARACTER(1));
	SET_STRING_ELT(typePrecision, 0, mkChar("precision"));

	PROTECT(result = maternPoints(
			points,
			params,
			typePrecision));

	//	result = crossprod(LinvQab) + result
	// blas DSYRK https://www.math.utah.edu/software/lapack/lapack-blas/dsyrk.html
	F77_NAME(dsyrk)(
			"L","T", &Ncol, &Nrow,
			&one, REAL(LinvQab), &Nrow,
			&one, REAL(GET_SLOT(result, install("x"))), &Ncol
			);

	UNPROTECT(2);
	return result;
}
