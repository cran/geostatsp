#include"geostatsp.h"


void computeBoxCox(
		double *obsCov,
		// number of observations, number of datasets
		const int *N, // Nobs, Nrep
		double *boxcox, //   3 by Nrep
		// c 1: boxcox par; c 2: sum log(Y), c 3: two log jacobian
		const int boxcoxType
		// 0= do nothing
		// 1= do box-cox
		// 2= put y in 1st column and log(y) in 2nd column of obscov
		// 3= same as 2 but 2nd column of Y and
		//    cols 2,3 of boxcox already computed
		// 4= everything's pre-computed
) {
//
	int D, Dbc, Nend;
	const int Nobs= N[0], Nrep= N[1];
	double *pLogY, *pRep, bcHere, bcEps, sumLogY;
	if( (boxcoxType==0) | (boxcoxType==4) ){
		return;
	}

	bcEps = 0.0005;

	if(boxcoxType == 1){
		pLogY = obsCov; // logs go in first column
		Nend=-1;
	} else {
		pLogY = &obsCov[Nobs]; // logs go in 2nd col
		Nend = 1;
	}

	if(boxcoxType < 3){ // haven't precomputed log Y , sum log Y
		sumLogY = 0.0;
		for(D=0;D<Nobs;++D) {
			pLogY[D] = log(obsCov[D]);
			sumLogY += pLogY[D];
		}
		// two log jacobian
		// to be subtracted from likelihood
		for(D=0;D<Nrep;++D){
			boxcox[D+Nrep] = sumLogY;
			boxcox[D+2*Nrep] = -2*(boxcox[D]-1)*sumLogY;
		}
	} else {
		sumLogY = boxcox[Nrep+1];
		for(D=2;D<Nrep;++D){
			boxcox[D+Nrep] = sumLogY;
			boxcox[D+2*Nrep] = -2*(boxcox[D]-1)*sumLogY;
		}
	}

	Dbc = Nrep-1;
	while(Dbc>Nend){
		pRep = &obsCov[Dbc*Nobs];
		bcHere = boxcox[Dbc];

		if(fabs(bcHere-1) < bcEps) {
			// bc is 1, no transform
			for(D=0;D<Nobs;++D) {
				pRep[D] = obsCov[D];
			}
			// set jacobian to zero
			boxcox[Dbc+2*Nrep] = 0.0;
		} else if(fabs(bcHere) > bcEps) {
			for(D=0;D<Nobs;++D) {
				pRep[D] = (exp(bcHere*pLogY[D]) - 1) /
					bcHere;
			}
		} else {
			for(D=0;D<Nobs;++D) {
				pRep[D] = pLogY[D];
			}
		}
		Dbc--;
	}
} // end box-cox


void maternLogLGivenChol(
		double *obsCov,
		const int *N,  // Nobs, Nrep, Ncov,
		const double *cholVariance,
		double *totalSsq, // an Nrep by 2 matrix
		// first column Y Vinv Y
		// second column Y Vinv X beta
		double *betaHat, // an Ncov by Nrep matrix
		double *varBetaHat, // an Ncov by Ncov by Nrep array
		double *determinants, // detVarHalf, detCholCovInvXcrossHalf
		double *LxLy) {

	const int *Nobs= &N[0], *Nrep= &N[1], *Ncov= &N[2];
	int D, Ncol, infoCholCovInvXcross, infoInvCholCovInvXcross,
		oneInt;
	double zero, one, *pCov, *cholCovInvXY;
	double *cholCovInvXcross;

	oneInt = 1;
	one=1.0;
	zero = 0.0;


	Ncol = *Ncov + *Nrep;


	cholCovInvXcross = varBetaHat;
	cholCovInvXY = obsCov;
	// cholCovInvXY = cholCovMat^{-1} %*% cbind(obs, covariates)


	// solve L x = b
	//      left, lower, not transposed, not diagonal
	F77_NAME(dtrsm)(
			"L", "L", "N","N",
			Nobs, &Ncol,
			&one,
			cholVariance,
			Nobs,
			cholCovInvXY,
			Nobs);

	//  cholCovInvXcross = Matrix::crossprod(cholCovInvX)
// transpose A, don't transpose B,
	pCov = &cholCovInvXY[*Nrep*(*Nobs)];
	//C :=
	//      alpha*op( A )*op( B ) + beta*C,
	F77_NAME(dgemm)("T", "N",
			// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
	  		Ncov, Ncov, Nobs,
			// alpha
			&one,
			// A, nrow(A)
	  		pCov, Nobs,
			// B, nrow(B)
	  		pCov, Nobs,
			// beta
	  		&zero	,
			// C, nrow(c)
			cholCovInvXcross, Ncov);

	//cholCovInvXcrossInv =       Matrix::solve(cholCovInvXcross)
	//detCholCovInvXcross = Matrix::determinant(cholCovInvXcross)$modulus
	//    A = L  * L**T,  if UPLO = 'L'
	// lower or upper, dim, A, nrow, info
	F77_CALL(dpotrf)("L", Ncov, cholCovInvXcross, Ncov, &infoCholCovInvXcross);
	// cholCovInvXcross is now cholesky of cholCovInvXcross
	determinants[1]=0.0;  // the log determinant
	for(D = 0; D < *Ncov; D++)
		determinants[1] += log(cholCovInvXcross[D*(*Ncov)+D]);

	//  DPOTRI computes the inverse of a real symmetric positive definite
	//  matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
	//  computed by DPOTRF.
	// L or U, dim, L, ncol
	F77_NAME(dpotri)("L", Ncov,
			cholCovInvXcross, Ncov,
			&infoInvCholCovInvXcross);
	// cholCovInvXcross is now cholCovInvXcrossInv

	//betaHat = as.vector(
	//      cholCovInvXcrossInv %*%
	//	 Matrix::crossprod(cholCovInvX, cholCovInvY))

	// LxLy=crossprod(cholCovInvX, cholCovInvY)
	// = t(L^(-1) X) L^(-1) Y
	F77_NAME(dgemm)(
			//	op(A), op(B),
			"T", "N",
			// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
	  		Ncov, Nrep, Nobs,
			// alpha
			&one,
			// A, nrow(A)
	  		pCov, Nobs,
			// B, nrow(B)
			cholCovInvXY, Nobs,
			// beta
	  		&zero,
			// C, nrow(c)
			LxLy, Ncov);
	// LxLy is Ncov by Nrep

	// betaHat = cholCovInvXcrossInv %*% LxLy

	F77_NAME(dsymm)(
			// Left or Right, lower ur upper
			"L", "L",
			// nrows of A, ncol ob(B)
	  		Ncov, Nrep,
			// alpha
			&one,
			// A, nrow(A)
			cholCovInvXcross, Ncov,
			// B, nrow(B)
			LxLy, Ncov,
			// beta
	  		&zero,
			// C, nrow(c)
			betaHat, Ncov);


	// totalSsq Y Vinv Y and
	// Y Vinv X betaHat  = Ly' Lx betaHat


	for(D=0;D<*Nrep;++D) {
		totalSsq[D] = F77_NAME(ddot)(
				Nobs,
				&cholCovInvXY[*Nobs*D], &oneInt,
				&cholCovInvXY[*Nobs*D], &oneInt
				);

	// LxLy is Ncov by Nrep
		totalSsq[*Nrep+D] =
		F77_NAME(ddot)(
	  		Ncov,
			&LxLy[*Ncov*D], &oneInt,
			&betaHat[*Ncov*D], &oneInt);
	}

}

// add addToDiag to varMat then compute likelihood
void maternLogLGivenVarU(
		double *varMat, // variance matrix
		const double *varDiag, // new entries for diagnoal
		double *obsCov, // Y, X
		const int *N, // Nobs, Nrep, Ncov,
		double *totalSsq, // length Nrep
		double *betaHat, // an Ncov by Nrep matrix
		double *varBetaHat, // an Ncov by Ncov by Nrep array
		double *determinants // detVarHalf, detCholCovInvXcrossHalf
		) {

	int D, infoCholVarmat;
	double *LxLy;

	for(D=0;D<N[0];++D){
		// diagonals
		varMat[D*N[0]+D] = *varDiag;
	}

	F77_CALL(dpotrf)("L", N, varMat, N, &infoCholVarmat);

	determinants[0]=0.0;  // the log determinant
	for(D = 0; D < N[0]; D++)
		determinants[0] += log(varMat[D*N[0]+D]);

	LxLy = (double *) calloc(N[1]*(N[2]),sizeof(double));

	maternLogLGivenChol(
			obsCov,
			N,
			varMat,
			totalSsq,
			betaHat, // an Ncov by Nrep matrix
			varBetaHat, // an Ncov by Ncov matrix
			determinants, LxLy
			);
	free(LxLy);
}

// computes the -2 log likelihood
// if there is more than one Y supplied
// with the minimum is the last element
void logLfromComponents(
		const int *N,
		// matrix with three rows, boxcox par, sum log L, and jacobian
		const double *boxcox,
		const int boxcoxType,
		double *totalSsq,// length 2*Nrep,
		// on entry, Y Vinv Y in first column
		//   and Y Vinv X' betaHat in second
		// on exit logL in first column
		//  and totalVarHat in second column
		const double *determinants,
		const int *Ltype
		// 0=ml, var estimated
		// 1=reml, var estimated
		//2=ml, var fixed
		// 3=reml, var fixed
){
	const int Nrep = N[1], Ncov= N[2];
	int Nadj, D;
	double *logL, detOne, Lstart, *totalVarHat, ssqXY;

	logL = totalSsq; // Y Vinv Y
	totalVarHat = &totalSsq[Nrep];

	if( (*Ltype==1) | (*Ltype == 3) ){// reml
		Nadj = N[0]-Ncov;
		detOne = determinants[1];
	} else {  // ml
		Nadj = N[0];
		detOne=0.0;// don't add detCholCovInvXcrossHalf
	}

	Lstart =
	      2 * ( Nadj * M_LN_SQRT_2PI +
	    		determinants[0] + detOne);

	if(*Ltype < 2 ){// var estimated
		Lstart += Nadj;//-Nadj*log(Nadj));
		for(D=0;D<Nrep;++D) {
			ssqXY = (logL[D] - totalVarHat[D])/Nadj;
			logL[D] = Lstart + Nadj*log(ssqXY);
			totalVarHat[D] = ssqXY;
		}
	} else { // var fixed
		for(D=0;D<Nrep;++D) {
			ssqXY = logL[D] - totalVarHat[D];
			logL[D] = Lstart + ssqXY;
			totalVarHat[D] = 1; // should be 1
		}
	}
	if(boxcoxType){
		// two log jacobians in 3rd row
		for(D=0;D<Nrep;++D) logL[D] += boxcox[D+(2*Nrep)];
	}

}


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
) {

	double determinants[2], *corMat, *LxLy, *totalSsq;
	int maternType=2, zero=0, oneI=1;
//////
	computeBoxCox(obsCov,
		N,
		boxcox,
		*boxcoxType);

	corMat = (double *) calloc(N[0]*N[0],sizeof(double));
	LxLy = (double *) calloc(N[1]*(N[2]),sizeof(double));
	totalSsq = (double *) calloc(2*(N[1]),sizeof(double));

	
	/////////
	maternForL(
		xcoord, ycoord,
		N,corMat,
		param,aniso,
		&zero,// don't ignore nugget
		&maternType,//chol of variance matrix
		determinants);

	////////
	maternLogLGivenChol(
		obsCov,
		N,  // Nobs, Nrep, Ncov,
		corMat,
		totalSsq, // a Nrep by 2 matrix
		betaHat, // an Ncov by Nrep matrix
		varBetaHat, // an Ncov by Ncov by Nrep array
		determinants, // detVarHalf, detCholCovInvXcrossHalf
		LxLy);

	free(corMat);
	free(LxLy);

 //////////
	logLfromComponents(
			N,boxcox,
			*boxcoxType,
			totalSsq,
			determinants,
			Ltype
	);

	F77_NAME(dcopy)(&N[1],totalSsq, &oneI,logL, &oneI);
	F77_NAME(dcopy)(&N[1],&totalSsq[N[1]], &oneI,  totalVarHat, &oneI);
	free(totalSsq);

	obsCov[0] = determinants[0];
	obsCov[1] = determinants[1];

	*Ltype = maternType;
}


