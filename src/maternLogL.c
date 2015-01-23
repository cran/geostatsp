#include"geostatsp.h"


void computeBoxCox(
		double *obsCov,
		// number of observations, number of datasets
		const int *N, // Nobs, Nrep
		double *boxcox, //  Nrep by 3
		// c 1: boxcox par; c 2: sum log(Y), c 3: two log jacobian
		int *boxcoxType
		// 0= do nothing
		// 1= do box-cox
		// 2= put y in 1st column and log(y) in 2nd column of obscov
		// 3= same as 2 but 2nd column of Y and
		//    rows 2,3 of boxcox already computed
		// 4= everything's pre-computed
) {
//
	int D, Dbc, Nend;
	const int Nobs= N[0], Nrep= N[1];
	double *pLogY, *pRep, bcHere, bcEps, *sumLogY, *twoLogJacobian;
	if(*boxcoxType==0 | *boxcoxType==4){
		return;
	}

	bcEps = 0.01;
	sumLogY = &boxcox[Nrep];
	twoLogJacobian = &boxcox[Nrep*2];

	if(*boxcoxType == 1){
		pLogY = obsCov; // logs go in first column
		Nend=-1;
	} else {
		pLogY = &obsCov[Nobs]; // logs go in 2nd col
		Nend = 1;
	}
	if(*boxcoxType < 3){
		*sumLogY = 0.0;
		for(D=0;D<Nobs;++D) {
			pLogY[D] = log(obsCov[D]);
			*sumLogY += pLogY[D];
		}
		for(D=0;D<Nrep;++D){
			sumLogY[D] = sumLogY[0];
			twoLogJacobian[D] = 2*
					(boxcox[D]-1)*sumLogY[D];
		}
	}

	Dbc = Nrep-1;
	while(Dbc>Nend){
		pRep = &obsCov[Dbc*(Nobs)];
		bcHere = boxcox[Dbc];
		if(abs(bcHere) > bcEps) {
			for(D=0;D<Nobs;++D) {
				pRep[D] = (exp(bcHere*pLogY[D]) - 1) /
					bcHere;
			}
		} else {
			for(D=0;D<Nobs;++D) {
				pRep[D] = pLogY[D];
			}

		}
	}
	*boxcoxType = 4;
} // end box-cox


void maternLogLGivenChol(
		double *obsCov,
		const int *N,  // Nobs, Nrep, Ncov,
		const double *cholVariance,
		double *totalSsq, // a 1 by Nrep matrix
		double *betaHat, // an Ncov by Nrep matrix
		double *varBetaHat, // an Ncov by Ncov by Nrep array
		double *determinants // detVarHalf, detCholCovInvXcrossHalf
		) {

	const int *Nobs= &N[0], *Nrep= &N[1], *Ncov= &N[2];
	int D, Ncol, infoCholCovInvXcross, infoInvCholCovInvXcross,
		NobsRep, NobsCovObs, oneInt;
	double zero, minusone, one, *pCov, *cholCovInvXY;
	double *cholCovInvXcross, detLx, *LxLy, detCholCovInvXcrossHalf;

	detCholCovInvXcrossHalf = determinants[1];
	oneInt = 1;
	one=1.0;
	minusone = -1.0;
	zero = 0.0;


	Ncol = *Ncov + *Nrep;
	NobsCovObs = *Nobs*Ncol;
	NobsRep = *Nobs*(*Nrep);

	cholCovInvXcross = varBetaHat;
	cholCovInvXY = obsCov;
	LxLy = (double *) calloc(*Ncov*(*Nrep),sizeof(double));

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
	detCholCovInvXcrossHalf=0;  // the log determinant
	for(D = 0; D < *Ncov; D++)
		detCholCovInvXcrossHalf += log(cholCovInvXcross[D*(*Ncov)+D]);

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

	// resids = obsCov[,1] - as.vector(obsCov[,-1] %*% betaHat)
	//  cholCovInvResid = Matrix::solve(cholCovMat, resids)

	// C=y, beta = 1, alpha = -1, A=X, B=betaHat
	F77_NAME(dgemm)("N", "N",
			// nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
	  		Nobs, Nrep, Ncov,
			// alpha
			&minusone,
			// A, nrow(A)
			pCov, Nobs,
			// B, nrow(B)
			betaHat, Ncov,
			// beta
	  		&one,
			// C, nrow(c)
			cholCovInvXY, Nobs);

	//  totalSsq = as.vector(Matrix::crossprod(cholCovInvResid))
	//F77_NAME(ddot)(const int *n, const double *dx, const int *incx,
	//	       const double *dy, const int *incy);
	for(D=0;D<*Nrep;++D) {
		totalSsq[D] = F77_NAME(ddot)(
				Nobs,
				cholCovInvXY, &oneInt,
				cholCovInvXY, &oneInt
				);
	}


	free(LxLy);
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

	int D, infoCholVarmat, zeroI=0;

	for(D=0;D<N[0];++D){
		// diagonals
		varMat[D*N[0]+D] = *varDiag;
	}

	F77_CALL(dpotrf)("L", N, varMat, N, &infoCholVarmat);

	determinants[0]=0.0;  // the log determinant
	for(D = 0; D < N[0]; D++)
		determinants[0] += log(varMat[D*N[0]+D]);

	maternLogLGivenChol(
			obsCov,
			N,
			varMat,
			totalSsq,
			betaHat, // an Ncov by Nrep matrix
			varBetaHat, // an Ncov by Ncov matrix
			determinants
			);
}

// computes the -2 log likelihood
// if there is more than one Y supplied
// with the minimum is the last element
void logLfromComponents(
		const int *N,
		const double *boxcox,
		const int *boxcoxType,
		double *totalSsq,// length Nrep+1, logL on exit
		double *totalVarHat,// length Nrep
		// can be the same as totalSsq
		// in which case totalVarHat is lost
		const double *determinants,
		const int *Ltype
		// 0=ml, var estimated
		// 1=reml, var estimated
		//2=ml, var fixed
		// 3=reml, var fixed
){
	const int Nrep = N[1], Ncov= N[2];
	int Nadj, D;
	const double *pBoxCox;
	double *pLogL, *logL, detOne, Lstart;

	logL = totalSsq;

	if(*Ltype==1 | *Ltype == 3){// reml
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
		Lstart += Nadj;
		for(D=0;D<Nrep;++D) {
			totalVarHat[D] = totalSsq[D]/Nadj;
			logL[D] = Lstart + Nadj * log(totalVarHat[D]);
		}
	} else { // var fixed
		for(D=0;D<Nrep;++D) {
			totalVarHat[D] = 1.0;
			logL[D] = Lstart + totalSsq[D];
		}
	}
	if(*boxcoxType){
		if(*boxcoxType != 4)
			warning(
"Box-Cox log jacobians should be computed prior to calling logLfromComponents"
					);
		pBoxCox= &boxcox[Nrep*2]; // two log jacobians in 3rd row
		for(D=0;D<Nrep;++D) logL[D] += pBoxCox[D];
	}

	// find the minimum lf all the likelihoods calculated
	pLogL = &logL[Nrep];
	*pLogL = logL[0];
	for(D=1;D<Nrep;++D){
		*pLogL = fmin2(*pLogL, logL[D]);
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
		double *logL,// length N[1] + 1
		// last element is the minimum
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

	double determinants[2], *pBoxCox, *corMat;
	int maternType=2, zero=0;

	computeBoxCox(obsCov,
		N,
		boxcox,
		boxcoxType);

	corMat = (double *) calloc(N[0]*N[0],sizeof(double));

	maternForL(
			xcoord, ycoord,
		N,corMat,
		param,aniso,
		&zero,// don't ignore nugget
		&maternType,//chol of variance matrix
		determinants);

	maternLogLGivenChol(
		obsCov,
		N,  // Nobs, Nrep, Ncov,
		corMat,
		logL, // a 1 by Nrep matrix
		betaHat, // an Ncov by Nrep matrix
		varBetaHat, // an Ncov by Ncov by Nrep array
		determinants // detVarHalf, detCholCovInvXcrossHalf
		);

	free(corMat);

	logLfromComponents(
			N,boxcox,boxcoxType,
			logL,totalVarHat,
			determinants,
			Ltype
	);
	*Ltype = maternType;
}


