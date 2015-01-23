
// gradient

#include"geostatsp.h"

void maternLogLmultipleParams(
		const double *xcoord, const double *ycoord,
		const double *param, // 6 by Npar
		// nugget, variance,
        // range, shape,
        // anisoRatio, ansioAngleRadians
		const int *aniso,
		double *obsCov,
		const int *N, // Nobs, Nrep, Ncov, Npar, NsameVar
		// the first NsameVar columns of param
		// only nugget is different
		double *boxcox,
		int *boxcoxType,
		double *logL,// matrix, Nrep + 1 by Npar
		// last row is the minimum
		const int *Ltype
		// 0=ml, var estimated
		// 1=reml, var estimated
		//2=ml, var fixed
		// 3=reml, var fixed
) {

	double *corMat,*corMatCopy, *obsCovCopy;
	const double *dParam;
	double varDiag,determinants[2];
	double *pLogL, *betaHat, *varBetaHat;
	const int NrowParam=6, Nrep=N[1], Ncov=N[2], Npar = N[3], NsameVar=N[4];
	int Dpar,Nsq, NobsCov;
	int NrepP1;
	int zeroI=0,oneI=1, maternType;

	NrepP1 = Nrep+1;
	Nsq = N[0]*N[0];
	NobsCov = N[0]*(N[1]+N[2]);

	computeBoxCox(
			obsCov,
			N,
			boxcox,
			boxcoxType);

	obsCovCopy = (double *) calloc(NobsCov,sizeof(double));
	corMat = (double *) calloc(Nsq,sizeof(double));
	betaHat = (double *) calloc(Ncov*Nrep,sizeof(double));
	varBetaHat = (double *) calloc(Ncov*Ncov,sizeof(double));

	// matern matrix without nugget
	if(NsameVar){

		maternType=1;
		maternForL(xcoord, ycoord,N,corMat,
			param,aniso,
			&oneI,// ignore nugget
			&maternType,//variance matrix
			determinants);

		corMatCopy = (double *) calloc(Nsq,sizeof(double));

		for(Dpar=NsameVar-1;Dpar>=0;--Dpar){

			if(Dpar){ // copy the variance matrix
			F77_NAME(dcopy)(&Nsq,
					corMat, &oneI,
					corMatCopy, &oneI);
			} else {
				// this is the last time we'll need
				// this variance matrix
				free(corMatCopy);
				corMatCopy = corMat;
			}
			F77_NAME(dcopy)(&NobsCov,
					obsCov, &oneI,
							obsCovCopy, &oneI);

			dParam=&param[NrowParam*Dpar];
			pLogL = &logL[NrepP1*Dpar];

			varDiag = dParam[0]+dParam[1];

			maternLogLGivenVarU(
					corMatCopy, // variance matrix
					&varDiag, // new entries for diagnoal
					obsCovCopy,
					N,
					pLogL, // length Nrep
					betaHat, // an Ncov by Nrep matrix
					varBetaHat, // an Ncov by Ncov by Nrep array
					determinants // detVarHalf, detCholCovInvXcrossHalf
					);

			logLfromComponents(
					N,boxcox,boxcoxType,
					pLogL,pLogL,
					determinants,
					Ltype
			);
		}
	}

	for(Dpar=NsameVar;Dpar<Npar;++Dpar){
		F77_NAME(dcopy)(&NobsCov,
				obsCov, &oneI,
				obsCovCopy, &oneI);

		dParam=&param[NrowParam*Dpar];
		pLogL = &logL[NrepP1*Dpar];

		maternType=2;
		maternForL(xcoord, ycoord,
				N,corMat,
				dParam,aniso,
			&zeroI,// don't ignore nugget
			&maternType,//chol of variance matrix
			determinants);

		maternLogLGivenChol(
				obsCovCopy,
			N,  // Nobs, Nrep, Ncov,
			corMat,
			pLogL, // a 1 by Nrep matrix
			betaHat, // an Ncov by Nrep matrix
			varBetaHat, // an Ncov by Ncov by Nrep array
			determinants // detVarHalf, detCholCovInvXcrossHalf
			);

		logLfromComponents(
				N,boxcox,boxcoxType,
				pLogL,pLogL,
				determinants,
				Ltype
		);

	}
	free(corMat);
	free(obsCovCopy);
	free(varBetaHat);
	free(betaHat);

}
