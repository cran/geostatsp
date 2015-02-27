#include"geostatsp.h"

int *SparamOpt, *limTypeOpt; // number of params, followed by indices of params to be optimized
double *paramOpt;
// 0 nugget, 1 variance,
// 2 range, 3 shape,
// 4 anisoRatio, 5 ansioAngleRadians,
// 6 boxcox
double *lower, *upper, *parscale, *ndeps;

double *obsCovOpt, *obsCovCopy; // pointer to the data
double *obsForBoxcoxOpt; // pointer three column matrix of Y, logY and boxcox Y

const double *xcoordOpt, *ycoordOpt;
double *corMatOpt, determinants[2];
double boxcoxParamOpt[9] = {1,0,-9,-9,-9,-9,-9,-9,-9};
int anisoOpt, LtypeOpt, boxcoxTypeOpt;
int Nopt[3];// Nobs, 1, Ncov
int NboxcoxOpt[2]; // Nobs, 3

double *totalSsqOpt;
int NforBoxCoxOpt[3];
double  *betaHatOpt, *varBetaHatOpt;
double *LxLyOpt;

int verboseOpt=0;

// objective function for minimizing -2 logL
double maternLogLObj(
		int junk,
		double *paramArg,
		void *ex
		) {
	int zero=0, oneI=1,  maternType = 2, Dparam;
	int NforCopy;
	// copy param to fullParam

	NforCopy = Nopt[0]*(Nopt[1]+Nopt[2]);

	for(Dparam=0;Dparam<SparamOpt[0];++Dparam){
		paramOpt[SparamOpt[Dparam+1]]=parscale[Dparam]*paramArg[Dparam];
	}

	boxcoxParamOpt[2] = paramOpt[6];
	computeBoxCox(
		obsForBoxcoxOpt,
		NboxcoxOpt,
		boxcoxParamOpt,
		boxcoxTypeOpt
	);

	// make a copy of the data
	F77_NAME(dcopy)(&NforCopy,
			obsCovOpt, &oneI,
			obsCovCopy, &oneI);

	maternForL(
		xcoordOpt, ycoordOpt,
		Nopt,corMatOpt,
		paramOpt,
		&anisoOpt,
		&zero,// don't ignore nugget
		&maternType,//chol of variance matrix
		determinants);

	maternLogLGivenChol(
		obsCovCopy,
		Nopt,  // Nobs, 1, Ncov,
		corMatOpt,
		totalSsqOpt, // a Nrep  by 2 matrix
		betaHatOpt, // an Ncov by Nrep matrix
		varBetaHatOpt, // an Ncov by Ncov by Nrep array
		determinants, // detVarHalf, detCholCovInvXcrossHalf
		LxLyOpt);

	logLfromComponents(
				Nopt,
				&boxcoxParamOpt[6],
				1,
				totalSsqOpt,
				determinants,
				&LtypeOpt
		);
//	Rprintf("\n p ");
//	for(Dparam=0;Dparam<SparamOpt[0];++Dparam)
//		Rprintf(" %f ", paramArg[Dparam] );
//	Rprintf("\n pf ");
//	for(Dparam=0;Dparam<7;++Dparam)
//		Rprintf(" %f ", paramOpt[Dparam] );
//	Rprintf("l %f ", logL[0]);

	if(ISNAN(totalSsqOpt[0])){
			Rprintf("\n p ");
			for(Dparam=0;Dparam<SparamOpt[0];++Dparam)
				Rprintf(" %f ", paramArg[Dparam] );
			Rprintf("\n pf ");
			for(Dparam=0;Dparam<7;++Dparam)
				Rprintf(" %f ", paramOpt[Dparam] );
			Rprintf("\nb %f ", boxcoxParamOpt[8]);
			Rprintf("d %f %f\n", determinants[0], determinants[1]	);

			Rprintf("l %f \n", totalSsqOpt[0]);
	}

	return totalSsqOpt[0];
}

void maternLogLgr(
		int junk,
		double *paramArg,
		double *result,
		void *ex
) {
	int Dpar, oneI=1, Nparam;
	double *parHere, deltaPar, *fullGr;

	R_CheckUserInterrupt();

	Nparam = SparamOpt[0];

	parHere = (double *) calloc(Nparam,sizeof(double));
	fullGr = (double *) calloc(6*Nparam,sizeof(double));


	if(verboseOpt) {
		Rprintf("\nGr npars=%d\nopt scale ", Nparam);
		for(Dpar=0;Dpar < Nparam;++Dpar){
			Rprintf("%f ", paramArg[Dpar]);
		}
		Rprintf("\nnatural scale ");
		for(Dpar=0;Dpar < Nparam;++Dpar){
			Rprintf("%f ", paramArg[Dpar]* parscale[Dpar]);
		}
		Rprintf("\n");
	}

	for(Dpar=0;Dpar < Nparam;++Dpar){
		deltaPar = ndeps[Dpar];

		if(verboseOpt){
	  		Rprintf("p%d=%f delta=%f bnd=%d lb=%f ub=%f\n",
					Dpar,
					paramArg[Dpar],
					deltaPar,
					limTypeOpt[Dpar],
					lower[Dpar], upper[Dpar]);
		}
		//lower
		F77_NAME(dcopy)(&Nparam,
				paramArg, &oneI,
				parHere, &oneI);

		parHere[Dpar]= paramArg[Dpar]-deltaPar;
		if( (limTypeOpt[Dpar] == 1) | (limTypeOpt[Dpar] == 2)){
			parHere[Dpar] = fmax(
								parHere[Dpar],
								lower[Dpar]
						);
		}
		fullGr[Dpar+1*Nparam]=parHere[Dpar];

		fullGr[Dpar+2*Nparam] = maternLogLObj(
				junk,parHere, ex);
		if(verboseOpt){
  		Rprintf("lp=%f lf=%f ",
  				parHere[Dpar],
				fullGr[Dpar+2*Nparam]);
		}

		//upper
		parHere[Dpar]=paramArg[Dpar]+deltaPar;
		if( (limTypeOpt[Dpar] == 3) | (limTypeOpt[Dpar] == 2) ){
			parHere[Dpar] = fmin(
					parHere[Dpar],
					upper[Dpar]
			);
		}
		fullGr[Dpar+3*Nparam]=parHere[Dpar];

		fullGr[Dpar+4*Nparam] = maternLogLObj(
				junk,parHere, ex);
		fullGr[Dpar+5*Nparam] =fullGr[Dpar+4*Nparam] -
				fullGr[Dpar+2*Nparam];
		result[Dpar] = fullGr[Dpar+5*Nparam]/
				(fullGr[Dpar+3*Nparam] - fullGr[Dpar+1*Nparam]);

		if(verboseOpt){
  		Rprintf("up=%f uf=%f gr=%f\n",
  				parHere[Dpar],
				fullGr[Dpar+4*Nparam],
				result[Dpar] );
		}
	}

	if(verboseOpt)	Rprintf("\n");

	free(parHere);
	free(fullGr);

}


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
) {

	double *paramArg, result, *resultGr;
	int oneI=1,Dparam, DparamForOpt, optimFail;
	int fncount, grcount;
	void *optimEx = &oneI;
	char themsg[100];

	// assign the global variables (which end in Opt)
	xcoordOpt=xcoord;
	ycoordOpt=ycoord;

	Nopt[0]	= N[0];
	Nopt[1]	= 1;
	Nopt[2]	= N[2];

	paramOpt=fullParam;
	anisoOpt=*aniso;


	SparamOpt = (int *) calloc(7,sizeof(int));

	// create a vector of indices of parameters to estimate
	DparamForOpt = 1;
	for(Dparam=0;Dparam<7;++Dparam){
		if(Sparam[Dparam]){
			SparamOpt[DparamForOpt] = Dparam;
			++DparamForOpt;
		}
	}
	// first entry is number of parameters to be estimated
	SparamOpt[0]=DparamForOpt-1;

	paramArg = (double *) calloc(SparamOpt[0],sizeof(double));


	// create vector of values for parameters to be estimated
	for(Dparam=0;Dparam<SparamOpt[0];++Dparam){
		paramArg[Dparam] = fullParam[SparamOpt[Dparam+1]];
	}

	// for optim
	lower = parLim;
	upper = &parLim[SparamOpt[0]];
	parscale = &parLim[SparamOpt[0]*2];
	ndeps = &parLim[SparamOpt[0]*3];
	limTypeOpt = limType;


	// prepare boxcox stuff
	// first two columns of obscov are Y and log(Y)
	// likelihood not computed for them.
	obsCovOpt = &obsCov[ 2*N[0] ];
	// but they are used for computing box-cox transform
	obsForBoxcoxOpt = obsCov;
	// are we doing box-cox?


	boxcoxParamOpt[0] = 1;
	boxcoxParamOpt[1] = 0;
	boxcoxParamOpt[2] = fullParam[6];
	NboxcoxOpt[0] = N[0];
	NboxcoxOpt[1] = 3;
	// put log(y) in position 2
	computeBoxCox(
			obsCov,
			NboxcoxOpt,
			boxcoxParamOpt,
			2
	);

	if(Sparam[6]){ // yes, it's being optimized
		boxcoxTypeOpt=3;
	} else { // not being optimized
		boxcoxTypeOpt=4;
	}


// allocate memory

	betaHatOpt= (double *) calloc(N[2],sizeof(double));
	varBetaHatOpt = (double *) calloc(N[2]*N[2],sizeof(double));
	LtypeOpt = *Ltype;
	obsCovCopy = (double *) calloc(N[0]*(N[1]+N[2]),sizeof(double));
	// enough memory for covariates and one vector of observations
	corMatOpt = (double *) calloc(N[0]*N[0],sizeof(double));
	LxLyOpt = (double *) calloc(N[1]*N[2],sizeof(double));
	totalSsqOpt = (double *) calloc(2*N[1],sizeof(double));


#ifdef UNDEF
	resultGr = (double *) calloc(SparamOpt[0]*6,sizeof(double));

	result = maternLogLObj(junk,paramArg, optimEx);
	scalarsF[0] = result;
	scalarsF[1] = totalVarHatOpt[0];
	scalarsF[2] =maternLogLObj(junk,paramArg, optimEx);

	maternLogLgr(
			junk,
			paramArg,
			resultGr,
			optimEx);

	for(Dparam=0;Dparam<SparamOpt[0]*6;++Dparam){
		parLim[Dparam] = resultGr[Dparam];
	}
	free(resultGr);
# endif

//#ifdef UNDEF

	if(scalarsInt[0]){
		verboseOpt=1;
	} else {
		verboseOpt = 0;
	}

// scale everything by parscale
	for(Dparam=0;Dparam<SparamOpt[0];++Dparam){
		paramArg[Dparam] = paramArg[Dparam]/parscale[Dparam];
		lower[Dparam] = lower[Dparam]/parscale[Dparam];
		upper[Dparam] = upper[Dparam]/parscale[Dparam];
	}

	lbfgsb(
			SparamOpt[0],//int n,
			scalarsInt[4],//int lmm,
			paramArg,//double *x,
			lower,
		    upper,
			limType,//int *nbd,
			&result,//double *Fmin,
			maternLogLObj,//optimfn fn,
			maternLogLgr,
			&optimFail,//int *fail,
			optimEx,//void *ex,
			scalarsF[6],//double factr,
			scalarsF[7],//double pgtol,
			&fncount,
			&grcount,
			scalarsInt[1],//int maxit,
			themsg,//char *msg,
			scalarsInt[0],//int trace,
			scalarsInt[2]//int nREP
		  );


// run once again to ensure the global objects
	// reflect the optimal parameters
	resultGr = (double *) calloc(SparamOpt[0],sizeof(double));
	maternLogLgr(
			*N,
			paramArg,
			resultGr,
			optimEx
	);

	maternLogLObj(*N,paramArg, optimEx);

	for(Dparam=0;Dparam<SparamOpt[0];++Dparam){
		paramArg[Dparam] = paramArg[Dparam]*parscale[Dparam];
	}

	strcpy(*msg, themsg);

	scalarsInt[0] = optimFail;
	scalarsInt[1] = fncount;
	scalarsInt[2] = grcount;
	scalarsF[0] = result;
	scalarsF[1] = totalSsqOpt[1];
	scalarsF[2] = boxcoxParamOpt[6];
	scalarsF[3] = boxcoxParamOpt[7];
	scalarsF[4] = boxcoxParamOpt[8];
	scalarsF[5] = determinants[0];
	scalarsF[6] = determinants[1];


	// make a copy of betahat
	F77_NAME(dcopy)(&N[2],
			betaHatOpt, &oneI,
			parLim, &oneI);
	// and varBetaHat
	Dparam = N[2]*N[2];
	F77_NAME(dcopy)(&Dparam,
			varBetaHatOpt, &oneI,
			&parLim[N[2]], &oneI);
	// and gradient
	F77_NAME(dcopy)(&SparamOpt[0],
			resultGr, &oneI,
			&parLim[N[2]+Dparam], &oneI);


//# endif
	free(resultGr);
	free(corMatOpt);
	free(obsCovCopy);
	free(LxLyOpt);
	free(paramArg);
	free(SparamOpt);
	free(totalSsqOpt);
	free(betaHatOpt);
	free(varBetaHatOpt);
}

