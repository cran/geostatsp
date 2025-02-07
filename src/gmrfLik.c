/*
 * needs modified matrix package
 * svn checkout svn://svn.r-forge.r-project.org/svnroot/matrix/pkg/Matrix Matrix
 R CMD build --no-build-vignettes Matrix
 no vignettes  because I don't have texi2pdf
 */

#include"geostatsp.h"
//#include<R.h>
//#include<Rmath.h>
//#include<R_ext/Lapack.h>
//#include<R_ext/Applic.h>
//#include<R_ext/Print.h>
//#include<R_ext/Utils.h>
//#include<R_ext/Rdynload.h>
//#include<Matrix.h> it's in matrix_stubs.c
//#include<Matrix_stubs.c>

double Brent_fmin(double ax, double bx, double (*f)(double, void *),
                  void *info, double tol);


/*
	a function for interfacing to the Matrix package
	hopefully this will be in Matrix_stubs.c soon
 */


/*
 *  global variables, so that an optimizer can be build eventually
 */
CHM_SP Q;
CHM_FR L;
CHM_DN obsCovRot, Lx;
CHM_DN YwkL, EwkL, YwkD, EwkD; // workspaces
cholmod_common c;
double *YXVYXglobal,  *YXYXglobal, *YrepAdd, *betaHatGlobal;
double *copyLx;
int Nxy, Nobs, Ncov, Nrep, Nxysq;
int Ltype;
int DxisqTausq, NxisqTausq;
const double logTwoPi = 1.837877066409345339082; //format(log(2*pi), digits=22)
const double onePlusLogTwoPi = 2.837877066409345339082; //format(log(1+2*pi), digits=22)
/*
 * compute sums of squares from cross products
 *
 */
void ssqFromXprod(
    double *YXVinvYX, // N by N
    double *detXVinvX,
    const int N, // Nxy
    const int Nrep,
    double *copyLx
){

  int oneI=1, infoCholXX, infoInvXX, Ncov;
  int D;
  double	oneD=1.0, moneD = -1.0,zeroD=0.0;
  double *xvx;

  /// copy of LyLx
  Ncov = N*Nrep;
  F77_CALL(dcopy)(&Ncov,
      YXVinvYX, &oneI,
      copyLx, &oneI);

  // xvinvx submatrix
  xvx = &YXVinvYX[N*Nrep+Nrep];
  Ncov = N-Nrep;

  //  cholesky X Vinv X
  F77_CALL(dpotrf)("L",
      &Ncov, xvx,
      &N, // Ncov by Ncov submatrix of N by N matrix
      &infoCholXX FCONE);

  // determinant for reml
  *detXVinvX  = 0.0;
  for(D=0;D<Ncov;++D){
      *detXVinvX  += log(xvx[D*N+D]);
  }
  *detXVinvX *= 2;

  // invert X Vinv X
  F77_NAME(dpotri)("L",
      &Ncov,
      xvx,
      &N,
      &infoInvXX
      FCONE);

  // put beta hat in last row
  // (first Nrep rows still have LyLx)
  // C= A B, A=xvx, B=LxLy
  F77_NAME(dsymm)(
      "L", "L", // A on left, A in lower
      &Ncov, &Nrep, // C has Nrep rows, Ncov columns
      &oneD,
      xvx, &N,// XVinvX^(-1)
      &copyLx[Nrep],&N, // Lx
      &zeroD,
      &YXVinvYX[Nrep], &N //betahat, ldc
      FCONE FCONE
  );

  // sum of squares
  //  blasBeta  C     alpha A     B
  //     (1)  LyLy + (-1)  LxLy beta
  F77_NAME(dgemm)(
      //	op(A), op(B),
      "T", "N",
      // nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
      &Nrep, &Nrep, &Ncov,
      // alpha
      &moneD,
      // A, lda(A), betahat
      &copyLx[Nrep], &N,
      // B, lda(B), LxLy
      &YXVinvYX[Nrep], &N, // betaHat
      // blasBeta
      &oneD,
      // C, nrow(c)
      YXVinvYX, &N
      FCONE FCONE);

}




double logLoneNuggetMoreArguments(
    double xisqTausq,
    double *DYXVYX,
    double *DbetaHat,
    double *determinant,
    double *determinantForReml,
    double *m2logL,
    double *m2logReL
){
  // Q, L, Lx, copyLx, c, obsCovRot, YXYX, YrepAdd are global
  // and int Nxy, Nobs, Nxysq, Nrep;
  // as are YwkL, EwkL, YwkD, EwkD; // workspaces

  double minusXisqTausq, oneD=1.0;
  double varHatMl, varHatReml;
  double result, *tempPointer;
  int oneI=1, D;

  double xisqTausqTwo[2] = {xisqTausq, 0.0}; // because cholmod wants vector of length 2

  /*
   * V = xisq Q^(-1) + tausq I
   *   =  xisq ( Q^(-1) + (tausq/xisq) I)
   *   =  tausq ( (xisq/tausq) Q^(-1) + I)
   *  = tausq Q^(-1) ( (xisq/tausq) I + Q)
   *  log | V | = N log(tau^2) − log |Q| +  log |(xisq/tausq) I + Q|
   *  woodbury
   *  \left(A+UCV \right)^{-1} =
   *    A^{-1} -
   *    A^{-1}U \left(C^{-1}+VA^{-1}U \right)^{-1} VA^{-1}
   *  C =  Q^{-1},  U = V = I, A = (tausq/xisq) I
   *  V^{-1} = (xisq/tausq) I -
   *    (xisq/tausq)^2 ( Q + (xisq/tausq) I )^{-1}
   * L Lt = Q + (xi^2/tau^2) I
   * V^(-1) =  (xisq/tausq) (I - (xisq/tausq) (L Lt)^(-1) )
   */
  // V^(-1) = (1/tau^2) * Q (L Lt)^(-1)
  // log | V | =N log(tau^2) − log |Q| + 2 log |L|


  // factorize beta*I+A or beta*I+A’*A
  M_cholmod_factorize_p(
      Q,
      xisqTausqTwo, // beta
      (int*)NULL, 0 /*fsize*/,
      L, // resulting factorization
      &c // common
  );

  *determinant = M_chm_factor_ldetL2(L) - *determinant;

  //Lx = L^(-1) obsCov
  M_cholmod_solve2(
      CHOLMOD_L,
      L,
      obsCovRot,
      &Lx,
      &YwkL, &EwkL,
      &c);

  // DYXVYX =  t(YX) (YX)
  // copy of YXYX
  F77_CALL(dcopy)(
      &Nxysq,
      YXYXglobal, // in
      &oneI,
      DYXVYX, // out
      &oneI);

  minusXisqTausq = - xisqTausq;//(xisqTausq*xisqTausq);

  // cross product
  // want   YXVYX -
  //   xisqTausq t(Lx) Lx

  // DYXVYX = - xisqTausq * t(Lx) Lx
  // C := alpha*op( A )*op( B ) + beta*C,
  F77_NAME(dgemm)(
      //	op(A), op(B),
      "T", "N",
      // nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
      &Nxy, &Nxy, &Nobs,
      // alpha
      &minusXisqTausq,
      // A, nrow(A)
      Lx->x, &Nobs,
      // B, nrow(B)
      Lx->x, &Nobs,
      // beta
      &oneD,
      // C, nrow(c)
      DYXVYX, &Nxy
      FCONE FCONE);

  // betaHat =  t(YX) Q^(-1) (YX)
  // copy YXVYX into beta hat
  F77_CALL(dcopy)(
      &Nxysq,
      DYXVYX, // in
      &oneI,
      DbetaHat, // out
      &oneI);

  // compute beta hat and t(resid) %*% Vinv %*% resid
  ssqFromXprod(
      DbetaHat,
      &determinantForReml[0],
      Nxy, Nrep, copyLx);

  for(D=0;D<Nrep;++D){

      determinant[D] = determinant[0];

      determinantForReml[D] = determinantForReml[0];

      // mlrepAdd
      varHatMl = DbetaHat[D*Nxy+D]/Nobs;
      m2logL[D] = Nobs*(onePlusLogTwoPi +log(varHatMl) ) +
          *determinant - YrepAdd[D];

      // reml
      varHatReml = DbetaHat[D*Nxy+D]/(Nobs-Ncov);
      m2logReL[D] = (Nobs-Ncov)*(onePlusLogTwoPi +log(varHatMl) ) +
          *determinant - *determinantForReml - YrepAdd[D];
  }

  if(Ltype){
      tempPointer = m2logReL;
  } else {
      tempPointer = m2logL;
  }

  // find minimum element
  // doesnt yet work
  // R_max_col(tempPointer,&oneI, &Nrep,&D,&oneI);
  D=0;

  result = tempPointer[D];
  return result;

}

/*
 * logL given xisqTausq
// needs global variables
// Q, L, c, detTwo
// obsCovRot, Lx, YwkL, EwkL, DLx, YwkD, EwkD
// YXVYXglobal, YXYX, Nxy, Nobs,
 */

double logLoneNugget(double xisqTausq, void *nothing){

  return
      logLoneNuggetMoreArguments(
          xisqTausq,
          &YXVYXglobal[DxisqTausq*Nxysq],
          &betaHatGlobal[DxisqTausq*Nxysq],
          &YXVYXglobal[Nxysq*NxisqTausq + DxisqTausq*Nrep],
          &YXVYXglobal[Nxysq*NxisqTausq + Nrep*NxisqTausq+ DxisqTausq*Nrep],
          &YXVYXglobal[Nxysq*NxisqTausq + 2*Nrep*NxisqTausq+ DxisqTausq*Nrep],
          &YXVYXglobal[Nxysq*NxisqTausq + 3*Nrep*NxisqTausq + DxisqTausq*Nrep]
      );

}

// logL for calling form Brent_fmin
// this needs YXVYXstart as well as the previous
// function's global arguments
double logLoneLogNugget(double logXisqTausq, void* nothing){

  double result;

  if(DxisqTausq>=NxisqTausq){
      DxisqTausq=1;
  }

  result = logLoneNugget(exp(logXisqTausq), nothing);
  DxisqTausq++;

  return(result);
}

/*
 * callable function from R
 */
SEXP gmrfLik(
    SEXP QR,
    SEXP obsCovR,
    SEXP xisqTausq, // first element must be zero
    SEXP reml,
    SEXP YrepAddR,// Jacobian, length used for Nrep
    SEXP optParam  // set first element -999 for verbose
){

  int Drep, dooptim, verbose=0; // length(xisqTausq)
  int oneI=1;
  double	oneD[2]={1.0,0.0}, zeroD[2]={0.0, 0.0};
  double optTol, optMin, optMax; // default interval for optimizer
  int NxisqMax = 100, Nxyvarmat, twoNxyvarmat; // number of xisqTausq's to retain when optimizing
  double *YXVYX, *YXYX, *determinant, *determinantForReml;
  double *betaHat;
  double *m2logL, *m2logReL, *varHatMl, *varHatReml, *resultXisqTausq;
  double nothing=0.0;
  SEXP resultR;
  CHM_DN obsCov;

  Ltype=*INTEGER(reml); // set to 1 for reml

  // default interval for optimizer
  optMin = REAL(optParam)[0];
  optMax = REAL(optParam)[1];
  optTol = REAL(optParam)[2];// tolerence for fmin_bren

  Nrep =LENGTH(YrepAddR);
  YrepAdd = REAL(YrepAddR);

  Nobs = INTEGER(GET_DIM(obsCovR))[0];
  Nxy = INTEGER(GET_DIM(obsCovR))[1];
  Ncov = Nxy - Nrep;
  Nxysq = Nxy*Nxy;

  NxisqTausq = LENGTH(xisqTausq);
  // if length zero, do optimization
  dooptim=!NxisqTausq;

  Nxyvarmat = Nxysq*NxisqTausq;
  twoNxyvarmat = 2*Nxyvarmat;

  if(dooptim){
      NxisqTausq = NxisqMax;
  } else {
      if(REAL(xisqTausq)[0]>0.0) {
          warning(
              "first element of xisqTausq is %f, must be zero\n",
              REAL(xisqTausq)[0]);
      }
      if(optMin < -998.9) {
          verbose = 1;
          Rprintf("d Nxixqtauxq %d, Nobs %d, Nrep %d, Ncov %d\n",
                  NxisqTausq, Nobs, Nrep, Ncov);
      }
  }

  resultR = PROTECT(allocVector(REALSXP,
                                twoNxyvarmat + 7*Nrep*NxisqTausq + Nxysq));

  for(Drep = 0; Drep < LENGTH(resultR);++Drep){
      REAL(resultR)[Drep] = NA_REAL;
  }

  /*
   * if NxisqTausq=2, Ny=2, Nx=1, Nxysq=9, Nxyvarmat =18
   *  result has length 2*2*(2+1)*(2+1) +7*2*2= 54+28=64
   *  result[0-17] = t(obscov) Q obscov
   *  result[18-35] = betaHat, varBetaHat
   *  result[39] = determinant
   *  result[43] = detReml
   *  result[47] = m2logL
   *  result[51] = m2logLreml
   *  result[55] = varhatMl
   *  result[59] = varhatReml
   *  result[63] = resultxisqTausq
   *  result[64-71] = t(YX) %*% YX
   */

  YXVYX = &REAL(resultR)[0];
  YXVYXglobal =  YXVYX;

  betaHat = &REAL(resultR)[Nxyvarmat];
  betaHatGlobal = betaHat;

  YXYX = &REAL(resultR)[twoNxyvarmat + 7*Nrep*NxisqTausq];
  YXYXglobal = YXYX;

  // determinant of V
  determinant = &REAL(resultR)[twoNxyvarmat];
  // determinant of X Vinv X
  determinantForReml =
      &REAL(resultR)[twoNxyvarmat + Nrep*NxisqTausq];
  m2logL = &REAL(resultR)[twoNxyvarmat + 2*Nrep*NxisqTausq];
  m2logReL = &REAL(resultR)[twoNxyvarmat + 3*Nrep*NxisqTausq];
  varHatMl = &REAL(resultR)[twoNxyvarmat + 4*Nrep*NxisqTausq];
  varHatReml = &REAL(resultR)[twoNxyvarmat + 5*Nrep*NxisqTausq];
  resultXisqTausq =
      &REAL(resultR)[twoNxyvarmat + 6*Nrep*NxisqTausq];

//  copyLx = (double *) calloc(Nxy*Nrep,sizeof(double));
copyLx = (double *) R_alloc(Nxy*Nrep, sizeof(double));


  Q = AS_CHM_SP(QR);
  obsCov = AS_CHM_DN(obsCovR);
  M_R_cholmod_start(&c);

  // get some stuff ready
  if(verbose) {
      Rprintf("starting\n");
  }

  // YXYX cross product of data
  F77_NAME(dgemm)(
      //	op(A), op(B),
      "T", "N",
      // nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
      &Nxy, &Nxy, &Nobs,
      // alpha
      oneD,
      // A, nrow(A)
      obsCov->x, &Nobs,
      // B, nrow(B)
      obsCov->x, &Nobs,
      // beta
      zeroD,
      // C, nrow(&c)
      YXYX, &Nxy
      FCONE FCONE);


  /* .. now allocate Lx .. */
  Lx = M_cholmod_copy_dense(obsCov,&c);


  // for YX Vinv YX
  // Lx is Q YX
  M_cholmod_sdmult(
      Q,
      0, oneD, zeroD, // transpose, scale, scale
      obsCov,Lx,// in, out
      &c);
  // now compute t(YX) Lx
  if(verbose) {
      Rprintf("put t(obscov) Q obscov in result\n");
  }
  F77_NAME(dgemm)(
      //	op(A), op(B),
      "T", "N",
      // nrows of op(A), ncol ob(B), ncol op(A) = nrow(opB)
      &Nxy, &Nxy, &Nobs,
      // alpha
      oneD,
      // A, nrow(A)
      obsCov->x, &Nobs,
      // B, nrow(B)
      Lx->x, &Nobs,
      // beta
      zeroD,
      // C, nrow(c)
      YXVYX, &Nxy
      FCONE FCONE);
  if(verbose) {
      Rprintf("first cholesky\n");
  }
  // Q = P' L D L' P
  L = M_cholmod_analyze(Q, &c);
  M_cholmod_factorize(Q,L, &c);

  obsCovRot = M_cholmod_solve(CHOLMOD_P, L,obsCov,&c);


  // determinant
  determinant[0] = M_chm_factor_ldetL2(L);
  resultXisqTausq[0]= NA_REAL;

  if(verbose) {
      Rprintf("likelihood without nugget, det %f\n",
              determinant[0]);
  }

  // betaHat =  t(YX) Q^(-1) (YX)
  // copy of YXVYX
  F77_CALL(dcopy)(
      &Nxysq,
      YXVYX, // in
      &oneI,
      betaHat, // out
      &oneI);

  ssqFromXprod(
      betaHat, // N by N
      determinantForReml,
      Nxy, Nrep,
      copyLx);

  for(Drep=0;Drep<Nrep;++Drep){

      determinant[Drep] = determinant[0];
      determinantForReml[Drep] = determinantForReml[0];

      varHatMl[Drep] = betaHat[Drep*Nxy+Drep]/Nobs;
      varHatReml[Drep] = betaHat[Drep*Nxy+Drep]/(Nobs-Ncov);

      // mlrepAdd
      m2logL[Drep] = Nobs*(onePlusLogTwoPi + log(varHatMl[Drep])) -
          *determinant -
          YrepAdd[Drep];

      // reml
      m2logReL[Drep] = (Nobs-Ncov)*(log(varHatReml[Drep]) + onePlusLogTwoPi) -
          *determinant - *determinantForReml -
          YrepAdd[Drep];

  }



  //	Rprintf("done zero %d\n", dooptim);


  if(dooptim){
      //		Rprintf(" opt %f %f %f\n", optMin, optMax, optTol);

      DxisqTausq=1;
      // do optimizer
      Brent_fmin(
          optMin, optMax,
          logLoneLogNugget,
          &nothing, optTol);



      NxisqTausq = DxisqTausq;

  } 	else { // not optimizing

      for(DxisqTausq=1;DxisqTausq < NxisqTausq;++DxisqTausq){
          if(verbose) {
              Rprintf("DxisqTausq %d, xisqTausq %f\n",
                      DxisqTausq, REAL(xisqTausq)[DxisqTausq]);
          }

          determinant[DxisqTausq*Nrep] = determinant[0];

          logLoneNuggetMoreArguments(
              REAL(xisqTausq)[DxisqTausq],
              &YXVYX[DxisqTausq*Nxysq], // data cross prod
              &betaHat[DxisqTausq*Nxysq], // betas
              &determinant[DxisqTausq*Nrep],
              &determinantForReml[DxisqTausq*Nrep],
              &m2logL[DxisqTausq*Nrep],
              &m2logReL[DxisqTausq*Nrep]
          );
      }
  }

  for(DxisqTausq=1;DxisqTausq < NxisqTausq;++DxisqTausq){
      for(Drep=0;Drep<Nrep;++Drep){

          varHatMl[DxisqTausq*Nrep + Drep] =
              betaHat[DxisqTausq*Nxysq+Drep*Nxy+Drep]/Nobs;

          varHatReml[DxisqTausq*Nrep + Drep] =
              betaHat[DxisqTausq*Nxysq+Drep*Nxy+Drep]/(Nobs-Ncov);

          resultXisqTausq[DxisqTausq*Nrep + Drep] =
              REAL(xisqTausq)[DxisqTausq];
      }

  }

//  determinant[DxisqTausq*Nrep+Drep] =
//      determinant[DxisqTausq*Nrep];

//  determinantForReml[DxisqTausq*Nrep+Drep]=
//      determinantForReml[DxisqTausq*Nrep];

//  m2logL[DxisqTausq*Nrep+Drep]  -=
//      determinant[0];

//  m2logReL[DxisqTausq*Nrep+Drep] -=
//      determinant[0];



  M_cholmod_free_factor(&L, &c);
  M_cholmod_free_dense(&obsCovRot, &c);

  M_cholmod_free_dense(&Lx, &c);


  // don't free Q because it's from an R object
  //	M_cholmod_free_sparse(&Q, &c);

  // don't free obsCov because it's from an R object
  //	M_cholmod_free_dense(&obsCov, &c);

//  free(copyLx);
  //	free(YXYX);
  M_cholmod_free_dense(&YwkL, &c);
  M_cholmod_free_dense(&YwkD, &c);
  M_cholmod_free_dense(&EwkL, &c);
  M_cholmod_free_dense(&EwkD, &c);

  M_cholmod_finish(&c);

  UNPROTECT(1);
  return resultR;
}


