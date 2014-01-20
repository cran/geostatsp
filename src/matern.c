#include<R.h>
#include<Rmath.h>

void maternArasterBpoints(double *Axmin, double *Axres, int *AxN,
		double *Aymax, double *Ayres, int *AyN,
		double *Bx, double *By, int *BN,
		double *result,
		double  *range, double*shape, double *variance,
		double *anisoRatio, double *anisoAngleRadians) {

int DB, DAx, DAy, AyN2, AxN2, BN2;
int Dindex,Ncell;
double distCellRight[2], distCellDown[2], distTopLeft[2], distRowHead[2];
double distTopLeftR[2], distHere[2];
double costheta, sintheta, anisoRatioSq;
double xscale, varscale,  thisx;
int nb, Nzeros;
double *bk, alpha,truncate;


AyN2 = *AyN;
AxN2 = *AxN;
Ncell = AyN2*AxN2;
BN2 = *BN;
*Axmin += *Axres/2; // add half a cells size, assign each cell
*Aymax -= *Ayres/2; // it's centroid rather than it's top left corner


costheta = cos(*anisoAngleRadians);
sintheta = sin(*anisoAngleRadians);
anisoRatioSq = (*anisoRatio)*(*anisoRatio);

distCellRight[0] = costheta *(*Axres);
distCellRight[1] = sintheta * (*Axres);

distCellDown[0] =  sintheta * (*Ayres);
distCellDown[1] =  - costheta * (*Ayres);


xscale = sqrt(8 * (*shape)) / *range;
varscale =  log(*variance)  - lgammafn(*shape ) -  (*shape -1)*M_LN2;

truncate = *variance*1e-06; // count a zero if var < truncate

alpha = *shape;
// code stolen from R's src/nmath/bessel_k.c
	nb = 1+ (int)floor(alpha);/* nb-1 <= |alpha| < nb */
	bk = (double *) calloc(nb, sizeof(double));


Nzeros=0;
//#pragma omp parallel for private(distTopLeft,distTopLeftR,distRowHead,DAy,distHere,thisx,Dindex)
for(DB=0;DB<BN2;++DB){ // loop through points
	Dindex = DB*Ncell;
	distTopLeft[0]= (Bx[DB]-*Axmin); // distance from point DB to
	distTopLeft[1]= (By[DB]-*Aymax); //   top left corner of raster

	distTopLeftR[0]= costheta * distTopLeft[0] - sintheta * distTopLeft[1];
	distTopLeftR[1] = sintheta *distTopLeft[0] + costheta * distTopLeft[1];

	distRowHead[0] = distTopLeftR[0]; // distance to leftmost cell of row DAy
	distRowHead[1] = distTopLeftR[1];

	for(DAy=0;DAy<AyN2;++DAy){ // loop through y of raster

		distHere[0] = distRowHead[0]; // dist to cell DAx Day
		distHere[1] = distRowHead[1];
		for(DAx=0;DAx<AxN2;++DAx){ // loop through x of raster

			thisx =  sqrt(distHere[0]*distHere[0] +
	      			distHere[1]*distHere[1]/anisoRatioSq)*xscale;

// if thiex is nan assume it's infinity
			if(isnan(thisx)) {
				if(isinf(xscale)) {
// range is probably zero.
					// if distance is zero set result to variance
					if(distHere[0]*distHere[0] +
			      			distHere[1]*distHere[1] < truncate){
						result[Dindex]= *variance;
					}
				} else {
// range is finite, distance must be zero
					result[Dindex] = 0;
				}
			} else {
			result[Dindex] = exp(varscale + *shape * log(thisx) )*
    				bessel_k_ex(thisx, alpha, 1.0, bk);
			}

			if(isnan(result[Dindex]))  {
				// assume distance is very small
				if(thisx < 1) {
					result[Dindex]= *variance;
				} else {
					result[Dindex]= 0;
				}
			}


    		if(result[Dindex]  <  truncate) ++Nzeros;

			++Dindex;
			distHere[0] -= distCellRight[0];
			distHere[1] -= distCellRight[1];
		}

		distRowHead[0] -=distCellDown[0];
		distRowHead[1] -=distCellDown[1];

	}

}
*BN = Nzeros;
*range = xscale;
*shape=varscale;
*anisoRatio = anisoRatioSq;
free(bk);

}

void maternAniso(double *x, double *y, int *N,
		double *result,
		double  *range, double*shape, double *variance,
		double *anisoRatio, double *anisoAngleRadians) {

	int Drow, Dcol, Nm1, Dcolp1, N2;
	int Dindex;

	double xscale, varscale,  thisx;
	double anisoRatioSq, dist[2], distRotate[2], costheta, sintheta;

    int nb, Nzeros;
    double *bk, alpha,truncate;

    costheta = cos(*anisoAngleRadians);
    sintheta = sin(*anisoAngleRadians);


    anisoRatioSq = (*anisoRatio)*(*anisoRatio);

	xscale = sqrt(8 * (*shape)) / *range;
	varscale =  log(*variance)  - lgammafn(*shape ) -  (*shape -1)*M_LN2;

    truncate = *variance*1e-06; // count a zero if var < truncate

	alpha = *shape;
	// code stolen from R's src/nmath/bessel_k.c
		nb = 1+ (int)floor(alpha);/* nb-1 <= |alpha| < nb */
		bk = (double *) calloc(nb, sizeof(double));


    Nm1 = *N-1;
    N2 = *N;
    Dindex = 0;
    Nzeros = 0;

    for(Dcol=0;Dcol < Nm1;++Dcol) {
    	Dcolp1 = Dcol + 1;
    	Dindex += Dcol;

    	for(Drow=Dcolp1;Drow < N2; ++Drow) {
    		Dindex += 1;

    		dist[0] = x[Dcol] - x[Drow];
    		dist[1] = y[Dcol] - y[Drow];

// rotate anticlockwise by aniso.angle.radians
    	// distRotate =  ( cos(theta)  -sin(theta)  )  dist
    	//               ( sin(theta)   cos(theta)  )
    		distRotate[0] = costheta *dist[0] - sintheta * dist[1];
    		distRotate[1] = sintheta *dist[0] + costheta * dist[1];

    		thisx =  sqrt(distRotate[0]*distRotate[0] +
      			distRotate[1]*distRotate[1]/anisoRatioSq)*xscale;

			if(isnan(thisx)) {
				if(isinf(xscale)) {
	// range is probably zero.
					// if distance is zero set result to variance
					if(distRotate[0]*distRotate[0] +
							distRotate[1]*distRotate[1] < truncate){
						result[Dindex]= *variance;
					}
				} else {
	// range is finite, distance must be zero
						result[Dindex] = 0;
				}
			} else { // thisx not nan
    		result[Dindex] = exp(varscale + *shape * log(thisx) )*
    				bessel_k_ex(thisx, alpha, 1.0, bk);
			}

			if(isnan(result[Dindex]))  {
					// assume distance is very small
					if(thisx < 1) {
						result[Dindex]= *variance;
					} else {
						result[Dindex]= 0;
					}
				}



    		if(result[Dindex]  <  truncate) ++Nzeros;

    	}
		Dindex += 1;
    }
	*N = Nzeros;
    free(bk);
}


void matern(double *distance, int *N,
		double *range, double *shape, double *variance) {

	int D, N2;
	double xscale, varscale,  thisx;

    int nb,  Nzeros;
    double *bk, alpha,truncate;

    truncate = *variance*1e-06; // count a zero if var < truncate
    Nzeros = 0;

	alpha = *shape;

// code stolen from R's src/nmath/bessel_k.c
	nb = 1+ (int)floor(alpha);/* nb-1 <= |alpha| < nb */

	bk = (double *) calloc(nb, sizeof(double));

	N2 = *N;// for some reason need D to be int, not long.
// evaluate the matern!

	/*
	xscale = abs(x)*(sqrt(8*param["shape"])/ param["range"])
	result = ( param["variance"]/(gamma(param["shape"])* 2^(param["shape"]-1)  ) ) *
			( xscale^param["shape"] *
				besselK(xscale , param["shape"]) )
*/

	xscale = sqrt(8 * (*shape)) / *range;
	varscale =  log(*variance)  - lgammafn(*shape ) -  (*shape -1)*M_LN2;
// #ifdef SUPPORT_OPENMP
// #pragma omp parallel for private(thisx)
// #endif
	for(D=0; D < N2; D++) {
		thisx = fabs(distance[D])*xscale;

		if(isnan(thisx)) {
//			warning("%f %f", thisx, xscale);
			if(isinf(xscale)) {
// range is probably zero.
				// if distance is zero set result to variance
				if(fabs(distance[D]) < truncate){
					distance[D]= *variance;
				}
			} else {
// range is finite, distance must be zero
				distance[D] = 0;
			}
		} else { // thisx not nan
			distance[D] = exp(varscale + *shape * log(thisx) )*
					bessel_k_ex(thisx, alpha, 1.0, bk);
		}
		if(isnan(distance[D])) {
			// assume distance is very small
			if(thisx < 1) {
				distance[D]= *variance;
			} else {
				distance[D]= 0;
			}
		}


		if(distance[D] <  truncate) ++Nzeros;
	}
	*range = xscale;
	*shape=varscale;
	*N = Nzeros;

    free(bk);

}


