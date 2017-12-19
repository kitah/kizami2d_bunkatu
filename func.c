#include <stdlib.h>
#include <math.h>
#include "func.h"

//#define EPS  1E-8
//#define SGMC  0.87
 
// PSF
int hsamp_gauss(int nf, double sgm, double **h){
  float sum;
  int n, nh = (nf - 1) / 2, m;

  /***** ガウス関数をサンプリングした h *****/  
  if (fabs(sgm) < EPS) {
    for(n = -nh ; n <= nh ; n++) {
      for(m = -nh ; m <= nh ; m++) {
	h[n][m] = 0.0;
      }
    }
    h[0][0] = 1.0;
  } else {
    sum = 0.0;
    for(n = -nh ; n <= nh ; n++) {
      for(m = -nh ; m <= nh ; m++) {
	h[n][m] = exp(-(double)(n * n + m * m) / 2.0  / (sgm * sgm + SGMC * SGMC) );
	sum += h[n][m];
      }
    }
    for(n = -nh ; n <= nh ; n++) {
      for(m = -nh ; m <= nh ; m++) {
	h[n][m] /= sum;
      }
    }
  }

  return 0;
}

/*--------------------------------------------------------------------------*/  
double gauss(void){
  double rnd(void);

  return rnd() + rnd() + rnd() + rnd() + rnd() + rnd() +
    rnd() + rnd() + rnd() + rnd() + rnd() + rnd() - 6.0;
}

/*--------------------------------------------------------------------------*/  
double rnd(void){
  return (double)random() / (double) RAND_MAX;
}


