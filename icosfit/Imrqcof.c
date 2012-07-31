#include "ICOSfit.h"
#include <stdio.h>
#include <setjmp.h>
#include "nortlib.h" // for debugging only
#include "nrutil.h"

extern jmp_buf Fit_buf;

#if ! defined( __QNXNTO__ )
  int isnanf( float f ) {
    union {
      float f;
      long int li;
    } u;
    u.f = f;
    return ( ( u.li & 0x7F800000 ) == 0x7F800000 );
  }
#endif

// Here I'm trying to override the object members alpha and beta
// with method arguments.
void fitdata::mrqcof( float *av, float **alpha, float *beta ) {
  int i,j,k,l,m;
  float wt,sig2i,dy;

  for (j=1;j<=mfit;j++) {
    for (k=1;k<=j;k++) alpha[j][k]=0.0;
    beta[j]=0.0;
  }
  chisq=0.0;
  for (i=1;i<=npts;i++) {
    func->evaluate( x[i], av );
    if ( isnanf( func->value ) ) {
      nl_error( 2, "evaluate x[%d]=%lf returned NaN", i, x[i] );
      longjmp(Fit_buf, 1);
      // throw 1;
    }
    dy = y[i] - func->value;
    for (j=0; j < ma; j++ ) {
      if ( isnanf( func->params[j].dyda ) ) {
        nl_error( 2, "evaluate x[%d]=%lf dy/da[%d] returned NaN", i, x[i], j );
        longjmp(Fit_buf, 1);
        // throw 1;
      }
      dyda[j+1] = func->params[j].dyda;
    }
    sig2i = 1.0 / (sig[i]*sig[i]);
    for (j=0,l=1;l<=ma;l++) {
      if (ia[l]) {
        wt=dyda[l]*sig2i;
        for (j++,k=0,m=1;m<=l;m++)
          if (ia[m]) alpha[j][++k] += wt*dyda[m];
        beta[j] += dy*wt;
      }
    }
    chisq += dy*dy*sig2i;
  }
  for (j=2;j<=mfit;j++)
    for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
}
