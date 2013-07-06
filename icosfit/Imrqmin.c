#include "config.h"
#include <stdio.h>
#include "ICOSfit.h"
#include "nrutil.h"

int fitdata::mrqmin() {
  void covsrt(float **covar, int ma, int ia[], int mfit);
  void gaussj(float **a, int n, float **b, int m);
  int j,k,l;

  if (alamda < 0.0) {
    for (mfit=0,j=1;j<=ma;j++)
      if (ia[j]) mfit++;
    if (mfit > mf_size ) {
      if ( mf_size != 0 )
        free_matrix(oneda,1,mf_size,1,1);
      mf_size = mfit;
      oneda=matrix(1,mfit,1,1);
    }
    alamda=0.001;
    mrqcof( a, alpha, beta );
    for (j=1;j<=ma;j++) atry[j]=a[j];
    ochisq = chisq;
  }
  // ochisq = chisq;
  for (j=1;j<=mfit;j++) {
    for (k=1;k<=mfit;k++) {covar[j][k]=alpha[j][k]; }
    covar[j][j]=alpha[j][j]*(1.0+(alamda));
    oneda[j][1]=beta[j];
  }
  if (verbose & 32 ) {
    fprintf(stderr, "In fitdata::mrqmin\n" );
    print_matrix( covar, "covar", mfit, mfit );
    print_matrix( oneda, "oneda", mfit, 1 );
  }
  gaussj(covar,mfit,oneda,1);
  if (verbose & 32 ) {
    fprintf(stderr, "In fitdata::mrqmin after gaussj:\n" );
    print_matrix( covar, "covar", mfit, mfit );
    print_matrix( oneda, "oneda", mfit, 1 );
  }
  for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
  if (alamda == 0.0) {
    covsrt(covar,ma,ia,mfit);
    covsrt(alpha,ma,ia,mfit);
    return 0;
  }
  for (j=0,l=1;l<=ma;l++)
    if (ia[l]) {
      atry[l]=a[l]+da[++j];
    }
  if ( adjust_params( atry) ) return 1;
  mrqcof( atry, covar, da );
  // if (chisq == ochisq ) ochisq = chisq/.9995;
  if (chisq <= ochisq) {
    alamda *= 0.1;
    // ochisq=chisq;
    for (j=1;j<=mfit;j++) {
      for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
      beta[j]=da[j];
    }
    for (l=1;l<=ma;l++) a[l]=atry[l];
  } else {
    alamda *= 10.0;
    // chisq=ochisq;
  }
  return 0;
}
