#ifndef ICOSfit_H_INCLUDED
#define ICOSfit_H_INCLUDED

#include "config.h"
#include "nortlib.h"
#include "funceval.h"
#include "ptread.h"

#define ICOSFIT_VERSION "2.12"
#define ICOSFIT_VERSION_DATE "12/23/2012"

class fitdata {
  public:
	int BaseStart;
	int BaseEnd;
	int SignalStart;
	int SignalEnd;
	int Start;
	int End;
	int FitBaseline;
	int npts;
	int npts_vec;
	int ma;
	int *ia;
	float *a, *x, *y, *sig;
	float *a_save;
	float **covar, **alpha;
	float chisq, ochisq;
	float alamda;
	func_evaluator *func;
	func_evaluator *base;
	func_abs *absorb;
	PTfile *PTf;
	ICOSfile *IFile;
	int verbose;
	FILE *vfp;
	mlf_def_t *vmlf;
	static const int n_input_params;
	static const int ScanNum_col;
	static const int dFN_col;
    fitdata( PTfile *ptf, ICOSfile *IF,
       func_evaluator *f, func_evaluator *baseline, func_abs *abs );
    void handle_restart( const char *ofname );
    int fit();
    void write();
    void lwrite( FILE *ofp, FILE *vofp, int fileno );
    int mrqmin( );
    void mrqcof( float *av, float **alpha, float *beta );
    int adjust_params( float *av );
  private:
	int mfit, mf_size;
	float *atry,*beta,*da,**oneda, *dyda;
};

// Located in build.cc
extern fitdata *build_func();

// Located in loadlin.cc
extern func_abs *load_absorb( char *filename );

// Located in fitfunc.cc
void print_matrix( float **mat, const char *name, int nrow, int ncol );
void print_vector( float *vec, const char *name, int ncol );

#endif
