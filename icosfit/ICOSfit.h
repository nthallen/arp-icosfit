#ifndef ICOSfit_H_INCLUDED
#define ICOSfit_H_INCLUDED

#include "config.h"
#include "nortlib.h"
#include "funceval.h"
#include "ptread.h"

#define ICOSFIT_VERSION "3.00alpha" FLOAT_STYLE
#define ICOSFIT_VERSION_DATE "01/03/2018"

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
    ICOS_Float *a, *x, *y, *sig;
    ICOS_Float *a_save;
    ICOS_Float **covar, **alpha;
    ICOS_Float chisq, ochisq;
    ICOS_Float alamda;
    func_evaluator *func;
    func_evaluator *base;
    func_abs *absorb;
    PTfile *PTf;
    ICOSfile *IFile;
    int verbose;
    FILE *vfp;
    mlf_def_t *vmlf;
    static int n_input_params;
    static const int ScanNum_col;
    fitdata( PTfile *ptf, ICOSfile *IF,
       func_evaluator *f, func_evaluator *baseline, func_abs *abs );
    void handle_restart( const char *ofname );
    int fit();
    void write();
    void lwrite( FILE *ofp, FILE *vofp, int fileno );
    int mrqmin( );
    void mrqcof( ICOS_Float *av, ICOS_Float **alpha, ICOS_Float *beta );
    int adjust_params( ICOS_Float *av );
  private:
    int mfit, mf_size;
    ICOS_Float *atry,*beta,*da,**oneda, *dyda;
};

// Located in build.cc
extern fitdata *build_func();

// Located in loadlin.cc
extern func_abs *load_absorb( char *filename );

// Located in fitfunc.cc
void print_matrix( ICOS_Float **mat, const char *name, int nrow, int ncol );
void print_vector( ICOS_Float *vec, const char *name, int ncol );

#endif
