#ifndef PTREAD_H_INCLUDED
#define PTREAD_H_INCLUDED
#include <stdio.h>
#include "mlf.h"
#include "f_vector.h"

class PTfile {
  public:
    FILE *fp;
    double time;
    double P, T;
    double cal_flow, inlet_flow;
    unsigned long ScanNum;
    int RORIS, RateS;
    unsigned long next_ScanNum;

    int readline();
    void backup();
    void calc_wndata();
    PTfile( const char *fname );
  private:
    int format;
    int n_vars;
    long int last_file_pos;
    double Etln_params[8];
};

class ICOSfile {
  public:
    ICOSfile(const char *fbase, const char *obase, int bin );
    int read(unsigned long int fileno);
    FILE *writefp();
    int wn_sample( float wn );

    mlf_def_t *mlf;
    mlf_def_t *omlf;
    FILE *ofp;
    f_vector *sdata;
    f_vector *fdata;
    static f_vector *bdata;
    static f_vector *wndata;
    static float nu_F0;
    static int dFN;
    int binary;
    static const int mindatasize;
};

#endif
