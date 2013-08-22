#ifndef FUNCEVAL_H_INCLUDED
#define FUNCEVAL_H_INCLUDED

#include <stdio.h>
#include <vector>

struct parameter {
  int index;
  int *ia;
  float init;
  float dyda;
  float low;
  float high;
  float prev;
};

class func_evaluator {
  public:
    func_evaluator( const char *name, int n_params = 0 ); 
    virtual void evaluate( float x, float *a );
    void init( float *a, int *ia );
    void init(float *a, int p1);
    virtual void init(float *a);
    virtual int adjust_params( float alamda, float P, float T, float *a );
    void append_func( func_evaluator *newfunc );
    inline void link_param( int src, func_evaluator *child, int dest ) {
      child->params[dest].index = params[src].index;
      child->params[dest].ia = params[src].ia;
    }
    inline float get_param( float *a, int idx ) {
      return a[params[idx].index];
    }
    inline float set_param( float *a, int idx, float value ) {
      a[params[idx].index] = value;
      return value;
    }
    inline void fix_param( int idx ) { *params[idx].ia = 0; }
    inline void float_param( int idx ) { *params[idx].ia = 1; }
    inline int param_fixed( int idx ) { return *params[idx].ia == 0; }
    void clamp_param_high( float *a, int idx );
    void clamp_param_low( float *a, int idx );
    void clamp_param_highlow( float *a, int idx );
    virtual int line_check(int include, float& start, float& end,
                            float P, float T, float *a);
    virtual int skew_samples();
    virtual void dump_params(float *a, int indent);
    void print_indent( FILE *fp, int indent );

    int n_params;
    parameter *params;
    float value;
    const char *name;
    func_evaluator *parent;
    func_evaluator *first;
    func_evaluator *last;
    func_evaluator *next;
};

// func_aggregate and its derived classes assume all the
// sub-functions have independent parameters. Hence total
// number of parameters is just the sum of the number of
// parameters of the children.
class func_aggregate : public func_evaluator {
  public:
    inline func_aggregate() : func_evaluator("aggregate", 0) {}
    inline func_aggregate(const char*sname) : func_evaluator(sname,0) {}
    void append_func( func_evaluator * );
};
class func_product : public func_aggregate {
  public:
    inline func_product() : func_aggregate("product") {}
    void evaluate(float x, float *a);
};
class func_sum : public func_aggregate {
  public:
    inline func_sum() : func_aggregate("sum") {}
    inline func_sum(char *name) : func_aggregate(name) {}
    void evaluate(float x, float *a);
    void evaluate(float x, float *a, int i);
};

class QTdata {
  public:
    QTdata(int isotopologue);
    ~QTdata();
    double evaluate(double T);
    int isotop;
  private:
    int Tmin, Tmax, dT;
    static const int Tref = 296;
    static const int QTBUFSIZE = 80;
    std::vector<double> QT;
};

//--------------------------------------------------------
// func_line objects share location and width parameters.
// Derived classes include gaussian, lorentzian and voigt,
// which each add additional parameters.
//
// nu - an approximation of the line wavenumber
// nu0 - static offset to extend precision of wavenumbers
//       Set when the first line is defined.
// nu1 - nu-nu0 (or a more precise version thereof)
// nu_P - nu1+delta*P/760. : pressure-shifted wavenumber less nu0
//--------------------------------------------------------
class func_line : public func_evaluator {
  public:
    func_line( const char *name, int np, int mol, int iso,
      double nu, double S, double Gair, double E, double n,
      double delta, unsigned int ipos, double threshold, int fix_w,
      int fix_fp );
    ~func_line();
    int adjust_params( float alamda, float P, float T, float *a );
    static const int l_idx, w_idx, n_idx;
    static int n_lines;
    int line_number;
    int fixed; // 0 = free, 1 = fixed
    int fix_finepos; // 0 = free-ish, 1 = fixed
    int fix_width; // 0 = free-ish, 1 = fixed
    inline func_line *lnext() { return (func_line *)next; }
    void init(float *a);
    virtual float line_start(float *a);
    virtual float line_end(float *a);
    virtual void line_fix();
    virtual void line_float();
    int line_check(int include, float& start, float& end, float P, float T, float *a);
    void print_config(FILE *fp);
    virtual void print_intermediates(FILE *fp);
    //--------------------------------------------------
    int isotopomer;
    double nu;
    float nu1, S, G_air, E, n_air, delta;
    unsigned int ipos; // was loc...
    float S_thresh;
    float molwt;
    float Ks, nu_P, Corr_Tref;
    float prev_numdens;
    float prev_ged;
    int rolledback;
    QTdata *QT;
    static double nu0;
    static const double DRTPI; // 1/SQRT(pi)
    static const double Tref; // 296. K
    static const double C2; // 1.4388 cm K second radiation constant hc/k
};

// func_abs has 1 common parameter, which is nu_F0
// In Release 2.2, we inherit nu_F0 from func_skew.
// Each line then has it's own dnu that func_abs owns
// and its own 4 parameters.
class func_abs : public func_evaluator {
  public:
    inline func_abs() : func_evaluator("abs") { n_params = 1; }
    void append_func( func_line *newfunc );
    void init(float *a);
    void evaluate(float x, float *a);
    inline func_line *lfirst() { return (func_line *)first; }
    int adjust_params(float alamda, float P, float T, float *a);
    void print_config(FILE *fp);
    void print_intermediates(FILE *fp);
    void fix_linepos(int linenum);
    void float_linepos(int linenum);
    void dump_params(float *a, int indent);
};
typedef func_abs *func_abs_p;
typedef func_line *func_line_p;
inline func_abs_p new_func_abs() { return new func_abs(); }
inline func_abs_p abs_append( func_abs_p abs, func_line_p line ) {
  abs->append_func( line );
  return abs;
}

class gaussian : public func_line {
  public:
    inline gaussian( int mol, int iso,
          double nu_in, double S_in, double G_air_in, double E_in,
          double n_in, double delta_in, int ipos_in, double threshold, int fix_w ) :
       func_line( "gaussian", 3, mol, iso, nu_in, S_in, G_air_in, E_in,
           n_in, delta_in, ipos_in, threshold, fix_w, 1 ) {};
    void evaluate(float x, float *a);
};
class lorentzian : public func_line {
  public:
    inline lorentzian( int mol, int iso,
          double nu_in, double S_in, double G_air_in, double E_in,
          double n_in, double delta_in, int ipos_in, double threshold, int fix_w ) :
       func_line( "lorentzian", 3, mol, iso, nu_in, S_in, G_air_in, E_in,
           n_in, delta_in, ipos_in, threshold, fix_w, 1 ) {};
    void evaluate(float x, float *a);
};
// Based on R.J.Wells' functions
class voigt : public func_line {
  public:
    voigt( int mol, int iso,
          double nu_in, double S_in, double G_air_in, double E_in,
          double n_in, double delta_in, int ipos_in, double threshold,
          int fix_dw, int fix_lw, int fix_fp );
    void init(float *a);
    void evaluate(float x, float *a);
    float line_width(float*a);
    float line_start(float*a);
    float line_end(float*a);
    static const int gl_idx; // 3
    int adjust_params( float alamda, float P, float T, float *a );
    void line_fix();
    void line_float();
    void dump_params(float *a, int indent);
    void print_intermediates(FILE *fp);

  private:
    float prev_gl;
    float prev_y;
    float X, Y, K;
    float YQ, XLIMA, XLIMB, XLIMC, XLIM4;
    int RGB, RGC, RGD;
    float A0, B1, C0, C2, D0, D1, D2, E0, E2, E4, F1, F3, F5;
    float G0, G2, G4, G6, H0, H2, H4, H6, P0, P2, P4, P6, P8;
    float Q1, Q3, Q5, Q7, R0, R2, W0, W2, W4, Z0, Z2, Z4, Z6, Z8;
    int fix_lwidth;
};
inline func_line_p new_voigt( int mol, int iso,
    double nu, double S, double Gair, double E, double n,
    double delta, int ipos, double threshold,
    int fix_dw, int fix_lw, int fix_fp ) {
  return new voigt( mol, iso, nu, S, Gair, E, n, delta, ipos, threshold,
                    fix_dw, fix_lw, fix_fp );
}

class func_quad : public func_evaluator {
  public:
    func_quad( float q, float l, float c );
    void evaluate(float x, float *a);
    static const int q_idx, l_idx, c_idx; // 0, 1, 2
};


// baseline functions. func_base is a virtual base class
// (perhaps not technically?)
class func_base : public func_evaluator {
  public:
    inline func_base() : func_evaluator("base",0) {}
    inline func_base(const char *name) : func_evaluator(name,0) {}
    int uses_nu_F0;
};

// This is the old func_base for SVDs as a function of x
// The file format is ICOS standard binary with the first
// element of each column being initial parameter settings
class func_base_svdx : public func_base {
  public:
    func_base_svdx( const char *filename );
    void evaluate( float x, float *a );
  private:
    float **baseline;
    int n_pts;
};

// This is a function supporting vectors as a function of
// nu as well as polynomials as a function of x or nu
// The file format is:
// i*4  0  To distinguish from func_base_svdx format
// i*4  1  To identify func_base_ptbnu format
// f*8 polynomial scale factor
// f*8 nu0 smallest value of nu
// f*8 dnu increment between samples
// i*2 n_vectors
// i*2 npts per vector
// i*2 number of polynomial coefficients (degree+1)
// i*2 polynomial of x == 0, nu == 1
// f*4 X n_vectors initial parameter values
// f*4 X n polynomial coefficients initial parameter values
// f*4 X npts X n_vectors: vector data in column order

// Can use global SignalRegion to determine range of x
class func_base_ptbnu : public func_base {
  public:
    func_base_ptbnu( const char *filename );
    void evaluate( float x, float *a );
    void init(float *a);
  private:
    struct {
      double poly_scale;
      double nu0;
      double dnu;
      unsigned short n_vectors;
      unsigned short n_pts;
      unsigned short poly_coeffs;
      unsigned short poly_of_nu;
    } cfg;
    float **vectors;
    float **dvdnu;
    float **polyvecs;
};

class func_base_input : public func_base {
  public:
    func_base_input( func_base *base );
    void evaluate( float x, float *a );
    void init(float *a);
  private:
};

extern func_base *pick_base_type( const char *filename );


//-------------
// func_skew.c
//-------------
class skew_data {
  public:
    skew_data();
    void set_n_params(int n_gp, int n_ap);
    float g;
    float *dg; // allocate n_params
    float *da; // allocate n_params
    float gN;
    float Power;
    int initialized;
};

// func_skew applies the cell skew function of ICOS. Its two
// children define the input power curve (base) and the
// intra-cavity absorption (abs).
// base and abs have mostly independent parameters, but if
// the base function uses nu_F0, that parameter is shared.
// The skew member is used as an M-element circular buffer.
class func_skew : public func_evaluator {
  public:
    func_skew(func_base *base, func_abs *abs);
    void init(float *a);
    void evaluate(float x, float *a);
    int skew_samples();
    void dump_params(float *a, int indent);
  private:
    float N;
    int M;
    float R2, R2N, P_scale;
    skew_data *skew; // We will have M of these
    float prev_x;
    int skewidx;
    int n_base_params, n_abs_params;
    func_base *basep;
    func_abs *absp;
};

// func_noskew calculates absorption for a simple multi-pass
// cell. Its two children define the input power curve (base)
// and the intra-cavity absorption (abs).
// base and abs have mostly independent parameters, but if
// the base function uses nu_F0, that parameter is shared.
class func_noskew : public func_evaluator {
  public:
    func_noskew(func_base *base, func_abs *abs);
    void init(float *a);
    void evaluate(float x, float *a);
    void dump_params(float *a, int indent);
  private:
    float N_Passes;
    int n_base_params, n_abs_params;
    func_base *basep;
    func_abs *absp;
};

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

#if defined( __QNXNTO__ )
  extern "C" {
    extern int isnanf(float) __attribute__((__const__));
  };
#else
  extern int isnanf( float f );
#endif


#endif
