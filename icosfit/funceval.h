#ifndef FUNCEVAL_H_INCLUDED
#define FUNCEVAL_H_INCLUDED

#include <stdio.h>
#include <vector>
#include "config.h"

class paramref {
  public:
    inline paramref(int arg_i, int param_i) {
      arg_num = arg_i;
      param_num = param_i;
    }
    int arg_num;
    int param_num;
};

class parameter {
  public:
    inline parameter(int idx) : index(idx) {}
    int index;
    ICOS_Float dyda;
    std::vector<paramref> refs;
};

class func_evaluator {
  public:
    func_evaluator(const char *name, bool indexed = false, int index = 0); 
    virtual void evaluate( ICOS_Float x, ICOS_Float *a );
    void init( ICOS_Float *a, int *ia );
    virtual void init(ICOS_Float *a);
    void append_func(func_evaluator *newfunc);
    virtual void adopted(func_evaluator *new_parent);
    virtual void fix_float_param(bool float_it);
    virtual bool param_fixed();
    inline bool param_fixed(int i) { return args[i]->param_fixed(); }
    inline void fix_param() { fix_float_param(false); }
    inline void fix_param(int i) { args[i]->fix_param(); }
    inline void float_param() { fix_float_param(true); }
    inline void float_param(int i) { args[i]->float_param(); }
    virtual ICOS_Float get_param(ICOS_Float *a);
    inline ICOS_Float get_param( ICOS_Float *a, int idx ) {
      return args[idx]->get_param(a);
    }
    virtual ICOS_Float set_param(ICOS_Float *a, ICOS_Float value);
    inline ICOS_Float set_param(ICOS_Float *a, int idx, ICOS_Float value) {
      return args[idx]->set_param(a, value);
    }
    virtual int adjust_params(ICOS_Float alamda, ICOS_Float P, ICOS_Float T, ICOS_Float *a);
    virtual int line_check(int include, ICOS_Float& start, ICOS_Float& end,
                            ICOS_Float P, ICOS_Float T, ICOS_Float *a);
    virtual int skew_samples();
    virtual void dump_params(ICOS_Float *a, int indent);
    void print_indent( FILE *fp, int indent );

    // int n_params;
    std::vector<func_evaluator*> args;
    std::vector<parameter> params;
    // parameter *params;
    ICOS_Float value;
    const char *name;
    func_evaluator *parent;
    // func_evaluator *first;
    // func_evaluator *last;
    // func_evaluator *next;
  protected:
    static std::vector<func_evaluator*> evaluation_order;
};

class func_parameter : public func_evaluator {
  public:
    func_parameter(const char *name, ICOS_Float init_value,
                   bool indexed = false, int index = 0);
    void init(ICOS_Float *a);
    void fix_float_param(bool float_it);
    bool param_fixed();
    void evaluate( ICOS_Float x, ICOS_Float *a );
    void dump_params(ICOS_Float *a, int indent);
    static inline void set_ia(int *ia) { func_parameter::ia = ia; }
    
    // Functions we might need in some form, possibly renamed
    // inline void link_param( int src, func_evaluator *child, int dest ) {
      // child->params[dest].index = params[src].index;
      // child->params[dest].ia = params[src].ia;
    // }
    inline ICOS_Float get_param(ICOS_Float *a) { return a[params[0].index]; }
    inline ICOS_Float set_param(ICOS_Float *a, ICOS_Float value) {
      a[params[0].index] = value;
      return value;
    }
    // inline void fix_param( int idx ) { ia[params[idx].index] = 0; }
    // inline void float_param( int idx ) { ia[params[idx].index] = 1; }
    // inline int param_fixed( int idx ) { return ia[params[idx].index] == 0; }
    void clamp_param_high( ICOS_Float *a, int idx );
    void clamp_param_low( ICOS_Float *a, int idx );
    void clamp_param_highlow( ICOS_Float *a, int idx );
    // End of possible functions
    
    int index; ///< Parameter's global index
  protected:
    bool is_fixed();
    ICOS_Float init_val; ///< Initialization value
    // ICOS_Float low;
    // ICOS_Float high;
    // ICOS_Float prev;
  private:
    static int *ia; ///< The 1-based vector of parameter fix/float settings
    static int n_parameters;
};

// func_aggregate and its derived classes assume all the
// sub-functions have independent parameters. Hence total
// number of parameters is just the sum of the number of
// parameters of the children.
// class func_aggregate : public func_evaluator {
  // public:
    // inline func_aggregate(const char*sname = "aggregate") : func_evaluator(sname) {}
    // void append_func( func_evaluator * );
// };
// class func_product : public func_aggregate {
  // public:
    // inline func_product() : func_aggregate("product") {}
    // void evaluate(ICOS_Float x, ICOS_Float *a);
// };
// class func_sum : public func_aggregate {
  // public:
    // inline func_sum() : func_aggregate("sum") {}
    // inline func_sum(char *name) : func_aggregate(name) {}
    // void evaluate(ICOS_Float x, ICOS_Float *a);
    // void evaluate(ICOS_Float x, ICOS_Float *a, int i);
// };

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
    func_line( const char *name, int mol, int iso,
      double nu, double S, double Gair, double E, double n,
      double delta, unsigned int ipos, double threshold, int fix_w,
      int fix_fp );
    ~func_line();
    int adjust_params( ICOS_Float alamda, ICOS_Float P, ICOS_Float T, ICOS_Float *a );
    static const int l_idx, w_idx, n_idx;
    static int n_lines;
    int line_number; ///< 1-based indexing
    int fixed; ///< 0 = free, 1 = fixed
    int fix_finepos; ///< 0 = free-ish, 1 = fixed
    int fix_width; ///< 0 = free-ish, 1 = fixed
    // ### haven't figured out how to deal with this
    // inline func_line *lnext() { return (func_line *)next; }
    void init(ICOS_Float *a);
    virtual ICOS_Float line_start(ICOS_Float *a);
    virtual ICOS_Float line_end(ICOS_Float *a);
    virtual void line_fix();
    virtual void line_float();
    int line_check(int include, ICOS_Float& start, ICOS_Float& end,
            ICOS_Float P, ICOS_Float T, ICOS_Float *a);
    void print_config(FILE *fp);
    virtual void print_intermediates(FILE *fp);
    //--------------------------------------------------
    int isotopomer;
    double nu;
    ICOS_Float nu1, S, G_air, E, n_air, delta;
    unsigned int ipos; // was loc...
    ICOS_Float S_thresh;
    ICOS_Float molwt;
    ICOS_Float Ks, nu_P, Corr_Tref;
    ICOS_Float prev_numdens;
    ICOS_Float prev_ged;
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
// ### Just review this comment:
// ### In Release 3.0, each line has both nu_F0 and dnu.
// ### Lines own dnu, nu_F0 is owned by func_abs, or maybe func_skew
class func_abs : public func_evaluator {
  public:
    func_abs();
    void append_func( func_line *newfunc );
    void init(ICOS_Float *a);
    void evaluate(ICOS_Float x, ICOS_Float *a);
    // inline func_line *lfirst() { return (func_line *)first; }
    int adjust_params(ICOS_Float alamda, ICOS_Float P, ICOS_Float T, ICOS_Float *a);
    void print_config(FILE *fp);
    void print_intermediates(FILE *fp);
    void fix_linepos(int linenum);
    void float_linepos(int linenum);
    void dump_params(ICOS_Float *a, int indent);
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
       func_line( "gaussian", mol, iso, nu_in, S_in, G_air_in, E_in,
           n_in, delta_in, ipos_in, threshold, fix_w, 1 ) {};
    void evaluate(ICOS_Float x, ICOS_Float *a);
};
class lorentzian : public func_line {
  public:
    inline lorentzian( int mol, int iso,
          double nu_in, double S_in, double G_air_in, double E_in,
          double n_in, double delta_in, int ipos_in, double threshold, int fix_w ) :
       func_line( "lorentzian", mol, iso, nu_in, S_in, G_air_in, E_in,
           n_in, delta_in, ipos_in, threshold, fix_w, 1 ) {};
    void evaluate(ICOS_Float x, ICOS_Float *a);
};
// Based on R.J.Wells' functions
class voigt : public func_line {
  public:
    voigt( int mol, int iso,
          double nu_in, double S_in, double G_air_in, double E_in,
          double n_in, double delta_in, int ipos_in, double threshold,
          int fix_dw, int fix_lw, int fix_fp );
    void init(ICOS_Float *a);
    void evaluate(ICOS_Float x, ICOS_Float *a);
    ICOS_Float line_width(ICOS_Float*a);
    ICOS_Float line_start(ICOS_Float*a);
    ICOS_Float line_end(ICOS_Float*a);
    static const int gl_idx; // 3
    int adjust_params( ICOS_Float alamda, ICOS_Float P, ICOS_Float T, ICOS_Float *a );
    void line_fix();
    void line_float();
    void dump_params(ICOS_Float *a, int indent);
    void print_intermediates(FILE *fp);

  private:
    ICOS_Float prev_gl;
    ICOS_Float prev_y;
    ICOS_Float X, Y, K;
    ICOS_Float YQ, XLIMA, XLIMB, XLIMC, XLIM4;
    int RGB, RGC, RGD;
    ICOS_Float A0, B1, C0, C2, D0, D1, D2, E0, E2, E4, F1, F3, F5;
    ICOS_Float G0, G2, G4, G6, H0, H2, H4, H6, P0, P2, P4, P6, P8;
    ICOS_Float Q1, Q3, Q5, Q7, R0, R2, W0, W2, W4, Z0, Z2, Z4, Z6, Z8;
    int fix_lwidth;
};
inline func_line_p new_voigt( int mol, int iso,
    double nu, double S, double Gair, double E, double n,
    double delta, int ipos, double threshold,
    int fix_dw, int fix_lw, int fix_fp ) {
  return new voigt( mol, iso, nu, S, Gair, E, n, delta, ipos, threshold,
                    fix_dw, fix_lw, fix_fp );
}

// class func_quad : public func_evaluator {
  // public:
    // func_quad( ICOS_Float q, ICOS_Float l, ICOS_Float c );
    // void evaluate(ICOS_Float x, ICOS_Float *a);
    // static const int q_idx, l_idx, c_idx; // 0, 1, 2
// };


// baseline functions. func_base is a virtual base class
// (perhaps not technically?)
class func_base : public func_evaluator {
  public:
    inline func_base(const char *name = "base") : func_evaluator(name) {}
    int uses_nu_F0; ///< 0 or 1
};

// This is the old func_base for SVDs as a function of x
// The file format is ICOS standard binary with the first
// element of each column being initial parameter settings
class func_base_svdx : public func_base {
  public:
    func_base_svdx( const char *filename );
    void evaluate( ICOS_Float x, ICOS_Float *a );
  private:
    ICOS_Float **baseline;
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
    void evaluate( ICOS_Float x, ICOS_Float *a );
    void init(ICOS_Float *a);
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
    ICOS_Float **vectors;
    ICOS_Float **dvdnu;
    ICOS_Float **polyvecs;
};

class func_base_input : public func_base {
  public:
    func_base_input( func_base *base );
    void evaluate( ICOS_Float x, ICOS_Float *a );
    void init(ICOS_Float *a);
  private:
};

extern func_base *pick_base_type(const char *filename, func_parameter *nu_F0);


//-------------
// func_skew.c
//-------------
class skew_data {
  public:
    skew_data();
    void set_n_params(int n_gp, int n_ap);
    ICOS_Float g;
    ICOS_Float *dg; // allocate n_params
    ICOS_Float *da; // allocate n_params
    ICOS_Float gN;
    ICOS_Float Power;
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
    void init(ICOS_Float *a);
    void evaluate(ICOS_Float x, ICOS_Float *a);
    int skew_samples();
    void dump_params(ICOS_Float *a, int indent);
  private:
    ICOS_Float N;
    int M;
    ICOS_Float R2, R2N, P_scale;
    skew_data *skew; // We will have M of these
    ICOS_Float prev_x;
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
    void init(ICOS_Float *a);
    void evaluate(ICOS_Float x, ICOS_Float *a);
    void dump_params(ICOS_Float *a, int indent);
  private:
    ICOS_Float N_Passes;
    int n_base_params, n_abs_params;
    func_base *basep;
    func_abs *absp;
};

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

#endif
