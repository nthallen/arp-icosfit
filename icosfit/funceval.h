#ifndef FUNCEVAL_H_INCLUDED
#define FUNCEVAL_H_INCLUDED

#include <stdio.h>
#include <vector>
#include <stdint.h>
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
    inline parameter(int idx, const char *iname = 0)
      : index(idx), name(iname) {}
    int index;
    const char *name;
    ICOS_Float dyda;
    std::vector<paramref> refs;
};

class func_evaluator;
class func_line;

class argref {
  public:
    inline argref(func_evaluator *arg_in, unsigned int rn) :
      arg(arg_in), refnum(rn), dyda(0) {}
    func_evaluator *arg;
    unsigned int refnum;
    ICOS_Float dyda; ///< Partial with respect to this argument
};

class evaluation_order {
  public:
    std::vector<func_evaluator*> order;
    void set(func_evaluator *func, bool top = true, bool clear = false);
    void set_pre_order(func_evaluator *func, bool top = true, bool clear = false);
    void set_children(func_evaluator *func);
    void add(func_evaluator *func);
    void evaluate(ICOS_Float x, ICOS_Float *a);
    void pre_eval(ICOS_Float x, ICOS_Float *a);
    void init(ICOS_Float *a);
    void dump();
    int adjust_params(ICOS_Float alamda,
      ICOS_Float P, ICOS_Float T, ICOS_Float *a);
};

class func_evaluator {
  public:
    func_evaluator(const char *name, bool indexed = false, int index = 0); 
    //static void evaluate_all(std::vector<func_evaluator*> &order,
    //    ICOS_Float x, ICOS_Float *a);
    virtual void evaluate(ICOS_Float x, ICOS_Float *a);
    static void pre_eval_all(ICOS_Float x, ICOS_Float *a);
    virtual void pre_eval(ICOS_Float x, ICOS_Float *a);
    virtual void evaluate_partials();
    void init( ICOS_Float *a, int *ia );
    virtual void init(ICOS_Float *a);
    void append_func(func_evaluator *newfunc);
    virtual unsigned int adopted(func_evaluator *new_parent);
    virtual void fix_float_param(bool float_it, unsigned int refnum);
    virtual bool param_fixed();
    inline bool param_fixed(int i) {
      return args[i].arg->param_fixed(); }
    inline void fix_param(int i) {
      args[i].arg->fix_float_param(false,args[i].refnum); }
    inline void float_param(int i) {
      args[i].arg->fix_float_param(true,args[i].refnum); }
    inline ICOS_Float get_arg( ICOS_Float *a, int idx ) {
      return args[idx].arg->value;
    }
    virtual ICOS_Float set_param(ICOS_Float *a, ICOS_Float value);
    inline ICOS_Float set_param(ICOS_Float *a, int idx, ICOS_Float value) {
      return args[idx].arg->set_param(a, value);
    }
    virtual int adjust_params(ICOS_Float alamda, ICOS_Float P,
      ICOS_Float T, ICOS_Float *a);
    virtual func_line *is_line();
    virtual int line_check(int include, ICOS_Float& start, ICOS_Float& end,
                            ICOS_Float P, ICOS_Float T, ICOS_Float *a);
    virtual int skew_samples();
    virtual void dump_params();
    virtual void dump_partials();
    void print_indent( FILE *fp, int indent );
    virtual void print_config(FILE *fp);
    virtual void print_intermediates(FILE *fp);
    virtual void output_params(FILE *ofp, bool fixed);
    // void set_evaluation_order(std::vector<func_evaluator*> &order,
    //    bool top = true, bool clear = false);

    // int n_params;
    std::vector<argref> args;
    std::vector<parameter> params;
    // parameter *params;
    ICOS_Float value;
    const char *name;
    func_evaluator *parent;
    unsigned int n_references;
    // func_evaluator *first;
    // func_evaluator *last;
    // func_evaluator *next;
    static evaluation_order global_evaluation_order;
    static evaluation_order pre_evaluation_order;
    static evaluation_order dump_evaluation_order;
    bool added_to_eval;
};

class func_parameter : public func_evaluator {
  public:
    func_parameter(const char *name, ICOS_Float init_value,
                   bool indexed = false, int index = 0);
    void init(ICOS_Float *a);
    void fix_float_param(bool float_it, unsigned int refnum);
    bool param_fixed();
    ICOS_Float set_param(ICOS_Float *a, ICOS_Float value);
    void evaluate( ICOS_Float x, ICOS_Float *a );
    void evaluate_partials();
    void dump_params();
    static inline void set_ia(int *ia) { func_parameter::ia = ia; }
    void output_params(FILE *ofp, bool fixed);
   
    int index; ///< Parameter's global index
    static int n_parameters;
  protected:
    ICOS_Float init_val; ///< Initialization value
    uint32_t refs_float; ///< Bit-mapped
  private:
    static int *ia; ///< The 1-based vector of parameter fix/float settings
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

/**
 * func_line objects share location and width parameters.
 * Derived classes include gaussian, lorentzian and voigt,
 * which each add additional parameters.
 *
 * nu - an approximation of the line wavenumber
 * nu0 - static offset to extend precision of wavenumbers
 *       Set when the first line is defined.
 * nu1 - nu-nu0 (or a more precise version thereof)
 * nu_P - nu1+delta*P/760. : pressure-shifted wavenumber less nu0
 */
class func_line : public func_evaluator {
  public:
    func_line( const char *name, int mol, int iso,
      double nu, double S, double Gair, double E, double n,
      double delta, unsigned int ipos, double threshold, int fix_w,
      int fix_fp, func_parameter *N = 0 );
    ~func_line();
    unsigned int adopted(func_evaluator *new_parent);
    int adjust_params( ICOS_Float alamda, ICOS_Float P,
        ICOS_Float T, ICOS_Float *a );
    func_line *is_line();
    static const int dnu_idx, w_idx, n_idx;
    int nu_F0_idx;
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
    void print_intermediates(FILE *fp);
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

/**
 * func_abs has 1 common parameter, which is nu_F0
 * In Release 2.2, we inherit nu_F0 from func_skew.
 * Each line then has it's own dnu that func_abs owns
 * and its own 4 parameters.
 * In Release 3.0, each line has both nu_F0 and dnu.
 * Lines own dnu, nu_F0 is owned by func_abs.
 */
class func_abs : public func_evaluator {
  public:
    func_abs();
    void init(ICOS_Float *a);
    void evaluate(ICOS_Float x, ICOS_Float *a);
    int adjust_params(ICOS_Float alamda, ICOS_Float P, ICOS_Float T, ICOS_Float *a);
    void print_config(FILE *fp);
    void print_intermediates(FILE *fp);
};
typedef func_abs *func_abs_p;
typedef func_line *func_line_p;
typedef func_parameter *func_parameter_p;
inline func_abs_p new_func_abs() { return new func_abs(); }
inline func_abs_p abs_append( func_abs_p abs, func_line_p line ) {
  abs->append_func( line );
  return abs;
}
func_parameter_p Nparameter(const char *group);

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
          int fix_dw, int fix_lw, int fix_fp, func_parameter *N = 0 );
    void init(ICOS_Float *a);
    void evaluate(ICOS_Float x, ICOS_Float *a);
    ICOS_Float line_width(ICOS_Float*a);
    ICOS_Float line_start(ICOS_Float*a);
    ICOS_Float line_end(ICOS_Float*a);
    static const int gl_idx; // 3
    int adjust_params( ICOS_Float alamda, ICOS_Float P, ICOS_Float T, ICOS_Float *a );
    void line_fix();
    void line_float();
    //void dump_params(ICOS_Float *a, int indent);
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
    int fix_dw, int fix_lw, int fix_fp, func_parameter *NP ) {
  return new voigt( mol, iso, nu, S, Gair, E, n, delta, ipos, threshold,
                    fix_dw, fix_lw, fix_fp, NP);
}

/** baseline functions. func_base is a virtual base class
 */
class func_base : public func_evaluator {
  public:
    inline func_base(const char *name = "base") :
      func_evaluator(name), uses_nu_F0(0) {}
    int uses_nu_F0; ///< 0 or 1
};

/**
 * This is the old func_base for SVDs as a function of x
 * The file format is ICOS standard binary with the first
 * element of each column being initial parameter settings
 */
class func_base_svdx : public func_base {
  public:
    func_base_svdx( const char *filename );
    void evaluate( ICOS_Float x, ICOS_Float *a );
  private:
    ICOS_Float **baseline;
    int n_pts;
};

/**
 * This is a function supporting vectors as a function of
 * nu as well as polynomials as a function of x or nu
 * The file format is:
 * i*4  0  To distinguish from func_base_svdx format
 * i*4  1  To identify func_base_ptbnu format
 * f*8 polynomial scale factor
 * f*8 nu0 smallest value of nu
 * f*8 dnu increment between samples
 * i*2 n_vectors
 * i*2 npts per vector
 * i*2 number of polynomial coefficients (degree+1)
 * i*2 polynomial of x == 0, nu == 1
 * f*4 X n_vectors initial parameter values
 * f*4 X n polynomial coefficients initial parameter values
 * f*4 X npts X n_vectors: vector data in column order

 * Can use global SignalRegion to determine range of x
 */
class func_base_ptbnu : public func_base {
  public:
    func_base_ptbnu(const char *filename, func_evaluator *nu_F0);
    void evaluate(ICOS_Float x, ICOS_Float *a);
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

extern func_base *pick_base_type(const char *filename, func_evaluator *nu_F0);


//-------------
// func_skew.c
//-------------
class func_beta : public func_evaluator {
  public:
    func_beta(func_abs *abs);
    void evaluate(ICOS_Float x, ICOS_Float *a);
};

class func_gamma : public func_evaluator {
  public:
    func_gamma(func_beta *beta);
    void evaluate(ICOS_Float x, ICOS_Float *a);
  private:
    ICOS_Float R2;
};

class func_epsilon : public func_evaluator {
  public:
    func_epsilon(func_gamma *gamma);
    void evaluate(ICOS_Float x, ICOS_Float *a);
  private:
    ICOS_Float N;
};

class func_delta : public func_evaluator {
  public:
    func_delta(func_epsilon *epsilon, func_gamma *gamma);
    void evaluate(ICOS_Float x, ICOS_Float *a);
};

class func_g : public func_evaluator {
  public:
    func_g(func_base *base, func_beta *beta, func_delta *delta);
    void evaluate(ICOS_Float x, ICOS_Float *a);
};

class skew_data {
  public:
    skew_data();
    void set_n_params(int n_gp, int n_epsp);
    ICOS_Float g;
    ICOS_Float *dg; // allocate n_params
    ICOS_Float eps;
    ICOS_Float *deps; // allocate n_params
    ICOS_Float Power;
    int initialized;
    int n;
};

/**
 * func_skew applies the cell skew function of ICOS. Its two
 * children define the input power curve (base) and the
 * intra-cavity absorption (abs).
 * base and abs have mostly independent parameters, but if
 * the base function uses nu_F0, that parameter is shared.
 * The skew member is used as an M-element circular buffer.
 */
class func_skew : public func_evaluator {
  public:
    func_skew(func_g *g, func_epsilon *epsilon);
    void init(ICOS_Float *a);
    void pre_eval(ICOS_Float x, ICOS_Float *a);
    void evaluate(ICOS_Float x, ICOS_Float *a);
    void evaluate_partials();
    int skew_samples();
    void output_params(FILE *ofp, bool fixed);
  private:
    void sub_eval(ICOS_Float x, ICOS_Float *a);
    ICOS_Float N;
    int M;
    ICOS_Float R2, R2N, P_scale;
    skew_data *skew; // We will have M of these
    int *depsi; // map from g parameters to eps parameters
    int skewidx;
    func_evaluator *basep;
    evaluation_order skew_eval_order;
};

/**
 * func_noskew calculates absorption for a simple multi-pass
 * cell. Its two children define the input power curve (base)
 * and the intra-cavity absorption (abs).
 * base and abs have mostly independent parameters, but if
 * the base function uses nu_F0, that parameter is shared.
 */
class func_noskew : public func_evaluator {
  public:
    func_noskew(func_base *base, func_abs *abs);
    void evaluate(ICOS_Float x, ICOS_Float *a);
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
