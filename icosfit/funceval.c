#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nortlib.h"
#include "funceval.h"
#include "global.h"
#include "ptread.h"

func_evaluator::func_evaluator( const char *sname, int n ) {
  n_params = n;
  name = sname;
  if ( n_params == 0 ) params = 0;
  else params = new parameter[n_params];
  first = last = next = parent = 0;
  // printf("func_evaluator( %s, %d );\n", sname, n );
}

void func_evaluator::append_func( func_evaluator *newfunc ) {
  // printf( "func_evaluator::append_func(%s, %s)\n", name, newfunc->name );
  if ( this->last == 0 ) {
    this->first = this->last = newfunc;
  } else {
    this->last->next = newfunc;
    this->last = newfunc;
  }
  newfunc->parent = this;
}

// func_evaluator::init() initializes the a and ia pointers in
// the children's params list. The object's own pointers were
// either initialized by their parent or by the init(a,ia)
// call. This base version assumes all the children have
// independent parameters and distributes its own parameters
// sequentially to the children, starting at the specified
// index.
// When init() is called, we are guaranteed that params have
// been allocated for this object, but we may still need to
// allocate params for children of this object.
void func_evaluator::init(float *a, int p1) {
  func_evaluator *child;
  int p2;


  // printf( "func_evaluator::init(%s, %d); n_params=%d\n", name, p1, n_params );
  for ( p2 = 0; p2 < n_params; p2++ )
    a[params[p2].index] = params[p2].init;
  for ( child = first; child != 0; child = child->next ) {
    if ( p1 + child->n_params > n_params ) {
      fprintf( stderr, "Too many child params: n_params = %d\n", n_params );
      exit(1);
    }
    if ( child->params == 0 )
      child->params = new parameter[child->n_params];
    for ( p2 = 0; p2 < child->n_params; p2++ ) {
      link_param( p1, child, p2 );
      p1++;
    }
    child->init(a);
  }
}

void func_evaluator::init(float *a) {
  init( a, 0 );
}

// a and ia are the 1-based vectors from mrqmin. i.e.
// their 0th element is unused.
// init is called once when the fitdata object is
// created, and it recursively calls init(a) on all
// its children. The entire function structure has
// been assembled at this point, so the init(a)
// methods are a good place to do initializations
// that need to assume that.
void func_evaluator::init( float *a, int *ia ) {
  int i;

  // printf( "func_evaluator::init( %s, a, ia );\n", name );
  if ( params == 0 ) params = new parameter[n_params];
  for ( i = 0; i < n_params; i++ ) {
    params[i].index = i+1;
    params[i].ia = &ia[i+1];
  }
  init(a); // initialize the children
}

// evaluate is used as a helper function to the
// routines in the derived classes that do the
// real work. This routine just evaluates all the
// children so their results are available to
// their parents.
void func_evaluator::evaluate(float x, float *a) {
  func_evaluator *child;
  for ( child = first; child != 0; child = child->next ) {
    child->evaluate(x, a);
  }
}

int func_evaluator::adjust_params( float alamda, float P, float T, float *a ) {
  func_evaluator *child;
  int rv = 0;

  for ( child = first; child != 0; child = child->next ) {
    if ( child->adjust_params( alamda, P, T, a ) )
      rv = 1;
  }
  return rv;
}

/* clamp_param_high() is a method to be called from adjust_params.
   It requires that the specified parameter has a high value set.
   It will not take any action if the parameter is fixed.
*/
void func_evaluator::clamp_param_high( float *a, int idx ) {
  if ( ! param_fixed(idx) ) {
    float val = get_param(a,idx);
    float limit = params[idx].high;
    if ( val > limit ) {
      val = (params[idx].prev + limit)/2;
      set_param( a, idx, val );
    }
    params[idx].prev = val;
  }
}

/* clamp_param_low() is a method to be called from adjust_params.
   It requires that the specified parameter has a low value set.
   It will not take any action if the parameter is fixed.
*/
void func_evaluator::clamp_param_low( float *a, int idx ) {
  if ( ! param_fixed(idx) ) {
    float val = get_param(a,idx);
    float limit = params[idx].low;
    if ( val < limit ) {
      val = (params[idx].prev + limit)/2;
      set_param( a, idx, val );
    }
    params[idx].prev = val;
  }
}

/* clamp_param_highlow() is a method to be called from adjust_params.
   It requires that the specified parameter has a high value set.
   It will not take any action if the parameter is fixed.
*/
void func_evaluator::clamp_param_highlow( float *a, int idx ) {
  if ( ! param_fixed(idx) ) {
    float val = get_param(a,idx);
    float limit = params[idx].high;
    if ( val > limit ) {
      val = (params[idx].prev + limit)/2;
      set_param( a, idx, val );
    } else {
      limit = params[idx].low;
      if ( val < limit ) {
        val = (params[idx].prev + limit)/2;
        set_param( a, idx, val );
      }
    }
    params[idx].prev = val;
  }
}


//---------------------------------------------------------
// aggregate
//---------------------------------------------------------
void func_aggregate::append_func( func_evaluator *newfunc ) {
  // printf( "func_aggregate::append_func(%s, %s)\n", name, newfunc->name );
  func_evaluator::append_func(newfunc);
  n_params += newfunc->n_params;
}

//---------------------------------------------------------
// sum
//  Allow inheritance with a number of fixed parameters.
//---------------------------------------------------------
void func_sum::evaluate(float x, float *a) {
  evaluate( x, a, 0 );
}

void func_sum::evaluate(float x, float *a, int i) {
  int j;
  func_evaluator *child;

  func_evaluator::evaluate(x, a);
  value = 0.;
  for ( child = first; child != 0; child = child->next ) {
    value += child->value;
    for ( j = 0; j < child->n_params; j++ )
      params[i++].dyda = child->params[j].dyda;
  }
}

//---------------------------------------------------------
// product
//---------------------------------------------------------
void func_product::evaluate(float x, float *a) {
  int i = 0, j;
  func_evaluator *child, *child2;

  func_evaluator::evaluate(x, a);
  value = 1.;
  for ( child = first; child != 0; child = child->next )
    value *= child->value;
  for ( child = first; child != 0; child = child->next ) {
    float prod;
    if ( child->value == 0. ) {
      prod = 1.;
      for ( child2 = first; child2 != 0; child2 = child2->next )
        if ( child2 != child ) prod *= child2->value;
    } else prod = value / child->value;
    for ( j = 0; j < child->n_params; j++ )
      params[i++].dyda = prod * child->params[j].dyda;
  }
}

static int get_molwt( int isotopomer ) {
  switch (isotopomer) {
    case 11: return 18; // H_2O
    case 12: return 20; // H_2{}^{18}O
    case 13: return 19; // H_2{}^{17}O
    case 14: return 19; // HDO
    case 21: return 44; // CO2
    case 22: return 45; // CO_2{}^{13}C
    case 23: return 46; // CO_2{}^{18}O
    case 24: return 45; // CO_2{}^{17}O
    case 41: return 44; // N2O
    case 42: return 45; // N^{15}NO
    case 43: return 45; // NN^{15}O
    case 61: return 16; // CH_4
    case 62: return 17; // C13H4
    case 63: return 18; // CH3D
    default:
      nl_error( 3,
        "Uncatalogued isotopomer '%d': Edit funceval.c get_molwt()",
         isotopomer );
      return 0; // Never reached
  }
}


//---------------------------------------------------------
// func_line object has at least 3 parameters
//   l_idx = 0: Fine location in cm-1
//   w_idx = 1: Doppler E-folding halfwidth in cm-1
//   n_idx = 2: Number Density in molecules/cm3
// These parameters are all initialized to zero.
//
// There are also a raft of fixed parameters based on HITRAN
// data for individual lines:
// molecule
// isotope
// nu = wavenumber
// S = spectral line intensity
// G_air = Gamma_air
// E = Lower state energy of the transition
// n = coefficient of temerature dependence or air-broadened hw
// delta = air-broadened pressure shift
// ipos = sample number of line center at start
//---------------------------------------------------------
const int func_line::l_idx = 0;
const int func_line::w_idx = 1;
const int func_line::n_idx = 2;
const double func_line::DRTPI = 0.5641895835477563; // 1/SQRT(pi)
const double func_line::Tref = 296.; // K
const double func_line::C2 = 1.4388; // cm K second radiation constant hc/k
double func_line::nu0 = 0.;
int func_line::n_lines = 0;

func_line::func_line( const char *name, int np, int mol, int iso,
          double nu_in, double S_in, double G_air_in, double E_in,
          double n_in, double delta_in, unsigned int ipos_in, double threshold,
          int fix_w, int fix_fp ) :
  func_evaluator( name, np ) {
  params[l_idx].init = 0.;
  params[w_idx].init = 1.;
  params[n_idx].init = 0.;
  line_number = ++n_lines;
  fixed = 0;
  fix_finepos = fix_fp;
  fix_width = fix_w;
  prev_numdens = 0.;
  prev_ged = 1.;
  rolledback = 0;
  isotopomer = mol*10 + iso;
  molwt = get_molwt(isotopomer);
  if ( nu0 == 0. ) nu0 = floor(nu_in);
  nu = nu_in;
  nu1 = nu_in - nu0;
  S = S_in;
  G_air = G_air_in;
  E = E_in;
  n_air = n_in;
  delta = delta_in;
  S_thresh = threshold;
  ipos = ipos_in;
  Corr_Tref = 1/(exp(-C2 * E / Tref ) * (1-exp(-C2*nu/Tref)));
}

void func_line::init(float *a) {
  if ( fix_width ) fix_param( w_idx );
  func_evaluator::init(a);
}

void func_line::print_config( FILE *fp ) {
  fprintf( fp,
    "  %d %d %.6lf %.4e %.4f %.4f %.2f %.6f %d\n",
    isotopomer/10, isotopomer%10, nu1+nu0, S, G_air, E,
    n_air, delta, n_params );
}

int func_line::adjust_params( float alamda, float P, float T, float *a ) {
  // Eliminated a check for drifting. Taken care of in func_abs.
  if ( alamda < -1.5 ) {
    double Spt = S * pow(Tref/T, 1.5) * exp(-C2*E/T) * (1-exp(-C2*nu/T))
            * Corr_Tref;
    Ks = Spt * GlobalData.CavityLength * DRTPI; 
    nu_P = nu1 + delta * (1 - P/760.);
    rolledback = 0;
  }
  float numdens = get_param( a, n_idx );
  // if ( numdens < 0. ) numdens = set_param( a, n_idx, prev_numdens/2. );
  // prev_numdens = numdens;
  float gamma_ed;
  // const float ln2 = log(2.);
  if ( param_fixed( w_idx ) ) {
    if ( alamda < 0 ) {
      // The following formula comes from Liz's old lecture notes
      // (Webster?)
      // Value is stored in wavenumbers. It must be
      // multiplied by s_wn (param sc_idx) to get samples.
      // This will probably have to be adjusted
      // to take into account the laser line width
      // gamma_ed = set_param( a, w_idx,
      //  3.581e-7 * nu * sqrt(T/(ln2*molwt)));
      gamma_ed = set_param( a, w_idx,
        4.30213e-7 * nu * sqrt(T/molwt));
      prev_ged = gamma_ed;
    } else {
      gamma_ed = get_param(a, w_idx);
    }
  } else {
    const float min_ged = 1e-3;
    gamma_ed = get_param(a, w_idx);
    if ( gamma_ed <= min_ged )
      gamma_ed = set_param(a, w_idx, (prev_ged-min_ged)/2 + min_ged );
    prev_ged = gamma_ed;
  }
  if ( alamda < 0 && param_fixed( l_idx ) ) {
    // I'm pretty sure this is handled by func_abs::adjust_params
    // set_param( a, l_idx, 0. );
  }
  if ( fixed ) {
    if ( alamda == 0 ) {
      float strength = gamma_ed > 0 ? Ks * numdens / gamma_ed : 0.;
      if ( strength > S_thresh * 4. ) {
        if ( rolledback < 2 ) {
          nl_error( 0, "Floating line %d (strength %g)",
                      line_number, strength );
          line_float();
          return 1;
        } else nl_error( 0, "NOT re-floating line %d",
                  line_number );
      }
    }
  } else {
    float strength = gamma_ed > 0 ? Ks * numdens / gamma_ed : 0.;
    if ( strength <= S_thresh ) {
      nl_error( 0, "Fixing line %d (strength %g)",
                    line_number, strength );
      line_fix();
      if (alamda >= 0 ) rolledback++;
      return 1;
    }
  }
  return 0;
}

float func_line::line_start(float *a) {
  return (nu_P - GlobalData.RightLineMarginMultiplier*get_param(a, w_idx));
}
float func_line::line_end(float *a) {
  return (nu_P + GlobalData.LeftLineMarginMultiplier*get_param(a, w_idx));
}

int func_evaluator::line_check(int include, float& start, float& end,
                float P, float T, float *a) {
  if ( first != 0 && first->line_check( include, start, end, P, T, a ) )
    return 1;
  if ( next != 0 && next->line_check( include, start, end, P, T, a ) )
    return 1;
  return 0;
}

int func_evaluator::skew_samples() {
  int rv = 0;
  int frv;
  if ( first != 0 ) {
    frv = first->skew_samples();
    if ( frv > rv ) rv = frv;
  }
  if ( next != 0 ) {
    frv = next->skew_samples();
    if ( frv > rv ) rv = frv;
  }
  return rv;
}

// dump_params() is a diagnostic tool that is invoked
// at times of failure. The default version simply lists
// all the parameters and their values without recursing
// to children, but overrides can delegate to children.
void func_evaluator::dump_params(float *a, int indent) {
  int i;
  print_indent( stderr, indent );
  fprintf( stderr, "Parameters for '%s':\n", name );
  indent += 2;
  for ( i = 0; i < n_params; i++ ) {
    int indx = params[i].index;
    print_indent( stderr, indent );
    fprintf( stderr, "[%2d] %g\n", indx, a[indx] );
  }
}

void func_evaluator::print_indent( FILE *fp, int indent ) {
  while (indent-- > 0) fputc(' ', fp );
}

// line_check(include, start, end, P, T, a );
// operates in two passes. First, include is set to 0 to
// indicate the 'exclude' step. A line is excluded if it
// hits the start or end boundary, and if so, the boundaries
// are moved in to also exclude any significant portion
// of this line. If the boundaries need to be moved, we
// return 1, and another exclude pass is indicated. This
// is necessary to deal with overlapping lines.
//
// Refinement to the exclusion: If a line hits the boundary,
// we first try to re-evaluate it's width by fixing the line
// and calling adjust_params. If that helps, then its threshold
// needs to be raised. If it doesn't help, then turn it off.
//
// Once a final set of lines has been determined, include
// is set to 1 to indicate the 'include' step. Here the
// boundaries are expanded to include all the lines which
// are still 'on'.
int func_line::line_check(int include, float& start, float& end,
                    float P, float T, float *a ) {
  float ls = line_start(a);
  float le = line_end(a);
  func_line *next = lnext();
  if ( ! include ) {
    int rv = -1;
    if ( ! fixed && ( ls < start || le > end ) ) {
      float save_thresh = S_thresh;
      line_fix();
      S_thresh = Ks * get_param(a, n_idx)*2/get_param(a, w_idx);
      adjust_params( -1, P, T, a );
      ls = line_start(a);
      le = line_end(a);
      if ( ls >= start && le <= end ) {
        nl_error(0, "Raised threshold on line %d near boundary",
                      line_number );
      } else {
        nl_error( 0, "Fixing line %d (%.4f,%.4f)",
                          line_number, ls, le );
        S_thresh = save_thresh;
      }
    }
    if ( le < start || ls > end ) {
      float lem = le+GlobalData.LeftLineMargin;
      float lsm = ls-GlobalData.RightLineMargin;
      rv = 0;
      if ( ls < start && lem > start ) {
        start = lem; rv = 1;
        if ( GlobalData.Verbosity & 2 )
          nl_error( 0, "Exclude: Updated start to %.4f", start );
      }
      if ( le > end && lsm < end ) {
        end = lsm; rv = 1;
        if ( GlobalData.Verbosity & 2 )
          nl_error( 0, "Exclude: Updated end to %.4f", end );
      }
      if ( ! param_fixed(n_idx) ) {
        nl_error( 0, "Turning off line %d (%.4f,%.4f)",
                          line_number, ls, le );
        fix_param(n_idx);
        set_param( a, n_idx, 0. );
        line_fix();
      }
      if ( rv != 0 ) return rv;
    }
    if ( next != 0 && next->line_check( include, start, end, P, T, a ) )
      return 1;
    if ( rv < 0 && param_fixed(n_idx) ) {
      nl_error( 0, "Turning on line %d (%.4f,%.4f)",
                          line_number, ls, le );
      float_param(n_idx);
      // We don't actually line_float() until the fit raises the
      // number density high enough.
    }
  } else {
    if ( next != 0 && next->line_check( include, start, end, P, T, a ) )
      return 1;
    if ( ! param_fixed(n_idx) ) {
      if ( start == 0. || ls-GlobalData.RightLineMargin < start )
        start = ls-GlobalData.RightLineMargin;
      if ( end == 0. || le+GlobalData.LeftLineMargin > end )
        end = le+GlobalData.LeftLineMargin;
    }
  }
  return 0;
}

void func_line::line_fix() {
  func_abs *p = (func_abs *)parent;
  fix_param(l_idx);
  p->fix_linepos(line_number);
  fix_param(w_idx);
  fixed = 1;
}

void func_line::line_float() {
  func_abs *p = (func_abs *)parent;
  // float_param(l_idx);
  if ( fix_finepos == 0 ) p->float_linepos(line_number);
  if ( fix_width == 0 ) float_param(w_idx);
  fixed = 0;
}


//---------------------------------------------------------
// func_abs: in funcabs.cc
//---------------------------------------------------------


//---------------------------------------------------------
// gaussian
//---------------------------------------------------------
// Evaluation of this will get out of hand if
// w becomes small, or more specifically if
// w^2 becomes small
// This form of the gaussian is normalized in
// various ways presumably to match the lorentzian.
// ### This needs to be fixed for cm-1 scale
void gaussian::evaluate(float x, float *a) {
  static const float four_log_2 = 4 * log(2.);
  static const float fl2_pi = sqrt(four_log_2/M_PI);
  float dnu = a[params[l_idx].index];
  float w = a[params[w_idx].index];
  float s = a[params[n_idx].index];
  assert( s >= 0. );
  float w2 = w * w;
  float v1 = fl2_pi * s / w;
  int xx = int(x);
  float diff = ICOSfile::wndata->data[xx] - nu_P + dnu;
  float diff2 = diff * diff;
  float v2 = exp( - diff2 * four_log_2 / w2 );
  value = v1 * v2;
  params[n_idx].dyda = fl2_pi * v2 / w;
  params[l_idx].dyda = 2*v1*v2*diff*four_log_2/w2;
  float v3 = fl2_pi * s * v2 / w2;
  params[w_idx].dyda = -(v3/w) + v3 * 2 * diff2 * four_log_2 / w2;
}

//---------------------------------------------------------
// lorentzian
//---------------------------------------------------------
// ### This needs to be fixed for wavenumber scale
void lorentzian::evaluate(float x, float *a) {
  float dnu = a[params[l_idx].index];
  float w = a[params[w_idx].index];
  float s = a[params[n_idx].index];
  int xx = int(x);
  float diff = ICOSfile::wndata->data[xx] - nu_P + dnu;
  float v1 = M_PI * ( 4 * diff *diff + w * w );
  float v12 = v1 * v1;
  value = 2 * w * s / v1;
  params[n_idx].dyda = 2 * w / v1;
  params[l_idx].dyda = 16*M_PI*w*s*diff / v12;
  params[w_idx].dyda = (2*s*v1 - 4*M_PI*s*w*w)/v12;
}

//----------------------------------------------------------
// voigt: in humdev.cc
//----------------------------------------------------------

//----------------------------------------------------------
// func_quad: second-order polynomial fit
//----------------------------------------------------------
const int func_quad::q_idx = 0, func_quad::l_idx = 1, func_quad::c_idx = 2;

func_quad::func_quad(float q, float l, float c) : func_evaluator("quad",3) {
  params[q_idx].init = q;
  params[l_idx].init = l;
  params[c_idx].init = c;
}

void func_quad::evaluate(float x, float *a) {
  float q = a[params[q_idx].index];
  float l = a[params[l_idx].index];
  float c = a[params[c_idx].index];
  value = q*x*x + l*x + c;
  params[q_idx].dyda = x*x;
  params[l_idx].dyda = x;
  params[c_idx].dyda = 1;
}
