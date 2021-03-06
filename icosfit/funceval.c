#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <set>
#include "nortlib.h"
#include "funceval.h"
#include "global.h"
#include "ptread.h"

evaluation_order func_evaluator::global_evaluation_order;
evaluation_order func_evaluator::pre_evaluation_order;
evaluation_order func_evaluator::dump_evaluation_order;

void evaluation_order::set(func_evaluator*func, bool top, bool clear) {
  if (top) {
    set(func, false, true);
  }
  std::vector<argref>::iterator arg;
  for (arg = func->args.begin(); arg != func->args.end(); ++arg) {
    set(arg->arg, false, clear);
  }
  if (clear) {
    func->added_to_eval = false;
  } else if (!func->added_to_eval) {
    add(func);
    func->added_to_eval = true;
  }
}

void evaluation_order::set_pre_order(func_evaluator*func, bool top, bool clear) {
  if (top) {
    set(func, false, true);
  }
  if (clear) {
    func->added_to_eval = false;
  } else if (!func->added_to_eval) {
    add(func);
    func->added_to_eval = true;
  }
  std::vector<argref>::iterator arg;
  for (arg = func->args.begin(); arg != func->args.end(); ++arg) {
    set_pre_order(arg->arg, false, clear);
  }
}

void evaluation_order::set_children(func_evaluator*func) {
  set(func, false, true);
  func->added_to_eval = true;
  set(func, false, false);
}

void evaluation_order::add(func_evaluator*func) {
  order.push_back(func);
}

void evaluation_order::evaluate(ICOS_Float x, ICOS_Float *a) {
  std::vector<func_evaluator*>::iterator func;
  for (func = order.begin(); func != order.end(); ++func) {
    (*func)->evaluate(x, a);
    (*func)->evaluate_partials();
  }
}

void evaluation_order::pre_eval(ICOS_Float x, ICOS_Float *a) {
  std::vector<func_evaluator*>::iterator func;
  for (func = order.begin(); func != order.end(); ++func) {
    (*func)->pre_eval(x, a);
  }
}

void evaluation_order::init(ICOS_Float *a) {
  std::vector<func_evaluator*>::iterator func;
  for (func = order.begin(); func != order.end(); ++func) {
    (*func)->init(a);
  }
}

int evaluation_order::adjust_params(ICOS_Float alamda, ICOS_Float P,
      ICOS_Float T, ICOS_Float *a) {
  int rv = 0;
  std::vector<func_evaluator*>::iterator func;
  for (func = order.begin(); func != order.end(); ++func) {
    if ((*func)->adjust_params(alamda, P, T, a))
      rv = 1;
  }
  return rv;
}

void evaluation_order::dump() {
  fprintf(stderr, "Evaluation State:\n");
  std::vector<func_evaluator*>::iterator func;
  for (func = order.begin(); func != order.end(); ++func) {
    (*func)->dump_params();
    (*func)->dump_partials();
  }
}

/**
 * Evaluates the partial derivatives with respect to
 * all relevant parameters. func_evaluator::evaluate()
 * is responsible for calculating partials with respect
 * to the arguments, which may be functions of other
 * arguments and parameters. This finishes the calculation
 * with respect to the underlying parameters.
 */
void func_evaluator::evaluate_partials() {
  std::vector<parameter>::iterator p;
  for (p = params.begin(); p != params.end(); ++p) {
    p->dyda = 0;
    std::vector<paramref>::iterator ref;
    for (ref = p->refs.begin(); ref != p->refs.end(); ++ref) {
      p->dyda += 
        args[ref->arg_num].dyda *
        args[ref->arg_num].arg->params[ref->param_num].dyda;
    }
  }
}

func_evaluator::func_evaluator(const char *sname, bool indexed, int idx) {
  if (indexed) {
    char buf[80];
    snprintf(buf, 80, "%s[%d]", sname, idx);
    name = nl_strdup(buf);
  } else {
    name = sname;
  }
  parent = 0;
  value = 0.;
  n_references = 0;
  added_to_eval = false;
}

/**
 * @param newfunc The new argument, appended to args vector.
 * @param owner if true, newfunc's parent is set
 */
void func_evaluator::append_func( func_evaluator *newfunc) {
  args.push_back(argref(newfunc, newfunc->adopted(this)));
}

/**
 * @param new_parent The called object is now a child of the new_parent object
 * By overloading this function, children can gain access to needed
 * information from the parents. Specifically, line types can access the
 * parent func_abs's nu_F0 parameter.
 */
unsigned int func_evaluator::adopted(func_evaluator *new_parent) {
  if (parent == 0) {
    parent = new_parent;
  }
  return n_references++;
}

/**
 * @param a The initial parameter value vector
 * Besides finalizing some important internal statistics,
 * init(a) also defines the evaluation order.
 */
void func_evaluator::init(ICOS_Float *a) {
  std::vector<argref>::iterator child;
  std::set<int> pidx;
  
  for (child = args.begin(); child != args.end(); ++child) {
    // child->arg->init(a);
    // child now has params defined
    std::vector<parameter>::iterator cp;
    for (cp = child->arg->params.begin(); cp < child->arg->params.end(); ++cp) {
      pidx.insert(cp->index);
    }
  }
  // Now we know exactly how many and which parameters we depend on
  std::set<int>::iterator ipidx;
  int i = 0;
  for (ipidx = pidx.begin(); ipidx != pidx.end(); ++ipidx, ++i) {
    params.push_back(parameter(*ipidx));
  }
  // Now go back through the children and record where the
  // child references appear
  unsigned int argi, argpi, pi;
  for (argi = 0; argi < args.size(); ++argi) {
    for (argpi = 0; argpi < args[argi].arg->params.size(); ++argpi) {
      for (pi = 0; pi < params.size(); ++pi) {
        if (params[pi].index == args[argi].arg->params[argpi].index) {
          params[pi].name = args[argi].arg->params[argpi].name;
          params[pi].refs.push_back(paramref(argi, argpi));
          break;
        }
      }
    }
  }
  for (child = args.begin(); child != args.end(); ++child) {
    std::vector<parameter>::iterator cp;
    for (cp = child->arg->params.begin(); cp != child->arg->params.end(); ++cp) {
      pidx.insert(cp->index);
    }
  }
}

/**
  @param a vector of paramter values
  @param ia vector of parameter enable flags
  
  a and ia are the 1-based vectors from mrqmin. i.e.
  their 0th element is unused.
  init is called once when the fitdata object is
  created, and it init(a) on the global_evaluation_order.
  The entire function structure has
  been assembled at this point, so the init(a)
  methods are a good place to do initializations
  that need to assume that.
 */
void func_evaluator::init( ICOS_Float *a, int *ia ) {
  // printf( "func_evaluator::init( %s, a, ia );\n", name );
  func_parameter::set_ia(ia);
  global_evaluation_order.init(a);
}

/**
 * @param float true to float, false to fix parameter.
 * This call is currently illegal unless object is a
 * func_parameter.
 */
void func_evaluator::fix_float_param(bool float_it, unsigned int refnum) {
  nl_error(3, "Illegal attempt to %s a non-parameter",
      float_it ? "float" : "fix");
}

/**
 * @return true if the parameter is fixed.
 * This call is currently illegal unless object is a
 * func_parameter.
 */
bool func_evaluator::param_fixed() {
  nl_error(3, "Illegal attempt to query fixed state of a non-parameter");
  return false;
}

ICOS_Float func_evaluator::set_param(ICOS_Float *a, ICOS_Float value) {
  nl_error(3, "Illegal attempt to set value of a non-parameter");
  return 0.;
}

/**
 * No-op method for virtual override by func_skew.
 */
void func_evaluator::pre_eval(ICOS_Float x, ICOS_Float *a) {}

/**
 * With version 3, evaluate() no longer recurses, since
 * the execution order is determined during init().
 * Sub classes will override this method and are tasked
 * with calculating value and the partial derivative with
 * respect to each of their arguments (args).
 */
void func_evaluator::evaluate(ICOS_Float x, ICOS_Float *a) {}

/**
 * With version 3, adjust_params() no longer recurses, since
 * the execution order is determined during init().
 * Sub classes will override this method and are tasked
 * with calculating value and the partial derivative with
 * respect to each of their arguments (args).
 * @return non-zero if a parameter value has been changed.
 */
int func_evaluator::adjust_params( ICOS_Float alamda, ICOS_Float P, ICOS_Float T, ICOS_Float *a ) {
  return 0;
}

/**
  @param ofp FILE pointer to ICOSsum.dat
  @param fixed true to output fixed/floating bool values, false to output parameter values.
  Default will recurse through arguments, but we
  will override for efficiency and consistency as
  necessary.
  Specific overrides are necessary for:
    - func_skew: base, abs
    - func_parameter: Just the value
 */
void func_evaluator::output_params(FILE *ofp, bool fixed) {
  std::vector<argref>::iterator child;
  for (child = args.begin(); child != args.end(); ++child) {
    child->arg->output_params(ofp, fixed);
  }
}

/**
 * @return zero for non-line classes
 */
func_line *func_evaluator::is_line() { return 0; }

int func_parameter::n_parameters = 0;
int *func_parameter::ia = 0;

func_parameter::func_parameter(const char *name, ICOS_Float init_value,
        bool indexed, int idx) : func_evaluator(name,indexed,idx) {
  index = ++n_parameters;
  init_val = init_value;
  refs_float = 0;
}

void func_parameter::init(ICOS_Float *a) {
  func_evaluator::init(a); // This should do nothing
  params.push_back(parameter(index, name));
  params.back().dyda = 1.0;
  a[index] = init_val;
  // ia is initialized as all floating, so we will
  // do the same on a per-reference basis. This means
  // references must ask to fix.
  refs_float = (1<<n_references)-1;
  func_evaluator::init(a);
}

/**
 * @param float_it true to float the parameter
 * @param refnum The reference number
 * If any reference wants the parameter to float, it will float.
 * All references must request a parameter to be fixed before it is
 * actually fixed.
 * Note that parameters cannot be fixed or floated until after init()
 * (or at least until init() has run on all dependents)
 */
void func_parameter::fix_float_param(bool float_it, unsigned int refnum) {
  if (float_it) {
    refs_float |= 1 << refnum;
    ia[params[0].index] = 1;
  } else {
    refs_float &= ~(1U << refnum);
    if (refs_float == 0)
    ia[params[0].index] = 0;
  }
}

bool func_parameter::param_fixed() {
  return ia[params[0].index] == 0;
}

ICOS_Float func_parameter::set_param(ICOS_Float *a, ICOS_Float value) {
  this->value = a[params[0].index] = value;
  return value;
}

void func_parameter::evaluate( ICOS_Float x, ICOS_Float *a ){
  value = a[index];
}

/**
 * Does nothing. Partial is always 1.0
 */
void func_parameter::evaluate_partials() {
}

void func_parameter::output_params(FILE *ofp, bool fixed) {
  if (fixed) {
    fprintf(ofp, " %d", param_fixed() ? 0 : 1);
  } else {
    fprintf(ofp, " %13.7" FMT_E, value);
  }
}

func_parameter_p Nparameter(const char *group) {
  char pname[40];
  snprintf(pname, 40, "N[%s]", group);
  return new func_parameter(nl_strdup(pname), 0);
}

static ICOS_Float get_molwt( int isotopomer ) {
  switch (isotopomer) {
    case 11: return 18.010565; // H_2O
    case 12: return 20.014811; // H_2{}^{18}O
    case 13: return 19.01478;  // H_2{}^{17}O
    case 14: return 19.01674;  // HDO
    case 15: return 21.020985; // HD^{18}O
    case 16: return 20.020956; // HD^{17}O
    case 21: return 43.989830; // CO2
    case 22: return 44.993183; // CO_2{}^{13}C
    case 23: return 45.994076; // CO_2{}^{18}O
    case 24: return 44.994045; // CO_2{}^{17}O
    case 41: return 44.001060; // N2O
    case 42: return 44.998096; // N^{15}NO
    case 43: return 44.998096; // NN^{15}O
    case 61: return 16.031300; // CH_4
    case 62: return 17.034655; // C13H4
    case 63: return 17.037476; // CH3D
    case 81: return 29.997990; // NO
    case 101: return 45.992905; // NO2
    case 131: return 17.002741; // OH
    case 151: return 35.976677; //H35Cl
    case 152: return 37.973728; //H37Cl
    case 153: return 36.983; // D35Cl
    case 154: return 38.98; // D37Cl
    case 201: return 30.010565; // H2C0
    case 202: return 31.013920; // H213CO
    case 203: return 32.014812; // H2C18O
    default:
      nl_error( 3,
        "Uncatalogued isotopomer '%d': Edit funceval.c get_molwt()",
         isotopomer );
      return 0; // Never reached
  }
}


/**
  func_line object has at least 3 parameters
    dnu_idx = 0: Fine location in cm-1
    w_idx = 1: Doppler E-folding halfwidth in cm-1
    n_idx = 2: Number Density in molecules/cm3
  These parameters are all initialized to zero.
 
  There are also a raft of fixed parameters based on HITRAN
  data for individual lines:
  molecule
  isotope
  nu = wavenumber
  S = spectral line intensity
  G_air = Gamma_air
  E = Lower state energy of the transition
  n = coefficient of temerature dependence or air-broadened hw
  delta = air-broadened pressure shift
  ipos = sample number of line center at start
 */
const int func_line::dnu_idx = 0;
const int func_line::w_idx = 1;
const int func_line::n_idx = 2;
const double func_line::DRTPI = 0.5641895835477563; // 1/SQRT(pi)
const double func_line::Tref = 296.; // K
const double func_line::C2 = 1.4388; // cm K second radiation constant hc/k
double func_line::nu0 = 0.;
int func_line::n_lines = 0;

func_line::func_line( const char *name, int mol, int iso,
          double nu_in, double S_in, double G_air_in, double E_in,
          double n_in, double delta_in, unsigned int ipos_in, double threshold,
          int fix_w, int fix_fp, func_parameter *N ) :
    func_evaluator(name, true, ++n_lines) {
  line_number = n_lines;
  append_func(new func_parameter("dnu", 0., true, line_number));
  append_func(new func_parameter("gd", 1., true, line_number));
  append_func(N ? N : new func_parameter("N", 0., true, line_number));
  fixed = 0;
  fix_finepos = fix_fp;
  fix_width = fix_w;
  prev_numdens = 0.;
  prev_ged = 1.;
  rolledback = 0;
  isotopomer = mol*10 + iso;
  QT = 0;
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
  nu_F0_idx = -1;
}

func_line::~func_line() {
  delete(QT);
}

unsigned int func_line::adopted(func_evaluator *new_parent) {
  append_func(new_parent->args[0].arg);
  nu_F0_idx = args.size()-1;
  return func_evaluator::adopted(new_parent);
}

void func_line::init(ICOS_Float *a) {
  func_evaluator::init(a);
  if (fix_width) fix_param( w_idx );
  if (QT == 0)
    QT = new QTdata(isotopomer);
}

void func_line::print_config( FILE *fp ) {
  fprintf( fp,
    "  %d %d %.6lf %.4" FMT_E
    " %.4" FMT_F " %.4" FMT_F
    " %.2" FMT_F " %.6" FMT_F " %lu\n",
    isotopomer/10, isotopomer%10, nu1+nu0, S, G_air, E,
    n_air, delta, params.size() );
}

void func_line::print_intermediates(FILE *fp) {}

/**
 * @param alamda The Levenberg-Marquardt lambda parameter
 * @param P Current pressure in Torr
 * @param T Current temperature in Kelvin
 * @param a Pointer to paramater values vector
 *
 * According to the description of the Levenberg-Marquardt Method
 * in Numerical Recipes in C, alamda takes on a few special values.
 * Values less than zero indicate initialization. alamda is also set
 * to zero for the final computation of the covariance matrix.
 * Otherwise it takes values greater than zero. In adjust_params(),
 * we use alamda primarily to identify initialization steps. In
 * particular, alamda == -2 on the first initialization of each
 * fit and then set to -1 at the beginning of each iteration
 * step.
 */
int func_line::adjust_params(ICOS_Float alamda,
          ICOS_Float P, ICOS_Float T, ICOS_Float *a) {
  // Eliminated a check for drifting. Taken care of in func_abs.
  if ( alamda < -1.5 ) { // very first initialization
    double Spt = S * QT->evaluate(T) * exp(-C2*E/T) * (1-exp(-C2*nu/T))
            * Corr_Tref;
    Ks = Spt * GlobalData.CavityLength * DRTPI; 
    nu_P = nu1 + delta * P/760.;
    rolledback = 0;
  }
  ICOS_Float numdens = get_arg( a, n_idx );
  // Negative number densities, although physically nonsensical,
  // are important for statistical purposes. If we arbitrarily
  // force the fit to return values >= 0, that pushes the mean
  // values above zero even when no absorption is present.
  ICOS_Float gamma_ed;
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
      gamma_ed = get_arg(a, w_idx);
    }
  } else {
    const ICOS_Float min_ged = 1e-3;
    gamma_ed = get_arg(a, w_idx);
    if ( gamma_ed <= min_ged )
      gamma_ed = set_param(a, w_idx, (prev_ged-min_ged)/2 + min_ged );
    prev_ged = gamma_ed;
  }
  if ( fixed ) {
    if ( alamda == 0 ) {
      ICOS_Float strength = gamma_ed > 0 ? Ks * numdens / gamma_ed : 0.;
      if ( strength > S_thresh * 4. ) {
        if ( rolledback < 2 ) {
          nl_error( 0, "Floating line %d (strength %" FMT_G ")",
                      line_number, strength );
          line_float();
          return 1;
        } else nl_error( 0, "NOT re-floating line %d",
                  line_number );
      }
    }
  } else {
    ICOS_Float strength = gamma_ed > 0 ? Ks * numdens / gamma_ed : 0.;
    if ( strength <= S_thresh ) {
      nl_error( 0, "Fixing line %d (strength %" FMT_G ")",
                    line_number, strength );
      line_fix();
      if (alamda >= 0 ) rolledback++;
      return 1;
    }
  }
  if ( ! param_fixed(dnu_idx)) {
    ICOS_Float dnu = get_arg(a, dnu_idx);
    if ( fabs(dnu) > GlobalData.TolerableDrift ) {
      nl_error( 0,
        "Increasing threshold for wandering line at %.4" FMT_F
        ", dnu = %.4" FMT_F, nu, dnu );
      line_fix();
      S_thresh =
        Ks * get_arg( a, n_idx ) /
          get_arg(a, w_idx );
      return 1;
    }
  }
  return 0;
}

func_line *func_line::is_line() { return this; }

ICOS_Float func_line::line_start(ICOS_Float *a) {
  return (nu_P - GlobalData.RightLineMarginMultiplier*get_arg(a, w_idx));
}
ICOS_Float func_line::line_end(ICOS_Float *a) {
  return (nu_P + GlobalData.LeftLineMarginMultiplier*get_arg(a, w_idx));
}

int func_evaluator::line_check(int include, ICOS_Float& start, ICOS_Float& end,
                ICOS_Float P, ICOS_Float T, ICOS_Float *a) {

  std::vector<argref>::iterator arg;
  for (arg = args.begin(); arg != args.end(); ++arg) {
    if (arg->arg->line_check(include, start, end, P, T, a))
      return 1;
  }
  return 0;
}

int func_evaluator::skew_samples() {
  int rv = 0;
  int frv;
  std::vector<argref>::iterator arg;
  for (arg = args.begin(); arg != args.end(); ++arg) {
    frv = arg->arg->skew_samples();
    if ( frv > rv ) rv = frv;
  }
  return rv;
}

/**
  dump_params() is a diagnostic tool that is invoked
  at times of failure. The default version simply lists
  all the parameters and their values without recursing
  to children, but overrides can delegate to children.
 */
void func_evaluator::dump_params() {
  // print_indent( stderr, indent );
  const char *comma = "";
  int len = fprintf( stderr, "  %s(", name );
  std::vector<argref>::iterator arg;
  for (arg = args.begin(); arg != args.end(); ++arg) {
    if (len+strlen(arg->arg->name) > 50) {
      fprintf(stderr, "%s\n      ", comma);
      len = 6;
      comma = "";
    }
    len += fprintf(stderr, "%s%s", comma, arg->arg->name);
    comma = ",";
  }
  fprintf(stderr, ") = %" FMT_G "\n", value);
  for (arg = args.begin(); arg != args.end(); ++arg) {
    fprintf(stderr, "    d(%s)/d(%s) = %" FMT_G "\n",
      name, arg->arg->name, arg->dyda);
    comma = ",";
  }
}

void func_evaluator::dump_partials() {
  fprintf(stderr, "    Partials with respect to parameters for %s():\n", name);
  std::vector<parameter>::iterator param;
  for (param = params.begin(); param != params.end(); ++param) {
    fprintf(stderr, "      /d[%d] (%s) = %" FMT_G "\n", param->index,
      param->name, param->dyda);
  }
}

void func_parameter::dump_params() {
  // print_indent( stderr, indent );
  fprintf( stderr, "  %s (a[%d]) = %" FMT_G "\n", name, params[0].index, value);
}

void func_evaluator::print_indent( FILE *fp, int indent ) {
  while (indent-- > 0) fputc(' ', fp );
}

void func_evaluator::print_config(FILE *fp) {}
void func_evaluator::print_intermediates(FILE *fp) {}

/**
 * line_check(include, start, end, P, T, a );
 * operates in three passes. First, include is set to 0 to
 * indicate the 'exclude' step. A line is excluded if it
 * hits the start or end boundary, and if so, the boundaries
 * are moved in to also exclude any significant portion
 * of this line. If the boundaries need to be moved, we
 * return 1, and another exclude pass is indicated. This
 * is necessary to deal with overlapping lines.
 *
 * Refinement to the exclusion: If a line hits the boundary,
 * we first try to re-evaluate it's width by fixing the line
 * and calling adjust_params. If that helps, then its threshold
 * needs to be raised. If it doesn't help, then turn it off.
 *
 * Next, include is set to 1. Any lines that were previously
 * excluded that now fall within the sample range are re-enabled.
 *
 * Once a final set of lines has been determined, include
 * is set to 2 to indicate the 'include' step. Here the
 * boundaries are expanded to include all the lines which
 * are still 'on'.
 */
int func_line::line_check(int include, ICOS_Float& start, ICOS_Float& end,
                    ICOS_Float P, ICOS_Float T, ICOS_Float *a ) {
  ICOS_Float ls = line_start(a);
  ICOS_Float le = line_end(a);
  // func_line *next = lnext();
  if ( include == 0 ) {
    int rv = -1;
    if ( ! fixed && ( ls < start || le > end ) ) {
      ICOS_Float save_thresh = S_thresh;
      line_fix();
      S_thresh = Ks * get_arg(a, n_idx)*2/get_arg(a, w_idx);
      adjust_params( -1, P, T, a );
      ls = line_start(a);
      le = line_end(a);
      if ( ls >= start && le <= end ) {
        nl_error(0, "Raised threshold on line %d near boundary",
                      line_number );
      } else {
        nl_error( 0, "Fixing line %d (%.4" FMT_F ",%.4" FMT_F ")",
                          line_number, ls, le );
        S_thresh = save_thresh;
      }
    }
    if ( le < start || ls > end ) {
      ICOS_Float lem = le+GlobalData.LeftLineMargin;
      ICOS_Float lsm = ls-GlobalData.RightLineMargin;
      rv = 0;
      if ( ls < start && lem - GlobalData.LineMarginHysteresis > start ) {
        start = lem; rv = 1;
        if ( GlobalData.Verbosity & 2 )
          nl_error( 0, "Exclude: Updated start to %.4" FMT_F, start );
      }
      if ( le > end && lsm + GlobalData.LineMarginHysteresis < end ) {
        end = lsm; rv = 1;
        if ( GlobalData.Verbosity & 2 )
          nl_error( 0, "Exclude: Updated end to %.4" FMT_F, end );
      }
      if ( ! param_fixed(n_idx) ) {
        nl_error( 0, "Turning off line %d (%.4" FMT_F ",%.4" FMT_F ")",
                          line_number, ls, le );
        fix_param(n_idx);
        set_param(a, n_idx, 0.); // ### This will need to change with shared N
        line_fix();
      }
      if ( rv != 0 ) return rv;
    }
  }
  if (include == 1) {
    if ( ls >= start && le <= end && param_fixed(n_idx) ) {
      nl_error( 0, "Turning on line %d (%.4" FMT_F ",%.4" FMT_F ")",
                          line_number, ls, le );
      float_param(n_idx);
      // We don't actually line_float() until the fit raises the
      // number density high enough.
    }
  }
  if (include == 2) {
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
  fix_param(dnu_idx);
  fix_param(nu_F0_idx);
  fix_param(w_idx);
  fixed = 1;
}

void func_line::line_float() {
  float_param(nu_F0_idx);
  if ( fix_finepos == 0 ) {
    float_param(dnu_idx);
  }
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
void gaussian::evaluate(ICOS_Float x, ICOS_Float *a) {
  static const ICOS_Float four_log_2 = 4 * log(2.);
  static const ICOS_Float fl2_pi = sqrt(four_log_2/M_PI);
  ICOS_Float dnu = a[params[dnu_idx].index] +
    a[params[nu_F0_idx].index];
  ICOS_Float w = a[params[w_idx].index];
  ICOS_Float s = a[params[n_idx].index];
  assert( s >= 0. );
  ICOS_Float w2 = w * w;
  ICOS_Float v1 = fl2_pi * s / w;
  int xx = int(x);
  ICOS_Float diff = ICOSfile::wndata->data[xx] - nu_P + dnu;
  ICOS_Float diff2 = diff * diff;
  ICOS_Float v2 = exp( - diff2 * four_log_2 / w2 );
  value = v1 * v2;
  args[n_idx].dyda = fl2_pi * v2 / w;
  args[dnu_idx].dyda = args[nu_F0_idx].dyda =
    2*v1*v2*diff*four_log_2/w2;
  ICOS_Float v3 = fl2_pi * s * v2 / w2;
  args[w_idx].dyda = -(v3/w) + v3 * 2 * diff2 * four_log_2 / w2;
}

//---------------------------------------------------------
// lorentzian
//---------------------------------------------------------
// ### This needs to be fixed for wavenumber scale
void lorentzian::evaluate(ICOS_Float x, ICOS_Float *a) {
  ICOS_Float dnu = a[params[dnu_idx].index] +
    a[params[nu_F0_idx].index];
  ICOS_Float w = a[params[w_idx].index];
  ICOS_Float s = a[params[n_idx].index];
  int xx = int(x);
  ICOS_Float diff = ICOSfile::wndata->data[xx] - nu_P + dnu;
  ICOS_Float v1 = M_PI * ( 4 * diff *diff + w * w );
  ICOS_Float v12 = v1 * v1;
  value = 2 * w * s / v1;
  args[n_idx].dyda = 2 * w / v1;
  args[dnu_idx].dyda = args[nu_F0_idx].dyda =
    16*M_PI*w*s*diff / v12;
  args[w_idx].dyda = (2*s*v1 - 4*M_PI*s*w*w)/v12;
}

//----------------------------------------------------------
// voigt: in humdev.cc
//----------------------------------------------------------
