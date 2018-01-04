#include <math.h>
#include <stdlib.h>
#include "funceval.h"
#include "ptread.h"
#include "global.h"
#include "nortlib.h"

func_beta::func_beta(func_abs *abs)
    : func_evaluator("beta") {
  append_func(abs);
}

void func_beta::evaluate(ICOS_Float x, ICOS_Float *a) {
  value = exp(-args[0].arg->value);
  args[0].dyda = -value;
  if ( isnan(value) ) {
    nl_error(2,"Absorb(%.0lf) is NaN", x);
  }
}

func_gamma::func_gamma(func_beta *beta)
    : func_evaluator("gamma") {
  ICOS_Float R = 1 - GlobalData.MirrorLoss;
  R2 = R*R;
  append_func(beta);
}

void func_gamma::evaluate(ICOS_Float x, ICOS_Float *a) {
  value = R2 * args[0].arg->value * args[0].arg->value;
  args[0].dyda = 2*R2*args[0].arg->value;
}

func_epsilon::func_epsilon(func_gamma *gamma)
    : func_evaluator("epsilon") {
  append_func(gamma);
  const ICOS_Float c = 2.99792458e10; // cm/s
  N = c/(2*GlobalData.CavityLength*GlobalData.SampleRate);
}

void func_epsilon::evaluate(ICOS_Float x, ICOS_Float *a) {
  value = pow(args[0].arg->value, N);
  args[0].dyda = N*pow(value, N-1);
}

func_delta::func_delta(func_epsilon *epsilon, func_gamma *gamma)
    : func_evaluator("delta") {
  append_func(epsilon);
  append_func(gamma);
}

void func_delta::evaluate(ICOS_Float x, ICOS_Float *a) {
  ICOS_Float gm1 = 1/(args[1].arg->value - 1);
  value = (args[0].arg->value - 1)*gm1;
  args[0].dyda = 1*gm1;
  args[1].dyda = (1-args[0].arg->value) * gm1 * gm1;
}

func_g::func_g(func_base *base, func_beta *beta, func_delta *delta)
    : func_evaluator("g0") {
  append_func(base);
  append_func(beta);
  append_func(delta);
} 
     
void func_g::evaluate(ICOS_Float x, ICOS_Float *a) {
  ICOS_Float P = args[0].arg->value;
  ICOS_Float beta = args[1].arg->value;
  ICOS_Float delta = args[2].arg->value;
  if ( isnan(P) ) {
    nl_error(2,"Base(%.0lf) is NaN", x);
  }
  value = P*beta*delta;
  args[0].dyda = beta*delta;
  args[1].dyda = P*delta;
  args[2].dyda = P*beta;
}

//-----------------------------------------------
// To recast func_skew to accept a base function
// that depends on nu_F0, we should not derive
// from func_aggregate, since that 
// assumes all params are independent.
//-----------------------------------------------

skew_data::skew_data() {
  dg = deps = 0;
  g = eps = 0.;
  Power = 0.;
  initialized = 0;
  n = 0;
}

void skew_data::set_n_params(int n_gp, int n_epsp ) {
  dg = new ICOS_Float[n_gp];
  deps = new ICOS_Float[n_epsp];
}

func_skew::func_skew( func_g *g, func_epsilon *eps ) :
    func_evaluator("skew") {
  append_func(g);
  append_func(eps);
  basep = g->args[0].arg;
  
  // if ( basep->uses_nu_F0 != 0 && basep->uses_nu_F0 != 1)
    // nl_error(4,"uses_nu_F0 must be 0 or 1");
  // n_base_params = base->n_params - basep->uses_nu_F0;
  // n_abs_params = abs->n_params;
  // n_params = n_base_params + n_abs_params;
  // initialize the structures here?
  const ICOS_Float c = 2.99792458e10; // cm/s
  N = c/(2*GlobalData.CavityLength*GlobalData.SampleRate);
  ICOS_Float R = 1 - GlobalData.MirrorLoss;
  R2 = R*R;
  R2N = pow(R2,N);
  P_scale = (1-R2N)/(1-R2);
  M = (int) ceil(log(GlobalData.SkewTolerance)/(2*N*log(R)));
  skew = new skew_data[M];
  skew_eval_order.set_children(this);
  pre_evaluation_order.add(this);
}

// Assign the first parameters to base's parameters
// Then if base->uses_nu_F0, link abs's first parameter
// to our first (which is base's first, which is nu_F0)
void func_skew::init(ICOS_Float *a) {
  func_evaluator::init(a);
  for (int i = 0; i < M; i++ )
    skew[i].set_n_params( args[0].arg->params.size(),
      args[1].arg->params.size());
}
  // func_evaluator *child;
  // int p1 = 0;
  // int p2;

  // // printf( "func_skew::init(%s, %d); n_params=%d\n", name, p1, n_params );
  // for ( p2 = 0; p2 < n_params; p2++ )
    // a[params[p2].index] = params[p2].init;
  // for ( child = first; child != 0; child = child->next ) {
    // if ( child->params == 0 )
      // child->params = new parameter[child->n_params];
    // if ( child == first && basep->uses_nu_F0 ) {
      // link_param( n_base_params, child, 0 );
      // p2 = 1;
    // } else {
      // p2 = 0;
    // }
    // if ( p1 + child->n_params - p2 > n_params ) {
      // fprintf( stderr, "Too many child params: n_params = %d\n", n_params );
      // exit(1);
    // }
    // for ( ; p2 < child->n_params; p2++ ) {
      // link_param( p1, child, p2 );
      // p1++;
    // }
    // child->init(a);
  // }
// }

int func_skew::skew_samples() {
  int rv = func_evaluator::skew_samples();
  if ( M > rv ) rv = M;
  return rv;
}

void func_skew::pre_eval(ICOS_Float x, ICOS_Float *a) {
  if (GlobalData.PTE_MirrorLoss_col) {
    ICOS_Float R = 1 - GlobalData.input.MirrorLoss;
    R2 = R*R;
    R2N = pow(R2,N);
    P_scale = (1-R2N)/(1-R2);
  }
  for ( int i = 0; i < M; i++ )
    skew[i].initialized = 0;
  skewidx = 0;
  for (ICOS_Float xi = x - M + 1.; xi < x; ++xi) {
    skew_eval_order.evaluate(xi, a);
    sub_eval(xi, a);
  }
}

/**
 * Reinitializes skew[skewidx] with n=0 and iterates
 * all the other elements that have been initialized
 */
void func_skew::sub_eval(ICOS_Float x, ICOS_Float *a) {
  int i;
  func_evaluator *arg;
  for (i = 0; i < M; ++i) {
    skew_data *curskew = &skew[i];
    if (i == skewidx) {
      curskew->n = 0;
      arg = args[0].arg;
      curskew->g = arg->value;
      // base is arg[0] of func_g, which is arg[0] of
      // func_skew
      curskew->Power = basep->value;
      for (unsigned j = 0; j < arg->params.size(); ++j) {
        curskew->dg[j] = arg->params[j].dyda;
      }
      arg = args[1].arg;
      curskew->eps = arg->value;
      for (unsigned j = 0; j < arg->params.size(); ++j) {
        curskew->deps[j] = arg->params[j].dyda * curskew->g;
      }
      curskew->initialized = 1;
    } else if (skew[i].initialized) {
      arg = args[0].arg;
      curskew->g *= curskew->eps;
      curskew->Power *= R2N;
      for (unsigned j = 0; j < arg->params.size(); ++j) {
        curskew->dg[j] *= curskew->eps;
      }
      if (++curskew->n > 1) {
        arg = args[1].arg;
        for (unsigned j = 0; j < arg->params.size(); ++j) {
          curskew->deps[j] *= curskew->eps;
        }
      }
    }
  }
  if (++skewidx == M)
    skewidx = 0;
}

// We depend on the fact that the functions are evaluated
// for monotonically increasing values of x. We know we are
// starting over if x decreases. We also know that x is
// essentially integral.
// Note: The calculation of 'Power' is not part of the fit
// per se. It is performed only to generate the effective
// skewed baseline for diagnostic output. The baseline
// function determines the input power, but for useful
// comparison with the raw data, we need to estimate the
// output power in the absence of absorption.
void func_skew::evaluate(ICOS_Float x, ICOS_Float *a) {
  sub_eval(x, a);
  value = 0;
  ICOS_Float Pwr = 0;
  for (int i = 0; i < M; i++) {
    skew_data *cs = &skew[i];
    value += cs->g;
    Pwr += cs->Power;
  }
  basep->value = Pwr * P_scale; // For diagnostic output
}

void func_skew::evaluate_partials() {
  std::vector<parameter>::iterator p;
  for (p = params.begin(); p != params.end(); ++p) {
    p->dyda = 0;
  }
  for (int i = 0; i < M; i++) {
    skew_data *cs = &skew[i];
    for (p = params.begin(); p != params.end(); ++p) {
      std::vector<paramref>::iterator ref;
      for (ref = p->refs.begin(); ref != p->refs.end(); ++ref) {
        if (ref->arg_num == 0) {
          p->dyda += cs->dg[ref->param_num];
        } else {
          p->dyda += cs->n*cs->deps[ref->param_num];
        }
      }
    }
  }
}

// void func_skew::dump_params(ICOS_Float *a, int indent) {
  // print_indent( stderr, indent );
  // fprintf( stderr, "Parameters for '%s':\n", name );
  // indent += 2;
  // basep->dump_params( a, indent );
  // absp->dump_params( a, indent );
// }
