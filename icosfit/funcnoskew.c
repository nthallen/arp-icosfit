/* funcnoskew.c */
#include <math.h>
#include <stdlib.h>
#include "funceval.h"
#include "global.h"
#include "nortlib.h"

// func_noskew calculates absorption for a simple multi-pass
// cell. Its two children define the input power curve (base)
// and the intra-cavity absorption (abs).
// base and abs have mostly independent parameters, but if
// the base function uses nu_F0, that parameter is shared.
// value = base * exp(- NPasses * abs)
// d/dbase = exp(-NPasses*abs)
// d/dabs = base * (-NPasses)*exp(- NPasses * abs)
func_noskew::func_noskew(func_base *base, func_abs *abs) :
    func_evaluator("noskew") {
  append_func(base);
  append_func(abs);
  basep = base;
  absp = abs;
  if ( basep->uses_nu_F0 != 0 && basep->uses_nu_F0 != 1)
    nl_error(4,"uses_nu_F0 must be 0 or 1");
  // n_base_params = base->n_params - basep->uses_nu_F0;
  // n_abs_params = abs->n_params;
  // n_params = n_base_params + n_abs_params;
  N_Passes = GlobalData.N_Passes +
    GlobalData.CavityFixedLength/GlobalData.CavityLength;
}

// Identical to func_skew::init()
// nu_F0 comes after n_base_params in our list
// If base->uses_nu_F0, it's first param is linked
// to abs' first, which is our n_base_params.
// void func_noskew::init(ICOS_Float *a) {
  // func_evaluator *child;
  // int p1 = 0;
  // int p2;

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

void func_noskew::evaluate(ICOS_Float x, ICOS_Float *a) {
  // ### int i;
  // func_evaluator::evaluate( x, a ); // evaluate base and abs
  ICOS_Float P = basep->value;
  ICOS_Float Abs = absp->value;
  if ( isnan(P) )
    nl_error(2,"Base(%.0lf) is NaN", x);
  if ( isnan(Abs) )
    nl_error(2,"Absorb(%.0lf) is NaN", x);
  ICOS_Float eNabs = exp( -N_Passes * Abs );
  value = P * eNabs;
  if ( isnan(value) )
    nl_error(2,"noskew(%.0lf) is NaN", x);
  ICOS_Float NPeNabs = - N_Passes * value;
  args[0].dyda = eNabs;
  args[1].dyda = P * NPeNabs;

  // for ( i = 0; i < n_base_params; i++ )
    // params[i].dyda = eNabs *
	// basep->params[i+basep->uses_nu_F0].dyda;
  // for ( i = 0; i < n_abs_params; i++ )
    // params[i+n_base_params].dyda = NPeNabs * absp->params[i].dyda;
  // if ( basep->uses_nu_F0 )
    // params[n_base_params].dyda += eNabs * basep->params[0].dyda;
}

// void func_noskew::dump_params(ICOS_Float *a, int indent) {
  // print_indent( stderr, indent );
  // fprintf( stderr, "Parameters for '%s':\n", name );
  // indent += 2;
  // basep->dump_params( a, indent );
  // absp->dump_params( a, indent );
// }
