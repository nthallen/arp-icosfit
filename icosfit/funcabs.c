#include <math.h>
#include "ICOSfit.h"
#include "global.h"

void func_abs::print_config( FILE *fp ) {
  func_line *child;
  fprintf( fp, "CavLen = %.1f;\n", GlobalData.CavityLength );
  fprintf( fp, "n_abs_params = 1;\nn_abs_line_params = 1;\n" );
  fprintf( fp, "lines = [\n" );
  for ( child = lfirst(); child != 0; child = child->lnext() ) {
    child->print_config( fp );
  }
  fprintf( fp, "];\n" );
}

// static float interpolate( float x0, float x1, float y0, float y1, float x ) {
//   return (x-x0)*(y1-y0)/(x1-x0) + y0;
// }

// Assumes vals is increasing strictly monotonically.
// Returns the interpolated or extrapolated index of the value
// ------Currently unused-------------------------------------
// static float binary_lookup( float *vals, int n_vals, float value ) {
//   int low = 0, high = n_vals-1;
//   if ( value < vals[0] )
//     return interpolate( vals[0], vals[1], 0., 1., value );
//   if ( value > vals[high] )
//     return interpolate( vals[high-1], vals[high], high-1., high, value );
//   while ( high > low+1 ) {
//     int mid = (high+low)/2;
//     if ( value < vals[mid] )
//       high = mid;
//     else low = mid;
//   }
//   return interpolate( vals[low], vals[high], low, high, value );
// }

int func_abs::adjust_params( float alamda, float P, float T, float *a ) {
  // New: Call adjust_params for all children, since samples/cm-1
  // won't be changing.
  // static float max_drift_per_sec = 30; // 8
  int rv = 0;
  int lines_floating = 0;
  func_line *child;
  for ( child = lfirst(); child != 0; child = child->lnext() ) {
    if ( child->adjust_params( alamda, P, T, a ) )
      rv = 1;
    if ( ! child->fixed ) lines_floating++;
  }
  if ( lines_floating == 0 && ! param_fixed(0) ) {
    nl_error( 0, "Fixing nu_F0" );
    fix_param(0);
    rv = 1;
  }
  if ( alamda == 0 && lines_floating > 0 && param_fixed(0) ) {
    nl_error( 0, "Floating nu_F0" );
    float_param(0);
    rv = 1;
  }
  float nu_F0 = get_param(a,0);
  if (alamda < 0 && GlobalData.input.nu_F0 != 0.) {
    nu_F0 = GlobalData.input.nu_F0;
    set_param( a, 0, nu_F0 );
  }
  // This is true only on the very first fit (if then)
  if ( nu_F0 == 0. && GlobalData.input.nu_F0 == 0. ) {
    int n_lines = 0;
    for ( child = lfirst(); child != 0; child = child->lnext() ) {
      if ( child->ipos >= GlobalData.SignalRegion[0] &&
           child->ipos <= GlobalData.SignalRegion[1] ) {
        nu_F0 += child->nu_P - ICOSfile::wndata->data[child->ipos];
        n_lines++;
      }
    }
    if ( n_lines > 0 ) nu_F0 /= n_lines;
    else nl_error( 3, "No valid line starting positions recorded" );
    set_param( a, 0, nu_F0 );
  }
  // Should I ever see the case nu_F0 == 0 and GlobalData.input.nu_F0 != 0.?
  
  if ( rv != 0 ) return rv;
  // nu_F0 += GlobalData.input.nu_F0; // No longer required
  { int p1 = 1;
    for ( child = lfirst(); child != 0; child = child->lnext() ) {
      float dnu = get_param( a, p1 );
      child->set_param( a, child->l_idx, nu_F0 + dnu );
      if ( ! param_fixed(p1)) {
        if ( fabs(dnu) > GlobalData.TolerableDrift ) {
          nl_error( 0,
            "Increasing threshold for wandering line at %.4f, "
            "dnu = %.4f", child->nu, dnu );
          child->line_fix();
          child->S_thresh =
            child->Ks * child->get_param( a, child->n_idx ) /
              child->get_param(a, child->w_idx );
          // child->adjust_params( alamda, P, T, a );
          rv = 1;
        }
      }
      p1 += 1 + child->n_params;
    }
  }
  return rv;
}

void func_abs::fix_linepos( int linenum ) {
  fix_param( linenum*5 - 4 );
}

void func_abs::float_linepos( int linenum ) {
  float_param( linenum*5 - 4 );
}

void func_abs::append_func( func_line *newfunc ) {
  func_evaluator::append_func(newfunc);
  n_params += 1 + newfunc->n_params;
}

// Too much of this is copied from func_evaluator::init.
// This goes to the general issue of more clearly characterizing
// what a parameter is. The problem is for func_abs we are
// abandoning the assumption of independence.
void func_abs::init(float *a) {
  func_evaluator *child;
  int p1, p2;


  // printf( "func_abs::init(%s); n_params=%d\n", name, n_params );
  // set_param(a, 0, 0.); // nu_F0 (oops, wrong init)
  params[0].init = 0.; // nu_F0
  fix_param(0);
  p1 = 1; // skip nu_F0
  for ( child = first; child != 0; child = child->next ) {
    params[p1].init = 0.;
    fix_param(p1);
    p1++;
    if ( p1 + child->n_params > n_params ) {
      nl_error( 4, "Too many child params: n_params = %d\n", n_params );
    }
    if ( child->params == 0 )
      child->params = new parameter[child->n_params];
    for ( p2 = 0; p2 < child->n_params; p2++ ) {
      link_param( p1, child, p2 );
      p1++;
    }
  }
  for ( p2 = 0; p2 < n_params; p2++ )
    a[params[p2].index] = params[p2].init;
  for ( child = first; child != 0; child = child->next ) {
    child->init(a);
    child->fix_param(0);
  }
}

void func_abs::evaluate(float x, float *a) {
  func_evaluator *child;
  int p1, i;
  
  func_evaluator::evaluate(x, a);
  params[0].dyda = 0.;
  value = 0;
  p1 = 1;
  for ( child = first; child != 0; child = child->next ) {
    value += child->value;
    params[p1++].dyda = child->params[0].dyda;
    params[0].dyda += child->params[0].dyda;
    for ( i = 0; i < child->n_params; i++ ) {
      params[p1++].dyda = child->params[i].dyda;
    }
  }
}

void func_abs::dump_params(float *a, int indent) {
  func_evaluator *child;
  int indx, p1;

  print_indent( stderr, indent );
  fprintf( stderr, "Parameters for '%s':\n", name );
  indent += 2;
  indx = params[0].index;
  print_indent( stderr, indent );
  fprintf( stderr, "[%2d] nu_F0: %g + nu0(%lg) = %lg\n",
    indx, a[indx], func_line::nu0, a[indx]+func_line::nu0);
  indent += 2;
  p1 = 1;
  for ( child = first; child != 0; child = child->next ) {
    print_indent( stderr, indent );
    indx = params[p1].index;
    fprintf( stderr, "[%2d] dnu: %g\n", indx, a[indx] );
    child->dump_params( a, indent );
    p1 += 1 + child->n_params;
  }
}
