#include <math.h>
#include "ICOSfit.h"
#include "global.h"

func_abs::func_abs() : func_evaluator("abs") {
  append_func(new func_parameter("nu_F0", 0.));
}

void func_abs::print_config(FILE *fp) {
  std::vector<func_evaluator>::iterator *child;
  func_line *line;
  
  fprintf( fp, "CavLen = %.1" FMT_F ";\n", GlobalData.CavityLength );
  fprintf( fp, "n_abs_params = 1;\nn_abs_line_params = 1;\n" );
  fprintf( fp, "lines = [\n" );
  // Skip the first arg, which is nu_F0
  for (child = ++args.begin(); child != args.end(); ++child) {
    (*child)->print_config( fp );
  }
  fprintf( fp, "];\n" );
}

void func_abs::print_intermediates(FILE *fp) {
  func_line *child;
  for ( child = lfirst(); child != 0; child = child->lnext() ) {
    child->print_intermediates(fp);
  }
}

// static ICOS_Float interpolate( ICOS_Float x0, ICOS_Float x1, ICOS_Float y0, ICOS_Float y1, ICOS_Float x ) {
//   return (x-x0)*(y1-y0)/(x1-x0) + y0;
// }

// Assumes vals is increasing strictly monotonically.
// Returns the interpolated or extrapolated index of the value
// ------Currently unused-------------------------------------
// static ICOS_Float binary_lookup( ICOS_Float *vals, int n_vals, ICOS_Float value ) {
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

int func_abs::adjust_params( ICOS_Float alamda, ICOS_Float P, ICOS_Float T, ICOS_Float *a ) {
  // New: Call adjust_params for all children, since samples/cm-1
  // won't be changing.
  // static ICOS_Float max_drift_per_sec = 30; // 8
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
  ICOS_Float nu_F0 = get_param(a,0);
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
      ICOS_Float dnu = get_param( a, p1 );
      child->set_param( a, child->dnu_idx, nu_F0 + dnu );
      if ( ! param_fixed(p1)) {
        if ( fabs(dnu) > GlobalData.TolerableDrift ) {
          nl_error( 0,
            "Increasing threshold for wandering line at %.4" FMT_F
	    ", dnu = %.4" FMT_F, child->nu, dnu );
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

// void func_abs::fix_linepos( int linenum ) {
  // fix_param( linenum*5 - 4 );
// }

// void func_abs::float_linepos( int linenum ) {
  // float_param( linenum*5 - 4 );
// }

// void func_abs::append_func( func_line *newfunc ) {
  // func_evaluator::append_func(newfunc);
  // n_params += 1 + newfunc->n_params;
// }

/**
 * @param a The global parameter value vector.
 * We mark nu_F0 as fixed so that it floats only when
 * there are lines which are strong enough to float.
 */
void func_abs::init(ICOS_Float *a) {
  func_evaluator::init();
  fix_param(0); // nu_F0
}

void func_abs::evaluate(ICOS_Float x, ICOS_Float *a) {
  std::vector<argref>::iterator child;
  int p1, i;
  
  // func_evaluator::evaluate(x, a);
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

void func_abs::dump_params(ICOS_Float *a, int indent) {
  func_evaluator *child;
  int indx, p1;

  print_indent( stderr, indent );
  fprintf( stderr, "Parameters for '%s':\n", name );
  indent += 2;
  indx = params[0].index;
  print_indent( stderr, indent );
  fprintf( stderr, "[%2d] nu_F0: %" FMT_G " + nu0(%lg) = %lg\n",
    indx, a[indx], func_line::nu0, a[indx]+func_line::nu0);
  indent += 2;
  p1 = 1;
  for ( child = first; child != 0; child = child->next ) {
    print_indent( stderr, indent );
    indx = params[p1].index;
    fprintf( stderr, "[%2d] dnu: %" FMT_G "\n", indx, a[indx] );
    child->dump_params( a, indent );
    p1 += 1 + child->n_params;
  }
}
