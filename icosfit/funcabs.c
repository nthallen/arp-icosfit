#include <math.h>
#include "ICOSfit.h"
#include "global.h"

func_abs::func_abs() : func_evaluator("abs") {
  append_func(new func_parameter("nu_F0", 0.));
}

void func_abs::print_config(FILE *fp) {
  std::vector<argref>::iterator child;
  
  fprintf( fp, "CavLen = %.1" FMT_F ";\n", GlobalData.CavityLength );
  fprintf( fp, "n_abs_params = 1;\nn_abs_line_params = 0;\n" );
  fprintf( fp, "lines = [\n" );
  // Skip the first arg, which is nu_F0
  for (child = args.begin(); child != args.end(); ++child) {
    child->arg->print_config( fp );
  }
  fprintf( fp, "];\n" );
}

void func_abs::print_intermediates(FILE *fp) {
  std::vector<argref>::iterator child;
  for (child = args.begin(); child != args.end(); ++child) {
    child->arg->print_intermediates(fp);
  }
}

int func_abs::adjust_params( ICOS_Float alamda, ICOS_Float P, ICOS_Float T, ICOS_Float *a ) {
  int rv = 0;
  std::vector<argref>::iterator child;
  for (child = args.begin(); child != args.end(); ++child) {
    if ( child->arg->adjust_params( alamda, P, T, a ) )
      rv = 1;
  }
  ICOS_Float nu_F0 = get_arg(a,0);
  if (alamda < 0 && GlobalData.input.nu_F0 != 0.) {
    nu_F0 = GlobalData.input.nu_F0;
    set_param( a, 0, nu_F0 );
  }
  // This is true only on the very first fit (if then)
  // ### Can this all be moved to func_line::adjust_params()?
  // ### No, I don't think so, since it is assimilating data from
  // ### all the lines. Why not init?
  if ( nu_F0 == 0. && GlobalData.input.nu_F0 == 0. ) {
    int n_lines = 0;
    for (child = args.begin(); child != args.end(); ++child) {
      func_line *line = child->arg->is_line();
      if (line) {
        if ( line->ipos >= GlobalData.SignalRegion[0] &&
             line->ipos <= GlobalData.SignalRegion[1] ) {
          nu_F0 += line->nu_P - ICOSfile::wndata->data[line->ipos];
          n_lines++;
        }
      }
    }
    if ( n_lines > 0 ) nu_F0 /= n_lines;
    else nl_error( 3, "No valid line starting positions recorded" );
    set_param( a, 0, nu_F0 );
  }
  // Should I ever see the case nu_F0 == 0 and GlobalData.input.nu_F0 != 0.?
  
  if ( rv != 0 ) return rv;
  return rv;
}

/**
 * @param a The global parameter value vector.
 * We mark nu_F0 as fixed so that it floats only when
 * there are lines which are strong enough to float.
 */
void func_abs::init(ICOS_Float *a) {
  func_evaluator::init(a);
  fix_param(0); // nu_F0
}

void func_abs::evaluate(ICOS_Float x, ICOS_Float *a) {
  std::vector<argref>::iterator child;
  
  args[0].dyda = 0.;
  value = 0;
  for (child = ++args.begin(); child != args.end(); ++child) {
    value += child->arg->value;
    child->dyda = 1;
  }
}
