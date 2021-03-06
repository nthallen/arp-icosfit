#include "ICOSfit.h"
#include <setjmp.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#define NRANSI
#include "nrutil.h"
#include "nr.h" /* for lfit */
#include "global.h"
#include "clp.h"

static fitdata *crntfit;
jmp_buf Fit_buf;

void BaseFitFunction( ICOS_Float x, ICOS_Float *afunc, int ma ) {
  int i;
  func_evaluator *bfunc = crntfit->base;
  for ( i = 0; i < ma; i++ ) afunc[i+1] = 0.;
  bfunc->evaluate( x, afunc );
  for ( i = 0; i < ma; i++ ) afunc[i+1] = bfunc->params[i].dyda;
}

fitdata::fitdata( PTfile *ptf, ICOSfile *IF,
   func_evaluator *f, func_evaluator *baseline, func_abs *abs ) {
  int i;
  
  PTf = ptf;
  IFile = IF;
  verbose = GlobalData.Verbosity;
  vfp = 0;
  vmlf = 0;
  if ( verbose & 8 ) {
    vmlf = mlf_init( 3, 60, 1, GlobalData.OutputDir, "dir", NULL );
  }
  BaseStart = GlobalData.BackgroundRegion[0];
  BaseEnd = GlobalData.BackgroundRegion[1];
  func = f;
  base = baseline;
  absorb = abs;
  FitBaseline = 1;
  SignalStart = GlobalData.SignalRegion[0] + func->skew_samples();
  SignalEnd = GlobalData.SignalRegion[1];
  
  Start = End = 0;
  npts = npts_vec = 0;
  x = y = sig = 0;
  ma = func_parameter::n_parameters;
  a    = vector(1,ma);
  a_save = vector(1,ma);
  atry = vector(1,ma);
  beta = vector(1,ma);
  da = vector(1,ma);
  dyda = vector(1, ma);
  ia   = ivector(1,ma);
  if ( a == 0 || ia == 0 )
    nl_error( 3, "Out of memory in fitdata::fitdata" );
  for ( i = 1; i <= ma; i++ ) ia[i] = 1;
  func_evaluator::global_evaluation_order.set(f);
  func_evaluator::dump_evaluation_order.set_pre_order(f);
  f->init( a, ia );
  covar = matrix( 1, ma, 1, ma );
  alpha = matrix( 1, ma, 1, ma );
  mfit = 0;
  mf_size = 0;
}

#define RESTART_BUFSIZE 4096

/**
  Resets internal parameters to match a specific start condition
  based on a previous run.
 */
void fitdata::handle_restart( const char *ofname ) {
  // unsigned int prev_ScanNum = 0;
  if ( RestartAt != NoKey)
    GlobalData.RestartAt = GetClpValue(RestartAt, 0);;
  if ( GlobalData.RestartAt > 0 ) {
    char bakname[PATH_MAX];
    FILE *ifp = fopen( ofname, "r" );
    if ( ifp == 0 )
      nl_error( 3, "Unable to read output file %s for restart", ofname );

    for (;;) {
      if ( PTf->readline() ) {
        if (PTf->ScanNum == GlobalData.RestartAt) {
          PTf->backup();
          break;
        // } else {
          // prev_ScanNum = PTf->ScanNum;
        }
      } else {
        nl_error( 3, "Did not find ScanNum %d in PTE file", GlobalData.RestartAt );
      }
    }
    
    if ( GlobalData.PreserveOutput == 1 ) {
      snprintf( bakname, PATH_MAX-1, "%s.%d", ofname, GlobalData.RestartAt );
      bakname[PATH_MAX-1] = '\0';
      IFile->ofp = fopen( bakname, "w" );
      if ( IFile->ofp == 0 )
        nl_error( 3, "Unable to open output file %s", bakname );
    } else {
      snprintf( bakname, PATH_MAX-1, "%s.bak", ofname );
      bakname[PATH_MAX-1] = '\0';
      unlink( bakname );
      rename( ofname, bakname );
      IFile->ofp = fopen( ofname, "w" );
      if ( IFile->ofp == 0 )
        nl_error( 3, "Unable to reopen output file %s", ofname );
    }

    { char buf[RESTART_BUFSIZE];
      unsigned int ScanNum = 0;
      while ( fgets( buf, RESTART_BUFSIZE, ifp ) != 0 ) {
        int col;
        char *p = buf;
        while ( isspace( *p ) ) p++;
        for ( col = 1; col < ScanNum_col; col++ ) {
          while ( ! isspace(*p) && *p != '\0' ) p++;
          while ( isspace(*p) ) p++;
        }
        if ( *p == '\0' ) break;
        ScanNum = strtoul( p, &p, 10 );
        if (ScanNum >= GlobalData.RestartAt) { // What if ScanNum is decreasing?
          break;
        } else {
          for ( ++col; col <= n_input_params; ++col ) {
            while ( isspace(*p) ) p++;
            while ( ! isspace(*p) && *p != '\0' ) p++;
          }
          // double nu_F0d = strtod( p, &p );
          // ICOSfile::nu_F0 = nu_F0d - func_line::nu0;
          // ICOSfile::dFN = (int) strtoul( p, &p, 10 );
          int i;
          for ( i = 1; i <= ma; ++i ) {
            a[i] = strtod( p, &p );
            if ( ! isspace(*p) ) break;
          }
          if ( i <= ma )
            nl_error(3, "Syntax error reading Scan %u from %s", ScanNum, ofname);
          for ( i = 1; i <= ma && *p != '\0'; i++ ) {
            while ( isspace(*p) ) p++;
            if ( *p == '0' || *p == '1' ) ia[i] = *p++ - '0';
            else nl_error( 3, "Expected 0 or 1 during Restart" );
          }
          assert(i > ma);
          fprintf( IFile->ofp, "%s", buf );
        }
      }
      fclose( ifp );
      IFile->read( ScanNum ); // To initialize wndata
      { std::vector<argref>::iterator arg;
        for (arg = absorb->args.begin(); arg != absorb->args.end(); ++arg) {
          func_line *line = arg->arg->is_line();
          if ( line && line->param_fixed( line->dnu_idx ) ) {
            line->fixed = 1;
            if ( line->param_fixed( line->n_idx ) )
              nl_error( 0, "Line at %.4" FMT_F " is off", line->nu );
            else
              nl_error( 0, "Line at %.4" FMT_F " is fixed", line->nu );
          }
        }
      }
      return;
    }
    nl_error( 3, "Unexpected loop exit in handle_restart");
  } else {
    IFile->ofp = fopen( ofname, "a" );
    if ( IFile->ofp == 0 )
      nl_error( 3, "Unable to open output file %s", ofname );
  }
}

static FILE *pathopen( const char *path, const char *format, int fileno ) {
  static const int POBUFSIZE = 256;
  char buf[POBUFSIZE];
  int n = snprintf( buf, POBUFSIZE-1, "%s/", path );
  if ( n < 0 || n >= POBUFSIZE-1 ) return 0;
  n = snprintf( buf+n, POBUFSIZE-n-1, format, fileno );
  if ( n < 0 || n >= POBUFSIZE-1 ) return 0;
  FILE *fp = fopen( buf, "w" );
  return fp;
}

int fitdata::adjust_params( ICOS_Float *av ) {
  return func_evaluator::global_evaluation_order.adjust_params(alamda,
    PTf->P, PTf->T, av );
}

void print_matrix( ICOS_Float **mat, const char *name, int nrow, int ncol ) {
  int row, col;
  for ( row = 1; row <= nrow; row++ ) {
    fprintf( stderr, "%s[%d][] = ", name, row );
    int cols = 0;
    for ( col = 1; col <= ncol; col++ ) {
      if ( cols > 500 ) {
        fprintf( stderr, "\n  " );
        cols = 0;
      }
      fprintf( stderr, "  %12.4" FMT_G, mat[row][col] );
      cols += 14;
    }
    fprintf( stderr, "\n" );
  }
}

void print_vector( ICOS_Float *vec, const char *name, int ncol ) {
  int col;
  fprintf( stderr, "%s[] = ", name );
  int cols = 0;
  for ( col = 1; col <= ncol; col++ ) {
    if ( cols > 500 ) {
      fprintf( stderr, "\n  " );
      cols = 0;
    }
    fprintf( stderr, "\t%12.4" FMT_G, vec[col] );
    cols += 14;
  }
  fprintf( stderr, "\n" );
}

/**
  Handles the complete fit for one raw scan file.
 */
int fitdata::fit( ) {
  ICOS_Float *yin = IFile->sdata->data;
  int i;

  // store the object for reference from lower funcs
  // This is solely for the baseline fit, and hence
  // should go away. The linear fit routine should also
  // be objectified.
  crntfit = this;

  alamda=-2;
  while ( adjust_params( a ) != 0 ) alamda = -1;

  // Evaluation of range of the input over which we should
  // fit. This should not be limited to the absorb func.
  // It should be a virtual function of func_evaluator
  Start = End = 0;
  {
    // Note: Wavenumber decreases with sample number
    ICOS_Float nu_F0 = absorb->get_arg( a, 0 );
    ICOS_Float wnStart = IFile->wndata->data[SignalEnd] + nu_F0;
    ICOS_Float wnEnd = IFile->wndata->data[SignalStart] + nu_F0;
    while ( func->line_check( 0, wnStart, wnEnd, PTf->P, PTf->T, a ) != 0 );
    func->line_check( 1, wnStart, wnEnd, PTf->P, PTf->T, a );
    ICOS_Float EwnStart = 0., EwnEnd = 0.;
    func->line_check( 2, EwnStart, EwnEnd, PTf->P, PTf->T, a );
    if ( EwnStart != 0. && EwnEnd != 0. ) {
      if ( EwnStart > wnStart ) wnStart = EwnStart;
      if ( EwnEnd < wnEnd ) wnEnd = EwnEnd;
      Start = IFile->wn_sample( wnEnd - nu_F0 );
      End = IFile->wn_sample( wnStart - nu_F0 );
    }
  }
  if ( End == 0 && Start == 0 ) {
    // nl_error( 3, "All lines have been excluded!" );
    End = SignalEnd;
    Start = SignalStart;
  }
  assert( End > Start );
  npts = End - Start + 1;
  if ( npts > npts_vec ) {
    if ( npts_vec > 0 ) {
      free_vector(x,1,npts_vec);
      free_vector(sig,1,npts_vec);
    }
    npts_vec = npts;
    x = vector(1,npts_vec);
    sig =  vector(1,npts_vec);
    if ( x == 0 || sig == 0 )
      nl_error(3,"Out of memory resizing in fitdata::fit" );
  }

  for ( i = 1; i <= npts; i++ ) x[i] = i + Start - 1;
  y = yin + Start - 1;

  for ( i = 1; i <= ma; i++ ) a_save[i] = a[i];

  // And this should be done with exceptions:
  int counter = 0;
  int err_val = setjmp(Fit_buf);
  if ( err_val == 0 ) {
    if ( FitBaseline != 0 ) {
      // do fit for baseline avoiding any of the lines
      // This approach assumes the baseline function's parameters
      // are the first ones. Also assumes that the baseline function
      // by itself approximates the baseline, which is to say that
      // the function by which it is multiplied is very close to 1
      // away from the lines.
      //
      // As implemented, I only do this rudimentary fit on the
      // first line, so its purpose is simply to establish
      // initial values for the quad fit parameters. As such,
      // it is a kluge, and should be replaced by a proper
      // initialization from the config file. Let the Matlab routine
      // figure out what the value should be.
      //
      // In reality, I hope never to use the quadratic baseline
      // again, so I don't see a value in investing effort to clean
      // up its implementation, but in the short term, I probably
      // need the code to support comparisons. In the mid term, we
      // could get the same result as a quadratic fit by providing a
      // quadratic baseline file.
      
      // initialize sig to avoid the lines
      ICOS_Float nu_F0 = absorb->get_arg( a, 0 );
      for ( i = 1; i <= npts; i++ ) {
        sig[i] = 1;
        std::vector<argref>::iterator child;
        
        for (child = absorb->args.begin(); child != absorb->args.end(); ++child) {
          func_line *line = child->arg->is_line();
          if (line) {
            if ( x[i] > IFile->wn_sample(line->line_end(a) - nu_F0) &&
                 x[i] < IFile->wn_sample(line->line_start(a) - nu_F0) )
              sig[i] = 0;
          }
        }
      } 
      lfit(x,y,sig,npts,a,ia,base->params.size(),covar,&chisq, BaseFitFunction);
      FitBaseline = 0;
    }

    // Now initialize sig properly
    for ( i = 1; i <= npts; i++ ) sig[i] = GlobalData.Sigma;

    int converging = 0;
    int vctr = 0;
    ochisq = -1;
    chisq = 0;

    if ( verbose & 8 ) {
      if ( vfp != 0 ) fclose( vfp );
      vfp = 0;
      mlf_set_index( vmlf, IFile->mlf->index );
      if ( mlf_next_dir( vmlf ) ) {
        vfp = pathopen( vmlf->fpath, "ICOSsum.out", 0 );
      }
    }
    
    for ( counter = 0; counter < GlobalData.MaxIterations; counter++) {
      if ( verbose & 8 ) {
        FILE *vvfp = pathopen( vmlf->fpath, "%04d.dat", vctr );
        this->lwrite( vfp, vvfp, vctr++ );
      }
      if (verbose & 4) {
        func_evaluator::dump_evaluation_order.dump();
      }
      int mrqrv = mrqmin();
      if ( mrqrv ) {
        for ( i = 1; i <= ma; i++ ) a[i] = a_save[i];
        nl_error( 0,
          "Parameters rolled back, Retrying after disabling lines: %d",
          counter );
        // Now set the width of the newly disabled lines
        alamda = -1;
        adjust_params( a );
      } else {
        if ( verbose & 2 ) {
          fprintf( stderr, "%ld:%d:%d: chisq:%" FMT_G ", alamda:%" FMT_G,
              PTf->ScanNum, counter, converging, chisq, alamda );
          if ( ochisq != 0 ) {
            fprintf( stderr, " ochisq:%lg, dchisq:%lg",
               ochisq, (ochisq-chisq)/ochisq );
          }
          fprintf(stderr, "\n");
        }

        // These termination conditions are rather arbitrary. Plenty
        // of room here for tweaking.
        assert( ochisq >= 0 );
        if ( chisq <= ochisq ) {
          if ((ochisq-chisq)/ochisq <= GlobalData.ConvergenceStep ) {
            if ( ++converging >= GlobalData.ConvergenceCount ) {
              alamda=0.0;
              if ( adjust_params( a ) ) {
                nl_error( 0, "Retrying after enabling lines: %d", counter );
                for ( i = 1; i <= ma; i++ ) a_save[i] = a[i];
                alamda = -1; // Need to reinitialize mrqmin for new free params
                converging = 0;
              } else {
                mrqmin();
                return 1;
              }
            }
          } else {
            converging = 0;
          }
          ochisq = chisq;
        }
      }
    }
    nl_error( 1, "%s: Failed to converge", IFile->mlf->fpath );
    alamda = 0.0;
    mrqmin();
    return 1;
  } else {
    nl_error(1, "Failure after %d iterations", counter);
    func_evaluator::dump_evaluation_order.dump();
    return 0;
  }
}

int fitdata::n_input_params = 6;
const int fitdata::ScanNum_col = 1;

/**
  @param ofp Output FILE pointer for ICOSsum.dat
  @param vofp Output FILE pointer for verbose scan data if requested
  @param fileno The scan number
  Write the fit results. This function is used for the final fit
  as well as intermediate results (verbose & 8). It writes a single
  line to the ICOSsum.dat file and, if verbose & 1,
  also generates a separate verbose scan fit file.
 */
void fitdata::lwrite( FILE *ofp, FILE *vofp, int fileno ) {
  int i;
  if ( vofp != 0 ) {
    // Write verbose output file
    jmp_buf Fit_buf_save;
    memcpy( Fit_buf_save, Fit_buf, sizeof(Fit_buf) );
    if ( setjmp(Fit_buf) == 0 ) {
      func_evaluator::pre_evaluation_order.pre_eval(x[1], a);
      for ( i = 1; i <= npts; i++ ) {
        ICOS_Float yfit;
        func_evaluator::global_evaluation_order.evaluate(x[i], a);
        yfit = func->value;
        fprintf( vofp, "%12.6" FMT_E " %14.8" FMT_E " %12.6" FMT_E
          " %12.6" FMT_E " %12.6" FMT_E " %12.6" FMT_E,
          x[i], ICOSfile::wndata->data[i+Start-1],
          y[i], yfit, base->value, absorb->value );
        if (verbose & 128)
          absorb->print_intermediates(vofp);
        if (verbose & 16) {
          for (unsigned j = 0; j < func->params.size(); j++ ) {
            if (ia[func->params[j].index]) {
              fprintf( vofp, " %12.6" FMT_E, func->params[j].dyda );
            } else {
              fprintf( vofp, " 0" );
            }
          }
        }
        fprintf( vofp, "\n" );
      }
    } else {
      nl_error( 3, "unexpected longjmp during re-evaluation" );
    }
    memcpy( Fit_buf, Fit_buf_save, sizeof(Fit_buf_save) );
    fclose(vofp);
  }
  if ( ofp != 0 ) {
    // Write line of output to ICOSsum.dat
    int mfit = 0, n_i_p = 6;
    for ( i = 1; i <= ma; i++ ) {
      if ( ia[i] != 0 ) mfit++;
    }
    assert( ScanNum_col == 1 );
    fprintf( ofp, "%6d %6.2lf %6.2lf %12.5le %d %d",
      fileno, PTf->P, PTf->T, chisq/(End-Start-mfit+1),
      Start, End);
    { std::vector<argref>::iterator child;
      for (child = absorb->args.begin(); child != absorb->args.end(); ++child) {
        func_line *line = child->arg->is_line();
        if (line) {
          fprintf( ofp, " %d %12.5" FMT_E, line->fixed, line->S_thresh );
          n_i_p += 2;
        }
      }
    }
    assert(n_i_p == n_input_params);
    // This is not a simple recursive dump of parameters, as
    // that could be unnecessarily dependent on function build
    // order and structure. Instead, specific subclasses will
    // impose order as necessary.
    func->output_params(ofp, false);
    func->output_params(ofp, true);
    if (verbose & 128) {
      std::vector<argref>::iterator child;
      for (child = absorb->args.begin(); child != absorb->args.end(); ++child) {
        func_line *line = child->arg->is_line();
        if (line) {
          fprintf(ofp, " %13.7le %13.7le", line->Corr_Tref, line->Ks);
        }
      }
    }
    fprintf( ofp, "\n" );
    fflush( ofp );
  }
}

void fitdata::write() {
  FILE *fp = (verbose&1) ? IFile->writefp() : 0;
  this->lwrite( IFile->ofp, fp, PTf->ScanNum );
}
