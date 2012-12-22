#include "ICOSfit.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include "clp.h"
#include "mlf.h"
#include "global.h"

int verbose = 8;
/* verbose & 1 => output fits
   verbose & 2 => output gaussj [apparently no longer in use but see 32]
                  output line_check changes
   verbose & 4 => output baseline fit [apparently no longer in use]
                  output baseline vectors after reading in (func_base)
   verbose & 8 => output partial fits
   verbose & 16 => output derivatives in output files
   verbose & 32 => output ochisq, alamda and covariance matrix on each iteration
   verbose & 64 => output fringe positions in fit_fringes()
*/
void ICOS_init() {
  if (ShowVersion) {
    printf("icosfit version %s, %s\n", ICOSFIT_VERSION, ICOSFIT_VERSION_DATE);
    exit(0);
  }
  #if HAVE_LIBMALLOC_G
	mallopt(MALLOC_CKACCESS, 1);
	mallopt(MALLOC_FILLAREA, 1);
  #endif
}

static const char *output_filename( const char *name ) {
  static char fname[PATH_MAX];

  if ( name[0] == '/' ) return name;
  snprintf( fname, PATH_MAX-1, "%s/%s", GlobalData.OutputDir, name );
  fname[PATH_MAX-1] = '\0';
  return fname;
}

void ICOS_main() {
  fitdata *fitspecs = build_func();
  if ( GlobalData.LogFile != 0 ) {
    const char *fname;
    char pipename[PATH_MAX+14];
    FILE *fp;

    fname = output_filename( GlobalData.LogFile );
    snprintf( pipename, PATH_MAX+13, "/usr/bin/tee -a %s", fname );
    pipename[PATH_MAX+13] = '\0';
    fp = popen( pipename, "w" );
    if ( fp == 0 ) nl_error( 3, "Unable to create pipe to %s", pipename );
    if ( dup2( fileno(fp), 1 ) == -1 )
       nl_error( 3, "Unable to dup2: %s", strerror(errno) );
    if ( dup2( fileno(fp), 2 ) == -1 )
      nl_error( 3, "Unable to dup stdout to stderr: %s", strerror(errno) );
    // fclose(fp);
    if ( GlobalData.RestartAt <= 0 )
	  fprintf( stderr, "ICOSfit Start\n" );
	else fprintf( stderr, "\nICOSfit Restart at %d\n",
	   GlobalData.RestartAt );
  }
  while ( fitspecs->PTf->readline() != 0 ) {
    if ( ( GlobalData.ScanNumRange[0] == 0 ||
		   fitspecs->PTf->ScanNum >= GlobalData.ScanNumRange[0] ) &&
		 ( GlobalData.ScanNumRange[1] == 0 ||
		   fitspecs->PTf->ScanNum <= GlobalData.ScanNumRange[1] ) ) {
	  if ( fitspecs->IFile->read( fitspecs->PTf->ScanNum ) ) {
		if ( GlobalData.PTformat == 2 ) fitspecs->PTf->calc_wndata();
		if ( fitspecs->fit() != 0 ) {
		  fitspecs->write();
		  fprintf(stderr, "Successfully fit %lu: chisq = %f\n",
				   fitspecs->IFile->mlf->index,
				   fitspecs->chisq );
		} else {
		  fprintf( stderr, "Failed to fit %lu\n", fitspecs->IFile->mlf->index );
		  exit(1);
		}
	  }
	}
  }
  #if HAVE_LIBMALLOC_G
    fprintf( stderr, "Checking Heap\n" );
    malloc_dump_unreferenced( 2, 1 );
  #endif
}

fitdata *build_func() {
  PTfile *ptf = new PTfile( GlobalData.PTFile );
  ICOSfile *IF = new ICOSfile( GlobalData.ICOSdir,
     GlobalData.OutputDir, GlobalData.binary );

  func_evaluator *func;
  func_base *base;
  if ( GlobalData.BaselineFile == 0 )
    nl_error(3, "BaselineFile is now required" );
  base = pick_base_type( GlobalData.BaselineFile );
  func_abs *abs = GlobalData.absorb;
  if ( GlobalData.N_Passes > 0 ) {
    func = new func_noskew( base, abs );
    nl_error( 0, "Using func_noskew()" );
  } else {
    if ( GlobalData.SampleRate == 0 )
      nl_error(3, "SampleRate required for skew calculation" );
    func = new func_skew( base, abs );
    nl_error( 0, "Using func_skew()" );
  }
  fitdata *fd = new fitdata( ptf, IF, func, base, abs );
  
  // This disables the func_quad code in fitdata::fit
  // This is kluge and will be eliminated soon.
  fd->FitBaseline = ( GlobalData.BaselineFile == 0 );

  fd->handle_restart( output_filename( GlobalData.OutputFile ) );

  { const char *fnamep;
    FILE *fp;

    assert( GlobalData.MFile != 0 );
    fnamep = output_filename( GlobalData.MFile );
    fp = fopen( fnamep, "w" );
    if ( fp == 0 ) nl_error( 3, "Unable to open MFile '%s'", fnamep );
    assert( abs != 0 && abs->first != 0 && abs->first->params != 0 );
    fprintf( fp,
      "%% ICOS configuration data\n"
      "n_input_params = %d;\n"
      "n_base_params = %d;\n"
      "binary = %d;\n"
      "nu0 = %.0f;\n",
      fd->n_input_params,
      abs->params[0].index - 1,
      GlobalData.binary,
      func_line::nu0 );
    fprintf( fp, "BaselineFile = '%s';\n",
      GlobalData.BaselineFile ? GlobalData.BaselineFile : "" );
    fprintf( fp, "PTEfile = '%s';\n",
      (GlobalData.PTformat == 2 && GlobalData.PTFile) ?
	GlobalData.PTFile : "" );
    abs->print_config( fp );
    fclose(fp);
  }
  return fd;
}
