#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include "ICOSfit.h"
#include "global.h"

// func_base is a func_evaluator that uses a canonical
// baseline shape derived from zero air for the baseline.
// The baseline data is written in my ICOS-standard binary
// file format, which uses two 32-bit unsigned ints at
// the start to specify rows and columns, and then records
// the remaining data as floats. For the baseline, the first
// element of each column is the parameter initialization.

func_base_svdx::func_base_svdx( const char *filename ) :
    func_base( "func_base_svdx" ) {
  uses_nu_F0 = 0;
  FILE *fp = fopen( filename, "r" );
  if ( fp == 0 )
    nl_error( 3, "Unable to open baseline file %s", filename );
  icos_hdr_t header[2];
  if ( fread_swap32( &header, sizeof(icos_hdr_t), 2, fp ) != 2 )
    nl_error( 3, "%s: Error reading header: %s", filename,
      strerror(errno) );
  if ( header[0] <= 1 || header[1] <= 0 )
    nl_error( 3, "%s: Invalid header: %ld, %ld", filename,
              header[0], header[1] );
  n_params = header[1];
  params = new parameter[n_params];
  n_pts = header[0]-1;
  baseline = new float *[header[1]+1];
  if (baseline == 0) nl_error(3, "Out of memory in func_base" );
  unsigned int i;
  for ( i = 0; i < header[1]; i++ ) {
    /* The first element of each column is the initial parameter */
    baseline[i] = new float[header[0]];
    if ( baseline[i] == 0 )
      nl_error(3, "Out of memory in func_base");
    if ( fread_swap32( baseline[i], sizeof(float), header[0], fp )
          != header[0] )
      nl_error( 3, "%s: Error reading baseline: %s", filename,
        strerror(errno) );
    params[i].init = baseline[i][0];
  }
  fclose(fp);
}

void func_base_svdx::evaluate( float x, float *a ) {
  int i;
  int ix = (int)x;
  if ( ix < 1 || ix > n_pts )
    nl_error( 3,
      "x out of range in func_base::evaluate: %d", ix );
  value = 0.;
  for ( i = 0; i < n_params; i++ ) {
    float ai = get_param( a, i );
    float bix = baseline[i][ix];
    value += ai * bix;
    params[i].dyda = bix;
  }
}

func_base_ptbnu::func_base_ptbnu( const char *filename ) :
    func_base( "func_base_ptbnu" ) {
  uses_nu_F0 = 0;
  // The file format is specified in funceval.h
  FILE *fp = fopen( filename, "r" );
  icos_hdr_t header[2];
  if ( fread_swap32( &header, sizeof(icos_hdr_t), 2, fp ) != 2 )
    nl_error( 3, "%s: Error reading header: %s", filename,
      strerror(errno) );
  if (header[0]!=0 || header[1]!=1)
    nl_error(3,"Input file is not in ptbnu format: '%s'", filename );
  if ( sizeof(cfg) != 32 )
    nl_error(4, "ASSERT FAILURE: func_base_ptbnu sizeof(cfg) != 32 (%d)", sizeof(cfg) );
  if (fread(&cfg, sizeof(cfg), 1, fp) != 1)
    nl_error(3, "%s: Error reading cfg header: %s", filename,
      strerror(errno));
  #ifdef USE_BIG_ENDIAN
    cfg.poly_scale = bswap_64(cfg.poly_scale);
    cfg.nu0 = bswap_64(cfg.nu0);
    cfg.dnu = bswap_64(cfg.dnu);
    cfg.n_vectors = bswap_16(cfg.n_vectors);
    cfg.n_pts = bswap_16(cfg.n_pts);
    cfg.poly_coeffs = bswap_16(cfg.poly_coeffs);
    cfg.poly_of_nu = bswap_16(cfg.poly_of_nu);
  #endif
  if ( cfg.n_vectors > 0 && cfg.n_pts < 2 )
    nl_error(3, "%s: n_pts too small: %d", filename, cfg.n_pts );
  uses_nu_F0 = (cfg.n_vectors | cfg.poly_of_nu) ? 1 : 0;
  n_params = uses_nu_F0 + cfg.n_vectors + cfg.poly_coeffs;
  params = new parameter[n_params];
  // Read in initial parameter values
  int i;
  for ( i = 0; i < cfg.n_vectors; i++) {
    if ( fread_swap32( &params[i+uses_nu_F0].init, sizeof(float), 1, fp ) != 1 )
      nl_error( 3, "%s: Error reading vector param init: %s", filename,
        strerror(errno));
  }
  for ( i = 0; i < cfg.poly_coeffs; i++ ) {
    if ( fread_swap32( &params[uses_nu_F0+cfg.n_vectors+i].init, sizeof(float), 1, fp )
            != 1 )
      nl_error( 3, "%s: Error reading polynomial param init: %s", filename,
        strerror(errno));
  }
  if ( cfg.n_vectors ) {
    vectors = new float *[cfg.n_vectors];
    if (vectors == 0) nl_error(3, "Out of memory in func_base_ptbnu" );
    for ( i = 0; i < cfg.n_vectors; i++ ) {
      vectors[i] = new float[cfg.n_pts];
      if ( fread_swap32(vectors[i], sizeof(float), cfg.n_pts, fp) != cfg.n_pts )
        nl_error( 3, "%s: Error reading ptb vector: %s", filename,
          strerror(errno));
    }
  }
  fclose(fp);
  if ( cfg.n_vectors ) {
    dvdnu = new float *[cfg.n_vectors];
    if ( dvdnu == 0) nl_error(3, "Out of memory allocating dvdnu" );
    for ( i = 0; i < cfg.n_vectors; i++ ) {
      float *vi, *di;
      di = dvdnu[i] = new float[cfg.n_pts];
      vi = vectors[i];
      for (int j = 0; j < cfg.n_pts-1; j++ )
        di[j] = (vi[j+1]-vi[j])/cfg.dnu;
      di[cfg.n_pts-1] = di[cfg.n_pts-2];
    }
    if ( GlobalData.Verbosity & 4 ) {
      fprintf( stderr, "Baseline Nu0 = %.4lf  dnu = %.5lf\n", cfg.nu0, cfg.dnu );
      fprintf( stderr, "Baseline Vectors = [\n" );
      for ( int j = 0; j < cfg.n_pts; j++ ) {
        for ( i = 0; i < cfg.n_vectors; i++ ) {
          fprintf(stderr, "  %.5e", vectors[i][j] );
        }
        fprintf( stderr, "\n" );
      }
      fprintf( stderr, "];\n" );
      fprintf( stderr, "Baseline dVdnu = [\n" );
      for ( int j = 0; j < cfg.n_pts; j++ ) {
        for ( i = 0; i < cfg.n_vectors; i++ ) {
          fprintf(stderr, "  %.5e", dvdnu[i][j] );
        }
        fprintf( stderr, "\n" );
      }
      fprintf( stderr, "];\n" );
    }
  }
}

void func_base_ptbnu::init( float *a ) {
  int p2;
  for ( p2 = 0; p2 < n_params; p2++ )
    a[params[p2].index] = params[p2].init;

  cfg.nu0 -= func_line::nu0;
  // Now setup polynomial stuff
  if ( cfg.poly_coeffs > 0 && ! cfg.poly_of_nu ) {
    int nx = GlobalData.SignalRegion[1]+1;
    polyvecs = new float *[cfg.poly_coeffs-1]; // don't bother with constant
    if ( polyvecs == 0 ) nl_error(3, "Out of memory in func_base_ptbnu::init" );
    for ( int i = 0; i < cfg.poly_coeffs-1; i++ ) {
      polyvecs[i] = new float[nx];
      if ( polyvecs[i] == 0 ) nl_error( 3, "Out of memory in func_base_ptbnu::init" );
    }
    for ( int j = 0; j < nx; j++ ) {
      float x = j/cfg.poly_scale;
      float power = x;
      for ( int i = 0; i < cfg.poly_coeffs-1; i++ ) {
        polyvecs[i][j] = power;
        power *= x;
      }
    }
  }
}

// given x and parameters a, calculate value and
// params[].dyda
void func_base_ptbnu::evaluate( float x, float *a ) {
  int ix = int(x);
  float nu = 0.;

  value = 0;
  if ( uses_nu_F0 ) {
    float nu_F0 = get_param(a,0);
    nu = ICOSfile::wndata->data[ix] + nu_F0;
    params[0].dyda = 0;
  }
  if ( cfg.n_vectors ) {
    float bins = (nu-cfg.nu0)/cfg.dnu;
    if ( bins < 0. || bins >= cfg.n_pts )
      nl_error( 3, "Input nu (%.2f) out of range in func_base_ptbnu::evaluate", nu );
    int nui = (int) floor(bins);
    float fbin = bins - nui; // fraction of a bin
    for ( int i = 0; i < cfg.n_vectors; i++ ) {
      float ai = get_param(a, i+uses_nu_F0);
      float dvdnui = dvdnu[i][nui];
      float vnui = vectors[i][nui] + fbin * cfg.dnu * dvdnui;
      value += ai * vnui;
      params[0].dyda += ai * dvdnui;
      params[i+uses_nu_F0].dyda = vnui;
    }
  }
  // Now for the polynomials
  if ( cfg.poly_of_nu ) {
    float nupower = 1;
    float prevpower = 1;
    for ( int i = 0; i <= cfg.poly_coeffs; i++ ) {
      int pi = i + uses_nu_F0 + cfg.n_vectors;
      float ai = get_param(a,pi);
      value += ai * nupower;
      params[pi].dyda = nupower;
      params[0].dyda += i*ai*prevpower;
      prevpower = nupower;
      nupower *= nu;
    }
  } else {
    value += get_param(a,cfg.n_vectors+uses_nu_F0); // Constant
    params[cfg.n_vectors+uses_nu_F0].dyda = 1;
    for ( int i = 0; i < cfg.poly_coeffs-1; i++ ) {
      int pi = uses_nu_F0 + cfg.n_vectors + 1 + i;
      float xpower = polyvecs[i][ix];
      value += get_param(a,pi) * xpower;
      params[pi].dyda = xpower;
    }
  }
}

func_base_input::func_base_input( func_base *base ) :
            func_base("func_base_input") {
  append_func(base);
  n_params = base->n_params + 1;
  uses_nu_F0 = base->uses_nu_F0;
}

void func_base_input::init(float *a) {
  int p1, p2;
  a[params[uses_nu_F0].index] = 1;
  if ( first->params == 0 )
    first->params = new parameter[first->n_params];
  if ( uses_nu_F0 )
    link_param( 0, first, 0 );
  for ( p1 = 1+uses_nu_F0, p2 = uses_nu_F0; p2 < first->n_params; p2++ ) {
    link_param( p1, first, p2 );
    p1++;
  }
  first->init(a);
}

void func_base_input::evaluate( float x, float *a ) {
  int ix = int(x);
  func_evaluator::evaluate( x, a ); // evaluate base
  value = a[params[uses_nu_F0].index]*ICOSfile::bdata->data[ix] + first->value;
  params[uses_nu_F0].dyda = ICOSfile::bdata->data[ix];
  if ( uses_nu_F0 )
    params[0].dyda = first->params[0].dyda;
  int i;
  for ( i = 1+uses_nu_F0; i < n_params; i++ )
    params[i].dyda = first->params[i-1].dyda;
}

func_base *pick_base_type( const char *filename ) {
  func_base *base;
  FILE *fp = fopen( filename, "r" );
  if ( fp == 0 )
    nl_error( 3, "Unable to open baseline file %s", filename );
  icos_hdr_t header[2];
  if ( fread_swap32( &header, sizeof(icos_hdr_t), 2, fp ) != 2 )
    nl_error( 3, "%s: Error reading header: %s", filename,
      strerror(errno) );
  fclose(fp);
  if (header[0]) {
    base = new func_base_svdx(filename);
  } else if (header[1] == 1) {
    base = new func_base_ptbnu(filename);
  } else nl_error( 3, "Unrecognized baseline file format: %s", filename );
  if ( GlobalData.BaselineInput )
    base = new func_base_input( base );
  return base;
}
