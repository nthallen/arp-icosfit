#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include "ptread.h"
#include "nortlib.h"
#include "global.h"
#include "mlf.h"

PTfile::PTfile( const char *fname ) {
  fp = fopen( fname, "r" );
  ScanNum = next_ScanNum = 0;
  if ( fp == 0 )
    nl_error( nl_response, "Unable to open input file '%s'", fname );
  format = GlobalData.PTformat;
  last_file_pos = -1;
  switch ( format ) {
    case 0: n_vars = 12; break;
    case 1: n_vars = 7; break;
    case 2: n_vars = 11; break;
    default: nl_error( nl_response, "Unknown format code: %d", format );
  }
}

const int MYBUFSIZE = 256;
const int MAX_VARS = 12;

int PTfile::readline() {
  if ( fp == 0 ) return 0;
  if ( ScanNum < next_ScanNum ) {
    // This does not include the case where ScanNum==0, because then
    // next_ScanNum == 0 as well.
    // Note that if format != 2 (PTEfile) next_ScanNum will always be 0
    // As we are eliminating other PT formats, this can be ignored.
    ScanNum++;
    return 1;
  }
  for (;;) {
    char buf[MYBUFSIZE], *p, *ep;
    double data[MAX_VARS];
    int i;

    last_file_pos = ftell(fp);
    if ( fgets( buf, MYBUFSIZE, fp ) == 0 ) {
      fclose(fp);
      fp = 0;
      return 0;
    }
    for ( p = buf, i = 0; i <= n_vars; i++ ) {
      data[i] = strtod( p, &ep );
      if ( i == n_vars ) {
        GlobalData.input.nu_F0 = (p == ep) ? 0. : data[i];
      } else if ( p == ep ) {
        nl_error( 2, "Invalid number of parameters in PTFile\n" );
        fclose(fp);
        fp = 0;
        return 0;
      }
      p = ep;
    }
    if ( format == 2 ) {
      time = 0.;
      ScanNum = int(data[0]);
      P = data[1];
      T = data[2];
      for ( i = 0; i < 8; i++ ) Etln_params[i] = data[i+3];
      return 1;
    } else {
      time = data[0];
      P = data[1];
      if ( format == 0 ) {
        T = 273.15 + ( data[2] + data[3] + data[4] + data[5] ) / 4.;
        next_ScanNum = int(data[6]);
        cal_flow = data[8];
        inlet_flow = data[9];
        RORIS = int(data[10]);
        RateS = int(data[11]);
      } else {
        T = data[2];
        next_ScanNum = int(data[3]);
        cal_flow = data[4];
        inlet_flow = data[5];
        RORIS = int(data[6]);
        RateS = 0;
      }
      if ( T < 249. ) T = 273.15 + GlobalData.DefaultTemp;
      if ( ScanNum != next_ScanNum ) {
        if ( RORIS == GlobalData.QCLI_Wave ) {
          ScanNum = (ScanNum==0) ? next_ScanNum :
            ( (ScanNum<next_ScanNum) ? ScanNum+1 : ScanNum-1 );
          return 1;
        } else ScanNum = next_ScanNum;
      }
    }
  }
}

/**
 * Repositions read pointer to before the most recently read line.
 * Used during restart processing.
 */
void PTfile::backup() {
  if (last_file_pos >= 0)
    fseek(fp, last_file_pos, SEEK_SET);
}

void PTfile::calc_wndata() {
  int from = GlobalData.SignalRegion[0];
  int to = GlobalData.SignalRegion[1];
  int i;
  
  if ( GlobalData.PTformat != 2 )
    nl_error( 4, "calc_wndata called erroneously" );
  if ( ICOSfile::wndata == 0 ) ICOSfile::wndata = new f_vector( to, 1 );
  ICOSfile::wndata->check( to );
  for ( i = from; i <= to; i++ ) {
    double ii = (i - Etln_params[0] + 1) * 1e-3;
    double fn = Etln_params[1] + Etln_params[2]*ii + Etln_params[3]*ii*ii
      + Etln_params[4]*exp(-ii/Etln_params[5])
      + Etln_params[6]*exp(-ii/Etln_params[7]);
    ICOSfile::wndata->data[i] = -GlobalData.EtalonFSR * fn;
  }
  ICOSfile::wndata->n_data = to;
}

const int ICOSfile::mindatasize = 1024;

ICOSfile::ICOSfile( const char *fbase, const char *obase, int bin ) {
  binary = bin;
  mlf = mlf_init( 3, 60, 0, fbase, "dat", NULL );
  omlf = mlf_init( 3, 60, 1, obase, "dat", NULL );
  ofp = 0;
  sdata = new f_vector(mindatasize, 1);
  edata = new f_vector(mindatasize, 1);
  bdata = new f_vector(mindatasize, 1);
  fdata = new f_vector(100, 0);
}

f_vector *ICOSfile::bdata;
f_vector *ICOSfile::wndata;
int ICOSfile::dFN;
f_vector *wndebug;

#ifdef NOT_IMPLEMENTED
  static void err_throw( int except, int level, char *fmt, ... ) {
    va_list arg;

    va_start(arg, fmt);
    nl_verror(stderr, level, fmt, arg);
    va_end(arg);
    throw except;
  }
#endif

// returns 1 on success, 0 if there was an error
// For ASCII Processing:
//  fgets() reads in the line including the newline
//  strtod() updates the endptr (ep) to point to the char
//  immediately following the converted value. This should
//  satisfy isspace(*ep) for every successful conversion,
//  and since strtod swallows whitespace at the beginning,
//  it should be false for every unsuccessful conversion.
//
//  An ASCII file is deemed to have an etalon if there are
//  two numbers on a line. Clearly it would make sense
//  to make this determination once per file and avoid
//  testing each line.

int ICOSfile::read( unsigned long int fileno ) {
  FILE *fp;
  int has_etalon = 1;
  mlf_set_index( mlf, fileno );
  fp = mlf_next_file(mlf);
  if ( fp == 0 ) return 0;
  sdata->clear();
  edata->clear();
  fdata->clear();
  bdata->clear();
  if ( binary ) {
    icos_hdr_t header[2];
    if ( fread_swap32( header, sizeof(icos_hdr_t), 2, fp ) != 2 ) {
      nl_error( 2, "%s: Error reading header: %s", mlf->fpath,
        strerror(errno) );
      fclose(fp);
      return 0;
    }
    // Support for new SSP file format
    if (header[0] == 0x10006 && header[1] > 255) {
      unsigned long data[4];
      if (fread(data, sizeof(unsigned long), 4, fp) != 4) {
        nl_error( 2, "%s: Error reading SSP header: %s", mlf->fpath,
          strerror(errno) );
        fclose(fp);
        return 0;
      }
      header[0] = header[1]>>16;
      header[1] &= 0xFF;
    }
    if ( header[0] <= 0 || header[1] <= 0 || header[1] > 3 ) {
      nl_error( 2, "%s: Invalid header ( %ld, %ld )", mlf->fpath,
        header[0], header[1] );
      fclose(fp);
      return 0;
    }
    has_etalon = ( header[1] >= 2 );
    sdata->check(header[0]);
    sdata->n_data =
      fread_swap32( sdata->data+sdata->offset, sizeof(float), header[0], fp );
    if ( sdata->n_data != (int)header[0] ) {
      nl_error( 2, "%s: Error reading sdata: %s", strerror(errno) );
      fclose(fp);
      return 0;
    }
    if ( has_etalon ) {
      edata->check(header[0]);
      edata->n_data =
        fread_swap32( edata->data+edata->offset, sizeof(float), header[0], fp );
      if ( edata->n_data != (int)header[0] ) {
        nl_error( 2, "%s: Error reading edata: %s", strerror(errno) );
        fclose(fp);
        return 0;
      }
    }
    if ( header[1] >= 3 ) {
      bdata->check(header[0]);
      bdata->n_data =
        fread_swap32( bdata->data+bdata->offset, sizeof(float), header[0], fp );
      if ( bdata->n_data != (int)header[0] ) {
        nl_error( 2, "%s: Error reading bdata: %s", strerror(errno) );
        fclose(fp);
        return 0;
      }
    }
  } else {
    for (;;) {
      char buf[MYBUFSIZE], *ep;
      double value;

      if ( fgets( buf, MYBUFSIZE, fp ) == 0 ) {
        fclose(fp);
        break;
      }
      value = strtod( buf, &ep );
      if ( !isspace(*ep) ) {
        nl_error( 2, "%s:%d: No value read", mlf->fpath, sdata->n_data+1 );
        fclose(fp);
        return 0;
      }
      sdata->append(value);
      if ( has_etalon ) {
        value = strtod( ep, &ep );
        if ( !isspace(*ep) ) {
          if ( edata->n_data == 0 ) has_etalon = 0;
          else {
            nl_error( 2, "%s:%d: Expected second value",
              mlf->fpath, sdata->n_data );
            fclose(fp);
            return 0;
          }
        } else edata->append(value);
      }
    }
  }
  fclose(fp);
  if ( GlobalData.PTformat != 2 && edata->n_data > 0
          && fit_fringes(fileno) == 0 ) return 0;

  // This first baseline calculation is the zero baseline which
  // precedes the laser-on ICOS data. It should not be confused
  // with the baseline function which follows which attempts to
  // match the laser power function.
  if ( GlobalData.BackgroundRegion[0] <=
       GlobalData.BackgroundRegion[1] ) {
    float baseline = 0., *yin = sdata->data;
    unsigned int i;
    for ( i = GlobalData.BackgroundRegion[0];
          i <= GlobalData.BackgroundRegion[1];
          i++ ) baseline += yin[i];
    baseline /=
      GlobalData.BackgroundRegion[1]-GlobalData.BackgroundRegion[0]+1;
    for ( i = 1; (int)i <= sdata->n_data; i++ )
      yin[i] = yin[i] - baseline;
  }
  return 1;
}

// return 1 if there are no problems, otherwise
// return 0.
int ICOSfile::fit_fringes( unsigned long int fileno ) {
  if ( wndata == 0 ) wndata = new f_vector( edata->n_data, 1 );
  wndata->check( edata->n_data );
  wndebug = wndata;

  int i, j, rv = 1;
  int from = GlobalData.SignalRegion[0];
  int to = GlobalData.SignalRegion[1];

  if ( GlobalData.TuningRate != 0 ) {
    if (wndata->n_data == 0) {
      float wn = 0.;
      for ( i = from; i <= to; i++ ) {
        wndata->data[i] = wn;
        wn -= GlobalData.TuningRate;
      }
      wndata->n_data = to+1;
    }
  } else {
    // These initializations should be done
    // during the ICOSfile contruction
    // const int n = 3;
    // const float sx2 = 28; // sum(x^2) from -n to n
    // const float sx4 = 196; // sum(x^4)
    // const int n = 15;
    // const float sx2 = 2480; // sum(x^2) from -n to n
    // const float sx4 = 356624; // sum(x^4)
    // const int N = 2*n+1;
    // const float sx22 = sx2*sx2;
    // const float b2den = N*sx4-sx22; // 588
    static int n = 0, N;
    static float sx2, sx4, sx22, b2den;
    float prev_b2a = 0.;
    float prev_a1 = 0.;
    int new_n = int(floor(GlobalData.MinimumFringeSpacing/4));
    int decreasing = -1; // -1 means uninitialized, 0 = b2 is increasing, 1 = decreasing
    int filterwid[400]; // For debugging
    
    // Evaluate b2 for each set of points
    // b2(n-1) b2(n) b2(n+1)
    // if b2(n-1) >= b2(n) && b2(n) < b2(n+1)
    // && b2(n) < 0, then we've located a
    // fringe peak at n, and we want to set
    // the location of that peak to n-a1
    // or n+b1/(2*b2). This means we need
    // to know which way we're heading and
    
    if ( edata->n_data < from ) return 0;
    if ( edata->n_data < to ) to = edata->n_data;
    for ( i = from; i < to; i++ ) {
      if ( new_n != n ) {
        int ni;
        n = new_n;
        N = n*2 + 1;
        sx2 = sx4 = 0;
        for ( ni = -n; ni <= n; ni++ ) {
          float x2 = ni*ni;
          sx2 += x2;
          sx4 += x2*x2;
        }
        sx22 = sx2*sx2;
        b2den = N*sx4 - sx22;
      }
      if ( i+N >= to ) break;
          
      float sy = 0., syx = 0., syx2 = 0.;
      int j;
      for ( j = -n; j <= n; j++ ) {
        float y = edata->data[i+j+n];
        sy += y;
        syx += y*j;
        syx2 += y*j*j;
      }
      float b2a = (N*syx2 - sy*sx2);
      if ( decreasing < 0 ) decreasing = 0;
      else if ( b2a < prev_b2a ) {
        decreasing = 1;
        float b2 = b2a/b2den;
        float b1 = syx/sx2;
        prev_a1 = -b1/(2*b2);
      } else if ( decreasing > 0 && b2a > prev_b2a && prev_b2a < 0) {
        if ( prev_a1 > -n && prev_a1 < n ) {
          float new_fr = i + n - 1 + prev_a1;
          filterwid[fdata->n_data] = n;
          fdata->append( new_fr );
          int nd = fdata->n_data;
          if ( nd >= 3 ) {
            float dsfr = (new_fr - fdata->data[nd-2]) /
              ( fdata->data[nd-2] - fdata->data[nd-3] );
            if ( dsfr < GlobalData.DSFRLimits[0] ||
                 dsfr > GlobalData.DSFRLimits[1] ) {
              nl_error( 1,
                "%lu: Etalon dsfr out of range at "
                "fringe %d, sample %d: %f",
                 mlf->index, nd, i, dsfr );
              // return 0;
              rv = 0; break;
            }
          }
          if ( nd >= 2 ) {
            new_n = int(floor((fdata->data[nd-1] - fdata->data[nd-2])/4));
          }
          //for ( i += new_n; i < to; i++ ) {
          //	if ( edata->data[i] > edata->data[i-1] ) break;
          //}
          decreasing = -1;
        } else {
          decreasing = 0;
        }
      }
      prev_b2a = b2a;
    }
    if ( fdata->n_data < 2 )
      nl_error( 3, "Unable to locate more than one fringe" );
    if ( GlobalData.Verbosity & 64 ) {
      fprintf(stderr, " FR: %ld", fileno );
      for ( i = 0; i < fdata->n_data; i++ )
        fprintf(stderr, " %.2lf", fdata->data[i] );
      fprintf(stderr, "\n FF: %ld", fileno );
      for ( i = 0; i < fdata->n_data; i++ )
        fprintf(stderr, " %d", filterwid[i] );
      fprintf(stderr, "\n" );
    }
    if (rv == 0) return 0;
    if ( wndata->n_data > fdata->data[0] ) {
      dFN = -int(floor(wndata->data[ int(floor(fdata->data[0]+.5)) ]/
              GlobalData.EtalonFSR + .5));
    }
    j = 0;
    for ( i = GlobalData.SignalRegion[0]; i <= (int)GlobalData.SignalRegion[1]; j++ ) {
      int limit;
      double slope = -GlobalData.EtalonFSR/(fdata->data[j+1] - fdata->data[j]);
      double wn = -(j+dFN)*GlobalData.EtalonFSR + (i - fdata->data[j])*slope;
      if (j + 2 == fdata->n_data || GlobalData.SignalRegion[1] < fdata->data[j+1] )
        limit = GlobalData.SignalRegion[1];
      else limit = int(floor(fdata->data[j+1]));
      for (; i <= limit; i++ ) {
        wndata->data[i] = wn;
        wn += slope;
      }
    }
    wndata->n_data = GlobalData.SignalRegion[1];
  }
  return 1;
}

// Translate wavenumber back to sample number
// We will now assume that nu_F0 (the free parameter) has been
// subtracted from wn before the call, and hence wn is suitable
// for direct lookup in wndata.
int ICOSfile::wn_sample( float wn ) {
  // assert( nu_F0 > 0 );
  assert( wndata->n_data >= (int)GlobalData.SignalRegion[1] );
  // Note that wavenumber decreases with sample, so I use 'low'
  // to identify the lower wavenumber value, although it is the
  // higher index number
  // wn -= nu_F0;
  int low = GlobalData.SignalRegion[1];
  int high = GlobalData.SignalRegion[0];
  float wnlow = wndata->data[low];
  float wnhigh = wndata->data[high];
  if ( wn <= wnlow ) return low;
  if ( wn >= wnhigh ) return high;
  while ( low > high ) {
    int mid = int(low + (wn-wnlow)*(high-low)/(wnhigh-wnlow));
    assert( high <= mid && mid <= low );
    if ( mid == low ) return low;
    if ( mid == high ) return high;
    float wnmid = wndata->data[mid];
    if ( wn < wnmid ) {
      high = mid;
      wnhigh = wnmid;
    } else {
      low = mid;
      wnlow = wnmid;
    }
  }
  nl_error( 4, "wn_sample failed" );
  return 0;
}

FILE *ICOSfile::writefp() {
  FILE *fp;
  mlf_set_index( omlf, mlf->index );
  fp = mlf_next_file(omlf);
  return fp;
}

#ifdef USE_BIG_ENDIAN
 int fread_swap32( void *buf, size_t size, size_t count, FILE *fp ) {
   int rv = fread( buf, size, count, fp );
   icos_hdr_t *bptr = (icos_hdr_t *)buf;
   if (size != 4) nl_error( 4, "fread_swap32 requires size == 4" );
   int i;
   for ( i = 0; i < rv; i++ ) {
     bptr[i] = bswap_32(bptr[i]);
   }
   return rv;
 }
#endif
