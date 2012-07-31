#include <string.h>
#include "f_vector.h"
#include "nortlib.h"

// f_vector is a class of resizable vectors of floats.
// f_vector *fv = new f_vector( int min, int off );
//   min is the minimum number of floats for the first allocation.
//   as the vector grows, it will be in multiples of this size.
//   off is the index of the first element. It defaults to 0.
//   1 is the only common value.
// fv->append( 3.5 ); // Add 3.5 to the end of the vector
// fv->clear; // Set the number of points to 0
// fv->check( 7 ); // Make sure vector has room for 7 elements
// fv->data is the data vector
// fv->datasize is the size of the vector
// fv->n_data is the number of elements currently in use

f_vector::f_vector( int min, int off ) {
  data = 0;
  n_data = 0;
  datasize = 0;
  if ( min <= 0 || off < 0 )
    nl_error( 4, "Invalid minimum size or offset to f_vector" );
  offset = off;
  min_size = min+off;
}

void f_vector::clear() { n_data = 0; }

void f_vector::check( int size ) {
  float *newdata;
  if ( size+offset > datasize ) {
    if ( datasize == 0 )
      datasize = min_size;
    while ( datasize < size+offset && datasize >= 0 ) datasize *= 2;
    if ( datasize <= 0  )
      nl_error( 4, "datasize overflow in f_vector::check" );
    newdata = new float[datasize];
    if ( newdata == 0 ) nl_error( 4, "Out of memory in f_vector::check" );
    if ( data != 0 ) {
	  if ( n_data > 0 )
	    memcpy( newdata+offset, data+offset, n_data * sizeof(float) );
	  delete data;
    }
	data = newdata;
  }
}

void f_vector::append( float f ) {
  check(n_data+1);
  data[offset+n_data++] = f;
}

