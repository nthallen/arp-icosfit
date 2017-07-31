#ifndef F_VECTOR_H_INCLUDED
#define F_VECTOR_H_INCLUDED
#include "config.h"

class f_vector {
  public:
	ICOS_Float *data;
	int n_data;
	int datasize;
	int min_size;
	int offset;
	f_vector( int minsize, int offset=0 );
	void clear();
	void check( int size );
	void append( ICOS_Float f );
};

#endif
