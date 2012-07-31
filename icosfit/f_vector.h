class f_vector {
  public:
	float *data;
	int n_data;
	int datasize;
	int min_size;
	int offset;
	f_vector( int minsize, int offset=0 );
	void clear();
	void check( int size );
	void append( float f );
};

