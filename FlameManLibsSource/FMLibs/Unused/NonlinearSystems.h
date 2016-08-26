typedef void ( *NewtonFunc )( Double *x, Double **alpha, Double *bet );

void myludcmp( Double **a, int n, int *indx, Double *d );
void mylubksb( Double **a, int n, int *indx, Double *b );
void mnewt( int ntrial, Double *x, int n, Double tolx, Double tolf
		, NewtonFunc usrfun );
