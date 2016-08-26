/*
	svd.h: header file for "singular value decomposition" package
*/

#ifndef __svd__
#define __svd__

void svdcmp( double **a, int m, int n, double w[], double **v, double *tmp );
void svbksb( double **u, double w[], double **v, int m, 
	int n, double b[], double x[], double *tmp );
void svedit( double w[], int n, double tol );
double svcond( double w[], int n );
void svinverse( double **ainv, int n, double **u, double w[], double **v );
int svrank( double w[], int n, double tol );

#endif	/* __svd__ */
