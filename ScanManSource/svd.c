/*
	svd.c: "singular value decomposition" taken from numerical recipes
*/

#ifdef qUseDump
#pragma load "ScanMan.dump"
#else 
#include "ScanMan.h"
#endif


#ifdef applec
#include <SANE.h>
#endif

#include "svd.h"


#define SIGN(a,b)	copysign(b,a)

#ifdef HP		/* where is copysign? */

FPUType copysign( FPUType x, FPUType y )
{
	return (x >= 0) ? fabs(y) : -fabs(y);
}

#endif

/*static FPUType fmax( FPUType a, FPUType b )
{
	return (a > b) ? a : b;
}*/


static /*inline*/ FPUType pythag( FPUType a, FPUType b )
{
	register FPUType c;
	
	a = fabs( a );
	b = fabs( b );
	
	if ( a > b ) {
		c = b / a;
		c = a * sqrt(1.0 + c*c);
	} else if ( b != 0.0 ) {
		c = a / b;
		c = b * sqrt(1.0 + c*c);
	} else {
		c = 0.0;
	}
	
	return c;
}


/*void svdcmp( double **a, int m, int n, double w[], double **v, double *tmp )
{
	int		flag, i, its, j, jj, k, l, nm;
	double	*rv1 = NULL;
	FPUType	c, f, h, s, x, y, z;
	FPUType	anorm = 0.0, g = 0.0, scale = 0.0;

	if ( m < n ) FatalError("svdcmp(): You must augment A with extra zero rows");
	rv1 = (tmp) ? tmp : New1DArray( n );
	
	for ( i = 0; i < n; ++i ) {
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if ( i < m ) {
			for ( k = i; k < m; ++k ) scale += fabs( a[k][i] );
			if ( scale ) {
				for ( k = i; k < m; ++k ) {
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN( sqrt(s), f );
				h = f * g - s;
				a[i][i] = f - g;
				if ( i != n-1 ) {
					for ( j = l; j < n; ++j ) {
						for ( s = 0.0, k = i; k < m; ++k ) s += a[k][i] * a[k][j];
						f = s / h;
						for ( k = i; k < m; ++k ) a[k][j] += f * a[k][i];
					}
				}
				for ( k = i; k < m; ++k ) a[k][i] *= scale;
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if ( i < m && i != n-1 ) {
			for ( k = l; k < n; ++k ) scale += fabs(a[i][k]);
			if ( scale ) {
				for ( k = l; k < n; ++k ) {
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = -SIGN( sqrt(s), f );
				h = f * g - s;
				a[i][l] = f - g;
				for ( k = l; k < n; ++k ) rv1[k] = a[i][k] / h;
				if ( i != m-1 ) {
					for ( j = l; j < m; ++j ) {
						for ( s = 0.0, k = l; k < n; ++k ) s += a[j][k] * a[i][k];
						for ( k = l; k < n; ++k ) a[j][k] += s * rv1[k];
					}
				}
				for ( k = l; k < n; ++k ) a[i][k] *= scale;
			}
		}
		anorm = fmax( anorm, fabs(w[i])+fabs(rv1[i]) );
	}
	
	for ( i = n-1; i >= 0; --i ) {
		if ( i < n-1 ) {
			if ( g ) {
				for ( j = l; j < n; ++j )
					v[j][i] = (a[i][j] / a[i][l]) / g;
				for ( j = l; j < n; ++j ) {
					for ( s = 0.0, k = l; k < n; ++k ) s += a[i][k] * v[k][j];
					for ( k = l; k < n; ++k ) v[k][j] += s * v[k][i];
				}
			}
			for ( j = l; j < n; ++j) v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	
	for ( i = n-1; i >= 0 ; --i ) {
		l = i + 1;
		g = w[i];
		if ( i < n-1 )
			for ( j = l; j < n; ++j) a[i][j] = 0.0;
		if ( g ) {
			g = 1.0 / g;
			if ( i != n-1 ) {
				for ( j = l; j < n; ++j ) {
					for ( s = 0.0, k = l; k < m; ++k ) s += a[k][i] * a[k][j];
					f = (s / a[i][i]) * g;
					for ( k = i; k < m; ++k ) a[k][j] += f * a[k][i];
				}
			}
			for ( j = i; j < m; ++j ) a[j][i] *= g;
		} else {
			for ( j = i; j < m; ++j)  a[j][i] = 0.0;
		}
		a[i][i] += 1.0;
	}
	
	for ( k = n-1; k >= 0; --k ) {
		for ( its = 1;its <= 30; ++its ) {
			flag = 1;
			for ( l = k; l >= 0; --l ) {
				nm = l - 1;
				if ( fabs(rv1[l]) + anorm == anorm ) {
					flag = 0;
					break;
				}
				if ( fabs(w[nm]) + anorm == anorm ) break;
			}
			if ( flag ) {
				c = 0.0;
				s = 1.0;
				for ( i = l; i <= k; ++i ) {
					f = s * rv1[i];
					if ( fabs(f) + anorm != anorm ) {
						g = w[i];
						h = pythag( f, g );
						w[i] = h;
						h = 1.0 / h;
						c = g * h;
						s = (-f * h);
						for ( j = 0; j < m; ++j ) {
							y = a[j][nm];
							z = a[j][i];
							a[j][nm] = y * c + z * s;
							a[j][i] = z * c - y * s;
						}
					}
				}
			}
			z = w[k];
			if ( l == k ) {
				if ( z < 0.0 ) {
					w[k] = -z;
					for ( j = 0;j < n; ++j ) v[j][k] = (-v[j][k]);
				}
				break;
			}
			if ( its == 30 ) FatalError("No convergence in 30 SVDCMP iterations");
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag( f, 1.0 );
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;
			for ( j = l; j <= nm; ++j ) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;
				for ( jj = 0;jj < n; ++jj ) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z = pythag(f, h);
				w[j] = z;
				if (z) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = (c * g) + (s * y);
				x = (c * y) - (s * g);
				for ( jj = 0; jj < m; ++jj ) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	
	if ( !tmp ) Free1DArray( rv1 );
}
*/
#undef SIGN

/*
void svbksb( double **u, double w[], double **v, int m, 
			 int n, double b[], double x[], double *tmp )
{
	int		j, i;
	double	*ws;
	FPUType	s;
	
	ws = (tmp) ? tmp : New1DArray( n );

	for ( j = 0; j < n; ++j ) {
		s = 0.0;
		if ( w[j] ) {
			for ( i = 0; i < m; ++i ) s += u[i][j] * b[i];
			s /= w[j];
		}
		ws[j] = s;
	}
	for ( j = 0; j < n; ++j ) {
		s = 0.0;
		for ( i = 0; i < n; ++i ) s += v[j][i] * ws[i];
		x[j] = s;
	}
	
	if ( !tmp ) Free1DArray( ws );
}
*/

void svedit( double w[], int n, double tol )
{
	int		i;
	FPUType	threshold;
	
	threshold = w[0];
	for ( i = 1; i < n; ++i )
		if ( w[i] > threshold ) threshold = w[i];
	threshold *= tol;
	for ( i = 0; i < n; ++i )
		if ( w[i] < threshold ) w[i] = 0.0;
}


double svcond( double w[], int n )
{
	double	maxi, mini;
	
	MinMax( w, n, &maxi, &mini );
	return maxi / mini;
}


void svinverse( double **ainv, int n, double **u, double w[], double **v )
{
	int		i, j, k;
	FPUType	sum;
	
	for ( i = 0; i < n; ++i ) {
		for ( j = 0; j < n; ++j ) {
			sum = 0.0;
			for ( k = 0; k < n; ++k )
				sum += v[i][k] * u[j][k] / w[k];
			ainv[i][j] = sum;
		}
	}
}


int svrank( double w[], int n, double tol )
{
	int	i, rank = 0;
	
	for ( i = 0; i < n; ++i )
		if ( w[i] > tol ) ++rank;
	
	return rank;
}
