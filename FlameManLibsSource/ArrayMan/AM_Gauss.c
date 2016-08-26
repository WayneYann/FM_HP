  /***************************************************************
 ******************************************************************
****		Gaussian Elimination with partial pivoting			****
****			 and condition number estimation.				****
****		This file contains the factorization driver and 	****
****		condition number estimation routine geco(), the		****
****		factorization routine gefa(), and the routine		****
****					routine gesl().							****
 ******************************************************************
  ****************************************************************/

/* © Copyright 1989,90,91,92 by Josef Goettgens. All rights reserved. */

#include <math.h>
#include "ArrayManager.h"

#ifdef __MWERKS__
#include <fp.h>
#endif

static int gesl_T( MatrixPtr a, int *ipvt, Double b[]);


int geco( MatrixPtr a, int ipvt[], Double *rcond, Double z[] )
/*
  PURPOSE
	  geco() factors a real matrix by gaussian elimination
	  and estimates the condition of the matrix.

  REMARKS
	  If  rcond  is not needed, gefa() is slightly faster.
	  to solve	A*x = b , follow geco() by gesl().
	  To compute  inverse(A)*c , follow geco() by gesl().
	  To compute  determinant(A) , follow geco() by gedi().
	  To compute  inverse(A) , follow geco() by gedi().

  INPUT
		a		A pointer to the matrix structure. 
				See the definition in ArrayManager.h.

  OUTPUT
		a		A pointer to the matrix structure containing
				an upper triangular matrix and the multipliers
				which were used to obtain it.
				The factorization can be written  a = l*u  where
				l  is a product of permutation and unit lower
				triangular matrices and  u	is upper triangular.
		ipvt	An integer vector (of length a->cols) of pivot indices.
		rcond	A Double estimate of the reciprocal condition of  A .
				for the system	A*x = b , relative perturbations
				in	A  and	b  of size	epsilon  may cause
				relative perturbations in  x  of size  epsilon/rcond .
				If	rcond  is so small that the logical expression
							1.0 + rcond .eq. 1.0
				is true, then  a  may be singular to working
				precision.	In particular,	rcond  is zero	if
				exact singularity is detected or the estimate
				underflows.
		 z		A Double vector (of length a->cols) for a work vector
				whose contents are usually unimportant.  If  A	is
				close to a singular matrix, then  z  is an approx-
				imate null vector in the sense that:
						 norm(a*z) = rcond*norm(a)*norm(z) .

	RETURNS
				= -1  Matrix is not square.
				=  0  Normal return value.
				=  k  if  u(k,k) .eq. 0.0 .  This is not an error
					  condition for this subroutine, but it does
					  indicate that gesl() or gedi() will divide by zero
					  if called.  Use  rcond  in geco() for a reliable
					  indication of singularity.

	ROUTINES
		gefa(), blas sasum() and sdot(), copysign(), fabs();

	WARNINGS
		This routine uses the UN*X math library routines
		copysign() and fabs().
*/
{
	/*register*/ int		j;
	int					k;
	int					n = a->cols, info = 0;
	/*register*/ Double	s, **m = a->mat;
#ifdef applec
	extended			ek, anorm, ynorm;
#else
	Double				ek, anorm, ynorm;
#endif

	/* Compute 1-norm of A.												*/
	for ( j = 0, anorm = 0.0; j < n; ++j ) {
		s = asum( n, a->mat[j], 1 );
		anorm = (s > anorm) ? s : anorm;
	}

	info = gefa( a, ipvt );				/* Factor A.					*/

	/*
	 * rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
	 * estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
	 * trans(a)  is the transpose of a .  The components of  e	are
	 * chosen to cause maximum local growth in the elements of w  where
	 * trans(u)*w = e .  The vectors are frequently rescaled to avoid
	 * overflow.
	 */

	ek = 1.0;							/* solve trans(u)*w = e			*/
	for ( j = 0; j < n; ++j ) z[j] = 0.0;

	for ( k = 0; k < n; ++k ) {
		
#ifdef applec
		register extended zk = z[k];
#else
		register Double	zk = z[k];
#endif
		Double	wk, wkm, sm;
		int		kp1 = k + 1;

		m[k][k] = m[k][k];
		
		if ( zk != 0.0 ) ek = copysign( ek, -zk );
		if ( fabs( ek - zk ) > fabs( m[k][k] ) ) {
			s   = fabs( m[k][k] )/ fabs( ek - zk );
			scal( n, s, z, 1 );
			zk  = z[k];
			ek *= s;
		}
		wk  = ek - zk;
		wkm = -ek - zk;
		s   = fabs( wk );
		sm  = fabs( wkm );
		if ( m[k][k] == 0.0 ) {
			wk  = 1.0;
			wkm = 1.0;
		} else {
			wk  /= m[k][k];
			wkm /= m[k][k];
		}
		if( kp1 < n ) {
			for ( j = kp1; j < n; ++j ) {
				sm    = sm + fabs( z[j] + wkm * m[j][k] );
				z[j] += wk * m[j][k];
				s    += fabs( z[j] );
			}
			if( s < sm ) {
#ifdef applec
				register extended t = wkm - wk;
#else
				register Double t = wkm - wk;
#endif
				wk = wkm;
				for( j = kp1; j < n; ++j ) z[j] += t * m[j][k];
			}
		}
		z[k] = wk;
	}
	s = 1.0 / asum( n, z, 1 );
	scal( n, s, z, 1 );

	for( k = n - 1; k >= 0; --k ) {			/* Solve trans(L)*y = w.		*/
		register int		l;
#ifdef applec
		register extended t;
#else
		register Double t;
#endif
		if ( k < n-1 ) z[k] += dot( n-k-1, a->mat[k]+ k+1, 1, (z+k+1), 1 );
		if ( fabs( z[k] ) > 1.0 ) {
			s = 1.0 / fabs( z[k] );
			scal( n, s, z, 1 );
		}
		l    = ipvt[k];
		t    = z[l];
		z[l] = z[k];
		z[k] = t;
	}
	s = 1.0 / asum( n, z, 1 );
	scal( n, s, z, 1);

	ynorm = 1.0;

	for ( k = 0; k < n; ++k ) {				/* Solve L*v = y.				*/
		register int		l;
#ifdef applec
		register extended t;
#else
		register Double t;
#endif
		l    = ipvt[k];
		t    = z[l];
		z[l] = z[k];
		z[k] = t;
		if ( k < n-1 ) saxpy( n-k-1, t, (a->mat[k]+k+1), 1, (z+k+1), 1 );
		if ( fabs( z[k] ) > 1.0)  {
			s = 1.0 / fabs( z[k] );
			scal( n, s, z, 1 );
			ynorm *= s;
		}
	}
	s = 1.0 / asum( n, z, 1 );
	scal( n, s, z, 1 );
	ynorm *= s;

	for( k = n-1; k >= 0; --k ) {			/* Solve  U*z = v.					*/
#ifdef applec
		register extended t;
#else
		register Double t;
#endif
		if ( fabs( z[k] ) > fabs( m[k][k] ) ) {
			s = fabs( m[k][k] ) / fabs( z[k] );
			scal( n, s, z, 1 );
			ynorm *= s;
		}
		if ( m[k][k] == 0.0 )	z[k] = 1.0;
		else					z[k] /= m[k][k];
		t = -z[k];
		saxpy( k, t, a->mat[k], 1, z, 1 );
	}

	s = 1.0 / asum( n, z, 1 );				/* Make znorm = 1.0.				*/
	scal( n, s, z, 1 );
	ynorm *= s;

	*rcond = (anorm == 0.0) ? 0.0 : ynorm / anorm;

	return info;
}



int gefa( MatrixPtr a, int *ipvt )
/*
	PURPOSE
		gefa() factors a real matrix by Gaussian elimination.
	
	REMARKS
		gefa() is usually called by geco(), but it can be called
		directly with a saving in time if  rcond  is not needed.
		(time for geco) = (1 + 9/n)*(time for gefa) .
	
	INPUT
		m		A pointer to the matrix structure. 
				See the definition in ArrayManager.h.
	
	OUTPUT
		m		A pointer to the matrix structure containing
				an upper triangular matrix and the multipliers
				which were used to obtain it.
				The factorization can be written  m = l*u  where
				l  is a product of permutation and unit lower
				triangular matrices and  u	is upper triangular.
		ipvt	An integer vector (of length m->cols) of pivot indices.
	
	RETURNS
				= -2  The partition scheme of a is not kColumnPointers
				= -1  Matrix is not square.
				=  0  Normal return value.
				=  k  if  u(k,k) .eq. 0.0 .  This is not an error
					  condition for this subroutine, but it does
					  indicate that sgesl or sgedi will divide by zero
					  if called.  Use  rcond  in sgeco for a reliable
					  indication of singularity.
	
	ROUTINES
		iamax
*/
{
	/*register*/ int	  i, j;
	int 			  k, l;
	int				  n = a->cols;			/* Number of equations				*/
	int				  nm1 = n - 1;			/* Number of equations - 1			*/
	int				  info = 0;				/* Assume nothing will go wrong!	*/
	Double			  *akk = a->mat[0];		/* Pointer to first column			*/
	Double			  *alk;
	/*register*/ Double t, *mik;
	
	/* Do some error checking.													*/
	if ( a->partition != kColumnPointers )	return -2;
	if ( a->rows != a->cols )				return -1;

	/* Gaussian elimination with partial pivoting.								*/
	if ( n < 2 ) goto CLEAN_UP;
	
	for ( k = 0; k < nm1; ++k, ++ipvt ) {	/*	Loop over Diagonal				*/
		
		/* Find index of max elem in col below the diagonal (l = pivot index).	*/
		akk   = a->mat[k] + k;
		l	  = iamax( n - k, akk, 1 ) + k;
		*ipvt = l;
		
		/* Zero pivot implies this column already triangularized.				*/
		alk = a->mat[k] + l;
		if ( *alk == 0.0 ) {
			info = k;
			continue;
		}
		
		/* Interchange a(k,k) and a(l,k) if necessary.							*/
		if ( l != k ) {
			t	 = *alk;
			*alk = *akk;
			*akk = t;
		}
		
		/* Compute multipliers for this column.									*/
		t = -1.0 / (*akk);
		for ( i = k+1, mik = a->mat[k]; i < n; ++i ) mik[i] *= t;
		
		/* Column elimination with row indexing.								*/
		if ( l != k ) {
			/* Interchange a(k,j) and a(l,j) if necessary.						*/
			/*register*/ Double **m = a->mat;
			for ( j = k+1; j < n; ++j ) {
				t = m[j][k];					/* t = pelem(a,k,j);			*/
				m[j][k] = m[j][l];				/* pelem(a,k,j) = pelem(a,l,j);	*/
				m[j][l] = t;					/* pelem(a,l,j) = t;			*/
			}
		}
		for ( j = k+1; j < n; ++j ) {
			/*register*/ Double *aij = a->mat[j];
			Double **m = a->mat;
			t = m[j][k];						/* t = pelem(a,k,j);			*/
			for ( i = k+1, mik = a->mat[k]; i < n; ++i )
				aij[i] += t * mik[i];
		}
	}	/* End of for k loop */
	
  CLEAN_UP:
	*ipvt = nm1;
	if ( *akk == 0.0 ) info = n;
	return info;
}


int gesl( MatrixPtr a, int *ipvt, Double b[], int job )
/*
	PURPOSE
		SGESL solves the real system
		a * x = b  or  trans(a) * x = b
		using the factors computed by geco() or gefa().

	INPUT
		a		A pointer to the matrix structure containing the factored
				matrix.  See the definition of Matrix in ArrayManager.h.
		ipvt	The pivot vector (of length a->cols) from geco() or gefa().
		b		The right hand side vector (of length a->cols).
		job 	= 0 		to solve  a*x = b ,
				= nonzero	to solve  trans(a)*x = b  where
							trans(a)  is the transpose.

	OUTPUT
		b		The solution vector x.

	REMARKS
		Error condition:
		A division by zero will occur if the input factor contains a
		zero on the diagonal.  Technically this indicates singularity
		but it is often caused by improper arguments or improper
		setting of lda .  It will not occur if the subroutines are
		called correctly and if sgeco has set rcond .gt. 0.0
		or sgefa has set info .eq. 0 .
*/
{
	/*register*/ Double t, *mik;
	Double			 *akk;
	Double			 **aMat = a->mat;
	/*register*/ int	 i, k;
	int				 l, n = a->cols, nm1 = n - 1;
	
	if ( job ) return ( gesl_T( a, ipvt, b ) );
	
	/* job = 0 , solve	A * x = b.												*/
	
	for ( k = 0; k < nm1; ++k ) {	/* Forward elimination. Solve L*y = b.		*/
	
		akk = aMat[k] + k; 		/* akk points to a(k,k).					*/
		
		l = ipvt[k];				/* Interchange b[k] and b[l] if necessary.	*/
		t = b[l];
		if ( l != k ) {
			b[l] = b[k];
			b[k] = t;
		}
		for ( i = k+1, mik = aMat[k]; i < n; ++i )
			b[i] += t * mik[i];
	}
	
	for ( k = nm1; k >= 0; --k ) {	/* Back substitution.  Solve  U*x = y.		*/
		/*register*/ Double *uik = aMat[k];
		akk = uik + k;
		b[k] /= (*akk);
		for ( i = 0; i < k; ++i ) b[i] -= uik[i] * b[k];
	}
	
	return 0;
}


static int gesl_T( MatrixPtr a, int *ipvt, Double b[] )
{
	/*register*/ Double t, *mik;
	Double			 *akk;
	Double			 **aMat = a->mat;
	/*register*/ int	 i, k;
	int				 l, n = a->cols, nm1 = n - 1;

	/* job = nonzero.  Solve  trans(A) * x = b.									*/
	
	/* First solve trans(U)*y = b.												*/
	for ( k = 0; k < n; ++k ) {
		/*register*/ Double *uik = aMat[k];
		akk = uik + k;
		t = 0.0;
		for ( i = 0; i < k; ++i ) t += uik[i]*b[i];
		b[k] = (b[k] - t) / (*akk);
	}
	
	/* b now contains y.  Solve trans(L)*x = y.									*/
	for ( k = n-2; k >= 0; --k ) {
		mik = aMat[k];
		t = 0.0;
		for ( i = k+1; i < n; ++i ) t += mik[i] * b[i];
		b[k] += t;
		
		l	 = ipvt[k];		/* Interchange b(k) and b(ipvt(k)) if necessary.	*/
		if( l == k ) continue;
		t	 = b[l];
		b[l] = b[k];
		b[k] = t;
	}
	return 0;
}
