/*
	ArrayManager.c: source code for ArrayManager functions
	
	MPW Library:
		Build ArrayLib.o.make

	HP720 Libaray:
		Build makefile

	© Copyright 1989-1993 by Josef Goettgens. All rights reserved.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>

#include "alligator.h"
#include "ArrayManager.h"

/*#ifndef HP*/
/*#define DBGNEW(t)		((t*)allocate_dbg(sizeof(t),file,line))*/
/*#define DBGNEW2(n,t)	((t*)allocate_dbg((unsigned)((n)*sizeof(t)),file,line))*/
/*#define DBGDELETE(ptr)	alligator_delete_dbg((ptr),file,line)*/
/*#else*/
#define DBGNEW(t)		((t*)allocate(sizeof(t),file,line))
#define DBGNEW2(n,t)	((t*)allocate((unsigned)((n)*sizeof(t)),file,line))
#define DBGDELETE(ptr)	alligator_delete((ptr))
/*#endif*/

/*
	void ClearArray( Double *ptr, int len );
	
	ClearArray() sets the len elements of the 1-D array pointed to
	by ptr to zero.
*/

void ClearArray( register Double *ptr, register int len )
{
	const FPUType zero = 0.0;
	register int i;
	
	for ( i = 0; i < len; ++i ) *ptr++ = zero;
}

void Clear1DArray( Double *v, int len )
{
	ClearArray(v , len );
}


void Clear2DArray( Double **m, int rows, int cols )
{
	ClearArray( *m , rows * cols );
}


void Clear3DArray( Double ***t, int planes, int rows, int cols )
{
	ClearArray( **t , planes * rows * cols );
}


void ClearVector( VectorPtr v )
{
	ClearArray( v->vec , v->phys_len );
}


void ClearMatrix( MatrixPtr m )
{
	ClearArray( *m->mat , m->phys_rows * m->phys_cols );
}


void ClearTensor( TensorPtr t )
{
	ClearArray( **t->tensor , t->phys_planes * t->phys_rows * t->phys_cols );
}



void ClearIntArray( register int *ptr, register int len )
{
	register int i;
	
	for ( i = 0; i < len; ++i ) *ptr++ = 0;
}


MatrixPtr SubMatrix( int plane, TensorPtr t )
{
	if ( plane >= 0 && plane < t->planes ) {
		t->mat = t->tensor[plane];
	}
	return (MatrixPtr)t;
}


static int CheckMIndex( int row, int col, MatrixPtr m )
{
	if ( !m )	return kNoData;
	
	/*	Check logical dimensions.										*/
	if ( row < 0 || row >= m->rows )
		return kInvalidLogicalRowIndex;
	if ( col < 0 || col >= m->cols )
		return kInvalidLogicalColIndex;

	/*	Check physical dimensions.										*/
	if ( row >= m->phys_rows )
		return kInvalidPhysicalRowIndex;
	if ( col >= m->phys_cols )
		return kInvalidPhysicalColIndex;
	
	return kIndexOK;
}


static int CheckTIndex( int plane, int row, int column, TensorPtr t )
{
	if ( !t )	return kNoData;

	if ( plane < 0 || plane >= t->planes )
		return kInvalidLogicalPlaneIndex;
	if ( plane >= t->phys_planes )
		return kInvalidPhysicalPlaneIndex;
	return CheckMIndex( row, column, SubMatrix( plane, t) );
}


int GetMElement( Double *value, int row, int column, MatrixPtr m )
{
	int		err;
	
	err = CheckMIndex(row, column, m);
	if ( err == kIndexOK ) {
		switch ( m->partition ) {
			case kRowPointers:
				*value = m->mat[row][column];
				break;
			case kColumnPointers:
				*value = m->mat[column][row];
				break;
			default:
				err = kInvalidPartition;
		}
	}
	return err;
}


int PutMElement( Double value, int row, int column, MatrixPtr m )
{
	int		err;
	
	err = CheckMIndex(row, column, m);
	if ( err == kIndexOK ) {
		switch ( m->partition ) {
			case kRowPointers:
				m->mat[row][column] = value;
				break;
			case kColumnPointers:
				m->mat[column][row] = value;
				break;
			default:
				err = kInvalidPartition;
		}
	}
	return err;
}


int GetTElement( Double *value, int plane, int row, int column, TensorPtr t )
{
	int		err;
	Double	***a;
	
	err = CheckTIndex( plane, row, column, t );
	if ( err == kIndexOK ) {
		a = t->tensor;
		switch ( t->partition ) {
			case kRowPointers:
				*value = a[plane][row][column];
				break;
			case kColumnPointers:
				*value = a[plane][column][row];
				break;
			default:
				err = kInvalidPartition;
		}
	}
	return err;
}


int PutTElement( Double value, int plane, int row, int column, TensorPtr t )
{
	int		err;
	Double	***a;
	
	err = CheckTIndex( plane, row, column, t );
	if ( err == kIndexOK ) {
		a = t->tensor;
		switch ( t->partition ) {
			case kRowPointers:
				a[plane][row][column] = value;
				break;
			case kColumnPointers:
				a[plane][column][row] = value;
				break;
			default:
				err = kInvalidPartition;
		}
	}
	return err;
}



/*
	void MinMax( const Double *x, int n, Double *maxi, Double *mini );
	
	MinMax() computes the maximum and minimum value of the array x
	of length n;
*/

void MinMax( const Double *x, int n, Double *maxi, Double *mini )
{
	int		i;
	FPUType	ma, mi;
	
	ma = mi = *x++;
	for ( i = 1; i < n; ++i, ++x ) {
		if ( *x > ma ) ma = *x;
		if ( *x < mi ) mi = *x;
	}
	*maxi = (Double)ma;
	*mini = (Double)mi;
}



/*
	Linear Algegra
*/

/*	void ludcmp( Double **a, int n, int *indx, Double *d, Double *vv );

	Taken from "Numerical Recipes in C" ( pt ).
	
	Given an n by n matrix a[0..n-1][0..n-1], this function replaces it 
	by the LU decomposition of a rowwise permutation of itself. 
	
	Double *a must be a pointer to rows (as opposed to columns!).
	
	a  and  n  are input.  a  is output, arranged as in equation (2.3.14) 
	above; indx[0..n-1]  is an output vector which records the row 
	permutation effected by the partial pivoting;  d  is output as +-1 
	depending on whether the number of row interchanges was even or 
	odd, respectively. This routine is used in combination with 
	lubksb to solve linear equations or invert a matrix.
	
	Changes to original NR version:
	o	type float changed to Double
	o	array indexing starts at 0
	o	additional parameter vv for external memory supply,
		which stores implicit row scaling
										( pt )
*/

void ludcmp( Double **a, int n, int *indx, Double *d, Double *vv )
{
	int		i, imax, j, k;
	Double	big, dum, sum, temp;
	const Double	kTiny = 1.0e-12;

	*d = 1.0;						/* No row interchanges yet. */
	for ( i = 0; i < n; ++i ) {		/* Loop over rows to get the implicit scaling information. */
		big = 0.0;
		for ( j = 0; j < n; ++j ) {
			if ( ( temp=fabs( a[i][j] ) ) > big) big = temp;
		}
		if (big == 0.0) {
			fprintf( stderr, "equation %d singular\n", i );
			FATAL( "Singular matrix in ludcmp()" );	/* No nonzero largest element. */
		}
		vv[i] = 1.0 / big;	/* Save the scaling. */
	}
	for ( j = 0; j < n ; ++j ) {	/* This is the loop over columns of Crout's method. */
		for ( i = 0; i < j ; ++i ) {	
			sum = a[i][j];
			for ( k = 0; k < i; ++k ) 
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;	/* Initialize for the search for largest pivot element. */
		for ( i = j; i < n; ++i ) {
			sum = a[i][j];
			for ( k = 0; k < j; ++k )
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ( ( dum = vv[i] * fabs( sum ) ) >= big ) {	/* Is the figure of merit for the pivot better than the best so far? */
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {					/* Do we need to interchange rows? */
			for ( k = 0; k < n; ++k ) {		/* Yes, do so */
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
		*d = -(*d);			/* change the parity of d */
		vv[imax] = vv[j];	/* also interchange the scale factor. */
		}
		indx[j] = imax;
		if ( a[j][j] == 0.0 ) a[j][j] = kTiny;
		/*	
			If the pivot element is zero the matrix is singular (at least to 
			the precision of the algorithm). For some applications on singular 
			matrices, it is desirable to substitute kTiny for zero.
		*/
		if (j != n - 1) {	/* Now, finally, divide by the pivot element. */
			dum = 1.0 / a[j][j];
			for ( i = j + 1; i < n; ++i ) a[i][j] *= dum;
		}
	}	/* Go back for the next column in the reduction. */
}


/*	void lubksb( Double **a, int n, const int *indx, Double *b );

	Taken from "Numerical Recipes in C" (pt).
	
	Solves the set of  n  linear equations  A X = B.  Here a[0..n-1][0..n-1]
	is input, not as the matrix A but rather as its LU decomposition, 
	determined by the routine ludcmp().  indx[0..n-1] is input as the permutation
	vector returned by ludcmp.  b[0..n-1]  is input as the right-hand side vector
	B, and returns with the solution vector X.  a,  n,  and  indx  are not modified
	by this routine and can be left in place for successive calls with different
	right-hand sides  b. This routine takes into account the possibility that  b
	will begin with many zero elements, so it is efficient for use in matrix
	inversion.
*/

void lubksb( Double **a, int n, const int *indx, Double *b )
{
	int 	i, ii = 0, ip, j;
	Double	sum;
/*	
	When ii is set to a not negative value, it will become the index of 
	the first nonvanishing element of b. We now do the forward
	substitution, equation 2.3.6. The only new wrinkle is to 
	unscramble the permutation as we go.
*/
	for ( i = 0; i < n; ++i ) {	
		ip = indx[i];
		sum = b[ip];	
		b[ip] = b[i];
		if ( ii > -1 ) {
			for ( j = ii; j <= i - 1; ++j ) sum -= a[i][j] * b[j];
		}
		else if ( sum ) ii = i;	{
			/* 	A nonzero element was encountered, so from now on we will 
				have to do the sums in the loop above. */
			b[i] = sum;
		}
	}	
	for ( i = n - 1; i >= 0; --i ) {	/* Now we do the backsubstitution */
		sum = b[i];
		for ( j = i + 1; j < n; ++j ) 
			sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i]; /* Store a component of the solution vector x. */
	}	/* All done! */
}


/*
	void thomas( Double *a, Double *b, Double *c, Double *r, int n );
	
	Given a tridiagonal system, described by the sub-, main-, and super-
	diagonal vectors, a, b, c, respectively, and the rhs-vector r, thomas()
	solves this linear system and returns the solution in r.
*/

void thomas( Double *a, Double *b, Double *c, Double *r, int n )
{
	register int i;
	register Double ar, br;
	register Double *rPtr = r, *cPtr = c + 1;
	
	*c /= *b;
	*rPtr++ /= *b;
	for ( i = 1; i < n; ++i ) {
		ar = *(a+i);
		br = *(b+i);
		br -= ar * *(cPtr-1);
		*cPtr++ /= br;
		*rPtr = ( *rPtr - ar * *(rPtr-1) ) / br;
		++rPtr;
	}
	i = n - 1;
	r += i;
	rPtr = r - 1;
	c += i - 1;
	do {
		*rPtr-- -= *c-- * *r--;
	} while ( --i );
}



#define	MIN(a,b)	((a)<(b)?(a):(b))

void DoolittleC( Double *r, Double **b, int nrows, int ncols, int ihalfb, int kkk )
/*
	DoolittleC solves a banded linear systems using Doolittle's method,
	where the matrix elements are stored according to b[col][row].

	Asymmetric band matrix solver doctored to ignore zeroes in lu decomp.
	
	b:		lhs matrix, dimensioned (neq,2*ihalfb+1) in calling function
	r:		rhs vector, dimensioned neq in calling function
	ihalfb:	half-bandwith, bandwith = 2*ihalfb+1

	Solution returns in r

	kkk = 1:	performs lu decomposition destructively
	kkk = 2:	performs back substitution
	kkk = 3:	performs options 1 & 2
	
	This C version was written by Alexandra Kees.
*/
{
	int		i, j, kc, jc;
	int		k, kd, lim, mr;
	Double	pivot, help, sum;
		
	if ( kkk == kBackSubstitution ) goto ModifyR;
	
	/*	Triangularize matrix b using Doolittle method
	 */
	for ( k = 0; k < nrows-1; ++k ) {
		pivot = b[ihalfb][k];
		kc = ihalfb;
		for ( i = k+1, --kc; (i < nrows) && (kc >= 0); ++i, --kc ) {
			help = -b[kc][i] / pivot;
			if ( help == 0.0 ) continue;
			b[kc][i] = help;
			lim = kc + ihalfb;
			for ( j = kc+1; j <= lim; ++j ) {
				jc = ihalfb + j - kc;
				b[j][i] += help * b[jc][k];
			}
		}
	}
#	ifdef DEBUG
	for(i = 0; i < nrows; i++)
		for(j = 0; j < ncols; j++)
			fprintf(stdout, "%d\t%d\t\t%g\n", i, j, b[j][i]);
#	endif
	
	if ( kkk == kLUDecomposition ) return;
	
 ModifyR:
	for (k = 1; k < nrows; ++k){
		kc = ihalfb - k;
		kd = 0.0;
		if (kc < 0){
			kc = 0.0;
			kd = k - ihalfb;
		}
		sum = 0.0;
		for (i = kc; i < ihalfb; ++i){
			sum += b[i][k] * r[kd];
			++kd;
		}
		r[k] += sum;
	}
#	ifdef DEBUG
	fprintf(stdout, "\n###Modified r\n");
	for(i = 0; i < nrows; ++i)
		fprintf(stdout, "r[%d] = %g\n", i, r[i]);	
#	endif
	
	/* back solution */
	 
	r[nrows-1] = r[nrows-1] / b[ihalfb][nrows-1];
	for(k = 2; k <= nrows; ++k){
		kc = nrows - k;
		jc = kc;
		mr = MIN (ncols, ihalfb + k);
		sum = 0.0;
		for(j = ihalfb+1; j < mr; j++){
			jc++;
			sum += b[j][kc] * r[jc];
		}
		r[kc] = (r[kc] - sum) / b[ihalfb][kc];
	}
#	ifdef DEBUG
	fprintf(stdout, "\n\n###SOLUTION:\n");
	for(i = 0; i < nrows; ++i)
		fprintf(stdout, "%g\n",  r[i]);	
#	endif

}


void DoolittleR( Double *r, Double **b, int nrows, int ncols, int ihalfb, int kkk )
/*
	DoolittleR solves a banded linear systems using Doolittle's method,
	where the matrix elements are stored according to b[row][col].

	Asymmetric band matrix solver doctored to ignore zeroes in lu decomp.
	
	b:		lhs matrix, dimensioned (neq,2*ihalfb+1) in calling function
	r:		rhs vector, dimensioned neq in calling function
	ihalfb:	half-bandwith, bandwith = 2*ihalfb+1

	Solution returns in r

	kkk = 1:	performs lu decomposition destructively
	kkk = 2:	performs back substitution
	kkk = 3:	performs options 1 & 2
	
	This C version was written by Alexandra Kees.
*/
{
	int		i, j, kc, jc;
	int		k, kd, lim, mr;
	Double	pivot, help, sum;
		
	if ( kkk == kBackSubstitution ) goto ModifyR;
	
	/*
	 *	triangularize matrix b using Doolittle method
	 */
	for ( k = 0; k < nrows-1; ++k ) {
		pivot = b[k][ihalfb];
		kc = ihalfb;
		for ( i = k+1, --kc; (i < nrows) && (kc >= 0); ++i, --kc ) {
			help = -b[i][kc] / pivot;
			if ( help == 0.0 ) continue;
			b[i][kc] = help;
			lim = kc + ihalfb;
			for ( j = kc+1; j <= lim; ++j ) {
				jc = ihalfb + j - kc;
				b[i][j] += help * b[k][jc];
			}
		}
	}
#	ifdef DEBUG
	for(i = 0; i < nrows; i++)
		for(j = 0; j < ncols; j++)
			fprintf(stdout, "\n%d\t%d\t\t%g", i, j, b[i][j]);
#	endif
	
	if ( kkk == kLUDecomposition ) return;
	
 ModifyR:
	for (k = 1; k < nrows; ++k){
		kc = ihalfb - k;
		kd = 0.0;
		if (kc < 0){
			kc = 0.0;
			kd = k - ihalfb;
		}
		sum = 0.0;
		for (i = kc; i < ihalfb; ++i){
			sum += b[k][i] * r[kd];
			++kd;
		}
		r[k] += sum;
	}
#	ifdef DEBUG
	fprintf(stdout, "\n###Modified r\n");
	for(i = 0; i < nrows; ++i)
		fprintf(stdout, "\nr[%d] = %g", i, r[i]);	
#	endif
		
	/* back solution */
	 
	r[nrows-1] = r[nrows-1] / b[nrows-1][ihalfb];
	for(k = 2; k <= nrows; ++k){
		kc = nrows - k;
		jc = kc;
		mr = MIN (ncols, ihalfb + k);
		sum = 0.0;
		for(j = ihalfb+1; j < mr; j++){
			jc++;
			sum += b[kc][j] * r[jc];
		}
		r[kc] = (r[kc] - sum) / b[kc][ihalfb];
	}
#	ifdef DEBUG
	fprintf(stdout, "\n\n###SOLUTION:\n");
	for(i = 0; i < nrows; ++i)
		fprintf(stdout, "\n%g",  r[i]);	
#	endif

}

#undef MIN
