/*
	AM_BlockTriDiagSolver.c: Routines to solve linear block tri-diagonal systems
	
	History:
		02/25/93	calls to SpinCursor()
		09/19/92	calls to rc_dot now contain the correct physical offset (jg)
		04/24/92	first release (jg)

	© Copyright 1989,90,91,92 by Josef Goettgens. All rights reserved.
*/

#include "ArrayManager.h"

#ifdef applec
#include <CursorCtl.h>
extern int StandAlone;
#endif


static FPUType rc_dot(int n, const Double *r, int inc_r, 
						const Double *c, int inc_c)
{
	const int	m = n % 6;
	int			i;
	FPUType		sum = 0.0;
	
	if ( m ) {
		for ( i = 0; i < m; ++i ) {
			sum += *r * *c++;	r += inc_r;
		}
		if ( n < 6 ) return sum;
	}
	n /= 6;
	for ( i = 0; i < n; ++i ) {
		sum += *r * *c++;	r += inc_r;
		sum += *r * *c++;	r += inc_r;
		sum += *r * *c++;	r += inc_r;
		sum += *r * *c++;	r += inc_r;
		sum += *r * *c++;	r += inc_r;
		sum += *r * *c++;	r += inc_r;
	}
	return sum;
}

/*
	int decbt( TensorPtr a, TensorPtr b, TensorPtr c, int *ip );
	
	Block-tridiagonal matrix decomposition routine.
	
	Written by A. C. Hindmarsh.
	Latest revision.. november 10, 1983 (ach) 
	Reference.. ucid-30150
				Solution of block-tridiagonal systems of linear 
				algebraiequations
				A.C. Hindmarsh
				February 1977

	This C version was written by J. Goettgens (04/24/92).
	
	The input matrix contains three blocks of elements in each block-row, 
	including blocks in the [0][2] and [n-1][n-3] block positions.
	decbt uses block gauss elimination and subroutines gefa and gesl
	for solution of blocks.  Partial pivoting is done within 
	block-rows only.
	
	Note.. This version uses linpack routines gefa/gesl instead of
	dec/sol for solution of blocks, and it uses the blas routine dot
	for dot product calculations. 
	
	Input.. (The tensors know about their dimensions)
		a = m by m by n TensorPtr containing diagonal blocks. 
		b = m by m by n TensorPtr containing the super-diagonal blocks
			(in b[k] for k = 0,...,n-2) and the block in the [n-1][n-3]
			block position (in b[n-1]). 
		c = m by m by n TensorPtr containing the subdiagonal blocks 
			(in c[k] for k = 1,2,3,...,n-1) and the block in the
			[0][2] block position (in c[0]).
	   ip = integer array of length m*n for working storage.
	Note..
		m = order of each block 
			(= a->rows/cols = b->rows/cols = c->rows/cols).
		n = number of blocks in each direction of the matrix. 
			n must be 4 or more.  The complete matrix has order m*n.
			(n = a->planes = b->planes = c->planes)
	Output.. 
	a,b,c = Pointers to tensors containing the block lu decomposition
			of the input matrix. 
	   ip = m * n array of pivot information.
	Return value..
	      =  0  if no trouble occurred, or
		  = -2  invalid partition scheme (need kColumnPointers)
		  = -1  if the input value of m or n was illegal, or 
		  =  k  if a singular matrix was found in the k-th diagonal block.
	
	Use solbt to solve the associated linear system.
	
	External routines required.. dgefa and dgesl (from linpack) and 
	ddot (from the blas, or basic linear algebra package). 
*/

int decbt( TensorPtr a, TensorPtr b, TensorPtr c, int *ip )
{
	const int	m = a->rows;	/* Order of each block.							*/
	const int	mp = a->phys_rows;
	const int	n = a->planes;	/* Number of blocks in each direction of the
								   matrix.  n must be 4 or more.  The complete 
								   matrix has order m * n.						*/
	const int	nm1 = n - 1, nm2 = n - 2;
	
	int			err = 0;
	int			i, j, k, km1;
	MatrixPtr	sub_mat = NULL;
	Double		***at = a->tensor,
				***bt = b->tensor,
				***ct = c->tensor;
	
	/*	Do some error checking.													*/
	if ( m < 1 || n < 4 )	return -1;
	if ( a->partition != kColumnPointers )	return -2;
	
	/*	Process the first block-row.											*/
	sub_mat = SubMatrix( 0, a );
	if ( (err = gefa( sub_mat, ip )) != 0 )
		return 1;											/* error return		*/
	for ( j = 0; j < m; ++j ) {
		err = gesl( sub_mat, ip, bt[0][j], 0 );
		err = gesl( sub_mat, ip, ct[0][j], 0 );
	}
	
	/*	Adjust b[1][*][*].														*/
	{	register Double **mat = bt[1];
		for ( j = 0; j < m; ++j ) {
			for ( i = 0; i < m; ++i ) {
				mat[j][i] -= rc_dot( m, &ct[1][0][i], mp, ct[0][j], 1 );
			}
		}
	}
	
	/*	Main loop.  Process block-rows 2 to n-1.								*/
	{	register Double **mat;
		for ( k = 1, km1 = 0; k < nm1; ++k, ++km1 ) {
			mat = at[k];
			for ( j = 0; j < m; ++j ) {
				for ( i = 0; i < m; ++i ) {
					mat[j][i] -= rc_dot( m, &ct[k][0][i], mp, bt[km1][j], 1 );
				}
			}
			sub_mat = SubMatrix( k, a );
			if ( (err = gefa(sub_mat, ip+k*m)) != 0 )
				return k + 1; 									/* error return	*/
			for ( j = 0; j < m; ++j ) {
				err = gesl( sub_mat, ip+k*m, bt[k][j], 0 );
			}
#ifdef applec
			if (!StandAlone) SpinCursor(1);
#endif
		}	/* for ( k = 1; k < nm1; ++k )	*/
	}
	
	/*	Process last block-row and return.										*/
	{	register Double **mat = ct[nm1];
		for ( j = 0; j < m; ++j ) {
			for ( i = 0; i < m; ++i ) {
				mat[j][i] -= rc_dot( m, &bt[nm1][0][i], mp, bt[nm2-1][j], 1 );
			}
		}
		for ( j = 0, mat = at[nm1]; j < m; ++j ) {
			for ( i = 0; i < m; ++i ) {
				mat[j][i] -= rc_dot( m, &ct[nm1][0][i], mp, bt[nm2][j], 1 );
			}
		}
	}
	sub_mat = SubMatrix( nm1, a );
	if ( (err = gefa( sub_mat, ip+nm1*m )) != 0 )
		return n;										/* error return		*/
	
	return err;											/* success return	*/
}


/*
	int solbt( TensorPtr a, TensorPtr b, TensorPtr c, MatrixPtr y, int *ip );

	Solution of block-tridiagonal linear system. 
	Cefficient matrix must have been previously processed by decbt.
	a, b, c, and ip  must not have been changed since call to decbt. 
	Written by A. C. Hindmarsh.
	C version by J. Goettgens.
	
	Input.. 
	a,b,c = Pointers to m by m by n Tensors containing block lu
			decomposition of coefficient matrix from decbt.
	   ip = m by n integer array of pivot information from decbt. 
		y = Pointer to matrix of length m by n containg the right-hand
			side vector (treated as an m by n array here). 
	Output.. 
		y = solution vector, of length m * n. 
	
	External routines required.. gesl (linpack) and dot (blas). 
*/

int solbt( TensorPtr a, TensorPtr b, TensorPtr c, MatrixPtr y, int *ip )
{
	const int	n = a->planes;
	const int	m = a->rows;
	const int	mp = a->phys_rows;
	const int	nm1 = n - 1, nm2 = n - 2;
	
	int			i, k, km1;
	int			err = 0;
	Double		**mat = y->mat;
	
	/*	Do some error checking.													*/
	if ( y->partition != kColumnPointers )	return -2;
	
	/*	Forward solution sweep.													*/
	err = gesl( SubMatrix( 0, a ), ip, mat[0], 0 );
	for ( k = 1; k < nm1; ++k ) {
		km1 = k - 1;
		for ( i = 0; i < m; ++i ) {
			mat[k][i] -= rc_dot( m, &c->tensor[k][0][i], mp, mat[km1], 1 );
		}
		err = gesl( SubMatrix( k, a ), ip + k*m, mat[k], 0 );
#ifdef applec
		if (!StandAlone) SpinCursor(1);
#endif
	}
	
	for ( i = 0; i < m; ++i ) {
		mat[nm1][i] -= rc_dot( m, &c->tensor[nm1][0][i], mp, mat[nm2], 1 ) + 
			 		   rc_dot( m, &b->tensor[nm1][0][i], mp, mat[nm2-1], 1 );
	}
	err = gesl( SubMatrix( nm1, a ), ip + nm1*m, mat[nm1], 0 );
	
	/*	Backward solution sweep.												*/
	for ( k = nm2; k >= 0; --k ) {
		for ( i = 0; i < m; ++i ) {
			mat[k][i] -= rc_dot( m, &b->tensor[k][0][i], mp, mat[k+1], 1 );
		}
#ifdef applec
		if (!StandAlone) SpinCursor(1);
#endif
	}
	for ( i = 0; i < m; ++i ) {
		mat[0][i] -= rc_dot( m, &c->tensor[0][0][i], mp, mat[2], 1 );
	}
	
	return err;
}
