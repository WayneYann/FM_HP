/*
	IndependentEQ.c: Routines that deal with linearly independent
					 reaction equations.

	The SVD code was taken from the book "Numerical Recipes in C",
	W.H. Press, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
	Cambridge University Press (1988)
*/


#ifdef qUseDump
#pragma load "ScanMan.dump"
#else
#include "ScanMan.h"
#endif

#include "svd.h"
#include "Redux.h"


int IndependentEQ( MatrixPtr m )
{
	int			i, j, rows = m->rows, cols = m->cols, newRows;
	int			indepEQ = 0;
	double		**a = m->mat;
	double		**u = NULL, **v = NULL, *w = NULL, *vr1 = NULL;
	
	/*	Make sure we have at least as much rows as columns.						*/
	if ( cols > rows ) newRows = cols;
	else newRows = rows;

	u = New2DArray( newRows, cols );
	v = New2DArray( newRows, cols );
	w = New1DArray( newRows );

	for ( i = 0; i < rows; ++i ) {
		for ( j = 0; j < cols; ++j ) {
			u[i][j] = a[i][j];
		}
	}
	
	svdcmp( u, newRows, cols, w, v, NULL );
	
	indepEQ = svrank( w, newRows, kTiny );
	
	Free1DArray( w );
	Free2DArray( v );
	Free2DArray( u );
	
	return indepEQ;
}


MatrixPtr Inverse( MatrixPtr m )
{
	int			i, j, rows = m->rows, cols = m->cols, newRows;
	int			indepEQ = 0;
	double		**a = m->mat;
	double		**u = NULL, **v = NULL, *w = NULL, *vr1 = NULL;
	MatrixPtr	mInv = NULL;
	
	/*	Make sure we have at least as much rows as columns.						*/
	if ( cols > rows ) newRows = cols;
	else newRows = rows;

	u = New2DArray( newRows, cols );
	v = New2DArray( newRows, cols );
	w = New1DArray( newRows );

	for ( i = 0; i < rows; ++i ) {
		for ( j = 0; j < cols; ++j ) {
			u[i][j] = a[i][j];
		}
	}
		
	svdcmp( u, newRows, cols, w, v, NULL );
	
	if ( ( indepEQ = svrank( w, newRows, kTiny ) ) != newRows )	
		FatalError( "Unexpected singular matrix." );
			
	mInv = NewMatrix( cols, rows, kRowPointers );
	
	svinverse( mInv->mat, cols, u, w, v );

	Free1DArray( w );
	Free2DArray( v );
	Free2DArray( u );
	
	return mInv;
}


