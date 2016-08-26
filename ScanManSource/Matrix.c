/*
	Matrix.c
*/

#ifdef qUseDump
#pragma load "ScanMan.dump"
#else
#include "ScanMan.h"
#endif

#include "Redux.h"

void MatrixPrint( const MatrixPtr m, FILE *fp, char *label )
{
	int		i, j, row = m->rows, col = m->cols;
	Double	**a = m->mat;
	
	if ( !fp ) fp = stdout;
	
	if ( label ) fprintf( fp, "%s:\n", label );

	for ( i = 0; i < row; ++i ) {
		for ( j = 0; j < col; ++j ) {
			fprintf( fp, "\t%g", a[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fprintf( fp, "\n" );
}


static FPUType dotprod( int n, Double *x, int incx, Double *y, int incy )
{
	int		i;
	FPUType	sum = 0.0;
	
	for ( i = 0; i < n; ++i, x += incx, y += incy ) {
		sum += *x * *y;
	}
	return sum;
}


MatrixPtr Transpose( const MatrixPtr m )
{
	int			i, j, row = m->rows, col = m->cols;
	MatrixPtr	trans = NULL;
	Double		**ta = NULL, **ma = NULL;
	
	trans = NewMatrix( col, row, kRowPointers );
	
	ta = trans->mat, ma = m->mat;
	for ( i = 0; i < row; ++i ) {
		for ( j = 0; j < col; ++j ) {
			ta[j][i] = ma[i][j];
		}
	}
	
	return trans;
}


MatrixPtr MatrixMult( const MatrixPtr mat1, const MatrixPtr mat2 )
{
	int			i, j; 
	int			col1 = mat1->cols, col2 = mat2->cols;
	int			row1 = mat1->rows, row2 = mat2->rows;
	MatrixPtr	mul = NULL;
	Double		**a1 = NULL, **a2 = NULL, **a = NULL;
	
	if ( col1 != row2 ) 
		FatalError( "Illegal matrix dimensions for multiplication." );
		
	mul = NewMatrix( row1, col2, kRowPointers );
	
	a1 = mat1->mat, a2 = mat2->mat, a = mul->mat;
	for ( i = 0; i < row1; ++i ) {
		for ( j = 0; j < col2; ++j ) {
			a[i][j] = dotprod( col1, a1[i], 1, &a2[0][j], col2 );
		}
	}
	
	return mul;
}

