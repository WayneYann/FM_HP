/*
	AM_alloc.c: allocation/deallocation functions for ArrayManager

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

#ifdef qMemDebug
#define SKIPOMPARSE 1
#if SKIPOMPARSE
/*  OOps!  This file MUST not be compiled with qMemDebug defined. */
#undef qMemDebug
#endif
#undef SKIPOMPARSE
#endif

/*
	Double *New1DArray( int n )
	
	Allocate a 1-dimensional array with "n" elements.
	If memory allocation is successful, a pointer to the block
	of memory is returned, otherwise the program terminates
	with a call to FATAL. The array is initialized to zero.
*/

Double *_New1DArray( int n )
{
	return NEW2( n, Double );
}

Double *New1DArray_dbg( int n, const char *file, int line )
{
	return DBGNEW2( n, Double );
}


/*
	void Free1DArray( Double *ptr )
	
	Free the memory of a 1-dimensional array the memory of which
	was previously allocated with New1DArray().
*/

void _Free1DArray( Double *ptr )
{
	DELETE( ptr );
}

void Free1DArray_dbg( Double *ptr, const char *file, int line )
{
	DBGDELETE( ptr );
}


/*
	Double **New2DArray( int rows, int columns )

	New2DArray allocates a 2-dimensional array dimensioned  rows * columns.
	If memory allocation is successful, a pointer of type  Double **  is
	returned, otherwise the program terminates with a call to FATAL.
	The array should be indexed  a[row][col]. The array is initialized to zero.
	
	The memory layout for an array with name  a  having 10 rows and 7 columns
	is as follows:
	
		a[0] ->	* * * * * * *
		a[1] ->	* * * * * * *
		a[2] ->	* * * * * * *
		a[3] ->	* * * * * * *
		a[4] ->	* * * * * * *
		a[5] ->	* * * * * * *
		a[6] ->	* * * * * * *
		a[7] ->	* * * * * * *
		a[8] ->	* * * * * * *
		a[9] ->	* * * * * * *

	Essentially the matrix is interpreted as a collection of row vectors.
	By indexing the matrix with  a[column][row]  this matrix is interpreted
	as a collection of column vectors.
*/

Double **_New2DArray( int rows, int columns )
{
	const int kAlign = 8;
	char	*pool = NULL;
	Double	**p = NULL, *a = NULL;
	int		size, offset;
	int		i;

	offset = rows * sizeof( a );
	if ( offset % kAlign ) offset += kAlign - (offset % kAlign);
	size = rows * columns * sizeof( *a ) + offset;
	pool = NEW2( size, char );
	
	p = (Double **)pool;
	a = (Double *)(pool + offset);
	for ( i = 0; i < rows; ++i ) *(p + i) = a + i * columns;
	
	return p;
}

Double **New2DArray_dbg( int rows, int columns, const char *file, int line )
{
	const int kAlign = 8;
	char	*pool = NULL;
	Double	**p = NULL, *a = NULL;
	int		size, offset;
	int		i;

	offset = rows * sizeof( a );
	if ( offset % kAlign ) offset += kAlign - (offset % kAlign);
	size = rows * columns * sizeof( *a ) + offset;
	pool = DBGNEW2( size, char );
	
	p = (Double **)pool;
	a = (Double *)(pool + offset);
	for ( i = 0; i < rows; ++i ) *(p + i) = a + i * columns;
	
	return p;
}


/*
	void Free2DArray( Double **ptr )
	
	Free the memory of a 2-dimensional array, the memory of which
	was previously allocated with New2DArray().
*/

void _Free2DArray( Double **p )
{
	DELETE( p );
}


void Free2DArray_dbg( Double **p, const char *file, int line )
{
	DBGDELETE( p );
}



/*
	Double **Make2DArray( Double *a, int rows, int columns )
	
	Given a pointer to a sufficiently large block of memory, a, and the
	dimensions of a 2D array, rows and columns, this routine allocates
	a vector of length rows of pointers to Double, such that subsequently
	the memory which a points to can be accessed as a 2D array via the
	pointer that is returned by the function. 
*/

Double **_Make2DArray( Double *a, int rows, int columns )
{
	int		i;
	Double	**p;
	
	p = NEW2(rows, Double *);
	for ( i = 0; i < rows; ++i ) p[i] = a + i * columns;
	
	return p;
}

Double **Make2DArray_dbg( Double *a, int rows, int columns, 
	const char *file, int line )
{
	int		i;
	Double	**p;
	
	p = DBGNEW2(rows, Double *);
	for ( i = 0; i < rows; ++i ) p[i] = a + i * columns;
	
	return p;
}


/*
	void Release2DArray( Double **p )
	
	Given an array pointer p, previously allocated by Make2DArray(), this
	routine releases the memory p points to. Note that the memory of the
	array entries itself is not freed.
*/

void _Release2DArray( Double **p )
{
	DELETE( p );
}

void Release2DArray_dbg( Double **p, const char *file, int line )
{
	DBGDELETE( p );
}


/*
	Double ***New3DArray( int planes, int rows, int columns )
	
	Allocate a 3-dimensional array with dimensions columns, rows
	and planes.  The array elements are initialized to zero.
	The array should be indexed  a[plane][row][column] .
*/

Double ***_New3DArray( int planes, int rows, int columns )
{
	const int	kAlign = 8;
	int			k, i;					/* counter for planes and rows	*/
	int			offset;					/* offset of 1st array element	*/
	size_t		size;					/* total size of array			*/
	char		*pool  = NULL;			/* pointer to memory pool		*/
	Double		***a   = NULL,			/* pointer to 3D array			*/
				**p    = NULL,			/* pointer to start of planes	*/
				*elems = NULL;			/* pointer to start of elements	*/
	
	/* Allocate memory for array elements and all pointers.				*/
	offset = planes * (sizeof(Double **) + rows * sizeof(Double *));
	if ( offset % kAlign ) offset += kAlign - (offset % kAlign);
	size = planes * rows * columns * sizeof(Double) + offset;
	pool = NEW2(size, char);
	
	/* Fix pointers to start of planes and elements.					*/
	p = (Double **)(pool + planes * sizeof(Double **));
	elems = (Double *)(pool + offset);

	/* Fix pointer to 3D array.											*/
	a = (Double ***)pool;

	/* Fix pointers to planes.											*/
	for ( k = 0; k < planes; ++k ) *(a + k) = p + k * rows;

	/* Fix pointers to rows.											*/
	for ( k = 0; k < planes; ++k )
		for ( i = 0; i < rows; ++i )
			*(*(a + k) + i) = elems + columns * (k * rows + i);
	
	return a;
}

Double ***New3DArray_dbg( int planes, int rows, int columns, 
	const char *file, int line )
{
	const int	kAlign = 8;
	int			k, i;					/* counter for planes and rows	*/
	int			offset;					/* offset of 1st array element	*/
	size_t		size;					/* total size of array			*/
	char		*pool  = NULL;			/* pointer to memory pool		*/
	Double		***a   = NULL,			/* pointer to 3D array			*/
				**p    = NULL,			/* pointer to start of planes	*/
				*elems = NULL;			/* pointer to start of elements	*/
	
	/* Allocate memory for array elements and all pointers.				*/
	offset = planes * (sizeof(Double **) + rows * sizeof(Double *));
	if ( offset % kAlign ) offset += kAlign - (offset % kAlign);
	size = planes * rows * columns * sizeof(Double) + offset;
	pool = DBGNEW2(size, char);
	
	/* Fix pointers to start of planes and elements.					*/
	p = (Double **)(pool + planes * sizeof(Double **));
	elems = (Double *)(pool + offset);

	/* Fix pointer to 3D array.											*/
	a = (Double ***)pool;

	/* Fix pointers to planes.											*/
	for ( k = 0; k < planes; ++k ) *(a + k) = p + k * rows;

	/* Fix pointers to rows.											*/
	for ( k = 0; k < planes; ++k )
		for ( i = 0; i < rows; ++i )
			*(*(a + k) + i) = elems + columns * (k * rows + i);

	return a;
}



/*
	void Free3DArray( Double ***p, int planes )
	
	Free all memory of a 3D array, which was previously allocated
	by New3DArray().
*/

void _Free3DArray( Double ***a )
{
	DELETE( a );
}

void Free3DArray_dbg( Double ***a, const char *file, int line )
{
	DBGDELETE( a );
}



/*
	Double ***Make3DArray( Double *a, int planes, int rows, int columns )

	Given a pointer to a sufficiently dimensioned memory block  a, the
	dimensions of a 3D array, this functions returns (if successful)
	a pointer of type Double ***, i.e. the memory to which  a  points can
	subsequently be indexed as a 3-dimensional array.
*/

Double ***_Make3DArray( Double *v, int planes, int rows, int columns )
{
	int			k, i;					/* counter for planes and rows	*/
	size_t		size;					/* size of all pointers			*/
	char		*pool = NULL;			/* pointer to memory pool		*/
	Double		***a  = NULL,			/* pointer to 3D array			*/
				**p   = NULL;			/* pointer to start of planes	*/

	/* Allocate memory for pointers.									*/
	size = planes * (sizeof(Double **) + rows * (sizeof(Double *)));
	pool = NEW2( size, char );

	/* Fix pointers to start of planes.									*/
	p = (Double **)(pool + planes * sizeof(Double **));

	/* Fix pointer to 3D array.											*/
	a = (Double ***)pool;

	/* Fix pointers to planes.											*/
	for ( k = 0; k < planes; ++k ) *(a + k) = p + k * rows;
	
	/* Fix pointers to rows.											*/
	for ( k = 0; k < planes; ++k )
		for ( i = 0; i < rows; ++i )
			*(*(a + k) + i) = v + columns * (k * rows + i);
	
	return a;
}

Double ***Make3DArray_dbg( Double *v, int planes, int rows, int columns, 
	const char *file, int line )
{
	int			k, i;					/* counter for planes and rows	*/
	size_t		size;					/* size of all pointers			*/
	char		*pool = NULL;			/* pointer to memory pool		*/
	Double		***a  = NULL,			/* pointer to 3D array			*/
				**p   = NULL;			/* pointer to start of planes	*/

	/* Allocate memory for pointers.									*/
	size = planes * (sizeof(Double **) + rows * (sizeof(Double *)));
	pool = DBGNEW2( size, char );

	/* Fix pointers to start of planes.									*/
	p = (Double **)(pool + planes * sizeof(Double **));

	/* Fix pointer to 3D array.											*/
	a = (Double ***)pool;

	/* Fix pointers to planes.											*/
	for ( k = 0; k < planes; ++k ) *(a + k) = p + k * rows;
	
	/* Fix pointers to rows.											*/
	for ( k = 0; k < planes; ++k )
		for ( i = 0; i < rows; ++i )
			*(*(a + k) + i) = v + columns * (k * rows + i);
	
	return a;
}


/*
	void Release3DArray( Double ***p )
	
	Release the memory which was previously allocated by a call
	to Make3DArray(). The memory of the array entries is not freed.
*/

void _Release3DArray( Double ***p )
{
	DELETE( p );
}

void Release3DArray_dbg( Double ***p, const char *file, int line )
{
	DBGDELETE( p );
}


/*
	VectorPtr NewVector( int phys_len, int log_len, Double scale, int scaled )
	
	Given the physical length of a vector phys_len, this function returns a
	pointer to a Vector of physical and logical length phys_len.
*/

VectorPtr _NewVector( int phys_len )
{
	VectorPtr v = NEW(Vector);
	v->vec = New1DArray( phys_len );
	v->len = v->phys_len = phys_len;
	return v;
}

VectorPtr NewVector_dbg( int phys_len, const char *file, int line )
{
	VectorPtr v = DBGNEW(Vector);
	v->vec = New1DArray_dbg( phys_len, file, line );
	v->len = v->phys_len = phys_len;
	return v;
}


void _DisposeVector( VectorPtr v )
{
	DELETE( v->vec );
	DELETE( v );
}

void DisposeVector_dbg( VectorPtr v, const char *file, int line )
{
	DBGDELETE( v->vec );
	DBGDELETE( v );
}


MatrixPtr _NewMatrix( int phys_rows, int phys_cols, int partition )
{
	MatrixPtr m = NEW(Matrix);
	switch ( partition ) {
		case kRowPointers:
			m->mat = New2DArray( phys_rows, phys_cols );
			break;
		case kColumnPointers:
			m->mat = New2DArray( phys_cols, phys_rows );
			break;
		default:
			FATAL( "NewMatrix(): unknown partition type" );
	}
	
	m->rows = m->phys_rows = phys_rows;
	m->cols = m->phys_cols = phys_cols;
	m->partition = partition;
	
	return m;
}

MatrixPtr NewMatrix_dbg( int phys_rows, int phys_cols, int partition, 
	const char *file, int line )
{
	MatrixPtr m = DBGNEW(Matrix);
	switch ( partition ) {
		case kRowPointers:
			m->mat = New2DArray_dbg( phys_rows, phys_cols, file, line );
			break;
		case kColumnPointers:
			m->mat = New2DArray_dbg( phys_cols, phys_rows, file, line );
			break;
		default:
			FATAL( "NewMatrix(): unknown partition type" );
	}
	
	m->rows = m->phys_rows = phys_rows;
	m->cols = m->phys_cols = phys_cols;
	m->partition = partition;
	
	return m;
}


void _DisposeMatrix( MatrixPtr m )
{
	Free2DArray( m->mat );
	DELETE( m );
}

void DisposeMatrix_dbg( MatrixPtr m, const char *file, int line )
{
	Free2DArray_dbg( m->mat, file, line );
	DBGDELETE( m );
}


TensorPtr _NewTensor( int phys_planes, int phys_rows, int phys_cols,
					 int partition )
{
	TensorPtr t = NEW(Tensor);
	switch ( partition ) {
		case kRowPointers:
			t->tensor = New3DArray( phys_planes, phys_rows, phys_cols );
			break;
		case kColumnPointers:
			t->tensor = New3DArray( phys_planes, phys_cols, phys_rows );
			break;
	}
	t->rows = t->phys_rows = phys_rows;
	t->cols = t->phys_cols = phys_cols;
	t->planes = t->phys_planes = phys_planes;
	t->partition = partition;
	t->mat = NULL;					/* Don't point to anything right now.	*/
	return t;
}

TensorPtr NewTensor_dbg( int phys_planes, int phys_rows, int phys_cols,
					 int partition, const char *file, int line )
{
	TensorPtr t = DBGNEW(Tensor);
	switch ( partition ) {
		case kRowPointers:
			t->tensor = New3DArray_dbg( phys_planes, phys_rows, phys_cols, file, line );
			break;
		case kColumnPointers:
			t->tensor = New3DArray_dbg( phys_planes, phys_cols, phys_rows, file, line );
			break;
	}
	t->rows = t->phys_rows = phys_rows;
	t->cols = t->phys_cols = phys_cols;
	t->planes = t->phys_planes = phys_planes;
	t->partition = partition;
	t->mat = NULL;					/* Don't point to anything right now.	*/
	return t;
}


void _DisposeTensor( TensorPtr t )
{
	Free3DArray( t->tensor );
	DELETE( t );
}

void DisposeTensor_dbg( TensorPtr t, const char *file, int line )
{
	Free3DArray_dbg( t->tensor, file, line );
	DBGDELETE( t );
}



/*
	Integer arrays
*/

int *_New1DIntArray( int n )
{
	return NEW2(n, int);

}

int *New1DIntArray_dbg( int n, const char *file, int line )
{
	return DBGNEW2(n, int);

}


int **_New2DIntArray( int rows, int columns )
{
	char	*pool = NULL;
	int		**p = NULL, *a = NULL;
	int		i, size;

	size = rows * (columns * sizeof(*a) + sizeof(a));
	pool = NEW2(size, char);
	
	p = (int **)pool;
	a = (int *)(pool + rows * sizeof(a));
	for ( i = 0; i < rows; ++i ) *(p + i) = a + i * columns;
	
	return p;
}

int **New2DIntArray_dbg( int rows, int columns, const char *file, int line )
{
	char	*pool = NULL;
	int		**p = NULL, *a = NULL;
	int		i, size;

	size = rows * (columns * sizeof(*a) + sizeof(a));
	pool = DBGNEW2(size, char);
	
	p = (int **)pool;
	a = (int *)(pool + rows * sizeof(a));
	for ( i = 0; i < rows; ++i ) *(p + i) = a + i * columns;
	
	return p;
}


void _Free1DIntArray( int *iv )
{
	DELETE( iv );
}

void Free1DIntArray_dbg( int *iv, const char *file, int line )
{
	DBGDELETE( iv );
}


void _Free2DIntArray( int **im )
{
	DELETE( im );
}

void Free2DIntArray_dbg( int **im, const char *file, int line )
{
	DBGDELETE( im );
}


IntVectorPtr _NewIntVector( int phys_len )
{
	IntVectorPtr iv = NEW(IntVector);
	iv->vec = New1DIntArray( phys_len );
	iv->len = iv->phys_len = phys_len;
	
	return iv;
}

IntVectorPtr NewIntVector_dbg( int phys_len, const char *file, int line )
{
	IntVectorPtr iv = DBGNEW(IntVector);
	iv->vec = New1DIntArray_dbg( phys_len, file, line );
	iv->len = iv->phys_len = phys_len;
	
	return iv;
}


void _DisposeIntVector( IntVectorPtr iv )
{
	DELETE( iv->vec );
	DELETE( iv );
}

void DisposeIntVector_dbg( IntVectorPtr iv, const char *file, int line )
{
	DBGDELETE( iv->vec );
	DBGDELETE( iv );
}


IntMatrixPtr _NewIntMatrix( int phys_rows, int phys_cols, int partition )
{
	IntMatrixPtr m = NEW(IntMatrix);
	switch ( partition ) {
		case kRowPointers:
			m->mat = New2DIntArray( phys_rows, phys_cols );
			break;
		case kColumnPointers:
			m->mat = New2DIntArray( phys_cols, phys_rows );
			break;
		default:
			FATAL( "NewIntMatrix(): unknown partition type" );
	}
	
	m->rows = m->phys_rows = phys_rows;
	m->cols = m->phys_cols = phys_cols;
	m->partition = partition;
	
	return m;
}

IntMatrixPtr NewIntMatrix_dbg( int phys_rows, int phys_cols, int partition, 
	const char *file, int line )
{
	IntMatrixPtr m = DBGNEW(IntMatrix);
	switch ( partition ) {
		case kRowPointers:
			m->mat = New2DIntArray_dbg( phys_rows, phys_cols, file, line );
			break;
		case kColumnPointers:
			m->mat = New2DIntArray_dbg( phys_cols, phys_rows, file, line );
			break;
		default:
			FATAL( "NewIntMatrix(): unknown partition type" );
	}
	
	m->rows = m->phys_rows = phys_rows;
	m->cols = m->phys_cols = phys_cols;
	m->partition = partition;
	
	return m;
}


void _DisposeIntMatrix( IntMatrixPtr m )
{
	Free2DIntArray( m->mat );
	DELETE( m );
}

void DisposeIntMatrix_dbg( IntMatrixPtr m, const char *file, int line )
{
	Free2DIntArray_dbg( m->mat, file, line );
	DBGDELETE( m );
}

