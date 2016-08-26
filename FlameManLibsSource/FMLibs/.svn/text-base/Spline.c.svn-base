/*
	Spline.c
	4. Dezember 1991
	( pt )
	3. Oktober 1992
	external memory support and random/sequential access functions for simple splines

#	MPW Library:
	Begin
		C Spline.c -r -mc68020 -mc68881 -elems881 -s Spline
		Lib Spline.c.o -o "{myLibs}"Spline.lib
		C Spline.c -r -mc68020 -mc68881 -elems881 -s Spline -d qMemDebug -o Spline.dbg.c.o
		Lib Spline.dbg.c.o -o "{myLibs}"Spline.dbg.lib
		duplicate -y Spline.h "{myIncludes}"
	End ·· "{WorkSheet}"
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "alligator.h"
#include "ArrayManager.h"
#include "Spline.h"

#if defined (applec) || defined (powerc)
#include <CursorCtl.h>
#endif


#undef	DEBUG
#define	B_SPLINES


/* local prototypes */
static SplinePtr NewSpline( int nPoints, void *xStorage, void *storage );
static SplinePtr ComputeBSpline( Double *x, Double *f, int nPoints, Double sL, 
	Double sR, void *xStorage, void *storage );
static SplinePtr ComputeCSpline( Double *x, Double *f, int nPoints, Double sL, 
	Double sR, void *xStorage, void *storage );
static PSplinePtr NewPSpline( int nPoints, void *storage );
static Double ArcLength( PSplinePtr ps );


/*	NewSpline allocates the memory for a Spline spanning nPoints data points.
	If storage is not NULL, it is assumed to point to a sufficiently large 
	block to hold four or five Double arrays dimensioned nPoints.
	If xStorage is not NULL, no memory is allocated for the data points,
	as it is the case for PSplines which use their own set of data points.
*/
static SplinePtr NewSpline( int nPoints, void *xStorage, void *storage )
{
	SplinePtr	p;
	int			len;
	
	p = NEW( Spline );
	
	p->n = nPoints;
	p->xStorage = ( xStorage != NULL );
	p->storage = ( storage != NULL );
	len = ( xStorage ) ? 4 * nPoints : 5 * nPoints;
	
	if ( storage ) {		/* don't allocate for parametric splines */
		p->a = (Double *)storage;
		Clear1DArray( p->a, len );
	}
	else {	
		p->a = New1DArray( len );
	}
	p->b = p->a + nPoints;
	p->c = p->b + nPoints;
	p->d = p->c + nPoints;

	p->x = ( xStorage ) ? (Double *)xStorage : p->d + nPoints;

	return p;
}


/*	FreeSpline frees the memory previously allocated by a call to NewSpline,
	which is implicitly called by ComputeSimpleSpline. Call this function 
	after you are done with the Spline *p. 
*/
void FreeSpline( SplinePtr p )
{
	if ( !p ) return;
	if ( !p->storage ) Free1DArray( p->a );
	DELETE( p );
}


/*	ComputeBSpline computes a cubic spline for the nPoints data points given
	by x and f. Gradient BCs can be supplied by specifying appropriate
	values for the slope at the left end, sL, and the slope at the right end,
	sR. Alternatively, a zero curvature BC (natural spline) is applied at the
	left and/or right side, if the freeLeft and/or freeRight are set to TRUE. 
	A linear equation for the spline coefficients b[i] is solved.
*/
static SplinePtr SimpleBSpline( Double *x, Double *f, int nPoints, 
	int freeLeft, Double sL, int freeRight, Double sR, void *xStorage,
	void *storage )
{
	int			j, n = nPoints - 1;
	SplinePtr	s = NULL;
	Double		*a = NULL, *b = NULL, *c = NULL, *rhs = NULL;
	Double		hj, hjm;
	
	/*	allocate memory for the spline
	*/
	s = NewSpline( nPoints, xStorage, storage );
	
/*	Build matrix and rhs
*/
	a = s->a;		/* temporary fields for thomas algorithm */
	b = s->d;
	c = s->c;
	rhs = s->b;		/* rhs will contain the solution for the b's */
	
	if ( freeLeft ) {
		a[0] = 0.0;
		b[0] = 1.0;
		c[0] = 0.5;
		rhs[0] = 1.5 / (x[1] - x[0]) * (f[1] - f[0]);
	}
	else {
		a[0] = 0.0;
		b[0] = 1.0;
		c[0] = 0.0;
		rhs[0] = sL;
	}
	
	for ( j = 1; j < n; ++j ) {
		hj = x[j + 1] - x[j]; 
		hjm = x[j] - x[j - 1];
		b[j] = 2.0 * (hj + hjm);
		c[j] = hjm;
		a[j] = hj;
		rhs[j] = 3.0 * ( f[j] * (hj/hjm - hjm/hj) + 
			f[j+1] * hjm/hj - f[j-1] * hj/hjm );
#		if defined (applec) || defined (powerc)
		SpinCursor( 1 );
#		endif
	}
	
	if ( freeRight ) {
		a[n] = 0.5;
		b[n] = 1.0;
		c[n] = 0.0;
		rhs[n] = 1.5 / (x[n] - x[n-1]) * (f[n] - f[n-1]);
	}
	else {
		a[n] = 0.0;
		b[n] = 1.0;
		c[n] = 0.0;
		rhs[n] = sR;
	}
	
	thomas( a, b, c, rhs, nPoints );
	
	memcpy( s->x, x, nPoints * sizeof(Double) );
	memcpy( s->a, f, nPoints * sizeof(Double) );

#	ifdef DEBUG
	fprintf( stderr, "# equation for b-coefficients solved\n" );
#	endif
	
	if ( freeLeft ) 				/* natural spline */
		s->c[0] = 0.0;
	else {
		hj = x[1] - x[0];
		s->c[0] = ( 3.0 * ( f[1] - f[0] ) / hj - 
			( 2.0 * s->b[0] + s->b[1] ) ) / hj;
	}
	for ( j = 1; j < nPoints; ++j ) {
		hjm = x[j] - x[j - 1];
		s->c[j] = ( s->b[j] - s->b[j - 1] ) / hjm - s->c[j - 1];
#		if defined (applec) || defined (powerc)
		SpinCursor( 1 );
#		endif
	}
	if ( freeRight ) 
		s->c[n] = 0.0;				/* override c[n] in case of natural spline */

	for ( j = 0; j < n; ++j ) {
		hj = x[j + 1] - x[j];
		s->d[j] = ( s->c[j + 1] - s->c[j] ) / ( 3.0 * hj );
#		if defined (applec) || defined (powerc)
		SpinCursor( 1 );
#		endif
	}

	return s;
}


/*	ComputeBSpline computes a cubic spline for the nPoints data points given
	by x and f. Gradient BCs can be supplied by specifying appropriate
	values for the slope at the left end, sL, and the slope at the right end,
	sR. Alternatively, a zero curvature BC (natural spline) is applied at the
	left and/or right side, if the freeLeft and/or freeRight are set to TRUE. 
	A linear equation for the spline coefficients c[i] is solved.
*/
static SplinePtr SimpleCSpline( Double *x, Double *f, int nPoints, 
	int freeLeft, Double sL, int freeRight, Double sR, void *xStorage, 
	void *storage )
{
	int			j, n = nPoints-1;
	SplinePtr	s = NULL;
	Double		*a = NULL, *b = NULL, *c = NULL, *rhs = NULL;
	Double		hj, hjm;
	
	/*	allocate memory for the spline
	*/
	s = NewSpline( nPoints, xStorage, storage );
	
/*	Build matrix and rhs
*/
	a = s->a;		/* temporary fields for thomas algorithm */
	b = s->b;
	c = s->d;
	rhs = s->c;		/* rhs will contain the solution for the c's */
	
	hj = x[1] - x[0];

	if ( freeLeft ) {
		a[0] = 0.0;
		b[0] = 1.0;
		c[0] = 0.0;
		rhs[0] = 0.0;
	}
	else {
		a[0] = 0.0;
		b[0] = 2.0 * hj;
		c[0] = hj;
		rhs[0] = 3.0 * ( ( f[1] - f[0] ) / hj - sL );
	}

	for ( j = 1; j < n; ++j ) {
		hjm = hj;				/* = x[j] - x[j-1] */
		hj = x[j + 1] - x[j]; 
		
		a[j] = hjm;
		b[j] = 2.0 * (hj + hjm);
		c[j] = hj;
		rhs[j] = 3.0 * ( (f[j + 1] - f[j])/hj - (f[j] - f[j - 1])/hjm );
#		if defined (applec) || defined (powerc)
		SpinCursor( 1 );
#		endif
	}

	if ( freeRight ) {
		a[n] = 0.0;
		b[n] = 1.0;
		c[n] = 0.0;
		rhs[n] = 0.0;
	}
	else {
		hj = x[n] - x[n-1];
		a[n] = hj;
		b[n] = 2.0 * hj;
		c[n] = 0.0;
		rhs[n] = 3.0 * ( sR - ( f[n] - f[n-1] ) / hj );
	}
	
	thomas( a, b, c, rhs, nPoints );
	
	memcpy( s->x, x, nPoints * sizeof(Double) );
	memcpy( s->a, f, nPoints * sizeof(Double) );
	
	for ( j = nPoints - 2; j >= 0; --j ) {
		hj = x[j + 1] - x[j];
		s->b[j] = (f[j + 1] - f[j]) / hj - (s->c[j + 1] + 2.0 * s->c[j]) * hj / 3.0;
		s->d[j] = (s->c[j + 1] - s->c[j]) / ( 3.0 * hj );
#		if defined (applec) || defined (powerc)
		SpinCursor( 1 );
#		endif
	}
	
	return s;
}


/*	ComputeSimpleSpline is the user level call for the computation of a spline.
	It calls SimpleSpline with NULL as xStorage parameter, which assures 
	allocation for the data points field. If storage is NULL, memory allocation
	is handled by the package. The computed Spline data structure is returned.
	If bSpline is TRUE, ComputeBSpline() is called, ComputeCSpline() otherwise.
*/
SplinePtr ComputeSimpleSpline( Double *x, Double *f, int nPoints, 
	int freeLeft, Double sL, int freeRight, Double sR, void *storage, int bSpline )
{
	SplinePtr	p = NULL;
	
	if ( bSpline ) {
		p = SimpleBSpline( x, f, nPoints, freeLeft, sL, freeRight, sR, NULL, storage );
	}
	else {
		p = SimpleCSpline( x, f, nPoints, freeLeft, sL, freeRight, sR, NULL, storage );
	}
	return p;
}


/*	Locate employs a binary search algorithm to return the index in the 
	sorted array xx dimensioned n, for which xx[index] <= x <= x[index+1]. 
	For x out of range, either -1 or n-1 are returned, respectively.
*/
int Locate( const Double xx[], int n, Double x )
{
	int		ascnd, upper = n, middle, lower = -1;
	
	ascnd = ( xx[n-1] > xx[0] );			/* ordering scheme (ascending/descending)	*/
	
	while ( upper - lower > 1 ) {
		middle = ( upper + lower ) >> 1;	/* compute middle point						*/
		if ( ( x >= xx[middle] ) == ascnd )	/* bisection								*/
			lower = middle;
		else
			upper = middle;
	}

#	if defined (applec) || defined (powerc)
	SpinCursor( 1 );
#	endif

	return lower;
}


/*	SplineInterpolate interpolates nPoints data values f(x) at the 
	points given by x using the Spline *s. f(x) for x values out of
	the Spline's range are set to kBogus = -9.9999e99 and a warning 
	is issued. SplineInterpolate assumes the x values to be sorted 
	in ascending order. No memory is allocated for x or f, make sure
	to supply sufficient storage!

*/
void SplineInterpolate( SplinePtr s, Double *x, Double *f, int points )
{
	const Double	kBogus = -9.9999e99;
	const Double	kTiny = 1.e-12;
	int				i, j, k, nm1 = s->n - 1;
	int				done = FALSE, ascend;
	Double			dx;
	Double			*a = s->a, *b = s->b, *c = s->c, *d = s->d, *xs = s->x;
	
	ascend = ( x[points-1] > x[0] );
	k = ( ascend ) ? 0 : nm1;
	j = 0;
		
	while ( x[k] < xs[0] ) {	/* check for x out of range */
#		ifdef DEBUG
		fprintf( stderr, "### warning: point \"%e\" out of range.\n", x[k] );
#		endif
		f[k] = kBogus;
		++k;
		if ( k == points ) return;
#		if defined (applec) || defined (powerc)
		SpinCursor( 1 );
#		endif
	}
	
	for ( i = k; i < points; ++i ) {
		while ( xs[j+1] <= x[i] ) {		/* find the correct intervall */
			if ( ++j >= nm1 ) {
				if ( fabs( x[i] - xs[j] ) < kTiny ) {	/* x[i] matches last spline point ± roundoff */
					f[i++] = a[j];		/* increment i */
				}
				else {
					fprintf( stderr, "### warning: end of splines reached!\n" );
				}
				done = TRUE;
				break;
			}
		}
		if ( done ) break;
		dx = x[i] - xs[j];
		f[i] = a[j] + dx * ( b[j] + dx * ( c[j] + dx * d[j] ) );
#		if defined (applec) || defined (powerc)
		SpinCursor( 1 );
#		endif
	}
	
	for ( k = i; k < points; ++k ) {	/* set the remaining values */
#		ifdef DEBUG
		fprintf( stderr, "### warning: point \"%e\" out of range.\n", x[k] );
#		endif
		f[k] = kBogus;
#		if defined (applec) || defined (powerc)
		SpinCursor( 1 );
#		endif
	}
}


/*	SplineRandomAccess interpolates the Spline s at the point x using the 
	function Locate() to find the interval. The interpolated result is stored in
	*value. On success, the spline's coordinate index is returned, otherwise -1.
*/
int SplineRandomAccess( SplinePtr s, Double x, Double *value )
{
	int		n = s->n, loc;
	Double	*xs = s->x;
	
	loc = Locate( xs, n, x );

#	ifdef DEBUG
//	fprintf( stderr, "# at x = %g Locate() returned %d.\n", x, loc );
#	endif

	if ( loc < 0 || loc > n-1 ) {
#		ifdef DEBUG
		fprintf( stderr, "# warning: point \"%e\" out of range in "
			"SplineRandomAccess().\n", x );
#		endif
		return -1;	
	}
	
	x -= xs[loc];
	*value = s->a[loc] + x * ( s->b[loc] + x * ( s->c[loc] + x * s->d[loc] ) );

	return loc;
}


/*	SplineSequentialAccess interpolates the Spline *s at the point x. If x is outside
	of the last evaluation interval, the function Locate() is employed to find the new
	the interval. The interpolated value is stored in value. On success, the spline's
	coordinate index is returned, -1 otherwise.
*/
int SplineSequentialAccess( SplinePtr s, Double x, Double *value )
{
	static int		loc = -1;
	Double			*xs = s->x;
	
	if ( loc >= 0 ) {
		if ( ( x >= xs[loc] ) && ( x <= xs[loc+1] ) ) {
			x -= xs[loc];
			*value = s->a[loc] + x * ( s->b[loc] + x * ( s->c[loc] + x * s->d[loc] ) );
			return loc;
		}
	}
	loc = SplineRandomAccess( s, x, value );
	if ( loc < 0 ) {
#		ifdef DEBUG
		fprintf( stderr, "# warning: point \"%e\" out of range in "
			"SplineSequentialAccess().\n", x );
#		endif
	}
	else {
		x -= xs[loc];
		*value = s->a[loc] + x * ( s->b[loc] + x * ( s->c[loc] + x * s->d[loc] ) );
	}
	return loc;
}


int SplineDerivative( SplinePtr s, Double x, Double *deriv )
{
	int		n = s->n, loc;
	Double	*xs = s->x;
	
	loc = Locate( xs, n, x );

#	ifdef DEBUG
//	fprintf( stderr, "# at x = %g Locate() returned %d.\n", x, loc );
#	endif

	if ( loc < 0 || loc > n-1 ) {
#		ifdef DEBUG
		fprintf( stderr, "# warning: point \"%e\" out of range in "
			"SplineDerivative().\n", x );
#		endif
		return -1;
	}
	
	x -= xs[loc];
	*deriv = s->b[loc] + x * ( 2.0 * s->c[loc] + 3.0 * x * s->d[loc] );

	return loc;
}


/*	PrintSpline prints the Spline pointed to by s to the file specified by fp,
	using the format string format to display the spline coefficients. If format 
	is NULL, "%g" is used, if fp is NULL, stdout is used.
*/
void PrintSpline( SplinePtr s, char *format, FILE *fp )
{
	int		i, n = s->n;
	char	formatString[80] = {0};
	
	if ( format ) {
		sprintf( formatString, "%%d\t%s\t%s\t%s\t%s\t%s\n", 
			format, format, format, format, format );
	}
	else {
		sprintf( formatString, "%%d\t%%g\t%%g\t%%g\t%%g\t%%g\n" );
	}
	
	
	if ( !fp ) fp = stdout;

	fprintf( fp, "# spline spanning %d points with %s memory supply.\n", 
		s->n, (s->storage) ? "external" : "internal" );

	fprintf( fp, "i\tx[i]\ta[i]\tb[i]\tc[i]\td[i]\n" );
	for ( i = 0; i < n; ++i ) {
		fprintf( fp, formatString, i, s->x[i], s->a[i], s->b[i], s->c[i], s->d[i] );
#			if defined (applec) || defined (powerc)
			SpinCursor( 1 );
#			endif
	}
}


/*	Writes the spline sp in binary format to the stream fp
*/
void WriteSimpleSpline( SplinePtr sp, FILE *fp )
{
	const char	*error = "WriteSimpleSpline() failed.";

	if ( fwrite( sp, sizeof(*sp), 1, fp ) != 1 ) FATAL( error );

	if ( fwrite( sp->a, sizeof(*sp->a), sp->n, fp ) != sp->n ) FATAL( error );
	if ( fwrite( sp->b, sizeof(*sp->b), sp->n, fp ) != sp->n ) FATAL( error );
	if ( fwrite( sp->c, sizeof(*sp->c), sp->n, fp ) != sp->n ) FATAL( error );
	if ( fwrite( sp->d, sizeof(*sp->d), sp->n, fp ) != sp->n ) FATAL( error );
	if ( fwrite( sp->x, sizeof(*sp->x), sp->n, fp ) != sp->n ) FATAL( error );
}


/*	Reads a spline in binary format from the stream fp, and returns the 
	spline. If storage is not NULL, it is assumed to point to a sufficiently
	large block to hold at least size bytes. 
*/
SplinePtr ReadSimpleSpline( FILE *fp, void *storage, int size )
{
	int			n, sizeA;
	Spline		mySpline = {0};
	SplinePtr	ptr = NULL;
	const char	*error = "ReadSimpleSpline() failed.";
	
	if ( fread( &mySpline, sizeof(mySpline), 1, fp ) != 1 ) FATAL( error );
	
	sizeA = sizeof(*mySpline.a);
	n = mySpline.n;
	if ( size - 5 * n * sizeA < 0 ) {
		ptr = NewSpline( n, NULL, NULL );		/* allocate your own memory */
	}
	else {
		ptr = NewSpline( n, NULL, storage );	/* storage is large enough */
	}

	if ( fread( ptr->a, sizeof(*ptr->a), ptr->n, fp ) != ptr->n ) 
		FATAL( error );
	if ( fread( ptr->b, sizeof(*ptr->b), ptr->n, fp ) != ptr->n )
		FATAL( error );
	if ( fread( ptr->c, sizeof(*ptr->c), ptr->n, fp ) != ptr->n ) 
		FATAL( error );
	if ( fread( ptr->d, sizeof(*ptr->d), ptr->n, fp ) != ptr->n ) 
		FATAL( error );
	if ( fread( ptr->x, sizeof(*ptr->x), ptr->n, fp ) != ptr->n ) 
		FATAL( error );
		
	return ptr;
}


/*	NewPSpline allocates memory for a parametric spline, PSpline spanning 
	nPoints points. The memory for the individual splines is set up when
	the Splines are computed in ComputePSpline. If storage is NULL, memory
	allocation is handled by the package, otherwise storage has to be large
	enough to hold at least 9 arrays of Double dimensioned nPoints.
*/
static PSplinePtr NewPSpline( int nPoints, void *storage )
{
	PSplinePtr	ps = NULL;
	
	ps = NEW( PSpline );
	
	ps->n = nPoints;
	ps->xs = NULL;			/* these are set up by a call to SimpleSpline */
	ps->ys = NULL;
	if ( storage ) {
		ps->storage = TRUE;
		ps->u = (Double *)storage;
		Clear1DArray( ps->u, 9 * nPoints );
	}
	else {
		ps->u = New1DArray( 9 * nPoints );
	}
	
	return ps;
}


/*	FreePSpline frees the memory of the PSpline *ps which was previously 
	allocated with a call to NewPSpline or ComputePSpline
*/
void FreePSpline( PSplinePtr ps )
{
	if ( !ps ) return;
	
	if ( !ps->storage ) Free1DArray( ps->u );
	FreeSpline( ps->ys );
	FreeSpline( ps->xs );
	
	DELETE( ps );
}


/*	ArcLength employs Simpson's itegration to compute the arc length of 
	the parametric spline ps.
*/
static Double ArcLength( PSplinePtr ps )
{
	const Double	kSim[] = {1.0/3.0, 4.0/3.0, 1.0/3.0};	/* coefficients for simpson's rule */
	int				i, j, n = ps->n-1;
	SplinePtr		xs = ps->xs, ys = ps->ys;
	Double			*u = ps->u;
	Double			len = 0.0, oldLen = 0.0, dx, dy, u0, u1, du, p;
	
	for ( i = 0; i < n; ++i ) {
		u0 = u[i], u1 = u[i + 1], du = 0.5 * ( u1 - u0 );
		for ( j = 0; j < 3; ++j ) {
			p = (Double)j * du;		/* = u - u[i] */

			dx = xs->b[i] + 2.0 * p * ( xs->c[i] + 1.5 * p * xs->d[i] );
			dy = ys->b[i] + 2.0 * p * ( ys->c[i] + 1.5 * p * ys->d[i] );
			
			len += kSim[j] * sqrt( dx*dx + dy*dy ) * du;
		}
		u[i] = oldLen;	
		oldLen = len;
#		if defined (applec) || defined (powerc)
		SpinCursor( 1 );
#		endif
	}
	u[n] = len;

#	ifdef DEBUG
//	fprintf( stderr, "# arc length = %.10g\n", len );
#	endif

	return len;
}


/*	ComputePSpline allocates a PSpline data structure and fits
	a parametric spline to the two dimensional set of data points x and y.
	Starting from a fit using the accumulated chord length as a parameter, 
	the spline is recomputed using it's arc length as parameter, until the
	change in arc length between two iteration steps is less than kTol. The 
	computed PSpline data structure is returned. If NULL is passed as storage
	parameter, memory allcoation is handled by the package. Otherwise storage 
	is assumed to be large enough to hold 9 Double arrays dimensioned nPoints.
*/
PSplinePtr ComputePSpline( Double *x, Double *y, int nPoints, void *storage, int bSpline )
{
	const int		kMaxIter = 20;
	const Double	kTol = 1.e-6;
	int				i, counter = kMaxIter;
	PSplinePtr		ps = NULL;
	Double			*u = NULL;
	Double			*ws1 = NULL, *ws2 = NULL;
	Double			dx, dy;
	Double			len = 0.0, oldLen = 0.0;
	
	ps = NewPSpline( nPoints, storage );
	
	ws1 = ps->u + nPoints;
	ws2 = ps->u + 5 * nPoints;
	
	u = ps->u;
	u[0] = 0.0;
	for ( i = 1; i < nPoints; ++i ) {
		dx = x[i] - x[i - 1], dy = y[i] - y[i - 1];
		u[i] = u[i - 1] + sqrt( dx * dx + dy * dy );
#		if defined (applec) || defined (powerc)
		SpinCursor( 1 );
#		endif
	}
	len = u[nPoints-1] - u[0];
#	ifdef DEBUG
//	fprintf( stderr, "# first arc length estimate: %g.\n", len );
#	endif
	
	if ( bSpline ) {
		ps->xs = SimpleBSpline( u, x, nPoints, TRUE, 0.0, TRUE, 0.0, ps->u, ws1 );
		ps->ys = SimpleBSpline( u, y, nPoints, TRUE, 0.0, TRUE, 0.0, ps->u, ws2 );
	}
	else {
		ps->xs = SimpleCSpline( u, x, nPoints, TRUE, 0.0, TRUE, 0.0, ps->u, ws1 );
		ps->ys = SimpleCSpline( u, y, nPoints, TRUE, 0.0, TRUE, 0.0, ps->u, ws2 );
	}
	len = ArcLength( ps );

	while ( fabs( len - oldLen) > kTol * len ) {	/* recompute using the spline's arc length */
		/*	resample arc length
		*/
		
		/*	free the simple splines.
		*/
		FreeSpline( ps->xs );
		FreeSpline( ps->ys );

		if ( bSpline ) {
			ps->xs = SimpleBSpline( u, x, nPoints, TRUE, 0.0, TRUE, 0.0, ps->u, ws1 );
			ps->ys = SimpleBSpline( u, y, nPoints, TRUE, 0.0, TRUE, 0.0, ps->u, ws2 );
		}
		else {
			ps->xs = SimpleCSpline( u, x, nPoints, TRUE, 0.0, TRUE, 0.0, ps->u, ws1 );
			ps->ys = SimpleCSpline( u, y, nPoints, TRUE, 0.0, TRUE, 0.0, ps->u, ws2 );
		}

		oldLen = len;
		len = ArcLength( ps );

		if ( !--counter ) {	/* check number of iterations */
			fprintf( stderr, "### Too many iterations in ComputePSpline().\n" );
			break;
		}
	}
	ps->arcLength = len;
	
	return ps;
}


/*	PrintPSpline prints the PSpline pointed to by ps to the file fp, 
	using the format string format to display the spline coefficients. 
	If format is NULL, "%g" is used, if fp is NULL, stdout is used.
*/
void PrintPSpline( PSplinePtr ps, char *format, FILE *fp )
{
	if ( !fp ) fp = stdout;

	fprintf( fp, "# parametric spline spanning %d points with %s memory supply.\n", 
		ps->n, (ps->storage) ? "external" : "internal" );
	fprintf( fp, "# arc length = %g\n", ps->arcLength );
	
	fprintf( fp, "# x-spline:\n" );
	PrintSpline( ps->xs, format, fp );
	
	fprintf( fp, "# y-spline:\n" );
	PrintSpline( ps->ys, format, fp );
}


Double PSplineArcLength( PSplinePtr ps )
{
	return ps->arcLength;
}


/*	PSplineInterpolate interpolates a two dimensional set of data values
	using the PSpline *ps. The nPoints interpolated points are returned 
	in x and y, which are assumed to be large enough. If storage is NULL, 
	workspace is allocated and freed by PSplineInterpolate(), otherwise 
	make sure that storage is large enough to hold an array of Double
	dimensioned points. 
*/
void PSplineInterpolate( PSplinePtr ps, Double *x, Double *y, int points, void *storage )
{
	int			i;
	Double		dl, *l = (Double *)storage;
	
	if ( !storage ) l = New1DArray( points );
	else			Clear1DArray( l, points );
	
	dl = ps->arcLength / (Double)( points - 1 );
	for ( i = 0; i < points; ++i ) {
		l[i] = i * dl;
#		if defined (applec) || defined (powerc)
		SpinCursor( 1 );
#		endif
	}
	
	SplineInterpolate( ps->xs, l, x, points );
	SplineInterpolate( ps->ys, l, y, points );
	
	if ( !storage ) Free1DArray( l );
}


