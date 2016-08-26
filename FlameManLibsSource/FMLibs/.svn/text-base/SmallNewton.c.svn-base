#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "alligator.h"
#include "ArrayManager.h"
#include "SmallNewton.h"

#undef NDEBUG
#include <assert.h>

#undef NUMDEBUG

#ifndef TRUE
#define TRUE	1
#define FALSE	0
#endif

#undef DEBUG
#ifdef DEBUG
static FILE		*gDebug = NULL;
#endif


void warning_fl( const char *errorString, const char *file, int line )
{
	fprintf( stderr, "File %s; Line %d;\t# %s.\n", file, line, errorString );
}


static Double Norm( VectorPtr vec )
{
	int		i, n = vec->len;
	Double	*x = vec->vec;
	Double	norm = 0.0, tmp;
	
	i = n;
	while ( i-- ) {
		tmp	= *x++;
		norm += tmp * tmp;
	}
	
	return sqrt( norm ) / (Double)n;
}


static int NumDeriv( const VectorPtr vec, MatrixPtr m, NewtonInfoPtr ni, void *object )
{
	int				i, j, ret, n = vec->len;
	Double			*x = vec->vec;
	Double			**jac = m->mat;
	NewtonFunc		func = ni->func;
	VectorPtr		fVec = ni->ws, f0Vec = ni->relErr;
	Double			*f = fVec->vec, *f0 = f0Vec->vec;
	const Double	kDelta = 1.e-4;
	const Double	kDeltaAdd = 1.e-12;
	const Double	kThreshold = 1.0e-20;
	Double			xOrig, delta, big, temp;
	
	if ( n != m->cols || n != m->rows ) FATAL( "incompatible matrix dimensions" );
	
	if ( ret = (*func)( vec, f0Vec, object ) ) return ret;
	
	ClearMatrix( m );
	
	for ( i = 0; i < n; ++i ) {
	  /*	  fprintf(stderr, "Newton:i %d\tx %g\n",i,x[i]);*/
	  if ( fabs( xOrig = x[i] ) < kThreshold ) /*continue*/;
		
		delta = xOrig * kDelta + kDeltaAdd;
		/*		x[i] *= ( 1.0 + kDelta );*/
		x[i] += delta;
		delta = x[i]-xOrig;
		

		if ( ret = (*func)( vec, fVec, object ) ) {
			x[i] = xOrig;			/* reset x to original values! */
			return ret;
		}
#ifdef NUMDEBUG
			  fprintf(stderr, "Newton:i %d\tx %g\tf0 %g\tf %g\n",i,x[i],f0[i],f[i]);
#endif
		x[i] = xOrig;
		
		for ( j = 0; j < n; ++j ) {
			jac[j][i] = ( f[j] - f0[j] ) / delta;
		}
	}
	
	for ( j = 0; j < n; ++j ) {		/* Loop over rows */
		big = 0.0;
#ifdef NUMDEBUG
				fprintf(stderr, "j %d\tDeriv: %g\tDelta %g\n", j, jac[j][j],delta);
#endif
		for ( i = 0; i < n; ++i ) {
			if ( ( temp=fabs( jac[j][i] ) ) > big ) big = temp;
		}
		if ( big == 0.0 ) {
			return -j;
		}
	}

	return 0;
}


void ResetNewtonScale( NewtonInfoPtr ni, Double val )
{
	int				n = ni->n;
	Double			*scale = ni->scale->vec;
	const Double	kMinScale = (val > ni->minScale) ? val : ni->minScale; 
	
#	ifdef DEBUG
	fputs( "# resetting newton scale.\n", gDebug );
#	endif
	
/*
	while ( n-- ) {
		*scale = 1.0;
		++scale;
	}
*/
	while ( n-- ) {
		*scale = kMinScale;
		++scale;
	}

}


static void NewtonDefaults( NewtonInfoPtr ni )
{
	int				n = ni->n;
	const Double	kTol = 1.0e-10, kMinScale = 1.0e-6;
	
	SetNewtonFuncs( ni, NULL, NumDeriv, Norm );
	
	ni->maxSteps = 10;
	ni->modified = 1;
	ni->converged = FALSE;
	
	ni->tol = kTol;
	ni->minScale = kMinScale;
	
	ResetNewtonScale( ni, 0.0 );
}


NewtonInfoPtr NewNewtonInfo( int n, VectorPtr xIn )
{
	NewtonInfoPtr	ni = NULL;
	
#ifdef DEBUG
	if ( !( gDebug = fopen( "Newton.diagnostics", "w" ) ) )
		FATAL( "error opening file" );
# ifdef applec
	fsetfileinfo( "Newton.diagnostics", 'MPS ', 'TEXT' );
# endif
#endif

	ni = NEW( NewtonInfo );
	
	assert( n > 0 );
	ni->n = n;
	
	ni->index = NewIntVector( n );
/*	ni->x = NewVector( n );*/
	if ( xIn ) {
		if ( xIn->len == n ) {
			ni->x = xIn;
			ni->freeX = 0;
		}
		else {
			fprintf( stderr, "#error: incompatible vector length" );
			exit(2);
		}
	}
	else {
		ni->x = NewVector( n );
		ni->freeX = 1;
	}
	ni->inc = NewVector( n );
	ni->f = NewVector( n );
	ni->ws = NewVector( n );
	ni->scale = NewVector( n );
	ni->relErr = NewVector( n );
	ni->jac = NewMatrix( n, n, kRowPointers );
	
	NewtonDefaults( ni );
	
	return ni;
}


void FreeNewtonInfo( NewtonInfoPtr ni )
{
	if ( !ni ) return;
	
	DisposeMatrix( ni->jac );
	DisposeVector( ni->relErr );
	DisposeVector( ni->scale );
	DisposeVector( ni->ws );
	DisposeVector( ni->f );
	DisposeVector( ni->inc );
/*	DisposeVector( ni->x );*/
	DisposeIntVector( ni->index );
	
#	ifdef DEBUG
	if ( gDebug ) {
		fclose( gDebug );
		gDebug = NULL;
	}
#	endif

	DELETE( ni );
}


void SetNewtonFuncs( NewtonInfoPtr ni, int (*func)( const VectorPtr vec, VectorPtr f, void *object ), NewtonDeriv deriv, 
	NewtonNorm norm )
{
	ni->func = func;
	if ( deriv ) ni->deriv = deriv; else ni->deriv = NumDeriv;
	ni->norm = ( norm ) ? norm : Norm;
}


void SetNewtonScale( NewtonInfoPtr ni, VectorPtr scale )
{
	int				i, n = scale->len;
	const Double	*scl = scale->vec;
	Double			*newtonScale = ni->scale->vec;
	const Double	kMinScale = ni->minScale;
	
	if ( n != ni->n ) {
		WARNING( "incompatible dimensions" );
		return;
	}
	
	/*	set newton scaling to scaling from Vector scale, except for 
		scale[i] < kMinScale.
	*/
	for ( i = 0; i < n; ++i ) {
		if ( scl[i] > kMinScale ) newtonScale[i] = scl[i];
		/*
		newtonScale[i] = ( scl[i] > kMinScale ) ? scl[i] : 1.0;
		*/
	}
	
#	if 0
	{
		AMPrintOptions	prnt;
		
		DefaultAMPOpts( &prnt );
		prnt.title = "newton scale";
		PrintVector( ni->scale, &prnt, stderr );
	}
#	endif
}


int NewtonSolve( NewtonInfoPtr ni, int resetScale, int newJacobian, void *object )
{
	int				i, ret, n = ni->n;
	int				step = 0, modified = ni->modified;
	int				*index = ni->index->vec;
	VectorPtr		fVec = ni->f, xVec = ni->x;
	VectorPtr		incVec = ni->inc;
	VectorPtr		relErrVec = ni->relErr;
	MatrixPtr		jacMat = ni->jac;
	Double			**jac = jacMat->mat;
	Double			*x = xVec->vec, *f = fVec->vec;
	Double			*inc = incVec->vec;
	Double			*scale = ni->scale->vec, *relErr = relErrVec->vec;
	Double			d;
	NewtonFunc		func = ni->func;
	NewtonDeriv		deriv = ni->deriv;
	NewtonNorm		norm = ni->norm;
	Double			xNorm, lastXNorm, resNorm, lastResNorm;
	const Double	kMaxResRaise = 1.0e3;
	const Double	kTol = ni->tol;

#	ifdef DEBUG
	const char		*stepType = NULL;
	fputs( "step\tinc norm\tres norm\tstep type\n", gDebug );
#	endif
	newJacobian = newJacobian;
	
	if ( !func ) FATAL( "missing function" );
	if ( modified < 1 ) modified = ni->modified = 1;
	
	/*
	SetNewtonScale( ni, ni->x );
	*/
	if ( resetScale ) ResetNewtonScale( ni, 0.0 );

	ni->converged = FALSE;
	
	xNorm = lastXNorm = resNorm = lastResNorm = 1.0;

	if ( ret = (*func)( xVec, fVec, object ) ) return ret;
	for ( i = 0; i < n; ++i ) {
		relErr[i] = f[i] / scale[i];
	}
	resNorm = (*norm)( relErrVec );
	
	while ( step++ < ni->maxSteps ) {

#		ifdef DEBUG
		stepType = "F";
#		endif
		if ( !(step % modified) || xNorm > 0.25 * lastXNorm ) {
#				ifdef DEBUG
				stepType = "C";
#				endif
				if ( ret = (*deriv)( xVec, jacMat, ni, object ) ) return ret;
				ludcmp( jac, n, index, &d, inc );
		}
		for ( i = 0; i < n; ++i ) inc[i] = -f[i];
		
		lubksb( jac, n, index, inc );

		for ( i = 0; i < n; ++i ) {
			/*
			if ( ( d = fabs( x[i] ) ) > scale[i] ) scale[i] = d;
			*/
			relErr[i] = inc[i] / scale[i];
		}
		lastXNorm = xNorm;
		xNorm = (*norm)( relErrVec );

		for ( i = 0; i < n; ++i ) x[i] += inc[i];	/* update solution */
		
#ifdef DEBUG
		fprintf( gDebug, "%d\t%e\t%e\t%s\n", step, xNorm, resNorm, stepType );
		/*
		for ( i = 0; i < n; ++i ) {
			fprintf( gDebug, "\t(%.10g\t%.10g\t%.10g\t%.10g)\n", 
				x[i], inc[i], relErr[i], scale[i] );
		}
		*/
#endif

		if ( ret = (*func)( xVec, fVec, object ) ) return ret;
		for ( i = 0; i < n; ++i ) {
			relErr[i] = f[i] / scale[i];
		}
		lastResNorm = resNorm;
		resNorm = (*norm)( relErrVec );
		
		if ( (xNorm < kTol /*&& step > 1*/ ) || resNorm == 0.0 ) {
			ni->converged = TRUE;
			/*SetNewtonScale( ni, ni->x );*/
			/*	update scaling vector;
			*/
			for ( i = 0; i < n; ++i ) {
				if ( ( d = fabs( x[i] ) ) > scale[i] ) scale[i] = d;
				/*scale[i] = fabs( x[i] );*/
			}
			break;
		}
		/*
		if (  resNorm > lastResNorm * kMaxResRaise && resNorm > 1.e-15 ) {
			fprintf( stderr, "# residual %g explodes.\n", resNorm );
			return -1;
		}
		*/		
	}
	ni->step = step;
	
	return ( ni->converged ) ? 0 : -1;
}


void PrintNewtonInfo( NewtonInfoPtr ni, FILE *fp, const char *label )
{
	if ( !fp ) fp = stdout;
	
	fprintf( fp, "\n# NewtonInfo %s:\n", (label) ? label: "" );
	fprintf( fp, "# system of %d equations.\n# max number of steps %d, modified steps %d.\n", 
		ni->n, ni->maxSteps, ni->modified );
	fprintf( fp, "# required tolerance %g, minimum scale %g\n", ni->tol, ni->minScale );
	fprintf( fp, "# %s derivative.\n", (ni->deriv) ? "user supplied" : "numerical" );
}

