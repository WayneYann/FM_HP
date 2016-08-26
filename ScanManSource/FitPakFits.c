/*
	Fits.c: package to do linear and non-linear fits

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef applec
pascal unsigned long TickCount(void) = 0xA975;
#endif

#include "FitPak.h"

#undef DEBUG


void svdfit( DataSetPtr data_set, ModelPtr model, Double **u, Double **v,
			 Double w[], BasisFunctions funcs )
{
	const Double	tol = 1.0e-8;
	int				ma = model->m,
					ndata = data_set->items,
					j, i;
	Double			wmax, tmp, thresh, sum, chi2, *b;
	Double			afunc[kMaxParameters];
	DataPointPtr	ptr;
	
	/*	Do some error checking.
	*/
	if ( model->m > kMaxParameters-1 ) nrerror( "SVDFIT: Too many parameters" );
	
	/*	Obtain temporary workspace for sigma weighted function values, b[].
	*/
	b = vector(1, ndata);
	
	for ( ptr = data_set->data, i = 1; i <= ndata; ++i, ++ptr ) {
		(*funcs)( ptr->x, afunc, ma );
		tmp = 1.0 / ptr->sigma;
		for ( j = 1; j <= ma; ++j ) u[i][j] = afunc[j] * tmp;
		b[i] = ptr->y * tmp;
	}
	
	/*	Perform singular value decomposition.
	*/
	svdcmp( u, ndata, ma, w, v );
	
	/*	Edit the singular values.
	*/
	for ( wmax = 0.0, j = 1; j <= ma; ++j ) if (w[j] > wmax) wmax = w[j];
	thresh = tol * wmax;
	for (j = 1; j <= ma; ++j) if (w[j] < thresh) w[j] = 0.0;

	/*	Solve for the coefficients a[].
	*/
	svbksb( u, w, v, ndata, ma, b, model->a );
	
	/*	Compute chi^2.
	*/
	for ( ptr = data_set->data, chi2 = 0.0, i = 1; i <= ndata; ++i, ++ptr ) {
		(*funcs)( ptr->x, afunc, ma );
		for ( sum = 0.0, j = 1; j <= ma; ++j ) sum += model->a[j] * afunc[j];
		tmp = (ptr->y - sum) / ptr->sigma;
		chi2 += tmp * tmp;
	}
	model->chi_squared = chi2;
	
	/*	Clean up.	*/
	free_vector(b, 1, ndata);
}


void svdvar( Double **v, int ma, Double w[], Double **cvm )
{
	int		k, j, i;
	Double	sum, *wti;

	wti = vector( 1, ma );
	
	for ( i = 1; i <= ma; ++i ) {
		wti[i] = 0.0;
		if ( w[i] ) wti[i] = 1.0 / (w[i] * w[i]);
	}
	
	for ( i = 1; i <= ma; ++i) {
		for (j = 1; j <= i; ++j ) {
			for (sum = 0.0, k = 1; k <= ma; ++k)
				sum += v[i][k] * v[j][k] * wti[k];
			cvm[j][i] = cvm[i][j] = sum;
		}
	}
	
	free_vector(wti, 1, ma);
}



void mrqmin( DataSetPtr data_set, ModelPtr model, Double **covar,
			 Double **alpha, ModelFunction func )
{
	static Double ochisq, atry[kMaxParameters], da[kMaxParameters],
				   beta[kMaxParameters];
	static int	   *iw = NULL;
	int			   k, kk, j, ihit, *a_list = model->list;
	Double		   *a = model->a, lambda = model->lambda;
	
	if ( model->m > kMaxParameters-1 ) nrerror( "MRQMIN: Too many parameters" );

	if ( lambda < 0.0 ) {
		iw = ivector( 1, 3 * model->m_fit );
		kk = model->m_fit + 1;
		for ( j = 1; j <= model->m; ++j ) {
			ihit = 0;
			for ( k = 1; k <= model->m_fit; ++k )
				if ( a_list[k] == j ) ++ihit;
			if ( ihit == 0 )	a_list[kk++] = j;
			else if (ihit > 1)	nrerror( "Bad a_list permutation in MRQMIN-1" );
		}
		if ( kk != model->m + 1 ) nrerror( "Bad a_list permutation in MRQMIN-2" );
		model->lambda = lambda = 0.001;
		mrqcof( data_set, model, alpha, beta, func );
		ochisq = model->chi_squared;
	}

	for ( j = 1; j <= model->m_fit; ++j ) {
		for ( k = 1; k <= model->m_fit; ++k ) covar[j][k] = alpha[j][k];
		covar[j][j] = alpha[j][j] * (1.0 + lambda);
		da[j] = beta[j];
	}
	
#ifdef DEBUG
	fprintf( stdout, "\n" );
	for ( j = 1; j <= model->m_fit; ++j ) {
		for ( k = 1; k <= model->m_fit; ++k )
			fprintf( stdout, "%12g", covar[j][k] );
		fprintf( stdout, "    %12g\n", da[j] );
	}
#endif
	gaussj1( covar, model->m_fit, da, iw );
	
	if ( lambda == 0.0 ) {
		covsrt( covar, model->m, a_list, model->m_fit );
		free_ivector( iw, 1, 3 * model->m_fit );
		return;
	}
	
	for ( j = 1; j <= model->m; ++j ) atry[j] = a[j];	/* optimize!!! */
	for ( j = 1; j <= model->m_fit; ++j )
		atry[a_list[j]] = a[a_list[j]] + da[j];
	
	model->a = atry;
	mrqcof( data_set, model, covar, da, func );
	model->a = a;
	
	if ( model->chi_squared < ochisq ) {				/* step was successful */
		model->lambda *= 0.1;
		ochisq = model->chi_squared;
		for (j = 1; j <= model->m_fit; ++j ) {
			for ( k = 1; k <= model->m_fit; ++k ) alpha[j][k] = covar[j][k];
			beta[j] = da[j];
			a[a_list[j]] = atry[a_list[j]];
		}
	} else {											/* failure, increase lambda */
		model->lambda *= 10.0;
		model->chi_squared = ochisq;
	}
}


/*
	void mrqcof( Double x[], Double y[], Double sig[], int ndata, Double a[], int ma,
				 int lista[], int mfit, Double **alpha, Double beta[], Double *chisq,
				 void (*funcs)(Double, Double *, Double *, Double *, int) )
	
	Used by mrqmin() to evaluate the linearized fitting matrix **alpha and
	vector beta[]:
					2						   2   2
		     1 ¶ chi					   1  ¶ chi
		§  = - ÑÑÑÑÑÑÑ			alpha   =  - ÑÑÑÑÑÑÑÑ
		 k   2  ¶ a					 kl    2 ¶ a ¶ a
		           k							k	l
*/


void mrqcof( DataSetPtr data_set, ModelPtr model, Double **alpha,
			 Double beta[], ModelFunction func )
{
	int				k, j, i,
					ndata = data_set->items,
					mfit = model->m_fit,
					*lista = model->list;
	Double			ymod, wt, sig2i, dy, chisq,
					dyda[kMaxParameters];
	DataPointPtr	ptr = NULL;

	for ( j = 1; j <= mfit; ++j ) {
		for ( k = 1; k <= j; ++k ) alpha[j][k] = 0.0;
		beta[j] = 0.0;
	}
	for ( ptr = data_set->data, chisq = 0.0, i = 1; i <= ndata; ++i, ++ptr ) {
		(*func)( ptr->x, model->a, &ymod, dyda, model->m );
		sig2i = 1.0 / (ptr->sigma * ptr->sigma);
		dy = ptr->y - ymod;
		for ( j = 1; j <= mfit; ++j ) {
			wt = dyda[lista[j]] * sig2i;
			for ( k = 1; k <= j; ++k ) alpha[j][k] += wt * dyda[lista[k]];
			beta[j] += dy * wt;
		}
		chisq += dy * dy * sig2i;
	}
	model->chi_squared = chisq;
	for ( j = 2; j <= mfit; ++j )
		for ( k = 1; k <= j-1; ++k ) alpha[k][j] = alpha[j][k];
}



/*
	void covsrt( Double **covar, int ma, int lista[], int mfit );
	
	Given the covariance matrix **covar of a fit of ma total parameters,
	and their ordering, lista[], repack the covariance matrix to the true
	order of the parameters.  Elements associated with fixed parameters
	will be zero.
*/

void covsrt( Double **covar, int ma, int lista[], int mfit )
{
	int		i, j;
	Double	swap;

	for ( j = 1; j<ma; ++j )
		for ( i = j+1; i <= ma; ++i ) covar[i][j] = 0.0;
	for ( i = 1; i < mfit; ++i )
		for ( j = i+1; j <= mfit; ++j ) {
			if ( lista[j] > lista[i] )
				covar[lista[j]][lista[i]] = covar[i][j];
			else
				covar[lista[i]][lista[j]] = covar[i][j];
		}
	swap = covar[1][1];
	for (j = 1; j <= ma; ++j ) {
		covar[1][j] = covar[j][j];
		covar[j][j] = 0.0;
	}
	covar[lista[1]][lista[1]] = swap;
	for (j = 2; j <= mfit; ++j ) covar[lista[j]][lista[j]] = covar[1][j];
	for (j = 2; j <= ma; ++j )
		for (i = 1; i <= j-1; ++i ) covar[i][j] = covar[j][i];
}


static void progress_info( int step, ModelPtr model, FILE *info )
{
	int		i;
	
	if ( !info ) return;
	
	if ( model->type == kLinear )
		fprintf( info, "\nLinear Fit:  chi-squared: %g\n", model->chi_squared );
	else
		fprintf( info, "\nIteration #%2d:  %s %g    %s %9.2e\n", step, 
			"chi-squared:", model->chi_squared, "lambda:", model->lambda );

	for ( i = 1; i <= model->m; ++i ) fprintf( info, "        a[%2d]", i );
	fprintf( info, "\n" );
	for (i = 1; i <= model->m; ++i ) fprintf( info, " %12g", model->a[i] );
	fprintf( info, "\n" );
	if ( info == stdout ) fflush( stdout );
}

/*
	void NonLinearFit( DataSetPtr data_set, ModelPtr model, 
					   ModelFunction func, FILE *info );
	
	Levenberg-Marquardt method, attempting to reduce the value chi^2 of a fit between
	a DataSet data_set and a nonlinear function, pointed to by func, depending on 
	model->m coefficients model->a.  The array model->list numbers the parameters
	model->a such that the first model->m_fit elements correspond to values actually
	being adjusted; the remaining m - m_fit parameters are held fixed at their value.
	
	The program returns the best-fit values for the ma fit parameters a, and chi^2.
	
	The type ModelFunction is defined as
	
		typedef void (*ModelFunction)( Double *x, Double a[], Double *y,
									   Double dyda[], int n );
	
	where x points to the independant variable vector, a[] is the array of
	fitting parameters, *y points to the function value, dyda[] is  the
	array of derivatives with respect to a, and n the number of parameters
	to fit.

	This is a template for a nonlinear fitting function:

	void myFunc( Double *x, Double a[], Double *y, Double *dyda[], int n );
	{
		*y = some function of x[1], x[2], ..., and a[1], a[2], ... a[n].
		dyda[1] = the derivative of y with respect to a[1]
		dyda[2] = the derivative of y with respect to a[2]
		...
		dyda[n] = the derivative of y with respect to a[n]
	}
	
	If info is a valid pointer to a FILE, some progress information is
	printed to this stream.  If info is the NULL pointer, nothing is
	printed.
	
	See also: DataSet, Model, DataPoint, ModelFunction
	
	CAUTION!
	Currently FitPak uses the ArrayManager package, but FitPak vectors and
	matrices are not compatible with ArrayManager functions!  All arrays are
	1 based and not 0 based as in ArrayManager.
*/

void NonLinearFit( DataSetPtr data_set, ModelPtr model, ModelFunction func, FILE *info )
{
	int 	i, step = 1,				/* counts the number of iterations */
			test = 0;				/* test flag */
	Double	ochisq,					/* old value of chi^2, used to monitor convergence */
			*u = NULL,				/* alias for uncertainties */
			**covar = NULL,			/* covariance matrix */
			**alpha = NULL;			/* curvature matrix */

#ifdef applec
	int 	ticks;
	ticks = TickCount();
#endif
	covar = matrix( 1, model->m, 1, model->m );
	alpha = matrix( 1, model->m, 1, model->m );
	
	model->lambda = -1.0;
	mrqmin( data_set, model, covar, alpha, func );
	
	while ( test < 3 ) {
	
		/*	Print progress information.
		*/
		if ( info ) progress_info( step, model, info );
		
		/*	Update step and the old value of chi^2.
		*/
		++step;
		ochisq = model->chi_squared;
		
		mrqmin( data_set, model, covar, alpha, func );
		
		/*	Check stop criterion.
		*/
		if ( model->chi_squared >= ochisq )					test = 0;	/* step failed */
		else if ( fabs(ochisq - model->chi_squared) < 0.1 )	++test;
	}
	
	model->lambda = 0.0;
	mrqmin( data_set, model, covar, alpha, func );
	u = model->uncertainties;
	for ( i = 1; i <=  model->m; ++i ) u[i] =  sqrt( covar[i][i] );

#ifdef applec
	ticks = TickCount() - ticks;
#endif
	
	if ( info ) {
		progress_info( step, model, info );
		fprintf( info, "\nUncertainties:\n" );
		for ( i = 1; i <=  model->m; ++i ) fprintf( info, " %12g", u[i] );
		fprintf( info, "\n" );
#ifdef applec
		fprintf( info, "\nElapsed time: %d ticks\n", ticks );
#endif
		if ( info == stdout ) fflush( stdout );
	}
	
	/*	Clean up.
	*/
	free_matrix( alpha, 1, model->m, 1, model->m );
	free_matrix( covar, 1, model->m, 1, model->m) ;
}


/*
	void LinearFit( DataSetPtr data_set, ModelPtr model,
					BasisFunctions funcs, FILE *info );
	
	Given a data set, described by a pointer to a DataSet struct use chi^2
	minimization to determine model->m coefficients model->a of the fitting
	function y = · a[i] * basis[i].  Here we solve the fitting equations using
	singular value decomposition of the data_set->items by model->m
	matrix.  The program returns values for the model->m fit parameters model->a,
	and chi^2.  The user supplies a pointer to a function funcs that returns
	model->m basis functions evaluated at model->x in an array.
	
	The type BasisFunctions is defined as
	
		typedef void (*BasisFunctions)(Double *x, Double *basis, int n );
	
	where x is a pointer to the independant variable vector, basis a pointer
	to the vector of basis functions, and n the number of basis functions.
	
	This is a template for the function that evaluates the basis functions:

	void myBasis( Double *x, Double p[], int np )
	{
		p[1] = some functions of x[1], x[2], ...
		p[2] = ...
		...
		p[np] = ...
	}
	
	If info is a valid pointer to a FILE, some progress information is
	printed to this stream.  If info is the NULL pointer, nothing is
	printed.
	
	See also: DataSet, DataPoint, Model, BasisFunctions
	
	CAUTION!
	Currently FitPak uses the ArrayManager package, but FitPak vectors and
	matrices are not compatible with ArrayManager functions!  All arrays are
	1 based and not 0 based as in ArrayManager.
*/

void LinearFit( DataSetPtr data_set, ModelPtr model, BasisFunctions funcs, FILE *info )
{
	int			i;
	Double		**cvm = NULL, **u = NULL, **v = NULL, *w = NULL,
				*unc = model->uncertainties;

#ifdef applec
	int			ticks;
	ticks = TickCount();
#endif
	
	u = matrix( 1, data_set->items, 1, model->m );
	cvm = matrix( 1, model->m, 1, model->m );
	v = matrix( 1, model->m, 1, model->m );
	w = vector( 1, model->m );
	
	svdfit( data_set, model, u, v, w, funcs );
	svdvar( v, model->m, w, cvm );
	for ( i = 1; i <= model->m; ++i ) unc[i] = sqrt( cvm[i][i] );
	
#ifdef applec
	ticks = TickCount() - ticks;
#endif

	if ( info ) {
		progress_info( 1, model, info );
		fprintf( info, "\nUncertainties:\n" );
		for ( i = 1; i <=  model->m; ++i ) fprintf( info, " %12g", unc[i] );
		fprintf( info, "\n" );
#ifdef applec
		fprintf( info, "\nElapsed time: %d ticks\n", ticks );
#endif
		if ( info == stdout ) fflush( stdout );
	}

	/*	Clean up.
	*/
	free_vector( w, 1, model->m );
	free_matrix( v, 1, model->m, 1, model->m );
	free_matrix( cvm, 1, model->m, 1, model->m );
	free_matrix( u, 1, data_set->items, 1, model->m );
}
