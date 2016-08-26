#define FOOL_SOFTBENCH( x ) 

/*
	TofZ.c: This file contains the functions that compute the Temperature
			as a function of the mixture fraction.  
	
	Revision history:
		7/22/91		first version by jg
		7/28/91		implementation of the derivative function (jg)
					wrong slope at Z = § corrected (jg)
*/


#undef DEBUG


#ifdef DEBUG
#include <stdio.h>
#endif

#include "ArrayManager.h"
#include "TofZ.h"

#ifndef FALSE
#define FALSE	0
#define	TRUE	1
#endif


/*
 *	GLOBALS (file scope)
 */

static int			tz_init = FALSE;	/* signals whether coefficients are valid	*/
static Double		c0, c1, c2, c3;		/* coefficients of 3rd degree polynomial	*/
static Double		m1, m2;				/* left and right slope 					*/
static TofZParams	theParams = { 300.0, 2100.0, 0.055, 0.0275, 0.11 };

static void compute_coefficients( void );



static void compute_coefficients( void )
{
	Double a, b, a2, b2, t1, t11, t2, t3, t4;

	/*
		Left and right slope of the original function T(Z)
	*/
	m1 = (theParams.Tmax - theParams.Tu) / theParams.Zst;
	m2 = (theParams.Tmax - theParams.Tu) / (1.0 - theParams.Zst);
	
	/*
		Set up local variables
	*/
	a = theParams.alpha;
	b = theParams.beta;
	a2 = a * a;
	b2 = b * b;
	t1 = a - b;
	t11 = 1.0 / (t1 * t1 * t1);		/* 1 / (a - b)^3 */
	t2 = 1.0 - b;
	t3 = a + 2.0 * b;
	t4 = a * b;
	
	/*
		Computation of c0
	*/
	c0  = 2 * t4 * t4 * m1;
	c0 += a2 * (t1 - 2.0 * b * t2) * m2;
	c0 *= t11;
	c0 += theParams.Tu;

	/*
		Computation of c1
	*/
	c1  = -b * (4.0 * a2 + t4 + b2) * m1;
	c1 -= a * (a2 - 6.0 * b + t4 + 4.0 * b2) * m2;
	c1 *= t11;

	/*
		Computation of c2
	*/
	c2  = 2.0 * (a2 + t4 + b2) * m1;
	c2 += (-3.0 * (a + b) + 2.0 * (a2 + t4 + b2)) * m2;
	c2 *= t11;
	
	/*
		Computation of c3
	*/
	c3  = -(a + b) * m1;
	c3 -= (a + b - 2.0) * m2;
	c3 *= t11;
	
	/*
		Finally set the flag that signals a valid set of
		coefficients to true.
	*/
	tz_init = TRUE;

#	ifdef DEBUG
	fprintf( stderr, "\nTofZ coefficients:\n" );
	fprintf( stderr, "\tm1 = %g\n", m1 );
	fprintf( stderr, "\tm2 = %g\n", m2 );
	fprintf( stderr, "\tc0 = %g\n", c0 );
	fprintf( stderr, "\tc1 = %g\n", c1 );
	fprintf( stderr, "\tc2 = %g\n", c2 );
	fprintf( stderr, "\tc3 = %g\n", c3 );
#	endif
}


void set_coefficients( TofZParamsPtr ptr )
{
	theParams = *ptr;
	compute_coefficients();

#	ifdef DEBUG
	fprintf( stderr, "\nTofZParams:\n" );
	fprintf( stderr, "\tTu    = %g K\n", theParams.Tu    );
	fprintf( stderr, "\tTmax  = %g K\n", theParams.Tmax  );
	fprintf( stderr, "\tZst   = %g\n",   theParams.Zst   );
	fprintf( stderr, "\talpha = %g\n",   theParams.alpha );
	fprintf( stderr, "\tbeta  = %g\n",   theParams.beta  );
#	endif
}


void get_coefficients( TofZParamsPtr ptr )
{
	*ptr = theParams;
}


Double temperature( Double z )
{
	Double	t;				/* holds the temperature in [K] on exit */
	
	compute_coefficients();
	
	/* patch z if out of bounds */
	if ( z < 0.0 )	z = 0.0;
	if ( z > 1.0 )	z = 1.0;
	
	if ( z < theParams.alpha )
		t = theParams.Tu + m1 * z;
	else if ( z < theParams.beta )
		t = ((c3 * z + c2) * z + c1) * z + c0;
	else
		t = theParams.Tu + m2 * (1.0 - z);
	
	return t;
}


Double dTdZ( Double z )
{
	Double	dt;				/* holds the temperature in [K] on exit */
	
	if ( !tz_init ) compute_coefficients();
	
	/* patch z if out of bounds */
	if ( z < 0.0 )	z = 0.0;
	if ( z > 1.0 )	z = 1.0;
	
	if ( z < theParams.alpha )
		dt = m1;
	else if ( z < theParams.beta )
		dt = (3.0 * c3 * z + 2.0 * c2) * z + c1;
	else
		dt = -m2;
	
	return dt;
}
