  /***************************************************************
  *****************************************************************
 *******************************************************************
*****								*****
*****				BLAS				*****
*****		Basic Linear Algebra Subroutines		*****
*****	      Written in the C Programming Language.		*****
*****								*****
*****	Functions include:					*****
*****	isamax, sasum, saxpy, scopy, sdot, snrm2,		*****
*****								*****
*****   In addition a few other routines are included:		*****
*****	vexopy, vfill						*****
*****								*****
*****	If your 3M library does not have the copysign function	*****
*****   then compile this file with -DCOPYSIGN and one will be	*****
*****   be supplied.						*****
 *******************************************************************
  *****************************************************************
   ***************************************************************/

/* © Copyright 1989,90,91,92 by Josef Goettgens. All rights reserved. */

#include <math.h>
#include <float.h>
#include "ArrayManager.h"


int iamax( int n, const Double *sx, int incx )
/*
    PURPOSE
        Finds the index of element having max. absolute value.

    INPUT
	n	Number of elements to check.
	sx	Vector to be checked.
	incx	Every incx-th element is checked.

*/
{
  register int i, istmp = 0;
  register FPUType smax = 0.0;

  if ( n <= 1 ) return( istmp );

  if ( incx != 1 ) {

    /* Code for increment not equal to 1.		*/

    if ( incx < 0 ) sx = sx + ((-n+1)*incx + 1);
    istmp = 0;
    smax  = fabs( *sx );
    sx += incx;
    for ( i = 1; i < n; ++i, sx += incx ) 
      if( fabs( *sx ) > smax ) {
	istmp = i;
	smax  = fabs( *sx );
      }

  } else {
  
    /* Code for increment equal to 1.			*/
  
    smax  = fabs(*sx++);
    for ( i = 1; i < n; ++i, ++sx ) 
      if ( fabs( *sx ) > smax ) { 
        istmp = i;
        smax  = fabs( *sx );
      }
  }
  return( istmp );
}
 

Double asum(int n, const Double *sx, int incx)
/*
  PURPOSE
      Returns sum of magnitudes of Double precision SX.
      sasum = sum from 0 to n-1 of  ABS(SX(1+I*INCX))

  INPUT
    n		Number of elements to multiply.
    sx		Pointer to float vector to take abs sum of.
    incx	Storage incrament for sx.

  RETURNS
    sasum       Double variable with the result.

  WARNINGS
    This routine uses the UN*X math library function fabs().

*/
{
  register int i, m;
  register FPUType sum = 0.0;
    
    if( n <= 0 ) return sum;
    
    if ( incx == 1 ) {	 	/* Code for increments equal to 1.	*/

        /* Clean-up loop so remaining vector length is a multiple of 6.	*/
	m = n % 6;
	if ( m ) {
	    for ( i = 0; i < m; ++i ) sum += fabs( sx[i] );
	    if ( n < 6 ) return( sum );
	}
	
	for ( i = m, sx += m; i < n; i += 6 ) {
	    sum += fabs( *sx++ ) + fabs( *sx++ ) + fabs( *sx++ ) + 
	      	   fabs( *sx++ ) + fabs( *sx++ ) + fabs( *sx++ );
	}

    } else {			/* Code for increments not equal to 1.	*/
        
	m = n * incx;
	for ( i = 0; i < m; i += incx ) sum += fabs( sx[i] );
    }

    return sum;
}


void saxpy(int n, Double sa, Double *sx, int incx, Double *sy, int incy)
/*
  PURPOSE
    Vector times a scalar plus a vector.  sy = sy + sa*sx.

  INPUT
    n		Number of elements to multiply.
    sa		Scalar to multiply by (note that this is a Double).
    sx		Pointer to float vector to scale.
    incx	Storage incrament for sx.
    sy		Pointer to float vector to add.
    incy	Storage incrament for sy.

  OUTPUT
    sy		sy = sy + sa*sx
*/
{
    register int i;
    register FPUType ssa = sa;

    if( n<=0 || ssa==0.0 ) return;
    
    if( incx == incy ) {
	if( incx == 1 ) {

	    /* Both increments = 1 */
	    for( i=0; i<n; i++ ) 
	      sy[i] += ssa*sx[i];
	    return;
	}
	if( incx>0 ) {

	    /* Equal, positive, non-unit increments. */
	    for( i=0; i<n; i++,sx+=incx,sy+=incx )
	      *sy += ssa*(*sx);
	    return;
	}
    }

    /* Unequal or negative increments. */
    if( incx < 0 ) sx += ((-n+1)*incx + 1);
    if( incy < 0 ) sy += ((-n+1)*incy + 1);
    for( i=0; i<n; i++,sx+=incx,sy+=incy ) 
      *sy += ssa*(*sx);
}


void saxpyx(int n, Double sa, Double *sx, int incx, Double *sy, int incy)
/*
  PURPOSE
    Vector times a scalar plus a vector.  sx = sy + sa*sx.

  INPUT
    n		Number of elements to multiply.
    sa		Scalar to multiply by (note that this is a Double).
    sx		Pointer to float vector to scale.
    incx	Storage incrament for sx.
    sy		Pointer to float vector to add.
    incy	Storage incrament for sy.

  OUTPUT
    sx		sx = sy + sa*sx
*/
{
    register int i;
    register FPUType ssa = sa;

    if( n<=0 || ssa==0.0 ) return;
    
    if( incx == incy ) {
	if( incx == 1 ) {

	    /* Both increments = 1 */
	    for( i=0; i<n; i++ ) 
	      sx[i] = sy[i] + ssa*sx[i];
	    return;
	}
	if( incx>0 ) {

	    /* Equal, positive, non-unit increments. */
	    for( i=0; i<n; i++, sx+=incx, sy+=incx)
	      *sx = *sy + ssa*(*sx);
	    return;
	}
    }

    /* Unequal or negative increments. */
    if( incx < 0 ) sx += ((-n+1)*incx + 1);
    if( incy < 0 ) sy += ((-n+1)*incy + 1);
    for( i=0; i<n; i++,sx+=incx,sy+=incy ) 
      *sx = *sy + ssa*(*sx);
}


void copy(int n, Double *sx, int incx, Double *sy, int incy)
/*
    PURPOSE
        Copies vector sx into vector sy.
 
    INPUT
        n    Number of components to copy.
	sx   Source vector
	incx Index increment for sx.
        incy Index increment for sy.
 
    OUTPUT
        sy   Destination vector.
*/
{
    register int i;

    if( n<1  ) return;
    if( incx == incy ) {
	if( incx == 1 ) {

	    /* Both increments = 1 */
	    for( i=0; i<n; i++ )
	      sy[i] = sx[i];
	    return;
	}
	if( incx > 0 ) {

	    /* Equal, positive, non-unit increments. */
	    for( i=0; i<n; i++, sx+=incx, sy+=incx)
	      *sy = *sx;
	    return;
	}
    }

    /* Non-equal or negative increments. */
    if( incx < 0 ) sx += ((-n+1)*incx + 1);
    if( incy < 0 ) sy += ((-n+1)*incy + 1);
    for( i=0; i<n; i++,sx+=incx,sy+=incy ) 
      (*sx) = (*sy);
    return;
}


Double dot(int n, const Double *sx, int incx, const Double *sy, int incy)
/*
    PURPOSE
        Forms the dot product of a vector.

    INPUT
        n       Number of elements to sum.
        sx      Address of first element of x vector.
        incx    Incrament for the x vector.
        sy      Address of first element of y vector.
        incy    incrament for the y vector.

    OUPUT
        sdot    Dot product x and y.  Double returned
		due to `C' language features.
*/
{
    register int i;
    register FPUType temp = 0.0;

    if ( n < 1 ) return temp;
    
    if ( incx == incy ) {
    
	if ( incx == 1 ) {

	    /* Both increments = 1				*/
	    register int m = n % 6;
	    
	    if ( m ) {
	      for ( i = 0; i < m; ++i ) temp += *sx++ * *sy++;
	      if ( n < 6 ) return temp;
	    }
	    n /= 6;		/* Number of remaining loops.	*/
	    for ( i = 0; i < n; ++i ) {
/* following is poor style */
/*	      temp += *sx++ * *sy++ + *sx++ * *sy++ + *sx++ * *sy++ +*/
/*	      	    + *sx++ * *sy++ + *sx++ * *sy++ + *sx++ * *sy++;*/
	      temp += sx[0] * sy[0] + sx[1] * sy[1] + sx[2] * sy[2] +
	      	    + sx[3] * sy[3] + sx[4] * sy[4] + sx[5] * sy[5];
	    }
	    
	    return temp;
	}
	if ( incx > 0 ) {
	    /* Equal, positive, non-unit increments.		*/
	    for( i = 0; i < n; ++i, sx += incx, sy += incx)
	      temp += (*sx) * (*sy);
	    return temp;
	}
    }

    /* Unequal or negative increments.				*/
    if ( incx < 0 ) sx += ((-n+1)*incx + 1);
    if ( incy < 0 ) sy += ((-n+1)*incy + 1);
    for ( i = 0; i < n; ++i, sx += incx, sy += incy ) 
      temp += (*sx) * (*sy);
    return temp;
}


Double norm( Double *x, int n )
{
    FPUType	sum1 = 0.0, sum2 = 0.0;
    Double	*ptr1 = x + n - 1, *ptr2;
    
    if ( n % 2 ) {
	sum1 = *ptr1 * *ptr1;
	--ptr1;
    }
    ptr2 = ptr1 - 1;
    while ( ptr1 >= x ) {
	sum1 += *ptr1 * *ptr1;
	sum2 += *ptr2 * *ptr2;
	ptr1 -= 2, ptr2 -= 2;
    }
    return sqrt( sum1 + sum2 );
}



Double nrm2(int n, const Double *sx, int incx)
/*
    PURPOSE
        Computes the Euclidean norm of sx while being
	very careful of distructive underflow and overflow.

    INPUT
        n       Number of elements to use.
        sx      Address of first element of x vector.
        incx    Incrament for the x vector (>0).

    OUPUT
        snrm2   Euclidean norm of sx.  Returns Double
		due to `C' language features.
    REMARKS
        This algorithm proceeds in four steps.
	1) scan zero components.
	2) do phase 2 when component is near underflow.
*/
{
  register int i;
  static FPUType cutlo, cuthi;
  FPUType sum = 0.0e0, hitst;
  FPUType xmax;

  if( n<1 || incx<1 ) return( sum );

  /* Calculate near underflow */
  if( cutlo == 0.0 ) cutlo = sqrt( DBL_MIN/DBL_EPSILON );
  /* Calculate near  overflow */
  if( cuthi == 0.0 ) cuthi = sqrt( DBL_MAX );
  hitst = cuthi / (Double)n;
  i     = 0;

  /* Zero Sum. */
  while( *sx == 0.0 && i<n ) {
    i++;
    sx += incx;
  }
  if( i>=n ) return( sum );

START:
  if( fabs( *sx ) > cutlo ) {
    for( ; i<n; i++, sx+=incx ) {	/* Loop over elements. */
      if( fabs( *sx ) > hitst ) goto GOT_LARGE;  
      sum += (*sx) * (*sx);
    }
    sum = sqrt( sum );
    return( sum );			/* Sum completed normaly. */
  }
  else {				/* Small sum prepare for phase 2. */
    xmax  = fabs( *sx );
    sx += incx;
    i++;
    sum += 1.0;
    for( ; i<n; i++, sx+=incx ) {
      if( fabs( *sx ) > cutlo ) {	/* Got normal elem.  Rescale and process. */
	sum = (sum*xmax)*xmax;
	goto START;
      }
      if( fabs( *sx ) > xmax ) {
	sum  = 1.0 + sum*(xmax /(*sx))*(xmax /(*sx));
	xmax = fabs( *sx );
	continue;
      }
      sum += ((*sx)/xmax)*((*sx)/xmax);
    }
    return( (Double)xmax*sqrt( sum ) );
  }					/* End of small sum. */

 GOT_LARGE:
  sum  = 1.0 + (sum/(*sx))/(*sx);	/* Rescale and process. */
  xmax = fabs( *sx );
  sx   += incx;
  i++;
  for( ; i<n; i++, sx+=incx ) {
    if( fabs( *sx ) > xmax ) {
      sum  = 1.0 + sum*(xmax /(*sx))*(xmax /(*sx));
      xmax = fabs( *sx );
      continue;
    }
    sum += ((*sx)/xmax)*((*sx)/xmax);
  }
  return( (Double)xmax*sqrt( sum ) );	/* End of small sum. */
}					/* End of ---SDOT--- */


Double r1mach(void)
/***********************************************************************
 ****   This routine computes the unit roundoff for DOUBLE precision 
 ****	of the machine.  This is defined as the smallest positive 
 ****	machine number u such that  1.0 + u .ne. 1.0
 ****	Returns a Double due to `C' language features.
 **********************************************************************/
{
#ifdef DBL_EPSILON
    return DBL_EPSILON;
#else
    Double u = 1.0e0, _comp_;
 
    do {
        u *= 0.5e0;
        _comp_ = 1.0e0 + u;
    }
    while( _comp_ != 1.0e0 );
    return( (Double)u*2.0e0 );
#endif
} /*-------------------- end of function r1mach ------------------------*/


#if 0
int min0(int n, int a, int b, int c, int d, int e, int f, int g, 
	int h, int i, int j, int k, int l, int m, int o, int p)
{
    int mt;
 
    if( n < 1 || n > 15 ) return( -1 );
    mt = a;
    if( n == 1 ) return( mt );
 
    if( mt > b ) mt = b;
    if( n == 2 ) return( mt );
 
    if( mt > c ) mt = c;
    if( n == 3 ) return( mt );
 
    if( mt > d ) mt = d;
    if( n == 4 ) return( mt );
 
    if( mt > e ) mt = e;
    if( n == 5 ) return( mt );
 
    if( mt > f ) mt = f;
    if( n == 6 ) return( mt );
 
    if( mt > g ) mt = g;
    if( n == 7 ) return( mt );
 
    if( mt > h ) mt = h;
    if( n == 8 ) return( mt );
 
    if( mt > i ) mt = i;
    if( n == 9 ) return( mt );
 
    if( mt > j  ) mt = j;
    if( n == 10 ) return( mt );
 
    if( mt > k  ) mt = k;
    if( n == 11 ) return( mt );
 
    if( mt > l  ) mt = l;
    if( n == 12 ) return( mt );
 
    if( mt > m ) mt = m;
    if( n == 13 ) return( mt );
 
    if( mt > o  ) mt = o;
    if( n == 14 ) return( mt );
 
    if( mt > p  ) mt = p;
    return( mt );
}
#endif


void scal(int n, Double sa, Double *sx, int incx)
/*
    PURPOSE
        Scales a vector by a constant.
 
    INPUT
        n    Number of components to scale.
        sa   Scale value (note that this is a Double).
        sx   Vector to scale.
        incx Every incx-th element of sx will be scaled.
 
    OUTPUT
        sx   Scaled vector.
*/
{
    register int i;
    register FPUType ssa = sa;
 
    if( n < 1 ) return;

    /* Code for increment not equal to 1.*/
    if( incx != 1 ) {
	if( incx < 0 ) sx += (-n+1)*incx;
	for( i=0; i<n; i++, sx+=incx )
	  *sx *= ssa;
        return;
    }

    /*  Code for unit increment. */
    for( i=0; i<n; i++ ) 
      sx[i] *= ssa;
    return;
}


void vexopy(int n, Double *v, Double *x, Double *y, int itype)
/*
  Purpose:
    To operate on the vectors x and y.

  Input:
    n		Number of elements to scale.
    x		First operand vector.
    y		Second operand vector.
    itype	Type of operation to perform:
		itype = 1 => '+'
		itype = 2 => '-'

  Output:
    v		Result vector of x op y.
*/
{
    register int i;

    if( n<1 ) return;

    if( itype == 1 )		/* ADDITION. */
      for( i=0; i<n; i++ )
	v[i] = x[i] + y[i];
    else			/* SUBTRACTION. */
      for( i=0; i<n; i++ )
	v[i] = x[i] - y[i];
}


void vfill(int n, Double *v, Double val)
/*
  Purpose
    To fill the DOUBLE vector v with the value val.
    Make sure that if you pass a value for val that you cast
    it to Double, viz.,
		vfill( 100, vector, (Double)0.0 );
*/
{
    register int i;
    register FPUType vval = val;

    if ( n < 1 ) return;
    for ( i = 0; i < n; ++i ) *v++ = vval;
}


#ifdef COPYSIGN
Double copysign(register Double x, register Double y)
/*
  PURPOSE
     This routine is recomended by the IEEE standard 754.  Most C 
     3M libraries contain this function.  If yours does not then
     compile this file with the -DCOPYSIGN option and the following
     code will be generated.  This routine copies the sign from y
     onto x.

  INPUT
     x		Thing to have it's sign changed.
     y		Variable to supply the sign.

  RETURNS
     copysign	Returns x with the sign of y (the fortran intrinsic
     		sign(x,y)).

  MPW
     #include <sane.h>
*/
{
    if( x * y < 0.0 )	return( -x );
    else		return( x );
}
#endif
