/*  
		ludcmp.c		Given an n x n matrix a[0...n-1][0...n-1], this routine
						replaces it by the  LU decomposition of a rowwise
						permutation of itself. 
						
						In combination with lubksb this routine solves linear
						equations or invert a matrix.
						
						input 	-	matrix a
						output	-	a 
									indx[0...n-1] records the row permutation
									effected by the partial pivoting
									d is +1 or -1 depending on whether the
									number of row interchanges was even or odd,
									respectively, this pointer is input for
									lubksb
					
		cf.		Numerical Recipes in C


*/


#include "ArrayManager.h"
#include "NonlinearSystems.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define TINY 1.0e-20;

void myludcmp( Double **a, int n, int *indx, Double *d )
{
	int i,imax,j,k;
	Double big,dum,sum,temp;
	Double *vv;
	
	vv = (Double *)malloc ( n * sizeof(*vv) );
	memset (vv, 0, n);
	/* vv=vector(1,n); */
	*d=1.0;
	for (i = 0;i < n; i++) {
		big=0.0;
		for (j = 0;j < n; ++j)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
		vv[i]=1.0/big;
	}
	for (j = 0; j < n; ++j) {
		for (i = 0;i < j; i++) {
			sum=a[i][j];
			for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i = j;i < n; i++) {
			sum=a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k = 0;k < n; k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	/* free_vector(vv,1,n); */
	free (vv);

}	/*  ludcmp	*/

#undef TINY


/*  

			lubksb.c		Solves the set of n linear equations A x = b.
							Here a[0...n-1][0...n-1] is input, not as the
							matrix A rather as its LU decomposition, determined by
							the routine ludcmp.
							
							input:
							a				-	matrix a decompesed with ludcmp
							n				-	number of equations
							b[0...n-1]		-	right hand side of equations
							indx[0...n-1]	-	output of ludcmp
							
							output:
							a				-	not modified
							b				-	solution vector x

*/

void mylubksb( Double **a, int n, int *indx, Double *b )
{
	int i, ii = -1, ip, j;
	Double sum;

	for (i = 0;i < n; i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii >= 0)
			for (j = ii; j <= i-1; ++j) sum -= a[i][j]*b[j];
		else if (sum) ii = i;
		b[i]=sum;
	}
	for (i = n-1; i >= 0; i--) {
		sum=b[i];
		for (j = i+1; j <  n; ++j) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}	/*  lubksb	*/




void mnewt( int ntrial, Double *x, int n, Double tolx, Double tolf 
		, NewtonFunc usrfun )
{
	int k,i,*indx;
	Double errx,errf,d,*bet,**alpha;

	indx = New1DIntArray( n );
	bet = New1DArray( n );
	alpha = New2DArray( n, n );
	for (k=0;k<ntrial;k++) {
		usrfun(x,alpha,bet);
		errf=0.0;
		for (i=0;i<n;i++) errf += fabs(bet[i]);
		if (errf <= tolf) break;
		myludcmp(alpha,n,indx,&d);
		mylubksb(alpha,n,indx,bet);
		errx=0.0;
		for (i=0;i<n;i++) {
			errx += fabs(bet[i]);
			x[i] += bet[i];
		}
		if (errx <= tolx) break;
	}
	Free2DArray( alpha );
	Free1DArray( bet );
	Free1DIntArray( indx );
	return;
}

#undef FREERETURN
