/*
	util.c: utility routines
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FitPak.h"

#undef DEBUG


void nrerror( const char *error_text )
{
	fprintf(stderr,"\nRun-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to shell...\n");
	exit(2);
}


Double *vector( int nl, int nh )
{
	int		i, n = nh-nl+1;
	Double *v;

	v = (Double *)malloc( (unsigned) (nh-nl+1)*sizeof(Double) );
	if (!v) nrerror("allocation failure in vector()");
	for ( i = 0; i < n; ++i ) v[i] = 0.0;
	
	return v-nl;
}


int *ivector(int nl, int nh )
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}


Double **matrix( int nrl, int nrh, int ncl, int nch )
{
	int i;
	Double **m;

	m=(Double **) malloc((unsigned) (nrh-nrl+1)*sizeof(Double*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(Double *) malloc((unsigned) (nch-ncl+1)*sizeof(Double));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}



Double **submatrix( Double **a, int oldrl, int oldrh, int oldcl, 
				   int oldch, int newrl, int newcl )
{
#	ifdef applec
#	pragma unused(oldch)
#	endif

	int i,j;
	Double **m;

	m=(Double **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(Double*));
	if (!m) nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	return m;
}



void free_vector( Double *v, int nl, int nh )
{
#	ifdef applec
#	pragma unused(nh)
#	endif

	free( (v+nl) );
}


void free_ivector( int *v, int nl, int nh )
{
#	ifdef applec
#	pragma unused(nh)
#	endif

	free( (v+nl) );
}


void free_matrix( Double **m, int nrl, int nrh, int ncl, int nch )
{
#	ifdef applec
#	pragma unused(nch)
#	endif

	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}



void free_submatrix(Double **b, int nrl, int nrh, int ncl, int nch )
{
#	ifdef applec
#	pragma unused(nrh,ncl,nch)
#	endif

	free( (b+nrl) );
}



Double **convert_matrix( Double *a, int nrl, int nrh, int ncl, int nch )
{
	int i,j,nrow,ncol;
	Double **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	m = (Double **) malloc((unsigned) (nrow)*sizeof(Double*));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	return m;
}



void free_convert_matrix( Double **b, int nrl, int nrh, int ncl, int nch )
{
#	ifdef applec
#	pragma unused(nrh,ncl,nch)
#	endif

	free( (b+nrl) );
}


void fpoly( Double x, Double p[], int np )
{
	int j;

	p[1] = 1.0;
	for ( j = 2; j <= np; ++j ) p[j] = p[j-1] * x;
}


void fleg( Double x, Double pl[], int nl )
{
	int		j;
	Double	twox, f2, f1, d;

	pl[1] = 1.0;
	pl[2] = x;
	if ( nl > 2 ) {
		twox = 2.0 * x;
		f2 = x;
		d = 1.0;
		for ( j = 3; j <= nl; ++j ) {
			f1 = d;
			d += 1.0;
			f2 += twox;
			pl[j] = (f2 * pl[j-1] - f1 * pl[j-2]) / d;
		}
	}
}


DataSetPtr NewDataSet( int num_of_points, int x_dim )
{
	char			*mem_error = "NewDataSet: allocation failed";
	int				i;
	Double			*v = NULL;
	DataSetPtr		data_set = NULL;
	DataPointPtr	ptr = NULL;
	
	if ( !(data_set = (DataSetPtr)malloc( sizeof(DataSet) )) )
		nrerror( mem_error );
	
	data_set->max_items = num_of_points;
	data_set->items = num_of_points;
	data_set->dim = x_dim;

	data_set->data = (DataPointPtr)calloc( num_of_points, sizeof(DataPoint) );
	if ( !data_set->data ) nrerror( mem_error );

	v = vector( 1, num_of_points * x_dim );
	if ( !v ) nrerror( mem_error );

	/*	Fix pointers to x values for each DataPoint record and set the
		standard deviation to a default value of 1.
	*/
	for ( ptr = data_set->data, i = 0; i < num_of_points; ++i, ++ptr ) {
		ptr->x = v;
		v += x_dim;
		ptr->sigma = 1.0;
	}
			
#	ifdef DEBUG
	fprintf( stderr, "\n# DataSet @ 0x%X\n", data_set );
	fprintf( stderr, "#   max_items = %d\n", data_set->max_items );
	fprintf( stderr, "#   items = %d\n", data_set->items );
	fprintf( stderr, "#   dim = %d\n", data_set->dim );
	fprintf( stderr, "#   sizeof(DataPoint) = %d\n\n", sizeof(DataPoint) );
#	endif
	
	return data_set;
}


void FreeDataSet( DataSetPtr data_set )
{
	free_vector( data_set->data->x, 1, data_set->dim );
	free( data_set->data );
	free( data_set );
}


ModelPtr NewModel( int m, ModelType type )
{
	int			i;
	char		*mem_error = "NewModel: allocation failed";
	ModelPtr	model = NULL;
	
	if ( !(model = (ModelPtr)malloc( sizeof(Model) )) ) nrerror( mem_error );
	model->m_fit = model->m = m;
	model->chi_squared = 0.0;
	model->lambda = -1.0;
	model->type = type;
	model->a = vector( 1, model->m );
	model->uncertainties = vector( 1, model->m );
	model->list = ivector( 1, model->m );
	
	model->m_fit = model->m;	/* initially assume that we use all parameters */
	for ( i = 1; i <= m; ++i ) model->list[i] = i;
	
	return model;
}


void FreeModel( ModelPtr model )
{
	free_ivector( model->list, 1, model->m );
	free_vector( model->uncertainties, 1, model->m );
	free_vector( model->a, 1, model->m );
	free( model );
}
