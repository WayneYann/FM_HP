/*
	FitPak.h: header file for FitPak functions
	
	Most of the stuff is taken from the book  Numerical Recipes in C  by
	W.H. Press, B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling.
*/

#ifndef __FitPak__
#define __FitPak__

#ifdef applec
#pragma once
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "ArrayManager.h"

#ifndef FALSE
#define FALSE	0
#define TRUE	1
#endif

/*
	Type definitions
*/

typedef struct DataPoint {
	Double		y;				/* dependant variable */
	Double		sigma;			/* standard deviation of the data point */
	Double		*x;				/* independant variable vector */
} DataPoint, *DataPointPtr;

typedef struct DataSet {
	int			max_items;		/* number of data points (physical dimension) */
	int				items;		/* number of records used (logical dimension) */
	int				dim;		/* dimension of independant variable vector */
	DataPointPtr	data;		/* pointer to data points */
} DataSet, *DataSetPtr;

typedef enum ModelType { kLinear, kNonLinear } ModelType;

typedef struct Model {
	int			m;				/* dimension of basis */
	int			m_fit;			/* number of parameters to fit (<= m) */
	int			*list;			/* list of parameters to fit */
	Double		*a;				/* vector of coefficients */
	Double		*uncertainties;	/* vector of uncertainties */
	Double		chi_squared;	/* chi squared */
	Double		lambda;			/* fudge factor (Levenberg-Marquardt method) */
	ModelType	type;			/* type of fitting function (kLinear or kNonLinear) */
} Model, *ModelPtr;


typedef void (*ModelFunction)(Double *, Double *, Double *, Double *, int);
typedef void (*BasisFunctions)(Double *, Double *, int);

/*
	Prototypes
*/

void nrerror( const char *error_text );
Double *vector( int nl, int nh );
int *ivector(int nl, int nh );
Double **matrix( int nrl, int nrh, int ncl, int nch );
Double **submatrix( Double **a, int oldrl, int oldrh, int oldcl, 
				   int oldch, int newrl, int newcl );
void free_vector( Double *v, int nl, int nh );
void free_ivector( int *v, int nl, int nh );
void free_matrix( Double **m, int nrl, int nrh, int ncl, int nch );
void free_submatrix(Double **b, int nrl, int nrh, int ncl, int nch );
Double **convert_matrix( Double *a, int nrl, int nrh, int ncl, int nch );
void free_convert_matrix( Double **b, int nrl, int nrh, int ncl, int nch );
void fpoly( Double x, Double p[], int np );
void fleg( Double x, Double pl[], int nl );

DataSetPtr NewDataSet( int num_of_points, int x_dim );
void FreeDataSet( DataSetPtr );
ModelPtr NewModel( int m, ModelType type );
void FreeModel( ModelPtr );

void gaussj( Double **a, int n, Double **b, int m );
void gaussj1( Double **a, int n, Double *b, int *iw );
void ludcmp_old( Double **a, int n, int *indx, Double *d );
void lubksb_old( Double **a, int n, int *indx, Double b[] );
void mprove( Double **a, Double **alud, int n, int indx[], Double b[], Double x[] );
void svdcmp( Double **a, int m, int n, Double *w, Double **v );
void svbksb( Double **u, Double w[], Double **v, int m, int n, Double b[], Double x[] );

void svdfit( DataSetPtr data_set, ModelPtr model, Double **u, Double **v,
			  Double w[], void (*funcs)(Double *, Double *, int) );
void svdvar( Double **v, int ma, Double w[], Double **cvm );
void mrqmin( DataSetPtr data_set, ModelPtr model, Double **covar,
			 Double **alpha, ModelFunction func );
void mrqcof( DataSetPtr data_set, ModelPtr model, Double **alpha,
			 Double beta[], ModelFunction func );
void covsrt( Double **covar, int ma, int lista[], int mfit );

Double ran0( int *idum );
Double ran1( int *idum );
Double expdev( int *idum );
Double gasdev( int *idum );
Double gamdev( int ia, int *idum );

void fpoly( Double x, Double p[], int np );
void fleg( Double x, Double pl[], int nl );

/*
	High Level interface
*/

void LinearFit( DataSetPtr data_set, ModelPtr model, BasisFunctions funcs, FILE *info );
void NonLinearFit( DataSetPtr data_set, ModelPtr model, ModelFunction func, FILE *info );

#define kMaxParameters	25

#ifdef __cplusplus
}
#endif

#endif	/* __FitPak__ */
