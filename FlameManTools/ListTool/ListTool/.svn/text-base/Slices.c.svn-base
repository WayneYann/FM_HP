/*
	Slices.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "alligator.h"
#include "ArrayManager.h"
#include "List.h"
#include "dfsd.h"
#include "WSS.h"
#include "Slices.h"

static ListPtr	gSliceList = NULL;

typedef struct Slice {
	int		n;
	Double	*x;
	Double	*data;
	char	*name;
} Slice, *SlicePtr;


typedef struct SliceSet {
	ListPtr		sliceList;
	int			maxSlices;
	int			numSlices;
	float		*timeSteps;
	char		*path;
} SliceSet, *SliceSetPtr;


static SliceSet		gSlices;
static WSSPtr		gWSS = NULL;

static const char	*kSuffix = ".bin";


void InitSlices( char *path, int maxTimeSteps, int gridPoints )
{
	memset( &gSlices, 0, sizeof(gSlices) );
		
	gSlices.sliceList = NewList( "slices" );
	gSlices.numSlices = 0;
	gSlices.maxSlices = maxTimeSteps+1;
	gSlices.path = path;
	
	gWSS = NewWSS( 100000 + (maxTimeSteps+1)*(gridPoints+1) * sizeof(float) );

	gSlices.timeSteps = WSSPush( gWSS, gSlices.maxSlices * sizeof(gSlices.timeSteps) );
}


void AddSlice( int n, Double *x, Double *data, char *name )
{
	SlicePtr	slice = NEW( Slice );
	char		fname[255];
	
	slice->n = n;
	slice->x = x;
	slice->data = data;
	slice->name = name;
	
	if ( !gSlices.sliceList ) FATAL( "slices not initialized" );
	AddItem( gSlices.sliceList, slice, kListTail );
	
	sprintf( fname, "%s%s%s", gSlices.path, name, kSuffix );
	remove( fname );
}


static void BinWriteSlice( void *ptr, void *aux )
{
#ifdef applec
#pragma unused( aux )
#endif
	SlicePtr	hdf = (SlicePtr)ptr;
	int			i, n = hdf->n;
	char		fname[255];
	FILE		*fp = NULL;
	float		*data = WSSPush( gWSS, n * sizeof(*data) );
	float		*x = WSSPush( gWSS, n * sizeof(*x) );

	if ( !data || !x ) FATAL( "workspace exhausted" );
	
	for ( i = 0; i < n; ++i ) {
		data[i] = (float)hdf->data[i], x[i] = (float)hdf->x[i];
	}
	
	sprintf( fname, "%s%s%s", gSlices.path, hdf->name, kSuffix );	
	if ( !( fp = fopen( fname, "ab" ) ) ) FATAL( "couldn't open binary file" );
	if ( fwrite( data, sizeof(*data), n, fp ) != n ) FATAL( "couldn't write data" );
	fclose( fp );
	sprintf( fname, "%s%s_x%s", gSlices.path, hdf->name, kSuffix );	
	if ( !( fp = fopen( fname, "ab" ) ) ) FATAL( "couldn't open binary file" );
	if ( fwrite( x, sizeof(*x), n, fp ) != n ) FATAL( "couldn't write x" );
	fclose( fp );
	
	WSSPop( gWSS );
	WSSPop( gWSS );
}


int WriteSlices( Double time )
{
	gSlices.timeSteps[gSlices.numSlices++] = (float)time;
	if ( gSlices.numSlices >= gSlices.maxSlices ) {
		fputs( "# too many slices, no output.\n", stderr );
		--gSlices.numSlices;
		return -1;
	}
	ForAll( gSlices.sliceList, BinWriteSlice, NULL, kListHead );
	CheckAlligatorMemory();

	return 0;
}


static void CheckHDF( int ret, const char *routine )
{
	if ( ret < 0 ) {
		fprintf( stderr, "# HDF error: %s returned %d\n", routine, ret );
		exit( 2 );
	}
}


static void SaveHDF( int nx, int ny, float *x, float *y, float *data, 
	char *name, char *path )
{
	int32	dims[2];
	int		ret;
	char	fname[255];
	
	dims[0] = ny, dims[1] = nx;

	ret = DFSDsetdims( 2, dims );
	CheckHDF( ret, "DFSDsetdims" );
		
	ret = DFSDsetdatastrs( name, "x", "", "" );
	CheckHDF( ret, "DFSDsetdatastrs" );
	
	if ( y ) {
		ret = DFSDsetdimscale( 1, dims[0], y );
		CheckHDF( ret, "DFSDsetdimscale(1)" );
	}
	if ( x ) {
		ret = DFSDsetdimscale( 2, dims[1], x );
		CheckHDF( ret, "DFSDsetdimscale(2)" );
	}
	
	sprintf( fname, "%s%s.hdf", path, name );
	ret = DFSDputdata( fname, 2, dims, data );
	CheckHDF( ret, "DFSDputdata" );
}


static void CreateSDS( void *ptr, void *aux )
{
#ifdef applec
#pragma unused( aux )
#endif
	SlicePtr	hdf = (SlicePtr)ptr;
	const int	n = hdf->n, nData = n * gSlices.numSlices;
	int			i;
	float		*data = WSSPush( gWSS, nData * sizeof(*data) );
	float		*col = WSSPush( gWSS, n * sizeof(*data) );
	Double		*x = hdf->x;
	char		fname[255];
	FILE		*fp = NULL;
		
	if ( !data || !col ) FATAL( "work space stack exhausted" );
	
	for ( i = 0; i < n; ++i ) col[i] = (float)x[i];

	sprintf( fname, "%s%s%s", gSlices.path, hdf->name, kSuffix );
	if ( !( fp = fopen( fname, "rb" ) ) ) FATAL( "error opening binary file" );
	if ( fread( data, sizeof(*data), nData, fp ) != nData ) {
		FATAL( "error reading data" );
	}
	fclose( fp );
	remove( fname );
	/*	remove the temporary coordinate file...	*/
	sprintf( fname, "%s%s_x%s", gSlices.path, hdf->name, kSuffix );
	remove( fname );

	SaveHDF( n, gSlices.numSlices, col, gSlices.timeSteps, data, hdf->name, gSlices.path );
	
#ifdef qDebug
	fprintf( stderr, "# %d bytes of free memory in work space stack\n", WSSFreeMem( gWSS ) );
#endif	

	WSSPop( gWSS );
	WSSPop( gWSS );
}


static void CreateGnuPlotOutput( void *ptr, void *aux )
{
#ifdef applec
#pragma unused( aux )
#endif
	SlicePtr	hdf = (SlicePtr)ptr;
	const int	n = hdf->n, nData = n * gSlices.numSlices;
	int			i, j;
	float		*data = WSSPush( gWSS, nData * sizeof(*data) );
	float		*x = WSSPush( gWSS, nData * sizeof(*x) );
	char		fname[255];
	FILE		*fp = NULL;
		
	sprintf( fname, "%s%s%s", gSlices.path, hdf->name, kSuffix );
	if ( !( fp = fopen( fname, "rb" ) ) ) FATAL( "error opening binary file" );
	if ( fread( data, sizeof(*data), nData, fp ) != nData ) {
		FATAL( "error reading data" );
	}
	fclose( fp );
	remove( fname );
	sprintf( fname, "%s%s_x%s", gSlices.path, hdf->name, kSuffix );
	if ( !( fp = fopen( fname, "rb" ) ) ) FATAL( "error opening binary file" );
	if ( fread( x, sizeof(*x), nData, fp ) != nData ) {
		FATAL( "error reading data" );
	}
	fclose( fp );
	remove( fname );

	sprintf( fname, "%s%s.gp", gSlices.path, hdf->name );
	if ( !( fp = fopen( fname, "w" ) ) ) FATAL( "error opening gnuplot file" );
	fprintf( fp, "# x\ty\t%s\n", hdf->name );
	for ( j = 0; j < gSlices.numSlices; ++j ) {
		for ( i = 0; i < n; ++i ) {
			fprintf( fp, "%g\t%g\t%g\n", 
				x[i + j * n], gSlices.timeSteps[j], data[i + j * n] );
		}
		fputc( '\n', fp );
	}
	fclose( fp );

	WSSPop( gWSS );
	WSSPop( gWSS );
	
#ifdef qDebug
	fprintf( stderr, "# %d bytes of free memory in work space stack\n", WSSFreeMem( gWSS ) );
#endif	
}


void CleanupSlices( int text )
{
	if ( !gSlices.numSlices ) {		/*	nothing has been written...	*/
		fputs( "# no slices written to disk => no output\n", stderr );
	} else {
		ForAll( gSlices.sliceList, (text) ? CreateGnuPlotOutput : CreateSDS, NULL, kListHead );
	}
	
	DeleteList( gSlices.sliceList );
	WSSPop( gWSS );
	WSSTrace( gWSS, TRUE );
	DisposeWSS( gWSS );
}


