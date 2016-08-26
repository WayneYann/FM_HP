/*
	AM_Printing: printing functions for the ArrayManager
	
	Revision history:
		05/01/92	removal of internal printing flags
		04/19/92	first release

	© Copyright 1989,90,91,92 by Josef Goettgens. All rights reserved.
	
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ArrayManager.h"

static void print_a_row( const Double *v, int inc, int len,
		const AMPrintOptionsPtr prnt, FILE *fp );
static void print_an_int_row( const int *v, int inc, int len, 
							const AMPrintOptionsPtr prnt, FILE *fp );
static void print_a_plane( Double **m, int rows, int cols, 
		int partition, AMPrintOptionsPtr prnt, FILE *fp );
static void print_an_int_plane( int **m, int rows, int cols, 
						int partition, AMPrintOptionsPtr prnt, FILE *fp );
static void Adorn( int theChar, int len, FILE *fp );
static void PrintTitle( AMPrintOptionsPtr prnt, FILE *fp );
static void PrintLabel( const char *label, int count, int adorn, FILE *fp );


static AMPrintOptions gAMPOpts = {
	80, 		/* line width	*/
	"%g",		/* fp format	*/
	"%d",		/* int format	*/
	" ",		/* separator	*/
	NULL,		/* title		*/
	NULL,		/* column label	*/
	"Row %d",	/* row label	*/
	"Plane %d",	/* plane label	*/
};

/*
	void AccessError( int err, int plane, int row, int col, char *label );
	
	AccesError() prints an error string according to the error code
	returned by on of the access functions  Get¥Element, or Put¥Element.
	If no error occured, AccesError() does nothing.  This function does
	not terminate the program in case of an error.
*/

void AccessError( int err, int plane, int row, int col, char *label )
{
	char	*invalid = "### %s: invalid %s %s index %d\n",
			*log = "logical",
			*phys = "physical";
			
	switch ( err ) {
		case kIndexOK:
			/* do nothing */
			break;
		case kInvalidPartition:
			fprintf( stderr, "### %s: invalid paritition type\n", label );
			break;
		case kNoData:
			fprintf( stderr, "### %s: NULL pointer detected\n", label );
			break;
		case kInvalidLogicalRowIndex:
			fprintf( stderr, invalid, label, log, "row", row );
			break;
		case kInvalidLogicalColIndex:
			fprintf( stderr, invalid, label, log, "column", col );
			break;
		case kInvalidLogicalPlaneIndex:
			fprintf( stderr, invalid, label, log, "plane", plane );
			break;
		case kInvalidPhysicalRowIndex:
			fprintf( stderr, invalid, label, phys, "row", row );
			break;
		case kInvalidPhysicalColIndex:
			fprintf( stderr, invalid, label, phys, "column", col );
			break;
		case kInvalidPhysicalPlaneIndex:
			fprintf( stderr, invalid, label, phys, "plane", plane );
			break;
		default:
			FATAL( "Invalid error code in AccessError()" );
	}
}



/*
	char *AMInfo( char * )
	
	Write information about this library into info. The required
	memory is allocated by this routine and may be released with
	a call to DELETE().
*/

char *AMInfo( char *info )
{
	
	if ( info == NULL ) info = NEW2(256, char);
	
	sprintf( info, "ArrayManager, vs.2.2.0 (%s)\n  sizeof(Double) = %d bytes\n", 
		__DATE__, (int)sizeof(Double) );
	strcat( info, "  Compiler Flags:" );
#	ifdef MAC 
	strcat( info, " MAC" );
#	endif
#	ifdef THINK_C
	strcat( info, " THINK_C" );
#	endif
#	ifdef __SC__
	strcat( info, " __SC__" );
#	endif
#	ifdef applec
	strcat( info, " applec" );
#	endif
#	ifdef macintosh
	strcat( info, " macintosh" );
#	endif
#	ifdef HP
	strcat( info, " HP");
#	endif
	strcat( info, "\n" );
	
	return info;
}


AMPrintOptionsPtr NewAMPrintOptions( void )
{
	AMPrintOptionsPtr	prnt = NEW(AMPrintOptions);
	
	DefaultAMPOpts( prnt );
	return prnt;
}


void FreeAMPrintOptions( AMPrintOptionsPtr prnt )
{
	DELETE( prnt );
}


void DefaultAMPOpts( AMPrintOptionsPtr prnt )
{
#if 0
	prnt->lineWidth = 80;
	prnt->format = "%g";
	prnt->intFormat = "%d";
	prnt->sep = " ";
	prnt->title = NULL;
	prnt->colLabel = NULL;
	prnt->rowLabel = "Row #%d";
	prnt->planeLabel = "Plane #%d";
#endif
	*prnt = gAMPOpts;
}


static void Adorn( int theChar, int len, FILE *fp )
{
	int		i;
	
	for ( i = 0; i < len; ++i ) putc( theChar, fp );
	putc( '\n', fp );
}



static void PrintTitle( AMPrintOptionsPtr prnt, FILE *fp )
{
	if ( prnt->title ) {
		putc( '\n', fp );
		fprintf( fp, "%s:\n", prnt->title );
		Adorn( '=', strlen(prnt->title) + 1, fp );
	}
}


static void PrintLabel( const char *label, int count, int adorn, FILE *fp )
{
	char	buffer[80];
	
	if ( label ) {
		sprintf( buffer, label, count );
		fputs( buffer, fp );
		putc( '\n', fp );
		if (adorn)	Adorn( adorn, strlen(buffer), fp );
	}
}


static void print_a_row( const Double *v, int inc, int len, 
							const AMPrintOptionsPtr prnt, FILE *fp )
{
	char	buffer[80], fmt[40];
	int		j, num_len, sep_len;
	int		old_len = 32767;					/* To enforce no sep. on 1st number.*/
	
	sep_len = strlen( prnt->sep );
	sprintf( fmt, "%s%s", prnt->sep, prnt->format );
	
	for ( j = 0; j < len; ++j, v += inc ) {
	
		sprintf( buffer, fmt, *v );
		num_len = strlen( buffer );
		if ( old_len + num_len <= prnt->lineWidth ) {
			old_len += num_len;
			fputs( buffer, fp );
		} else {
			if ( j ) putc( '\n', fp );			/* No new line first time around.	*/
			old_len = num_len - sep_len;		/* No separator here.				*/
			fputs( buffer + sep_len, fp );
		}
		
	}
	putc( '\n', fp );							/* Always terminate current line.	*/
}


static void print_an_int_row( const int *v, int inc, int len, 
							const AMPrintOptionsPtr prnt, FILE *fp )
{
	char	buffer[80], fmt[40];
	int		j, num_len, sep_len;
	int		old_len = 32767;					/* To enforce no sep. on 1st number.*/
	
	sep_len = strlen( prnt->sep );
	sprintf( fmt, "%s%s", prnt->sep, prnt->format );
	
	for ( j = 0; j < len; ++j, v += inc ) {
	
		sprintf( buffer, fmt, *v );
		num_len = strlen( buffer );
		if ( old_len + num_len <= prnt->lineWidth ) {
			old_len += num_len;
			fputs( buffer, fp );
		} else {
			if ( j ) putc( '\n', fp );			/* No new line first time around.	*/
			old_len = num_len - sep_len;		/* No separator here.				*/
			fputs( buffer + sep_len, fp );
		}
		
	}
	putc( '\n', fp );							/* Always terminate current line.	*/
}


void Print1DArray( Double *v, int len, AMPrintOptionsPtr prnt, FILE *fp )
{
	if ( prnt == NULL ) prnt = &gAMPOpts;
	PrintTitle( prnt, fp );	
	print_a_row( v, 1, len, prnt, fp );
}


void Print1DIntArray( int *v, int len, AMPrintOptionsPtr prnt, FILE *fp )
{
	if ( prnt == NULL ) prnt = &gAMPOpts;
	PrintTitle( prnt, fp );	
	print_an_int_row( v, 1, len, prnt, fp );
}


static void print_a_plane( Double **m, int rows, int cols, 
						int partition, AMPrintOptionsPtr prnt, FILE *fp )
{
	int		i;								/* Row counter.						*/
	int		inc;							/* Increment between row entries.	*/
	Double	*v = NULL;						/* Points to column vectors.		*/
	
	switch ( partition ) {
		case kRowPointers:		inc = 1;			break;
		case kColumnPointers:	inc = m[1] - m[0];	break;
		default:
			FATAL( "print_a_plane(): Invalid partition type detected" );
	}

	for ( i = 0; i < rows; ++i ) {			/* Print all rows.					*/
		PrintLabel( prnt->rowLabel, i, 0, fp );
		v = (partition == kRowPointers) ? *(m + i) : *m + i;
		print_a_row( v, inc, cols, prnt, fp );
	}
}


static void print_an_int_plane( int **m, int rows, int cols, 
						int partition, AMPrintOptionsPtr prnt, FILE *fp )
{
	int		i;								/* Row counter.						*/
	int		inc;							/* Increment between row entries.	*/
	int		*v = NULL;						/* Points to column vectors.		*/
	
	switch ( partition ) {
		case kRowPointers:		inc = 1;			break;
		case kColumnPointers:	inc = m[1] - m[0];	break;
		default:
			FATAL( "print_an_int_plane(): Invalid partition type detected" );
	}

	for ( i = 0; i < rows; ++i ) {			/* Print all rows.					*/
		PrintLabel( prnt->rowLabel, i, 0, fp );
		v = (partition == kRowPointers) ? *(m + i) : *m + i;
		print_an_int_row( v, inc, cols, prnt, fp );
	}
}


void Print2DArray( Double **m, int rows, int cols, int partition, 
				   AMPrintOptionsPtr prnt, FILE *fp )
{
	if ( prnt == NULL ) prnt = &gAMPOpts;
	PrintTitle( prnt, fp );
	print_a_plane( m, rows, cols, partition, prnt, fp );
}


void Print2DIntArray( int **m, int rows, int cols, int partition, 
				   AMPrintOptionsPtr prnt, FILE *fp )
{
	if ( prnt == NULL ) prnt = &gAMPOpts;
	PrintTitle( prnt, fp );
	print_an_int_plane( m, rows, cols, partition, prnt, fp );
}


void Print3DArray( Double ***t, int planes, int rows, int cols,
				int partition, AMPrintOptionsPtr prnt, FILE *fp )
{
	int		k;										/* Counter for planes.		*/
	
	if ( prnt == NULL ) prnt = &gAMPOpts;
	PrintTitle( prnt, fp );
	
	for ( k = 0; k < planes; ++k ) {				/* Print all planes.		*/
		PrintLabel( prnt->planeLabel, k, '-', fp );
		print_a_plane( *(t + k), rows, cols, partition, prnt, fp );
	}

}


void PrintVector( VectorPtr v, AMPrintOptionsPtr prnt, FILE *fp )
{
	if ( prnt == NULL ) prnt = &gAMPOpts;
	Print1DArray( v->vec, v->len, prnt, fp );
}

void PrintIntVector( IntVectorPtr v, AMPrintOptionsPtr prnt, FILE *fp )
{
	if ( prnt == NULL ) prnt = &gAMPOpts;
	Print1DIntArray( v->vec, v->len, prnt, fp );
}


void PrintMatrix( MatrixPtr m, AMPrintOptionsPtr prnt, FILE *fp )
{
	if ( prnt == NULL ) prnt = &gAMPOpts;
	Print2DArray( m->mat, m->rows, m->cols, m->partition, prnt, fp );
}

void PrintIntMatrix( IntMatrixPtr m, AMPrintOptionsPtr prnt, FILE *fp )
{
	if ( prnt == NULL ) prnt = &gAMPOpts;
	Print2DIntArray( m->mat, m->rows, m->cols, m->partition, prnt, fp );
}


void PrintTensor( TensorPtr t, AMPrintOptionsPtr prnt, FILE *fp )
{
	if ( prnt == NULL ) prnt = &gAMPOpts;
	Print3DArray( t->tensor, t->planes, t->rows, t->cols, 
				  t->partition, prnt, fp );
}

