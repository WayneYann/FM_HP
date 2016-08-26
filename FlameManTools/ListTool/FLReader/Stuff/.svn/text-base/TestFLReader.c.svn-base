/*
	TestFLReader.c
	Test program for the Flamelet Reader package
	© Josef Gšttgens, Peter Terhoeven, Ian Herwono, 1993
	
	Version 1.b2		07/21/93
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "alligator.h"
#include "ArrayManager.h"
#include "List.h"

#include "FLReader.h"


typedef struct Options {
	int			change;
	char		*input;
	char		*output;
	char		*pattern;
	char		*delSymbol;
} Options, *OptionsPtr;


static OptionsPtr NewOptions( void )
{
	OptionsPtr		opts = NEW( Options );
	
	return opts;
}


static void FreeOptions( OptionsPtr opts )
{
	if ( !opts ) return;
	
	DELETE( opts->input );
	DELETE( opts->output );
	DELETE( opts->pattern );
	DELETE( opts->delSymbol );

	DELETE( opts );
}


static Double Time( void )
{
	return clock() / (Double)CLOCKS_PER_SEC;
}


static void HandleCommandLine( int argc, char *argv[], OptionsPtr opts )
{
	int				c, errflg = 0;
	extern char	*optarg;
	const char		*usage = 	"# Usage: %s -h[elp] | "
								"-i <input> [-o <output>] "
								"[-p <value to detach> "
								"-d <symbol to delete>]\n";
	
	while ( ( c = getopt( argc, argv, "chi:o:p:d:" ) ) != EOF ) {
		switch ( c ) {
			case 'i':
				opts->input = NEW2( strlen(optarg)+1, char );
				strcpy( opts->input, optarg );
				break;
			case 'o':
				opts->output = NEW2( strlen(optarg)+1, char );
				strcpy( opts->output, optarg );
				break;
			case 'p':
				opts->pattern = NEW2( strlen(optarg)+1, char );
				strcpy( opts->pattern, optarg );
				break;
			case 'd': 
				opts->delSymbol = NEW2( strlen(optarg)+1, char );
				strcpy( opts->delSymbol, optarg );
				break;
			case 'c':
				opts->change = !opts->change;
				break;
			case 'h':
				fprintf( stderr, usage, argv[0] );
				FreeOptions( opts );
				exit( 0 );
				break;
			case '?':
				errflg++;
				break;
		 }
	}

	if ( errflg || !opts->input ) {
		 fprintf(stderr, usage, argv[0] );
		 FreeOptions( opts );
		 exit( -1 );
	}
}


static void myDeleteDetachedSymbolItems( void *ptr, void *aux )
{
#	ifdef applec
#	pragma unused( aux )
#	endif
	FLRSymbolPtr	sym = (FLRSymbolPtr)ptr;
	
	switch ( sym->type ) {
		case kFLRArray:
			DisposeVector( sym->val.v );
			break;
		case kFLRString:
			DELETE( sym->val.s );
			break;
		default:
			/* nothing is allocated */
			break;
	}
}


int main( int argc, char *argv[] )
{
	int				status = 0;
	FILE			*out = NULL, *in = NULL;
	ListPtr			list = NULL;
	OptionsPtr		opts = NULL;
	Double			t;
	
	opts = NewOptions();
	HandleCommandLine( argc, argv, opts );
	
	t = Time();
	InitFLReader(); 				/* initialize the flamelet reader */

	if ( !( in = fopen( opts->input, "r" ) ) ) FATAL( "error opening input file" );

	status = ReadFlamelet( in );
	if (!status) FATAL( "error reading flamelet file" );
	
	fclose( in );
#	ifdef applec
	fsetfileinfo( opts->input, 'MPS ', 'TEXT' );
#	endif

	if ( opts->change )	{
		if ( !( status = ChangeString( "date", "MŠrz 1993" ) ) )
			fputs( "# change string failed.\n", stderr );
		else {
			char	string[128];
			if ( !(status = GetString( "date", string ) ) )
				FATAL( "GetString failed" );
			fprintf( stderr, "  date is now \"%s\"\n", string );
		}
		if ( !( status = ChangeNumber( "pressure", -1, "bar" ) ) )
			fputs( "# change number failed.\n", stderr );
		else {
			Double	val;
			if ( !(status = GetNumber( "pressure", &val, NULL ) ) )
				FATAL( "GetNumber failed" );
			fprintf( stderr, "  pressure is now %g\n", val );
		}
	}

	if ( opts->delSymbol ) {
		/*	handle delete symbol
		*/
		if ( !( status = DeleteSymbol( opts->delSymbol ) ) ) 
			fprintf( stderr, "# delete symbol \"%s\" failed.\n", opts->delSymbol );
		else {
			fprintf( stderr, "# symbol \"%s\" deleted.\n", opts->delSymbol );
		}
	}

	if ( opts->output ) {
		if ( !( out = fopen( opts->output, "w" ) ) ) FATAL( "error opening output file" );
#		ifdef applec
		fsetfileinfo( opts->output, 'MPS ', 'TEXT' );
#		endif
		WriteFlameletFile( out, kFLRAll );  /* write out the symbollist in flamelet format */
		fclose( out );
	}
	
	if ( opts->pattern ) {
		list = GetPatternedItems( opts->pattern ); 
		if ( list ) {
			ForAll( list, DetachItem, NULL, kListHead );
			fputs( "detached items: \n", stdout );
			ForAll( list, PrintSymbol, stdout, kListHead );
			ForAll( list, myDeleteDetachedSymbolItems, NULL, kListHead );

			RemoveList( list );
		} else {
			fprintf( stderr, "# no symbols matching the pattern \"%s\" found.\n", 
				opts->pattern );
		}
	}
	
	CleanupFLReader();  /* clean up the flamelet reader (delete the global lists) */

	t = Time() - t;
	fprintf( stderr, "# time to process flamelet file: %.2f seconds.\n", t );
	
	FreeOptions( opts );
	
	return 0;
}
