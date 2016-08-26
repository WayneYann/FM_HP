/*
	ListTool.c
	Flamelet Listing Tool based on the Flamelet Reader package
	© Josef Gšttgens, Peter Terhoeven, Ian Herwono, 1993
	
	version 1.1§		4/13/93
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>

#include "alligator.h"
#include "ArrayManager.h"
#include "List.h"

#include "FLReader.h"

#ifdef applec
#include <CursorCtl.h>
#endif


#define PUT_NL(xx)		fputc('\n',xx)

typedef struct Options {
	char			*suffix;		/* output file suffix									*/
	char			*repFile;		/* report file: max/min values of all Vectors			*/
	char			*symFile;		/* check for occurence of symbols in symFile in header	*/
	char			*symOutFile;	/* output file for header symbols matching symbol file	*/
	char			*coordOutFile;	/* output file for fixed coordinate output file			*/
	char			*coordName;		/* name of coordinate									*/
	int				progress;		/* some progress information							*/
	int 			printHeader;	/* print the header to stdout							*/
	int 			printTrailer;	/* print the trailer to stdout							*/
	int				help;			/* print help information to stdout						*/
	int 			time;			/* print timing information								*/
	int				slicesAsText;	/* print slices as text									*/
	int				stripMF;		/* strip leading "massfraction"/"molefraction"			*/
	unsigned int	creator, type;	/* creator/type information used by fsetfileinfo()		*/
	char			creatorType[9];	/* creator/type for output file (e.g. "MPS TEXT")		*/
	Double			loc;			/* location for fixed coordinate output					*/
	ListPtr			flameletFiles;	/* list of flamelet files to process					*/
	ListPtr			symbolList;		/* list of symbols to process */
	ListPtr			reportSymbolList;		/* list of symbols to report							*/
} Options, *OptionsPtr;


typedef struct ReportSymbol {
	char	*name;
	Double	*data;
} ReportSymbol, *ReportSymbolPtr;


static const char	kDefaultCreatorType[] = "QKPTCGTX";
static const char	kSuffix[] = ".kg";
static VectorPtr	gZ = NULL;
static VectorPtr	gData = NULL;

#define	kBufLen			255
#define	Progress( s ) 	if (opts->progress) fprintf( stderr, "# %s...\n",(s) )


static Double Time( void )
{
/*	return clock() / (Double)CLOCKS_PER_SEC;*/
	return 1;
}


static void CopyVector( VectorPtr dest, VectorPtr src )
{
	if ( dest->phys_len < src->len ) FATAL( "dest too small" );
	
	dest->len = src->len;
	memcpy( dest->vec, src->vec, src->len*sizeof(*dest->vec) );
}


static void Warning (const char *format, ... )
{
	va_list ap;
	int result;

	fputs( "# Warning: ", stderr );
	va_start( ap, format );
	result = vfprintf( stderr, format, ap );
	va_end( ap );
	fputs( ".\n", stderr );
}


static void StripID( char *name, const char *id )
{
	const char	*strings[] = {	"massfraction-",
								"molefraction-",
								"molarfraction-",
								"concentration-",
								NULL
							};
	const char	*replace[] = {	"Y-",
								"X-",
								"X-",
								"C-",
								NULL
							};
	int			i, len, saveChar;
	const char	**ptr = strings;
	
	for ( i = 0; strings[i]; ++i ) { 
		len = strlen( strings[i] );
		saveChar = name[len];
		name[len] = '\0';
		if ( FLRCompare( name, strings[i] ) == 0 ) {
			sprintf( name, "%s%s", replace[i], id + len );
			break;
		} else {
			name[len] = saveChar;
		}
	}
	
}


static int IsFlameletFile(FILE *fp)
{
	/*	int IsFlameletFile(FILE* fp);												*/
	/*																				*/
	/*	This function tries to find out whether fp points to a flamelet file.		*/
	/*	If the word "Header" appears as the first word within the first 19			*/
	/*	characters of the file, then it is assumed that fp points to a flamelet		*/
	/*	file; otherwise the file is signalled as not being a flamelet file.			*/
	/*																				*/
	int result = FALSE;
	char buffer[20];
	char *ptr = buffer;

	if ( fread(buffer, sizeof(char), 19, fp) == 19 ) {
		buffer[19] = '\0';
		while ( *ptr && isspace (*ptr) )
			++ptr;
		if ( *ptr == 'H' || *ptr == 'h' ) {
			++ptr;
			if ( *ptr == 'E' || *ptr == 'e' ) {
				++ptr;
				if ( *ptr == 'A' || *ptr == 'a' ) {
					++ptr;
					if ( *ptr == 'D' || *ptr == 'd' ) {
						++ptr;
						if ( *ptr == 'E' || *ptr == 'e' ) {
							++ptr;
							if ( *ptr == 'R' || *ptr == 'r' ) {
								result = TRUE;
							}
						}
					}
				}
			}
		}
	}
	return result;
}


static char *CleanString( char *buffer )
{
	char	*bptr = buffer, *ptr = NULL;
	
	/*	remove leading spaces
	*/
	while ( *bptr && isspace( *bptr ) ) ++bptr;
	
	/*	check for comment and empty lines
	*/
	if ( *bptr == '\0' || *bptr == '#' ) return NULL;

	/*	save start of symbol
	*/
	ptr = bptr;	
	while ( *bptr && !isspace( *bptr )  ) ++bptr;

	/*	terminate symbol
	*/	
	*bptr = '\0';

	return ptr;
}


static ReportSymbolPtr NewReportSymbol( const char *name, int len )
{
	ReportSymbolPtr	ptr = NEW( ReportSymbol );
	
	ptr->name = NEW2( strlen(name) + 1, char );
	strcpy( ptr->name, name );
	ptr->data = NEW2( len, Double );
	
	return ptr;
}


static void FreeReportSymbol( ReportSymbolPtr ptr )
{
	DELETE( ptr->data );
	DELETE( ptr->name );
	DELETE( ptr );
}


static void delete_report_symbols( void *ptr, void *aux )
{
#ifdef applec
#pragma unused(aux)
#endif
	FreeReportSymbol( (ReportSymbolPtr)ptr );
}


static OptionsPtr NewOptions( void ) 
{
	OptionsPtr opts = NEW( Options );
	
#	ifdef applec
	strcpy( opts->creatorType, kDefaultCreatorType );
	opts->creator = 'QKPT';
	opts->type = 'CGTX';
#	endif

	opts->flameletFiles = NewList( "Flamelets" );
	opts->symbolList = NewList( "symbols" );
	opts->reportSymbolList = NewList( "report symbols" );

	return opts;
}


static void DeleteOptions( OptionsPtr opts )
{
	if ( !opts ) return;
	
	DELETE( opts->suffix );
	DELETE( opts->repFile );
	DELETE( opts->symFile );
	DELETE( opts->symOutFile );	
	DELETE( opts->coordOutFile );	
	DELETE( opts->coordName );	
	DeleteList( opts->flameletFiles );
	DeleteList( opts->symbolList );
	/*	first remove data from lis, since DELETE() does not work...
	*/
	ForAll( opts->reportSymbolList, delete_report_symbols, NULL, kListHead );
	RemoveList( opts->reportSymbolList );
	
	DELETE( opts );
}

static void PrintOptions( OptionsPtr opts, FILE *fp )
{
	const char	*ptr = NULL;
	
	if (opts->suffix)		fprintf(fp,"# output suffix  : \".%s\"\n", opts->suffix );
	if (opts->repFile)		fprintf(fp,"# report filename  : \"%s\"\n", opts->repFile );
	if (opts->symFile)		fprintf(fp,"# symbol filename  : \"%s\"\n", opts->symFile );

#	ifdef applec
	ptr = (const char *)( &opts->creator );
	fprintf( fp, "# creator '%c%c%c%c', ", ptr[0], ptr[1], ptr[2], ptr[3] );
	ptr = (const char *)( &opts->type );
	fprintf( fp, "type '%c%c%c%c'\n", ptr[0], ptr[1], ptr[2], ptr[3] );
#	endif
}


static void PrintHelp( const char *program, FILE *fp )
{
	if ( !fp ) fp = stdout;
	
	fprintf( fp, "# usage: "
		"%s -h | %s [options] <file1> [<file2> ...]\n\n", program, program );
	fputs( "options : \n", fp );
	fputs( "  -h                  # print help\n", fp );
	fputs( "  -o <output suffix>  # suffix for output file containing vectors\n", fp );
	fputs( "  -p                  # progress Flag\n", fp );
	fputs( "  -r <report file>    # print out max/min values of vectors to report file\n", fp );
	fputs( "  -s <symbol file>    # print out the symbols in header to symbol file\n", fp );
	fputs( "  -t                  # print timing information\n", fp );
	fputs( "  -H                  # print out header to stdout\n", fp );
	fputs( "  -M                  # strip leading \"massfraction\"/\"molefraction\".\n", fp );
	fputs( "  -T                  # print out trailer to stdout\n", fp );
	fputs( "  -S <id[,coordinate]># create hdf-file for <id> over (time,<coordinate>) from flamelet files\n", fp );
	fputs( "  -g                  # create gnuplot compatible file instead of hdf file\n", fp );
	fputs( "  -l <name>,<value>   # reports all vector information at location <value> of coordinate <name>\n", fp );
	fputs( "  -X \"string\"         # set type and creator for output file, ", fp );
	fprintf( fp, "e.g. \"MPS TEXT\", default: \"%s\".\n", kDefaultCreatorType );
	
}


static void CollectSymbols( ListPtr list, FILE *inFile )
{
	char			buffer[kBufLen+1];
	const char		*bptr = NULL;
	char			*name = NULL;

	while ( fgets( buffer, kBufLen, inFile ) ) {
		if ( !( bptr = CleanString( buffer ) ) ) continue;
		
		name = NEW2( strlen(bptr) + 1, char );
		strcpy( name, bptr );
/*		if ( opts->stripMF ) {
			StripID( name, bptr );
		}
*/
#ifdef qDebug
		fprintf( stderr, "# added \"%s\" to the symbol list\n", bptr );
#endif
		AddItem( list, name, kListTail );
		
#		ifdef applec
		RotateCursor( -32 );
#		endif
	}
}


static void HandleCommandLine( int argc, char *argv[], OptionsPtr opts )
{
	int 		c, errflg = 0;
	char		*ptr = NULL, *tempString = NULL;
	extern char	*optarg;
	extern int	optind;
	
	while ( ( c = getopt( argc, argv, "hptgHMTl:o:r:s:S:X:" ) ) != EOF ) {
		switch ( c ) {
			case 'o':
				if ( !strlen( optarg ) ) {		/*	prevent an empty suffix	*/
					Warning( "empty suffix ignored, using \"%s\" instead", kSuffix );
					break;
				}
				opts->suffix = NEW2( strlen(optarg)+1, char );
				strcpy( opts->suffix,optarg );
				break;
			case 'r':
				opts->repFile = NEW2( strlen(optarg)+1, char );
				strcpy( opts->repFile, optarg );
				break;
			case 's':
				opts->symFile = NEW2( strlen(optarg)+1, char );
				strcpy( opts->symFile, optarg );
				break;
			case 'H':
				opts->printHeader = !opts->printHeader;
				break;
			case 'p':
				opts->progress = !opts->progress;
				break;
			case 'g':
				opts->slicesAsText = !opts->slicesAsText;
				break;
			case 'T':
				opts->printTrailer = !opts->printTrailer;
				break;
			case 'h':
				opts->help = !opts->help;
				break;
			case 't':
				opts->time = !opts->time;
				break;
			case 'M':
				opts->stripMF = !opts->stripMF;
				break;
			case 'l':
				tempString = NEW2( strlen(optarg)+1, char );
				strcpy( tempString, optarg );
				ptr = strtok( tempString, "," );
				if ( !ptr ) FATAL( "error in -l argument" );
#ifdef DEBUG
				fprintf( stderr, "# coordinate = \"%s\"\n", ptr );
#endif
				opts->coordName = NEW2( strlen(ptr)+1, char );
				strcpy( opts->coordName, ptr );
				ptr = strtok( NULL, "," );
				if ( !ptr ) FATAL( "error in -l argument" );
#ifdef DEBUG
				fprintf( stderr, "# location = \"%s\"\n", ptr );
#endif
				opts->loc = atof( ptr );
				if ( tempString ) {
					DELETE( tempString );
					tempString = NULL;
				}
#ifdef DEBUG
				fprintf( stderr, "# report %s = %g\n", opts->coordName, opts->loc );
#endif
				break;
			case 'X':
				if ( strlen(optarg) >= 8 ) {
					strncpy( opts->creatorType, optarg, 8 );
					opts->creatorType[8] = '\0';	/* terminate string correctly */
				}
				break;
			default:
				++errflg;
		}
	}
	/*	Build list of input data files.												*/
	for ( ; optind < argc; optind++ ) {
		FILE	*fp = NULL;
		char	*name = NULL;
		
		if ( (fp = fopen(argv[optind], "r")) != NULL ) {
			if ( IsFlameletFile(fp) ) {
				name = NEW2( strlen(argv[optind]) + 1, char );
				strcpy( name, argv[optind] );
				AddItem( opts->flameletFiles, name, kListTail );
			} else {
				Warning("\"%s\" is not a flamelet file", argv[optind]);
			}
			fclose(fp);
		} else {
			Warning("Couldn't access flamelet file \"%s\"", argv[optind]);
		}
	}	

	if ( errflg || opts->help ) {
		PrintHelp( argv[0], stdout );
		DeleteOptions( opts );
		exit( (errflg) ? 1 : 0 );
	}
	
	/*	do some error checking
	*/
	Progress( "checking input files" );
	if ( !ItemsInList( opts->flameletFiles ) ) {
		PrintHelp( argv[0], stdout );
		DeleteOptions( opts );
		FATAL( "flamelet input file required !!" );
	}
	
	/*	if neither suffix nor report file are specified, create output file only,
		if both are specified, ignore report file.
	*/
	if ( !opts->suffix && !opts->repFile ) {
		opts->suffix = NEW2( strlen(kSuffix) + 1, char );
		strcpy( opts->suffix, kSuffix );
	}
	if ( opts->suffix && opts->repFile ) {
		fputs( "# either output file or report file option is possible, "
			   "only output file is generated.\n", stderr );
		DELETE( opts->repFile );
		opts->repFile = NULL;
	}
	
	if ( opts->progress ) PrintOptions( opts, stderr );

#	ifdef applec
	if ( strlen( opts->creatorType ) == 8 ) {
		memcpy( &opts->creator, &opts->creatorType[0], sizeof(opts->creator) );
		memcpy( &opts->type, &opts->creatorType[4], sizeof(opts->type) );
	}
	else {
		fprintf( stderr, "# illegal creator/type combination '%s', using '%s'.\n", 
			opts->creatorType, kDefaultCreatorType );
	}
#	endif
	
	if ( opts->symFile ) {
		const char	*symOutSuffix = ".out";
		FILE		*symin = NULL;
		
		if ( !( symin = fopen( opts->symFile, "r" ) ) ) {
			FATAL( "couldn't open symbol file" );
		}
		
		opts->symOutFile = NEW2( strlen(opts->symFile) + strlen(symOutSuffix) + 1, char );
		sprintf( opts->symOutFile, "%s%s", opts->symFile, symOutSuffix );
		CollectSymbols( opts->symbolList, symin );
		
		fclose( symin );
	}
	if ( opts->coordName ) {
		opts->coordOutFile = NEW2( strlen(opts->coordName) + 4, char );
		sprintf( opts->coordOutFile, "%s.kg", opts->coordName );
	}

}


void PrintArrays( ListPtr VectorList, const char *fname, OptionsPtr opts )
{
	int 			n = VectorList->fNumItems;
	Double			**ptr = NEW2( n, Double * );
	FLRSymbolPtr	sym = NULL;
	LinkPtr			temp = NULL;
	int	  			i = 0, j = 0, len = 0;
	char			sep;
	FILE			*fp = NULL;
	char			*name = NEW2( 256, char );
		
	if ( !( fp = fopen( fname,"w" ) ) ) FATAL( "error opening output file" );
	fputc( '*', fp );
	
	sep = '\n';
	for( i = 0, temp = VectorList->fHead; temp != NULL; ++i, temp = temp->fNext ) {
		sym = temp->fItem;

		ptr[i] = sym->val.v->vec;
		if ( len ) {
			if ( len != sym->val.v->len ) FATAL( "inconsistent array dimensions" );
		}
		else {
			len = sym->val.v->len;
		}
		
		strcpy( name, sym->id );
		if ( opts->stripMF ) {
			StripID( name, sym->id );
		}
		
		fprintf( fp,"%c%s", sep, name );
		if ( sym->unit && strcmp( sym->unit, "1" ) != 0 ) 
			fprintf( fp, " [%s]", sym->unit );

		sep = '\t';
	}
	
	for ( j = 0; j < len; ++j ) {
		fprintf( fp,"\n%g", ptr[0][j] );
		for ( i = 1; i < n; ++i ) {
			fprintf( fp,"\t%g", ptr[i][j] );
		}
#		ifdef applec
		RotateCursor( -16*j );
#		endif
	}
	PUT_NL( fp );
	fclose( fp );
	
	DELETE( ptr );
	DELETE( name );
	
#	ifdef applec
	fsetfileinfo( fname, opts->creator, opts->type );
#	endif
}


static void PrintHeader( FILE *fp ) 
{
	WriteFlameletFile( fp, kFLRHeader );
}


static int FindTrailer( FILE *fp )
{
	int		c;
	
#	define	nextchar( c )	if ( ( c = getc( fp ) ) == EOF ) return FALSE
	
	while ( ( c = getc( fp ) ) != EOF ) {
		if ( c == 't' || c == 'T' ) {
			nextchar( c );
			if ( c == 'r' || c == 'R' ) {
				nextchar( c );
				if ( c == 'a' || c == 'A' ) {
					nextchar( c );
					if ( c == 'i' || c == 'I' ) {
						nextchar( c );
						if ( c == 'l' || c == 'L' ) {
							nextchar( c );
							if ( c == 'e' || c == 'E' ) {
								nextchar( c );
								if ( c == 'r' || c == 'R' ) return TRUE;
							}
						} 
					}
				}
			}
		}
		if (c == '"' ) {
			while ( ( c = getc( fp ) ) != EOF ) {
				if ( c == '"' ) break;
			}
		}
	}
#	undef nextchar

	return FALSE;
}


static void PrintTrailer( FILE *in, FILE *out )
{
	if ( !FindTrailer( in ) ) {
		fputs( "# trailer not found.\n", stderr );
		return;
	} else {
		int		c = 0;
		
		if ( !out ) out = stdout;
		fputs( "\nTrailer", out );
		while ( ( c = getc( in ) ) != EOF ) putc( c, out );
		PUT_NL( out );
	}
}


static void HandleSymbol( void *ptr, void *aux )
{
#ifdef applec
#pragma unused( aux )
#endif
	const char		*id = (const char *)ptr;
	char			*name = NULL, *unit = NULL;
	VectorPtr		vec = NULL;
	int				found = 0;
	Double			val;
	
	name = NEW2( kBufLen+1, char );
	unit = NEW2( kBufLen+1, char );

	found = 0;
	if ( GetString( id, name ) ) {
		++found;
#		ifdef qDebug
		fprintf( stdout, "# found string \"%s\": \"%s\".\n", id, name );
#		endif
	}
	if ( GetNumber( id, &val, unit ) ) {
		++found;
#		ifdef qDebug
		fprintf( stdout, "# found number \"%s\": %g [%s]\n", id, val, (unit) ? unit : "" );
#		endif
	}
	if ( GetArray( id, &vec, unit ) ) {
		++found;
#		ifdef qDebug
		fprintf( stdout, "# found array \"%s\": %d elements stored at 0x%X, unit [%s]\n", 
			id, vec->len, vec, (unit) ? unit : "" );
#		endif
	}
	
	if ( !found ) {
		fprintf( stderr, "# symbol \"%s\" not found.\n", id );
	}
	
#	ifdef applec
	RotateCursor( -32 );
#	endif

	DELETE( name );
	DELETE( unit );
}


static Double maxDblValue( const Double *f, int n )
{
	Double	m = *f;
	
	++f;
	while ( --n ) {		/* pre-increment since we already did the first element */
		if ( m < *f ) m = *f;
		++f;
	}
	
	return m;
}


static int find_report_symbol( void *ptr, void *aux )
{
	return ( FLRCompare( ((ReportSymbolPtr)ptr)->name, (const char *)aux ) == 0 );
}


/*	Locate employs a binary search algorithm to return the index in the 
	sorted array xx dimensioned n, for which xx[index] <= x <= x[index+1]. 
	For x out of range, either -1 or n-1 are returned, respectively.
*/
static int locate( const Double xx[], int n, Double x )
{
	int		ascnd, upper = n, middle, lower = -1;
	
	ascnd = ( xx[n-1] > xx[0] );			/* ordering scheme (ascending/descending)	*/
	
	while ( upper - lower > 1 ) {
		middle = ( upper + lower ) >> 1;	/* compute middle point						*/
		if ( ( x >= xx[middle] ) == ascnd )	/* bisection								*/
			lower = middle;
		else
			upper = middle;
	}

#	ifdef applec
	RotateCursor( 32 );
#	endif

	return lower;
}


static Double FindLoc( Double val, VectorPtr coord, VectorPtr data )
{
	const int	n = coord->len;
	int			loc = -1;
	Double		*x = coord->vec, *y = data->vec;
	
	loc = locate( x, n, val );
	if ( loc < 0 || loc >= n-1 ) {
		FATAL( "location out of bounds" );
	}
	
	return y[loc-1] + (val - x[loc-1]) / (x[loc] - x[loc-1]) * (y[loc] - y[loc-1]);
}


static void CollectReport( OptionsPtr opts, ListPtr list, int len, int count )
{
	LinkPtr			link = NULL;
	FLRSymbolPtr	sym = NULL;
	ReportSymbolPtr	rep = NULL;
	VectorPtr		indepCoord = NULL;
	ListPtr			symList = opts->symbolList;
	char			*name = NEW2( 255, char );
	Double			val = 0.0;
	
	ForAll( symList, HandleSymbol, NULL, kListHead );
	
	for ( link = symList->fHead; link != NULL; link = link->fNext ) {
		const char	*id = (char *)link->fItem;
		
		if ( GetNumber( id, &val, NULL ) ) {
			if ( (rep = (ReportSymbolPtr)FindItem( list, (void *)id, find_report_symbol ) ) == NULL ) {
				rep = NewReportSymbol( id, len );
#ifdef qDebug
				fprintf( stderr, "# symbol \"%s\" added to report symbol list.\n", id );
#endif
				AddItem( list, rep, kListTail );
			} else {
#ifdef qDebug
				fprintf( stderr, "# symbol \"%s\" found in report symbol list.\n", id );			
#endif
			}
			rep->data[count] = val;
		}
	}
	
	if ( opts->coordName ) {
		for ( link = gFLRArrayList->fHead; link != NULL; link = link->fNext ) {
			sym = (FLRSymbolPtr)link->fItem;
			if ( strcmp( sym->id, opts->coordName ) == 0 ) {
				indepCoord = sym->val.v;
				break;
			}
		}
		if ( !indepCoord ) {
			fprintf( stderr, "# independent coordinate \"%s\" not found in flamelet file\n", opts->coordName ); FATAL( "" );
		}
	}
#ifdef DEBUG
	fprintf( stderr, "# independent coordinate \"%s\" found in flamelet file\n", opts->coordName );
#endif
	
	if ( opts->repFile || indepCoord ) {
		for ( link = gFLRArrayList->fHead; link != NULL; link = link->fNext ) {
			sym = (FLRSymbolPtr)link->fItem;
			strcpy( name, sym->id );
			if ( opts->stripMF ) {
				StripID( name, sym->id );
			}
			if ( (rep = (ReportSymbolPtr)FindItem( list, name, find_report_symbol ) ) == NULL ) {
				rep = NewReportSymbol( name, len );
#ifdef qDebug
				fprintf( stderr, "# array \"%s\" added to report symbol list.\n", name );
#endif
				AddItem( list, rep, kListTail );
			} else {
#ifdef qDebug
				fprintf( stderr, "# array \"%s\" found in report symbol list.\n", name );			
#endif
			}
			if ( indepCoord ) {
				rep->data[count] = FindLoc( opts->loc, indepCoord, sym->val.v );
			} else {
				rep->data[count] = maxDblValue( sym->val.v->vec, sym->val.v->len );
			}
#			ifdef applec
			RotateCursor( -32 );
#			endif
		}
	}
	
	DELETE( name );
}


void Process( void *ptr, void *aux )
{
	int			status = 0, ok = 0;
	FILE		*fpout = NULL, *fpin = NULL;
	Double		t;
	const char	*inFile = (const char *)ptr;
	char		*outFile = NULL;
	static int	count = 0;
	
	OptionsPtr	opts = (OptionsPtr)aux;

	t = Time();
	Progress("initializing Flamelet Reader");
	InitFLReader();
	fprintf( stderr, "\tprocessing \"%s\"...\n", inFile );
	if ( !( fpin = fopen( inFile, "r" ) ) ) FATAL("can't open input file");
	
	Progress( "reading flamelet file" );
	status = ReadFlamelet( fpin );
	if ( !status ) {
		FATAL( "error reading flamelet file" );
	}
	fclose( fpin );
	
	if ( opts->suffix ) {
		outFile = NEW2( strlen(inFile) + strlen(opts->suffix) + 1, char );
		sprintf( outFile, "%s%s", inFile, opts->suffix );
		Progress( "writing output file" );
		PrintArrays( gFLRArrayList, outFile, opts );
		DELETE( outFile );
	}
	if ( opts->printHeader ) {
		Progress( "printing header section to stdout" );
		PrintHeader( stdout );
	}
	
	if ( opts->printTrailer ) {
		if ( !( fpin = fopen( inFile, "r" ) ) ) FATAL("can't open flamelet file");
		Progress( "printing trailer section to stdout" );
		PrintTrailer( fpin, stdout );
		fclose( fpin );
	}
	
	/*	collect symbol and/or report information.
		if a report file is requested, symbols are prepended to the report file
	*/
	if ( opts->repFile || ItemsInList( opts->symbolList ) || opts->coordName ) {
		const int	n = ItemsInList(opts->flameletFiles);
		
		Progress( "processing symfile/generating report" );
		CollectReport( opts, opts->reportSymbolList, n, count );
	}
	
	
	Progress( "cleaning up FlameletReader" );
	CleanupFLReader();
	++count;
	
	t = Time() - t;
	if ( opts->time ) fprintf( stderr, "# time to process flamelet file : %.2f seconds.\n", t );
}


static void print_report_symbols_header( void *ptr, void *aux )
{
	ReportSymbolPtr	sym = (ReportSymbolPtr)ptr;
	FILE			*fp = (FILE *)aux;

	fprintf( fp, "%s\t", sym->name );
}


static int gCount = 0;

static void print_report_symbols( void *ptr, void *aux )
{
	ReportSymbolPtr	sym = (ReportSymbolPtr)ptr;
	FILE			*fp = (FILE *)aux;
	
	fprintf( fp, "%g\t", sym->data[gCount] );
}


void GenerateReport( ListPtr list, const char *fname, int len )
{
	FILE	*fp = NULL;
	int		i;
	
	if ( !(fp = fopen( fname, "w" ) ) ) FATAL( "error opening report file" );
	
	fputs( "*\n", fp );
	ForAll( list, print_report_symbols_header, fp, kListHead );
	PUT_NL( fp );
	for ( i = 0; i < len; ++i ) {
		gCount = i;
		ForAll( list, print_report_symbols, fp, kListHead );
		PUT_NL( fp );
	}
	
	fclose( fp );
	
#ifdef applec
	fsetfileinfo( fname, 'QKPT', 'CGTX' );
#endif
}

int main( int argc, char *argv[] )
{
	OptionsPtr	opts = NULL;
	int			n = 0;

#	ifdef applec
	InitCursorCtl( NULL );
	RotateCursor( 32 );
#	endif

	opts = NewOptions();
	
	HandleCommandLine( argc, argv, opts);
	
	n = ItemsInList(opts->flameletFiles);
	ForAll( opts->flameletFiles, Process, opts, kListHead );
	if ( ItemsInList(opts->reportSymbolList) ) {
		GenerateReport( opts->reportSymbolList, 
			(opts->repFile) ? opts->repFile : 
			((opts->symOutFile) ? opts->symOutFile : opts->coordOutFile), n );
	}

	DeleteOptions( opts );

	return 0;
}

