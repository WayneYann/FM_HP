/*
	CreateBinFile.c: program to create a binary data bank of thermodynamic
		and transport properties.
*/

#ifdef applec
#include <CursorCtl.h>
pascal unsigned long TickCount(void) = 0xA975; /* from Events.h */
#endif

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ArrayManager.h"

#include "ThermoProperties.h"

#undef DEBUG

/*#ifdef applec*/
#define	DUMPSTRUCTS
/*#elif HP*/
/*#define	DUMPSTRUCTS*/
/*#else*/
/*#undef	DUMPSTRUCTS*/
/*#endif*/

#ifndef FALSE
#define FALSE	0
#define TRUE	1
#endif

#include "ThermoProperties.h"

#define kMaxLen		512			/* max. length of an input line/output record */


static void FullNameForDataFiles( char *fullName, const char *fileName );

/*
	Globals
*/
char *usage = "# Usage: %s -i infile -o outfile [-m molfile] [-h] [-p] [-d]\n";
int	 gProgress = FALSE;				/* print some progress info to stderr	*/


void FatalError( char *errorString )
{
	fprintf( stderr, "### %s\n", errorString );
	exit( 2 );
}


void DoHelp( void )
{
	fprintf( stderr, "# Format of the input file:\n" );
	fprintf( stderr, "#\tName  a_i  b_i  M eps/k  sigma\n" );
}


void ReadWrite( FILE *inFile, FILE *outFile, FILE *molFile )
{
	int		line_count = 0;			/* counts the number of lines			*/
	int		items_read;				/* counts the number of items parsed	*/
	char	buffer[kMaxLen];		/* holds one line/recort at a time		*/
	char	molName[kMaxLen];		/* holds the molspec name				*/
	char	*format = "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf";
	char	*formtherm1 = "%s";
/*	char	*formtherm2 = "%15lf%15lf%15lf%15lf%15lf";*/
/*	char	*formtherm4 = "%15lf%15lf%15lf%15lf";*/
	char	*formtherm2 = "%lf%lf%lf%lf%lf";
	char	*formtherm4 = "%lf%lf%lf%lf";
	size_t	size;					/* holds the number of bytes written to disk */
	SpeciesRecord sr;				/* holds the info about a species		*/
	Double	eps_over_k;
	Double	dummyDouble;
	int		dummyInt;
/**go** atomic weights as defines including F and BR added */    
	const Double	AWT_H = 1.00797;
	const Double	AWT_C = 12.01115;
	const Double	AWT_N = 14.00670;
	const Double	AWT_O = 15.99940; 
	const Double	AWT_F = 18.99840;
	const Double	AWT_BR = 79.90900;  
	const Double	AWT_AR = 39.94400;
	const Double    AWT_HE = 4.00;
#	ifndef DUMPSTRUCTS
	char	*ptr;					/* current insertion position of buffer	*/
	static int first = TRUE;
#	endif

	while ( fgets( buffer, kMaxLen, inFile ) ) {

		++line_count;
#		ifdef applec
		RotateCursor( 16 * line_count );
#		endif
	
		if ( (*buffer == '#') || (*buffer == '\n') )
			continue;				/* skip comment and empty lines			*/
		
		/*
			Parse buffer
		*/
		memset( sr.name, 0, kNameLen );
		if ( molFile ) {
			items_read = sscanf( buffer, formtherm1, sr.name );

			if ( strcmp( "THERMO", sr.name ) == 0 ) {
				fgets( buffer, kMaxLen, inFile );
				/*fprintf( stderr, "THERMO found next line is: '%s'\n", buffer );*/
				continue;
			}
			else if ( strcmp( "END", sr.name ) == 0 ) {
				/*fprintf( stderr, "comment found: '%s'\n", buffer );*/
				continue;
			}
			else if ( *buffer == '!' ) {
				/*fprintf( stderr, "comment found: '%s'\n", buffer );*/
				continue;
			}

			fgets( buffer, kMaxLen, inFile );
			items_read += sscanf( buffer, formtherm2
					, &sr.hot[0], &sr.hot[1], &sr.hot[2], &sr.hot[3], &sr.hot[4] );
			fgets( buffer, kMaxLen, inFile );
			items_read += sscanf( buffer, formtherm2
					, &sr.hot[5], &sr.hot[6], &sr.cold[0], &sr.cold[1], &sr.cold[2] );
			fgets( buffer, kMaxLen, inFile );
			items_read += sscanf( buffer, formtherm4
					, &sr.cold[3], &sr.cold[4], &sr.cold[5], &sr.cold[6] );

/*			fprintf( stderr, "species %s found in infile\n", sr.name );*/
/*			fprintf( stderr, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n"*/
/*				, sr.hot[0], sr.hot[1], sr.hot[2], sr.hot[3], sr.hot[4], sr.hot[5], sr.hot[6] );*/
/*			fprintf( stderr, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n"*/
/*				, sr.cold[0], sr.cold[1], sr.cold[2], sr.cold[3], sr.cold[4], sr.cold[5], sr.cold[6] );*/
			sr.M = -1.0;
			eps_over_k = -1.0;
			sr.sigma = -1.0;
			rewind(molFile);
			while ( fgets( buffer, kMaxLen, molFile ) ) {
				if ( *buffer == '!' ) {
					continue;
				}
				sscanf( buffer, "%s", molName );
				if ( strcmp( sr.name, molName ) == 0 ) {
					sscanf( buffer, "%s%d%lf%lf%lf%lf%lf"
						, molName, &dummyInt, &eps_over_k, &sr.sigma
						, &dummyDouble, &dummyDouble, &dummyDouble );
/*					fprintf( stderr, "species %s found in moldat\n", molName );*/
					break;
				}
			}
/*			fprintf( stderr, "M = %g\tepsoverk = %g\tsigma = %g\n"*/
/*				, sr.M, eps_over_k, sr.sigma );*/
		}
		else {
			items_read = sscanf( buffer, format, sr.name,
				&sr.hot[0], &sr.hot[1], &sr.hot[2], &sr.hot[3], &sr.hot[4],
				&sr.hot[5], &sr.hot[6], &sr.cold[0], &sr.cold[1], &sr.cold[2],
				&sr.cold[3], &sr.cold[4], &sr.cold[5], &sr.cold[6],
				&sr.M, &eps_over_k, &sr.sigma );
			if ( items_read != 18 ) {
				fprintf( stderr, "Line %4d   #  Not a valid record\n", line_count );
				continue;				/* try next line						*/
			}
		}
		sr.k_over_eps = 1.0 / eps_over_k;	/* save the inverse value		*/
		
		
		/*
			Compute mu_coeff and the number of C, H, N, or O atoms...
		*/

/**go** 'ComputeMuCoeff' moved, otherwise MuCoeff is computed with M=-1.0 if M is not provided
    ComputeMuCoeff( &sr );
**go**/

/*		ComputeComposition( &sr );*/
		/*
			Print a message, if the record doesn't contain a valid mol mass,
			compute the molar mass
		*/
/*		if ( sr.M < 0.0 ) {*/
/*			fprintf( stdout, "# \"%s\":\tmolar mass unknown\n", sr.name );*/
			/*continue;	*/			/* don't save it in the binary file		*/

/**go** new calc of M with defines and atoms F and BR added */
/*			sr.M = (Double) sr.composition_elems.C  * AWT_C  +*/
/*				(Double) sr.composition_elems.H  * AWT_H  +*/
/*				(Double) sr.composition_elems.N  * AWT_N  +*/
/*				(Double) sr.composition_elems.O  * AWT_O  +*/
/*				(Double) sr.composition_elems.F  * AWT_F  +*/
/*				(Double) sr.composition_elems.BR * AWT_BR +*/
/*				(Double) sr.composition_elems.AR * AWT_AR;*/
/**/
/*      fprintf( stdout, "  -> calculated molar mass = %lf\n", sr.M );*/

/**go** old calculation of M **/
/*			sr.M = (Double) sr.composition_elems.C * 12. +*/
/*				(Double) sr.composition_elems.H * 1.008 +*/
/*				(Double) sr.composition_elems.N * 14.01 +*/
/*				(Double) sr.composition_elems.O * 16.000;*/
	
/**go**/
/*fprintf( stderr, "%-10s %16.8E\n", sr.name, sr.M );*/
/*		}*/

/**go** new position of 'ComputeMuCoeff' */		
/*		ComputeMuCoeff( &sr );*/
		/* fill up unused fields */
		sr.mu_coeff = -1.0;		
		sr.D_coeff = sr.omega_coeff = NULL;
		sr.id = -1;
		
		
		/*
			...and dump everything to disc
		*/
        UpperString(sr.name);
#		ifdef DUMPSTRUCTS

		if ( ( size = fwrite( &sr, sizeof(sr), 1, outFile ) ) != 1 ) {
			fprintf( stderr, "# fwrite() returned %d instead of expected 1.\n",
                (int)size );
		}

#		else		/* user buffer */

		ptr = buffer;
		/* fill up the buffer */
		memcpy( ptr, sr.name, kNameLen );
		ptr += kNameLen;
		memcpy( ptr, &sr.M, sizeof(Double) );
		ptr += sizeof(Double);
		memcpy( ptr, &sr.composition_elems, sizeof(Composition) );
		ptr += sizeof(Composition);
		memcpy( ptr, sr.hot, 7 * sizeof(Double) );
		ptr += 7 * sizeof(Double);
		memcpy( ptr, sr.cold, 7 * sizeof(Double) );
		ptr += 7 * sizeof(Double);
		memcpy( ptr, &sr.k_over_eps, sizeof(Double) );
		ptr += sizeof(Double);
		memcpy( ptr, &sr.sigma, sizeof(Double) );
		ptr += sizeof(Double);
		memcpy( ptr, &sr.mu_coeff, sizeof(Double) );
		ptr += sizeof(Double);
		memcpy( ptr, &sr.D_coeff, sizeof(Double *) );
		ptr += sizeof(Double *);
		memcpy( ptr, &sr.omega_coeff, sizeof(Double *) );
		ptr += sizeof(Double *);
		memcpy( ptr, &sr.id, sizeof(int) );
		ptr += sizeof(int);
		
		/* and write it to disk */
		size = fwrite( buffer, 1, ptr - buffer, outFile );
		if ( size != ptr - buffer ) {
			fprintf( stderr, "# fwrite() returned %d instead of expected %d.\n",
				size, ptr - buffer );
		}
		
		if ( first ) {
			first = FALSE;
			/*
				Try to find out whether structs can be dumped to disk.
			*/
			if ( memcmp( buffer, &sr, ptr - buffer ) == 0 )
				fprintf( stderr, "\n#   *** Structs may be dumped ***\n\n" );
			else
				fprintf( stderr, "\n#   !!! Do NOT dump structs !!!\n\n" );
		}

#		endif	/* DUMPSTRUCTS */

	}	/* while */
	
}

void UpperString( char *string )
{
  while ( *string ) {
	*string = toupper( *string );
	++string;
  }
}

int main( int argc, char *argv[] )
{
	char	infile[256],			/* input file name (mandatory)				*/
			molfile[256],			/* input file name for moldat     			*/
			outfile[256];			/* output file name (mandatory)				*/
	int		inFlag = FALSE,			/* TRUE, if input file seen					*/
			molFlag = FALSE,		/* TRUE, if moldat file seen				*/
			outFlag = FALSE,		/* TRUE, if output file seen				*/
			timingFlag = FALSE;		/* print timing info to stderr				*/
 	int  	parms;					/* counts the number of parameters			*/
	FILE	*inFp = NULL,			/* pointer to input file					*/
			*molFp = NULL,			/* pointer to moldat file					*/
			*outFp = NULL;			/* pointer to output file					*/
	int		outputToMyData = FALSE;	/* TRUE, if output should be written to myData 	*/
	char	fullName[256];			/* full name of outputFile ( including path ) 	*/
#	ifdef applec
	unsigned long	time;			/* counts the number of elapsed ticks	*/
#	endif	/* applec */

#	ifdef applec
	time = TickCount();
	InitCursorCtl( NULL );
	RotateCursor( 0 );
#	endif	/* applec */

	/*
		Parse command line arguments
	*/
	for (parms = 1; parms < argc; parms++) {
#		ifdef DEBUG
		fprintf(stderr,"argv[%d] = %s\n", parms, argv[parms]);
#		endif
		if (*(argv[parms]) != '-') continue;
		else if (tolower(*(argv[parms]+1)) == 'i') {
			if(sscanf(argv[parms+1], "%s", infile) != 1) {
				fprintf(stderr, usage, argv[0]);
				return ( 1 );
			}
			inFlag = TRUE;
		}
		else if (tolower(*(argv[parms]+1)) == 'm') {
			if(sscanf(argv[parms+1], "%s", molfile) != 1) {
				fprintf(stderr, usage, argv[0]);
				return ( 1 );
			}
			molFlag = TRUE;
		}
		else if (tolower(*(argv[parms]+1)) == 'o') {
			if(sscanf(argv[parms+1], "%s", outfile) != 1) {
				fprintf(stderr, usage, argv[0]);
				return ( 1 );
			}
			outFlag = TRUE;
		}
		else if (*(argv[parms]+1) == 'h') {
			fprintf(stderr, usage, argv[0]);
			DoHelp();
		}
		else if (*(argv[parms]+1) == 'd') {
			outputToMyData = TRUE;
		}
		else if (tolower(*(argv[parms]+1)) == 'p') 
			gProgress = !gProgress;
		else {
			fprintf(stderr, "### %s - \"%s\"is not an option.\n", argv[0], argv[parms]);
			fprintf(stderr, usage, argv[0]);
			return ( 1 );
		}
	}
		
	if (!inFlag) {
		fprintf(stderr, usage, argv[0]);
		FatalError("Input file is missing");
	}
	if (!outFlag) {
		fprintf(stderr, usage, argv[0]);
		FatalError("Output file is missing");
	}
	if (!(inFp = fopen(infile, "r"))) {
		fprintf(stderr, "### Couldn't open input data file \"%s\".\n", infile);
		exit(2);
	}
	if (molFlag) {
		if (!(molFp = fopen(molfile, "r"))) {
			fprintf(stderr, "### Couldn't open input data file \"%s\".\n", molfile);
			exit(2);
		}
	}
	
	if ( outputToMyData ) {
		FullNameForDataFiles( fullName, outfile );
	}
	else {
		strcpy( fullName, outfile );
	}
	if (!(outFp = fopen(fullName, "wb"))) {
		fprintf(stderr, "### Couldn't open output file \"%s\".\n", outfile);
		exit(2);
	}

	if ( gProgress ) {
		fprintf(stderr, "#\tInput file: \"%s\"\n", infile);
		fprintf(stderr, "#\tOutput file: \"%s\"\n", fullName);
	}
	
	ReadWrite( inFp, outFp, molFp );
	
	/*
		Clean up
	*/
	fclose( inFp );
	fclose( outFp );

	return 0;
}

static void FullNameForDataFiles( char *fullName, const char *fileName )
{
	int		i = 0;
	char	*pathPtr = NULL;
	
	fullName[0] = '\0';
	pathPtr = getenv( "myData" );

	if ( pathPtr ) {
#		ifdef applec
		const char	kSeparator = ':';
#		else
		const char	kSeparator = '/';
#		endif

		while ( *pathPtr ) fullName[i++] = *pathPtr++;
		if ( fullName[i-1] != kSeparator ) fullName[i++] = kSeparator;
		fullName[i] = '\0';
	}
	strcat ( fullName, fileName );
}
