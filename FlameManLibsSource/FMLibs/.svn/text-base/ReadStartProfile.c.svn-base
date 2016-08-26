#define FOOL_SOFTBENCH( x ) 

/*
	ReadStartProfile.c
*/

#include "ListTool.h"

#undef TIMETEST

#undef DEBUG

/* Globals */
char		*gStream = NULL;	/* pointer to header(-tokens) */
char		*gLabels = NULL;	/* pointer to vector identifiers */
int			gNumOfParas = 0,	/* number of assignment statements */
			gRows = 0,			/* the length of a vector */
			gVectors = 0;		/* total number of vectors */
doublePtr	gData = NULL;		/* pointer to floating point data */
parameter	*gParas = NULL;			/* pointer to header contents */

void CleanReadStartProfile(void)
{
	if (gStream)			free(gStream);
	if (gLabels)			free(gLabels);
	if (gParas)				free(gParas);
	if (gData)				free(gData);
}


void printHeader(char *fname)
{
	int		i;
	FILE	*fp = stderr;
	
	fprintf( fp, "\nFile: %s\n", fname);
	fprintf( fp, "%-28s%-14s%s\n", "Identifier", "Context", "Value" );
	for ( i = 0; i < 90; ++i ) putc( '-', fp ); putc( '\n', fp );
	
	for ( i = 0; i < gNumOfParas; ++i ) {
		fprintf( fp, "%-28s%-14s", gParas[i].identifier, gParas[i].context );
		if ( gParas[i].tag == kTextType )
			fprintf( fp, "%s\n", gParas[i].what.string );
		else
			fprintf( fp, "%-12.5g [%s]\n",
				gParas[i].what.quantity.value, gParas[i].what.quantity.unit );
	}
	putc( '\n', fp );
}

double GetVariable( const char *varName )
{
	int i;
	for ( i = 0; i < gNumOfParas; ++i ) {
		if ( strcmp( gParas[i].identifier, varName ) == 0 ) {
			if ( gParas[i].tag == kPhysicalQuantity ) {
				return gParas[i].what.quantity.value;
			}
			else {
				fprintf( stderr, "the variable %s has no value\n", gParas[i].identifier );
				exit(2);
			}
		}
	}
	fprintf( stderr, "can't find variable %s\n", varName );
	exit(2);
	
	return 0.0;
}

parameter *GetParameter( const char *identifier )
{
	int			i;
	parameter	*ptr = NULL;
	
	for ( i = 0; i < gNumOfParas; ++i ) {
		if ( streq(identifier, gParas[i].identifier) ) {
			ptr = &gParas[i];
			break;
		}
	}
	return ptr;
}


void OutOfRange( const char *id, double value )
{
	fprintf( stderr, "# %s out of range: %g\n", id, value );
}


/*
int GetHeaderInfo( DataSetPtr ds )
{
	const char	*not_found = "Couldn't find a required parameter";
	int			i, len;
	parameter	*theParameter = NULL;
	
	if ( !(theParameter = GetParameter( "burningvelocity" )) )
		FatalError( not_found );
	ds->burningVelocity = theParameter->what.quantity.value;
	if ( ds->burningVelocity < 5.0 ) return -1;
	
	if ( !(theParameter = GetParameter( "fuel" )) )
		FatalError( not_found );
	ds->fuel = theParameter->what.string;
	
	if ( ds->fuel[0] == '"' ) {	
		len = strlen(ds->fuel);
		for ( i = 0; i < len; ++i ) ds->fuel[i] = ds->fuel[i+1];
		ds->fuel[len-2] = '\0';
	}
	for ( i = CH4; i <= HYDROGEN; ++i ) {
		if ( strcmp( ds->fuel, kFuelName[i] ) == NULL ) break;
	}
	if ( i > HYDROGEN ) FatalError( "Cannot handle this fuel" );
	ds->fuelID = i;

	if ( !(theParameter = GetParameter( "pressure" )) )
		FatalError( not_found );
	ds->pressure = theParameter->what.quantity.value;

	if ( !(theParameter = GetParameter( "fuel-air-equivalence-ratio" )) )
		FatalError( not_found );
	ds->phi = theParameter->what.quantity.value;
	if ( ds->phi > 1.05 ) return -1;

	if ( !(theParameter = GetParameter( "temperature" )) )
		FatalError( not_found );
	ds->Tu = theParameter->what.quantity.value;

	if ( !(theParameter = GetParameter( "gridpoints" )) )
		FatalError( not_found );
	ds->gridPoints = (int)theParameter->what.quantity.value;

	ds->gridPoints = gRows;
	
	CheckDataSet( ds );
	
	return 0;
}
*/



int mainz( int argc, char *argv[] )
{	
	char	*usage = "# Usage: %s -i infile -o outfile [-s] [-t] [-p]\n";
	FILE	*infp;
	char	*infile = NULL,			/* input file name (mandatory) */
			*outfile = NULL;		/* optional output file */
	int		inFlag = FALSE,			/* signals whether -i argument was succesfully scanned */
			outFlag = FALSE,		/* TRUE if using an output file */
			timingFlag = FALSE,		/* print timing info to stderr */
			scalarFlag = FALSE;		/* generate only scalar info */
 	int  	parms, time;
	StartProfile	sp;
	
#	if defined (applec) || defined (powerc)
	time = TickCount();
	InitCursorCtl(nil);
	RotateCursor(0);
#	endif
	
	for ( parms = 1; parms < argc; ++parms ) {
#		ifdef DEBUG
		fprintf(stderr,"argv[%d] = %s\n", parms, argv[parms]);
#		endif
		if (*(argv[parms]) != '-') continue;
		else if (tolower(*(argv[parms]+1)) == 'i') {
			infile = argv[parms+1];
			inFlag = TRUE;
		}
		else {
			fprintf(stderr, "### %s - \"%s\"is not an option.\n", argv[0], argv[parms]);
			fprintf(stderr, usage, argv[0]);
			return ( 1 );
		}
	}
		
	if ( !inFlag ) {
		fprintf(stderr, usage, argv[0]);
		FatalError("Input file is missing");
	}
	if (!(infp = fopen(infile, "r"))) {
		fprintf(stderr, "### Error: Couldn't open input data file \"%s\".\n", infile);
		exit(2);
	}	
	
	
	ReadStartProfiles( &sp, infp );
	
	return ( 0 );
}

StartProfilePtr ReadStartProfiles( StartProfilePtr sp, FILE *infp )
{	
	char	outfile[] = "StartProfiles.out";	/* debugging output file */
#ifdef DEBUG 
	FILE	*outFp;
#endif
#	if defined (applec) || defined (powerc)
 	int  	time;
		
	time = TickCount();
	InitCursorCtl(nil);
	RotateCursor(0);
#	endif
	

	if (!(gParas = (parameter *)calloc(kMaxParameters, sizeof(parameter))))
		FatalError("Out of memory");
	if (!(gLabels = (char *)calloc(kMaxLabelChars, sizeof(char))))
		FatalError("Memory allocation for labels failed");
	if (!(gStream = (char *)malloc(kMaxHeader * sizeof(char))))
		FatalError("Out of memory");

	parseHeader(gStream, infp);				/* is always done */

#ifdef DEBUG 
	printHeader( "" );
#endif

	parseBody( infp );							/* is always done */

	sp->labels = gLabels;
	sp->data = gData;
	sp->gridPoints = gRows;
	sp->variables = gVectors;
/*	fprintf( stderr, "********param '%s' = %g\n", GetParameter( "pressure" )->identifier, GetParameter( "pressure" )->what.quantity.value );*/
	
#ifdef DEBUG 
	fprintf(stderr, "#\t%d vectors of length %d\n", gVectors, gRows);
#endif
	
	
/*	ds->TStar = New1DArray( ds->gridPoints );
	ds->x_over_lF = New1DArray( ds->gridPoints );
	
	//	 set up x, xStar, T 
	ds->x = &gData[0];
	ds->xStar = &gData[ds->gridPoints];
	ds->T = &gData[2*ds->gridPoints];
*/

#ifdef DEBUG 
	if (!(outFp = fopen(outfile, "w"))) {
		fprintf(stderr, "### Error: Couldn't open output file \"%s\".\n", outfile);
		exit(2);
	}
	saveData( outFp );
	fclose( outFp );
#endif
	/*CleanUp();*/
	
#if defined (applec) || defined (powerc)
#ifdef TIMETEST 
	time = TickCount() - time;
	fprintf(stderr, "#\t%g Seconds\n", (double)time / 60.0);
#endif
#endif
	
	return ( sp );
}

int GetStartProfileGridPoints( FILE *infp )
{	
	char	outfile[] = "StartProfiles.out";	/* debugging output file */
 	int  	time;
	int		gridPoints;
		
#	if defined (applec) || defined (powerc)
	time = TickCount();
	InitCursorCtl(nil);
	RotateCursor(0);
#	endif
	

	if (!(gParas = (parameter *)calloc(kMaxParameters, sizeof(parameter))))
		FatalError("Out of memory");
	if (!(gLabels = (char *)calloc(kMaxLabelChars, sizeof(char))))
		FatalError("Memory allocation for labels failed");
	if (!(gStream = (char *)malloc(kMaxHeader * sizeof(char))))
		FatalError("Out of memory");

	parseHeader(gStream, infp);				/* is always done */
	gridPoints = GetVariable( "gridpoints" );
#ifdef DEBUG 
	printHeader( "" );
#endif
	CleanReadStartProfile();
	
#ifdef TIMETEST 
	time = TickCount() - time;
	fprintf(stderr, "#\t%g Seconds\n", (double)time / 60.0);
#endif
	
	return gridPoints;
}
