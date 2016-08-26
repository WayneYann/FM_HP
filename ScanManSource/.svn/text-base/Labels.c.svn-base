/*
	Labels.c: 11/11/93
*/

#include "ScanMan.h"

static void EnumerateSpecies( FILE *fp );
static void EnumerateReactions( FILE *fp );
static void EnumerateThirdBody( FILE *fp );
static void EnumerateArrCoeffsF77( FILE *fp );
static void EmitSpecies( void *ptr, void *aux );
static void EmitReactions( void *ptr, void *aux );
static void EmitThirdBody( void *ptr, void *aux );
static void EmitArrCoeffF77A( void *ptr, void *aux );
static void EmitArrCoeffF77N( void *ptr, void *aux );
static void EmitArrCoeffF77E( void *ptr, void *aux );
static void ListSpeciesF77( FILE *fp );
static void ListReactionsF77( FILE *fp );
static void ListThirdBodyF77( FILE *fp );
static void EnumerateSpeciesF77( FILE *fp );
static void EnumerateReactionsF77( FILE *fp );
static void EnumerateThirdBodyF77( FILE *fp );
static void EmitSpeciesF77( void *ptr, void *aux );
static void EmitReactionsF77( void *ptr, void *aux );
static void EmitThirdBodyF77( void *ptr, void *aux );
static void EmitSpeciesNumF77( void *ptr, void *aux );
static void EmitReactionsNumF77( void *ptr, void *aux );
static void EmitThirdBodyNumF77( void *ptr, void *aux );
static void Emit( FILE *fp, const char *s );
static void EnumerateArrCoeffsF90( FILE *fp );
static void EmitArrCoeffF90A( void *ptr, void *aux );
static void EmitArrCoeffF90N( void *ptr, void *aux );
static void EmitArrCoeffF90E( void *ptr, void *aux );
static void ListSpeciesF90( FILE *fp );
static void ListReactionsF90( FILE *fp );
static void ListThirdBodyF90( FILE *fp );
static void EnumerateSpeciesF90( FILE *fp );
static void EnumerateReactionsF90( FILE *fp );
static void EnumerateThirdBodyF90( FILE *fp );
static void EmitSpeciesF90( void *ptr, void *aux );
static void EmitReactionsF90( void *ptr, void *aux );
static void EmitThirdBodyF90( void *ptr, void *aux );
static void EmitSpeciesNumF90( void *ptr, void *aux );
static void EmitReactionsNumF90( void *ptr, void *aux );
static void EmitThirdBodyNumF90( void *ptr, void *aux );


static int				gSpeciesNo;
static int				gReactionNo;
static int				gThirdBodyNo;
static const int	gNOfElementsLine = 2;
static const int	gNOfElementsBlock = 30;
const char		*gPrefixSteadyState = "s";
const char		*gPrefixComputed = "s";
const char		*gPrefixReactions = "r";
const char		*gPrefixThirdBody = "m";
const char			*gPrefixSteadyStateF77 = "S";
const char			*gPrefixComputedF77 = "S";
const char			*gPrefixReactionsF77 = "R";
const char			*gPrefixThirdBodyF77 = "M";


void SaveLabels( void )
{
	FILE	*fp;
	char	fname[100];
	
	sprintf( fname, "%s.h", gOptions->base );
	if ( !( fp = fopen( fname, "w") ) ) { 
		fprintf( stderr, "# SaveLabels: unable to open file %s\n", fname );
		exit(2);
	}
	fprintf( fp, "#ifndef __MAGICFILE__\n" );
	fprintf( fp, "#define __MAGICFILE__\n" );
	fprintf( fp, "#define MECHANISM \"%s\"\n", gOptions->base );
	EnumerateSpecies( fp );
	EnumerateReactions( fp );
	EnumerateThirdBody( fp );
	
	fprintf( fp, "#endif\n" );
	fclose( fp );
}


static void EnumerateSpecies( FILE *fp )
{
	gSpeciesNo	= 0;
	Emit( fp, "\ntypedef enum SpeciesLabel {\n\n" );
	Emit( fp, "\t/* Computed species s.. */\n" );
	Emit( fp, "\t/* Steady-state species ss.. */\n" );
	Iterate( gSpeciesList, EmitSpecies, fp );
	Emit( fp, "\tsEnd\n" );
	Emit( fp, "} SpeciesLabel;\n\n" );
}


static void EnumerateReactions( FILE *fp )
{
	gReactionNo	= 0;
	Emit( fp, "\ntypedef enum ReactionLabel {\n" );
	Emit( fp, "\t/* Reactions */\n" );
	Iterate( gReactionList, EmitReactions, fp );
	Emit( fp, "\t/* PAHReactions */\n" );
	Iterate( gPAHReactionList, EmitReactions, fp );
	Emit( fp, "\t/* SootReactions */\n" );
	Iterate( gSootReactionList, EmitReactions, fp );
	Emit( fp, "\trEnd\n" );
	Emit( fp, "} ReactionLabel;\n\n" );
}


static void EnumerateThirdBody( FILE *fp )
{
	gThirdBodyNo = 0;
	Emit( fp, "\ntypedef enum TirdBodyLabel {\n\n" );
	Iterate( gUsedThirdBodyList, EmitThirdBody, fp );
	Emit( fp, "\tmEnd\n" );
	Emit( fp, "} TirdBodyLabel;\n\n" );
}


char *CSymbolLeadNum( const char *src )
{
	static char	dest[80];
	char 			*ptr = dest;
	
	switch( *src ) {
		case ',':
		case '.':
			*ptr++ = 'D';
		case ';':
		case ':':
			/*	remove these characters	*/
			break;
		case '1':
			*ptr++ = 'P';
			break;
		case '2':
			*ptr++ = 'S';
			break;
		case '3':
			*ptr++ = 'T';
			break;
		case '4':
			*ptr++ = 'Q';
			break;
		case '5':
			*ptr++ = 'P';
			break;
		case '-':
			*ptr++ = 'X';
			break;
		case '*':
			*ptr++ = 'Y';
			break;
		default:
			*ptr++ = *src;
	}
	while ( *++src ) {
		switch( *src ) {
			case ',':
			case '.':
				*ptr++ = 'D';
			case ';':
			case ':':
				/*	remove these characters	*/
				break;
			case '-':
				*ptr++ = 'X';
				break;
			case '*':
				*ptr++ = 'Y';
				break;
			default:
				*ptr++ = *src;
		}
	}
	*ptr = '\0';

	return dest;
}

char *CSymbol( const char *src )
{
	static char	dest[80];
	char 			*ptr = dest;
	
	--src;
	while ( *++src ) {
		switch( *src ) {
			case ',':
			case '.':
				*ptr++ = 'P';
			case ';':
			case ':':
				/*	remove these characters	*/
				break;
			case '-':
				*ptr++ = 'X';
				break;
			case '*':
				*ptr++ = 'Y';
				break;
			default:
				*ptr++ = *src;
		}
	}
	*ptr = '\0';

	return dest;
}


static void EmitSpecies( void *ptr, void *aux )
{
	SpeciesPtr	sPtr = (SpeciesPtr)ptr;
	FILE		*fp = (FILE *)aux;
	
	if ( sPtr->isSteadyState == FALSE ) {
		fprintf( fp, "\t%s%s = %d,\n", gPrefixComputed, CSymbol(sPtr->name), gSpeciesNo++ );
	}
	else if ( sPtr->isSteadyState == TRUE ) {
		fprintf( fp, "\t%s%s = %d,\n", gPrefixSteadyState, CSymbol(sPtr->name), gSpeciesNo++ );
	}
	else {
		fprintf( fp, "# %s: unknown steady state flag \n", CSymbol(sPtr->name) );
	}
}


static void EmitReactions( void *ptr, void *aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	
	fprintf( fp, "\t%s%s = %d,\n", gPrefixReactions, CSymbol(rPtr->label), gReactionNo++ );
}


static void EmitThirdBody( void *ptr, void *aux )
{
	String			rPtr = (String)ptr;
	FILE			*fp = (FILE *)aux;
	
	fprintf( fp, "\t%s%s = %d,\n", gPrefixThirdBody, CSymbol(rPtr), gThirdBodyNo++ );
}


static void Emit( FILE *fp, const char *s )
{
#ifdef applec
	pascal void SpinCursor(short increment);

	SpinCursor( 1 );
#endif

	if (fp && s) {
		fputs(s, fp);
	}
	if ( ferror(fp) ) {
		perror( "C code generator" );
		exit( 2 );
	}
}

void SaveLabelsF77( void )
{
	FILE	*fp;
	char	fname[100];
	
	sprintf( fname, "%sF.h", gOptions->base );
	if ( !( fp = fopen( fname, "w") ) ) { 
		fprintf( stderr, "# SaveLabelsF77: unable to open file %s\n", fname );
		exit(2);
	}
	
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C ======= %s =======\n", fname );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );

	ListSpeciesF77( fp );
	ListReactionsF77( fp );
	ListThirdBodyF77( fp );
	EnumerateSpeciesF77( fp );
	EnumerateReactionsF77( fp );
	EnumerateThirdBodyF77( fp );
	fprintf( fp, "\n\n      DOUBLE PRECISION A(REND),N(REND),E(REND)\n" );
	fprintf( fp, "\n\n      INTEGER IH\n" );
	EnumerateArrCoeffsF77( fp );
	fprintf( fp, "\n" );

	fclose( fp );
}

static void ListSpeciesF77( FILE *fp )
{
	gSpeciesNo	= 0;
	Iterate( gSpeciesList, EmitSpeciesF77, fp );
	Emit( fp, "\n     &    , SEND" );
}


static void ListReactionsF77( FILE *fp )
{
	gReactionNo	= 0;
	Iterate( gReactionList, EmitReactionsF77, fp );
/*	Iterate( gPAHReactionList, EmitReactionsF77, fp );*/
/*	Iterate( gSootReactionList, EmitReactionsF77, fp );*/
	Emit( fp, "\n     &    , REND" );
}


static void ListThirdBodyF77( FILE *fp )
{
	gThirdBodyNo = 0;
	Iterate( gUsedThirdBodyList, EmitThirdBodyF77, fp );
	Emit( fp, "\n     &    , MEND" );
}

static void EnumerateSpeciesF77( FILE *fp )
{
	gSpeciesNo	= 0;
	Iterate( gSpeciesList, EmitSpeciesNumF77, fp );
	fprintf( fp, "\n     &    , SEND = %d )", 1 + gSpeciesNo );
}


static void EnumerateArrCoeffsF77( FILE *fp )
{
	gReactionNo	= 0;
	Iterate( gReactionList, EmitArrCoeffF77A, fp );
	fprintf( fp, "\n     &    /" );
	gReactionNo	= 0;
	Iterate( gReactionList, EmitArrCoeffF77N, fp );
	fprintf( fp, "\n     &    /" );
	gReactionNo	= 0;
	Iterate( gReactionList, EmitArrCoeffF77E, fp );
	fprintf( fp, "\n     &    /" );
}


static void EnumerateReactionsF77( FILE *fp )
{
	gReactionNo	= 0;
	Iterate( gReactionList, EmitReactionsNumF77, fp );
/*	Iterate( gPAHReactionList, EmitReactionsNumF77, fp );*/
/*	Iterate( gSootReactionList, EmitReactionsNumF77, fp );*/
	fprintf( fp, "\n     &    , REND = %d )", 1 + gReactionNo );
}


static void EnumerateThirdBodyF77( FILE *fp )
{
	gThirdBodyNo = 0;
	Iterate( gUsedThirdBodyList, EmitThirdBodyNumF77, fp );
	fprintf( fp, "\n     &    , MEND = %d )", 1 + gThirdBodyNo );
}

static void EmitSpeciesF77( void *ptr, void *aux )
{
	SpeciesPtr	sPtr = (SpeciesPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gSpeciesNo % gNOfElementsLine == 0 ) {
		if ( gSpeciesNo % gNOfElementsBlock == 0 ) {
			fprintf( fp, "\n      INTEGER" );
			separator = ' ';
		}
		else {
			fprintf( fp, "\n     &    " );
		}
	}

	if ( sPtr->isSteadyState == FALSE ) {
		fprintf( fp, "%c %s%s", separator, gPrefixComputedF77, CSymbol(sPtr->name) );
	}
	else if ( sPtr->isSteadyState == TRUE ) {
		fprintf( fp, "%c %s%s", separator, gPrefixSteadyStateF77, CSymbol(sPtr->name) );
	}
	else {
		fprintf( fp, "\n# %s: unknown steady state flag \n", CSymbol(sPtr->name) );
	}
	++gSpeciesNo;
}


static void EmitReactionsF77( void *ptr, void *aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gReactionNo % gNOfElementsLine == 0 ) {
		if ( gReactionNo % gNOfElementsBlock == 0 ) {
			fprintf( fp, "\n      INTEGER" );
			separator = ' ';
		}
		else {
			fprintf( fp, "\n     &    " );
		}
	}
	
	fprintf( fp, "%c %s%s", separator, gPrefixReactionsF77, CSymbol( rPtr->label ));
	++gReactionNo;
}


static void EmitThirdBodyF77( void *ptr, void *aux )
{
	String		rPtr = (String)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gThirdBodyNo % gNOfElementsLine == 0 ) {
		if ( gThirdBodyNo % gNOfElementsBlock == 0 ) {
			fprintf( fp, "\n      INTEGER" );
			separator = ' ';
		}
		else {
			fprintf( fp, "\n     &    " );
		}
	}
	
	fprintf( fp, "%c %s%s", separator, gPrefixThirdBodyF77, CSymbol(rPtr) );
	++gThirdBodyNo;
}


static void EmitSpeciesNumF77( void *ptr, void *aux )
{
	SpeciesPtr	sPtr = (SpeciesPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gSpeciesNo % gNOfElementsLine == 0 ) {
		if ( gSpeciesNo % gNOfElementsBlock == 0 ) {
			if ( gSpeciesNo != 0 ) {
				fprintf( fp, " )" );
			}
			fprintf( fp, "\n      PARAMETER (" );
			separator = ' ';
		}
		else {
			fprintf( fp, "\n     &    " );
		}
	}
	
	if ( sPtr->isSteadyState == FALSE ) {
		fprintf( fp, "%c %s%s = %d", separator, gPrefixComputedF77, CSymbol(sPtr->name), 1 + gSpeciesNo++ );
	}
	else if ( sPtr->isSteadyState == TRUE ) {
		fprintf( fp, "%c %s%s = %d", separator, gPrefixSteadyStateF77, CSymbol(sPtr->name), 1 + gSpeciesNo++ );
	}
	else {
		fprintf( fp, "\n# %s: unknown steady state flag \n", CSymbol(sPtr->name) );
	}
}


static void EmitReactionsNumF77( void *ptr, void *aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gReactionNo % gNOfElementsLine == 0 ) {
		if ( gReactionNo % gNOfElementsBlock == 0 ) {
			if ( gReactionNo != 0 ) {
				fprintf( fp, " )" );
			}
			fprintf( fp, "\n      PARAMETER (" );
			separator = ' ';
		}
		else {
			fprintf( fp, "\n     &    " );
		}
	}
	
	fprintf( fp, "%c %s%s = %d", separator, gPrefixReactionsF77, CSymbol(rPtr->label), 1 + gReactionNo++ );
}

static void EmitArrCoeffF77A( void *ptr, void *aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	char		buffer[32];
	const int	parsPerLine = 2;
	Double		a = 0.0, n = 0.0, e = 0.0;
	
	if ( gReactionNo % parsPerLine == 0 ) {
		if ( gReactionNo % gNOfElementsBlock == 0 ) {
			if ( gReactionNo != 0 ) {
				fprintf( fp, " /" );
			}
			fprintf( fp, "\n      DATA (A(IH),IH=%d,%d)  /", gReactionNo+1
						, MIN( gReactionNo+gNOfElementsBlock, gCounter->reactions ) );
			separator = ' ';
		}
		else {
			fprintf( fp, "\n     &    " );
		}
	}
	
	PrintArrCoeffsF77( rPtr, kK, &a, &n, &e );
	ConvertCDoubleToF77Double( a, buffer );
	fprintf( fp, "%c %s", separator, buffer );
	++gReactionNo;
}

static void EmitArrCoeffF77N( void *ptr, void *aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	const int	parsPerLine = 2;
	Double		a = 0.0, n = 0.0, e = 0.0;
	
	if ( gReactionNo % parsPerLine == 0 ) {
		if ( gReactionNo % gNOfElementsBlock == 0 ) {
			if ( gReactionNo != 0 ) {
				fprintf( fp, " /" );
			}
			fprintf( fp, "\n      DATA (N(IH),IH=%d,%d)  /", gReactionNo+1
						, MIN( gReactionNo+gNOfElementsBlock, gCounter->reactions ) );
			separator = ' ';
		}
		else {
			fprintf( fp, "\n     &    " );
		}
	}
	
	PrintArrCoeffsF77( rPtr, kK, &a, &n, &e );
	fprintf( fp, "%c %G", separator, n );
	++gReactionNo;
}

static void EmitArrCoeffF77E( void *ptr, void *aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	const int	parsPerLine = 2;
	Double		a = 0.0, n = 0.0, e = 0.0;
	
	if ( gReactionNo % parsPerLine == 0 ) {
		if ( gReactionNo % gNOfElementsBlock == 0 ) {
			if ( gReactionNo != 0 ) {
				fprintf( fp, " /" );
			}
			fprintf( fp, "\n      DATA (E(IH),IH=%d,%d)  /", gReactionNo+1
						, MIN( gReactionNo+gNOfElementsBlock, gCounter->reactions ) );
			separator = ' ';
		}
		else {
			fprintf( fp, "\n     &    " );
		}
	}
	
	PrintArrCoeffsF77( rPtr, kK, &a, &n, &e );
	fprintf( fp, "%c %.10G", separator, e );
	++gReactionNo;
}


static void EmitThirdBodyNumF77( void *ptr, void *aux )
{
	String			rPtr = (String)ptr;
	FILE			*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gThirdBodyNo % gNOfElementsLine == 0 ) {
		if ( gThirdBodyNo % gNOfElementsBlock == 0 ) {
			if ( gThirdBodyNo != 0 ) {
				fprintf( fp, " )" );
			}
			fprintf( fp, "\n      PARAMETER (" );
			separator = ' ';
		}
		else {
			fprintf( fp, "\n     &    " );
		}
	}
	
	fprintf( fp, "%c %s%s = %d", separator, gPrefixThirdBodyF77, CSymbol(rPtr), 1 + gThirdBodyNo++ );
}

/*-----------------------------*/
/* FORTRAN 90 Routines         */
/*-----------------------------*/

void SaveLabelsF90( void )
{
	FILE	*fp;
	char	fname[100];
	
	sprintf( fname, "%sF90.h", gOptions->base );
	if ( !( fp = fopen( fname, "w") ) ) { 
		fprintf( stderr, "# SaveLabelsF90: unable to open file %s\n", fname );
		exit(2);
	}
	
	fprintf( fp, "!-------------------------------------------------------------\n" );
	fprintf( fp, "! ======= %s =======\n", fname );
	fprintf( fp, "!-------------------------------------------------------------\n" );

/* 	ListSpeciesF90( fp ); */
/* 	ListReactionsF90( fp ); */
/* 	ListThirdBodyF90( fp ); */
	EnumerateSpeciesF90( fp );
	EnumerateReactionsF90( fp );
	EnumerateThirdBodyF90( fp );
	fprintf( fp, "\n\ninteger, parameter :: DP=kind(1.0d0)\n");
	fprintf( fp, "\n\n\treal(DP) :: A(REND),N(REND),E(REND)\n" );
	fprintf( fp, "\n\n\tinteger :: IH\n" );
	EnumerateArrCoeffsF90( fp );
	fprintf( fp, "\n" );

	fclose( fp );
}

static void ListSpeciesF90( FILE *fp )
{
	gSpeciesNo	= 0;
	Iterate( gSpeciesList, EmitSpeciesF90, fp );
	Emit( fp, " & \n     , SEND" );
}


static void ListReactionsF90( FILE *fp )
{
	gReactionNo	= 0;
	Iterate( gReactionList, EmitReactionsF90, fp );
	Emit( fp, " &  \n     , REND" );
}


static void ListThirdBodyF90( FILE *fp )
{
	gThirdBodyNo = 0;
	Iterate( gUsedThirdBodyList, EmitThirdBodyF90, fp );
	Emit( fp, " &  \n     , MEND" );
}

static void EnumerateSpeciesF90( FILE *fp )
{
	gSpeciesNo	= 0;
	Iterate( gSpeciesList, EmitSpeciesNumF90, fp );
	fprintf( fp, " & \n     , SEND = %d ", 1 + gSpeciesNo );
}


static void EnumerateArrCoeffsF90( FILE *fp )
{
	gReactionNo	= 0;
	Iterate( gReactionList, EmitArrCoeffF90A, fp );
	fprintf( fp, " & \n     /" );
	gReactionNo	= 0;
	Iterate( gReactionList, EmitArrCoeffF90N, fp );
	fprintf( fp, " & \n     /" );
	gReactionNo	= 0;
	Iterate( gReactionList, EmitArrCoeffF90E, fp );
	fprintf( fp, " & \n     /" );
}


static void EnumerateReactionsF90( FILE *fp )
{
	gReactionNo	= 0;
	Iterate( gReactionList, EmitReactionsNumF90, fp );
	fprintf( fp, " & \n     , REND = %d ", 1 + gReactionNo );
}


static void EnumerateThirdBodyF90( FILE *fp )
{
	gThirdBodyNo = 0;
	Iterate( gUsedThirdBodyList, EmitThirdBodyNumF90, fp );
	fprintf( fp, " & \n     , MEND = %d ", 1 + gThirdBodyNo );
}

static void EmitSpeciesF90( void *ptr, void *aux )
{
	SpeciesPtr	sPtr = (SpeciesPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gSpeciesNo % gNOfElementsLine == 0 ) {
		if ( gSpeciesNo % gNOfElementsBlock == 0 ) {
			fprintf( fp, "\n     integer :: " );
			separator = ' ';
		}
		else {
			fprintf( fp, " & \n     " );
		}
	}

	if ( sPtr->isSteadyState == FALSE ) {
		fprintf( fp, "%c %s%s", separator, gPrefixComputedF77, CSymbol(sPtr->name) );
	}
	else if ( sPtr->isSteadyState == TRUE ) {
		fprintf( fp, "%c %s%s", separator, gPrefixSteadyStateF77, CSymbol(sPtr->name) );
	}
	else {
		fprintf( fp, "\n# %s: unknown steady state flag \n", CSymbol(sPtr->name) );
	}
	++gSpeciesNo;
}


static void EmitReactionsF90( void *ptr, void *aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gReactionNo % gNOfElementsLine == 0 ) {
		if ( gReactionNo % gNOfElementsBlock == 0 ) {
			fprintf( fp, "\n     integer :: " );
			separator = ' ';
		}
		else {
			fprintf( fp, " & \n     " );
		}
	}
	
	fprintf( fp, "%c %s%s", separator, gPrefixReactionsF77, CSymbol( rPtr->label ));
	++gReactionNo;
}


static void EmitThirdBodyF90( void *ptr, void *aux )
{
	String		rPtr = (String)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gThirdBodyNo % gNOfElementsLine == 0 ) {
		if ( gThirdBodyNo % gNOfElementsBlock == 0 ) {
			fprintf( fp, "\n     integer :: " );
			separator = ' ';
		}
		else {
			fprintf( fp, " & \n     " );
		}
	}
	
	fprintf( fp, "%c %s%s", separator, gPrefixThirdBodyF77, CSymbol(rPtr) );
	++gThirdBodyNo;
}


static void EmitSpeciesNumF90( void *ptr, void *aux )
{
	SpeciesPtr	sPtr = (SpeciesPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gSpeciesNo % gNOfElementsLine == 0 ) {
		if ( gSpeciesNo % gNOfElementsBlock == 0 ) {
			if ( gSpeciesNo != 0 ) {
				fprintf( fp, " " );
			}
			fprintf( fp, "\n     integer, parameter :: " );
			separator = ' ';
		}
		else {
			fprintf( fp, " & \n     " );
		}
	}
	
	if ( sPtr->isSteadyState == FALSE ) {
		fprintf( fp, "%c %s%s = %d", separator, gPrefixComputedF77, CSymbol(sPtr->name), 1 + gSpeciesNo++ );
	}
	else if ( sPtr->isSteadyState == TRUE ) {
		fprintf( fp, "%c %s%s = %d", separator, gPrefixSteadyStateF77, CSymbol(sPtr->name), 1 + gSpeciesNo++ );
	}
	else {
		fprintf( fp, "\n# %s: unknown steady state flag \n", CSymbol(sPtr->name) );
	}
}


static void EmitReactionsNumF90( void *ptr, void *aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gReactionNo % gNOfElementsLine == 0 ) {
		if ( gReactionNo % gNOfElementsBlock == 0 ) {
			if ( gReactionNo != 0 ) {
				fprintf( fp, " " );
			}
			fprintf( fp, "\n      integer, parameter :: " );
			separator = ' ';
		}
		else {
			fprintf( fp, "  & \n     " );
		}
	}
	
	fprintf( fp, "%c %s%s = %d", separator, gPrefixReactionsF77, CSymbol(rPtr->label), 1 + gReactionNo++ );
}

static void EmitArrCoeffF90A( void *ptr, void *aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	char		buffer[32];
	const int	parsPerLine = 2;
	Double		a = 0.0, n = 0.0, e = 0.0;
	
	if ( gReactionNo % parsPerLine == 0 ) {
		if ( gReactionNo % gNOfElementsBlock == 0 ) {
			if ( gReactionNo != 0 ) {
				fprintf( fp, " /" );
			}
			fprintf( fp, "\n      DATA (A(IH),IH=%d,%d)  /", gReactionNo+1
						, MIN( gReactionNo+gNOfElementsBlock, gCounter->reactions ) );
			separator = ' ';
		}
		else {
			fprintf( fp, " & \n     " );
		}
	}
	
	PrintArrCoeffsF77( rPtr, kK, &a, &n, &e );
	ConvertCDoubleToF77Double( a, buffer );
	fprintf( fp, "%c %s", separator, buffer );
	++gReactionNo;
}

static void EmitArrCoeffF90N( void *ptr, void *aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	const int	parsPerLine = 2;
	Double		a = 0.0, n = 0.0, e = 0.0;
	
	if ( gReactionNo % parsPerLine == 0 ) {
		if ( gReactionNo % gNOfElementsBlock == 0 ) {
			if ( gReactionNo != 0 ) {
				fprintf( fp, " /" );
			}
			fprintf( fp, "\n      DATA (N(IH),IH=%d,%d)  /", gReactionNo+1
						, MIN( gReactionNo+gNOfElementsBlock, gCounter->reactions ) );
			separator = ' ';
		}
		else {
			fprintf( fp, " & \n     " );
		}
	}
	
	PrintArrCoeffsF77( rPtr, kK, &a, &n, &e );
	fprintf( fp, "%c %G", separator, n );
	++gReactionNo;
}

static void EmitArrCoeffF90E( void *ptr, void *aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	char		separator = ',';
	const int	parsPerLine = 2;
	Double		a = 0.0, n = 0.0, e = 0.0;
	
	if ( gReactionNo % parsPerLine == 0 ) {
		if ( gReactionNo % gNOfElementsBlock == 0 ) {
			if ( gReactionNo != 0 ) {
				fprintf( fp, " /" );
			}
			fprintf( fp, "\n      DATA (E(IH),IH=%d,%d)  /", gReactionNo+1
						, MIN( gReactionNo+gNOfElementsBlock, gCounter->reactions ) );
			separator = ' ';
		}
		else {
			fprintf( fp, " & \n      " );
		}
	}
	
	PrintArrCoeffsF77( rPtr, kK, &a, &n, &e );
	fprintf( fp, "%c %.10G", separator, e );
	++gReactionNo;
}


static void EmitThirdBodyNumF90( void *ptr, void *aux )
{
	String			rPtr = (String)ptr;
	FILE			*fp = (FILE *)aux;
	char		separator = ',';
	
	if ( gThirdBodyNo % gNOfElementsLine == 0 ) {
		if ( gThirdBodyNo % gNOfElementsBlock == 0 ) {
			if ( gThirdBodyNo != 0 ) {
				fprintf( fp, " " );
			}
			fprintf( fp, "\n     integer, parameter :: " );
			separator = ' ';
		}
		else {
			fprintf( fp, " & \n     " );
		}
	}
	
	fprintf( fp, "%c %s%s = %d", separator, gPrefixThirdBodyF77, CSymbol(rPtr), 1 + gThirdBodyNo++ );
}


char *ShortSymbol( const char *src, char *buffer )
{
	char				*isomere;

	strcpy( buffer, src );
	if ( isomere = strrchr( buffer, '-' ) ) {
		if ( ( isomere - buffer ) > 3 || ( strlen( buffer ) > 8 && ( isomere - buffer ) > 3 ) ) {
			*isomere = '\0';
		}
	}
	
	return buffer;
}
