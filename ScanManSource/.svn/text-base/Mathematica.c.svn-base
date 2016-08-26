/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 *
 *	Mathematica.[ch]
 *				Functions related to Mathematica output.
 */

#ifndef __ArrayManager__
#include "ArrayManager.h"
#endif

#ifdef qUseDump
#pragma load "ScanMan.dump"
#else
#include "ScanMan.h"
#endif

#undef NDEBUG
#include <assert.h>
#include "Mathematica.h"

#undef qDebug
#define qSpecial

/*
	GLOBALS
*/
static const FPUType	kR = 8.3147;				/*	in		J / (mole K)		*/
static MatrixPtr		gNu;		/*	pointer to matrix of stoichiometric coeffs	*/


int MathematicaNumber( char *m, FPUType x )
{
	char	buffer[40];
	char	*src = buffer;
	char	*dest = m;
	char	negative = FALSE;
	const FPUType	kTiny = 1.0e-12;
	
	if ( fabs(x) < kTiny ) {						/*	We don't want tiny numbers.	*/
		
		m[0] = '0';
		m[1] = '\0';

	} else {
	
		sprintf( buffer, "%g", x );
	
		while ( *src && toupper(*src) != 'E' )
			*dest++ = *src++;						/*	Copy the mantissa.			*/
		if ( *src++ ) {
			*dest++ = '*';
			*dest++ = '1';
			*dest++ = '0';
			*dest++ = '^';
			if ( *src == '-' ) {
				++src;
				*dest++ = '(';
				*dest++ = '-';
				negative = TRUE;
			} else if ( *src == '+' )				/*	Skip '+' sign in exponent.	*/
				++src;
			while ( *src ) *dest++ = *src++;		/*	Copy the exponent.			*/
			if ( negative ) *dest++ = ')';
		}
		*dest = '\0';
	}
	
	return 0;
}


static void MathematicaExpression( String math, const String fc )
{
	const char	*ptr = fc - 1;
	const char	*la = fc;							/*	Look-ahead pointer			*/
	char		*dest = math;
	char 		tmp[2] = { '\0', '\0' };
	int			paren_open_expected = FALSE;
	int			paren_close_expected = 0;
	
	while ( *++ptr ) {
		switch ( *ptr ) {
			case 'T':
				la = ptr + 1;
				if ( isspace(*la) || isdigit(*la) )
					*dest++ = 't';
				else {
					tmp[0] = *la;
					if ( strpbrk( tmp, "().,;+-*/" ) )	*dest++ = 't';
					else								*dest++ = *ptr;
				}
				break;
			case 'e':
				la = ptr + 1;
				if ( *la == 'x' && *(la + 1) == 'p' ) {
					*dest++ = 'E';
					paren_open_expected = TRUE;
				} else {
					*dest++ = *ptr;
				}
				break;
			case '(':
				if ( paren_open_expected ) {
					paren_open_expected = FALSE;
					paren_close_expected = 1;
					*dest++ = '[';
				} else {
					if (paren_close_expected) ++paren_close_expected;
					*dest++ = *ptr;
				}
				break;
			case ')':
				if ( paren_close_expected ) {
					if ( paren_close_expected == 1 ) {
						paren_close_expected = 0;
						*dest++ = ']';
					} else {
						--paren_close_expected;
						*dest++ = *ptr;
					}
				} else {
					*dest++ = *ptr;
				}
				break;
			default:
				*dest++ = *ptr;
		}
	}	/*	while ( *++ptr )	*/
	*dest = '\0';
}


void MathematicaSymbol( char dest[80], const char *src )
{
	char *ptr = dest;
	
	--src;
	while ( *++src ) {
		switch( *src ) {
			case ',':
			case '.':
			case ';':
			case ':':
				/*	remove these characters	*/
				break;
			case '-':
				*ptr++ = 'X';
				break;
			default:
				*ptr++ = tolower( *src );
		}
	}
	*ptr = '\0';
}


static void SpeciesID( void *ptr, void *aux )
{
	char		buffer[80];
	SpeciesPtr	sPtr = (SpeciesPtr)ptr;
	FILE		*fp = aux;
	
	MathematicaSymbol( buffer, sPtr->name );
	
	fprintf( fp, "%s = %d;\n", buffer, sPtr->number + 1 );
#ifdef qSpecial
	if ( strcmp( buffer, "nXc3h7" ) == 0 )
		fprintf( fp, "c3h7 = %d;\n", sPtr->number + 1 );
#endif
}


static void MolarWeight( void *ptr, void *aux )
{
	char		buffer[80];
	SpeciesPtr	sPtr = (SpeciesPtr)ptr;
	FILE		*fp = aux;
	
	MathematicaSymbol( buffer, sPtr->name );	
	fprintf( fp, "mw[[%s]] = %.3f;\n", buffer, sPtr->molarMass );
}


static void InsertCoefficients( Pointer ptr, Pointer aux )
{
	int			i;
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	MatrixPtr	nu = (MatrixPtr)aux;
	
	for ( i = 0; i < rPtr->numberOfSpecies; ++i ) {
		PutMElement( rPtr->speciesCoeff[i], rPtr->speciesNumber[i], rPtr->id, nu );
	}
}


static void StoichiometricCoefficients( MatrixPtr nu )
{
	int		i, n;
	Double	*ptr = NULL;
	
	Iterate( gReactionList, InsertCoefficients, nu );

	/*	Currently the convention is, that the speciesCoeffs are positive
		for the left sides and negative for the right sides.  Therefore
		the matrix of stoichiometric coefficients has to change sign to
		be right.
	*/
	ptr = nu->mat[0];
	n = nu->rows * nu->cols;
	for ( i = 0; i < n; ++i, ++ptr )
		if ( *ptr != 0.0 ) *ptr *= -1.0;

#ifdef qDebug
	{	AMPrintOptions	prnt;
		
		DefaultAMPOpts( &prnt );
		prnt.title = "Stoichiometric coefficients";
		prnt.rowLabel = NULL;
		prnt.format = "%2g";
		prnt.lineWidth = 132;
		PrintMatrix( nu, &prnt, stdout );
	}
#endif
}


static void ReactionLabels( Pointer ptr, Pointer aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = aux;
	
	fprintf( fp, "r%s = %d;", rPtr->label, rPtr->id + 1 );
	if ( rPtr->withThirdBody || rPtr->withLindemann )
		fputs( " (*", fp );
	if ( rPtr->withThirdBody ) {
		fprintf( fp, " inert=%d", rPtr->thirdBodyNumber );
	}
	if ( rPtr->withLindemann ) {
		fprintf( fp, " L=%d", rPtr->lindemannNumber );
	}
	if ( rPtr->withThirdBody || rPtr->withLindemann )
		fputs( " *)", fp );
	fputs( "\n", fp );
}


int LookupSpecies( Pointer ptr, Pointer find )
{
	int		item = *(int *)find;

	return ( ((SpeciesPtr)ptr)->number == item );
}


static void ThirdBodyInfo( Pointer ptr, Pointer aux )
{
	ThirdBodyPtr	tPtr = (ThirdBodyPtr)ptr;
	FILE			*fp = (FILE *)aux;
	int				i, n = tPtr->speciesNumber->len;
	IntVectorPtr	id = tPtr->speciesNumber;
	VectorPtr		coeff = tPtr->speciesCoeff;
	SpeciesPtr		species = NULL;
	
	n = tPtr->speciesNumber->len;
	if ( n != coeff->len )
		Warning( "Array lengths are not equal in ThirdBodyInfo()" );
	
	fprintf( fp, "eff%s = { ", tPtr->name );
	for ( i = 0; i < n; ++i ) {
		fprintf( fp, "%.2f", coeff->vec[i] );
		species = FindListItem( gSpeciesList, &id->vec[i], LookupSpecies );
		if ( species == NULL )
			FatalError( "Couldn't find required species" );
		if ( i != species->number )
			Warning( "Species ID & loop counter don't match in ThirdBodyInfo()" );
#ifdef qDebug
		fprintf( fp, " (* %s *)", species->name );
#endif
		fputs( (i < n-1) ? ", " : " };\n", fp );
		if ( (i+1) % 5 == 0  && i != n-1) fputs( "\n\t", fp );
	}
}


static void ThirdBodyConcentration( Pointer ptr, Pointer aux )
{
	ThirdBodyPtr	tPtr = (ThirdBodyPtr)ptr;
	FILE			*fp = (FILE *)aux;

	fprintf( fp, "\t%s = eff%s . conc;\n", tPtr->name, tPtr->name );
}

static void ThirdBodyList( Pointer ptr, Pointer aux )
{
	ThirdBodyPtr	tPtr = (ThirdBodyPtr)ptr;
	FILE			*fp = (FILE *)aux;

	fprintf( fp, ", %s", tPtr->name );
}


static void PrintArrheniusFunc( FPUType A, FPUType n, FPUType E, FILE *fp )
{
	char	buffer[40];

	MathematicaNumber( buffer, A );
	fputs( buffer, fp );
	if ( n > 0.0 )			fprintf( fp, " * t^%g", n );
	else if ( n < 0.0 )	fprintf( fp, " * t^(%g)", n );
	if ( E != 0.0 ) {
		MathematicaNumber( buffer, - 1.0e3 * E / kR );
		fprintf( fp, " * E^(%s / t)", buffer );
	}
}


static void RateCoefficients( Pointer ptr, Pointer aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = (FILE *)aux;
	
	fprintf( fp, "rk[%d][t_] := ", rPtr->id + 1 );

	if ( rPtr->withLindemann ) {
		
		String	fcString;
		String	out;
		
		fputs( "Module[ { kInfty, k0, phi, kL, fc, n, f },\n", fp );
		
		  fprintf( fp, "\t\t(* Reaction \"%s\" *)\n", rPtr->label );

		  fputs( "\tk0 = ", fp );
		  PrintArrheniusFunc( rPtr->a, rPtr->n, rPtr->e, fp );
		  fputs( ";\n", fp );
		
		  fputs( "\tkInfty = ", fp );
		  PrintArrheniusFunc( rPtr->aInf, rPtr->nInf, rPtr->eInf, fp );
		  fputs( ";\n", fp );
		  
		  fputs( "\tphi = (k0 / kInfty) (pressure / (10^6 r t)); (* [M] in mole/cm^3 *)\n", fp );
		  fputs( "\tkL = phi / (1.0 + phi);\n", fp );
		  
		  fcString = ListItem( gFcContentsList, rPtr->lindemannNumber );
		  if ( fcString == NULL ) FatalError( "Couldn't find required string for Fc" );
		  out = (char *)malloc( 2 * strlen( fcString ) + 1 );
		  if ( out == NULL ) FatalError( "memory exhausted (RateCoefficients)" );
		  memset( out, 0, 2 * strlen( fcString ) + 1 );
		  MathematicaExpression( out, fcString );
		  fprintf( fp, "\tfc = %s;\n", out );
		  free( out );
		  
		  fputs( "\tn = 0.75 - 1.27 * Log[10,fc];\n", fp );
		  fputs( "\tphi = Log[10,phi] / n;\n",		  fp );
		  fputs( "\tf = fc^(1.0/(1.0 + phi^2));\n",   fp );
		  
		  fputs( "\tf * kInfty * kL\n]", fp );
		
	} else if ( rPtr->withThirdBody ) {
		PrintArrheniusFunc( rPtr->a, rPtr->n, rPtr->e, fp );
	} else {												/*	Normal case			*/
		PrintArrheniusFunc( rPtr->a, rPtr->n, rPtr->e, fp );
	}
	fputs( "\n", fp );
}


static void ReactionRates( Pointer ptr, Pointer aux )
{
	ReactionPtr	rPtr = (ReactionPtr)ptr;
	FILE		*fp = aux;
	int			err;
	int			i, num_species = gNu->rows;
	Double		val;

#ifdef qDebug
	fprintf( fp, "\t(* %3s *) w[[%d]] =",  rPtr->label, rPtr->id + 1 );
#else
	fprintf( fp, "\tw[[%d]] =", rPtr->id + 1 );
#endif
	for ( i = 0; i < num_species; ++i ) {
		err = GetMElement( &val, i, rPtr->id, gNu );
		assert( err == kIndexOK );
		if ( val < 0.0 ) {
			val = fabs(val);
			if ( val == 1.0 )	fprintf( fp, " conc[[%d]]", i + 1 );
			else				fprintf( fp, " conc[[%d]]^%g", i + 1, val );
		}
	}
	if ( rPtr->withThirdBody ) {
		ThirdBodyPtr tPtr = ListItem( gThirdBodyList, rPtr->thirdBodyNumber+1 );
		if ( tPtr == NULL ) FatalError( "Didn't find requested \"third body\" record" );
		fprintf( fp, " %s", tPtr->name );
	}
	fputs( ";\n", fp );
}


void WriteMathematicaFile( void )
{
	const char	*suffix = ".m";
	char		*fName = NULL;
	int			numSpecies, numReactions, numThirdBodies;
	int			i, j;										/*	Simple counters		*/
	Double		val;
	MatrixPtr	nu;
	FILE		*fp = NULL;
	
	fName = (char *)malloc( strlen(gOptions->base) + strlen(suffix) + 1 );
	if ( fName == NULL ) FatalError( "memory exhausted (WriteMathematicaFile)" );
	sprintf( fName, "%s%s", gOptions->base, suffix );
	
	fp = fopen( fName, "w" );
	if ( fp == NULL ) FatalError( "Couldn't open Mathematica file" );
	
	fputs( "Begin[ \"ScanMan`\" ]\n\n", fp );
	fputs( "r = 8.3147; (* J/(mole K) *)\n", fp );
	
	/*	Species Information.														*/
	fputs( "\n(* Species information *)\n", fp );
	numSpecies = CountItems( gSpeciesList );
	fprintf( fp, "numSpecies = %d;\n", numSpecies );
	Iterate( gSpeciesList, SpeciesID, fp );
	fputs( "\nmw = Table[ 0, {numSpecies} ];\n", fp );
	Iterate( gSpeciesList, MolarWeight, fp );
	fputs( "\ny = Table[ 0, {numSpecies} ];\n", fp );
	fputs( "x = Table[ 0, {numSpecies} ];\n", fp );
	
	/*	Third body information.														*/
	fputs( "\n(* Third body \"efficiency\" arrays *)\n", fp );
	numThirdBodies = CountItems( gThirdBodyList );
	fprintf( fp, "numThirdBodies = %d;\n", numThirdBodies );
	Iterate( gThirdBodyList, ThirdBodyInfo, fp );
	
	/*	Reaction information.														*/
	fputs( "\n(* Reaction information *)\n", fp );
	numReactions = CountItems( gReactionList );
	fprintf( fp, "numReactions = %d;\n", numReactions );
	
	Iterate( gReactionList, ReactionLabels, fp );
	
	/*	Rate constants.																*/
	fputs( "\n(* Rate coefficients *)\n", fp );
	fputs( "Remove[ rk ]\n", fp );
	Iterate( gReactionList, RateCoefficients, fp );
	
	nu = NewMatrix( numSpecies, numReactions, kColumnPointers );
	StoichiometricCoefficients( nu );
	
	/*	All coefficients.															*/
	fputs( "\nnu = { { ", fp );
	for ( i = 0; i < numSpecies; ++i ) {
		for ( j = 0; j < numReactions; ++j ) {
			GetMElement( &val, i, j, nu );
			fprintf( fp, "%g%s", val, (j < numReactions-1) ? ", " : " }" );
		}
		fputs( (i < numSpecies-1) ? ",\n { " : "\n", fp );
	}
	fputs( "};\n", fp );

	/*	"Production" coefficients.													*/
	fputs( "pnu = { { ", fp );
	for ( i = 0; i < numSpecies; ++i ) {
		for ( j = 0; j < numReactions; ++j ) {
			GetMElement( &val, i, j, nu );
			if ( val < 0.0 ) val = 0.0;
			fprintf( fp, "%g%s", val, (j < numReactions-1) ? ", " : " }" );
		}
		fputs( (i < numSpecies-1) ? ",\n { " : "\n", fp );
	}
	fputs( "};\n", fp );

	/*	"Consumption" coefficients.													*/
	fputs( "cnu = { { ", fp );
	for ( i = 0; i < numSpecies; ++i ) {
		for ( j = 0; j < numReactions; ++j ) {
			GetMElement( &val, i, j, nu );
			if ( val > 0.0 )		val = 0.0;
			else if ( val < 0.0 )	val *= -1.0;
			fprintf( fp, "%g%s", val, (j < numReactions-1) ? ", " : " }" );
		}
		fputs( (i < numSpecies-1) ? ",\n { " : "\n", fp );
	}
	fputs( "};\n", fp );

	/*	Reaction rates.																*/
	fputs( "\n(* Reaction rates *)\n", fp );
	fputs( "ReactionRates[ conc_List, k_List ] := Module[ { w", fp );
	Iterate( gThirdBodyList, ThirdBodyList, fp );
	fputs( " },\n", fp );
	Iterate( gThirdBodyList, ThirdBodyConcentration, fp );
		
	fputs( "\tw = Table[ 0, {Length[k]}];\n", fp );
	gNu = nu;
	Iterate( gReactionList, ReactionRates, fp );
	gNu = NULL;
	fputs( "\tw * k\n]\n", fp );
	
	DisposeMatrix( nu );
	
	fputs( "\nEndAdd[];\n", fp );
	fclose( fp );
#ifdef applec
	fsetfileinfo( fName, 'OMEG', 'TEXT' );
#endif
	free( fName );
}
