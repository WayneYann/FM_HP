/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 */

#ifdef qUseDump
#pragma load "ScanMan.dump"
#else
#include "ScanMan.h"
#endif


#undef DEBUGATOMCONTROL
#undef SHOWALLINPUT

void UpperString( char *string )
{
	while ( *string ) {
		*string = toupper( *string );
		++string;
	}
}

void LowerString( char *string )
{
	while ( *string ) {
		*string = tolower( *string );
		++string;
	}
}

char *MatchAtom( char *yytext, int yyleng, int mode )
{
	static char stack[4];
	static int currentSize = 0;
	static char atom[2];
	
#ifdef DEBUGATOMCONTROL
	fprintf(stderr, "start to match: '%s' in mode %d\n", yytext, mode );
#endif

	ShiftStack( yytext, yyleng, &currentSize, stack );
	while ( currentSize > 1 ) {
/*  first try to match the first two chars of the stack,  */
		atom[0] = stack[currentSize-1];
		atom[1] = stack[currentSize-2];
#ifdef DEBUGATOMCONTROL
		fprintf(stderr, "now i try to match the first two chars of the stack: '%s'\n", atom );
#endif
		if ( MatchItem( atom, &currentSize ) ) {  /*  first two chars matched  */
			
		}
		else {
/*  and if it doesn't work try to match only the first char of the stack  */
			atom[0] = stack[currentSize-1];
			atom[1] = '\0';
#ifdef DEBUGATOMCONTROL
			fprintf(stderr, "both chars failed, try only the first: '%s'\n", atom );
#endif
			if ( MatchItem( atom, &currentSize ) ) {  /*  first char matched  */
			
			}
			else {
				ScanError( "i can't match the atom '%s'", atom, 0 );
				currentSize = 0;
				return NULL;
			}
		}
	}
			
	/*  and if the possible atom is completely in the stack, 
		match the last char, that is now the first, because
		the stack has been shifted since the last match                    */
	if ( mode == 1 && currentSize > 0 ){
		atom[0] = stack[currentSize-1];
		atom[1] = '\0';
#ifdef DEBUGATOMCONTROL
		fprintf(stderr, "match the last char  of the possible atom: '%s'\n", atom );
#endif
		if ( MatchItem( atom, &currentSize ) ) {  /*  first char matched  */
		
		}
		else {
			ScanError( "i can't match the atom '%s'", atom, 0 );
			currentSize = 0;
			return NULL;
		}
	}
	return atom;
}

int ShiftStack( char *yytext, int yyleng, int *currentSize, char *stack )
{
	unsigned int i;
	/*  first shift   */
	for ( i = MIN(3, *currentSize); i > 0; --i) {
		stack[i] = stack[i-1];
	}
	stack[0] = toupper( yytext[yyleng-1] );
	++*currentSize;
	return *currentSize;
}

int MatchItem( char *atom, int *currentSize )
{	
	int 			matched = 0;
	LinkPtr 		matchedLink = NULL;
	AtomsPtr		atomLink = NULL;

	if ( ListIterator( gAllowAtomList, FindString, atom ) ) {
		matched = 1;
		*currentSize -= strlen( atom );  /*  reduce  */
		if ( !( matchedLink = ListIterator( gAtomsList, FindAtom, atom ) ) ) {
			atomLink = AddAtom( gAtomsList, atom );
			++gComposition->vec[atomLink->number];
#ifdef DEBUGATOMCONTROL
			fprintf(stderr, "i've found the allowed atom '%s' and added it\n", atom );
#endif
		}
		else {
			/* set number of atoms of this kind equal to one and correct later if there 
			   is a higher number found  */
			atomLink = (AtomsPtr)matchedLink->item;
			++gComposition->vec[atomLink->number];	
#ifdef DEBUGATOMCONTROL
			fprintf(stderr, "i've found the allowed atom '%s' that is already in the list\n", atom );
#endif
		}
	}

	return matched;
}

void SkipWhites( char * label )
{
	unsigned int i, j = 0, length = 0;
	length = strlen( label );
	for ( i = 0; i < strlen( label ); ++i ) {
		if ( isspace( label[i] ) ) {
			--length;
		}
		else {
			label[j++] = label[i];
		}
	}
	label[length] = '\0';
}

void ClearComposition( IntVectorPtr composition )
{
	int i;
	for ( i = 0; i < composition->len; ++i ) {
		composition->vec[i] = 0;
	}
}

void LineCounter( void )
{
	static int	lineNumber = 0;
	int i;
	int	len;
	
#ifdef applec
	if ( gNumberOfLines > lineNumber ) {
	    lineNumber = gNumberOfLines;
		SpinCursor( 32 );
	}
#endif

	if ( gOptions->progress ) {
		fprintf( stderr, "yytext = %s\n", ( char * )yytext );
	}
	for ( i = 0; i < yyleng; ++i ) {
		if ( yytext[i] == '\n' ) {
			++gNumberOfLines;
			gLineContents[0] = '\0';
		}
		else {
			if ( ( len = strlen( gLineContents ) ) < LINELENGTH ) {
				gLineContents[len] = yytext[i];
				gLineContents[len+1] = '\0';
			}
		}
	}
	if ( gOptions->progress ) {
		fprintf( stderr, "gLineContents = %s\n", gLineContents );
	}
}

void AdjustNumberOfLines( void )
{
	int i;
	
	gLineContents[strlen( gLineContents ) - yyleng] = '\0';
	for ( i = 0; i < yyleng; ++i ) {
		if ( yytext[i] == '\n' ) {
			--gNumberOfLines;
		}
	}
}

void FillThirdBodyArray( int speciesNumber, Double Mcoefficient )
{
	ThirdBodyPtr	thirdBody = NULL;
	
	thirdBody = ( ThirdBodyPtr ) gThirdBodyList->lastItem->item;
	if ( thirdBody->set->vec[speciesNumber] == 0 ) {
	  thirdBody->speciesNumber->vec[speciesNumber] = speciesNumber;
	  thirdBody->set->vec[speciesNumber] = 1;
	}
	else {
	  fprintf( stderr, "Warning: Thirdbody number %d has been set before to %g\n", speciesNumber, Mcoefficient );
	  thirdBody->speciesNumber->vec[speciesNumber] = speciesNumber;
	}
	thirdBody->speciesCoeff->vec[speciesNumber] = Mcoefficient;
}

void FillOtherThirdBodyArray( Double coefficient )
{
	ThirdBodyPtr	thirdBody = NULL;
	int i;
	
	thirdBody = ( ThirdBodyPtr ) gThirdBodyList->lastItem->item;
	for ( i = 0; i < thirdBody->speciesNumber->len; ++i ) {
	  /*		if ( thirdBody->speciesCoeff->vec[i] == 0.0 ) {*/
		if ( thirdBody->set->vec[i] == 0 ) {
			thirdBody->speciesNumber->vec[i] = i;
			thirdBody->set->vec[i] = 1;
			thirdBody->speciesCoeff->vec[i] = coefficient;
		}
	}
}

char *CopyString( char *to, char *from )
{
	if ( strlen( from ) < 256 ) {
		strcpy( to, from );
	}
	else {
		ScanError( "name '%s' too long\n# maximum size of a name is 255 chars", from, 0 );
		exit( 2 );
	}

	return to;
}

void ScanError( char *message, char *text, Double value )
{
	++yynerrs;
	gErrorOut = TRUE;

	fprintf( stderr, "\n# parse error:\n# " );  
	fprintf( stderr, "%s\n#\n# ", gLineContents );
	fprintf( stderr, message, text, value );
	fprintf( stderr, "\n#\n" );
	fprintf( stderr, "#-----------------------------------------------------------------\n" );
#ifdef applec
	fprintf( stderr, " File \"%s\"; Line %d\n", gOptions->inName, gNumberOfLines );
#else
	fprintf( stderr, " vi +%d %s\n", gNumberOfLines, gOptions->inName );
#endif
	fprintf( stderr, "#-----------------------------------------------------------------\n#\n\n\n" );
	if ( yynerrs > 10 ) {
		fprintf( stderr, "# make errors fewer\n" );
		fprintf( stderr, "#\n" );
		PrintLists();
		CleanUp();
		exit( 2 );
	}	
}
