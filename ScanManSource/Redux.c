/*
	Redux.c
	
	Reduction module for ScanMan
	
	2.12.92	(pt)
*/


#define SSNUMER
#undef	NDEBUG
#undef DEBUGCOLLFORW
#define NEWREDUX

#include <assert.h>


#ifdef qUseDump
#pragma load "ScanMan.dump"
#else
#include "ScanMan.h"
#endif

#include "Mathematica.h"
#include "Redux.h"


/*	Globals																		*/
static const char			*kOutOfMemory = "memory exhausted";
static ForwardReactionPtr	gForwardReacs = NULL;		/* forward reactions	*/
static ReactantsPtr		gReactants = NULL;			/* reactants (no inerts)*/
static ReductionPtr		gReduction = NULL;			/* reductions			*/
static GlobalReactionPtr	gGlobalReacs = NULL;		/* global reactions		*/
const FPUType				kTiny = 1.e-8;
static const int			kNameLen = 32;
static const int			kLineLength = 60;
static char				*gBuffer = NULL;			/* size = 4*kLineLength	*/
static char				*gFloatBuffer = NULL;		/* size =   kLineLength	*/

#undef DEBUG
#define PATCH

static ForwardReactionPtr NewForwardReaction( int nReacs )
{
	ForwardReactionPtr	ptr = NULL;
	
	ptr = (ForwardReactionPtr)malloc( sizeof(*ptr) );
	if ( !ptr ) FatalError( kOutOfMemory );
	memset( ptr, 0, sizeof(*ptr) );
	
	ptr->nReacs = nReacs;

	ptr->reac = (ReactionPtr *)malloc( nReacs * sizeof(*ptr->reac) );
	if ( !ptr->reac ) FatalError( kOutOfMemory );
	memset( ptr->reac, 0, nReacs * sizeof(*ptr->reac) );
	
	return ptr;
}


static void FreeForwardReaction( ForwardReactionPtr r )
{
	if ( !r ) return;

	free( r->reac );
	free( r );
}


static ReactantsPtr NewReactants( int nSpecs )
{
	ReactantsPtr	ptr = NULL;
	
	ptr = (ReactantsPtr)malloc( sizeof(*ptr) );
	if ( !ptr ) FatalError( kOutOfMemory );
	memset( ptr, 0, sizeof(*ptr) );
	
	ptr->nSpecs = nSpecs;
	
	ptr->spec = (SpeciesPtr *)malloc( nSpecs * sizeof(*ptr->spec) );
	if ( !ptr->spec ) FatalError( kOutOfMemory );
	memset( ptr->spec, 0, nSpecs * sizeof(*ptr->spec) );
	
	return ptr;
}


static void FreeReactants( ReactantsPtr r )
{
	if ( !r ) return;

	free( r->spec );
	free( r );
}


static ReductionPtr NewReduction( int nSpecies, int nReactions )
{
	int				n = nReactions + nSpecies;
	ReductionPtr	ptr = NULL;
	int				*i = NULL;
	
	ptr = (ReductionPtr)malloc( sizeof(*ptr) );
	if ( !ptr ) FatalError( kOutOfMemory );
	memset( ptr, 0, sizeof(*ptr) );
	
	ptr->nReactions = nReactions;
	ptr->nSpecies = nSpecies;

	ptr->reac = (int *)malloc( (nReactions + nSpecies) * sizeof(*ptr->reac) );
	if ( !ptr->reac ) FatalError( kOutOfMemory );
	
	ptr->spec = ptr->reac + nReactions;
	
	i = ptr->reac;
	while ( n-- ) {
		*i = -1;
		++i;
	}
	
	return ptr;
}


static void FreeReduction( ReductionPtr r )
{
	if ( !r ) return;

	free( r->reac );
	free( r );
}


static GlobalReactionPtr NewGlobalReaction( int nGlobalReacs, int nGlobalSpecies )
{
	GlobalReactionPtr	ptr = NULL;
	int					*i = NULL;
	int					n = nGlobalReacs + nGlobalSpecies;
	
	ptr = (GlobalReactionPtr)malloc( sizeof(*ptr) );
	if ( !ptr ) FatalError( kOutOfMemory );
	memset( ptr, 0, sizeof(*ptr) );
	
	ptr->nGlobalReacs = nGlobalReacs;
	ptr->nGlobalSpecies = nGlobalSpecies;

	ptr->spec = (int *)malloc( n * sizeof(*ptr->reac) );
	if ( !ptr->spec ) FatalError( kOutOfMemory );
	
	ptr->reac = ptr->spec + nGlobalSpecies;
	
	i = ptr->spec;
	while ( n-- ) {
		*i = -1;
		++i;
	}
	
	return ptr;
}


static void FreeGlobalReaction( GlobalReactionPtr r )
{
	if ( !r ) return;

	free( r->spec );
	free( r );
}


int Column( ReactionPtr reac )
{
	int			i, n = gForwardReacs->nForward;
	ReactionPtr	*ptr = gForwardReacs->reac;
	
	for ( i = 0; i < n; ++i ) {
		if ( reac == ptr[i] ) return i;
#		ifdef applec
		SpinCursor( 1 );
#		endif
	}
	
	return -1;
}


int Row( SpeciesPtr spec )
{
	/*	int Row( SpeciesPtr spec );
	
		Given a SpeciesPtr spec, this function returns the row
		index of the matrix of stoichiometric coefficients.
	*/
	int			i, n = gReactants->nReactants;
	SpeciesPtr	*ptr = gReactants->spec;
	
	for ( i = 0; i < n; ++i ) {
		if ( spec == ptr[i] ) return i;
#		ifdef applec
		SpinCursor( 1 );
#		endif
	}
	
	return -1;
}


static void RomanNumeral( int number, char *roman, int caps )
{
	int			i, j, len, digit, order = 0;
	char		*ptr = roman;
	char		string[5];
	const char	letter[4][2] = { {'I', 'V'}, {'X', 'L'}, {'C', 'D'}, {'M', '0'} };

	if ( number < 1 || number > 3999 ) {
		*ptr = '\0';
		return;
	}
	
	sprintf( string, "%d", number );
	len = strlen( string );

	order = len-1;

	for ( i = 0; i < len; ++i ) {						/* Start with the first digit.		*/
		digit = string[i];

		switch ( digit ) {
			case '1':
			case '2':
			case '3':
				for ( j = 0; j < digit - '0'; ++j ) {
					*ptr++ = letter[order][0];
				}
				break;
			case '4':
				*ptr++ = letter[order][0];
				*ptr++ = letter[order][1];
				break;
			case '5':
				*ptr++ = letter[order][1];
				break;
			case '6':
			case '7':
			case '8':
				*ptr++ = letter[order][1];
				for ( j = 0; j < digit - '5'; ++j ) {
					*ptr++ = letter[order][0];
				}
				break;
			case '9':
				*ptr++ = letter[order][0];
				*ptr++ = letter[order+1][0];
				break;
			case '0':
				break;
		}
		--order;
	}
	*ptr = '\0';										/* Terminate string.				*/

	if ( !caps ) LowerString( roman );
}


static int PTFindReaction( Pointer ptr, Pointer aux )
{	
	ReactionPtr reaction = (ReactionPtr)ptr;
	const char	*name = (const char *)aux;
	
	if ( strcmp( name, reaction->label ) == 0 ) {
		return TRUE;
	} 
	else {
		return FALSE;
	}
}


static int PTFindSpecies( Pointer ptr, Pointer aux )
{	
	SpeciesPtr spec = (SpeciesPtr)ptr;
	int			*id = (int *)aux;
	
	if ( *id == spec->number ) {
		return TRUE;
	} 
	else {
		return FALSE;
	}
}


static int PTFindThirdBody( Pointer ptr, Pointer aux )
{	
	ThirdBodyPtr	tb = (ThirdBodyPtr)ptr;
	int				*id = (int *)aux;
	
	if ( *id == tb->id ) {
		return TRUE;
	} 
	else {
		return FALSE;
	}
}


static void FixPointers( Pointer ptr, Pointer aux )
{
#	ifdef applec
#	pragma unused ( aux )
#	endif
	
	ReactionPtr	reac = (ReactionPtr)ptr;
	const char	*label = reac->label;
	int			len = strlen( label );
	char		*partner = NULL;
	
	switch ( label[len-1] ) {
		case 'f':
			/*	find backward reaction
			*/
			reac->forward = reac;

			partner = (char *)malloc( len+1 );
			if ( !partner ) FatalError( kOutOfMemory );
			strcpy( partner, label );
			partner[len-1] = 'b';
			
			if ( !( reac->backward = FindListItem( gReactionList, partner, PTFindReaction ) ) ) {
				char	errorString[128];
				
				sprintf( errorString, "no backward reaction for reaction \"%s\"", label );
				FatalError( errorString );
			}
			free( partner );
			break;
		case 'b':
			/*	find forward reaction
			*/
			reac->backward = reac;

			partner = (char *)malloc( len+1 );
			if ( !partner ) FatalError( kOutOfMemory );
			strcpy( partner, label );
			partner[len-1] = 'f';
			
			if ( !( reac->forward = FindListItem( gReactionList, partner, PTFindReaction ) ) ) {
				char	errorString[128];
				
				sprintf( errorString, "no forward reaction for reaction \"%s\"", label );
				FatalError( errorString );
			}
			free( partner );
			break;
		default:
			reac->forward = reac;
			reac->backward = NULL;		/* no backward reaction */
			break;
	}
}


static void PrintForwardBackward( Pointer ptr, Pointer aux )
{
	ReactionPtr	reac = (ReactionPtr)ptr;
	const char	*label = reac->label;
	const char	*fLabel = reac->forward->label;
	const char	*bLabel = ( reac->backward ) ? reac->backward->label : NULL;
	FILE		*out = (aux) ? (FILE *)aux : stdout;
	
	if ( bLabel ) {
		if ( strcmp( label, fLabel ) == 0 ) {
			fprintf( out, "reaction \"%s\": backward reaction labelled \"%s\".\n",
				label, bLabel );
		}
		else {
			fprintf( out, "reaction \"%s\": forward reaction labelled \"%s\".\n",
				label, fLabel );
		}
	}
	else {
		fprintf( out, "reaction \"%s\": no backward reaction.\n",
			label );
	}
}


void FixForwardBackwardsPtrs( void )
{
	int		numReactions = CountItems( gReactionList );
	
	Iterate( gReactionList, FixPointers, NULL );
	
#	ifdef DEBUG
	Iterate( gReactionList, PrintForwardBackward, stderr );
#	endif
}


static int NumberOfGlobalReactions( void )
{
	int		numSpecies = CountItems( gSpeciesList );
	int		numAtoms = CountItems( gAtomsList );
	
	if ( !gReduction ) return -1;		/* we need number of steady state assumptions */
	
	return numSpecies - gReduction->nReducs - numAtoms;
}


int NumberOfSteadyStateSpecies( void )
{
	return (gReduction == NULL) ? -1 : gReduction->nReducs;
}

static int IsCommentLine( const char *s )
{
	char		c;
	
	while ( c = *s++ ) {
		if ( isspace( c ) ) continue;
		if ( c == '#' )	return TRUE;		/* line comment */
		else			return FALSE;		/* something else */
	};

	return TRUE;							/* empty line */
}


static void SteadyStateError( const char *fName, int line, const char *errorString )
{
	fprintf( stderr, "File \"%s\"; Line %d; # %s\n", fName, line, errorString );
	exit( 2 );
}


static int FindSteadyStateReac( Pointer ptr, Pointer aux )
{
	const char	*string = (const char *)aux;
	const char	*label = ((ReactionPtr)ptr)->label;
	char		c;
	
	while ( c = toupper(*string++) ) {
		if ( c != toupper(*label++) ) return FALSE;
	}
	
	/*	strings match up to here, check for possible trailing 'f'
	*/
	if ( (c = tolower(*label)) == '\0' || c == 'f' ) return TRUE;
	
	return FALSE;
}


static int FindSteadyStateSpecies( Pointer ptr, Pointer aux )
{
	const char		*string = (const char *)aux;
	const char		*name = ((SpeciesPtr)ptr)->name;
	
	return ( strcmp( string, name ) == 0 );
}


static int ReactionIsEliminated( ReactionPtr rItem, ReductionPtr red ) 
{
	int		n = red->nReducs;
	int		*reac = red->reac;

	while ( n-- ) {
		if ( gForwardReacs->reac[*reac++] == rItem ) return TRUE;
	}
	
	return FALSE;
}


int SpeciesIsSteadyState( SpeciesPtr sItem )
{
	int		n = gReduction->nReducs;
	int		*spec = gReduction->spec;

	while ( n-- ) {
		if ( gReactants->spec[*spec++] == sItem ) return TRUE;
	}
	
	return FALSE;
}


static void CollectGlobalSpecies( GlobalReactionPtr globalReacs, ReductionPtr red )
{
#ifdef applec
#pragma unused(red)
#endif
	int			i, nReactants = gReactants->nReactants;
	SpeciesPtr	*spec = gReactants->spec;
	int			*s = globalReacs->spec;
	
	for ( i = 0; i < nReactants; ++i ) {
		if ( SpeciesIsSteadyState( spec[i] /*, red*/ ) ) continue;
#		ifdef DEBUG
		fprintf( stderr, "# species %s is global.\n", spec[i]->name );
#	 	endif
		*s = i;
		++s;
	}
}



static void ReadSteadyStates( const char *fName, MatrixPtr nu )
{
	const int	kNameLen = 80;
	int			c, n, len, lineCount = 0, row, col;
	int 		nReactions = CountItems( gReactionList );
	int 		nSpecies = CountItems( gSpeciesList );
	int			nReactants = gReactants->nReactants;
	int			nGlobal = 0;
	FILE		*fp = NULL;
	char		*buffer = NULL;
	char		*spec = NULL, *reac = NULL, *ptr = NULL, *bptr = NULL;
	ReactionPtr	rItem = NULL;
	SpeciesPtr	sItem = NULL;
	LinkPtr		link = NULL;
	const char	*kGlobal = "GLOBAL";				/* key word for global reactions	*/
	
#ifdef NEWREDUX
/*	if ( !( fp = fopen( fName, "r" ) ) ) {*/
		if ( !( buffer = (char *)malloc( LINELENGTH+1 + 2 * kNameLen ) ) ) 
			FatalError( kOutOfMemory );
		
		gReduction = NewReduction( nSpecies, nReactions );
		
		gReduction->nReducs = NumberOfItems( gSteadyStateList );
/*		fprintf( stderr, "steady state list has %d items\n", gReduction->nReducs );*/
		for ( n = 0; n < gReduction->nReducs; ++n ) {
			link = Find_n_th_Item( n+1, gSteadyStateList );
			if ( link ) {
				spec = (char *) link->item;
			}
			else {
				fprintf( stderr, "error in link number %d\n", n );
				exit( 2 );
			}
			if ( !( sItem = FindListItem( gSpeciesList, spec, FindSteadyStateSpecies ) ) ) {
				fprintf( stderr, "illegal species '%s'\n", spec );
				exit(2);
			}
			
		if ( ( row = Row( sItem ) ) < 0 )
			fprintf( stderr, "something wrong in Redux: illegal species '%s'\n", sItem->name );

			gReduction->spec[n] = row;
			gReduction->reac[n] = 0;
		}


		return;
/*	}*/
#else
	if ( !( fp = fopen( fName, "r" ) ) ) {
		FatalError( "couldn't open steady state file" );
	}
	if ( !( buffer = (char *)malloc( LINELENGTH+1 + 2 * kNameLen ) ) ) 
		FatalError( kOutOfMemory );
	reac = buffer + LINELENGTH+1;
	spec = reac + kNameLen;
	memset( buffer, 0,  LINELENGTH+1 + 2 * kNameLen );
	
	gReduction = NewReduction( nSpecies, nReactions );
		
#endif
	
	while ( buffer = fgets( buffer, LINELENGTH, fp ) ) {
		++lineCount;
		if ( IsCommentLine( buffer ) ) continue;		/* skip comment/empty lines */

		if ( sscanf( buffer, "%s%s", spec, reac ) != 2 ) {
			/*	check for context switch
			*/
			UpperString( spec );
			if ( strcmp( spec, kGlobal ) != 0 ) {
				SteadyStateError( fName, lineCount, "expected species and reaction on one line.");
			}
			else {
#				ifdef DEBUG
				fputs( "# now reading global reactions.\n", stderr );
#				endif
				break;
			}
		}
		
#		ifdef DEBUG
		fprintf( stderr, "# species \"%s\", reaction \"%s\".\n", spec, reac );
#		endif
		UpperString( reac );
		UpperString( spec );
		
		if ( *reac != 'W' ) {
			SteadyStateError( fName, lineCount, "illegal reaction." );
		}
		ptr = reac+1;	/* strip 'w' */
		
		/*	strip trailing 'F' or 'B'
		*/
		len = strlen( reac )-1;
		if ( ( c = reac[len] ) == 'F' || c == 'B' ) reac[len] = '\0';
		
		if ( !( rItem = FindListItem( gReactionList, ptr, FindSteadyStateReac ) ) ) {
			SteadyStateError( fName, lineCount, "illegal reaction." );
		}
#		ifdef DEBUG
		fprintf( stderr, "# found reaction %s.\n", rItem->label );
#		endif
		/*	check species
		*/
		ptr = spec;
		
		if ( !( sItem = FindListItem( gSpeciesList, ptr, FindSteadyStateSpecies ) ) ) {
			SteadyStateError( fName, lineCount, "illegal species." );
		}
#		ifdef DEBUG
		fprintf( stderr, "# found species %s.\n", sItem->name );
#		endif
		
		if ( ( row = Row( sItem ) ) < 0 )
			SteadyStateError( fName, lineCount, "illegal species." );
		if ( ( col = Column( rItem ) ) < 0 )
			SteadyStateError( fName, lineCount, "illegal reaction." );
		
#		ifdef DEBUG
		fprintf( stderr, "# row = %d, col = %d.\n", row, col );
#		endif
		
		if ( fabs( nu->mat[row][col] ) < kTiny ) 
			SteadyStateError( fName, lineCount, "species not in reaction.");
		
		n = gReduction->nReducs;
		gReduction->spec[n] = row;
		gReduction->reac[n] = col;
		++gReduction->nReducs;
	}
	
#	ifdef DEBUG
	fprintf( stderr, "# %d steady state assumptions found.\n", 
		gReduction->nReducs );
#	endif
#ifdef PATCH
	free( buffer );
	fclose ( fp );

	return;

#else

	nGlobal = NumberOfGlobalReactions();
	gGlobalReacs = NewGlobalReaction( nGlobal, nReactants - gReduction->nReducs );
	
	CollectGlobalSpecies( gGlobalReacs, gReduction );
	
	nGlobal = 0;
	while ( TRUE ) {
		
		if ( !(buffer = fgets( buffer, LINELENGTH, fp ) ) )
			SteadyStateError( fName, lineCount, "unexpected end of file." );

		++lineCount;
		if ( IsCommentLine( buffer ) ) continue;		/* skip comment/empty lines */
		
		UpperString( buffer );
		bptr = buffer;

		while ( isspace( *bptr ) ) ++bptr;
		while ( *bptr ) {
	
			ptr = reac;
			while ( isalnum( *bptr ) ) *ptr++ = *bptr++;
			*ptr = '\0';
			
			ptr = reac;
			if ( *ptr == 'W' ) ++ptr;
			
			/*	strip trailing 'F' or 'B'
			*/
			len = strlen( ptr ) - 1;
			if ( ( c = ptr[len] ) == 'F' || c == 'B' ) ptr[len] = '\0';
			
			if ( !( rItem = FindListItem( gReactionList, ptr, FindSteadyStateReac ) ) ) {
				SteadyStateError( fName, lineCount, "illegal global reaction." );
			}
			if ( ReactionIsEliminated( rItem, gReduction ) )
				SteadyStateError( fName, lineCount, "reaction is eliminated." );
			
#			ifdef DEBUG
			fprintf( stderr, "# global reaction \"%s\".\n", rItem->label );
#			endif
				
			gGlobalReacs->reac[nGlobal] = Column( rItem );
			++nGlobal;
			if ( nGlobal == gGlobalReacs->nGlobalReacs ) break;
			
			while ( isspace( *bptr ) ) ++bptr;
		}
		if ( nGlobal == gGlobalReacs->nGlobalReacs ) break;
	}
				
	free( buffer );
	fclose ( fp );
#endif
}


static void CollectForwardReacs( Pointer ptr, Pointer aux )
{
	ForwardReactionPtr	fr = (ForwardReactionPtr)aux;
	ReactionPtr			reac = (ReactionPtr)ptr;
	const char			*label = reac->label;
	const char			*forward = reac->forward->label;
	
	if ( strcmp( label, forward ) == 0 ) {
		fr->reac[fr->nForward] = reac;
		++fr->nForward;
#		ifdef DEBUGCOLLFORW
		fprintf( stderr, "# reaction labelled \"%s\" is forward reaction.\n", label );
#		endif
	}
}


int FindSpeciesInReaction( Pointer ptr, Pointer aux )
{
	ReactionPtr	reac = (ReactionPtr)ptr;
	SpeciesPtr	spec = (SpeciesPtr)aux;
	int			i, n = reac->numberOfSpecies;
	
	for ( i = 0; i < n; ++i ) {
		if ( reac->speciesNumber[i] == spec->number ) return TRUE;
#		ifdef applec
		SpinCursor( 1 );
#		endif
	}
	
	return FALSE;
}


static void CollectReactants( Pointer ptr, Pointer aux )
{
	ReactantsPtr	rptr = (ReactantsPtr)aux;
	SpeciesPtr		spec = (SpeciesPtr)ptr;
	int				nCur = rptr->nReactants;
	
	if ( FindListItem( gReactionList, spec, FindSpeciesInReaction ) ) {
		rptr->spec[rptr->nReactants] = spec;
		++rptr->nReactants;
#		ifdef DEBUG
		fprintf( stderr, "# species %s is reactant with id %d.\n", spec->name, spec->number );
	}
	else {
		fprintf( stderr, "# species %s is inert with id %d.\n", spec->name, spec->number );
#		endif
	}
}


static MatrixPtr StoichCoeffMatrix( void )
{
	int			fCount = 0, i, j, k, n;
	int			numReacs = CountItems( gReactionList );
	int			numSpecs = CountItems( gSpeciesList );
	int			nReactants = 0, nForward = 0;
	ReactionPtr	r = NULL;
	SpeciesPtr	*s = NULL;
	MatrixPtr	nu = NULL;
	Double		**mat = NULL;
	
	gForwardReacs = NewForwardReaction( numReacs );
	Iterate( gReactionList, CollectForwardReacs, gForwardReacs );

	gReactants = NewReactants( numSpecs );
	Iterate( gSpeciesList, CollectReactants, gReactants );
	
	nForward = gForwardReacs->nForward;
	nReactants = gReactants->nReactants;
	
#	ifdef DEBUG
	fprintf( stderr, "# %d forward reactions found.\n", nForward );
	fprintf( stderr, "# %d reactants found.\n", nReactants );
#	endif
	
	nu = NewMatrix( nReactants, nForward, kRowPointers );
	mat = nu->mat;
	
	for ( i = 0; i < nForward; ++i ) {
		r = gForwardReacs->reac[i];
		s = gReactants->spec;
		n = r->numberOfSpecies;
		for ( j = 0; j < n; ++j ) {
			int		found = FALSE;
			int		specNum = r->speciesNumber[j];
			
#			ifdef applec
			SpinCursor( 1 );
#			endif
			for ( k = 0; k < nReactants; ++k ) {
				if ( specNum == s[k]->number ) {
					mat[k][i] = r->speciesCoeff[j];
					found = TRUE;
					break;
				}
			}
			assert( found );
		}
	}

	for ( i = 0; i < nReactants; ++i ) {
#		ifdef applec
		SpinCursor( 1 );
#		endif
		for ( j = 0; j < nForward; ++j ) 
			if ( mat[i][j] ) mat[i][j] = -mat[i][j];		/* flip signs */
	}
	
#	ifdef DEBUG
	{
		AMPrintOptions	prnt;
		
		DefaultAMPOpts( &prnt );
		prnt.title = "matrix of stoichiometric coefficients";
		prnt.sep = "\t";
		prnt.rowLabel = NULL;
		PrintMatrix( nu, &prnt, stdout );
		fputs( "\n", stdout );
	}
#	endif

	return nu;
}



static void Reduce( MatrixPtr nu )
{
	int			rows = nu->rows, cols = nu->cols;
	int			rowR, colR, row, col;
	int			i, nReducs = gReduction->nReducs;
	FPUType		factor = 0.0;
	Double		**a = nu->mat;
	
	/*	Eliminate the rates of the reactions that were specified together 
		with the steady state assumptions.										*/
	for ( i = 0; i < nReducs; ++i ) {	
		rowR = gReduction->spec[i];
		colR = gReduction->reac[i];

		/*	Check for spurious zero pivots (jg, 30.7.92).						*/
		if ( fabs(a[rowR][colR]) < kTiny )
			/*fprintf( stderr, "#error: reaction %s doesn't consume species %s\n"
						, gForwardReacs);*/
			FatalError( "Zero pivot element encountered" );

		for ( row = 0; row < rows; ++row ) {
			if ( ( row == rowR ) || ( fabs( a[row][colR] ) < kTiny ) ) continue;
			factor = -a[row][colR] / a[rowR][colR];
			for ( col = 0; col < cols; ++col)
				a[row][col] += factor * a[rowR][col];
		}

#		ifdef applec
		SpinCursor( 1 );
#		endif
	}
	
#	ifdef DEBUG
	{
		AMPrintOptions	prnt;
		
		DefaultAMPOpts( &prnt );
		prnt.title = "reduced matrix of stoichiometric coefficients";
		prnt.sep = "\t";
		prnt.rowLabel = NULL;
		PrintMatrix( nu, &prnt, stdout );
		fputs( "\n", stdout );
	}
#	endif
}

static char *GetSSIndexName( const char *name, char *buffer )
{
	char *wBuff = buffer;
	
	do {
		if ( *name == '-' ) {
			*wBuff++ = 'X';
		} 
		else if ( *name == '*' ) {
			*wBuff++ = 'Y';
		}
		else {
			*wBuff++ = *name;
		}
	} while ( *name++ != '\0' );
	
	return buffer;
}

static void SteadyStateRelations( MatrixPtr nu )
{
	int			nForward = gForwardReacs->nForward;
	int			numReduced = nForward - gReduction->nReducs;
	int			numSteady = gReduction->nReducs;
	int			i, j, l, m;
	const int	/**reac = gReduction->reac,*/ *spec = gReduction->spec;
	Double		**mat = nu->mat;
	FILE		*fp = stdout;
	char		*bptr = gBuffer;
	Flag			first = TRUE;
	Flag			tbUsed = FALSE;
	LinkPtr			speciesLink = NULL;
	SpeciesPtr 		species = NULL;
	LinkPtr			thirdBodyLink = NULL;
	char			*thirdBody = NULL;
	char			buffer[128];
	FILE			*fpOut;

	sprintf( buffer, "%s.cp", gOptions->base );
	fpOut = fopen( buffer, "w" );

	fprintf( fpOut, "#include \"FlameMaster.h\"\n" );
	fprintf( fpOut, "#include \"%s.h\"\n\n", gOptions->base );
	fprintf( fpOut, "void ComputeSteadyStates( Double *k, Double *c, Double *M )\n" );
	fprintf( fpOut, "{\n" );

	for ( i = 0; i < numSteady; ++i ) {
		int				s = spec[i];
		
		fprintf( fpOut, "\tc[s%s] = (", GetSSIndexName( gReactants->spec[s]->name, buffer ) );
		for ( j = 0; j < nForward; ++j ) {
			ReactionPtr		r = gForwardReacs->reac[j];
			/*	interpretation of matrix rows == conservation equations		*/
			if ( fabs( mat[s][j] ) > kTiny  ) {
				if ( mat[s][j] > 0.0 ) {
					if ( mat[s][j] > 1.5 ) {
						fprintf( fpOut, "\n\t\t%s k[%s%s] * %3.1f", (first)?"":"+", gPrefixReactions, CSymbol(r->label), mat[s][j] );
					}
					else {
						fprintf( fpOut, "\n\t\t%s k[%s%s]", (first)?"":"+", gPrefixReactions, CSymbol(r->label) );
					}
					first = FALSE;
					for ( l = 0; l < r->numberOfSpecies; ++l ) {
						if ( r->speciesCoeff[l] > 0.0 ) {
							speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &r->speciesNumber[l] );
							if ( !speciesLink ) FatalError( "something wrong in PrintRedux" );
							species = ( SpeciesPtr ) speciesLink->item;
							for ( m = 0; m < ( int )( r->speciesCoeff[l] + 0.5 ); ++m ) {
								fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
							}
						}
					}
					if ( r->withThirdBody && !r->withLindemann ) {
						thirdBodyLink = Find_n_th_Item( r->thirdBodyNumber+1, gUsedThirdBodyList );
						thirdBody = ( char * ) thirdBodyLink->item;
						fprintf( fpOut, " * M[m%s]", thirdBody );
						tbUsed = TRUE;
					}
				}
				if ( r->backward ) {
					if ( mat[s][j] < 0.0 ) {
						if ( mat[s][j] < -1.5 ) {
							fprintf( fpOut, "\n\t\t%s k[%s%s] * %3.1f", (first)?"":"+", gPrefixReactions, CSymbol(r->backward->label), -mat[s][j] );
						}
						else {
							fprintf( fpOut, "\n\t\t%s k[%s%s]", (first)?"":"+", gPrefixReactions, CSymbol(r->backward->label) );
						}
						first = FALSE;
						for ( l = 0; l < r->backward->numberOfSpecies; ++l ) {
							if ( r->backward->speciesCoeff[l] > 0.0 ) {
								speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &r->backward->speciesNumber[l] );
								if ( !speciesLink ) FatalError( "something wrong in PrintRedux" );
								species = ( SpeciesPtr ) speciesLink->item;
								for ( m = 0; m < ( int )(r->backward->speciesCoeff[l] + 0.5 ); ++m ) {
									fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
								}
							}
						}
						if ( r->backward->withThirdBody && !r->backward->withLindemann ) {
							thirdBodyLink = Find_n_th_Item( r->backward->thirdBodyNumber+1, gUsedThirdBodyList );
							thirdBody = ( char * ) thirdBodyLink->item;
							fprintf( fpOut, " * M[m%s]", thirdBody );
							tbUsed = TRUE;
						}
					}
				}
			}
		}
		
		first = TRUE;
		fprintf( fpOut, " ) / ( CatchZero( " );
		
		for ( j = 0; j < nForward; ++j ) {
			ReactionPtr		r = gForwardReacs->reac[j];
			/*	interpretation of matrix rows == conservation equations		*/
			if ( fabs( mat[s][j] ) > kTiny  ) {
				if ( mat[s][j] < 0.0 ) {
					if ( mat[s][j] < -1.5 ) {
						fprintf( fpOut, "\n\t\t%s k[%s%s] * %3.1f", (first)?"":"+", gPrefixReactions, CSymbol(r->label), -mat[s][j] );
					}
					else {
						fprintf( fpOut, "\n\t\t%s k[%s%s]", (first)?"":"+", gPrefixReactions, CSymbol(r->label) );
					}
					first = FALSE;
					for ( l = 0; l < r->numberOfSpecies; ++l ) {
						if ( r->speciesCoeff[l] > 0.0 ) {
							speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &r->speciesNumber[l] );
							if ( !speciesLink ) FatalError( "something wrong in PrintRedux" );
							species = ( SpeciesPtr ) speciesLink->item;
							if ( strcmp( gReactants->spec[s]->name, species->name ) == 0 ) {
								for ( m = 1; m < ( int )( r->speciesCoeff[l] + 0.5 ); ++m ) {
									fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
								}
							}
							else {
								for ( m = 0; m < ( int )( r->speciesCoeff[l] + 0.5 ); ++m ) {
									fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
								}
							}
						}
					}
					if ( r->withThirdBody && !r->withLindemann ) {
						thirdBodyLink = Find_n_th_Item( r->thirdBodyNumber+1, gUsedThirdBodyList );
						thirdBody = ( char * ) thirdBodyLink->item;
						fprintf( fpOut, " * M[m%s]", thirdBody );
						tbUsed = TRUE;
					}
				}
				if ( r->backward ) {
					if ( mat[s][j] > 0.0 ) {
						if ( mat[s][j] > 1.5 ) {
							fprintf( fpOut, "\n\t\t%s k[%s%s] * %3.1f", (first)?"":"+", gPrefixReactions, CSymbol(r->backward->label), mat[s][j] );
						}
						else {
							fprintf( fpOut, "\n\t\t%s k[%s%s]", (first)?"":"+", gPrefixReactions, CSymbol(r->backward->label) );
						}
						first = FALSE;
						for ( l = 0; l < r->backward->numberOfSpecies; ++l ) {
							if ( r->backward->speciesCoeff[l] > 0.0 ) {
								speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &r->backward->speciesNumber[l] );
								if ( !speciesLink ) FatalError( "something wrong in PrintRedux" );
								species = ( SpeciesPtr ) speciesLink->item;
								if ( strcmp( gReactants->spec[s]->name, species->name ) == 0 ) {
									for ( m = 1; m < ( int )( r->backward->speciesCoeff[l] + 0.5 ); ++m ) {
										fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
									}
								}
								else {
									for ( m = 0; m < ( int )( r->backward->speciesCoeff[l] + 0.5 ); ++m ) {
										fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
									}
								}
							}
						}
						if ( r->backward->withThirdBody && !r->backward->withLindemann ) {
							thirdBodyLink = Find_n_th_Item( r->backward->thirdBodyNumber+1, gUsedThirdBodyList );
							thirdBody = ( char * ) thirdBodyLink->item;
							fprintf( fpOut, " * M[m%s]", thirdBody );
							tbUsed = TRUE;
						}
					}
				}
			}
		}
		first = TRUE;
		fprintf( fpOut, " ) );\n\n" );
	}
	if ( !tbUsed ) {
		fprintf( fpOut, "M = M;\n" );
	}
	fprintf( fpOut, "}\n\n" );
	
	fprintf( fpOut, "int SteadyStatesFunc( const VectorPtr /*x*/, VectorPtr fVec, void *object )\n" );
	fprintf( fpOut, "{\n" );
	fprintf( fpOut, "\tint\t\tspeciesIn;\n" );
	fprintf( fpOut, "\tDouble\t*c;\n" );
	fprintf( fpOut, "\tDouble\t*k;\n" );
	fprintf( fpOut, "\tDouble\t*M;\n" );
	fprintf( fpOut, "\tDouble\t*f = fVec->vec;\n" );
	fprintf( fpOut, "\tSteadyStateInfoPtr ssInfo = ( SteadyStateInfoPtr )object;\n" );
	fprintf( fpOut, "\tssInfo->GetSteadyStateInfo( &c, &k, &M, &speciesIn );\n\n" );

	first = TRUE;
	tbUsed = FALSE;
	for ( i = 0; i < numSteady; ++i ) {
		int				s = spec[i];
		
		fprintf( fpOut, "\tf[s%s-speciesIn] = c[s%s] - ("
						, GetSSIndexName( gReactants->spec[s]->name, buffer )
						, GetSSIndexName( gReactants->spec[s]->name, buffer ) );
		for ( j = 0; j < nForward; ++j ) {
			ReactionPtr		r = gForwardReacs->reac[j];
			/*	interpretation of matrix rows == conservation equations		*/
			if ( fabs( mat[s][j] ) > kTiny  ) {
				if ( mat[s][j] > 0.0 ) {
					if ( mat[s][j] > 1.5 ) {
						fprintf( fpOut, "\n\t\t%s k[%s%s] * %3.1f", (first)?"":"+", gPrefixReactions, CSymbol(r->label), mat[s][j] );
					}
					else {
						fprintf( fpOut, "\n\t\t%s k[%s%s]", (first)?"":"+", gPrefixReactions, CSymbol(r->label) );
					}
					first = FALSE;
					for ( l = 0; l < r->numberOfSpecies; ++l ) {
						if ( r->speciesCoeff[l] > 0.0 ) {
							speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &r->speciesNumber[l] );
							if ( !speciesLink ) FatalError( "something wrong in PrintRedux" );
							species = ( SpeciesPtr ) speciesLink->item;
							for ( m = 0; m < ( int )( r->speciesCoeff[l] + 0.5 ); ++m ) {
								fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
							}
						}
					}
					if ( r->withThirdBody && !r->withLindemann ) {
						thirdBodyLink = Find_n_th_Item( r->thirdBodyNumber+1, gUsedThirdBodyList );
						thirdBody = ( char * ) thirdBodyLink->item;
						fprintf( fpOut, " * M[m%s]", thirdBody );
						tbUsed = TRUE;
					}
				}
				if ( r->backward ) {
					if ( mat[s][j] < 0.0 ) {
						if ( mat[s][j] < -1.5 ) {
							fprintf( fpOut, "\n\t\t%s k[%s%s] * %3.1f", (first)?"":"+", gPrefixReactions, CSymbol(r->backward->label), -mat[s][j] );
						}
						else {
							fprintf( fpOut, "\n\t\t%s k[%s%s]", (first)?"":"+", gPrefixReactions, CSymbol(r->backward->label) );
						}
						first = FALSE;
						for ( l = 0; l < r->backward->numberOfSpecies; ++l ) {
							if ( r->backward->speciesCoeff[l] > 0.0 ) {
								speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &r->backward->speciesNumber[l] );
								if ( !speciesLink ) FatalError( "something wrong in PrintRedux" );
								species = ( SpeciesPtr ) speciesLink->item;
								for ( m = 0; m < ( int )(r->backward->speciesCoeff[l] + 0.5 ); ++m ) {
									fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
								}
							}
						}
						if ( r->backward->withThirdBody && !r->backward->withLindemann ) {
							thirdBodyLink = Find_n_th_Item( r->backward->thirdBodyNumber+1, gUsedThirdBodyList );
							thirdBody = ( char * ) thirdBodyLink->item;
							fprintf( fpOut, " * M[m%s]", thirdBody );
							tbUsed = TRUE;
						}
					}
				}
			}
		}
		
		first = TRUE;
		fprintf( fpOut, " ) / ( CatchZero( " );
		
		for ( j = 0; j < nForward; ++j ) {
			ReactionPtr		r = gForwardReacs->reac[j];
			/*	interpretation of matrix rows == conservation equations		*/
			if ( fabs( mat[s][j] ) > kTiny  ) {
				if ( mat[s][j] < 0.0 ) {
					if ( mat[s][j] < -1.5 ) {
						fprintf( fpOut, "\n\t\t%s k[%s%s] * %3.1f", (first)?"":"+", gPrefixReactions, CSymbol(r->label), -mat[s][j] );
					}
					else {
						fprintf( fpOut, "\n\t\t%s k[%s%s]", (first)?"":"+", gPrefixReactions, CSymbol(r->label) );
					}
					first = FALSE;
					for ( l = 0; l < r->numberOfSpecies; ++l ) {
						if ( r->speciesCoeff[l] > 0.0 ) {
							speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &r->speciesNumber[l] );
							if ( !speciesLink ) FatalError( "something wrong in PrintRedux" );
							species = ( SpeciesPtr ) speciesLink->item;
							if ( strcmp( gReactants->spec[s]->name, species->name ) == 0 ) {
								for ( m = 1; m < ( int )( r->speciesCoeff[l] + 0.5 ); ++m ) {
									fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
								}
							}
							else {
								for ( m = 0; m < ( int )( r->speciesCoeff[l] + 0.5 ); ++m ) {
									fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
								}
							}
						}
					}
					if ( r->withThirdBody && !r->withLindemann ) {
						thirdBodyLink = Find_n_th_Item( r->thirdBodyNumber+1, gUsedThirdBodyList );
						thirdBody = ( char * ) thirdBodyLink->item;
						fprintf( fpOut, " * M[m%s]", thirdBody );
						tbUsed = TRUE;
					}
				}
				if ( r->backward ) {
					if ( mat[s][j] > 0.0 ) {
						if ( mat[s][j] > 1.5 ) {
							fprintf( fpOut, "\n\t\t%s k[%s%s] * %3.1f", (first)?"":"+", gPrefixReactions, CSymbol(r->backward->label), mat[s][j] );
						}
						else {
							fprintf( fpOut, "\n\t\t%s k[%s%s]", (first)?"":"+", gPrefixReactions, CSymbol(r->backward->label) );
						}
						first = FALSE;
						for ( l = 0; l < r->backward->numberOfSpecies; ++l ) {
							if ( r->backward->speciesCoeff[l] > 0.0 ) {
								speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &r->backward->speciesNumber[l] );
								if ( !speciesLink ) FatalError( "something wrong in PrintRedux" );
								species = ( SpeciesPtr ) speciesLink->item;
								if ( strcmp( gReactants->spec[s]->name, species->name ) == 0 ) {
									for ( m = 1; m < ( int )( r->backward->speciesCoeff[l] + 0.5 ); ++m ) {
										fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
									}
								}
								else {
									for ( m = 0; m < ( int )( r->backward->speciesCoeff[l] + 0.5 ); ++m ) {
										fprintf( fpOut, " * c[s%s]", GetSSIndexName( species->name, buffer ) );
									}
								}
							}
						}
						if ( r->backward->withThirdBody && !r->backward->withLindemann ) {
							thirdBodyLink = Find_n_th_Item( r->backward->thirdBodyNumber+1, gUsedThirdBodyList );
							thirdBody = ( char * ) thirdBodyLink->item;
							fprintf( fpOut, " * M[m%s]", thirdBody );
							tbUsed = TRUE;
						}
					}
				}
			}
		}
		first = TRUE;
		fprintf( fpOut, " ) );\n\n" );
	}
	if ( !tbUsed ) {
		fprintf( fpOut, "\tM = M;\n" );
	}
	fprintf( fpOut, "\treturn 0;\n" );
	fprintf( fpOut, "}\n\n" );

	fprintf( fpOut, "#ifndef MECHANISM\n" );
	fprintf( fpOut, "#define MECHANISM \"\"\n" );
	fprintf( fpOut, "#endif\n\n" );
	fprintf( fpOut, "void TReaction::CheckSteadyStatesMech( const char *mechName )\n" );
	fprintf( fpOut, "{\n" );
	fprintf( fpOut, "\tchar	*name = new char[strlen( MECHANISM ) + 5];\n" );
	fprintf( fpOut, "\tsprintf( name, \"/%%s.pre\", MECHANISM );\n" );
	fprintf( fpOut, "\tif ( strstr( mechName, name ) == 0 ) {\n" );
	fprintf( fpOut, "\t\tcerr << \"#error: program linked with mechanism \" << MECHANISM << \".mech\"\n" );
	fprintf( fpOut, "\t\t\t\t\t<< \", input mechanism is \" << mechName << NEWL;\n" );
	fprintf( fpOut, "\t\texit(2);\n" );
	fprintf( fpOut, "\t}\n" );
	fprintf( fpOut, "\telse {\n" );
	fprintf( fpOut, "\t\tcerr << \"use reduced kinetics \" << MECHANISM << NEWL;\n" );
	fprintf( fpOut, "\t}\n" );
	fprintf( fpOut, "}\n" );

	fclose( fpOut );
}

static void GlobalReactions( MatrixPtr nu, MatrixPtr *glob, MatrixPtr *rate )
{
	int			i, j, row, col;
	int			nForward = gForwardReacs->nForward;
	int			nReactants = gReactants->nReactants;
	int			nReducs = gReduction->nReducs;
	int			nReacs = nForward - nReducs;
	int			numAtoms = CountItems( gAtomsList );
	int			globalReacs = gGlobalReacs->nGlobalReacs;
	int			globalSpecies = gGlobalReacs->nGlobalSpecies;
	ReactionPtr	reac = NULL;
	SpeciesPtr	spec = NULL;
	Double		**mat = nu->mat;
	char		roman[16];
	MatrixPtr	red = NewMatrix( globalSpecies, nReacs, kRowPointers );
	MatrixPtr	globInv = NULL, globT = NULL, globTglob = NULL, globTred = NULL;
	Double		val = 0;
	FILE		*fp = stdout;
	char		*bptr = NULL;
	char		tempbuff[128];
	
	(*glob) = NewMatrix( globalSpecies, globalReacs, kRowPointers );

	/*	assemble glob and red, where
		
		¥	glob is the matrix of stoichiometric coefficients for 
			the reduced mechanism written with global reaction rates,
			
		¥	red is the matrix of stoichiometric coefficients for 
			the reduced mechanism written with the original reaction rates.
	*/
	for ( i = 0; i < globalSpecies; ++i ) {
		row = gGlobalReacs->spec[i];

		for ( j = 0; j < globalReacs; ++j ) {
			col = gGlobalReacs->reac[j];
			(*glob)->mat[i][j] = mat[row][col];
		}
		
		j = 0;
		for ( col = 0; col < nForward; ++col ) {
			if ( ReactionIsEliminated( gForwardReacs->reac[col], gReduction ) ) continue;
			
			red->mat[i][j++] = mat[row][col];
		}
	}


	if ( IndependentEQ( (*glob) ) != globalReacs ) 
		FatalError( "# global reaction are linearly dependent" );
	
	/*	to compute the global reaction rates, we have to solve 
		the system of linear equations 
			
			L[X] = red . (w_1, w_2,..., w_n)T = glob . (w_I, w_II, w_III,...)
			
		since the matrix glob has more rows than columns, we multiply 
		the equation with glob's transpose, globT. Now we can compute
		the inverse of globT . glob and take the product of the equation
		with (globT . glob)^-1, giving the result
		
			(w_I, w_II, w_III,...) = rate . (w_1, w_2,...,w_n).
		
		where rate denotes ( (globT . glob)^-1 . globT . red ).
		
	*/
	globT 		= Transpose( (*glob) );
	globTglob 	= MatrixMult( globT, (*glob) );
	globInv 	= Inverse( globTglob );
	globTred 	= MatrixMult( globT, red );
	(*rate) 		= MatrixMult( globInv, globTred );
	
	DisposeMatrix( globT );
	DisposeMatrix( globTglob );
	DisposeMatrix( globInv );
	DisposeMatrix( globTred );

#	ifdef DEBUG
	{
		AMPrintOptions	prnt;
		
		DefaultAMPOpts( &prnt );
		prnt.title = "matrix of stoichiometric coefficients of global reactions";
		prnt.sep = "\t";
		prnt.rowLabel = NULL;
		PrintMatrix( (*glob), &prnt, stdout );
		fputs( "\n", stdout );
		prnt.title = "matrix of stoichiometric coefficients of reduced reactions";
		prnt.sep = "\t";
		prnt.rowLabel = NULL;
		PrintMatrix( red, &prnt, stdout );
		fputs( "\n", stdout );
		prnt.title = "rate matrix";
		prnt.sep = "\t";
		prnt.rowLabel = NULL;
		PrintMatrix( (*rate), &prnt, stdout );
		fputs( "\n", stdout );
	}
#	endif
	
	bptr = gBuffer;
	
	/*	Print the global reaction rates as linear combination of the 
		original reaction rates.
	*/
	fputs( "\n# global reaction rates\n", fp );
	for ( i = 0; i < globalReacs; ++i ) {		
		RomanNumeral( i+1, roman, TRUE );
		sprintf( tempbuff, "w%s =", roman );
		strcat( gBuffer, tempbuff );
		
		col = 0;
		for ( j = 0; j < nForward; ++j ) {
			ReactionPtr		r = gForwardReacs->reac[j];
			
			if ( ReactionIsEliminated( r, gReduction ) ) continue;
			
			if ( ( fabs( val = (*rate)->mat[i][col++] ) ) > kTiny ) {
				if ( bptr - gBuffer > kLineLength ) {
					fprintf( fp, "%s\n\t", gBuffer );
					bptr = gBuffer;
				}
				
				sprintf( tempbuff, " %+g ", val );
				strcat( gBuffer, tempbuff );
				if ( r->backward ) {
					sprintf( tempbuff, "(w%s - w%s)", r->label, r->backward->label );
					strcat( gBuffer, tempbuff );
				}
				else {
					sprintf( tempbuff, "w%s", r->label );
					strcat( gBuffer, tempbuff );
				}
			}
		}
		fprintf( fp, "%s\n", gBuffer );
		bptr = gBuffer;
	}
	
	/*	Print the conservation equations for the global species.
	*/
	bptr = gBuffer;
	fputs( "\n# global reactions\n", fp );
	for ( i = 0; i < globalSpecies; ++i ) {
		sprintf( tempbuff, "L[%s] =", gReactants->spec[ gGlobalReacs->spec[i] ]->name );
		strcat( gBuffer, tempbuff );
		for ( j = 0; j < globalReacs; ++j ) {
			if ( fabs( val = (*glob)->mat[i][j] ) > kTiny ) {
				if ( bptr - gBuffer > kLineLength ) {
					fprintf( fp, "%s\n\t", gBuffer );
					bptr = gBuffer;
				}
				
				RomanNumeral( j+1, roman, TRUE );
				sprintf( tempbuff, " %+g ", val );
				strcat( gBuffer, tempbuff );
				sprintf( tempbuff, "w%s", roman );
				strcat( gBuffer, tempbuff );
			}
		}
		fprintf( fp, "%s\n", gBuffer );
		bptr = gBuffer;
	}
	
	DisposeMatrix( red );
	
	fflush( fp );
}


static void MathSpeciesID( void *ptr, void *aux )
{
	SpeciesPtr	sPtr = (SpeciesPtr)ptr;
	FILE		*fp = aux;
	
	MathematicaSymbol( gBuffer, sPtr->name );
	UpperString( gBuffer );
	fprintf( fp, "s%s = %d;\n", gBuffer, sPtr->number + 1 );
}


static void MathThirdBodyID( void *ptr, void *aux )
{
	ThirdBodyPtr	tPtr = (ThirdBodyPtr)ptr;
	int				numSpecies = CountItems( gSpeciesList );
	FILE			*fp = aux;
	
	MathematicaSymbol( gBuffer, tPtr->name );
	UpperString( gBuffer );
	fprintf( fp, "s%s = %d;\n", gBuffer, tPtr->id + 1 + numSpecies );
}


static void MathMolarWeight( void *ptr, void *aux )
{
	SpeciesPtr	sPtr = (SpeciesPtr)ptr;
	FILE		*fp = aux;
	
	MathematicaSymbol( gBuffer, sPtr->name );	
	UpperString( gBuffer );
	fprintf( fp, "mw[[s%s]] = %.3f;\n", gBuffer, sPtr->molarMass );
}


static void MathReactionLabels( Pointer ptr, Pointer aux )
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


static void MathReactionRates( Pointer ptr, Pointer aux )
{
	ReactionPtr		reac = (ReactionPtr)ptr;
	const int		kBufferSize = 256;
	int				i, n = reac->numberOfSpecies, steadyStateCount = 0;
	int				*spec = reac->speciesNumber;
	Double			*coeff = reac->speciesCoeff;
	SpeciesPtr		sPtr = NULL;
	FILE			*fp = (FILE *)aux;
	char			*bptr = NULL;
	char			tempbuff[128];
	
	bptr = gBuffer;
	
	sprintf( tempbuff, "\tw[r%s] -> rk[r%s]", reac->label, reac->label );
	strcat( gBuffer, tempbuff );
		
	for ( i = 0; i < n; ++i ) {
		if ( coeff[i] < 0.0 ) continue;		/* skip products */
		
		if ( !( sPtr = FindListItem( gSpeciesList, &spec[i], PTFindSpecies ) ) )
			FatalError( "species not found" );
		
		if ( SpeciesIsSteadyState( sPtr /*, gReduction*/ ) ) ++steadyStateCount;
			
		if ( fabs(coeff[i] - 1.0) < kTiny ) {
			sprintf( tempbuff, " c[s%s]", sPtr->name );
			strcat( gBuffer, tempbuff );
		}
		else {
			sprintf( tempbuff, " c[s%s]^%g", sPtr->name, coeff[i] );
			strcat( gBuffer, tempbuff );
		}
	}
	
	if ( reac->withThirdBody ) {
		ThirdBodyPtr	tb = NULL;
		int				tbID = reac->thirdBodyNumber;
		if ( !( tb = FindListItem( gThirdBodyList, &tbID, 
			PTFindThirdBody ) ) ) FatalError( "third body not found" );
		sprintf( tempbuff, " c[s%s]", tb->name );
		strcat( gBuffer, tempbuff );
	}
	
	/*	only if there is at least one steady-state species in the reaction,
		we expand the reaction rate in terms of concentrations and rate 
		coefficient.
	*/
	if ( steadyStateCount ) {
		fprintf( fp, "%s,\n", gBuffer );
	}
	else {
		fprintf( fp, "(* %s, *)\n", gBuffer );
	}
}


static void MathPrintSteadyStates( MatrixPtr nu, FILE *fp )
{
	int			nForward = gForwardReacs->nForward;
	int			numReduced = nForward - gReduction->nReducs;
	int			numSteady = gReduction->nReducs;
	int			i, j, count = 0;
	const int	*reac = gReduction->reac, *spec = gReduction->spec;
	Double		**mat = nu->mat;
	char		*bptr = gBuffer;
	char		tempbuff[128];
	
	fputs( "\n(* steady state relations *)\n", fp );
	fprintf( fp, "numSteadyStates = %d;\n", numSteady );
	fputs( "eq = Table[ 0, {numSteadyStates} ];\n", fp );
	
	for ( i = 0; i < numSteady; ++i ) {
		int		s = spec[i];
		int		first = TRUE;
		Double	val = 0.0;
		
		sprintf( tempbuff, "(* L[s%s] == 0 *) eq[[%d]] = ",
			gReactants->spec[s]->name, i + 1  );
		strcat( gBuffer, tempbuff );
		for ( j = 0; j < nForward; ++j ) {
			ReactionPtr		r = gForwardReacs->reac[j];

			if ( fabs( val = mat[s][j] ) > kTiny  ) {
				if ( bptr - gBuffer > kLineLength ) {
					fprintf( fp, "%s\n\t", gBuffer );
					bptr = gBuffer;
				}
				
				if ( first ) {
					first = FALSE;
					if ( fabs( val ) - 1.0 < kTiny ) {
						sprintf( tempbuff, (val > 0.0) ? "" : " -" );
						strcat( gBuffer, tempbuff );
					}
					else if ( val > 0.0 ) {
						sprintf( tempbuff, " %g", val );
						strcat( gBuffer, tempbuff );
					}
					else {
						sprintf( tempbuff, " -%g", -val );
						strcat( gBuffer, tempbuff );
					}
				}
				else {
					if ( fabs( val ) - 1.0 < kTiny ) {
						sprintf( tempbuff, (val > 0.0) ? " +" : " -" );
						strcat( gBuffer, tempbuff );
					}
					else if ( val > 0.0 ) {
						sprintf( tempbuff, " + %g", val );
						strcat( gBuffer, tempbuff );
					}
					else {
						sprintf( tempbuff, " - %g", -val );
						strcat( gBuffer, tempbuff );
					}
				}
				if ( r->backward ) {
					sprintf( tempbuff, " (w[r%s] - w[r%s])", r->label, r->backward->label );
					strcat( gBuffer, tempbuff );
				}
				else {
					sprintf( tempbuff, " w[%s]", r->label );
					strcat( gBuffer, tempbuff );
				}
			}
		}
		fprintf( fp, "%s;\n", gBuffer );
		bptr = gBuffer;
	}
	
}


static void MathPrintGlobalRates( MatrixPtr rate, FILE *fp )
{
	int			i, j, col;
	int			nForward = gForwardReacs->nForward;
	int			globalReacs = gGlobalReacs->nGlobalReacs;
	int			globalSpecies = gGlobalReacs->nGlobalSpecies;
	Double		val;
	Double		**mat = rate->mat;
	char		*bptr = gBuffer;
	char			tempbuff[128];

	fputs( "\n(* global reaction rates *)\n", fp );
	fprintf( fp, "numGlobalReactions = %d;\n", globalReacs );
	fputs( "wGlobal = Table[ 0, {numGlobalReactions} ];\n", fp );
	
	for ( i = 0; i < globalReacs; ++i ) {		
		int		first = TRUE;

		sprintf( tempbuff, "wGlobal[[%d]] =", i+1 );
					strcat( gBuffer, tempbuff );
		
		col = 0;
		for ( j = 0; j < nForward; ++j ) {
			ReactionPtr		r = gForwardReacs->reac[j];
			
			if ( ReactionIsEliminated( r, gReduction ) ) continue;
			
			if ( ( fabs( val = mat[i][col++] ) ) > kTiny ) {
				if ( bptr - gBuffer > kLineLength ) {
					fprintf( fp, "%s\n\t", gBuffer );
					bptr = gBuffer;
				}

				if ( first ) {
					first = FALSE;
					if ( fabs( val ) - 1.0 < kTiny ) {
						sprintf( tempbuff, (val > 0.0) ? "" : " -" );
						strcat( gBuffer, tempbuff );
					}
					else if ( val > 0.0 ) {
						sprintf( tempbuff, " %g", val );
						strcat( gBuffer, tempbuff );
				}
					else {
						sprintf( tempbuff, " -%g", -val );					
						strcat( gBuffer, tempbuff );
					}
				}
				else {
					if ( fabs( val ) - 1.0 < kTiny ) {
						sprintf( tempbuff, (val > 0.0) ? " +" : " -" );
						strcat( gBuffer, tempbuff );
					}
					else if ( val > 0.0 ) {
						sprintf( tempbuff, " + %g", val );
						strcat( gBuffer, tempbuff );
					}
					else {
						sprintf( tempbuff, " - %g", -val );
						strcat( gBuffer, tempbuff );
					}
				}
				if ( r->backward ) {
					sprintf( tempbuff, " (w[r%s] - w[r%s])", r->label, r->backward->label );
					strcat( gBuffer, tempbuff );
				}
				else {
					sprintf( tempbuff, " w[r%s]", r->label );
					strcat( gBuffer, tempbuff );
				}
			}
		}
		fprintf( fp, "%s;\n", gBuffer );
		bptr = gBuffer;
	}
}
 

static void MathPrintGlobalReacs( MatrixPtr glob, FILE *fp )
{
	int			i, j;
	int			globalReacs = gGlobalReacs->nGlobalReacs;
	int			globalSpecies = gGlobalReacs->nGlobalSpecies;
	Double		val;
	Double		**mat = glob->mat;
	char		*bptr = gBuffer;
	char		tempbuff[128];

	fputs( "\n(* global species conservation equations *)\n", fp );
	fputs( "L = Table[ 0, {numSpecies} ];\n", fp );
	
	for ( i = 0; i < globalSpecies; ++i ) {
		int		first = TRUE;

		sprintf( tempbuff, "L[[s%s]] =", 
			gReactants->spec[ gGlobalReacs->spec[i] ]->name );
		strcat( gBuffer, tempbuff );
		
		for ( j = 0; j < globalReacs; ++j ) {
			if ( ( fabs( val = mat[i][j] ) ) > kTiny ) {
				if ( bptr - gBuffer > kLineLength ) {
					fprintf( fp, "%s\n\t", gBuffer );
					bptr = gBuffer;
				}

				if ( first ) {
					first = FALSE;
					if ( fabs( val ) - 1.0 < kTiny ) {
						sprintf( tempbuff, (val > 0.0) ? "" : " -" );
						strcat( gBuffer, tempbuff );
					}
					else if ( val > 0.0 ) {
						sprintf( tempbuff, " %g", val );
						strcat( gBuffer, tempbuff );
					}
					else {
						sprintf( tempbuff, " -%g", -val );					
						strcat( gBuffer, tempbuff );
					}
				}
				else {
					if ( fabs( val ) - 1.0 < kTiny ) {
						sprintf( tempbuff, (val > 0.0) ? " +" : " -" );
						strcat( gBuffer, tempbuff );
					}
					else if ( val > 0.0 ) {
						sprintf( tempbuff, " + %g", val );
						strcat( gBuffer, tempbuff );
					}
					else {
						sprintf( tempbuff, " - %g", -val );
						strcat( gBuffer, tempbuff );
					}
				}
				sprintf( tempbuff, " wGlobal[[%d]]", j+1 );
				strcat( gBuffer, tempbuff );
			}
		}
		fprintf( fp, "%s;\n", gBuffer );
		bptr = gBuffer;
	}
}


void MathPrintSSConc( FILE *fp )
{
	int			nReducs = gReduction->nReducs;
	int			i;
	
	fputs( "\n(* Steady state concentrations *)\n", fp );
	fputs( "steadyStateConcentrations = {\n", fp );
	for ( i = 0; i < nReducs; ++i ) {
		SpeciesPtr	spec = gReactants->spec[ gReduction->spec[i] ];
		
		MathematicaSymbol( gBuffer, spec->name );	
		UpperString( gBuffer );
		fprintf( fp, "\tc[s%s]%s", gBuffer, (i == nReducs - 1) ? "\n": ",\n" );
	}
	fputs( "};\n", fp );
}


static void ReduxMathPrint( const char *fileName, MatrixPtr nu, MatrixPtr glob,
	MatrixPtr rate )
{
	int			numSpecies = CountItems( gSpeciesList );
	int			numReactions = CountItems( gReactionList );
	int			numThirdBodies = CountItems( gThirdBodyList );
	FILE		*fp = NULL;
	
	if ( !( fp = fopen( fileName, "w" ) ) )
		FatalError( "couldn't open mathematica file" );
	
	/*	Species and third body information.
	*/
	fputs( "\n(* Species and third body information *)\n", fp );
	fprintf( fp, "numSpecies = %d;\n", numSpecies );
	Iterate( gSpeciesList, MathSpeciesID, fp );
	fprintf( fp, "\nnumThirdBodies = %d;\n", numThirdBodies );
	Iterate( gThirdBodyList, MathThirdBodyID, fp );
	
	/*	Molar weights.
	*/
	fputs( "\nmw = Table[ 0, {numSpecies} ];\n", fp );
	Iterate( gSpeciesList, MathMolarWeight, fp );
	
	/*	Reaction information.
	*/
	fputs( "\n(* Reaction information *)\n", fp );
	fprintf( fp, "numReactions = %d;\n", numReactions );
	Iterate( gReactionList, MathReactionLabels, fp );
	
	/*	Reaction rates.
	*/
	fputs( "\n(* Expansion rules for the reaction rates *)\n", fp );
	fprintf( fp, "numReactions = %d;\n", numReactions );
	fputs( "Array[ w, numReactions ];\n", fp );
	fputs( "Array[ rk, numReactions ];\n", fp );
	fputs( "Array[ c, numSpecies + numThirdBodies ];\n", fp );
	fputs( "wExpand = {\n", fp );
	Iterate( gReactionList, MathReactionRates, fp );
	fputs( "\t1 -> 1\n};\n", fp );			/* a patch */
	
	/*	Steady state relations
	*/
	MathPrintSteadyStates( nu, fp );
	
	/*	Global reactions	
	*/
	MathPrintGlobalRates( rate, fp );
	MathPrintGlobalReacs( glob, fp );
	
	/*	Steady state concentrations
	*/
	MathPrintSSConc( fp );
	
	fclose( fp );
	
	fprintf( stderr, "# reduced mechanism written to \"%s\".\n", fileName );
#	ifdef applec
	fsetfileinfo( fileName, 'OMEG', 'TEXT' );
#	endif
}


void ReduxCleanup( void )
{	
	free( gBuffer );
	FreeForwardReaction( gForwardReacs );
	FreeReactants( gReactants );
	FreeReduction( gReduction );
	FreeGlobalReaction( gGlobalReacs );
}


void ReduceMechanism( void )
{
	MatrixPtr	nu = NULL;
	MatrixPtr	glob = NULL;
	MatrixPtr	rate = NULL;

/*	char		mathFile[256];*/
	
	if ( !( gBuffer = (char *)malloc( 5 * kLineLength * sizeof(*gBuffer) ) )  )
		FatalError( kOutOfMemory );
	memset( gBuffer, 0, 5 * kLineLength * sizeof(*gBuffer) );
	gFloatBuffer = gBuffer + 4 * kLineLength;

	FixForwardBackwardsPtrs();

	nu = StoichCoeffMatrix();

	ReadSteadyStates( gOptions->steadyStates, nu );
	
	SteadyStateRelations( nu );

#ifndef NEWREDUX

	Reduce( nu );
#endif

/*	GlobalReactions( nu, &glob, &rate );*/

/*	strcpy( mathFile, gOptions->steadyStates );*/
/*	strcat( mathFile, ".m" );*/
/*	ReduxMathPrint( mathFile, nu, glob, rate );*/
	
	if ( glob ) DisposeMatrix( glob );
	if ( rate ) DisposeMatrix( rate );
	if ( nu ) DisposeMatrix( nu );
	
	ReduxCleanup();
}


/*	void RedMech( SteadyStateInfoPtr ss )

	Same as ReduceMechanism() except for printing routines and
	clean up functions.
*/
void RedMech( SteadyStateInfoPtr ss )
{	
	if ( !( gBuffer = (char *)malloc( 5 * kLineLength * sizeof(*gBuffer) ) )  )
		FatalError( kOutOfMemory );
	memset( gBuffer, 0, 5 * kLineLength * sizeof(*gBuffer) );
	gFloatBuffer = gBuffer + 4 * kLineLength;

	FixForwardBackwardsPtrs();

	/*	ss->nu contains only forward reactions and reacting species, i.e
		no N2, Ar, etc.														*/
	ss->nu = StoichCoeffMatrix();
	ReadSteadyStates( gOptions->steadyStates, ss->nu );
	Reduce( ss->nu );
	
	/*SteadyStateRelations( ss->nu );*/
	/*GlobalReactions( ss->nu, &ss->glob, &ss->rate );*/
}


ReactionPtr GetForwardReaction( int k )
{
	/*	ReactionPtr GetForwardReaction( int k );
	
		Given the column index (== reaction index) to a matrix of
		stoichiometric coefficients, this function returns a pointer
		to the corresponding Reaction record.
	*/
	return *(gForwardReacs->reac + k);
}
