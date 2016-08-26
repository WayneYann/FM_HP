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

#include "FitPak.h"

#define NVAR
#define INSERTBACKWARDS
#undef DEBUGEQCONST
#undef NEWEQUICONST

static void FitBasisFunction( Double *x, Double p[], int np );
static Double ExactKBack( Double T, ReactionPtr forw, Double aKonv, Double eKonv ); 
static Double FreeEnthalpy( Double T, SpeciesPtr species );
static Double Enthalpy( Double T, SpeciesPtr species );
static void NewFitBasisFunction( Double *x, Double p[], int np );

void CheckReactions( void )
{
	char	*error = NULL;
	LinkPtr link = 0;
		
	if ( link = ListIterator( gReactionList, CheckReactionZeroEntry, NULL ) ) {
		fprintf( stderr, "# warning: error message printed to stdout\n\n" );
		fprintf( stdout, "# warning: species appears on both sides of reaction\n\n" );
		PrintReaction( link, gReactionList );
		/*exit( 2 );*/
	}

	if ( !( AllOnesInTwo( gUsedThirdBodyList, gThirdBodyList, FindThirdBody, &error ) ) ) {
		fprintf( stderr, "# error: there is no definition for the coefficients of\n#        the thirdbody efficiency '%s'\n", error );
		exit( 2 );
	}

	ListIterator( gReactionList, SetOrderOfReaction, NULL );

	if ( gOptions->useForwardBackward ) {
		if ( CheckForwardBackward() ) {
			exit( 2 );
		}
	}
	
	if ( !gHeader->globalMechanism && CheckForDoubles() ) {
/*		exit( 2 );*/
	}

	if ( gOptions->useForwardBackward ) ComputeLindemann();

	if ( CheckBroadening() ) {
		exit( 2 );
	}

}

/*int CheckSpecies( void )
{
	
}
*/
int CheckForDoubles( void )
{
	LinkPtr 	link1 = NULL;
	LinkPtr 	link2 = NULL;
	ReactionPtr	reaction1 = NULL;
	ReactionPtr	reaction2 = NULL;
	int			err = FALSE;
	int 		i;
			
	for ( i = 0; i < NumberOfItems( gReactionList ); ++i ) {
		link1 = Find_n_th_Item( i+1, gReactionList );
		reaction1 = ( ReactionPtr )link1->item;
		if ( link2 = ListIterator( gReactionList, CompareDoubleReaction, reaction1 ) ) {
			reaction2 = ( ReactionPtr )link2->item;
			if ( strpbrk( reaction1->label, "." ) != NULL 
				&& strpbrk( reaction2->label, "." ) != NULL
				&& strcspn( reaction1->label, "." ) == strcspn( reaction2->label, "." )
				&& strncmp( reaction1->label, reaction2->label, strcspn( reaction1->label, "." ) ) == 0
				 && strcmp( reaction1->label, reaction2->label ) != 0 ) {
/*  				fprintf( stderr, "reaction %s and %s are equal\n", reaction1->label*/
/*						, reaction2->label );*/
			}
			else {
				fprintf( stderr, "# warning: error message printed to stdout\n\n" );
				fprintf( stdout, "# warning: doubly defined reaction\n\n" );
				PrintReaction( link1, gReactionList );
				fprintf( stdout, "# and\n\n" );
				PrintReaction( link2, gReactionList );
				reaction1->isDouble = TRUE;
				reaction2->isDouble = TRUE;
				if ( reaction1->hasForward ) {
				  reaction1->forward->isDouble = TRUE;
				}
				else if ( reaction1->hasBackward ) {
				  reaction1->backward->isDouble = TRUE;
				}

				if ( reaction2->hasForward ) {
				  reaction2->forward->isDouble = TRUE;
				}
				else if ( reaction2->hasBackward ) {
				  reaction2->backward->isDouble = TRUE;
				}
				err = TRUE;
			}
		}
	}
	return err;
}

int CheckStoichiometry( ReactionPtr reaction )
{
	LinkPtr				speciesLink = NULL;
	SpeciesPtr			species = NULL;
	static VectorPtr	sum = NULL;
	int i, j;
	
	if ( !sum ) {
		sum = NewVector( gAllowAtomList->items );
	}
	
	for ( i = 0; i < sum->len; ++i ) {
		sum->vec[i] = 0.0;
	}
	
	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[i] );
		species = ( SpeciesPtr )speciesLink->item;
		for ( j = 0; j < sum->len; ++j ) {
			sum->vec[j] += reaction->speciesCoeff[i] * species->composition->vec[j];
		}
	}
	
	for ( i = 0; i < sum->len; ++i ) {
		if ( sum->vec[i] != 0 ) {
			return 1;
		}
	}	
	
	return 0;
}

int CheckForwardBackward( void )
{
	LinkPtr 	link1 = NULL;
	LinkPtr 	link2 = NULL;
	ReactionPtr	reaction1 = NULL;
	ReactionPtr	reaction2 = NULL;
	ReactionPtr	back = NULL;
	int			err = FALSE;
	int 		i, j, appendTo;
	int			len;
	char		label[32];
	ListPtr		errorList = NULL;
	Flag        isForw;
	
	if ( !errorList ) {
		errorList = NewList( "errorList" );
	}
			
	for ( i = 0; i < NumberOfItems( gReactionList ); ++i ) {
		link1 = Find_n_th_Item( i+1, gReactionList );
		reaction1 = ( ReactionPtr )link1->item;
		
		strcpy( label, reaction1->label );
		len = strlen( label );
		
		if ( ListIterator( errorList, FindString, label ) ) {
			continue;
		}
		
		if ( label[len-1] == 'f' ) {
			label[len-1] = 'b';
			isForw = TRUE;
		}
		else if ( label[len-1] == 'b' ) {
			label[len-1] = 'f';
			isForw = FALSE;
		}
		else {
			continue;
		}

		if ( link2 = ListIterator( gReactionList, FindReaction, label ) ) {
			reaction2 = ( ReactionPtr )link2->item;
			/* turn signs of coefficients of reaction1 and use the function CompareReaction */
			for ( j = 0; j < reaction1->numberOfSpecies; ++j ) {
				reaction1->speciesCoeff[j] *= -1.0;
			}
			
			if ( !( CompareReaction( reaction1, reaction2 ) ) ) {
				fprintf( stderr, "# error: \n\n " );
				PrintReaction( link1, gReactionList );
				fprintf( stderr, "# and\n\n " );
				PrintReaction( link2, gReactionList );
				fprintf( stderr, "# are not corresponding forward and backward reactions\n\n " );
				err = TRUE;
				AddString( errorList, label );
			}
			else {
			  if ( isForw == TRUE ) {
				reaction1->backward = reaction2;
				reaction2->forward = reaction1;
			  }
			  else {
				reaction2->backward = reaction1;
				reaction1->forward = reaction2;
			  }
			}
			/* second turn of signs */
			for ( j = 0; j < reaction1->numberOfSpecies; ++j ) {
				reaction1->speciesCoeff[j] *= -1.0;
			}
			

		}
		else {
			appendTo = 0;
#ifdef INSERTBACKWARDS
			back = InsertReaction( gReactionList, label, i );
#else
			back = AddReaction( gReactionList, label );
#endif
			back->numberOfSpecies = reaction1->numberOfSpecies;
			for ( j = 0; j < back->numberOfSpecies; ++j ) {
				back->speciesCoeff[j] = - reaction1->speciesCoeff[reaction1->numberOfSpecies-j-1];
				back->speciesNumber[j] = reaction1->speciesNumber[reaction1->numberOfSpecies-j-1];
			}
			back->withLindemann = reaction1->withLindemann;
			back->lindemannNumber = reaction1->lindemannNumber;
			back->fca = reaction1->fca;
			back->fcb = reaction1->fcb;
			back->fcc = reaction1->fcc;
			back->fcTa = reaction1->fcTa;
			back->fcTb = reaction1->fcTb;
			back->fcTc = reaction1->fcTc;
			back->withThirdBody = reaction1->withThirdBody;
			back->thirdBodyNumber = reaction1->thirdBodyNumber;
			back->partialEquilibrium = reaction1->partialEquilibrium;
/* changed by hp 04.08.97 */
/*			back->backward = reaction1;*/
/*			reaction1->forward = back;*/
			if ( isForw == TRUE ) {
			  back->forward = reaction1;
			  reaction1->backward = back;
			}
			else {
			  back->backward = reaction1;
			  reaction1->forward = back;
			}
			ComputeRateConstants( reaction1, back );
			SetOneOrderOfReaction( back );
			back->rateFromEquilibrium = TRUE;
		}
	}
		
	return err;
}

int GetOrderOfReaction( ReactionPtr reac, Flag infRequired )
{
	int		ordReac = ( reac->withThirdBody ) ? 1: 0;
	int		i;
	
	for ( i = 0; i < reac->numberOfSpecies; ++i ) {
		if ( reac->speciesCoeff[i] > 0.0 ) {
			ordReac += ( int )reac->speciesCoeff[i];
		}
	}
	
	if ( reac->withLindemann ) {
		if ( reac->withThirdBody && infRequired ) {
			ordReac -= 1;
		}
		else if ( !reac->withThirdBody && !infRequired ) {
			ordReac += 1;
		}
	}
	
	return ordReac;
}

Double GetEConverter( ReactionPtr reac, Flag infRequired )
{
/*  compute eKonverter */

	int				orderOfReaction;
	LinkPtr			link = NULL;
	DimensionPtr	dim;
	Double			EKonv = 0.0;

	orderOfReaction = GetOrderOfReaction( reac, infRequired );

	if ( !( link = ListIterator( gUsedDimensionList, FindDimension, "E" ) ) ) {
		FatalError( "no units for E specified" );
	}
	dim = link->item;
	EKonv = dim->value;
	if ( dim->orderOfReaction ) {
		EKonv *= pow( dim->orderOfReacValue, ( Double )orderOfReaction );
	}
	if ( dim->tempExponent ) {
		EKonv *= pow( dim->tempExpValue, reac->n );
	}

	return EKonv;
}

Double GetAConverter( ReactionPtr reac, Flag infRequired )
{
/*  compute AKonverter for forward reaction  */

	int				orderOfReaction;
	LinkPtr			link = NULL;
	DimensionPtr	dim;
	Double			AKonv = 0.0;

	orderOfReaction = GetOrderOfReaction( reac, infRequired );

	if ( !( link = ListIterator( gUsedDimensionList, FindDimension, "A" ) ) ) {
		FatalError( "no units for A specified" );
	}
	dim = (DimensionPtr) link->item;
	AKonv = dim->value;
	if ( dim->orderOfReaction ) {
		AKonv *= pow( dim->orderOfReacValue, ( Double )orderOfReaction );
	}
	if ( dim->tempExponent ) {
		AKonv *= pow( dim->tempExpValue, reac->n );
	}

	return AKonv;
}

Double GetETokmoleConverter( ReactionPtr reac, Flag infRequired )
{
/*  compute converter conv: E [J/kmole] = E [J/mole] * conv  */

	const Double	moleM1TokmoleM1 = 1000.0;

	return moleM1TokmoleM1;
}

Double GetATokmoleConverter( ReactionPtr reac, Flag infRequired )
{
/*  compute converter conv: A [...,kmole] = A [...,mole] * conv  */

	const Double	moleM1TokmoleM1 = 1000.0;
	int				orderOfReaction;

	orderOfReaction = GetOrderOfReaction( reac, infRequired );

	return pow( moleM1TokmoleM1, ( Double )orderOfReaction-1.0 );
}

void ComputeRateConstants( ReactionPtr forw, ReactionPtr backw )
{
/* 	
 *	computes A and E from y = ln k_b - n * ln T = ln A - E/RT, 
 *	where A = exp( model->a[1] ) and E = a[2] * R
 */
/*#ifdef DEBUGEQCONST*/
/*	static int		init = 0;*/
/*	Double			*Kc = New1DArray( gReactionList->items );*/
/*	Double			*KcFit = New1DArray( gSpeciesList->items );*/
/*#endif*/
#ifdef NVAR
	const int		nConsts = 3;
#else
	const int		nConsts = 2;
#endif
	const int		nPoints = 100;
	const Double	TRef = 298.16;
	int				i;
	Double			dx, xStart = 500.0, xEnd = 2000.0; /*GB before xEnd=3000K*/
	DataSetPtr		data_set = NULL;
	ModelPtr		model = NULL;
	DataPointPtr	ptr = NULL;
	LinkPtr			link = NULL;
	Double			eKonverter;
	Double			aKonverter;
	Double			aBackKonv;
	int				ordBackReac;
	SpeciesPtr		species = NULL;
	Double			*pi_A = New1DArray( gSpeciesList->items );
	Double			*pi_B = New1DArray( gSpeciesList->items );
#ifdef NEWEQUICONST
	DimensionPtr	dimE;
	Double			AK, nK, EK;
	int				j;
#endif
	const Double	moleM1TokmoleM1 = 1000.0;
#ifdef DEBUGEQCONST
	char			fName[128];
	FILE			*fp;
	Double			kbEqui, kb, kbFit, locA, locE;
	sprintf( fName, "KC_%s.dout", forw->label );
	fp = fopen( fName, "w" );
	fprintf( fp , "*\n1000/T\tkb\tkbExact\tkbfit\n" );	
#endif
	data_set = NewDataSet( nPoints, 1 );
	model = NewModel( nConsts, kLinear );

#ifdef applec
	SpinCursor( -32 ); 
#endif

/*  compute Konverter */
	eKonverter = GetEConverter( forw, FALSE ) * GetETokmoleConverter( forw, FALSE );
	aKonverter = GetAConverter( forw, FALSE ) * GetATokmoleConverter( forw, FALSE );
	
/*  compute orderOfReaction for backward reaction */
	ordBackReac = GetOrderOfReaction( backw, FALSE );

/*  compute AKonverter for backward reaction  */
	backw->n = forw->n;
	aBackKonv = GetAConverter( backw, FALSE ) * GetATokmoleConverter( backw, FALSE );
	
#ifdef NEWEQUICONST
/* first fit pi_i */
	for ( j = 0; j < gSpeciesList->items; ++j ) {
		link = Find_n_th_Item( j+1, gSpeciesList );
		species = ( SpeciesPtr ) link->item;
		/* compute pi_j(T) store in ptr->y */
		/*  equidistant 1/T  */
		dx = ( 1.0 / xEnd - 1.0 / xStart ) / ( nPoints - 1.0 );
		for ( i = 0; i <  nPoints; ++i ) {
			ptr = &data_set->data[i];
			ptr->x[1] = 1.0 / ( 1.0 / xStart + dx * i );
			ptr->y = ( Enthalpy( TRef, species ) 
						- FreeEnthalpy( ptr->x[1], species ) )
					/ ( RGAS * ptr->x[1] );
			ptr->sigma = 1.0;
			if ( strcmp( species->name, "H" ) == NULL ) {
				fprintf( stderr, "%g\t%g\n", ptr->x[1], ptr->y );
			}
		}
		for (i = 1; i <= nConsts; ++i) model->a[i] = 0.0;
		
		LinearFit( data_set, model, NewFitBasisFunction, NULL );
		fprintf( stderr, "piA_%s = %g\tpiB_%s = %g\n"
				, species->name, model->a[1], species->name, model->a[2] );
		pi_A[j] = model->a[1];
		pi_B[j] = model->a[2];
	}
		exit(2);

	AK = nK = EK = 0.0;
	for ( i = 0; i < forw->numberOfSpecies; ++i ) {
		link = Find_n_th_Item( forw->speciesNumber[i]+1, gSpeciesList );
		species = ( SpeciesPtr ) link->item;
		AK += -forw->speciesCoeff[i] * pi_A[forw->speciesNumber[i]];
		nK += -forw->speciesCoeff[i] * pi_B[forw->speciesNumber[i]];
		EK += -forw->speciesCoeff[i] * Enthalpy( TRef, species );
	}
	AK = exp( AK );
	
	backw->a = forw->a / AK;
	backw->n = forw->n - nK;
	backw->e = forw->e - EK / eKonverter;
	if ( backw->e < -10.0 ) {
		fprintf( stderr, "#warning: activation energy of reaction %s is %f\n", backw->label, backw->e );
	}
#else
/*  equidistant 1/T  */
	dx = ( 1.0 / xEnd - 1.0 / xStart ) / ( nPoints - 1.0 );
	for ( i = 0; i <  nPoints; ++i ) {
		ptr = &data_set->data[i];
		ptr->x[1] = 1.0 / ( 1.0 / xStart + dx * i );
		ptr->y = ExactKBack( ptr->x[1], forw, aKonverter, eKonverter );
		ptr->sigma = 1.0;
	}
	
	for (i = 1; i <= nConsts; ++i) model->a[i] = 0.0;
	
	LinearFit( data_set, model, FitBasisFunction, NULL );
	
#ifdef NVAR
	backw->n = model->a[3];
#else
	backw->n = forw->n;
#endif
	backw->a = exp( model->a[1] ) / aBackKonv;
	backw->e = model->a[2] * RGAS / eKonverter;
	if ( backw->e < -10.0 ) {
		fprintf( stderr, "#warning: activation energy of reaction %s is %f\n", backw->label, backw->e );
	}
	
#ifdef DEBUGEQCONST
	for ( i = 0; i <  nPoints; ++i ) {
		ptr = &data_set->data[i];
		kbEqui = exp(ptr->y);
#ifndef NVAR
		kbEqui *= pow( ptr->x[1], forw->n );
#endif
		kbFit = backw->a * aBackKonv * pow( ptr->x[1], backw->n ) * exp( -backw->e * eKonverter / ( RGAS * ptr->x[1] ) );
		/*1b*/
		/*locA = 1.568e13;
		locE = 3.52;*/
		
		/*2b*/
		locA = 2.222e4;
		locE = 18.29;
		
		/*5b*/
		/*locA = 3.19e18;
		locE = 195.93;*/
		
		/*36b*/
		/*locA = 4.845e12;
		locE = 37.43;*/

		/*69f*/
		/*locA = 3.98e13;
		locE = 293.1;*/

		/*69b*/
		/*locA = 1.267e13;
		locE = 32.48;*/

		/*51b*/
		/*locA = 6.245e41;
		locE = 27.5;*/

		kb = locA * aBackKonv * pow( ptr->x[1], forw->n ) * exp( -locE * eKonverter/( RGAS * ptr->x[1] ) );
		fprintf( fp , "%g\t%g\t%g\t%g\n", 1000.0/ptr->x[1], kb, kbEqui, kbFit );	
	}
	fclose( fp );
#endif
#endif
	FreeModel( model );
	FreeDataSet( data_set );
	Free1DArray( pi_B );	
	Free1DArray( pi_A );	
}

static void NewFitBasisFunction( Double *x, Double p[], int np )
{
#ifdef applec
#pragma unused(np)
#endif
	Double	xx = x[1];

	p[1] = 1.0;
	p[2] = log( xx );
}

static void FitBasisFunction( Double *x, Double p[], int np )
{
#ifdef applec
#pragma unused(np)
#endif
	Double	xx = x[1];

	p[1] = 1.0;
	p[2] = -1.0 / xx;
#ifdef NVAR
	p[3] = log( xx );
#endif
}

static Double ExactKBack( Double T, ReactionPtr forw, Double aKonv, Double eKonv ) 
{
/*	this function returns ln k_b - n * ln( T ) = ln ( k_f / K_C ) - nf * ln( T ) */
	
	static const Double p0 = 1.0133e5;
	int			i;
	int			nSpecies = forw->numberOfSpecies;
	Double		sumNu = 0.0, sumNuMu = 0.0, lnKC;
	Double		RT = RGAS * T;
	Double		lnT = log( T );
	Double		lnKfMinNLnT;
	SpeciesPtr	species;
	
/*  compute lnKf -  n * ln T  */
	lnKfMinNLnT = log( forw->a * aKonv ) - forw->e * eKonv / RT;
#ifdef NVAR
	lnKfMinNLnT += forw->n * lnT;
#endif

/*  compute lnKc */
	for ( i = 0; i < nSpecies; ++i ) {
		sumNu += forw->speciesCoeff[i];
		species = ( SpeciesPtr )ListIterator( gSpeciesList, FindSpeciesNumber, &forw->speciesNumber[i] )->item;
		sumNuMu += forw->speciesCoeff[i] * FreeEnthalpy( T, species );
	}
	/* coefficients of products are defined positive */
	sumNu *= -1.0;
	
	sumNuMu *= -1.0;
	
	lnKC = -sumNu * log( RT / p0 ) - sumNuMu / RT;
	
	return lnKfMinNLnT - lnKC;
}

static Double Enthalpy( Double T, SpeciesPtr species )
{
	Double	*a = ( T > 1000.0 ) ? species->coeffHot : species->coeffCold;
	
	return RGAS * ( a[5] + T * ( a[0] + T * ( a[1] / 2.0 + T * ( a[2] / 3.0 
						+ T * ( a[3] / 4.0 + T * a[4] / 5.0 ) ) ) ) );
}

static Double FreeEnthalpy( Double T, SpeciesPtr species )
{
	Double	mu = 0.0;
	Double	*a = ( T > 1000.0 ) ? species->coeffHot : species->coeffCold;
	
	mu = a[0] * ( 1.0 - log( T ) ) + a[5] / T - a[6];
	mu -= 0.5 * T * ( a[1] + T * ( a[2] / 3.0 + T * ( a[3] / 6.0 + 0.1 * T * a[4] ) ) );
	mu *= RGAS * T;
	
	return mu;
}

/*int Checkbroadening( void )
{
	LinkPtr		nameLink = NULL;
	LinkPtr		contentsLink = NULL;
	LinkPtr		numberLink = NULL;
	LinkPtr		otherLink = NULL;
	String		name = NULL;
	String		contents = NULL;	
	String		otherContents = NULL;	
	String		checkBuffer = NULL;
	String		otherCheckBuffer = NULL;
	char		otherLabel[32];	
	int 		i;
	int			len;
	int 		number;
	Flag		err = FALSE;
	ListPtr		errorList = NULL;
	
	errorList = NewList( "errorList" );

	ListIterator( gFcNameList, PrintBroadening, gFcContentsList );	
	for ( i = 0; i < NumberOfItems( gFcNameList ); ++i ) {
		contentsLink = Find_n_th_Item( i+1, gFcContentsList );
		contents = ( String )contentsLink->item;
		nameLink = Find_n_th_Item( i+1, gFcNameList );
		name = ( String )nameLink->item;
		len = strlen( name );
		if ( name[len-1] == 'f' ||  name[len-1] == 'b' ) {
			strcpy( otherLabel, name );
			if ( name[len-1] == 'f' ) {
				otherLabel[len-1] = 'b';
			} 
			else {
				otherLabel[len-1] = 'f';
			}
		}
		if ( strcmp( contents, "" ) == 0 ) {
			if ( name[len-1] == 'f' ||  name[len-1] == 'b' ) {
				number = NumberOfListItem( gFcNameList, FindString, otherLabel );
				otherLink = Find_n_th_Item( number, gFcContentsList );
				otherContents = ( String )otherLink->item;
				if ( ( strcmp( otherContents, "" ) ) ) {
					contentsLink->item = otherLink->item;
					otherLink = Find_n_th_Item( number, gFcLineNumberList );
					numberLink = Find_n_th_Item( i+1, gFcLineNumberList );
					numberLink->item = otherLink->item;
					ListIterator( gFcNameList, PrintBroadening, gFcContentsList );	
				}
				else {
					AddString( errorList, name );
					if ( ListIterator( errorList, FindString, otherLabel ) ) {
						fprintf( stderr, "# you have to specify an equation for the broadening factor for the reaction labeled %s or %s\n", name, otherLabel );
					}
					err = TRUE;
				}
			}
			else {
				AddString( errorList, name );
				fprintf( stderr, "# you have to specify an equation for the broadening factor for the reaction labeled %s\n", name );
				err = TRUE;
			}
		}
		else {
			if ( name[len-1] == 'f' ||  name[len-1] == 'b' ) {
				number = NumberOfListItem( gFcNameList, FindString, otherLabel );
				otherLink = Find_n_th_Item( number, gFcContentsList );
				otherContents = ( String )otherLink->item;
				checkBuffer = ( String ) malloc( sizeof( contents ) + 1 );
				otherCheckBuffer = ( String ) malloc( sizeof( otherContents ) + 1 );
				strcpy( checkBuffer, contents );
				SkipWhites( checkBuffer );
				strcpy( otherCheckBuffer, contents );
				SkipWhites( otherCheckBuffer );
				if ( ( strcmp( checkBuffer, otherCheckBuffer ) ) != 0 ) {
					fprintf( stderr, "# forward and backward reaction should use the same broadening factor\n", name );
					fprintf( stderr, "# i will use\n\nFc = %s\n\n# for both, %s and %s\n", contents, name, otherLabel );
					otherLink->item = contentsLink->item;
					ListIterator( gFcNameList, PrintBroadening, gFcContentsList );	
				}
				free( checkBuffer );
				free( otherCheckBuffer );
			}
		}
	}
	FreeList( errorList, FreeString );
	
	return err;
}
*/

int CheckBroadening( void )
{
	if ( !( ListIterator( gReactionList, CheckOneBroadening, gFcNameList) ) ) {
		return 0;
	} 
	else {
		return 1;
	}
}

int CheckOneBroadening( LinkPtr link, void *var )
{
	ReactionPtr		reaction = link->item;
	ListPtr			labelList = var;
	String			label = NULL;
	int				len;
	Flag			err = FALSE;
	static ListPtr	errorList = NULL;
	
	if ( !errorList ) {
		errorList = NewList( "errorList" );
	}
	
	if ( reaction->aInf && !( reaction->fca || reaction->fcb || reaction->fcc ) ) {
		/* find reaction label */
		len = strlen( reaction->label );
		label = ( String ) malloc( len + 1 );
		strcpy( label, reaction->label );
		if ( gOptions->useForwardBackward ) {
			if ( label[len-1] == 'f' || label[len-1] == 'b' ) {
				label[len-1] = '\0';
			}
		}
		if ( !( ListIterator( labelList, FindString, label ) ) ) {
			if ( !( ListIterator( errorList, FindString, label ) ) ) {			
				AddString( errorList, label );
				if ( gOptions->useForwardBackward ) {
					fprintf( stderr, "# error: you have to specify a broadeningfactor for reaction %sf or %sb\n\n", label, label );
				}
				else {
					fprintf( stderr, "# error: you have to specify a broadeningfactor for reaction %s\n\n", label );
				}
			}
			err = TRUE;
		}
		free( label );
	}
	if ( link == gReactionList->lastItem ) {
		FreeList( errorList, FreeString );
		errorList = NULL;
	}
	
	return err;
}

void ComputeLindemann()
{
	LinkPtr 	link1 = NULL;
	LinkPtr 	link2 = NULL;
	ReactionPtr	reaction1 = NULL;
	ReactionPtr	reaction2 = NULL;
	int 		i;
	int			len;
	char		label[32];
	
	for ( i = 0; i < NumberOfItems( gReactionList ); ++i ) {
		link1 = Find_n_th_Item( i+1, gReactionList );
		reaction1 = ( ReactionPtr )link1->item;
		
		len = strlen( reaction1->label );
		if ( reaction1->label[len-1] == 'f' && reaction1->withLindemann == TRUE ) {
			strcpy( label, reaction1->label );
			label[len-1] = 'b';
			
			if ( !( link2 = ListIterator( gReactionList, FindReaction, label ) ) ) {
				fprintf( stderr, "# error: i can't find the reaction labeled %s which should be the redirected reaction to\n#\n#", label );
				PrintReaction( link1, gReactionList );
				exit( 2 );
			}
			
			reaction2 = ( ReactionPtr )link2->item;
			
			if ( reaction1->lindemannNumber != reaction2->lindemannNumber )   {
				if ( reaction1->lindemannNumber >= 0 ) {
					reaction2->lindemannNumber = reaction1->lindemannNumber;
				}
				else {
					reaction1->lindemannNumber = reaction2->lindemannNumber;
				}
			} 

			if ( reaction1->fca || reaction1->fcb || reaction1->fcc ) {
				reaction2->fca = reaction1->fca;
				reaction2->fcb = reaction1->fcb;
				reaction2->fcc = reaction1->fcc;
				reaction2->fcTa = reaction1->fcTa;
				reaction2->fcTb = reaction1->fcTb;
				reaction2->fcTc = reaction1->fcTc;
			}
			else if ( reaction2->fca || reaction2->fcb || reaction2->fcc ) {
				reaction1->fca = reaction2->fca;
				reaction1->fcb = reaction2->fcb;
				reaction1->fcc = reaction2->fcc;
				reaction1->fcTa = reaction2->fcTa;
				reaction1->fcTb = reaction2->fcTb;
				reaction1->fcTc = reaction2->fcTc;
			}

			if ( reaction1->a == 0.0 && reaction1->aInf != 0.0 && 
						reaction2->a != 0.0 && reaction2->aInf != 0.0) {
				reaction1->a = reaction1->aInf / reaction2->aInf * reaction2->a;
				reaction1->n = reaction1->nInf - reaction2->nInf + reaction2->n;
				reaction1->e = reaction1->eInf - reaction2->eInf + reaction2->e;
			}
			else if ( reaction1->a != 0.0 && reaction1->aInf == 0.0 && 
						reaction2->a != 0.0 && reaction2->aInf != 0.0) {
				reaction1->aInf = reaction1->a / reaction2->a * reaction2->aInf;
				reaction1->nInf = reaction1->n - reaction2->n + reaction2->nInf;
				reaction1->eInf = reaction1->e - reaction2->e + reaction2->eInf;
			}
			else if ( reaction1->a != 0.0 && reaction1->aInf != 0.0 && 
						reaction2->a == 0.0 && reaction2->aInf != 0.0) {
				reaction2->a = reaction2->aInf / reaction1->aInf * reaction1->a;
				reaction2->n = reaction2->nInf - reaction1->nInf + reaction1->n;
				reaction2->e = reaction2->eInf - reaction1->eInf + reaction1->e;
			}
			else if ( reaction1->a != 0.0 && reaction1->aInf != 0.0 && 
						reaction2->a != 0.0 && reaction2->aInf == 0.0) {
				reaction2->aInf = reaction2->a / reaction1->a * reaction1->aInf;
				reaction2->nInf = reaction2->n - reaction1->n + reaction1->nInf;
				reaction2->eInf = reaction2->e - reaction1->e + reaction1->eInf;
			}
			else if ( reaction1->a != 0.0 && reaction1->aInf != 0.0 && 
						reaction2->a != 0.0 && reaction2->aInf != 0.0) {
				continue;		
			}
			else {
				fprintf( stderr, "# error: not enough constants specified for the reactions\n\n" );
				PrintReaction( link1, gReactionList );
				fprintf( stderr, "# and\n\n " );
				PrintReaction( link2, gReactionList );
				exit( 2 );
			}
		}
	}		
}

int WriteFcFile( LinkPtr link, void *var )
{
	FILE		*fp = NULL;
	FILE		*fph = NULL;
	char		*funcName = NULL;
	char		cName[] = "BroadeningFactors.cp";
	char		hName[] = "BroadeningFactors.h";
	static Flag	firstOpen = TRUE;
	String		label = link->item;
	ListPtr 	list = var;
	LinkPtr		broadLink = NULL;
	String		broadening = NULL;
	LinkPtr		broadLink2 = NULL;
	String		broadening2 = NULL;
	int			number;
	int 		i;
	LinkPtr		numberLink = NULL;
	int			lineNumber;
	
	number = NumberOfListItem( gFcNameList, FindString, label );
	broadLink = Find_n_th_Item( number, list );
	broadening = ( String )broadLink->item;
	funcName = ( char * ) malloc( strlen( label ) + 3 );
	strcpy( funcName, "Fc" );
	strcat( funcName, label );
	numberLink = Find_n_th_Item( number, gFcLineNumberList );
	lineNumber = *( int * )numberLink->item;

	if ( firstOpen ) {

		fp = fopen( cName, "w" );
		fph = fopen( hName, "w" );

		fprintf( fph,	"#define MECHANISM \"%s\"\n", gOptions->base );
/*		fprintf( fph,	"#include <stdlib.h>\n"
						"#include <math.h>\n"
						"#include \"ArrayManager.h\"\n\n" );*/
		fprintf( fph,	"#include \"FlameMaster.h\"\n\n" );
		fprintf( fph, "/*	Mechanism file: \"%s\"	*/\n\n", gOptions->inName );
		fprintf( fph, "typedef Double (*BFFunction)(Double T);\n\n" );
		
		fprintf( fp,	"/*  The functions in this file compute the broadening factors.\n"
						" *  The names of the functions start with \"Fc\" followed by\n"
						" *  the label of the corresponding reaction.\n"
						" */\n\n"
						"#include\"BroadeningFactors.h\"\n\n" );
		
		fprintf( fp, "/*	Mechanism file: \"%s\"	*/\n\n", gOptions->inName );
		
		fprintf( fph, 	"/* prototypes */\n" );

		fprintf( fp, 	"BFFunction gBroadening[%d] = { ", gFcNameList->items );
		for ( i = 0; i < gFcNameList->items; ++i ) {
			broadLink2 = Find_n_th_Item( i+1, gFcNameList );
			broadening2 = ( String )broadLink2->item;
			if ( i ) {
				fprintf( fp, 	", " );
			}
			fprintf( fp, 	"Fc%s", broadening2 );
			fprintf( fph, 	"Double Fc%s( Double T );\n", broadening2 );
		}
		fprintf( fph, 	"Double FcErr( Double T );\n" );
		fprintf( fp, 	" };\n\n" );
		fprintf( fph, 	"\n\nextern BFFunction gBroadening[%d];\n", gFcNameList->items );

		fclose( fph );
		fprintf( fp, "#ifndef MECHANISM\n" );
		fprintf( fp, "#define MECHANISM \"\"\n" );
		fprintf( fp, "#endif\n\n" );
		fprintf( fp, "void TReaction::CheckBroadeningFactors( const char *mechName )\n" );
		fprintf( fp, "{\n" );
		fprintf( fp, "\tchar	*name = new char[strlen( MECHANISM ) + 6];\n" );
		fprintf( fp, "\tsprintf( name, \"/%%s.pre\", MECHANISM );\n" );
		fprintf( fp, "\tif ( strstr( mechName, name ) == NULL ) {\n" );
		fprintf( fp, "\t\tfor ( int i = 0; i < %d; ++i ) {\n", gFcNameList->items );
		fprintf( fp, "\t\t\tgBroadening[i] = FcErr;\n" );
		fprintf( fp, "\t\t}\n" );
		fprintf( fp, "\t}\n" );
		fprintf( fp, "}\n\n" );

		fprintf( fp, 	"Double FcErr( Double /*T*/ )\n" );
		fprintf( fp, 	"{\n" );
		fprintf( fp, 	"\tfprintf( stderr, \"#error: wrong broadening factors (%%s) linked to program\\n\", MECHANISM );\n" );
		fprintf( fp, 	"\texit( 2 );\n" );
		fprintf( fp, 	"\n\treturn 0;\n"
						"}\n\n" );

#ifdef applec
		fsetfileinfo( hName, 'MPS ','TEXT' );
		fsetfileinfo( cName, 'MPS ','TEXT' );
#endif
		firstOpen = FALSE;

	}
	else {
		fp = fopen( cName, "a" );
	}
	fprintf( fp, 	"Double %s( Double T )\n", funcName );
	fprintf( fp, 	"{\n" 
					"#line %d \"%s\"\n\n", lineNumber-1, gOptions->inName );
	if ( strcmp( broadening, "" ) == 0 ) {
		fprintf( fp, 	"\tfprintf( stderr, \"#error: no broadening function specified\\n\" );\n" );
		fprintf( fp, 	"\texit( 2 );\n" );
	}
	else {
		fprintf( fp, 	"\tT = %s;\n", broadening );
	}
	fprintf( fp, 	"\n\treturn T;\n"
					"}\n\n" );
	fclose( fp );
	return 0;
}

int CheckMolarMass( LinkPtr link, void *var )
{
	int				i;
	LinkPtr			atomLink = NULL;
	AtomsPtr 		atom = NULL;
	const Double	mC = 12.01;
	const Double	mH = 1.008;
	const Double	mN = 14.01;
	const Double	mO = 16.00;
	const Double	mAr = 39.948;
	const Double	mBr = 79.90900;
	const Double	mF = 18.99840;
	const Double	mCl = 35.453;
	const Double	mHe = 4.00;
	Double			molarMass = 0.0;
	int 			err = FALSE;
	static Flag		headerOut = FALSE;
	
	SpeciesPtr	species = link->item;
	var = var;
	
	for ( i = 0, atomLink = gAtomsList->firstItem; 
			i < gAtomsList->items; ++i, atomLink = atomLink->next ) {
		atom = (AtomsPtr)atomLink->item;
		
		if ( strcmp( atom->name, "C" ) == 0 ) {
			molarMass += species->composition->vec[i] * mC;
		}
		else if ( strcmp( atom->name, "H" ) == 0 ) {
			molarMass += species->composition->vec[i] * mH;
		}
		else if ( strcmp( atom->name, "N" ) == 0 ) {
			molarMass += species->composition->vec[i] * mN;
		}
		else if ( strcmp( atom->name, "O" ) == 0 ) {
			molarMass += species->composition->vec[i] * mO;
		}
		else if ( strcmp( atom->name, "AR" ) == 0 ) {
			molarMass += species->composition->vec[i] * mAr;
		}
		else if ( strcmp( atom->name, "BR" ) == 0 ) {
			molarMass += species->composition->vec[i] * mBr;
		}
		else if ( strcmp( atom->name, "F" ) == 0 ) {
			molarMass += species->composition->vec[i] * mF;
		}
		else if ( strcmp( atom->name, "CL" ) == 0 ) {
			molarMass += species->composition->vec[i] * mCl;
		}
		else if ( strcmp( atom->name, "HE" ) == 0 ) {
			molarMass += species->composition->vec[i] * mHe;
		}
		else {
			if ( fabs( ( double ) species->composition->vec[i] ) > 1.0e-10 ) {
				/*fprintf( stderr, "#error: molar mass of atom %s unknown\n", atom->name );
				err = TRUE;*/
				molarMass = -1.0;
				break;
			}
		}
	}
		
	if ( fabs( ( molarMass - species->molarMass ) / molarMass ) > 1.0e-3 ) {
		if ( !headerOut ) {
/*			fprintf( stderr, "### Warning: input molar mass is different from computed "
							"molar mass for the following species:\n" );*/
			headerOut = TRUE;
		}
/*		fprintf( stderr, "\t%-25s:\tM_i = %7.3f\tM_c = %7.3f\t", species->name, species->molarMass, molarMass );*/
		if ( molarMass > 0.0 ) {
			species->molarMass = molarMass;
			/* compute mu_Coeff */
/*			species->muCoeff = 2.6693e-6 * sqrt( fabs(species->molarMass) ) / (species->sigma * species->sigma);*/
/*			fprintf( stderr, "M_c choosen\n" );*/
		}
		else {
/*			fprintf( stderr, "M_i choosen\n" );*/
		}
	}
	
	return 0;
}

int CheckSigmaEps( LinkPtr link, void *var )
{
	Double			molarMass = 0.0;
	int 			err = FALSE;
	static Flag		headerOut = FALSE;
	
	SpeciesPtr	species = link->item;
	var = var;
			
	if ( species->sigma <= 0.0 || species->k_over_eps <= 0.0 ) {
		if ( species->sigma <= 0.0 ) {
			species->sigma = 1.234 * pow( species->molarMass, 0.33 );
/*			species->muCoeff = 2.6693e-6 * sqrt( fabs(species->molarMass) ) / (species->sigma * species->sigma);*/
		}
		if ( species->k_over_eps <= 0.0 ) {
			species->k_over_eps = 1.0 / ( 37.15 * pow( species->molarMass, 0.58 ) );
/*			fprintf( stderr, "\tuse eps/k = %7.3f K", 1.0/species->k_over_eps );*/
		}
/*		fprintf( stderr, "\n" );*/
	}
	
/*	fprintf( stderr, "%s: M = %g\tk_over_eps = %g\tsigma = %g\tmuCoeff = %g\n",*/
/*			species->name, species->molarMass, species->k_over_eps, species->sigma, species->muCoeff );*/

	return 0;
}

int ComputeMuCoeff(LinkPtr link, void *var)
{ 
	Double			molarMass = 0.0;
	int 			err = FALSE;
	static Flag		headerOut = FALSE;
	
	SpeciesPtr	species = link->item;
	var = var;

/*	if ( species->muCoeff <= 0.) {*/
/*		ListIterator(gSpeciesList, CheckSigmaEps, NULL);*/
		species->muCoeff = 2.6693e-6 * sqrt( species->molarMass ) / (species->sigma * species->sigma);
/*		fprintf( stderr, "### Warning: muCoeff =\t%d \tsigma =\t%d \n ", species->muCoeff, species->sigma);*/
/*	}*/

	return 0;
}

