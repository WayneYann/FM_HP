
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

#define TROE

int FindAtom( LinkPtr link, void *var )
{
	char			*name = var;
	AtomsPtr		atom = link->item;
	
	if ( strcmp( name, atom->name) == 0 ) {
		return 1;
	} 
	else {
		return 0;
	}
}

int FindSpeciesName( LinkPtr link, void *var )
{
	char		*name = var;
	SpeciesPtr	species = link->item;
	
	if ( strcmp( name, species->name) == 0 ) {
		return 1;
	} 
	else {
		return 0;
	}
}

int FindSpeciesNumber( LinkPtr link, void *var )
{
	int			number = *( ( int * ) var );
	SpeciesPtr	species = link->item;
	
	if ( number == species->number ) {
		return 1;
	} 
	else {
		return 0;
	}
}

int FindSpecies( LinkPtr link, void *var )
{
	int			i;
	char		*name = var;
	char		iso1[32], iso2[32];
	int 		iso1Len, iso2Len;
	SpeciesPtr	species = link->item;
	
/*  test equality of molecular structure  */
	if ( strchr( name, '-' ) || strchr( species->name, '-' ) ) {
		if ( strchr( name, '-' ) ) {
			iso1Len = strcspn( name, "-");
			strncpy( iso1, name, iso1Len );
			iso1[iso1Len] = '\0';
		}
		else {
			return 0;
		}
		if ( strchr( species->name, '-' ) ) {
			iso2Len = strcspn( species->name, "-");
			strncpy( iso2, species->name, iso2Len );
			iso2[iso2Len] = '\0';
		}
		else {
			return 0;
		}
		UpperString( iso1 );
		UpperString( iso2 );
		if ( strcmp( iso1, iso2 ) ) {
			return 0;
		}
	}

/*  test equality of composition  */
	for ( i = 0; i < gComposition->len; ++i ) {
		if ( species->composition->vec[i] != gComposition->vec[i] ) {
			return 0;
		}
	}
	return 1;
}

int FindReaction( LinkPtr link, void *var )
{
	char		*name = var;
	ReactionPtr	reaction = link->item;
	
	if ( strcmp( name, reaction->label ) == 0 ) {
		return 1;
	} 
	else {
		return 0;
	}
}

int FindThirdBody( LinkPtr link, void *var )
{
	char			*name = var;
	ThirdBodyPtr	thirdBody = link->item;
	
	if ( strcmp( name, thirdBody->name) == 0 ) {
		return 1;
	} 
	else {
		return 0;
	}
}

int FindDimension( LinkPtr link, void *var )
{
	char				*name = var;
	DimensionPtr		dimension = link->item;
	
	if ( strcmp( name, dimension->name) == 0 ) {
		return 1;
	} 
	else {
		return 0;
	}
}

int FindParameter( LinkPtr link, void *var )
{
	char			*name = var;
	ParameterPtr	parameter = link->item;
	
	if ( strcmp( name, parameter->name) == 0 ) {
		return 1;
	} 
	else {
		return 0;
	}
}

int PrintAtoms( LinkPtr link, void *var )
{
	ListPtr 		list = var;
	AtomsPtr 		atom = link->item;

	if ( link == list->firstItem ) {
		ListInfo( list, stdout );	
		fprintf( stdout, "Used atoms are:\n");
	}

	if ( link->next ) {
		fprintf( stdout, "\tNo.%d: %s\n", atom->number, atom->name );
	}
	else {
		fprintf( stdout, "\tNo.%d: %s\n\n", atom->number, atom->name );
	}
	return 0;
}

int PrintSpeciesNames( LinkPtr link, void *var )
{
	ListPtr 		list = var;
	SpeciesPtr 		species = link->item;
	LinkPtr			atomLink = NULL;
	AtomsPtr 		atom = NULL;
	int 			i;
	Double			*dPtr;

	fprintf( stdout, "%s\n", species->name );

	return 0;
}



int PrintSpecies( LinkPtr link, void *var )
{
	ListPtr 		list = var;
	SpeciesPtr 		species = link->item;
	LinkPtr			atomLink = NULL;
	AtomsPtr 		atom = NULL;
	int 			i;
	Double			*dPtr;

	if ( link == list->firstItem ) {
		ListInfo( list, stdout );	
		fprintf( stdout, "Used species are:\n");
	}

	fprintf( stdout, "%s is species no. %d\n", species->name, species->number );
	if ( species->isSteadyState ) {
		fprintf( stdout, "\tspecies is steady state\n" );
	}
/*	for ( i = 0, atomLink = gAtomsList->firstItem; 
			i < gAtomsList->items; ++i, atomLink = atomLink->next ) {
		atom = (AtomsPtr)atomLink->item;
		if ( !i ) {
			fprintf( stdout, "\tComposition: %d %s", species->composition->vec[i], atom->name );
		}
		else {
			fprintf( stdout, "\t%d %s", species->composition->vec[i], atom->name );
		}
	}*/
	PrintSpeciesComposition( species, stdout );
	
	fprintf( stdout, "\tCoefficients of the nasa polynomials ( hot[7], cold[7] ): " );
	dPtr = species->coeffHot;
	fprintf( stdout, "%-14e\t%-14e\t%-14e\t%-14e\t%-14e\t%-14e\t%-14e\t", 
		dPtr[0], dPtr[1], dPtr[2], dPtr[3], dPtr[4], dPtr[5], dPtr[6] );
	dPtr = species->coeffCold;
	fprintf( stdout, "%-14e\t%-14e\t%-14e\t%-14e\t%-14e\t%-14e\t%-14e\n", 
		dPtr[0], dPtr[1], dPtr[2], dPtr[3], dPtr[4], dPtr[5], dPtr[6] );
	
	fprintf( stdout, "\tM = %-14e\tk/eps = %-14e\tsigma = %-14e\tmuCoeff = %-14e\n",
		species->molarMass, species->k_over_eps, species->sigma, species->muCoeff );
	
	fprintf( stdout, "\tconstant d_Coefficients are: " );
	for ( i = 0; i < species->dCoeff->len; ++i ) {
		fprintf( stdout, "%-14e\t", species->dCoeff->vec[i] );
	}	
	fprintf( stdout, "\n\tconstant omega_Coefficients are: " );
	for ( i = 0; i < species->omegaCoeff->len; ++i ) {
		fprintf( stdout, "%-14e\t", species->omegaCoeff->vec[i] );
	}	
	fprintf( stdout, "\n\n" );
	return 0;
}

void PrintSpeciesComposition( SpeciesPtr species, FILE *fp )
{
	int i;
	LinkPtr			atomLink = NULL;
	AtomsPtr 		atom = NULL;

	for ( i = 0, atomLink = gAtomsList->firstItem; 
			i < gAtomsList->items; ++i, atomLink = atomLink->next ) {
		atom = (AtomsPtr)atomLink->item;
		if ( !i ) {
			fprintf( fp, "\tComposition: %d %s", species->composition->vec[i], atom->name );
		}
		else {
			fprintf( fp, "\t%d %s", species->composition->vec[i], atom->name );
		}
	}
	fprintf( fp, "\n" );
}

int PrintReaction( LinkPtr link, void *var )
{
	ListPtr 		list = var;
	ReactionPtr 	reaction = link->item;
	LinkPtr			speciesLink = NULL;
	LinkPtr			thirdBodyLink = NULL;
	SpeciesPtr 		species = NULL;
	int 			i;
	Flag			first = TRUE;
	char			*thirdBody = NULL;

	if ( link == list->firstItem ) {
		ListInfo( list, stdout );	
	}

# ifdef applec
	SpinCursor( 1 );
#endif
	fprintf( stdout, "Reaction \"%s\" (id=%d) (no. %d):\n", reaction->label, reaction->id
		, NumberOfListItem( list, FindReaction, reaction->label )-1 );
	fprintf( stdout, "\tOrder of Reaction is %d\n", reaction->orderOfReaction );

/*  print left side of equation  */
	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[i] );
		if ( !speciesLink ) FatalError( "something wrong in PrintReaction()" );
		species = ( SpeciesPtr ) speciesLink->item;

		if ( reaction->speciesCoeff[i] > 0.0 ) {
			if ( first ) {
				if ( fabs( reaction->speciesCoeff[i] ) != 1.0 ) {
					fprintf( stdout, "\t%g %s", fabs( reaction->speciesCoeff[i] ), species->name );
					first = FALSE;
				}
				else {
					fprintf( stdout, "\t%s", species->name );
					first = FALSE;
				}
			}
			else {
				if ( fabs( reaction->speciesCoeff[i] ) != 1.0 ) {
					fprintf( stdout, " + %g %s", fabs( reaction->speciesCoeff[i] ), species->name );
				}
				else {
					fprintf( stdout, " + %s", species->name );
				}
			}
		}
	}

	if ( reaction->withThirdBody ) {
		if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
			thirdBody = ( char * ) thirdBodyLink->item;
			fprintf( stdout, " + %s", thirdBody );
		}
		else {
			fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
			exit( 2 );
		}
	}

/*  print assignment operator  */
	if ( reaction->partialEquilibrium ) {
		fprintf( stdout, " == " );
	}
	else {
		fprintf( stdout, " -> " );
	}

/*  print right side of equation  */
	first = TRUE;
	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[i] );
		if ( !speciesLink ) FatalError( "something wrong in PrintReaction()" );
		species = ( SpeciesPtr ) speciesLink->item;

		if ( reaction->speciesCoeff[i] < 0.0 ) {
			if ( first ) {
				if ( fabs( reaction->speciesCoeff[i] ) != 1.0 ) {
					fprintf( stdout, "%g %s", fabs( reaction->speciesCoeff[i] ), species->name );
					first = FALSE;
				}
				else {
					fprintf( stdout, "%s", species->name );
					first = FALSE;
				}
			}
			else {
				if ( fabs( reaction->speciesCoeff[i] ) != 1.0 ) {
					fprintf( stdout, " + %g %s", fabs( reaction->speciesCoeff[i] ), species->name );
				}
				else {
					fprintf( stdout, " + %s", species->name );
				}
			}
		}
	}

	if ( reaction->withThirdBody ) {
		fprintf( stdout, " + %s", thirdBody );
	}

	fprintf( stdout, "\n\n" );
	fprintf( stdout, "\tA = %g\n\tn = %g\n\tE = %g\n\n"
			, reaction->a
			, reaction->n, reaction->e );
	if ( reaction->withLindemann ) {
		fprintf( stdout, "\tA_inf = %g\n\tn_inf = %g\n\tE_inf = %g\n", reaction->aInf, reaction->nInf, reaction->eInf );
		fprintf( stdout, "\tthis reaction uses the broadeningfactor no. %d\n\n", reaction->lindemannNumber );
		fprintf( stdout, "\tfca = %g\tfcTa = %g\n\tfcb = %g\tfcTb = %g\n\tfcc = %g\tfcTc = %g\n\n"
						, reaction->fca, reaction->fcTa
						, reaction->fcb, reaction->fcTb
						, reaction->fcc, reaction->fcTc );
	}
		
	return 0;
}

int PrintMagicRateCoeff( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	const Double	moleM1TokmoleM1 = 1000.0;
	char			buff[128];
	LinkPtr			thirdBodyLink = NULL;
	char			*thirdBody = NULL;

	if ( !reaction->withLindemann ) {
		fprintf( fp, "\tk[%s%s] = " 
				, gPrefixReactions
				, CSymbol(reaction->label) );
		PrintOneRateCoeff( reaction, kK, fp );
	}
	else {
/* 		fprintf( fp, "\tDouble\tk%s0 = " , CSymbol(reaction->label) ); */
/* 		PrintOneRateCoeff( reaction, kK0, fp ); */
/* 		fprintf( fp, "\tDouble\tk%sInf = " , CSymbol(reaction->label) ); */
/* 		PrintOneRateCoeff( reaction, kKinf, fp ); */
		fprintf( fp, "\tkTroe0 = " );
		PrintOneRateCoeff( reaction, kK0, fp );
		fprintf( fp, "\tkTroeInf = " );
		PrintOneRateCoeff( reaction, kKinf, fp );
		if ( !reaction->hasForward ) {
			strcpy( buff, CSymbol(reaction->label) );
			if ( reaction->hasBackward ) {
				buff[strlen( buff )-1] = '\0';
			}
			fprintf( fp, "\tfcTroe = " );
			PrintBroadeningSource( reaction, buff, fp );
		}
		if ( reaction->withThirdBody ) {
			thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList );
			thirdBody = ( char * ) thirdBodyLink->item;
			fprintf( fp, "\tk[%s%s] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf\n"
				, gPrefixReactions
				 , CSymbol(reaction->label) );
			fprintf( fp, "\t\t\t\t, fcTroe, M[%s%s] );\n"
						, gPrefixThirdBody, CSymbol( thirdBody ) );
		}
		else {
			fprintf( fp, "\tk[%s%s] = GetLindRateCoeff( temp, pressure, kTroe0, kTroeInf\n"
				"\t\t\t\t, fcTroe, -1.0 );\n"
				, gPrefixReactions
				 , CSymbol(reaction->label) );
		}
	}

	return 0;
}


/* int PrintMagicRateCoeffF77( LinkPtr link, void *var ) */
/* { */
/* 	FILE 			*fp = ( FILE * )var; */
/* 	ReactionPtr 	reaction = link->item; */
/* 	const Double	moleM1TokmoleM1 = 1000.0; */
/* 	char			buff[128]; */
/* 	LinkPtr			thirdBodyLink = NULL; */
/* 	char			*thirdBody = NULL; */
/*  */
/* 	if ( !reaction->withLindemann ) { */
/* 		fprintf( fp, "      K(%s%s) = "  */
/* 				, gPrefixReactionsF77 */
/* 				, CSymbol(reaction->label) ); */
/* 		PrintOneRateCoeffF77( reaction, kK, fp ); */
/* 	} */
/* 	else { */
/* 		fprintf( fp, "      K%s0 = " , CSymbol(reaction->label) ); */
/* 		PrintOneRateCoeffF77( reaction, kK0, fp ); */
/* 		fprintf( fp, "      K%sINF = " , CSymbol(reaction->label) ); */
/* 		PrintOneRateCoeffF77( reaction, kKinf, fp ); */
/* 		strcpy( buff, CSymbol(reaction->label) ); */
/* 		if ( reaction->hasBackward || reaction->hasForward ) { */
/* 			buff[strlen( buff )-1] = '\0'; */
/* 		} */
/* 		if ( !reaction->hasForward ) { */
/* 			fprintf( fp, "      FC%s = " , buff ); */
/* 			PrintBroadeningSourceF77( reaction, buff, fp ); */
/* 		} */
/* 		if ( reaction->withThirdBody ) { */
/* 			thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ); */
/* 			thirdBody = ( char * ) thirdBodyLink->item; */
/* 			fprintf( fp, "      K(%s%s) = GETLINDRATECOEFF( TEMP, PRESSURE, K%s0, K%sINF\n" */
/* 				, gPrefixReactionsF77 */
/* 				 , CSymbol(reaction->label) */
/* 						  , CSymbol(reaction->label) */
/* 						  , CSymbol(reaction->label) ); */
/* 			fprintf( fp, "     &    , FC%s, M(%s%s) )\n" */
/* 				, buff, gPrefixThirdBody, CSymbol( thirdBody ) ); */
/* 		} */
/* 		else { */
/* 			fprintf( fp, "      K(%s%s) = GETLINDRATECOEFF( TEMP, PRESSURE, K%s0, K%sINF\n" */
/* 				"     &    , FC%s, CONCDEFAULT )\n" */
/* 				, gPrefixReactionsF77 */
/* 				 , CSymbol(reaction->label) */
/* 						  , CSymbol(reaction->label) */
/* 						  , CSymbol(reaction->label) */
/* 				, buff ); */
/* 		} */
/* 	} */
/*  */
/* 	return 0; */
/* } */

int PrintMagicRateCoeffF77( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	const Double	moleM1TokmoleM1 = 1000.0;
	char			buff[128];
	LinkPtr			thirdBodyLink = NULL;
	char			*thirdBody = NULL;

	if ( !reaction->withLindemann ) {
		fprintf( fp, "      K(%s%s) = " 
				, gPrefixReactionsF77
				, CSymbol(reaction->label) );
		PrintOneRateCoeffF77( reaction, kK, fp );
	}
	else {
/* This change requested by Varun Mittal */
		fprintf( fp, "      K0TROE = " );
		PrintOneRateCoeffF77( reaction, kK0, fp );
		fprintf( fp, "      KINFTROE = " );
		PrintOneRateCoeffF77( reaction, kKinf, fp );
		strcpy( buff, CSymbol(reaction->label) );
		if ( reaction->hasBackward || reaction->hasForward ) {
			buff[strlen( buff )-1] = '\0';
		}
		if ( !reaction->hasForward ) {
			fprintf( fp, "      FCTROE = " );
			PrintBroadeningSourceF77( reaction, buff, fp );
		}
		if ( reaction->withThirdBody ) {
			thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList );
			thirdBody = ( char * ) thirdBodyLink->item;
			fprintf( fp, "      K(%s%s) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE\n"
				, gPrefixReactionsF77
				, CSymbol(reaction->label) );
			fprintf( fp, "     &    , FCTROE, M(%s%s) )\n"
				, gPrefixThirdBody, CSymbol( thirdBody ) );
		}
		else {
			fprintf( fp, "      K(%s%s) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE\n"
				"     &    , FCTROE, CONCDEFAULT )\n"
				, gPrefixReactionsF77
				, CSymbol(reaction->label) );
/* 		fprintf( fp, "      K%s0 = " , CSymbol(reaction->label) ); */
/* 		PrintOneRateCoeffF77( reaction, kK0, fp ); */
/* 		fprintf( fp, "      K%sINF = " , CSymbol(reaction->label) ); */
/* 		PrintOneRateCoeffF77( reaction, kKinf, fp ); */
/* 		if ( !reaction->hasForward ) { */
/* 			strcpy( buff, CSymbol(reaction->label) ); */
/* 			if ( reaction->hasBackward ) { */
/* 				buff[strlen( buff )-1] = '\0'; */
/* 			} */
/* 			fprintf( fp, "      FC%s = " , buff ); */
/* 			PrintBroadeningSourceF77( reaction, buff, fp ); */
/* 		} */
/* 		if ( reaction->withThirdBody ) { */
/* 			thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ); */
/* 			thirdBody = ( char * ) thirdBodyLink->item; */
/* 			fprintf( fp, "      K(%s%s) = GETLINDRATECOEFF( TEMP, PRESSURE, K%s0, K%sINF\n" */
/* 				, gPrefixReactionsF77 */
/* 				 , CSymbol(reaction->label) */
/* 						  , CSymbol(reaction->label) */
/* 						  , CSymbol(reaction->label) ); */
/* 			fprintf( fp, "     &    , FC%s, M(%s%s) )\n" */
/* 				, buff, gPrefixThirdBody, CSymbol( thirdBody ) ); */
/* 		} */
/* 		else { */
/* 			fprintf( fp, "      K(%s%s) = GETLINDRATECOEFF( TEMP, PRESSURE, K%s0, K%sINF\n" */
/* 				"     &    , FC%s, CONCDEFAULT )\n" */
/* 				, gPrefixReactionsF77 */
/* 				 , CSymbol(reaction->label) */
/* 						  , CSymbol(reaction->label) */
/* 						  , CSymbol(reaction->label) */
/* 				, buff ); */
		}
	}

	return 0;
}

int PrintMagicRateCoeffF77Vectorize( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	const Double	moleM1TokmoleM1 = 1000.0;
	char			buff[128];
	LinkPtr			thirdBodyLink = NULL;
	char			*thirdBody = '\0';

	if ( reaction->withLindemann ) {
		fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);
		fprintf( fp, "      K0(K)=" );
		PrintOneRateCoeffF77Vectorize( reaction, kK0, fp );
		fprintf( fp, "      KINF(K)=" );
		PrintOneRateCoeffF77Vectorize( reaction, kKinf, fp );
		if ( !reaction->hasForward ) {
			fprintf( fp, "      FC(K)=" );
			PrintBroadeningSourceF77Vectorize( reaction, fp );
		}
		fprintf( fp, "%5d CONTINUE\n", gLabel);
		gLabel+=gLabelIncrement/2;


		fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);
		if( reaction->withThirdBody ) {/**/
			if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
				thirdBody = ( char * ) thirdBodyLink->item;
				fprintf( fp, "      CONC = M(K,%s%s)\n", gPrefixThirdBodyF77, CSymbol( thirdBody ) );
			}
		}
		else {
				fprintf( fp, "       CONC = PRESSURE / ( RGAS * TEMP(K) )\n" );
		}
		fprintf( fp, "      NTMP = 0.75D+00 - 1.27D+00 * DLOG10( FC(K) )\n"  );
		fprintf( fp, "      CCOEFF = - 0.4D+00 - 0.67D+00 * DLOG10( FC(K) )\n" );
		fprintf( fp, "      DCOEFF = 0.14D+00\n" );
		fprintf( fp, "      LGKNULL = DLOG10( K0(K) )\n"  );
		fprintf( fp, "      K0K = K0(K) * CONC / KINF(K)\n" );
		fprintf( fp, "      F = ( LGKNULL + CCOEFF )\n" );
		fprintf( fp, "     & / ( NTMP - DCOEFF * ( LGKNULL + CCOEFF ) )\n" );
		fprintf( fp, "      F =  FC(K) ** ( 1.0D+00 / ( F * F + 1.0D+00 ) )\n"  );
		fprintf( fp, "      KR(K,%s%s) = KINF(K) * F * K0K / ( 1.0D+00 + K0K )\n\n"
						, gPrefixReactionsF77
						, reaction->label );
		fprintf( fp, "%5d CONTINUE\n", gLabel);
		fprintf( fp, "C\nC\n" );
		
		gLabel+=gLabelIncrement/2;

/*		fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);*/
/*		fprintf( fp, "      CONC = PRESSURE / ( RGAS * TEMP(K) )\n" );*/
/*		fprintf( fp, "      NTMP = 0.75D+00 - 1.27D+00 * DLOG10( FC(K) )\n" );*/
/*		fprintf( fp, "      K0K = K0(K) * CONC / KINF(K)\n" );*/
/*		fprintf( fp, "      KLK = K0K / ( 1.0D+00 + K0K )\n" );*/
/*		fprintf( fp, "      F = DLOG10( K0K ) / NTMP\n" );*/
/*		fprintf( fp, "      F = FC(K) ** ( 1.0D+00 / ( F * F + 1.0D+00 ) )\n" );*/
/*		fprintf( fp, "      KR(K,%s%s) = KINF(K) * F * KLK\n"*/
/*						, gPrefixReactionsF77*/
/*						, reaction->label );*/
/*		fprintf( fp, "%5d CONTINUE\n", gLabel);*/
/*		fprintf( fp, "C\nC\n" );*/
/*		*/
/*		gLabel+=gLabelIncrement/2;*/

	}
	
	return 0;
}


int PrintMagicRateCoeffF77New( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	const Double	moleM1TokmoleM1 = 1000.0;
	char			buff[128];

	if ( reaction->withLindemann ) {
		fprintf( fp, "      K%s0 = " , reaction->label );
		PrintOneRateCoeffF77( reaction, kK0, fp );
		fprintf( fp, "      K%sINF = " , reaction->label );
		PrintOneRateCoeffF77( reaction, kKinf, fp );
		if ( !reaction->hasForward ) {
			strcpy( buff, reaction->label );
			if ( reaction->hasBackward ) {
				buff[strlen( buff )-1] = '\0';
			}
			fprintf( fp, "      FC%s = " , buff );
			PrintBroadeningSourceF77( reaction, buff, fp );
		}
		fprintf( fp, "      K(%s%s) = GETLINDRATECOEFF( TEMP, PRESSURE, K%s0, K%sINF\n"
				"     &    , FC%s )\n"
				, gPrefixReactionsF77
				, reaction->label
				, reaction->label
				, reaction->label
				, buff );
	}
	
	return 0;
}

int PrintMagicRateCoeffF90( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	const Double	moleM1TokmoleM1 = 1000.0;
	char			buff[128];
	LinkPtr			thirdBodyLink = NULL;
	char			*thirdBody = NULL;

	if ( !reaction->withLindemann ) {
		fprintf( fp, "      K(%s%s) = " 
				, gPrefixReactionsF77
				, CSymbol(reaction->label) );
		PrintOneRateCoeffF90( reaction, kK, fp );
	}
	else {
 		fprintf( fp, "      K0TROE = " ); 
 		PrintOneRateCoeffF90( reaction, kK0, fp ); 
 		fprintf( fp, "      KINFTROE = " ); 
 		PrintOneRateCoeffF90( reaction, kKinf, fp ); 
 		strcpy( buff, CSymbol(reaction->label) ); 
 		if ( reaction->hasBackward || reaction->hasForward ) { 
 			buff[strlen( buff )-1] = '\0'; 
 		} 
 		if ( !reaction->hasForward ) { 
 			fprintf( fp, "      FCTROE = " ); 
 			PrintBroadeningSourceF90( reaction, buff, fp ); 
 		} 
 		if ( reaction->withThirdBody ) { 
 			thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ); 
 			thirdBody = ( char * ) thirdBodyLink->item; 
 			fprintf( fp, "      K(%s%s) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &\n" 
 				, gPrefixReactionsF77 
 				, CSymbol(reaction->label) ); 
 			fprintf( fp, "\t   , FCTROE, M(%s%s) )\n" 
 				, gPrefixThirdBody, CSymbol( thirdBody ) ); 
 		} 
 		else { 
 			fprintf( fp, "      K(%s%s) = GETLINDRATECOEFF( TEMP, PRESSURE, K0TROE, KINFTROE &\n" 
 				"\t   , FCTROE, CONCDEFAULT )\n" 
 				, gPrefixReactionsF77 
 				, CSymbol(reaction->label) );
		}
	}

 
/* 			fprintf( fp, "      K%s0 = " , CSymbol(reaction->label) ); */
/* 			PrintOneRateCoeffF77( reaction, kK0, fp ); */
/* 			fprintf( fp, "      K%sINF = " , CSymbol(reaction->label) ); */
/* 			PrintOneRateCoeffF77( reaction, kKinf, fp ); */
/* 			if ( !reaction->hasForward ) { */
/* 			  strcpy( buff, CSymbol(reaction->label) ); */
/* 			  if ( reaction->hasBackward ) { */
/* 			    buff[strlen( buff )-1] = '\0'; */
/* 			  } */
/* 			  fprintf( fp, "      FC%s = " , buff ); */
/* 			  PrintBroadeningSourceF77( reaction, buff, fp ); */
/* 			} */
/* 			if ( reaction->withThirdBody ) { */
/* 			  thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ); */
/* 			  thirdBody = ( char * ) thirdBodyLink->item; */
/* 			  fprintf( fp, "      K(%s%s) = GETLINDRATECOEFF( TEMP, PRESSURE, K%s0, K%sINF\n" */
/* 				   , gPrefixReactionsF77 */
/* 				   , CSymbol(reaction->label) */
/* 				   , CSymbol(reaction->label) */
/* 				   , CSymbol(reaction->label) ); */
/* 			  fprintf( fp, "     &    , FC%s, M(%s%s) )\n" */
/* 				   , buff, gPrefixThirdBody, CSymbol( thirdBody ) ); */
/* 			} */
/* 			else { */
/* 			  fprintf( fp, "      K(%s%s) = GETLINDRATECOEFF( TEMP, PRESSURE, K%s0, K%sINF\n" */
/* 				   "     &    , FC%s, CONCDEFAULT )\n" */
/* 				   , gPrefixReactionsF77 */
/* 				   , CSymbol(reaction->label) */
/* 				   , CSymbol(reaction->label) */
/* 				   , CSymbol(reaction->label) */
/* 				   , buff ); */
/* 			} */

	return 0;
}

void PrintBroadeningSource( ReactionPtr reaction, char *redLabel, FILE *fp )
{
	Flag	first = TRUE;

	if ( reaction->fca || reaction->fcb || reaction->fcc ) {
		if ( reaction->fca ) {
			if ( !first ) {
				fprintf( fp, " + " );
			}
			else {
				first = FALSE;
			}
			fprintf( fp, "%g * exp( -temp / %g )"
				, reaction->fca, reaction->fcTa );
		}
		if ( reaction->fcb ) {
			if ( !first ) {
				fprintf( fp, " + " );
			}
			else {
				first = FALSE;
			}
			fprintf( fp, "%g * exp( -temp / %g )"
				, reaction->fcb, reaction->fcTb );
		}
		if ( reaction->fcc ) {
			if ( !first ) {
				fprintf( fp, " + " );
			}
			else {
				first = FALSE;
			}
			fprintf( fp, "%g * exp( -%g / temp )"
				, reaction->fcc, reaction->fcTc );
		}
		fprintf( fp, ";\n" );
	}
	else {
		fprintf( fp, "Fc%s( temp );\n", redLabel );
	}
}

void PrintBroadeningSourceF77Vectorize( ReactionPtr reaction, FILE *fp )
{
	Flag	first = TRUE;

	if ( reaction->fca || reaction->fcb || reaction->fcc ) {
		if ( reaction->fca ) {
			if ( !first ) {
				fprintf( fp, "+" );
			}
			else {
				first = FALSE;
			}
			fprintf( fp, "%G*EXP(-TEMP(K)/%g)"
				, reaction->fca, reaction->fcTa );
		}
		if ( reaction->fcb ) {
			if ( !first ) {
				fprintf( fp, "+" );
			}
			else {
				first = FALSE;
			}
			fprintf( fp, "%G*EXP(-TEMP(K)/%g)"
				, reaction->fcb, reaction->fcTb );
		}
		if ( reaction->fcc ) {
			if ( !first ) {
				fprintf( fp, "+" );
			}
			else {
				first = FALSE;
			}
			if ( reaction->fcTc ) {
				fprintf( fp, "%G*EXP(-%g/TEMP(K))"
					, reaction->fcc, reaction->fcTc );
			}
			else{
				fprintf( fp, "%G", reaction->fcc );
			}
		}
		fprintf( fp, "\n" );
	}
	else {
		fprintf( fp, "\n### WARNING: Function not available ###\n" );
	}
}

void PrintBroadeningSourceF77( ReactionPtr reaction, char *redLabel, FILE *fp )
{
	Flag	first = TRUE;

	if ( reaction->fca || reaction->fcb || reaction->fcc ) {
		if ( reaction->fca ) {
			if ( !first ) {
				fprintf( fp, "\n     &     + " );
			}
			else {
				first = FALSE;
			}
			fprintf( fp, "%G * EXP( -TEMP / %g )"
				, reaction->fca, reaction->fcTa );
		}
		if ( reaction->fcb ) {
			if ( !first ) {
				fprintf( fp, "\n     &     + " );
			}
			else {
				first = FALSE;
			}
			fprintf( fp, "%G * EXP( -TEMP / %g )"
				, reaction->fcb, reaction->fcTb );
		}
		if ( reaction->fcc ) {
			if ( !first ) {
				fprintf( fp, "\n     &     + " );
			}
			else {
				first = FALSE;
			}
			fprintf( fp, "%G * EXP( -%g / TEMP )"
				, reaction->fcc, reaction->fcTc );
		}
		fprintf( fp, "\n" );
	}
	else {
		fprintf( fp, "Fc%s( temp ) );\n", redLabel );
	}
}

void PrintBroadeningSourceF90( ReactionPtr reaction, char *redLabel, FILE *fp )
{
	Flag	first = TRUE;

	if ( reaction->fca || reaction->fcb || reaction->fcc ) {
		if ( reaction->fca ) {
			if ( !first ) {
				fprintf( fp, " &\n\t   + " );
			}
			else {
				first = FALSE;
			}
			fprintf( fp, "%G * EXP( -TEMP / %g )"
				, reaction->fca, reaction->fcTa );
		}
		if ( reaction->fcb ) {
			if ( !first ) {
				fprintf( fp, " &\n\t   + " );
			}
			else {
				first = FALSE;
			}
			fprintf( fp, "%G * EXP( -TEMP / %g )"
				, reaction->fcb, reaction->fcTb );
		}
		if ( reaction->fcc ) {
			if ( !first ) {
				fprintf( fp, " &\n\t   + " );
			}
			else {
				first = FALSE;
			}
			fprintf( fp, "%G * EXP( -%g / TEMP )"
				, reaction->fcc, reaction->fcTc );
		}
		fprintf( fp, "\n" );
	}
	else {
		fprintf( fp, "Fc%s( temp ) );\n", redLabel );
	}
}

int PrintDeclarationBroadeningFactorsF77( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	char			buff[128];

	if ( reaction->withLindemann ) {
	  fprintf( fp, "      DOUBLE PRECISION KINFTROE, K0TROE\n");
	  if ( !reaction->hasForward ) {
	    fprintf( fp, "      DOUBLE PRECISION FCTROE\n");
	    return 1;
	  }
	}
	return 0;
}

int PrintDeclarationBroadeningFactorsF90( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	char			buff[128];

	if ( reaction->withLindemann ) {
		fprintf( fp, "\n      real(DP) ::  KINFTROE, K0TROE\n");
		if ( !reaction->hasForward ) {
		  fprintf( fp, "      real(DP) ::  FCTROE\n");
		  return 1;
		}
	}

	return 0;
}


void PrintOneRateCoeff( ReactionPtr reaction, ReactionRateType rateType, FILE *fp )
{
	Double	a, n, e;

	switch ( rateType ) {

	case kK:
		a = reaction->a * GetAConverter( reaction, TRUE ) 
			* GetATokmoleConverter( reaction, TRUE );
		n = reaction->n;
		e = reaction->e * GetEConverter( reaction, TRUE ) 
			* GetETokmoleConverter( reaction, TRUE );
		break;

	case kK0:
		a = reaction->a * GetAConverter( reaction, FALSE ) 
			* GetATokmoleConverter( reaction, FALSE );
		n = reaction->n;
		e = reaction->e * GetEConverter( reaction, FALSE ) 
			* GetETokmoleConverter( reaction, FALSE );
		break;

	case kKinf:
		a = reaction->aInf * GetAConverter( reaction, TRUE ) 
			* GetATokmoleConverter( reaction, TRUE );
		n = reaction->nInf;
		e = reaction->eInf * GetEConverter( reaction, TRUE ) 
			* GetETokmoleConverter( reaction, TRUE );		
		break;

	default:
		fprintf( stderr, "#error: unknown ReactionRateType in PrintOneRateCoeffF77\n" );
		exit(2);
	}

	fprintf( fp, "%.10E", a );
	if ( n && e ) {
		fprintf( fp, " * exp( %g * lgt %s %.10G / rt )"
				, n
				, ( e > 0.0 ) ? "-" : "+" 
				, fabs( e ) );
	}
	else if ( n ) {
		fprintf( fp, " * exp( %g * lgt )", n );
	}
	else if ( e ) {
		fprintf( fp, " * exp( %s%.10G / rt )"
				, ( e > 0.0 ) ? "-" : "" 
				, fabs( e ) );
	}
	fprintf( fp, ";\n" );
}

void PrintOneRateCoeffF77Vectorize( ReactionPtr reaction, ReactionRateType rateType, FILE *fp )
{
	Double	a, n, e;
	char	buffer[32];

	switch ( rateType ) {

	case kK:
		a = reaction->a * GetAConverter( reaction, TRUE ) 
			* GetATokmoleConverter( reaction, TRUE );
		n = reaction->n;
		e = reaction->e * GetEConverter( reaction, TRUE ) 
			* GetETokmoleConverter( reaction, TRUE );
		break;

	case kK0:
		a = reaction->a * GetAConverter( reaction, FALSE ) 
			* GetATokmoleConverter( reaction, FALSE );
		n = reaction->n;
		e = reaction->e * GetEConverter( reaction, FALSE ) 
			* GetETokmoleConverter( reaction, FALSE );
		break;

	case kKinf:
		a = reaction->aInf * GetAConverter( reaction, TRUE ) 
			* GetATokmoleConverter( reaction, TRUE );
		n = reaction->nInf;
		e = reaction->eInf * GetEConverter( reaction, TRUE ) 
			* GetETokmoleConverter( reaction, TRUE );		
		break;

	default:
		fprintf( stderr, "#error: unknown ReactionRateType in PrintOneRateCoeffF77\n" );
		exit(2);
	}

	ConvertCDoubleToF77Double( a, buffer );
	fprintf( fp, "%s", buffer );
	if ( n && e ) {
/*		fprintf( fp, "     & * TEMP ** %G * DEXP( %s%.16G / RT )"	*/
		fprintf( fp, "*DEXP(%G*LT(K)%s%.10G*RTI(K))"
				, n
				, ( e > 0.0 ) ? "-" : "+" 
				, fabs( e ) );
	}
	else if ( n ) {
/*		fprintf( fp, "     & * TEMP ** %G ", n );	*/
		fprintf( fp, "*DEXP(%G*LT(K))", n );
	}
	else if ( e ) {
		fprintf( fp, "*DEXP(%s%.10G * RTI(K))"
				, ( e > 0.0 ) ? "-" : "" 
				, fabs( e ) );
	}
	fprintf( fp, "\n" );
}

void PrintOneRateCoeffF77( ReactionPtr reaction, ReactionRateType rateType, FILE *fp )
{
	Double	a, n, e;
	char	buffer[32];

	switch ( rateType ) {

	case kK:
		a = reaction->a * GetAConverter( reaction, TRUE ) 
			* GetATokmoleConverter( reaction, TRUE );
		n = reaction->n;
		e = reaction->e * GetEConverter( reaction, TRUE ) 
			* GetETokmoleConverter( reaction, TRUE );
		break;

	case kK0:
		a = reaction->a * GetAConverter( reaction, FALSE ) 
			* GetATokmoleConverter( reaction, FALSE );
		n = reaction->n;
		e = reaction->e * GetEConverter( reaction, FALSE ) 
			* GetETokmoleConverter( reaction, FALSE );
		break;

	case kKinf:
		a = reaction->aInf * GetAConverter( reaction, TRUE ) 
			* GetATokmoleConverter( reaction, TRUE );
		n = reaction->nInf;
		e = reaction->eInf * GetEConverter( reaction, TRUE ) 
			* GetETokmoleConverter( reaction, TRUE );		
		break;

	default:
		fprintf( stderr, "#error: unknown ReactionRateType in PrintOneRateCoeffF77\n" );
		exit(2);
	}

	ConvertCDoubleToF77Double( a, buffer );
	fprintf( fp, "%s", buffer );
	if ( n && e ) {
/*		fprintf( fp, "     & * TEMP ** %G * DEXP( %s%.16G / RT )"	*/
		fprintf( fp, "\n     & * DEXP(%G * LT %s %.10G / RT)"
				, n
				, ( e > 0.0 ) ? "-" : "+" 
				, fabs( e ) );
	}
	else if ( n ) {
/*		fprintf( fp, "     & * TEMP ** %G ", n );	*/
		fprintf( fp, " * DEXP(%G * LT)", n );
	}
	else if ( e ) {
		fprintf( fp, " * DEXP(%s%.10G / RT)"
				, ( e > 0.0 ) ? "-" : "" 
				, fabs( e ) );
	}
	fprintf( fp, "\n" );
}

void PrintOneRateCoeffF90( ReactionPtr reaction, ReactionRateType rateType, FILE *fp )
{
	Double	a, n, e;
	char	buffer[32];

	switch ( rateType ) {

	case kK:
		a = reaction->a * GetAConverter( reaction, TRUE ) 
			* GetATokmoleConverter( reaction, TRUE );
		n = reaction->n;
		e = reaction->e * GetEConverter( reaction, TRUE ) 
			* GetETokmoleConverter( reaction, TRUE );
		break;

	case kK0:
		a = reaction->a * GetAConverter( reaction, FALSE ) 
			* GetATokmoleConverter( reaction, FALSE );
		n = reaction->n;
		e = reaction->e * GetEConverter( reaction, FALSE ) 
			* GetETokmoleConverter( reaction, FALSE );
		break;

	case kKinf:
		a = reaction->aInf * GetAConverter( reaction, TRUE ) 
			* GetATokmoleConverter( reaction, TRUE );
		n = reaction->nInf;
		e = reaction->eInf * GetEConverter( reaction, TRUE ) 
			* GetETokmoleConverter( reaction, TRUE );		
		break;

	default:
		fprintf( stderr, "#error: unknown ReactionRateType in PrintOneRateCoeffF77\n" );
		exit(2);
	}

	ConvertCDoubleToF77Double( a, buffer );
	fprintf( fp, "%s", buffer );
	if ( n && e ) {
/*		fprintf( fp, "     & * TEMP ** %G * DEXP( %s%.16G / RT )"	*/
		fprintf( fp, " &\n\t   * exp(%G * LT %s %.10G / RT)"
				, n
				, ( e > 0.0 ) ? "-" : "+" 
				, fabs( e ) );
	}
	else if ( n ) {
/*		fprintf( fp, "     & * TEMP ** %G ", n );	*/
		fprintf( fp, " * exp(%G * LT)", n );
	}
	else if ( e ) {
		fprintf( fp, " * exp(%s%.10G / RT)"
				, ( e > 0.0 ) ? "-" : "" 
				, fabs( e ) );
	}
	fprintf( fp, "\n" );
}

void PrintArrCoeffsF77( ReactionPtr reaction, ReactionRateType rateType, Double *a,
						Double *n, Double *e )
{
	switch ( rateType ) {

	case kK:
		*a = reaction->a * GetAConverter( reaction, TRUE ) 
			* GetATokmoleConverter( reaction, TRUE );
		*n = reaction->n;
		*e = reaction->e * GetEConverter( reaction, TRUE ) 
			* GetETokmoleConverter( reaction, TRUE );
		break;

	case kK0:
		*a = reaction->a * GetAConverter( reaction, FALSE ) 
			* GetATokmoleConverter( reaction, FALSE );
		*n = reaction->n;
		*e = reaction->e * GetEConverter( reaction, FALSE ) 
			* GetETokmoleConverter( reaction, FALSE );
		break;

	case kKinf:
		*a = reaction->aInf * GetAConverter( reaction, TRUE ) 
			* GetATokmoleConverter( reaction, TRUE );
		*n = reaction->nInf;
		*e = reaction->eInf * GetEConverter( reaction, TRUE ) 
			* GetETokmoleConverter( reaction, TRUE );		
		break;

	default:
		fprintf( stderr, "#error: unknown ReactionRateType in PrintArrCoeffsF77\n" );
		exit(2);
	}
}

void ConvertCDoubleToF77Double( Double c, char *f ) 
{
/*
	function gets a Double in c of the form [-]m.ddddddE±xx
	and returns a string in f of the form "[-]m.ddddddD±xx"
*/
	char	*pos;
	
	sprintf( f, "%.10E", c );
	pos = strrchr( f, 'E' );
	*pos = 'D';
}

int PrintMagicThirdBody( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ListPtr 		list = var;
	ThirdBodyPtr 	thirdBody = link->item;
	LinkPtr			speciesLink = NULL;
	SpeciesPtr	 	species = NULL;
	int 			i, countOuts = 0;

	fprintf( fp, "\tM[%s%s] = ", gPrefixThirdBody, CSymbol( thirdBody->name ) );
	for ( i = 0; i < thirdBody->speciesNumber->len; ++i ) {
	    speciesLink = Find_n_th_Item( thirdBody->speciesNumber->vec[i]+1, gSpeciesList );
	    species = ( SpeciesPtr ) speciesLink->item;
	    if ( thirdBody->speciesCoeff->vec[i] && !species->isSteadyState ) {
	      if ( countOuts % 4 == 0 && countOuts != 0  ) {
		fprintf( fp, "\n\t\t" );
		countOuts = 0;
	      }
	      ++countOuts;
	      if ( i ) {
		fprintf( fp, " + " );
	      }
	      if ( thirdBody->speciesCoeff->vec[i] != 1.0 ) {
		fprintf( fp, "%g * c[%s%s]", thirdBody->speciesCoeff->vec[i]
			 , gPrefixComputed, CSymbol( species->name ) );
	      }
	      else {
		fprintf( fp, "c[%s%s]", gPrefixComputed, CSymbol( species->name ) );
	      }
	    }
	}
	fprintf( fp, ";\n\n" );

	return 0;
}

int PrintMagicThirdBodyF77Vectorize( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ListPtr 		list = var;
	ThirdBodyPtr 	thirdBody = link->item;
	LinkPtr			speciesLink = NULL;
	SpeciesPtr	 	species = NULL;
	int 			i;
	int 			countOuts = 0;
	const int		nOfElementsLine = 2;
	const int		nOfElementsBlock = 20;

	
	for ( i = 0; i < thirdBody->speciesNumber->len; ++i ) {
	    speciesLink = Find_n_th_Item( thirdBody->speciesNumber->vec[i]+1, gSpeciesList );
	    species = ( SpeciesPtr ) speciesLink->item;
	    if ( thirdBody->speciesCoeff->vec[i] && !species->isSteadyState ) {
	      if ( countOuts % nOfElementsLine == 0 ) {
		if ( countOuts % nOfElementsBlock == 0 ) {
		  fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);
		  fprintf( fp, "      M(K,%s%s) = ", gPrefixThirdBodyF77, CSymbol( thirdBody->name ) );
		  if ( countOuts != 0 ) {
		    fprintf( fp, "M(K,%s%s)", gPrefixThirdBodyF77, CSymbol( thirdBody->name ) );
		  }
		}
		else {
		  fprintf( fp, "     &    " );
		}
	      }
	      if ( countOuts != 0 && i ) {
		fprintf( fp, " + " );
	      }
	      if ( thirdBody->speciesCoeff->vec[i] != 1.0 ) {
		fprintf( fp, "%g * C(K,%s%s)", thirdBody->speciesCoeff->vec[i]
			 , gPrefixComputedF77, CSymbol( species->name ) );
	      }
	      else {
		fprintf( fp, "C(K,%s%s)", gPrefixComputedF77, CSymbol( species->name ) );
	      }
	      if ( (countOuts+1) % nOfElementsLine == 0 ) {
		fprintf( fp, "\n" );
	      }
	      if ( (countOuts+1) % nOfElementsBlock == 0 ) {
		fprintf( fp, "%5d CONTINUE\n", gLabel);
		gLabel+=gLabelIncrement;
	      }
	      ++countOuts;
	    }
	}
	if ( (countOuts+1) % nOfElementsLine == 0 ) {
	  fprintf( fp, "\n" );
	}
	if ( countOuts % nOfElementsBlock != 0 ) {
	  fprintf( fp, "%5d CONTINUE\n", gLabel);
	  gLabel+=gLabelIncrement;
	}
	return 0;
}

int PrintMagicThirdBodyF77( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ListPtr 		list = var;
	ThirdBodyPtr 	thirdBody = link->item;
	LinkPtr			speciesLink = NULL;
	SpeciesPtr	 	species = NULL;
	int 			i;
	int 			countOuts = 0;
	const int		nOfElementsLine = 2;
	const int		nOfElementsBlock = 30;

	
	for ( i = 0; i < thirdBody->speciesNumber->len; ++i ) {
		speciesLink = Find_n_th_Item( thirdBody->speciesNumber->vec[i]+1, gSpeciesList );
		species = ( SpeciesPtr ) speciesLink->item;
		if ( thirdBody->speciesCoeff->vec[i] && !species->isSteadyState ) {
			if ( countOuts % nOfElementsLine == 0 ) {
				if ( countOuts % nOfElementsBlock == 0 ) {
					fprintf( fp, "\n      M(%s%s) = ", gPrefixThirdBodyF77, CSymbol( thirdBody->name ) );
					if ( countOuts != 0 ) {
						fprintf( fp, "M(%s%s)", gPrefixThirdBodyF77, CSymbol( thirdBody->name ) );
					}
				}
				else {
					fprintf( fp, "\n     &    " );
				}
			}
			if ( countOuts != 0 && i ) {
				fprintf( fp, " + " );
			}
			if ( thirdBody->speciesCoeff->vec[i] != 1.0 ) {
				fprintf( fp, "%g * C(%s%s)", thirdBody->speciesCoeff->vec[i]
					, gPrefixComputedF77, CSymbol( species->name ) );
			}
			else {
				fprintf( fp, "C(%s%s)", gPrefixComputedF77, CSymbol( species->name ) );
			}
			++countOuts;
		}
	}
	return 0;
}

int PrintMagicThirdBodyF90( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ListPtr 		list = var;
	ThirdBodyPtr 	thirdBody = link->item;
	LinkPtr			speciesLink = NULL;
	SpeciesPtr	 	species = NULL;
	int 			i;
	int 			countOuts = 0;
	const int		nOfElementsLine = 2;
	const int		nOfElementsBlock = 30;

	
	for ( i = 0; i < thirdBody->speciesNumber->len; ++i ) {
		speciesLink = Find_n_th_Item( thirdBody->speciesNumber->vec[i]+1, gSpeciesList );
		species = ( SpeciesPtr ) speciesLink->item;
		if ( thirdBody->speciesCoeff->vec[i] && !species->isSteadyState ) {
			if ( countOuts % nOfElementsLine == 0 ) {
				if ( countOuts % nOfElementsBlock == 0 ) {
					fprintf( fp, "\n      M(%s%s) = ", gPrefixThirdBodyF77, CSymbol( thirdBody->name ) );
					if ( countOuts != 0 ) {
						fprintf( fp, "M(%s%s)", gPrefixThirdBodyF77, CSymbol( thirdBody->name ) );
					}
				}
				else {
					fprintf( fp, " &\n\t   " );
				}
			}
			if ( countOuts != 0 && i ) {
				fprintf( fp, " + " );
			}
			if ( thirdBody->speciesCoeff->vec[i] != 1.0 ) {
				fprintf( fp, "%g * C(%s%s)", thirdBody->speciesCoeff->vec[i]
					, gPrefixComputedF77, CSymbol( species->name ) );
			}
			else {
				fprintf( fp, "C(%s%s)", gPrefixComputedF77, CSymbol( species->name ) );
			}
			++countOuts;
		}
	}
	return 0;
}

int PrintMagicReacRate( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	LinkPtr			speciesLink = NULL;
	LinkPtr			thirdBodyLink = NULL;
	SpeciesPtr 		species = NULL;
	int 			i, j;
	char			*thirdBody = NULL;

	fprintf( fp, "\tw[%s%s] = k[%s%s]", gPrefixReactions, CSymbol( reaction->label )
					, gPrefixReactions, CSymbol( reaction->label ) );

	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[i] );
		if ( !speciesLink ) FatalError( "something wrong in PrintMagicReacRate()" );
		species = ( SpeciesPtr ) speciesLink->item;

		if ( reaction->speciesCoeff[i] > 0.0 ) {
			for ( j = 0; j < reaction->speciesCoeff[i]; ++j ) {
				fprintf( fp, " * c[%s%s]"
						, gPrefixComputed, CSymbol( species->name ) );
			}
		}
	}

	if ( reaction->withThirdBody && !reaction->withLindemann) {
		if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
			thirdBody = ( char * ) thirdBodyLink->item;
			fprintf( fp, " * M[%s%s]", gPrefixThirdBody, CSymbol( thirdBody ) );
		}
		else {
			fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
			exit( 2 );
		}
	}

	fprintf( fp, ";\n" );
		
	return 0;
}

int PrintMagicReacRateF77Vectorize( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	LinkPtr			speciesLink = NULL;
	LinkPtr			thirdBodyLink = NULL;
	SpeciesPtr 		species = NULL;
	int 			i, j;
	char			*thirdBody = NULL;
	static int		reactionCount	=	0;
	const int		nOfReactions = 10;

	if ( reactionCount % nOfReactions == 0 ) {
		fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);
	}
	fprintf( fp, "      W(K,%s%s)=KR(K,%s%s)", gPrefixReactionsF77, CSymbol( reaction->label )
					, gPrefixReactionsF77, CSymbol( reaction->label ) );

	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[i] );
		if ( !speciesLink ) FatalError( "something wrong in PrintMagicReacRate()" );
		species = ( SpeciesPtr ) speciesLink->item;

		if ( reaction->speciesCoeff[i] > 0.0 ) {
			for ( j = 0; j < reaction->speciesCoeff[i]; ++j ) {
				fprintf( fp, "*C(K,%s%s)"
						, gPrefixComputedF77, CSymbol( species->name ) );
			}
		}
	}

/* gb 23.12.99 */
/*	if ( reaction->withThirdBody ) {*/
	if ( reaction->withThirdBody && !reaction->withLindemann ) {
		if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
			thirdBody = ( char * ) thirdBodyLink->item;
			fprintf( fp, " * M(K,%s%s)",gPrefixThirdBodyF77, CSymbol( thirdBody ) );
		}
		else {
			fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
			exit( 2 );
		}
	}

	fprintf( fp, "\n" );
	if ( ( reactionCount + 1 ) % nOfReactions == 0 || link->next == NULL ) {
		fprintf( fp, "%5d CONTINUE\n", gLabel);
		gLabel+=gLabelIncrement;
	}
	++reactionCount;
		
	return 0;
}

int PrintMagicReacRateF77( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	LinkPtr			speciesLink = NULL;
	LinkPtr			thirdBodyLink = NULL;
	SpeciesPtr 		species = NULL;
	int 			i, j;
	char			*thirdBody = NULL;

	fprintf( fp, "      W(%s%s) = K(%s%s)", gPrefixReactionsF77, CSymbol( reaction->label )
					, gPrefixReactionsF77, CSymbol( reaction->label ) );

	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[i] );
		if ( !speciesLink ) FatalError( "something wrong in PrintMagicReacRate()" );
		species = ( SpeciesPtr ) speciesLink->item;

		if ( reaction->speciesCoeff[i] > 0.0 ) {
			for ( j = 0; j < reaction->speciesCoeff[i]; ++j ) {
				fprintf( fp, " * C(%s%s)"
						, gPrefixComputedF77, CSymbol( species->name ) );
			}
		}
	}

	if ( reaction->withThirdBody && !reaction->withLindemann ) {
		if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
			thirdBody = ( char * ) thirdBodyLink->item;
			fprintf( fp, " * M(%s%s)",gPrefixThirdBodyF77, CSymbol( thirdBody ) );
		}
		else {
			fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
			exit( 2 );
		}
	}

	fprintf( fp, "\n" );
		
	return 0;
}

int PrintMagicProdRate( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	SpeciesPtr 		species = link->item;
	ReactionPtr 	reaction = NULL;
	LinkPtr			reactionLink = NULL;
	int 			j;
	int 			countOuts = 0;
	Flag			first = TRUE;

	if ( species->isSteadyState ) {
		return 0;
	}

	fprintf( fp, "\tcdot[%s%s] = ", gPrefixComputed, CSymbol( species->name ) );

	for ( reactionLink = gReactionList->firstItem; reactionLink != NULL
					; reactionLink = reactionLink->next ) {
		reaction = ( ReactionPtr ) reactionLink->item;

		for ( j = 0; j < reaction->numberOfSpecies; ++j ) {
			if ( reaction->speciesNumber[j] == species->number && reaction->speciesCoeff[j] != 0.0 ) {
				/* check and set counter */
				if ( countOuts % 4 == 0 && countOuts != 0  ) {
					fprintf( fp, "\n\t\t" );
					countOuts = 0;
				}
				++countOuts;

				/* print sign */
				if ( !first && reaction->speciesCoeff[j] < 0.0 ) {
					fprintf( fp, "%s", ( countOuts == 1 ) ? "+ " : " + " );
				}
				else if ( reaction->speciesCoeff[j] > 0.0 ) {
					fprintf( fp, "%s", ( countOuts == 1 ) ? "- " : " - " );
				}
				first = FALSE;
				
				/* print coefficient */
				if ( fabs( reaction->speciesCoeff[j] ) != 1.0 ) {
					fprintf( fp, "%g * ", fabs( reaction->speciesCoeff[j] ) );
				}
				
				/* print reaction rate */
				fprintf( fp, "w[%s%s]"
						, gPrefixReactions, CSymbol( reaction->label ) );
				
			}
		}
	}
	
	if ( first ) {
		fprintf( fp, "0.0;\n\n" );
	}
	else {
		fprintf( fp, ";\n\n" );
	}
		
	return 0;
}

int PrintMagicProdRateF77Vectorize( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	SpeciesPtr 		species = link->item;
	ReactionPtr 	reaction = NULL;
	LinkPtr			reactionLink = NULL;
	int 			j;
	int 			countOuts = 0;
	const int		nOfElementsLine = 2;
	const int		nOfElementsBlock = 20;

	if ( species->isSteadyState ) {
		return 0;
	}

	for ( reactionLink = gReactionList->firstItem; reactionLink != NULL
					; reactionLink = reactionLink->next ) {
		reaction = ( ReactionPtr ) reactionLink->item;

		for ( j = 0; j < reaction->numberOfSpecies; ++j ) {
			if ( reaction->speciesNumber[j] == species->number && reaction->speciesCoeff[j] != 0.0 ) {
				/* check and set counter */
				if ( countOuts % nOfElementsLine == 0 ) {
					if ( countOuts % nOfElementsBlock == 0 ) {
						fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);
						fprintf( fp, "      CDOT(K,%s%s)=", gPrefixComputedF77, CSymbol( species->name ) );
						if ( countOuts != 0 ) {
							fprintf( fp, "CDOT(K,%s%s)", gPrefixComputedF77, CSymbol( species->name ) );
						}
					}
					else {
						fprintf( fp, "     &    " );
					}
				}

				/* print sign */
				if ( countOuts != 0 && reaction->speciesCoeff[j] < 0.0 ) {
					fprintf( fp, "+" );
				}
				else if ( reaction->speciesCoeff[j] > 0.0 ) {
					fprintf( fp, "%s", ( countOuts == 0 ) ? "-" : "-" );
				}
				
				/* print coefficient */
				
				if ( fabs( reaction->speciesCoeff[j] ) != 1.0 ) {
					fprintf( fp, "%G*", fabs( reaction->speciesCoeff[j] ) );
				}
				
				/* print reaction rate */
				fprintf( fp, "W(K,%s%s)"
						,gPrefixReactionsF77, CSymbol( reaction->label ) );
					
				if ( (countOuts+1) % nOfElementsLine == 0 ) {
					fprintf( fp, "\n" );
				}
				if ( (countOuts+1) % nOfElementsBlock == 0 ) {
					fprintf( fp, "%5d CONTINUE\n", gLabel);
					gLabel+=gLabelIncrement;
				}
				++countOuts;
			}
		}
	}
	if ( (countOuts+1) % nOfElementsLine == 0 ) {
		fprintf( fp, "\n" );
	}
	if ( (countOuts) % nOfElementsBlock != 0 ) {
		fprintf( fp, "%5d CONTINUE\n", gLabel);
		gLabel+=gLabelIncrement;
	}
	
	if ( countOuts == 0 ) {
		fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);
		fprintf( fp, "      CDOT(K,%s%s) = 0.0D+00", gPrefixComputedF77, CSymbol( species->name ) );
		fprintf( fp, "\n%5d CONTINUE\n", gLabel);
		gLabel+=gLabelIncrement;
	}
	return 0;
}

int PrintMagicProdRateF77( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	SpeciesPtr 		species = link->item;
	ReactionPtr 	reaction = NULL;
	LinkPtr			reactionLink = NULL;
	int 			j;
	int 			countOuts = 0;
	const int		nOfElementsLine = 3;
	const int		nOfElementsBlock = 45;

	if ( species->isSteadyState ) {
		return 0;
	}

	for ( reactionLink = gReactionList->firstItem; reactionLink != NULL
					; reactionLink = reactionLink->next ) {
		reaction = ( ReactionPtr ) reactionLink->item;

		for ( j = 0; j < reaction->numberOfSpecies; ++j ) {
			if ( reaction->speciesNumber[j] == species->number ) {
				/* check and set counter */
				if ( countOuts % nOfElementsLine == 0 ) {
					if ( countOuts % nOfElementsBlock == 0 ) {
						fprintf( fp, "\n      CDOT(%s%s) = ", gPrefixComputedF77, CSymbol( species->name ) );
						if ( countOuts != 0 ) {
							fprintf( fp, "CDOT(%s%s)", gPrefixComputedF77, CSymbol( species->name ) );
						}
					}
					else {
						fprintf( fp, "\n     &    " );
					}
				}

				/* print sign */
				if ( countOuts != 0 && reaction->speciesCoeff[j] < 0.0 ) {
					fprintf( fp, " + " );
				}
				else if ( reaction->speciesCoeff[j] > 0.0 ) {
					fprintf( fp, "%s", ( countOuts == 0 ) ? "- " : " - " );
				}
				
				/* print coefficient */
				
				if ( reaction->speciesCoeff[j] != 0.0 ) {
					if ( fabs( reaction->speciesCoeff[j] ) != 1.0 ) {
						fprintf( fp, "%G * ", fabs( reaction->speciesCoeff[j] ) );
					}
				
					/* print reaction rate */
					fprintf( fp, "W(%s%s)"
							,gPrefixReactionsF77, CSymbol( reaction->label ) );
					
					++countOuts;
				}
			}
		}
	}
	
	if ( countOuts == 0 ) {
		fprintf( fp, "\n      CDOT(%s%s) = 0.0D0", gPrefixComputedF77, CSymbol( species->name ) );
	}
	return 0;
}

int PrintMagicProdRateF90( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	SpeciesPtr 		species = link->item;
	ReactionPtr 	reaction = NULL;
	LinkPtr			reactionLink = NULL;
	int 			j;
	int 			countOuts = 0;
	const int		nOfElementsLine = 3;
	const int		nOfElementsBlock = 45;

	if ( species->isSteadyState ) {
		return 0;
	}

	for ( reactionLink = gReactionList->firstItem; reactionLink != NULL
					; reactionLink = reactionLink->next ) {
		reaction = ( ReactionPtr ) reactionLink->item;

		for ( j = 0; j < reaction->numberOfSpecies; ++j ) {
			if ( reaction->speciesNumber[j] == species->number ) {
				/* check and set counter */
				if ( countOuts % nOfElementsLine == 0 ) {
					if ( countOuts % nOfElementsBlock == 0 ) {
						fprintf( fp, "\n      CDOT(%s%s) = ", gPrefixComputedF77, CSymbol( species->name ) );
						if ( countOuts != 0 ) {
							fprintf( fp, "CDOT(%s%s)", gPrefixComputedF77, CSymbol( species->name ) );
						}
					}
					else {
						fprintf( fp, " & \n\t   " );
					}
				}

				/* print sign */
				if ( countOuts != 0 && reaction->speciesCoeff[j] < 0.0 ) {
					fprintf( fp, " + " );
				}
				else if ( reaction->speciesCoeff[j] > 0.0 ) {
					fprintf( fp, "%s", ( countOuts == 0 ) ? "- " : " - " );
				}
				
				/* print coefficient */
				
				if ( reaction->speciesCoeff[j] != 0.0 ) {
					if ( fabs( reaction->speciesCoeff[j] ) != 1.0 ) {
						fprintf( fp, "%G * ", fabs( reaction->speciesCoeff[j] ) );
					}
				
					/* print reaction rate */
					fprintf( fp, "W(%s%s)"
							,gPrefixReactionsF77, CSymbol( reaction->label ) );
					
					++countOuts;
				}
			}
		}
	}
	
	if ( countOuts == 0 ) {
		fprintf( fp, "\n      CDOT(%s%s) = 0.0_DP", gPrefixComputedF77, CSymbol( species->name ) );
	}
	return 0;
}

int PrintMagicMolarMass( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	SpeciesPtr 		species = link->item;

	fprintf( fp, "\tW[%s%s] = %15.8e;\n"
			, gPrefixComputed, CSymbol( species->name ), species->molarMass );

	return 0;
}

int PrintMagicMolarMassF77( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	SpeciesPtr 		species = link->item;

	fprintf( fp, "      MM(%s%s) = %15.8e\n"
			, gPrefixComputedF77, CSymbol( species->name ), species->molarMass );

	return 0;
}

int PrintMagicNames( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	SpeciesPtr 		species = link->item;

	fprintf( fp, "\tstrcpy( names[%s%s], \"%-20s\" );\n"
			, gPrefixComputed, CSymbol( species->name ), species->name );

	return 0;
}

int PrintMagicNamesF77( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	SpeciesPtr 		species = link->item;

	fprintf( fp, "      NAMES(%s%s)='%-20s'\n"
			, gPrefixComputedF77, CSymbol( species->name ), species->name );

	return 0;
}

int PrintOneThermoDataF77( LinkPtr link, void *var )
{
	struct THData {
		FILE	*fp;
		Flag	hot;
	} *thermoData = var;
	
	FILE 			*fp = thermoData->fp;
	SpeciesPtr 		species = link->item;
	Double			*a = ( thermoData->hot ) ? species->coeffHot : species->coeffCold;
	int 			j;

	if ( species->isSteadyState ) {
		return 0;
	}

	fprintf( fp, "      H(%s%s) = ", gPrefixComputedF77, CSymbol( species->name ) );
	fprintf( fp, "%15.8e * (\n"
					, RGAS / species->molarMass );
	fprintf( fp, "     &      T * ( %15.8e + T * ( %15.8e\n"
					, a[0], 0.5 * a[1] );
	fprintf( fp, "     &      + T * ( %15.8e + T * ( %15.8e\n"
					, a[2] / 3.0, 0.25 * a[3] );
	fprintf( fp, "     &      + T * ( %15.8e ) ) ) ) )%s %15.8e )\n"
					, 0.2 * a[4], (a[5]<0.0)?"":" +", a[5] );
	fprintf( fp, "      CP(%s%s) = ", gPrefixComputedF77, CSymbol( species->name ) );
	fprintf( fp, "%15.8e * (\n"
					, RGAS / species->molarMass );
	fprintf( fp, "     &      %15.8e + T * ( %15.8e \n"
					, a[0], a[1] );
	fprintf( fp, "     &      + T * ( %15.8e + T * ( %15.8e\n"
					, a[2], a[3] );
	fprintf( fp, "     &      + T * ( %15.8e ) ) ) ) )\n"
					, a[4] );

	return 0;
}

int PrintMagicThermoDataF77( FILE *fp )
{
	struct THData {
		FILE	*fp;
		Flag	hot;
	} thermoData;

	fprintf( fp, "      SUBROUTINE COMPTHERMODATA( H, CP, T )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS\n" );
	fprintf( fp, "C     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES\n" );
	fprintf( fp, "C     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.\n" );
	fprintf( fp, "C     THE PARAMETER H AND CP SHOULD PROVIDE WORKSPACE OF LENGTH %d\n"
					, gCounter->species - gCounter->steadyStates );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "      DOUBLE PRECISION H(%d), CP(%d), T\n",
		gCounter->species - gCounter->steadyStates, gCounter->species - gCounter->steadyStates );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IF (T.GT.1.0D3) THEN\n" );
	thermoData.fp = fp;
	thermoData.hot = TRUE;
	ListIterator( gSpeciesList, PrintOneThermoDataF77, &thermoData );
	fprintf( fp, "      ELSE\n" );
	thermoData.hot = FALSE;
	ListIterator( gSpeciesList, PrintOneThermoDataF77, &thermoData );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	
	return 0;
}

int PrintOneThermoDataF90( LinkPtr link, void *var )
{
	struct THData {
		FILE	*fp;
		Flag	hot;
	} *thermoData = var;
	
	FILE 			*fp = thermoData->fp;
	SpeciesPtr 		species = link->item;
	Double			*a = ( thermoData->hot ) ? species->coeffHot : species->coeffCold;
	int 			j;

	if ( species->isSteadyState ) {
		return 0;
	}

	fprintf( fp, "      H(%s%s) = ", gPrefixComputedF77, CSymbol( species->name ) );
	fprintf( fp, "%15.8e_DP * ( &\n"
					, RGAS / species->molarMass );
	fprintf( fp, "\t   T * ( %15.8e_DP + T * ( %15.8e_DP &\n"
					, a[0], 0.5 * a[1] );
	fprintf( fp, "\t   + T * ( %15.8e_DP + T * ( %15.8e_DP &\n"
					, a[2] / 3.0, 0.25 * a[3] );
	fprintf( fp, "\t   + T * ( %15.8e_DP ) ) ) ) )%s %15.8e_DP )\n"
					, 0.2 * a[4], (a[5]<0.0)?"":" +", a[5] );
	fprintf( fp, "      CP(%s%s) = ", gPrefixComputedF77, CSymbol( species->name ) );
	fprintf( fp, "%15.8e_DP * ( &\n"
					, RGAS / species->molarMass );
	fprintf( fp, "\t   %15.8e_DP + T * ( %15.8e_DP &\n"
					, a[0], a[1] );
	fprintf( fp, "\t   + T * ( %15.8e_DP + T * ( %15.8e_DP &\n"
					, a[2], a[3] );
	fprintf( fp, "\t   + T * ( %15.8e_DP ) ) ) ) )\n"
					, a[4] );

	return 0;
}

int PrintMagicThermoDataF90( FILE *fp )
{
	struct THData {
		FILE	*fp;
		Flag	hot;
	} thermoData;

	fprintf( fp, "      SUBROUTINE COMPTHERMODATA( H, CP, T )\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "!     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS\n" );
	fprintf( fp, "!     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES\n" );
	fprintf( fp, "!     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.\n" );
	fprintf( fp, "!     THE PARAMETER H AND CP SHOULD PROVIDE WORKSPACE OF LENGTH %d\n"
					, gCounter->species - gCounter->steadyStates );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "      implicit none\n" );
	fprintf( fp, "      include '%sF90.h'\n", gOptions->base );
	fprintf( fp, "      real(DP) :: H(%d), CP(%d), T\n",
		gCounter->species - gCounter->steadyStates, gCounter->species - gCounter->steadyStates );
	fprintf( fp, "\n" );
	fprintf( fp, "      IF (T.GT.1000.0_DP) THEN\n" );
	thermoData.fp = fp;
	thermoData.hot = TRUE;
	ListIterator( gSpeciesList, PrintOneThermoDataF90, &thermoData );
	fprintf( fp, "      ELSE\n" );
	thermoData.hot = FALSE;
	ListIterator( gSpeciesList, PrintOneThermoDataF90, &thermoData );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      END\n" );
	
	return 0;
}

int PrintMagicThermoData( FILE *fp )
{
	struct THData {
		FILE	*fp;
		Flag	hot;
	} thermoData;


	fprintf( fp, "void ComputeThermoData( double *h, double *cp, double T )\n" );
	fprintf( fp, "{\n" );
	fprintf( fp, "/*\n\tThis function computes enthalpy 'h' and heat capacity 'cp' as\n" );
	fprintf( fp, "\tfunction of temperature 'T' for all non steady state species\n" );
	fprintf( fp, "\tin units [J/kg] and [J/kg K], respectively.\n" );
	fprintf( fp, "\tThe parameter h and cp should provide workspace of length %d */\n\n\n"
					, gCounter->species - gCounter->steadyStates );
	fprintf( fp, "\tint i;\n" );

	fprintf( fp, "\tif ( T > 1000.0 ) {\n" );
	thermoData.fp = fp;
	thermoData.hot = TRUE;
	ListIterator( gSpeciesList, PrintOneThermoData, &thermoData );

	fprintf( fp, "\t}\n\telse if (T >= 299.999999 ) {\n" );
	thermoData.hot = FALSE;
	ListIterator( gSpeciesList, PrintOneThermoData, &thermoData );
	fprintf( fp, "\t}\n" );

	fprintf( fp, "\telse {\n" );
	fprintf( fp, "\t\tComputeThermoData( h, cp, 300.0 );\n" );
	fprintf( fp, "\t\tfor (i = 0; i < sEnd; i++) {\n" );
	fprintf( fp, "\t\t\th[i] = (T-300.)*cp[i] + h[i];\n" );
	fprintf( fp, "\t\t}\n" );
	fprintf( fp, "\t}\n" );
	fprintf( fp, "}\n" );
	
	return 0;
}

int PrintOneThermoData( LinkPtr link, void *var )
{
	struct THData {
		FILE	*fp;
		Flag	hot;
	} *thermoData = var;
	
	FILE 			*fp = thermoData->fp;
	SpeciesPtr 		species = link->item;
	Double			*a = ( thermoData->hot ) ? species->coeffHot : species->coeffCold;
	int 			j;

	if ( species->isSteadyState ) {
		return 0;
	}

	fprintf( fp, "\t\th[%s%s] = ", gPrefixComputed, CSymbol( species->name ) );
	fprintf( fp, "%15.8e * (\n"
					, RGAS / species->molarMass );
	fprintf( fp, "\t\t\tT * ( %15.8e + T * ( %15.8e\n"
					, a[0], 0.5 * a[1] );
	fprintf( fp, "\t\t\t+ T * ( %15.8e + T * ( %15.8e\n"
					, a[2] / 3.0, 0.25 * a[3] );
	fprintf( fp, "\t\t\t+ T * %15.8e ) ) ) )%s %15.8e );\n"
					, 0.2 * a[4], (a[5]<0.0)?"":" +", a[5] );
	fprintf( fp, "\t\tcp[%s%s] = ", gPrefixComputed, CSymbol( species->name ) );
	fprintf( fp, "%15.8e * (\n"
					, RGAS / species->molarMass );
	fprintf( fp, "\t\t\t%15.8e + T * ( %15.8e \n"
					, a[0], a[1] );
	fprintf( fp, "\t\t\t+ T * ( %15.8e + T * ( %15.8e\n"
					, a[2], a[3] );
	fprintf( fp, "\t\t\t+ T * %15.8e ) ) ) );\n"
					, a[4] );

	return 0;
}

int ReactionLabelsToUpper( LinkPtr link, void *aux )
{
	ReactionPtr 	reaction = link->item;
	UpperString( reaction->label );
	
	return 0;
}

int PrintOneThermoDataF77Vectorize( LinkPtr link, void *var )
{
	struct THData {
		FILE	*fp;
		Flag	hot;
	} *thermoData = var;
	
	FILE 			*fp = thermoData->fp;
	SpeciesPtr 		species = link->item;
	Double			*a = ( thermoData->hot ) ? species->coeffHot : species->coeffCold;
	char			*c = ( thermoData->hot ) ? "CHIGH" : "CLOW";
	int 			j;

	if ( species->isSteadyState ) {
		return 0;
	}

	fprintf( fp, "      %s(1,%s%s) = %15.8e\n", c, gPrefixComputedF77, CSymbol( species->name ), a[0] );
	fprintf( fp, "      %s(2,%s%s) = %15.8e\n", c, gPrefixComputedF77, CSymbol( species->name ), a[1] );
	fprintf( fp, "      %s(3,%s%s) = %15.8e\n", c, gPrefixComputedF77, CSymbol( species->name ), a[2] );
	fprintf( fp, "      %s(4,%s%s) = %15.8e\n", c, gPrefixComputedF77, CSymbol( species->name ), a[3] );
	fprintf( fp, "      %s(5,%s%s) = %15.8e\n", c, gPrefixComputedF77, CSymbol( species->name ), a[4] );
	fprintf( fp, "      %s(6,%s%s) = %15.8e\n", c, gPrefixComputedF77, CSymbol( species->name ), a[5] );
	fprintf( fp, "      %s(7,%s%s) = %15.8e\n", c, gPrefixComputedF77, CSymbol( species->name ), a[6] );
	fprintf( fp, "C\n" );

	return 0;
}

int PrintMagicNASAPolynomialsVectorize( FILE *fp )
{
	struct THData {
		FILE	*fp;
		Flag	hot;
	} thermoData;

	fprintf( fp, "      SUBROUTINE GETNASAPOLYNOMIALS(CLOW,CHIGH,NCOEF,NSPECS)\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     THIS SUBROUTINE RETURNS THE NASA-POLYNOMIALS OF ALL NON STEADY\n" );
	fprintf( fp, "C     STATE SPECIES. CLOW AND CHIGH CONTAIN THE COEFFICIENTS FOR \n" );
	fprintf( fp, "C     TEMPERATURES LOWER THAN 1000 K, OR HIGHER THAN 1000 K, \n" );
	fprintf( fp, "C     RESPECTIVELY.\n" );
	fprintf( fp, "C     IF A SPECIES IS IRRELEVANT FOR THE THERMO-BALANCE, ALL THE \n" );
	fprintf( fp, "C     COEFFICIENTS OF THIS SPECIES ARE EQUAL TO ZERO.\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "      INTEGER NCOEF,NSPECS\n" );
	fprintf( fp, "      DOUBLE PRECISION CLOW(7,%d), CHIGH(7,%d)\n",
		gCounter->species - gCounter->steadyStates, gCounter->species - gCounter->steadyStates );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IF (NCOEF.NE.7.OR.NSPECS.NE.%d) THEN\n", gCounter->species - gCounter->steadyStates );
	fprintf( fp, "        WRITE(*,*) '### ERROR in GETNASAPOLYNOMIALS:'\n" );
	fprintf( fp, "        WRITE(*,*) '    NCOEF  must be 7'\n" );
	fprintf( fp, "        WRITE(*,*) '    NSPECS must be %d'\n", gCounter->species - gCounter->steadyStates );
	fprintf( fp, "        STOP\n" );
	fprintf( fp, "      ENDIF\n" );
	fprintf( fp, "C\n" );
	thermoData.fp = fp;
	thermoData.hot = FALSE;
	ListIterator( gSpeciesList, PrintOneThermoDataF77Vectorize, &thermoData );
	fprintf( fp, "C\n" );
	thermoData.hot = TRUE;
	ListIterator( gSpeciesList, PrintOneThermoDataF77Vectorize, &thermoData );
	fprintf( fp, "      END\n" );
	
	return 0;
}

/*
int PrintMagicThermoDataF77Vectorize( FILE *fp )
{
	struct THData {
		FILE	*fp;
		Flag	hot;
	} thermoData;
	int label1,label2;
	gLabel=10;
	
	label1=gLabel;
	gLabel+=10;
	label2=gLabel;
	gLabel+=10;

	fprintf( fp, "      SUBROUTINE COMPTHERMODATA(H,CP,TEMP,CLOW,CHIGH,MOLMASS,NCOEF\n" );
	fprintf( fp, "     $                          ,NSPECS,NGRIDPOINTS,RGAS)\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS\n" );
	fprintf( fp, "C     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES\n" );
	fprintf( fp, "C     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INTEGER NCOEF,NGRIDPOINTS,NSPECS\n" );
	fprintf( fp, "      DOUBLE PRECISION TEMP(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION H(0:(NGRIDPOINTS+1),NSPECS)\n");
	fprintf( fp, "      DOUBLE PRECISION CP(0:(NGRIDPOINTS+1),NSPECS)\n");
	fprintf( fp, "      DOUBLE PRECISION MOLMASS(NSPECS)\n");
	fprintf( fp, "      DOUBLE PRECISION CLOW(NCOEF,NSPECS)\n");
	fprintf( fp, "      DOUBLE PRECISION CHIGH(NCOEF,NSPECS)\n");
	fprintf( fp, "      DOUBLE PRECISION RGAS\n");
	fprintf( fp, "C\n" );
	fprintf( fp, "      INTEGER I,K\n" );
	fprintf( fp, "      DOUBLE PRECISION C1,C2,C3,C4,C5,C6,RRMW,T1\n");
	fprintf( fp, "      DOUBLE PRECISION TLOW, THIGH\n");
	fprintf( fp, "      PARAMETER (TLOW=300.,THIGH=3000.)\n");
	fprintf( fp, "C\n" );
	fprintf( fp, "      DO %d I=1,NSPECS\n", label1);
	fprintf( fp, "        RRMW=RGAS/MOLMASS(i)\n");
	fprintf( fp, "        DO %d K=1,NGRIDPOINTS\n", label2);
	fprintf( fp, "          T1=TEMP(K)\n" );
	fprintf( fp, "          IF (T1.GT.1.0D+03) THEN\n" );
	fprintf( fp, "            C1=CHIGH(1,I)\n" );
	fprintf( fp, "            C2=CHIGH(2,I)\n" );
	fprintf( fp, "            C3=CHIGH(3,I)\n" );
	fprintf( fp, "            C4=CHIGH(4,I)\n" );
	fprintf( fp, "            C5=CHIGH(5,I)\n" );
	fprintf( fp, "            C6=CHIGH(6,I)\n" );
	fprintf( fp, "          ELSE\n" );
	fprintf( fp, "            C1=CLOW(1,I)\n" );
	fprintf( fp, "            C2=CLOW(2,I)\n" );
	fprintf( fp, "            C3=CLOW(3,I)\n" );
	fprintf( fp, "            C4=CLOW(4,I)\n" );
	fprintf( fp, "            C5=CLOW(5,I)\n" );
	fprintf( fp, "            C6=CLOW(6,I)\n" );
	fprintf( fp, "          END IF\n" );
	fprintf( fp, "          H(K,I)=((C1+(C2/2.0D+00+(C3/3.0D+00+(C4/4.0D+00\n" );
	fprintf( fp, "     $                +C5/5.0D+00*T1)*T1)*t1)*T1)*T1+C6)*RRMW\n" );
	fprintf( fp, "          CP(K,I)=(C1+(C2+(C3+(C4+C5*T1)*T1)*T1)*T1)*RRMW\n" );
	fprintf( fp, "%5d   CONTINUE\n", label2);
	fprintf( fp, "%5d CONTINUE\n", label1);
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	
	return 0;
}
*/
int PrintMagicThermoDataF77Vectorize( FILE *fp )
{
	struct THData {
		FILE	*fp;
		Flag	hot;
	} thermoData;
	int label1,label2,label3,label4,label5;
	gLabel=10;
	
	label1=gLabel;
	gLabel+=10;
	label2=gLabel;
	gLabel+=10;
	label3=gLabel;
	gLabel+=10;
	label4=gLabel;
	gLabel+=10;
	label5=gLabel;
	gLabel+=10;

	fprintf( fp, "      SUBROUTINE COMPTHERMODATA(H,CP,TEMP,CLOW,CHIGH,MOLMASS,NCOEF\n" );
	fprintf( fp, "     $                          ,NSPECS,NGRIDPOINTS,RGAS)\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     THIS FUNCTION COMPUTES ENTHALPY 'H' AND HEAT CAPACITY 'CP' AS\n" );
	fprintf( fp, "C     FUNCTION OF TEMPERATURE T FOR ALL NON STEADY STATE SPECIES\n" );
	fprintf( fp, "C     IN UNITS [J/KG] and [J/KG K], RESPECTIVELY.\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INTEGER NCOEF,NGRIDPOINTS,NSPECS\n" );
	fprintf( fp, "      DOUBLE PRECISION TEMP(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION H(0:(NGRIDPOINTS+1),NSPECS)\n");
	fprintf( fp, "      DOUBLE PRECISION CP(0:(NGRIDPOINTS+1),NSPECS)\n");
	fprintf( fp, "      DOUBLE PRECISION MOLMASS(NSPECS)\n");
	fprintf( fp, "      DOUBLE PRECISION CLOW(NCOEF,NSPECS)\n");
	fprintf( fp, "      DOUBLE PRECISION CHIGH(NCOEF,NSPECS)\n");
	fprintf( fp, "      DOUBLE PRECISION RGAS\n");
	fprintf( fp, "C\n" );
	fprintf( fp, "      INTEGER I,K\n" );
	fprintf( fp, "      DOUBLE PRECISION C1,C2,C3,C4,C5,C6,RRMW,T1\n");
	fprintf( fp, "      DOUBLE PRECISION TLOW, THIGH\n");
	fprintf( fp, "      PARAMETER (TLOW=300.,THIGH=3000.)\n");
	fprintf( fp, "C\n" );
	fprintf( fp, "      DO %d K=0,NGRIDPOINTS+1\n", label5 );
	fprintf( fp, "        T1=TEMP(K)\n" );
	fprintf( fp, "        IF(T1.LE.TLOW) THEN\n" );
	fprintf( fp, "          T1=TLOW\n" );
	fprintf( fp, "          DO %d I=1,NSPECS\n", label1 );
	fprintf( fp, "            RRMW=RGAS/MOLMASS(I)\n" );
	fprintf( fp, "            C1=CLOW(1,I)\n" );
	fprintf( fp, "            C2=CLOW(2,I)\n" );
	fprintf( fp, "            C3=CLOW(3,I)\n" );
	fprintf( fp, "            C4=CLOW(4,I)\n" );
	fprintf( fp, "            C5=CLOW(5,I)\n" );
	fprintf( fp, "            C6=CLOW(6,I)\n" );
	fprintf( fp, "            H(K,I)=((C1+(C2/2.0D+00+(C3/3.0D+00+(C4/4.0D+00\n" );
	fprintf( fp, "     $                 +C5/5.0D+00*T1)*T1)*t1)*T1)*T1+C6)*RRMW\n" );
	fprintf( fp, "            CP(K,I)=(C1+(C2+(C3+(C4+C5*T1)*T1)*T1)*T1)*RRMW\n" );
	fprintf( fp, "            H(K,I)=H(K,I)+CP(K,I)*(Temp(k)-TLOW)\n" );
	fprintf( fp, "%5d       CONTINUE\n", label1 );
	fprintf( fp, "        ELSE IF(T1.LE.THIGH) THEN\n" );
	fprintf( fp, "          IF(T1.LT.1.0D+03) THEN\n" );
	fprintf( fp, "            DO %d I=1,NSPECS\n", label2 );
	fprintf( fp, "              RRMW=RGAS/MOLMASS(I)\n" );
	fprintf( fp, "              C1=CLOW(1,I)\n" );
	fprintf( fp, "              C2=CLOW(2,I)\n" );
	fprintf( fp, "              C3=CLOW(3,I)\n" );
	fprintf( fp, "              C4=CLOW(4,I)\n" );
	fprintf( fp, "              C5=CLOW(5,I)\n" );
	fprintf( fp, "            C6=CLOW(6,I)\n" );
	fprintf( fp, "            H(K,I)=((C1+(C2/2.0D+00+(C3/3.0D+00+(C4/4.0D+00\n" );
	fprintf( fp, "     $                 +C5/5.0D+00*T1)*T1)*t1)*T1)*T1+C6)*RRMW\n" );
	fprintf( fp, "            CP(K,I)=(C1+(C2+(C3+(C4+C5*T1)*T1)*T1)*T1)*RRMW\n" );
	fprintf( fp, "%5d        CONTINUE\n", label2 );
	fprintf( fp, "        ELSE\n" );
	fprintf( fp, "          DO %d I=1,NSPECS\n", label3 );
	fprintf( fp, "            RRMW=RGAS/MOLMASS(I)\n" );
	fprintf( fp, "            C1=CHIGH(1,I)\n" );
	fprintf( fp, "            C2=CHIGH(2,I)\n" );
	fprintf( fp, "            C3=CHIGH(3,I)\n" );
	fprintf( fp, "            C4=CHIGH(4,I)\n" );
	fprintf( fp, "            C5=CHIGH(5,I)\n" );
	fprintf( fp, "            C6=CHIGH(6,I)\n" );
	fprintf( fp, "            H(K,I)=((C1+(C2/2.0D+00+(C3/3.0D+00+(C4/4.0D+00\n" );
	fprintf( fp, "     $                 +C5/5.0D+00*T1)*T1)*t1)*T1)*T1+C6)*RRMW\n" );
	fprintf( fp, "            CP(K,I)=(C1+(C2+(C3+(C4+C5*T1)*T1)*T1)*T1)*RRMW\n" );
	fprintf( fp, "%5d       CONTINUE\n", label3 );
	fprintf( fp, "        ENDIF\n" );
	fprintf( fp, "      ELSE\n" );
	fprintf( fp, "        T1=THIGH\n" );
	fprintf( fp, "        DO %d I=1,NSPECS\n", label4 );
	fprintf( fp, "          RRMW=RGAS/MOLMASS(I)\n" );
	fprintf( fp, "          C1=CHIGH(1,I)\n" );
	fprintf( fp, "          C2=CHIGH(2,I)\n" );
	fprintf( fp, "          C3=CHIGH(3,I)\n" );
	fprintf( fp, "        C4=CHIGH(4,I)\n" );
	fprintf( fp, "        C5=CHIGH(5,I)\n" );
	fprintf( fp, "        C6=CHIGH(6,I)\n" );
	fprintf( fp, "        H(K,I)=((C1+(C2/2.0D+00+(C3/3.0D+00+(C4/4.0D+00\n" );
	fprintf( fp, "     $             +C5/5.0D+00*T1)*T1)*t1)*T1)*T1+C6)*RRMW\n" );
	fprintf( fp, "        CP(K,I)=(C1+(C2+(C3+(C4+C5*T1)*T1)*T1)*T1)*RRMW\n" );
	fprintf( fp, "        H(K,I)=H(K,I)+CP(K,I)*(Temp(k)-THIGH)\n" );
	fprintf( fp, "%5d   CONTINUE\n", label4 );
	fprintf( fp, "      ENDIF\n" );
	fprintf( fp, "%5d   CONTINUE\n", label5 );
	fprintf( fp, "      END\n" );

	return 0;
}

int PrintRedMechProdRate( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	SpeciesPtr 		species = link->item;
	ReactionPtr 	reaction = NULL;
	LinkPtr			reactionLink = NULL;
	int 			j;
	int 			countOuts = 0;
	Flag			first = TRUE;

	fprintf( fp, "L[%s%s] = ", gPrefixComputed, CSymbol( species->name ) );

	for ( reactionLink = gReactionList->firstItem; reactionLink != NULL
					; reactionLink = reactionLink->next ) {
		reaction = ( ReactionPtr ) reactionLink->item;
		for ( j = 0; j < reaction->numberOfSpecies; ++j ) {
			if ( reaction->speciesNumber[j] == species->number ) {
				/* print sign */
				if ( !first && reaction->speciesCoeff[j] < 0.0 ) {
					fprintf( fp, "%s", " + " );
				}
				else if ( reaction->speciesCoeff[j] > 0.0 ) {
					fprintf( fp, "%s", ( first ) ? "- " : " - " );
				}
				first = FALSE;
				
				/* print coefficient */
				if ( fabs( reaction->speciesCoeff[j] ) != 1.0 ) {
					fprintf( fp, "%g * ", fabs( reaction->speciesCoeff[j] ) );
				}
				fprintf( fp, "w[%s%s]"
						, gPrefixReactions, CSymbol( reaction->label ) );
			}
		}
	}
	
	if ( first ) {
		fprintf( fp, ".;\n" );
	}
	else {
		fprintf( fp, ";\n" );
	}
		
	return 0;
}

int PrintRedMechReaction( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	LinkPtr			speciesLink = NULL;
	LinkPtr			thirdBodyLink = NULL;
	SpeciesPtr 		species = NULL;
	int 			i, j;
	char			*thirdBody = NULL;

	fprintf( fp, "w[%s%s] = k[%s%s]", gPrefixReactions, CSymbol( reaction->label )
					, gPrefixReactions, CSymbol( reaction->label ));

/*  print left side of equation  */
	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[i] );
		if ( !speciesLink ) FatalError( "something wrong in PrintReaction()" );
		species = ( SpeciesPtr ) speciesLink->item;

		if ( reaction->speciesCoeff[i] > 0.0 ) {
			for ( j = 0; j < reaction->speciesCoeff[i]; ++j ) {
				fprintf( fp, " * c[%s%s]"
						, gPrefixComputed, CSymbol( species->name ) );
			}
		}
	}

	if ( reaction->withThirdBody ) {
		if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
			thirdBody = ( char * ) thirdBodyLink->item;
			fprintf( fp, " * M[%s%s]", gPrefixThirdBody, CSymbol( thirdBody ) );
		}
		else {
			fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
			exit( 2 );
		}
	}

	fprintf( fp, ";\n" );
		
	return 0;
}

int PrintThirdBody( LinkPtr link, void *var )
{
	ListPtr 		list = var;
	ThirdBodyPtr 	thirdBody = link->item;
	LinkPtr			speciesLink = NULL;
	SpeciesPtr	 	species = NULL;
	int 			i, countOuts = 0;

	if ( link == list->firstItem ) {
		ListInfo( list, stdout );	
		fprintf( stdout, "The different third body efficiencies are:\n");
	}

	fprintf( stdout, "Nr.%d: [%s] = ", thirdBody->id, thirdBody->name );
	for ( i = 0; i < thirdBody->speciesNumber->len; ++i ) {
		if ( thirdBody->speciesNumber->vec[i] != -1 && thirdBody->speciesCoeff->vec[i] ) {
		  speciesLink = Find_n_th_Item( thirdBody->speciesNumber->vec[i]+1, gSpeciesList );
		  species = ( SpeciesPtr ) speciesLink->item;
			if ( !i ) {
				fprintf( stdout, "%g[%s] ", thirdBody->speciesCoeff->vec[i], species->name );
			}
			else {
				fprintf( stdout, "+ %g[%s] ", thirdBody->speciesCoeff->vec[i], species->name );
				if ( (countOuts+1) % 4 == 0 ) {
					fprintf( stdout, "\n\t\t\t\t" );
				}
			}
			++countOuts;
		}
	}
	fprintf( stdout, "\n\n" );

	return 0;
}

int PrintDimensionTable( LinkPtr link, void *var )
{
	ListPtr 		list = var;
	DimensionPtr 	dimension = link->item;

	if ( link == list->firstItem ) {
		ListInfo( list, stdout );	
		fprintf( stdout, "Allowed Dimensions are:\n\n");
		fprintf( stdout, "%-10s%-13s%-4s%-4s%-4s%-4s%-4s\n", "", "value", "kg","m","s", "K", "mole" );
	}

	fprintf( stdout, "%-10s%-13g%-4d%-4d%-4d%-4d%-4d\n", dimension->name, dimension->value, dimension->kg, dimension->m, dimension->s, dimension->K, dimension->mole );

	return 0;
}

int PrintDimension( LinkPtr link, void *var )
{
	ListPtr 		list = var;
	DimensionPtr 	dimension = link->item;
	Flag			first = TRUE;
	LinkPtr 		linker = NULL;

	if ( link == list->firstItem ) {
		ListInfo( list, stdout );	
		fprintf( stdout, "Used Dimensions are:\n\n");
	}

	fprintf( stdout, "The input value of '%s' has to be multiplied by %g", dimension->name, dimension->value );
	ListIterator( dimension->parameterValues, PrintParameterFactors, dimension->parameterValues );
	fprintf( stdout, ".\n\n" );
/*	ListIterator( dimension->parameterValues, PrintParameter, dimension->parameterValues );*/

	if ( dimension->orderOfReaction ) {
		fprintf( stdout, "The input value of '%s' has to be multiplied by %g^(order of reaction)\n", dimension->name, dimension->orderOfReacValue );
	}
	if ( dimension->tempExponent ) {
		fprintf( stdout, "The input value of '%s' has to be multiplied by %g^(exponent of temperature)\n", dimension->name, dimension->tempExpValue );
	}

	fprintf( stdout, "SI-units of '%s' are: [", dimension->name );

	if ( dimension->kg || dimension->kgParaCoeff ) {
		if ( !first ){
			fprintf( stdout, " *" );
		}
		fprintf( stdout, " kg" );
		if ( dimension->kg != 1 || dimension->kgParaCoeff != 0 ) {
			fprintf( stdout, "^(" );
			if ( dimension->kgParaCoeff != 0 ) {
				if ( dimension->kgParaCoeff != 1 ) {
					fprintf( stdout, "%d*%s", dimension->kgParaCoeff, dimension->kgParameter );
				}
				else {
					fprintf( stdout, "%s", dimension->kgParameter );
				}
				if ( dimension->kg > 0 ) {
					fprintf( stdout, "+%d)", dimension->kg );
				} else if ( dimension->kg < 0 ) {
					fprintf( stdout, "%d)", dimension->kg );
				} else {
					fprintf( stdout, ")" );
				}
			}
			else {
				fprintf( stdout, "%d)", dimension->kg );
			}
		}
		first = FALSE;
	}

	if ( dimension->m || dimension->mParaCoeff ) {
		if ( !first ){
			fprintf( stdout, " *" );
		}
		fprintf( stdout, " m" );
		if ( dimension->m != 1 || dimension->mParaCoeff != 0 ) {
			fprintf( stdout, "^(" );
			if ( dimension->mParaCoeff != 0 ) {
				if ( dimension->mParaCoeff != 1 ) {
					fprintf( stdout, "%d*%s", dimension->mParaCoeff, dimension->mParameter );
				}
				else {
					fprintf( stdout, "%s", dimension->mParameter );
				}
				if ( dimension->m > 0 ) {
					fprintf( stdout, "+%d)", dimension->m );
				} else if ( dimension->m < 0 ) {
					fprintf( stdout, "%d)", dimension->m );
				} else {
					fprintf( stdout, ")" );
				}
			}
			else {
				fprintf( stdout, "%d)", dimension->m );
			}
		}
		first = FALSE;
	}

	if ( dimension->s || dimension->sParaCoeff ) {
		if ( !first ){
			fprintf( stdout, " *" );
		}
		fprintf( stdout, " s" );
		if ( dimension->s != 1 || dimension->sParaCoeff != 0 ) {
			fprintf( stdout, "^(" );
			if ( dimension->sParaCoeff != 0 ) {
				if ( dimension->sParaCoeff != 1 ) {
					fprintf( stdout, "%d*%s", dimension->sParaCoeff, dimension->sParameter );
				}
				else {
					fprintf( stdout, "%s", dimension->sParameter );
				}
				if ( dimension->s > 0 ) {
					fprintf( stdout, "+%d)", dimension->s );
				} else if ( dimension->s < 0 ) {
					fprintf( stdout, "%d)", dimension->s );
				} else {
					fprintf( stdout, ")" );
				}
			}
			else {
				fprintf( stdout, "%d)", dimension->s );
			}
		}
		first = FALSE;
	}

	if ( dimension->K || dimension->KParaCoeff ) {
		if ( !first ){
			fprintf( stdout, " *" );
		}
		fprintf( stdout, " K" );
		if ( dimension->K != 1 || dimension->KParaCoeff != 0 ) {
			fprintf( stdout, "^(" );
			if ( dimension->KParaCoeff != 0 ) {
				if ( dimension->KParaCoeff != 1 ) {
					fprintf( stdout, "%d*%s", dimension->KParaCoeff, dimension->KParameter );
				}
				else {
					fprintf( stdout, "%s", dimension->KParameter );
				}
				if ( dimension->K > 0 ) {
					fprintf( stdout, "+%d)", dimension->K );
				} else if ( dimension->K < 0 ) {
					fprintf( stdout, "%d)", dimension->K );
				} else {
					fprintf( stdout, ")" );
				}
			}
			else {
				fprintf( stdout, "%d)", dimension->K );
			}
		}
		first = FALSE;
	}

	if ( dimension->mole || dimension->moleParaCoeff ) {
		if ( !first ){
			fprintf( stdout, " *" );
		}
		fprintf( stdout, " mole" );
		if ( dimension->mole != 1 || dimension->moleParaCoeff != 0 ) {
			fprintf( stdout, "^(" );
			if ( dimension->moleParaCoeff != 0 ) {
				if ( dimension->moleParaCoeff != 1 ) {
					fprintf( stdout, "%d*%s", dimension->moleParaCoeff, dimension->moleParameter );
				}
				else {
					fprintf( stdout, "%s", dimension->moleParameter );
				}
				if ( dimension->mole > 0 ) {
					fprintf( stdout, "+%d)", dimension->mole );
				} else if ( dimension->mole < 0 ) {
					fprintf( stdout, "%d)", dimension->mole );
				} else {
					fprintf( stdout, ")" );
				}
			}
			else {
				fprintf( stdout, "%d)", dimension->mole );
			}
		}
		first = FALSE;
	}

	fprintf( stdout, "]\n\n" );

	return 0;
}

int PrintParameter( LinkPtr link, void *var )
{
	ListPtr 		list = var;
	ParameterPtr 	parameter = link->item;

	if ( link == list->firstItem ) {
		ListInfo( list, stdout );	
		fprintf( stdout, "Used parameters are:\n");
	}

	if ( link->next ) {
		fprintf( stdout, "\tthe value of '%s' is %g\n", parameter->name, parameter->value );
	}
	else {
		fprintf( stdout, "\tthe value of '%s' is %g\n\n", parameter->name, parameter->value );
	}
	return 0;
}

int PrintParameterFactors( LinkPtr link, void *var )
{
	ListPtr 		list = var;
	ParameterPtr 	parameter = link->item;

	if ( parameter->value != 1.0 ) {
		fprintf( stdout, ", %g^%s", parameter->value, parameter->name );
	}

	return 0;
}

int PrintBroadening( LinkPtr link, void *var )
{
	ListPtr list = var;
	LinkPtr	contentsLink = NULL;
	LinkPtr	numberLink = NULL;
	String 	name = link->item;
	String	contents = NULL;
	int 	number = 0;
	int 	lineNumber;
	
	number = NumberOfListItem( gFcNameList, FindString, name );
	if ( number == 1 ) {
		ListInfo( list, stdout );	
		fprintf( stdout, "Specified Broadening functions are:\n");
	}
	contentsLink = Find_n_th_Item( number, list );
	contents = ( String )contentsLink->item;
	numberLink = Find_n_th_Item( number, gFcLineNumberList );
	lineNumber = *( int * )numberLink->item;
	fprintf( stdout, "File \"%s\"; Line %4d #    Fc%s = %s\n",
		gOptions->inName, lineNumber, name, contents );

	return 0;
}

int DumpAtoms( LinkPtr link, void *var )
{
	AtomsPtr 	atoms = ( AtomsPtr )link->item;
	FILE		*fp = gOptions->outFile;
	ListPtr		list = ( ListPtr )var;
	
	if ( !fwrite( atoms, sizeof( Atoms ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping struct of atoms\n");
		return 1;
	}
	strcpy( gbuffer, atoms->name );
	if ( !fwrite( gbuffer, gHeader->maxLenOfString, 1, fp ) ) {
		fprintf( stderr, "# error while dumping name of atom\n");
		return 1;
	}

	return 0;
}

int DumpSpecies( LinkPtr link, void *var )
{
	SpeciesPtr 	species = ( SpeciesPtr )link->item;
	FILE		*fp = gOptions->outFile;
	ListPtr		list = ( ListPtr )var;
		
	if ( !fwrite( species, sizeof( Species ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping struct of species\n");
		return 1;
	}
	strcpy( gbuffer, species->name );
	if ( !fwrite( gbuffer, gHeader->maxLenOfString, 1, fp ) ) {
		fprintf( stderr, "# error while dumping name of species\n");
		return 1;
	}
	if ( !fwrite( &species->composition->len, sizeof( int ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping length of compositionarray\n");
		return 1;
	}
	if ( fwrite( species->composition->vec, sizeof( int ), species->composition->len, fp ) != species->composition->len ) {
		fprintf( stderr, "# error while dumping compositionarray\n");
		return 1;
	}

	if ( !fwrite( &species->dCoeff->len, sizeof( int ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping length of dCoeffArray\n");
		return 1;
	}
	if ( fwrite( species->dCoeff->vec, sizeof( Double ), species->dCoeff->len, fp ) != species->dCoeff->len ) {
		fprintf( stderr, "# error while dumping dCoeffArray\n");
		return 1;
	}
	
	if ( !fwrite( &species->omegaCoeff->len, sizeof( int ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping length of omegaCoeffArray\n");
		return 1;
	}
	if ( fwrite( species->omegaCoeff->vec, sizeof( Double ), species->omegaCoeff->len, fp ) != species->omegaCoeff->len ) {
		fprintf( stderr, "# error while dumping omegaCoeffArray\n");
		return 1;
	}
	
	return 0;
}

int DumpReaction( LinkPtr link, void *var )
{
	ReactionPtr reaction = ( ReactionPtr )link->item;
	FILE		*fp = gOptions->outFile;
	ListPtr		list = ( ListPtr )var;
	
	if ( !fwrite( reaction, sizeof( Reaction ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping struct of reactions\n");
		return 1;
	}
	strcpy( gbuffer, reaction->label );
	if ( !fwrite( gbuffer, gHeader->maxLenOfString, 1, fp ) ) {
		fprintf( stderr, "# error while dumping label of reaction\n");
		return 1;
	}
	
	return 0;
}

int DumpThirdBody( LinkPtr link, void *var )
{
	ThirdBodyPtr 	thirdBody = ( ThirdBodyPtr )link->item;
	FILE			*fp = gOptions->outFile;
	ListPtr			list = ( ListPtr )var;
	int                     i;

	if ( !fwrite( thirdBody, sizeof( ThirdBody ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping struct of thirdbody efficienies\n");
		return 1;
	}
	strcpy( gbuffer, thirdBody->name );
	if ( !fwrite( gbuffer, gHeader->maxLenOfString, 1, fp ) ) {
		fprintf( stderr, "# error while dumping name of thirdbody efficieny\n");
		return 1;
	}

	/* thirbody species number has been initialized with -1. For those that are still -1, set to actual number and set coeff to zero */
	if ( !fwrite( &thirdBody->speciesNumber->len, sizeof( int ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping length of efficieny array\n");
		return 1;
	}
	if ( fwrite( thirdBody->speciesNumber->vec, sizeof( int ), thirdBody->speciesNumber->len, fp ) != thirdBody->speciesNumber->len ) {
		fprintf( stderr, "# error while dumping efficieny array\n");
		return 1;
	}

	if ( !fwrite( &thirdBody->speciesCoeff->len, sizeof( int ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping length of coefficients of efficieny array\n");
		return 1;
	}
	if ( fwrite( thirdBody->speciesCoeff->vec, sizeof( Double ), thirdBody->speciesCoeff->len, fp ) != thirdBody->speciesCoeff->len ) {
		fprintf( stderr, "# error while dumping coefficients of efficieny array\n");
		return 1;
	}

	return 0;
}

int DumpDimension( LinkPtr link, void *var )
{
	DimensionPtr 	dimension = ( DimensionPtr )link->item;
	FILE			*fp = gOptions->outFile;
	ListPtr			list = ( ListPtr )var;
	
	if ( !fwrite( dimension, sizeof( Dimension ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping struct of used units\n");
		return 1;
	}
	strcpy( gbuffer, dimension->name );
	if ( !fwrite( gbuffer, gHeader->maxLenOfString, 1, fp ) ) {
		fprintf( stderr, "# error while dumping name of used units\n");
		return 1;
	}

	return 0;
}

/*int CheckReactionZeroEntry( LinkPtr link, void *var )
{
	ReactionPtr 	reaction = link->item;
	int				i;

	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		if ( reaction->speciesCoeff[i] == 0.0 ) {
			return 1;
		}
	}

	return 0;
}
*/

int CheckReactionZeroEntry( LinkPtr link, void *var )
{
	ReactionPtr 	reaction = link->item;
	SpeciesPtr		species = NULL;
	ThirdBodyPtr	thirdBody = NULL;
	char			tbName[127];
	LinkPtr 		wLink;
	int				i;

	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		if ( reaction->speciesCoeff[i] == 0.0 ) {
			species = ( SpeciesPtr ) ListIterator( gSpeciesList, FindSpeciesNumber
								, &reaction->speciesNumber[i] )->item;
			sprintf( tbName, "M_%s", species->name );
			
			if ( !( wLink = ListIterator( gUsedThirdBodyList, FindString, tbName ) ) ) {
				AddString( gUsedThirdBodyList, tbName );
			}
	
			
			if ( !( ListIterator( gThirdBodyList, FindThirdBody, tbName ) ) ) {
				thirdBody = AddThirdBody( gThirdBodyList, tbName );
				thirdBody->speciesCoeff->vec[reaction->speciesNumber[i]] = 1.0;
			}

			reaction->withThirdBody = TRUE;
			reaction->thirdBodyNumber = NumberOfListItem( gUsedThirdBodyList, FindString, tbName )-1;
		}
	}

	return 0;
}

int CompareDoubleReaction( LinkPtr link, void *var )
{
	ReactionPtr		reaction1 = var;
	ReactionPtr 	reaction2 = link->item;
	static ListPtr	errorList = NULL;
	
	if ( !errorList ) {
		errorList = NewList( "errorList" );
	}

	if ( ListIterator( errorList, FindString, reaction2->label ) ) {
		return 0;
	}
	if ( reaction1 == reaction2 ) {
		return 0;
	}

	if ( CompareReaction( reaction1, reaction2 ) ) {
		AddString( errorList, reaction1->label );
		return 1;
	}

	return 0;
}

int CompareReaction( ReactionPtr reaction1, ReactionPtr reaction2 )
{
	int				equals = 0;
	int 			numberOfSpecies = reaction1->numberOfSpecies;
	int				i, j;

	if ( numberOfSpecies != reaction2->numberOfSpecies ) {
		return 0;
	}
	if ( reaction1->withThirdBody != reaction2->withThirdBody ) {
		return 0;
	}

	for ( i = 0; i < numberOfSpecies; ++i ) {
		for ( j = 0; j < numberOfSpecies; ++j ) {
			if ( reaction1->speciesNumber[i] == reaction2->speciesNumber[j] ) {
				if ( reaction1->speciesCoeff[i] == reaction2->speciesCoeff[j] ) {
					++equals;
					break;
				}
			}
		}
	}
	if ( equals == numberOfSpecies 
			&& reaction1->thirdBodyNumber == reaction2->thirdBodyNumber ) {
		return 1;
	}
	
	return 0;
}

Flag CompareThirdBody( LinkPtr link1, LinkPtr link2 )
{
	ThirdBodyPtr	tb1 = ( ThirdBodyPtr )link1->item;
	ThirdBodyPtr	tb2 = ( ThirdBodyPtr )link2->item;

	return ( tb1->id > tb2->id ) ? TRUE : FALSE;
}

void InitDimension( ListPtr list )
{
/*  the definition of 'value' is: 1 * [name] = value * [SI-units]   		*/

/*                      name		value			kg		m		s		K		mole  */	
	AddDimension( list, "kg", 		1.0, 			1, 		0, 		0, 		0, 		0 );
	AddDimension( list, "m", 		1.0, 			0, 		1, 		0, 		0, 		0 );
	AddDimension( list, "s", 		1.0, 			0, 		0, 		1, 		0, 		0 );
	AddDimension( list, "K", 		1.0, 			0, 		0, 		0, 		1, 		0 );
	AddDimension( list, "mole", 	1.0, 			0, 		0, 		0, 		0, 		1 );
	AddDimension( list, "mol", 		1.0, 			0, 		0, 		0, 		0, 		1 );

	AddDimension( list, "k", 		1.0e3, 			0, 		0, 		0, 		0, 		0 );
	AddDimension( list, "M", 		1.0e6, 			0, 		0, 		0, 		0, 		0 );
	AddDimension( list, "G", 		1.0e9, 			0, 		0, 		0, 		0, 		0 );
	AddDimension( list, "c", 		1.0e-2, 		0, 		0, 		0, 		0, 		0 );
	AddDimension( list, "milli", 	1.0e-3, 		0, 		0, 		0, 		0, 		0 );

	AddDimension( list, "ft", 		0.30480, 		0, 		1, 		0, 		0, 		0 );
	AddDimension( list, "lb_m", 	0.45359, 		1, 		0, 		0, 		0, 		0 );
	AddDimension( list, "lb_f", 	4.4482, 		1, 		1, 		-2, 	0, 		0 );
	AddDimension( list, "lb", 		0.45359, 		1, 		0, 		0, 		0, 		0 );
	AddDimension( list, "Lb", 		4.4482, 		1, 		1, 		-2, 	0, 		0 );
/* the following unit is the thermochemical calorie (cal), which is different to the
 * international steam table calorie (I. T. cal): 1.000654[cal] = 1[I. T. cal]
 */
	AddDimension( list, "cal", 		4.1840, 		1, 		2, 		-2, 	0, 		0 );
	AddDimension( list, "atm", 		1.0133e5,		1, 		-1, 	-2, 	0, 		0 );
	AddDimension( list, "h", 		60.0, 			0, 		0, 		1, 		0, 		0 );
	AddDimension( list, "J", 		1.0, 			1, 		2, 		-2, 	0, 		0 );
	AddDimension( list, "N", 		1.0, 			1, 		1, 		-2, 	0, 		0 );
	AddDimension( list, "W", 		1.0, 			1, 		2, 		-3, 	0, 		0 );
	AddDimension( list, "bar", 		1.0e5, 			1, 		-1, 	-2, 	0, 		0 );
}

void InitAllowedExponents( void )
{
	/* allowed parameters for the order of reaction */
	AddString( gAllowReacExpList, "n" );

	/* allowed parameters for the temperature exponent */
	AddString( gAllowTempExpList, "n_k" );
	AddString( gAllowTempExpList, "nk" );
}

int PrintThermoFile( LinkPtr link, void *var )
{
	int			i;
	char		buffer[128];
	struct CHEMFILES {
		FILE	*fpTherm;
		FILE	*fpTrans;
	} *chFiles = (struct CHEMFILES *) var;
	SpeciesPtr 	species = link->item;
	LinkPtr			atomLink = NULL;
	AtomsPtr 		atom = NULL;

/* write thermo data */
/* first line */
	fprintf( chFiles->fpTherm, "%-18s", CSymbolLeadNum( ShortSymbol( species->name, buffer ) ) );
	fprintf( chFiles->fpTherm, "%-6s", "000000" );
	for ( i = 0, atomLink = gAtomsList->firstItem; 
			i < gAtomsList->items; ++i, atomLink = atomLink->next ) {
		atom = (AtomsPtr)atomLink->item;
		fprintf( chFiles->fpTherm, "%-2s%3d", atom->name, species->composition->vec[i] );
	}
	fprintf( chFiles->fpTherm, "%1s", "G" );
	fprintf( chFiles->fpTherm, "%10.0f%10.0f%8.0f      1\n", 300.0, 5000.0, 1000.0 );

/* start with data */
	for ( i = 0; i < 5; ++i ) {
		fprintf( chFiles->fpTherm, "%15.8E", species->coeffHot[i] );
	}
	
	fprintf( chFiles->fpTherm, "    2\n" );
	
	for ( i = 5; i < 7; ++i ) {
		fprintf( chFiles->fpTherm, "%15.8E", species->coeffHot[i] );
	}
	for ( i = 0; i < 3; ++i ) {
		fprintf( chFiles->fpTherm, "%15.8E", species->coeffCold[i] );
	}
	 
	fprintf( chFiles->fpTherm, "    3\n" );
	
	for ( i = 3; i < 7; ++i ) {
		fprintf( chFiles->fpTherm, "%15.8E", species->coeffCold[i] );
	}

	fprintf( chFiles->fpTherm, "                   4\n" );


/* write transport data */
	fprintf( chFiles->fpTrans, "%-8s           0%10.3f%10.3f%10.3f%10.3f%10.3f\n"
		, CSymbolLeadNum( ShortSymbol( species->name, buffer ) )
		, 1.0 / species->k_over_eps, species->sigma, 0.0, 0.0, 0.0 );
	return 0;
}

int PrintChemkinReactions( LinkPtr link, void *var )
{
	FILE			*fp = (FILE *) var;
	ListPtr 		list = var;
	ReactionPtr 	reaction = link->item;
	LinkPtr			speciesLink = NULL;
	LinkPtr			thirdBodyLink = NULL;
	SpeciesPtr 		species = NULL;
	int 			i, dumint, tblen;
	Flag			first = TRUE;
	ThirdBodyPtr	thirdBody = NULL;
	char			*thirdBodyName;
	char			tmp[128];
	char			reac[128];
	char			shortsym[128];
	const Double	EKonvToCal = 1.0 / 4.184;

	if ( reaction->hasForward && reaction->forward != NULL ) {
/*	if ( reaction->hasForward ) {*/
		return 0;
	}

/*  print left side of equation  */
	reac[0] = '\0';
	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[i] );
		if ( !speciesLink ) FatalError( "something wrong in PrintReaction()" );
		species = ( SpeciesPtr ) speciesLink->item;

		if ( reaction->speciesCoeff[i] > 0.0 ) {
			if ( first ) {
				if ( fabs( reaction->speciesCoeff[i] ) != 1.0 ) {
					sprintf( tmp, "%g%s", fabs( reaction->speciesCoeff[i] )
							, CSymbolLeadNum( ShortSymbol( species->name, shortsym ) ) );
					strcat( reac, tmp );
					first = FALSE;
				}
				else {
					sprintf( tmp, "%s", CSymbolLeadNum( ShortSymbol( species->name, shortsym ) ) );
					strcat( reac, tmp );
					first = FALSE;
				}
			}
			else {
				if ( fabs( reaction->speciesCoeff[i] ) != 1.0 ) {
					sprintf( tmp, "+%g%s", fabs( reaction->speciesCoeff[i] ), CSymbolLeadNum( ShortSymbol( species->name, shortsym ) ) );
					strcat( reac, tmp );
				}
				else {
					sprintf( tmp, "+%s", CSymbolLeadNum( ShortSymbol( species->name, shortsym ) ) );
					strcat( reac, tmp );
				}
			}
		}
	}

	if ( reaction->withThirdBody && !reaction->withLindemann ) {
		sprintf( tmp, "+M" );
		strcat( reac, tmp );
	}
	if ( reaction->withLindemann ) {
		sprintf( tmp, "(+M)" );
		strcat( reac, tmp );
	}

/*  print assignment operator  */
	if ( ( reaction->hasBackward || reaction->hasForward ) ) {
		if ( reaction->forward == NULL && reaction->backward == NULL ) {
			sprintf( tmp, "=>" );
			strcat( reac, tmp );
		}
		else {
			sprintf( tmp, "<=>" );
			strcat( reac, tmp );
		}
	}
	else {
		sprintf( tmp, "=>" );
		strcat( reac, tmp );
	}

/*  print right side of equation  */
	first = TRUE;
	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[i] );
		if ( !speciesLink ) FatalError( "something wrong in PrintReaction()" );
		species = ( SpeciesPtr ) speciesLink->item;

		if ( reaction->speciesCoeff[i] < 0.0 ) {
			if ( first ) {
				if ( fabs( reaction->speciesCoeff[i] ) != 1.0 ) {
					sprintf( tmp, "%g%s", fabs( reaction->speciesCoeff[i] )
							, CSymbolLeadNum( ShortSymbol( species->name, shortsym ) )  );
					strcat( reac, tmp );
					first = FALSE;
				}
				else {
					sprintf( tmp, "%s", CSymbolLeadNum( ShortSymbol( species->name, shortsym ) ) );
					strcat( reac, tmp );
					first = FALSE;
				}
			}
			else {
				if ( fabs( reaction->speciesCoeff[i] ) != 1.0 ) {
					sprintf( tmp, "+%g%s", fabs( reaction->speciesCoeff[i] ), CSymbolLeadNum( ShortSymbol( species->name, shortsym ) ) );
					strcat( reac, tmp );
				}
				else {
					sprintf( tmp, "+%s", CSymbolLeadNum( ShortSymbol( species->name, shortsym ) ) );
					strcat( reac, tmp );
				}
			}
		}
	} 

	if ( reaction->withThirdBody && !reaction->withLindemann ) {
		sprintf( tmp, "+M" );
		strcat( reac, tmp );
	}
	if ( reaction->withLindemann ) {
		sprintf( tmp, "(+M)" );
		strcat( reac, tmp );
	}

	sprintf( tmp, " " );
	strcat( reac, tmp );
	dumint = 41 - strlen( reac );
	for ( i = 0; i < dumint; ++i ) {
		strcat( reac, tmp );
	}

	if ( reaction->withLindemann ) {
		sprintf( tmp, "%8.3e%9.3f%10.2f", reaction->aInf/* * GetAConverter( reaction, TRUE )*/
				, reaction->nInf
				, reaction ->eInf * GetEConverter( reaction, TRUE ) * EKonvToCal );
	}
	else {
		sprintf( tmp, "%8.3e%9.3f%10.2f", reaction->a/* * GetAConverter( reaction, FALSE )*/
				, reaction->n
				, reaction->e * GetEConverter( reaction, FALSE ) * EKonvToCal );
	}

	strcat( reac, tmp );

	fprintf( fp, "%s\n", reac );

	if ( reaction->withThirdBody ) {
		if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
			thirdBodyName = ( char * ) thirdBodyLink->item;
			thirdBodyLink = ListIterator( gThirdBodyList, FindThirdBody, thirdBodyName );
			if ( !thirdBodyLink ) FatalError( "something wrong in PrintChemkinReactions()" );
			thirdBody = ( ThirdBodyPtr ) thirdBodyLink->item;
		}
		else {
			fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
			exit( 2 );
		}
		tblen=0;
		for ( i = 0; i < thirdBody->speciesNumber->len; ++i ) {
			if ( thirdBody->speciesCoeff->vec[i] != 1.0 ) {
			  tblen += 7 + strlen(CSymbolLeadNum( ShortSymbol( species->name, shortsym ) ) );
			  if ( tblen > 80 ) {
				fprintf( fp, "\n" );
				tblen = 7 + strlen(CSymbolLeadNum( ShortSymbol( species->name, shortsym ) ) );
			  }
			  speciesLink = Find_n_th_Item( thirdBody->speciesNumber->vec[i]+1, gSpeciesList );
			  species = ( SpeciesPtr ) speciesLink->item;
			  fprintf( fp, "%s/%4.2f/ ", CSymbolLeadNum( ShortSymbol( species->name, shortsym ) )
					  , thirdBody->speciesCoeff->vec[i] );
			}
		}
		fprintf( fp, "\n" );
	}

	if ( reaction->withLindemann ) {
		fprintf( fp, "     LOW  /  " );
		fprintf( fp, "%8.3e%9.3f%10.2f /\n", reaction->a /** GetAConverter( reaction, FALSE )*/
				, reaction->n
				, reaction->e * GetEConverter( reaction, FALSE ) * EKonvToCal );
		if ( fabs(reaction->fcb - 1.0 + reaction->fca) < 1e-5 && fabs( reaction->fcc - 1.0 ) < 1e-5 ) { 
		  /* exact TROE form*/
			fprintf( fp, "     TROE/  " );
			fprintf( fp, "%6.4g %8.2f %9.2f %9.2f /\n", reaction->fcb, reaction->fcTa
						, reaction->fcTb, reaction->fcTc );
		}
		else if ( fabs(reaction->fcb - 1.0 + reaction->fca) < 1e-5 && fabs(reaction->fcc) < 1.0e-5 ) { 
		  /* exact TROE form without last term*/
			fprintf( fp, "!     FCCHECK/  " );
			fprintf( fp, "%6.4f %9.2f %6.4f %9.2f %6.4f %9.2f\n", reaction->fca, reaction->fcTa
					, reaction->fcb, reaction->fcTb, reaction->fcc, reaction->fcTc );
			fprintf( fp, "     TROE/  " );
			fprintf( fp, "%6.4g %8.2f %9.2f        /\n", reaction->fcb, reaction->fcTa
						, reaction->fcTb );
		}
		else if ( fabs(reaction->fca) < 1e-5 && fabs(reaction->fcb) < 1e-5 
				  && fabs(reaction->fcTc) < 1e-5) { 
			fprintf( fp, "!     FCCHECK/  " );
			fprintf( fp, "%6.4f %9.2f %6.4f %9.2f %6.4f %9.2f\n", reaction->fca, reaction->fcTa
					, reaction->fcb, reaction->fcTb, reaction->fcc, reaction->fcTc );
			fprintf( fp, "     TROE/  " );
			fprintf( fp, "%6.4g %8.2f %9.2f %9.2f /\n", reaction->fcc, 1.0, 1.0e7, 1.0e7 );
		}
		else if ( fabs(reaction->fca) < 1e-5 && fabs(reaction->fcb) < 1e-5) { 
			fprintf( fp, "!     FCCHECK/  " );
			fprintf( fp, "%6.4f %9.2f %6.4f %9.2f %6.4f %9.2f\n", reaction->fca, reaction->fcTa
					, reaction->fcb, reaction->fcTb, reaction->fcc, reaction->fcTc );
			fprintf( fp, "     TROE/  " );
			fprintf( fp, "%6.4g %8.2f %9.2f", 0.0, 1.0, 1.0 );
			if ( fabs(reaction->fcc) < 1e-5 ) {
			  fprintf( fp, " /\n" );
			}  
			else {
			  fprintf( fp, " %9.2g /\n", reaction->fcTc );
			}
		}
		else if ( fabs(reaction->fcb) < 1e-5 ) { 
			fprintf( fp, "!     FCCHECK/  " );
			fprintf( fp, "%6.4f %9.2f %6.4f %9.2f %6.4f %9.2f\n", reaction->fca, reaction->fcTa
					, reaction->fcb, reaction->fcTb, reaction->fcc, reaction->fcTc );
			fprintf( fp, "     TROE/  " );
			fprintf( fp, "%6.4g %8.2f %9.2f", reaction->fca, 1.0, reaction->fcTa );
			if ( fabs(reaction->fcc) < 1e-5 ) {
			  fprintf( fp, " /\n" );
			}  
			else {
			  fprintf( fp, " %9.2g /\n", reaction->fcTc );
			}
		}
		else {
			fprintf( fp, "!     FCCHECK/  " );
			fprintf( fp, "%6.4f %9.2f %6.4f %9.2f %6.4f %9.2f\n", reaction->fca, reaction->fcTa
					, reaction->fcb, reaction->fcTb, reaction->fcc, reaction->fcTc );
		}
	}
		
	if ( reaction->hasBackward && reaction->backward != NULL 
		&& reaction->backward->rateFromEquilibrium == FALSE ) {
		fprintf( fp, "     REV / " );
		fprintf( fp, "%8.3e%9.3f%10.2f /\n", reaction->backward->a
				, reaction->backward->n
				, reaction->backward->e * GetEConverter( reaction->backward, FALSE ) * EKonvToCal );
	}

	if ( reaction->isDouble == TRUE ) {
		fprintf( fp, "     DUPLICATE\n" );
	}
	
	return 0;
}

int PrintChemkinSpecies( LinkPtr link, void *var )
{
	FILE			*fp = (FILE *) var;
	static int		specinrow = 0;
	SpeciesPtr 		species = link->item;
	char			shortsym[128];

	if ( specinrow == 8 ) {
		fprintf( fp, "\n" );
		specinrow = 0;
	}
	fprintf( fp, "%-18s", CSymbolLeadNum( ShortSymbol( species->name, shortsym ) ) );
	++specinrow;
	
	return 0;
}

int PrintChemkinAtoms( LinkPtr link, void *var )
{
	FILE			*fp = (FILE *) var;
	AtomsPtr 		atom = link->item;

	fprintf( fp, "%-3s", atom->name );

	return 0;
}

int PrintMagicProdRateVectorize( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	LinkPtr			speciesLink = NULL;
	LinkPtr			thirdBodyLink = NULL;
	LinkPtr			backwLink = NULL;
	ReactionPtr		backw = NULL;
	SpeciesPtr 		species = NULL;
	int 			i, j, lenN, len;
	char			*thirdBody = NULL;
	String			label = NULL, labelN = NULL;
	char 			Wlabel[256];
	static int		reactionCount	=	0;
	const int		nOfReactions = 1;
	
	Wlabel[0] = '\0';
/* added to include backw reactions*/	
	if ( !reaction->hasForward ) {	/* it means reaction is not a backward reaction */
		if ( reaction->hasBackward ) {	/* it means there is a corresponding backward reaction */
			backwLink = link->next;
			backw = backwLink->item;
			lenN = strlen( backw->label );
			labelN = ( String ) malloc( lenN+1 );
			strcpy( labelN, backw->label );
			if ( labelN[lenN-1] == 'B' ) {
				labelN[lenN-1] = '\0';
				len = strlen( reaction->label ); 
				label = ( String ) malloc( len+1 );
				strcpy( label, reaction->label );
				label[len-1] = '\0';
				if ( strcmp( label, labelN ) == 0 ){
					backw = backwLink->item;
				}
				else {
					backw = NULL;
				}
			}
		}
	}
	else {
		return 0;
	}
/**/
	strcpy( Wlabel, reaction->label );
	if ( Wlabel[strlen( reaction->label )-1] == 'F' || Wlabel[strlen( reaction->label )-1] == 'B' ) {
		Wlabel[strlen( reaction->label )-1] = '\0';
	}
	if ( reactionCount % nOfReactions == 0 ) {
		fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);
	}
	fprintf( fp, "      W_R%s=KR(K,%s%s)", Wlabel, gPrefixReactionsF77, CSymbol( reaction->label ) );

	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[i] );
		if ( !speciesLink ) FatalError( "something wrong in PrintMagicReacRate()" );
		species = ( SpeciesPtr ) speciesLink->item;

		if ( reaction->speciesCoeff[i] > 0.0 ) {
			for ( j = 0; j < reaction->speciesCoeff[i]; ++j ) {
				fprintf( fp, "*C(K,%s%s)"
						, gPrefixComputedF77, CSymbol( species->name ) );
			}
		}
	}

	if ( reaction->withThirdBody && !reaction->withLindemann ) {
		if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
			thirdBody = ( char * ) thirdBodyLink->item;
			fprintf( fp, " * M(K,%s%s)",gPrefixThirdBodyF77, CSymbol( thirdBody ) );
		}
		else {
			fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
			exit( 2 );
		}
	}
/* added to include backw reactions*/	
	if ( backw ) {
		fprintf( fp, "\n     &    - KR(K,%s%s)", gPrefixReactionsF77, CSymbol( backw->label ) );
		for ( i = 0; i < backw->numberOfSpecies; ++i ) {
			speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &backw->speciesNumber[i] );
			if ( !speciesLink ) FatalError( "something wrong in PrintMagicReacRate()" );
			species = ( SpeciesPtr ) speciesLink->item;
	
			if ( backw->speciesCoeff[i] > 0.0 ) {
				for ( j = 0; j < backw->speciesCoeff[i]; ++j ) {
					fprintf( fp, "*C(K,%s%s)"
							, gPrefixComputedF77, CSymbol( species->name ) );
				}
			}
		}
	
		if ( backw->withThirdBody && !reaction->withLindemann ) {
			if ( thirdBodyLink = Find_n_th_Item( backw->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
				thirdBody = ( char * ) thirdBodyLink->item;
				fprintf( fp, " * M(K,%s%s)",gPrefixThirdBodyF77, CSymbol( thirdBody ) );
			}
			else {
				fprintf( stderr, "#something's wrong with the third body no.%d\n", backw->thirdBodyNumber );
				exit( 2 );
			}
		}
	
	}
/**/	
	fprintf( fp, "\n" );
	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {

				speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[i] );
				if ( !speciesLink ) FatalError( "something wrong in PrintMagicReacRate()" );
				species = ( SpeciesPtr ) speciesLink->item;
			    if ( !species->isSteadyState ) {
					fprintf( fp, "      CDOT(K,%s%s)=CDOT(K,%s%s)", gPrefixComputedF77, CSymbol( species->name ) 
					, gPrefixComputedF77, CSymbol( species->name ) );
					/* print sign */
					if ( reaction->speciesCoeff[i] < 0.0 ) {
						fprintf( fp, "+" );
					}
					else if ( reaction->speciesCoeff[i] > 0.0 ) {
						fprintf( fp, "%s", "-" );
					}
					
					/* print coefficient */
					
					if ( fabs( reaction->speciesCoeff[i] ) != 1.0 ) {
						fprintf( fp, "%G*", fabs( reaction->speciesCoeff[i] ) );
					}
					
					/* print reaction rate */
					fprintf( fp, "W_R%s\n", Wlabel );
				}
	}

	if ( ( reactionCount + 1 ) % nOfReactions == 0 || link->next == NULL ) {
		fprintf( fp, "%5d CONTINUE\n", gLabel);
		gLabel+=gLabelIncrement;
	}
	++reactionCount;
		
	return 0;
}

int PrintMagicMuCoeffF77( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	SpeciesPtr 		species = link->item;

	fprintf( fp, "      MUCOEFF(%s%s) = %15.8e\n"
			, gPrefixComputedF77, CSymbol( species->name ),
species->muCoeff );

	return 0;
}

int PrintMagicKOverEpsF77( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	SpeciesPtr 		species = link->item;

	fprintf( fp, "      KOVEREPS(%s%s) = %15.8e\n"
			, gPrefixComputedF77, CSymbol( species->name ),
species->k_over_eps );

	return 0;
}

int PrintDeclarattionOfW( LinkPtr link, void *var )
{
	FILE 			*fp = ( FILE * )var;
	ReactionPtr 	reaction = link->item;
	int 			len;
	String			label = NULL;
	char 			Wlabel[256];

	Wlabel[0] = '\0';
/* added to include backw reactions*/	
	if ( !reaction->hasForward ) {	/* it means reaction is not a backward reaction */
		if ( reaction->hasBackward ) {	/* it means there is a corresponding backward reaction */
			strcpy( Wlabel, reaction->label );
			if ( Wlabel[strlen( reaction->label )-1] == 'F' || Wlabel[strlen( reaction->label )-1] == 'B' ) {
				Wlabel[strlen( reaction->label )-1] = '\0';
			}
		}
		else{
			strcpy( Wlabel, reaction->label );	
		}
	}
	else {
		return 0;
	}
/**/
/*	strcpy( Wlabel, label );*/
/*	Wlabel[strlen( reaction->label )-1] = '\0';*/

	fprintf( fp, "      DOUBLE PRECISION W_R%s\n"
			, Wlabel );

	return 0;
}
