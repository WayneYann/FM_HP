/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 */

#include "ReactionParser.h"
#include "ReactionScan.tab.h"

#define NEWTREAT

#ifdef NEWTREAT
void PutSpeciesToReaction( ReactionPtr reaction, int speciesNumber, Double speciesCoeff )
{
	int j, found = 0;
	int	numberOfSpecies = reaction->numberOfSpecies;
	
	for ( j = 0; j < numberOfSpecies; ++j ) {
		if ( reaction->speciesNumber[j] == speciesNumber && reaction->speciesCoeff[j] > 0 ) {
			reaction->speciesCoeff[j] += speciesCoeff;
			found = 1;
			continue;
		}
	}
	if ( reaction->numberOfSpecies == MaxSpeciesOfReaction ) {
		fprintf( stderr, "###error: Too many species per reaction. Max allowed is %d. Increase MaxSpeciesOfReaction in ScanMan.h\n", MaxSpeciesOfReaction );
	}
	if ( !found ) {
		reaction->speciesNumber[numberOfSpecies] = speciesNumber;
		reaction->speciesCoeff[reaction->numberOfSpecies++] = speciesCoeff;
	}
}
#else
void PutSpeciesToReaction( ReactionPtr reaction, int speciesNumber, Double speciesCoeff )
{
	int j, found = 0;
	int	numberOfSpecies = reaction->numberOfSpecies;
	
	for ( j = 0; j < numberOfSpecies; ++j ) {
		if ( reaction->speciesNumber[j] == speciesNumber ) {
			if ( reaction->speciesCoeff[j] < 0 ) {
				fprintf( stderr, "#warning: species no. '%d' appears on both sides of reaction %s\n"
						, speciesNumber, reaction->label );
			}
			reaction->speciesCoeff[j] += speciesCoeff;
			found = 1;
			continue;
		}
	}
	if ( !found ) {
		reaction->speciesNumber[numberOfSpecies] = speciesNumber;
		reaction->speciesCoeff[reaction->numberOfSpecies++] = speciesCoeff;
	}
}
#endif

void yyerror( char *s )
{
	if ( !gErrorOut ) {
		gErrorOut = TRUE;
		fprintf( stderr, "\n# parse error:\n#\n" );  
		if ( strcmp( s, "parse error" ) ) {
			fprintf( stderr, "# %s\n#\n", s );
		}
		fprintf( stderr, "# %s\n#\n", gLineContents );
		switch (gtypeOfYylval) {
			case PTRDOUBLE:
				fprintf( stderr, "# the last thing i saw was: '%g'\n\n", *yylval.ptrdouble );
				break;
			case TYPDOUBLE:
				fprintf( stderr, "# the last thing i saw was: '%g'\n\n", yylval.typdouble );
				break;
			case TYPINT:
				fprintf( stderr, "# the last thing i saw was: '%d'\n\n", yylval.typint );
				break;
			case TYPSTRING:
				fprintf( stderr, "# the last thing i saw was: '%s'\n\n", yylval.typstring );
				break;
			case TYPREACTION: 
				fprintf( stderr, "# the last thing i saw was the reactionlabel: '%s'\n\n", yylval.typreaction->label );
				break;
			case TYPDIMENSION:
				fprintf( stderr, "# something wrong with the type of yylval ( typdimension )\n\n" );
				break;
			case TYPDIMENSIONPTR:
				fprintf( stderr, "# the last thing i saw was the factor of a unit with the value: '%g'\n\n", yylval.typdimensionptr->value );
				break;
			case TYPEXPONENT:
				fprintf( stderr, "# the last thing i saw was an exponent:\n" );
				fprintf( stderr, "#                        value: '%d'\n", yylval.typexponent.number );
				fprintf( stderr, "# coefficient of the parameter: '%d'\n", yylval.typexponent.paraCoeff );
				fprintf( stderr, "#                    parameter: '%s'\n\n", yylval.typexponent.parameter );
				break;
			default:
				fprintf( stderr, "# i can't find a valid type to print the last value i've matched\n\n" );
		}
		fprintf( stderr, "#-----------------------------------------------------------------\n" );
#ifdef applec
		fprintf( stderr, " File \"%s\"; Line %d\n", gOptions->inName, gNumberOfLines );
#else
		fprintf( stderr, " vi +%d %s\n", gNumberOfLines, gOptions->inName );
#endif
		fprintf( stderr, "#-----------------------------------------------------------------\n#\n\n\n" );
	
		if ( yynerrs > 10 ) {
			fprintf( stderr, "# make errors fewer" );
			fprintf( stderr, "#" );
			PrintLists();
			CleanUp();
			exit( 2 );
		}
	}
}

void ChangeSignsOfReactionCoeffs( ReactionPtr reaction, Flag arrowToRight ) 
{
	static Flag	rightArrow = 1;
	int	i;
	int	numberOfSpecies = reaction->numberOfSpecies;
	
	if ( rightArrow ) {
		for ( i = 0; i < numberOfSpecies; ++i ) {
			reaction->speciesCoeff[i] *= -1.0;
		}
	}
	rightArrow = arrowToRight;
}

int SetOneOrderOfReaction( ReactionPtr reaction )
{
	int i;
	
	if ( reaction->orderOfReaction == 0 ) {
		reaction->orderOfReaction = ( reaction->withThirdBody ) ? 1: 0;
		for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
			if ( reaction->speciesCoeff[i] > 0.0 ) {
				reaction->orderOfReaction += ( int )reaction->speciesCoeff[i];
			}
		}
		if ( reaction->withLindemann && !reaction->withThirdBody ) {
			reaction->orderOfReaction += 1;
		}
	}
	
	return 0;
}

int SetOrderOfReaction( LinkPtr link, void *var )
{
	ReactionPtr	reaction = link->item;
	
	SetOneOrderOfReaction( reaction );
	
	return 0;
}

int ShiftLetters( DimensionPtr dimension, char *dimName )
{
	LinkPtr			link = NULL;
	DimensionPtr 	allowedDim = NULL;
	
	memset( dimension, 0, sizeof(Dimension) );
	dimension->value = 1.0;
	dimension->parameterValues = NewList( "Parameters of units" );
	
	if ( link = ListIterator( gDimensionList, FindDimension, dimName ) ) {
		allowedDim = ( DimensionPtr )link->item;
		dimension->value *= allowedDim->value;
		dimension->kg += allowedDim->kg;
		dimension->m += allowedDim->m;
		dimension->s += allowedDim->s;
		dimension->K += allowedDim->K;
		dimension->mole += allowedDim->mole;
	}
	else {
		ScanError( "%s is not an allowed unit", dimName, 0 );
		return 1;
	}

/*	fprintf(stderr, "value of %s is %g\n", dimName, dimension->value );*/
	return 0;
}

int DivideDimension( DimensionPtr dimension, Dimension numerator, Dimension denominator )
{
	ListPtr 		listFirst = numerator.parameterValues;
	ListPtr 		listSecond = denominator.parameterValues;
	LinkPtr 		linkFirst = NULL;
	LinkPtr			linkSecond = NULL;
	ParameterPtr	parameterFirst = NULL;
	ParameterPtr	parameterSecond = NULL;

	dimension->value = numerator.value / denominator.value;

	dimension->kg = numerator.kg - denominator.kg;
	dimension->m = numerator.m - denominator.m;
	dimension->s = numerator.s - denominator.s;
	dimension->K = numerator.K - denominator.K;
	dimension->mole = numerator.mole - denominator.mole;

	for ( linkFirst = listFirst->firstItem; linkFirst != 0; linkFirst = linkFirst->next ) {
		parameterFirst = ( ParameterPtr )linkFirst->item;
		if ( linkSecond = ListIterator( listSecond, FindParameter, parameterFirst->name ) ) {
			parameterSecond = ( ParameterPtr )linkSecond->item;
			parameterFirst->value = parameterFirst->value / parameterSecond->value;
		}
	} 

	for ( linkSecond = listSecond->firstItem; linkSecond != 0; linkSecond = linkSecond->next ) {
		parameterSecond = ( ParameterPtr )linkSecond->item;
		if ( !( linkFirst = ListIterator( listFirst, FindParameter, parameterSecond->name ) ) ) {
			parameterFirst = AddParameter( listFirst, parameterSecond->name );
			parameterFirst->value = 1.0 / parameterSecond->value;
		}
	} 

	dimension->parameterValues = listFirst;

	if ( strcmp( numerator.kgParameter, denominator.kgParameter ) == 0 ) {
		strcpy( dimension->kgParameter, numerator.kgParameter );
	}
	else if ( ChooseOneParameter( dimension->kgParameter, numerator.kgParameter, numerator.kgParaCoeff, 
										denominator.kgParameter, denominator.kgParaCoeff ) ) {
		return 1;
	}
	dimension->kgParaCoeff = numerator.kgParaCoeff - denominator.kgParaCoeff;

	if ( strcmp( numerator.mParameter, denominator.mParameter ) == 0 ) {
		strcpy( dimension->mParameter, numerator.mParameter );
	}
	else if ( ChooseOneParameter( dimension->mParameter, numerator.mParameter, numerator.mParaCoeff, 
										denominator.mParameter, denominator.mParaCoeff ) ) {
		return 1;
	}
	dimension->mParaCoeff = numerator.mParaCoeff - denominator.mParaCoeff;

	if ( strcmp( numerator.sParameter, denominator.sParameter ) == 0 ) {
		strcpy( dimension->sParameter, numerator.sParameter );
	}
	else if ( ChooseOneParameter( dimension->sParameter, numerator.sParameter, numerator.sParaCoeff, 
										denominator.sParameter, denominator.sParaCoeff ) ) {
		return 1;
	}
	dimension->sParaCoeff = numerator.sParaCoeff - denominator.sParaCoeff;

	if ( strcmp( numerator.KParameter, denominator.KParameter ) == 0 ) {
		strcpy( dimension->KParameter, numerator.KParameter );
	}
	else if ( ChooseOneParameter( dimension->KParameter, numerator.KParameter, numerator.KParaCoeff, 
										denominator.KParameter, denominator.KParaCoeff ) ) {
		return 1;
	}
	dimension->KParaCoeff = numerator.KParaCoeff - denominator.KParaCoeff;

	if ( strcmp( numerator.moleParameter, denominator.moleParameter ) == 0 ) {
		strcpy( dimension->moleParameter, numerator.moleParameter );
	}
	else if ( ChooseOneParameter( dimension->moleParameter, numerator.moleParameter, numerator.moleParaCoeff, 
										denominator.moleParameter, denominator.moleParaCoeff ) ) {
		return 1;
	}
	dimension->moleParaCoeff = numerator.moleParaCoeff - denominator.moleParaCoeff;
	
#ifdef DEBUGPARSEUNITS
	ListIterator( dimension->parameterValues, PrintParameter, dimension->parameterValues );
#endif

	/* clean unused list listSecond */
	FreeList( listSecond, FreeParameter );

	return 0;
}

int MultiplyDimension( DimensionPtr dimension, Dimension first, Dimension second )
{
	ListPtr 		listFirst = first.parameterValues;
	ListPtr 		listSecond = second.parameterValues;
	LinkPtr 		linkFirst = NULL;
	LinkPtr			linkSecond = NULL;
	ParameterPtr	parameterFirst = NULL;
	ParameterPtr	parameterSecond = NULL;
	
	dimension->value = first.value * second.value;

	dimension->kg = first.kg + second.kg;
	dimension->m = first.m + second.m;
	dimension->s = first.s + second.s;
	dimension->K = first.K + second.K;
	dimension->mole = first.mole + second.mole;

	for ( linkFirst = listFirst->firstItem; linkFirst != 0; linkFirst = linkFirst->next ) {
		parameterFirst = ( ParameterPtr )linkFirst->item;
		if ( linkSecond = ListIterator( listSecond, FindParameter, parameterFirst->name ) ) {
			parameterSecond = ( ParameterPtr )linkSecond->item;
			parameterFirst->value = parameterFirst->value * parameterSecond->value;
		}
	} 

	for ( linkSecond = listSecond->firstItem; linkSecond != 0; linkSecond = linkSecond->next ) {
		parameterSecond = ( ParameterPtr )linkSecond->item;
		if ( !( linkFirst = ListIterator( listFirst, FindParameter, parameterSecond->name ) ) ) {
			parameterFirst = AddParameter( listFirst, parameterSecond->name );
			parameterFirst->value = parameterSecond->value;
		}
	} 

	dimension->parameterValues = listFirst;

	if ( strcmp( first.kgParameter, second.kgParameter ) == 0 ) {
		strcpy( dimension->kgParameter, first.kgParameter );
	}
	else if ( ChooseOneParameter( dimension->kgParameter, first.kgParameter, first.kgParaCoeff, 
										second.kgParameter, second.kgParaCoeff ) ) {
			return 1;
	} 
	dimension->kgParaCoeff = first.kgParaCoeff + second.kgParaCoeff;

	if ( strcmp( first.mParameter, second.mParameter ) == 0 ) {
		strcpy( dimension->mParameter, first.mParameter );
	}
	else if ( ChooseOneParameter( dimension->mParameter, first.mParameter, first.mParaCoeff, 
										second.mParameter, second.mParaCoeff ) ) {
		return 1;
	}
	dimension->mParaCoeff = first.mParaCoeff + second.mParaCoeff;

	if ( strcmp( first.sParameter, second.sParameter ) == 0 ) {
		strcpy( dimension->sParameter, first.sParameter );
	}
	else if ( ChooseOneParameter( dimension->sParameter, first.sParameter, first.sParaCoeff, 
										second.sParameter, second.sParaCoeff ) ) {
		return 1;
	}
	dimension->sParaCoeff = first.sParaCoeff + second.sParaCoeff;

	if ( strcmp( first.KParameter, second.KParameter ) == 0 ) {
		strcpy( dimension->KParameter, first.KParameter );
	}
	else if ( ChooseOneParameter( dimension->KParameter, first.KParameter, first.KParaCoeff, 
										second.KParameter, second.KParaCoeff ) ) {
		return 1;
	}
	dimension->KParaCoeff = first.KParaCoeff + second.KParaCoeff;

	if ( strcmp( first.moleParameter, second.moleParameter ) == 0 ) {
		strcpy( dimension->moleParameter, first.moleParameter );
	}
	else if ( ChooseOneParameter( dimension->moleParameter, first.moleParameter, first.moleParaCoeff, 
										second.moleParameter, second.moleParaCoeff ) ) {
		return 1;
	}
	dimension->moleParaCoeff = first.moleParaCoeff + second.moleParaCoeff;

#ifdef DEBUGPARSEUNITS
	ListIterator( dimension->parameterValues, PrintParameter, dimension->parameterValues );
#endif

	/* clean unused list listSecond */
	FreeList( listSecond, FreeParameter );

	return 0;
}

void CopyDimension( DimensionPtr dimensionPtr, Dimension dimension )
{
	dimensionPtr->value = dimension.value;

	dimensionPtr->kg = dimension.kg;
	dimensionPtr->m = dimension.m;
	dimensionPtr->s = dimension.s;
	dimensionPtr->K = dimension.K;
	dimensionPtr->mole = dimension.mole;
	
	dimensionPtr->kgParaCoeff = dimension.kgParaCoeff;
	dimensionPtr->mParaCoeff = dimension.mParaCoeff;
	dimensionPtr->sParaCoeff = dimension.sParaCoeff;
	dimensionPtr->KParaCoeff = dimension.KParaCoeff;
	dimensionPtr->moleParaCoeff = dimension.moleParaCoeff;
	
	strcpy( dimensionPtr->kgParameter, dimension.kgParameter );
	strcpy( dimensionPtr->mParameter, dimension.mParameter );
	strcpy( dimensionPtr->sParameter, dimension.sParameter );
	strcpy( dimensionPtr->KParameter, dimension.KParameter );
	strcpy( dimensionPtr->moleParameter, dimension.moleParameter );

	dimensionPtr->parameterValues = dimension.parameterValues;
}

int PowerDimension( DimensionPtr dimension, Dimension base, Exponent exp )
{
	LinkPtr link = NULL;
	ParameterPtr parameter = NULL;

	dimension->value = ( Double ) pow( base.value, ( Double ) exp.number );
	dimension->kg = base.kg * exp.number;
	dimension->m = base.m * exp.number;
	dimension->s = base.s * exp.number;
	dimension->K = base.K * exp.number;
	dimension->mole = base.mole * exp.number;

	if ( exp.paraCoeff ) {
		if ( link = ListIterator( dimension->parameterValues, FindParameter, exp.parameter ) ) {
			ScanError( "there are only linear exponents in terms of the parameter '%s' allowed", exp.parameter, 0 );
			return 1;
		}
		else {
			parameter = AddParameter( dimension->parameterValues, exp.parameter );
			parameter->value = pow( base.value, ( Double ) exp.paraCoeff );
		}
	}

	if ( ChooseOneParameter( dimension->kgParameter, base.kgParameter, base.kgParaCoeff, 
										exp.parameter, exp.paraCoeff ) ) {
		return 1;									
	}
	if ( ChooseOneParameter( dimension->mParameter, base.mParameter, base.mParaCoeff, 
										exp.parameter, exp.paraCoeff ) ) {
		return 1;									
	}
	if ( ChooseOneParameter( dimension->sParameter, base.sParameter, base.sParaCoeff, 
										exp.parameter, exp.paraCoeff ) ) {
		return 1;									
	}
	if ( ChooseOneParameter( dimension->KParameter, base.KParameter, base.KParaCoeff, 
										exp.parameter, exp.paraCoeff ) ) {
		return 1;									
	}
	if ( ChooseOneParameter( dimension->moleParameter, base.moleParameter, base.moleParaCoeff, 
										exp.parameter, exp.paraCoeff ) ) {
		return 1;									
	}

	dimension->kgParaCoeff = base.kgParaCoeff * exp.number + base.kg * exp.paraCoeff;
	dimension->mParaCoeff = base.mParaCoeff * exp.number + base.m * exp.paraCoeff;
	dimension->sParaCoeff = base.sParaCoeff * exp.number + base.s * exp.paraCoeff;
	dimension->KParaCoeff = base.KParaCoeff * exp.number + base.K * exp.paraCoeff;
	dimension->moleParaCoeff = base.moleParaCoeff * exp.number + base.mole * exp.paraCoeff;

#ifdef DEBUGPARSEUNITS
	ListIterator( dimension->parameterValues, PrintParameter, dimension->parameterValues );
#endif

	return 0;
}

int MuliplyExponent( ExponentPtr exp, Exponent first, Exponent second )
{
	exp->number = first.number * second.number;

	if ( ChooseOneParameter( exp->parameter, first.parameter, first.paraCoeff, 
										second.parameter, second.paraCoeff ) ) {
		return 1;
	}

	exp->paraCoeff = first.number * second.paraCoeff + first.paraCoeff * second.number;	

	return 0;
}

int AddExponent( ExponentPtr exp, Exponent first, Exponent second )
{
	exp->number = first.number + second.number;

	if ( strcmp( first.parameter, second.parameter ) == 0 ) {
		strcpy( exp->parameter, first.parameter );
	}
	else if ( ChooseOneParameter( exp->parameter, first.parameter, first.paraCoeff, 
										second.parameter, second.paraCoeff ) ) {
		return 1;
	}
	
	exp->paraCoeff = first.paraCoeff + second.paraCoeff;	

	return 0;
}

int DivideExponent( ExponentPtr exp, Exponent first, Exponent second )
{
	if ( !second.parameter[0] || !second.paraCoeff ) {
		exp->number = first.number / second.number;
		strcpy( exp->parameter, first.parameter );
		exp->paraCoeff = first.paraCoeff + second.number;	
	}
	else {
		ScanError("it is not valid to use a parameter '%s' as a denominator of an exponent", second.parameter, 0 );
		return 1;
	}
	
	return 0;
}

int SubtractExponent( ExponentPtr exp, Exponent first, Exponent second )
{
	exp->number = first.number - second.number;

	if ( strcmp( first.parameter, second.parameter ) == 0 ) {
		strcpy( exp->parameter, first.parameter );
	} else if ( ChooseOneParameter( exp->parameter, first.parameter, first.paraCoeff, 
										second.parameter, second.paraCoeff ) ) {
		return 1;
	}
	
	exp->paraCoeff = first.paraCoeff - second.paraCoeff;	
	
	return 0;
}

int ChooseOneParameter( char *exp, char *firstPara, int firstCoeff, char *secondPara, int secondCoeff )
{
	if ( firstCoeff && secondCoeff ) {
		strcat( firstPara, " and " );
		strcat( firstPara, secondPara );
		ScanError("there should be almost one parameter for each SI-unit, i found '%s'", firstPara, 0 );
		return 1;
	} 
	else if ( firstCoeff ) {
		strcpy( exp, firstPara );
	}
	else if ( secondCoeff ) {
		strcpy( exp, secondPara );
	}
	else {
		exp[0] = '\0';
	}

	return 0;
}

void ShiftExponent( ExponentPtr exp, Exponent letter )
{
	exp->number = letter.number;
	strcpy( exp->parameter, letter.parameter );
	exp->paraCoeff = letter.paraCoeff;
}

void SetExponentFlags( DimensionPtr dimension )
{
	int 			i;
	ListPtr 		list = dimension->parameterValues;
	LinkPtr 		link = NULL;
	ParameterPtr 	para = NULL;
	
	for ( i = 0; i < list->items; ++i ) {
		link = Find_n_th_Item( i+1, list );
		para = ( ParameterPtr )link->item;
		if ( ListIterator( gAllowTempExpList, FindString, para->name ) ) {
			dimension->tempExponent = TRUE;
			dimension->tempExpValue = para->value;
		}
		if ( ListIterator( gAllowReacExpList, FindString, para->name ) ) {
			dimension->orderOfReaction = TRUE;
			dimension->orderOfReacValue = para->value;
		}
	}

}

void WriteCompoOfReaction( ReactionPtr reaction )
{
	int			i;
	SpeciesPtr	species;

	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		species = ( SpeciesPtr ) ListIterator( gSpeciesList, FindSpeciesNumber
								, &reaction->speciesNumber[i] )->item;
		fprintf( stderr, "Species %s has the coefficient %g:\n\t"
				, species->name, reaction->speciesCoeff[i] );
		PrintSpeciesComposition( species, stderr );
	}
}
