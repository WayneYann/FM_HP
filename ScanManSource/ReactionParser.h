/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 */

#ifndef __ScanMan__
#include "ScanMan.h"
#endif

#undef DEBUGPARSESPECIES
#undef DEBUGPARSEUNITS

typedef struct Exponent {
	char	parameter[256];
	int		number;
	int		paraCoeff;
} Exponent, *ExponentPtr;

enum yylvalTypes { PTRDOUBLE, TYPDOUBLE, TYPINT, TYPSTRING, TYPREACTION, TYPDIMENSION,
					TYPDIMENSIONPTR, TYPEXPONENT }; 

void PutSpeciesToReaction( ReactionPtr reaction, int speciesNumber, Double speciesCoeff );

void ChangeSignsOfReactionCoeffs( ReactionPtr reaction, Flag arrowToRight );

void yyerror( char * );

int ShiftLetters( DimensionPtr dimension, char *dimName );

int DivideDimension( DimensionPtr dimension, Dimension numerator, Dimension denominator );

int MultiplyDimension( DimensionPtr dimension, Dimension first, Dimension second );

void CopyDimension( DimensionPtr dimensionPtr, Dimension dimension );

int PowerDimension( DimensionPtr dimension, Dimension base, Exponent exp );

int MuliplyExponent( ExponentPtr exp, Exponent first, Exponent second );

int AddExponent( ExponentPtr exp, Exponent first, Exponent second );

int DivideExponent( ExponentPtr exp, Exponent first, Exponent second );

int SubtractExponent( ExponentPtr exp, Exponent first, Exponent second );

int ChooseOneParameter( char *exp, char *firstPara, int firstCoeff, char *secondPara, int secondCoeff );

void ShiftExponent( ExponentPtr exp, Exponent letter );

void SetExponentFlags( DimensionPtr dimension );

int yylex( void );

void WriteCompoOfReaction( ReactionPtr reaction );

#define RIGHT 1
#define LEFT 0
