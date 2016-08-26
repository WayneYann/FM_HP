/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 */
 
#ifndef __ScanMan__
#define __ScanMan__

#define VERSION "3.3.9"

#include "ListMan.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef __ArrayManager__
#include "ArrayManager.h"
#endif
#ifdef applec
#include <CursorCtl.h>
#endif


#define MAX( a, b )     (( (a) > (b) ) ? (a) : (b))
#define MIN( a, b )     (( (a) < (b) ) ? (a) : (b))

/* [J / kmole K] */
#ifndef RGAS
#define RGAS 8314.34
#endif

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif
#define LINELENGTH 512
#define MaxSpeciesOfReaction 100

typedef int ReactionRateType;

typedef char *String;

typedef struct Options {
	char	*base;			/* base name for file names that are derived from	*/
	char	*inName;		/*	the input file name.							*/
	char	*outName;		/* binary output file								*/
	char	*thermoName;
	char	*steadyStates;	/* file containing steady state assumptions			*/
	FILE	*reactionInFile;
	FILE	*thermoInFile;
	FILE	*outFile;		/* binary output file								*/
	Flag	printAtoms;
	Flag	printSpecies;
	Flag	printReactions;
	Flag	printThird;
	Flag	printUnits;
	Flag	debugScanner;	/* debug parser flag in yydebug  */
	Flag	progress;
	Flag	useForwardBackward;
	Flag	latexTable;
	Flag	mathematica;
	Flag	outputWriter;
	Flag	outputWriterF77;
	Flag	outputWriterF90;
	Flag	vectorizeF77;
	Flag	printBroadening;
	Flag	noBinaryOutput;	/* if TRUE, no binary output file is generated (jg)	*/
	Flag	CCode;			/* if TRUE, generate C code for k and w (jg)		*/
	int		Z;				/* flag for exploring things (RESERVED FOR JG!)		*/
} Options, *OptionsPtr;

typedef struct Header {
	char	version[32];
	int		maxLenOfString;
	Flag	globalMechanism;
} Header, *HeaderPtr;

typedef struct Counter {
	int	atoms;
	int	species;
	int	pahSpecies;
	int	sootSpecies;
	int	steadyStates;
	int	reactions;
	int	pahReactions;
	int	sootReactions;
	int	thirdBodies;
	int dimensions;
} Counter, *CounterPtr;


extern ListPtr				gAtomsList;
extern ListPtr				gSpeciesList;
extern ListPtr				gReactionList;
extern ListPtr				gPAHReactionList;
extern ListPtr				gSootReactionList;
extern ListPtr				gAllowAtomList;
extern ListPtr				gUsedThirdBodyList;
extern ListPtr				gThirdBodyList;
extern ListPtr				gDimensionList;
extern ListPtr				gAllowTempExpList;
extern ListPtr				gAllowReacExpList;
extern ListPtr				gUsedDimensionList;
extern ListPtr				gCommentList;
extern ListPtr				gLabelFollowingCommentList;
extern ListPtr				gSteadyStateList;
extern IntVectorPtr			gComposition;
extern int					gNumberOfLines;
extern char					gLineContents[];
extern int					gtypeOfYylval;
extern int 					gErrorOut;
extern OptionsPtr			gOptions;
extern HeaderPtr			gHeader;
extern CounterPtr			gCounter;
extern char					*gbuffer;
extern ListFuncPtr			gSpeciesCompareFunc;

extern ListPtr				gFcNameList;
extern ListPtr				gFcContentsList;
extern ListPtr				gFcLineNumberList;

extern const char			*gPrefixSteadyState;
extern const char			*gPrefixComputed;
extern const char			*gPrefixReactions;
extern const char			*gPrefixThirdBody;
extern const char			*gPrefixSteadyStateF77;
extern const char			*gPrefixComputedF77;
extern const char			*gPrefixReactionsF77;
extern const char			*gPrefixThirdBodyF77;

extern int 					yyleng;
extern char			 		*yytext;

#ifdef CRAY
extern long					yynerrs;
extern long					yydebug;
#else
extern int					yynerrs;
extern int					yydebug;
#endif

typedef struct Atoms {
	char			*name;
	int				number;
} Atoms, *AtomsPtr;

typedef struct Species {
	char			*name;
	IntVectorPtr	composition;
	int				number;
	int				isSteadyState;	/* signals whether this is a steady state species		*/
	int				aux;			/* auxiliary field needed for reduced mechanisms		*/
	Double			molarMass;
	Double			coeffHot[7];
	Double			coeffCold[7];
	Double			k_over_eps;		/* temperature factor for evaluating sigma (= k/eps);
									   input is eps/k  (computed)							*/
	Double			sigma;
	Double			muCoeff;		/* constant factor for evaluation of mu in kg/(m*s)
									   (= 2.6693e-6 * sqrt(M)/sigma)  (computed)			*/
	VectorPtr		dCoeff;			/* constant factor for evaluation of D in m^2/s 
									   = 0.0188292*sqrt(1/Mi+1/Mj)/(.5*(sigma_i+sigma_j))^2
									   (computed)											*/
	VectorPtr		omegaCoeff;		/* constant factor for evaluation of omega_D 
									   (= sqrt(k_over_eps_i * k_over_eps_j))  (computed)	*/
} Species, *SpeciesPtr;

typedef struct Reaction {
	char		*label;
	int			id;					/* (jg) for internal purposes							*/
	int			numberOfSpecies;	/*  without thirdBody  									*/
	int			speciesNumber[MaxSpeciesOfReaction];
	Double		speciesCoeff[MaxSpeciesOfReaction]; /* speciesCoeff[products] < 0 */
	Double		a;
	Double		n;
	Double		e;
	Flag		withLindemann;
	Double		aInf;
	Double		nInf;
	Double		eInf;
	int			lindemannNumber;	/* the lowest value for lindemannNumber is 0			*/
	Double		fca;
	Double		fcb;
	Double		fcc;
	Double		fcTa;
	Double		fcTb;
	Double		fcTc;
	Flag		withThirdBody;
	int			thirdBodyNumber;	/* the lowest value for thirdBodyNumber is 0			*/
	Flag		partialEquilibrium;
	int			orderOfReaction;	
	Flag		hasForward;
	Flag		hasBackward;
	struct Reaction	*forward;		/*	pointer to forward reaction							*/
	struct Reaction	*backward;		/*	pointer to backward reaction, or NULL				*/
	Flag          isDouble;
	Flag          rateFromEquilibrium;
} Reaction, *ReactionPtr;

typedef struct ThirdBody {
	char			*name;
	int				id;
	IntVectorPtr	speciesNumber;
	VectorPtr		speciesCoeff;
        IntVectorPtr   set;  /* 0 means not set, 1 means set */
} ThirdBody, *ThirdBodyPtr;

typedef struct Dimension {
	char			*name;
	Double			value;
	int				kg;
	int				kgParaCoeff;
	char			kgParameter[10];
	int				m;
	int				mParaCoeff;
	char			mParameter[10];
	int				s;
	int				sParaCoeff;
	char			sParameter[10];
	int				K;
	int				KParaCoeff;
	char			KParameter[10];
	int				mole;
	int				moleParaCoeff;
	char			moleParameter[10];
	Flag			orderOfReaction;
	Double			orderOfReacValue;
	Flag			tempExponent;
	Double			tempExpValue;
	ListPtr			parameterValues;
} Dimension, *DimensionPtr;

typedef struct Parameter {
	char			*name;
	Double			value;
} Parameter, *ParameterPtr;

AtomsPtr NewAtoms	( char *name );
SpeciesPtr NewSpecies	( char *name );
ReactionPtr NewReaction	( char *label );
ThirdBodyPtr NewThirdBody( char *name );
DimensionPtr NewDimension( char *name, Double value, int kg, int m, 
	int s, int K, int mole );
ParameterPtr NewParameter( char *name );
int FindAtom	( LinkPtr link, void *var );
int FindSpeciesName( LinkPtr link, void *var );
int FindSpeciesNumber( LinkPtr link, void *var );
int FindSpecies	( LinkPtr link, void *var );
int FindReaction	( LinkPtr link, void *var );
int FindThirdBody( LinkPtr link, void *var );
int FindDimension( LinkPtr link, void *var );
int FindParameter( LinkPtr link, void *var );
ParameterPtr AddParameter( ListPtr list, char *name );
int PrintAtoms( LinkPtr link, void *var );
int PrintSpecies( LinkPtr link, void *var );
int PrintSpeciesNames( LinkPtr link, void *var );
void PrintSpeciesComposition( SpeciesPtr species, FILE *fp );
int PrintReaction( LinkPtr link, void *var );
int PrintThirdBody( LinkPtr link, void *var );
int PrintDimensionTable( LinkPtr link, void *var );
int PrintDimension( LinkPtr link, void *var );
int PrintParameter( LinkPtr link, void *var );
int PrintParameterFactors( LinkPtr link, void *var );
int PrintBroadening( LinkPtr link, void *var );
int PrintMagicRateCoeff( LinkPtr link, void *var );
#ifndef powerc
void UpperString	( char *string );
#endif
void LowerString	( char *string );
AtomsPtr AddAtom	( ListPtr list, char *name );
SpeciesPtr AddSpecies	( ListPtr list, char *name, IntVectorPtr composition );
ReactionPtr AddReaction	( ListPtr list, char *name );
void *AddReactionVoid( ListPtr list, void *var );
DimensionPtr AddDimension( ListPtr list, char *name, Double value, int kg, 
	int m, int s, int K, int mole );
ParameterPtr AddParameter( ListPtr list, char *name );
ThirdBodyPtr AddThirdBody( ListPtr list, char *name );
ReactionPtr InsertReaction( ListPtr list, char *name, int appendTo );
char *MatchAtom( char *yytext, int yyleng, int mode );
int ShiftStack( char *yytext, int yyleng, int *currentSize, char *stack );
int MatchItem( char *atom, int *currentSize );
void SkipWhites( char * label );
void ClearComposition( IntVectorPtr composition );
int FreeAtoms( LinkPtr link, void *var );
int FreeSpecies( LinkPtr link, void *var );
int FreeReaction( LinkPtr link, void *var );
int FreeThirdBody( LinkPtr link, void *var );
int FreeDimension( LinkPtr link, void *var );
int FreeParameter( LinkPtr link, void *var );
int FreeOneAtoms( LinkPtr link );
int FreeOneSpecies( LinkPtr link );
int FreeOneReaction( LinkPtr link );
int FreeOneThirdBody( LinkPtr link );
int FreeOneDimension( LinkPtr link );
int FreeOneParameter( LinkPtr link );
int DumpAtoms( LinkPtr link, void *var );
int DumpSpecies( LinkPtr link, void *var );
int DumpReaction( LinkPtr link, void *var );
int DumpThirdBody( LinkPtr link, void *var );
int DumpDimension( LinkPtr link, void *var );
void LineCounter( void );
void AdjustNumberOfLines( void );
void FillThirdBodyArray( int speciesNumber, Double Mcoefficient );
void FillOtherThirdBodyArray( Double coefficient );
char *CopyString( char *to, char *from );
void ScanError( char *message, char *text, Double value );
void InitDimension( ListPtr list );
int CheckReactionZeroEntry( LinkPtr link, void *var );
int CompareDoubleReaction( LinkPtr link, void *var );
int CompareReaction( ReactionPtr reaction1, ReactionPtr reaction2 );
int CheckStoichiometry( ReactionPtr reaction );
int SetOneOrderOfReaction( ReactionPtr reaction );
int CheckForwardBackward( void );
OptionsPtr NewOptions( void );
void FreeOptions( OptionsPtr command );
void PrintLists( void );
void CleanUp( void );
void InitAllowedExponents( void );
void ComputeLindemann();
int WriteFcFile( LinkPtr link, void *var );
void SetCounter( void );
void WriteOutput( void );
void ReadInput( void );
void CheckReactions( void );
void ComputeRateConstants( ReactionPtr forw, ReactionPtr backw );
int CheckForDoubles( void );
int CheckBroadening( void );
int CheckOneBroadening( LinkPtr link, void *var );
void WriteLatexTable( void );
HeaderPtr InitHeader( HeaderPtr header );
CounterPtr NewCounter( void );
void ReadThermoProperties( void );
void SetZeroProp(SpeciesPtr species);
int GetThermoProps( LinkPtr link, void *var );
int GetComposition( char *name, IntVectorPtr composition );
int MatchBuffer( ListPtr list, String buffer, int *currentLoc, String name, 
	IntVectorPtr composition );
void CreateSortedList( ListPtr charList, int number );
void ComputeOmegaDCoeff( SpeciesPtr species1, int j );
int ComputeDiffusionConstants( LinkPtr link, void *var );
IntVectorPtr *NewIntVectorArray( int number, int length );
void FreeIntVectorArray( IntVectorPtr *ptr, int number );
void FatalError( const char *s );
void Warning( const char *s );
int LookupSpecies( Pointer ptr, Pointer find );
void SortSpecies( void );
void SortSteadyStates( void );
int UpdateReactionList( LinkPtr link, void *var );
int UpdateSpeciesList( LinkPtr link, void *var );
int UpdateThirdBodyList( LinkPtr link, void *var );
void SaveLabels( void );
void SaveLabelsF77( void );
void SaveLabelsF90( void );
int SetOrderOfReaction( LinkPtr link, void *var );
int CheckMolarMass( LinkPtr link, void *var );
int CheckSigmaEps( LinkPtr link, void *var );
int ComputeMuCoeff( LinkPtr link, void *var );
int FindSpeciesInReaction( Pointer ptr, Pointer aux );
int GetOrderOfReaction( ReactionPtr reac, Flag infRequired );
Double GetEConverter( ReactionPtr reac, Flag infRequired );
Double GetAConverter( ReactionPtr reac, Flag infRequired );
Double GetETokmoleConverter( ReactionPtr reac, Flag infRequired );
Double GetATokmoleConverter( ReactionPtr reac, Flag infRequired );
void WriteRedMechReactions( void );
int PrintRedMechReaction( LinkPtr link, void *var );
int PrintRedMechProdRate( LinkPtr link, void *var );
int ReactionLabelsToUpper( LinkPtr link, void *aux );
void WriteMagicC( void );
void WriteMagicF77( void );
void WriteMagicF90( void );
void WriteMagicF77Vectorize( void );
void WriteScanManMech( void );
int PrintMagicRateCoeff( LinkPtr link, void *var );
int PrintMagicRateCoeffF77( LinkPtr link, void *var );
int PrintMagicRateCoeffF77Vectorize( LinkPtr link, void *var );
int PrintMagicRateCoeffF90( LinkPtr link, void *var );
int PrintMagicThirdBody( LinkPtr link, void *var );
int PrintMagicThirdBodyF77( LinkPtr link, void *var );
int PrintMagicThirdBodyF77Vectorize( LinkPtr link, void *var );
int PrintMagicThirdBodyF90( LinkPtr link, void *var );
int PrintMagicReacRate( LinkPtr link, void *var );
int PrintMagicReacRateF77( LinkPtr link, void *var );
int PrintMagicReacRateF77Vectorize( LinkPtr link, void *var );
int PrintMagicProdRate( LinkPtr link, void *var );
int PrintMagicProdRateF77( LinkPtr link, void *var );
int PrintMagicProdRateF77Vectorize( LinkPtr link, void *var );
int PrintMagicProdRateF90( LinkPtr link, void *var );
void PrintBroadeningSource( ReactionPtr reaction, char *redLabel, FILE *fp );
void PrintBroadeningSourceF77( ReactionPtr reaction, char *redLabel, FILE *fp );
void PrintBroadeningSourceF77Vectorize( ReactionPtr reaction, FILE *fp );
void PrintBroadeningSourceF90( ReactionPtr reaction, char *redLabel, FILE *fp );
int PrintDeclarationBroadeningFactorsF77( LinkPtr link, void *var );
int PrintDeclarationBroadeningFactorsF90( LinkPtr link, void *var );
void PrintOneRateCoeff( ReactionPtr reaction, ReactionRateType rateType, FILE *fp );
void PrintOneRateCoeffF77( ReactionPtr reaction, ReactionRateType rateType, FILE *fp );
void PrintOneRateCoeffF77Vectorize( ReactionPtr reaction, ReactionRateType rateType, FILE *fp );
void PrintOneRateCoeffF90( ReactionPtr reaction, ReactionRateType rateType, FILE *fp );
void ConvertCDoubleToF77Double( Double c, char *f );
int PrintMagicThermoDataF77( FILE *fp );
int PrintOneThermoDataF77( LinkPtr link, void *var );
int PrintMagicThermoDataF90( FILE *fp );
int PrintOneThermoDataF90( LinkPtr link, void *var );
int PrintMagicThermoDataF77Vectorize( FILE *fp );
int PrintMagicThermoData( FILE *fp );
int PrintOneThermoData( LinkPtr link, void *var );
int PrintOneThermoDataF77Vectorize( LinkPtr link, void *var );
int PrintMagicNamesF77( LinkPtr link, void *var );
int PrintMagicMolarMassF77( LinkPtr link, void *var );
int PrintMagicRateCoeffF77New( LinkPtr link, void *var );
void PrintArrCoeffsF77( ReactionPtr reaction, ReactionRateType rateType, Double *a,
						Double *n, Double *e );
int PrintMagicMolarMass( LinkPtr link, void *var );
int PrintMagicNames( LinkPtr link, void *var );
char *CSymbolLeadNum( const char *src );
char *CSymbol( const char *src );
Flag CompareThirdBody( LinkPtr link1, LinkPtr link2 );
char *GetSpeciesName( char *buff, char *name );
char *ShortSymbol( const char *src, char *buffer );
int PrintThermoFile( LinkPtr link, void *var );
void WriteChemkin( void );
int PrintChemkinAtoms( LinkPtr link, void *var );
int PrintChemkinSpecies( LinkPtr link, void *var );
int PrintChemkinReactions( LinkPtr link, void *var );
int PrintDeclarattionOfW( LinkPtr link, void *var );
int PrintMagicKOverEpsF77( LinkPtr link, void *var );
int PrintMagicMuCoeffF77( LinkPtr link, void *var );
int PrintMagicProdRateVectorize( LinkPtr link, void *var );
int PrintMagicNASAPolynomialsVectorize( FILE *fp );


/*	Flags for full kinetic or reduced chemistry models						*/
enum kChemFlag { kAll, kComputedSpecies, kSteadyStateSpecies };
enum ReactionRateType		{ kK, kK0, kKinf };
extern enum kChemFlag gChemFlag;

int	gLabel;
int	gLabelIncrement;

#ifdef qUseDump
#ifdef applec
#pragma dump "ScanMan.dump"
#endif
#endif

#endif
