/*
 *	ScanMan:	Written by Heinz Pitsch and Josef Goettgens
 *				Copyright 1991-92, Heinz Pitsch
 *				All rights reserved.
 */

#undef VECTOR
#undef SHOWTHERMOFILE
#undef DEBUGCOMPOSITION
#undef DEBUGCHOICE
#undef DEBUG
#undef qDebug
#define PETER
#undef JOSEF
#undef OLDWRONGLATEXTABLE
#undef REFERENCE
#undef ENUMLATEXTABLE
#define SUPPRESSCOMMENTS

#ifdef qUseDump
#pragma load "ScanMan.dump"
#else
#include "ScanMan.h"
#endif

#include "ReadThermoProps.h"
#ifdef JOSEF
#include "Mathematica.h"
#endif
#ifdef PETER
#include "Redux.h"
#endif

#define TROE

#include	<unistd.h>

#define SKIPOMPARSE 1

/*
	GLOBALS
*/

extern FILE *yyin;              /* input file of lex.yy.c    */
extern int yyparse( void );

ListPtr					gAtomsList = NULL;
ListPtr					gSpeciesList = NULL;
ListPtr					gReactionList = NULL;
ListPtr					gPAHReactionList = NULL;
ListPtr					gSootReactionList = NULL;
ListPtr 				gAllowAtomList = NULL;
ListPtr					gUsedThirdBodyList = NULL;
ListPtr					gThirdBodyList = NULL;
ListPtr					gDimensionList = NULL;
ListPtr					gAllowTempExpList = NULL;
ListPtr					gAllowReacExpList = NULL;
ListPtr					gUsedDimensionList = NULL;
ListPtr					gCommentList = NULL;
ListPtr					gLabelFollowingCommentList = NULL;
ListPtr					gFcNameList = NULL;
ListPtr					gFcContentsList = NULL;
ListPtr					gFcLineNumberList = NULL;
ListPtr					gSteadyStateList = NULL;

IntVectorPtr			gComposition = NULL;
OptionsPtr				gOptions = NULL;
HeaderPtr				gHeader = NULL;
CounterPtr				gCounter = NULL;
char					*gbuffer = NULL;  /* initialized in InitHeader */
ListFuncPtr				gSpeciesCompareFunc = FindSpecies;

/*  error reporting  */
char					*gProgram;
int						gNumberOfLines = 1;
char					gLineContents[LINELENGTH];
int						gtypeOfYylval = 0;
int						gErrorOut = FALSE;

extern int yyparse( void );

/*	Prototypes for local functions.											*/
static void ConflictWarning( const char *, const char *, const char * );
static void Initialize( void );
static void HandleCommandLine( int argc, char *argv[] );
static void Usage( int exit_code );
static char *DefaultThermoFile( const char* name );


int main( int argc, char *argv[] )
{
	gProgram = argv[0];
#	ifdef applec
	InitCursorCtl( NULL );
#	endif

	Initialize();
	HandleCommandLine( argc, argv );
	
	if ( gOptions->printUnits && !gOptions->inName ) {
	
		ListIterator( gDimensionList, PrintDimensionTable, gDimensionList );
		
	}
	else {
		if ( yyparse() ) {
			PrintLists();
			FatalError( "Parser signalled errors" );
		}
		if ( gErrorOut ) {
			PrintLists();
			FatalError( "Parser signalled errors" );
		}
/* hp */
		ReadThermoProperties();
		ListIterator( gSpeciesList, CheckMolarMass, NULL );
		ListIterator( gSpeciesList, CheckSigmaEps, NULL );
		ListIterator( gSpeciesList, ComputeMuCoeff, NULL );
		ListIterator( gSpeciesList, ComputeDiffusionConstants, gSpeciesList );
		CheckReactions();
		if ( gThirdBodyList->items ) {
			SortList( gThirdBodyList, CompareThirdBody );
		}
#ifdef PETER
		if ( gOptions->steadyStates ) {
			/*	Will only be called if  gOptions->CCode  is FALSE.				*/
			ReduceMechanism();
		}
#endif
		SortSteadyStates();		/* (mbo)	*/
		SaveLabels();
		if ( gOptions->latexTable ) WriteLatexTable();
/*		WriteScanManMech();*/
		PrintLists();

		/*	Take care of output.												*/
		if ( gOptions->steadyStates ){
			WriteRedMechReactions();
		}
		
		if ( gOptions->mathematica ){
			/*WriteMathematica();*/
		}
#ifdef JOSEF
		if ( gOptions->mathematica )
			WriteMathematicaFile();
#endif
		if ( gOptions->outFile && (gOptions->noBinaryOutput == FALSE) ) {
			ListIterator( gFcNameList, WriteFcFile, gFcContentsList );
			WriteOutput();
		}
		
		if ( gOptions->outputWriter ) {
			WriteMagicC();
		}

		if ( gOptions->outputWriterF77 ) {
			WriteMagicF77();
		}
		
		if ( gOptions->outputWriterF90 ) {
			WriteMagicF90();
		}
	
		if ( gOptions->vectorizeF77 ) {	/* gb */
			WriteMagicF77Vectorize();
		}

		WriteChemkin();


#ifdef JOSEF
#ifdef qCCode
		if ( gOptions->CCode ) {
			/*	If  gOptions->steadyStates  is TRUE as well, C output for 
				reduced mechanisms is generated.								*/
			GenerateCCode( (gOptions->steadyStates != NULL) );
		} else 
#endif
#endif
	}
	
	CleanUp();
	fprintf( stderr, "done\n" );
	return 0;
}

void WriteChemkin( void )
{
	char	fName[128];
	FILE	*fpMech;
	

	struct CHEMFILES {
		FILE	*fpTherm;
		FILE	*fpTrans;
	} chFiles;

	sprintf( fName, "%s.chmech", gOptions->base );
	if ( !( fpMech = fopen( fName, "w" ) ) ) {
		fprintf( stderr, "#error: unable to open file %s\n", fName );
		return;
	}
	sprintf( fName, "%s.chthermo", gOptions->base );
	if ( !( chFiles.fpTherm = fopen( fName, "w" ) ) ) {
		fprintf( stderr, "#error: unable to open file %s\n", fName );
		return;
	}
	sprintf( fName, "%s.chtrans", gOptions->base );
	if ( !( chFiles.fpTrans = fopen( fName, "w" ) ) ) {
		fprintf( stderr, "#error: unable to open file %s\n", fName );
		return;
	}

	fprintf( fpMech, "ELEMENTS\n" );
	ListIterator( gAtomsList, PrintChemkinAtoms, fpMech );
	fprintf( fpMech, "\nEND\n" );
	fprintf( fpMech, "SPECIES\n" );
	ListIterator( gSpeciesList, PrintChemkinSpecies, fpMech );
	fprintf( fpMech, "\nEND\n" );
	fprintf( fpMech, "REACTIONS\n" );
	ListIterator( gReactionList, PrintChemkinReactions, fpMech );
	fprintf( fpMech, "END\n" );


	fprintf( chFiles.fpTrans, "#SpecName      Dummy     eps/k     sigma     Dummy     Dummy     Dummy\n" );
	fprintf( chFiles.fpTherm, "THERMO\n" );
	fprintf( chFiles.fpTherm, "%10.0f%10.0f%10.0f\n", 300.0, 1000.0, 5000.0 );
	ListIterator( gSpeciesList, PrintThermoFile, &chFiles );
	
	fclose( fpMech );
	fclose( chFiles.fpTherm );
	fclose( chFiles.fpTrans );
}

void WriteMagicC( void )
{
	char	fName[128];
	FILE	*fp;
	
	sprintf( fName, "%sC.c", gOptions->base );
	if ( !( fp = fopen( fName, "w" ) ) ) {
		fprintf( stderr, "#error: unable to open file %s\n", fName );
		return;
	}

/*	fprintf( fp, "#include \"maininclude.h\"\n", gOptions->base );*/
	fprintf( fp, "#include \"%s.h\"\n\n", gOptions->base );

	fprintf( fp, "#include <stdio.h>\n" );
	fprintf( fp, "#include <math.h>\n" );
	fprintf( fp, "#include <stdlib.h>\n" );
	fprintf( fp, "#include <string.h>\n" );
	fprintf( fp, "#include <time.h>\n" );

	fprintf( fp, "static double GetLindRateCoeff( double temp, double pressure\n" );
	fprintf( fp, "\t\t\t, double k0, double kInf, double fc, double conc );\n\n" );

	fprintf( fp, "void ComputeProductionRates( double *cdot, double *w, double *k\n" );
	fprintf( fp, "\t\t\t, double *c, double *M, double temp, double pressure );\n" );

	fprintf( fp, "double MAX_C(double X1, double X2);\n" );

/* 	fprintf( fp, "void Calc_h_cp_i_T_300K( double *cp, double *h);\n\n" ); */
/* 	fprintf( fp, "\t\t\t, double k0, double kInf, double fc, double conc );\n\n" ); */
	fprintf( fp, "void ComputeProductionRates( double *cdot, double *w, double *k\n" );
	fprintf( fp, "\t\t\t, double *c, double *M, double temp, double pressure )\n" );
	fprintf( fp, "{\n" );
	fprintf( fp, "/*\n\tThis function computes rates of production cdot in [kmole/(m^3s)].\n" );
	fprintf( fp, "\tThe parameters w ( reaction rate ), k ( rate coefficient )\n" );
	fprintf( fp, "\tand M ( third body concentrations ) are just work space for this function.\n" );
	fprintf( fp, "\tc contains the concentrations of non steady state species in [kmole/m^3]\n" );
	fprintf( fp, "\tand is workspace for the steady state concentrations, which are computed\n" );
	fprintf( fp, "\tin this function.\n" );
	fprintf( fp, "\ttemp is the temperature in [K] and pressure is the pressure in [Pa].\n" );
	fprintf( fp, "\tCalled functions are 'GetLindRateCoeff', 'ComputeSteadyStates',\n" );
	fprintf( fp, "\t'CatchZero' and the functions that evaluate the broadening factors\n" );
	fprintf( fp, "\tof the Troe formulation of pressure dependant rate coefficients 'Fc*'\n" );
	fprintf( fp, "*/\n\n" );

	fprintf( fp, "\tint\tnSpec = %d;\n\tint\tnSpecIn = %d;\n",
		gCounter->species, gCounter->species - gCounter->steadyStates );

	fprintf( fp, "\tdouble	kTroe0, kTroeInf, fcTroe;\n" );
	fprintf( fp, "\tdouble	RGAS = 8314.34;\n" );
	fprintf( fp, "\tdouble	lgt = log( temp );\n" );
	fprintf( fp, "\tdouble	rt = RGAS * temp;\n\n" );
	ListIterator( gThirdBodyList, PrintMagicThirdBody, fp );
	fprintf( fp, "\n\n" );
	ListIterator( gReactionList, PrintMagicRateCoeff, fp );
	fprintf( fp, "\n\n" );
	if ( gCounter->steadyStates - gCounter->pahSpecies - gCounter->sootSpecies ) {
		fprintf( fp, "\tComputeSteadyStates( k, c, M );\n\n" );
		fprintf( fp, "\tdouble\tcTot = pressure / ( RGAS * temp );\n");
		fprintf( fp, "\tfor ( int iNow = nSpecIn; iNow < nSpec; ++iNow ) {\n");
		fprintf( fp, "\t\tif ( c[iNow] > cTot ) {\n");
		fprintf( fp, "\t\t\tc[iNow] = cTot;\n");
		fprintf( fp, "\t\t}\n");
		fprintf( fp, "\t}\n");
	}
	ListIterator( gReactionList, PrintMagicReacRate, fp );
	fprintf( fp, "\n\n" );
	ListIterator( gSpeciesList, PrintMagicProdRate, fp );
	fprintf( fp, "}\n\n" );
	
/* 	fprintf( fp, "double Calc_h_cp_i_T_300K( double *cp, double *h)\n" ); */
/* 	fprintf( fp, "{\n" ); */
/* 	fprintf( fp, "	ComputeThermoData( h, cp, 300.0 )\n" ); */
/* 	fprintf( fp, "	return;\n" ); */
/* 	fprintf( fp, "}\n\n" ); */

	fprintf( fp, "double GetLindRateCoeff( double temp, double pressure\n" );
	fprintf( fp, "				, double k0, double kInf\n" );
	fprintf( fp, "				, double fc, double conc )\n" );
	fprintf( fp, "{\n" );
	fprintf( fp, "	const double	R = 8314.34;   /* [J / kmole K] */\n" );
	fprintf( fp, "	double			Ntmp;\n" );
	fprintf( fp, "	double			kl;\n" );
	fprintf( fp, "	double			f;\n" );
/*	fprintf( fp, "	double			conc;\n" );*/
	fprintf( fp, "	double			cCoeff, dCoeff, log10kNull;\n\n" );
#ifdef TROE
	fprintf( fp, "	int				iTroe = 1;\n\n" );
#else
	fprintf( fp, "	int				iTroe = 0;\n\n" );
#endif
	fprintf( fp, "	if ( conc <= 0.0 ) {\n" );
	fprintf( fp, "		conc = pressure / ( R * temp );\n" );
	fprintf( fp, "	}\n" );


	fprintf( fp, "	Ntmp = 0.75 - 1.27 * log10( fc );\n" );
	fprintf( fp, "\tif ( iTroe ) {\n" );
	fprintf( fp, "\t\tcCoeff = - 0.4 - 0.67 * log10( fc );\n" );
	fprintf( fp, "\t\tdCoeff = 0.14;\n" );
	fprintf( fp, "\t\tk0 *= conc / MAX_C(kInf, 1.0e-60);\n" );
	fprintf( fp, "\t\tlog10kNull = log10( k0 );\n" );
	fprintf( fp, "\t\tf = ( log10kNull + cCoeff ) / ( Ntmp - dCoeff * ( log10kNull + cCoeff ) );\n" );
	fprintf( fp, "\t\tf = pow( fc, 1.0 / ( f * f + 1.0 ) );\n" );
	fprintf( fp, "\t\tkInf *= f * k0 / ( 1.0 + k0 );\n" );
	fprintf( fp, "\t}\n" );
	fprintf( fp, "\telse {\n" );
	fprintf( fp, "\t	k0 = k0 * conc / kInf;\n" );
	fprintf( fp, "\t	kl = k0 / ( 1.0 + k0 );\n" );
	fprintf( fp, "\t	f = log10( k0 ) / Ntmp;\n" );
	fprintf( fp, "\t	f = pow( fc, 1.0 / ( f * f + 1.0 ) );\n" );
	fprintf( fp, "\t	kInf = kInf * f * kl;\n\n" );
	fprintf( fp, "\t}\n" );
	fprintf( fp, "	return kInf;\n" );
	fprintf( fp, "\t\n" );

	fprintf( fp, "}\n\n" );

	fprintf( fp, "/*\n" );
	fprintf( fp, "double CatchZero( double a )\n" );
	fprintf( fp, "{\n" );
	fprintf( fp, "	return ( a == 0.0 ) ? 1.0e-20 : a;\n" );
	fprintf( fp, "}\n\n" );
	fprintf( fp, "*/\n" );

	fprintf( fp, "void GetMolarMass( double *W )\n" );
	fprintf( fp, "{\n" );
	ListIterator( gSpeciesList, PrintMagicMolarMass, fp );
	fprintf( fp, "}\n\n" );

	fprintf( fp, "void GetSpeciesNames( char **names )\n" );
	fprintf( fp, "{\n" );
	ListIterator( gSpeciesList, PrintMagicNames, fp );
	fprintf( fp, "}\n\n" );

	fprintf( fp, "\n\n" );

	PrintMagicThermoData( fp );

	fprintf( fp, "int GetNSpecies( void )\n" );
	fprintf( fp, "{\n" );
	fprintf( fp, "\treturn sEnd;\n" );
	fprintf( fp, "}\n\n" );

	fprintf( fp, "int GetNSpecs( void )\n" );
	fprintf( fp, "{\n" );
	fprintf( fp, "\treturn %i;\n",gCounter->species - gCounter->steadyStates);
	fprintf( fp, "}\n\n" );

	fprintf( fp, "int GetNReactions( void )\n" );
	fprintf( fp, "{\n" );
	fprintf( fp, "\treturn rEnd;\n" );
	fprintf( fp, "}\n\n" );

	fprintf( fp, "double MAX_C(double X1, double X2)\n" );
	fprintf( fp, "{\n" );
	fprintf( fp, "  return ( (X1 > X2) ? X1 : X2 );\n" );
	fprintf( fp, "}\n" );

	fclose( fp );
}

void WriteMagicF77Vectorize( void )
{
	char	fName[128];
	FILE	*fp;
	Flag	isReduced = ( gCounter->steadyStates - gCounter->pahSpecies - gCounter->sootSpecies 
				|| gHeader->globalMechanism ) ? TRUE : FALSE;
	
	sprintf( fName, "%sF.f", gOptions->base );
	if ( !( fp = fopen( fName, "w" ) ) ) {
		fprintf( stderr, "#error: unable to open file %s\n", fName );
		return;
	}
	ListIterator( gReactionList, ReactionLabelsToUpper, NULL );
	SaveLabelsF77();
	
	gLabel	=	10;
	
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C ======= %s =======\n", fName );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      SUBROUTINE PRODRATES(CDOT,W,KR,C,M,TEMP,PRESSURE,NGRIDPOINTS,\n" );
	fprintf( fp, "     &                     NREACTIONS,NSPECIES,NSPECS,NTHIRDBCONC,RGAS,\n" );
	fprintf( fp, "     &                     ISOTHERMAL)\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     THIS SUBROUTINE COMPUTES RATES OF PRODUCTION CDOT\n" );
	fprintf( fp, "C     IN [KMOLE/(M^3S)]. THE PARAMETERS W ( REACTION RATE ),\n" );
	fprintf( fp, "C     K ( RATE COEFFICIENT ) AND M ( THIRD BODY CONCENTRATIONS ) ARE\n" );
	fprintf( fp, "C     JUST WORK SPACE FOR THIS FUNCTION.\n" );
	fprintf( fp, "C     C CONTAINS THE CONCENTRATIONS OF NON STEADY STATE SPECIES IN\n" );
	fprintf( fp, "C     [KMOLE/M^3] AND IS WORKSPACE FOR THE STEADY STATE CONCENTRATIONS,\n" );
	fprintf( fp, "C     WHICH ARE COMPUTED IN THIS FUNCTION.\n" );
	fprintf( fp, "C     TEMP IS THE TEMPERATURE IN [K] AND\n" );
	fprintf( fp, "C     PRESSURE IS THE PRESSURE IN [PA].\n" );
	fprintf( fp, "C     CALLED FUNCTIONS ARE 'GETLINDRATECOEFF', 'COMPSTEADYSTATES',\n" );
	fprintf( fp, "C     'CTCHZERO'\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "      INCLUDE 'RIF_DIMENSION.h'\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INTEGER          NGRIDPOINTS, NREACTIONS, NSPECIES, NSPECS,\n" );
	fprintf( fp, "     &                 NTHIRDBCONC, ISOTHERMAL\n");
	fprintf( fp, "      DOUBLE PRECISION CDOT(0:(NGRIDPOINTS+1),NSPECS)\n");
	fprintf( fp, "      DOUBLE PRECISION W\n");
	ListIterator( gReactionList, PrintDeclarattionOfW, fp );
/*	fprintf( fp, "      DOUBLE PRECISION W(0:(NGRIDPOINTS+1),NREACTIONS)\n");*/
	fprintf( fp, "      DOUBLE PRECISION KR(0:(NGRIDPOINTS+1),NREACTIONS)\n");
	fprintf( fp, "      DOUBLE PRECISION C(0:(NGRIDPOINTS+1),NSPECIES)\n");
	fprintf( fp, "      DOUBLE PRECISION M(0:(NGRIDPOINTS+1),NTHIRDBCONC)\n");
	fprintf( fp, "      DOUBLE PRECISION TEMP(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION RGAS,PRESSURE\n");
#ifdef TROE
    	fprintf( fp, "      DOUBLE PRECISION NTMP,K0K,KLK,CONC,F,CCOEFF,DCOEFF,LGKNULL\n");
#else
    	fprintf( fp, "      DOUBLE PRECISION NTMP,K0K,KLK,CONC,F\n");
#endif
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INTEGER I,K\n" );
	fprintf( fp, "      DOUBLE PRECISION LT(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION RT(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION RTI(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION K0(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION KINF(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION FC(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      DO %d I=1,NSPECS\n", gLabel);
	fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);
	fprintf( fp, "        CDOT(K,I)=0.0D0\n");
	fprintf( fp, "%5d CONTINUE\n", gLabel);
	fprintf( fp, "C\n" );

	gLabel+=gLabelIncrement;
	ListIterator( gThirdBodyList, PrintMagicThirdBodyF77Vectorize, fp );
	fprintf( fp, "C\n" );

	gLabel+=gLabelIncrement;
	fprintf( fp, "      IF( ISOTHERMAL.NE.1 ) THEN\n" );
	fprintf( fp, "        DO %d K=1,NGRIDPOINTS\n", gLabel);
	fprintf( fp, "          LT(K)=DLOG(TEMP(K))\n");
	fprintf( fp, "          RT(K)=RGAS*TEMP(K)\n");
	fprintf( fp, "          RTI(K)=1.0D+00/RT(K)\n");
	fprintf( fp, "%5d   CONTINUE\n", gLabel);
	gLabel+=gLabelIncrement;
	fprintf( fp, "C\n" );
	fprintf( fp, "        DO %d I=1,NREACTIONS\n", gLabel);
	fprintf( fp, "        DO %d K=1,NGRIDPOINTS\n", gLabel);
	fprintf( fp, "          KR(K,I)=A(I)*DEXP(N(I)*LT(K)-E(I)*RTI(K))\n");
	fprintf( fp, "%5d   CONTINUE\n", gLabel);
	gLabel+=gLabelIncrement;
	fprintf( fp, "C\n" );
	ListIterator( gReactionList, PrintMagicRateCoeffF77Vectorize, fp );
	fprintf( fp, "C\n" );
	fprintf( fp, "      ENDIF\n" );
/*	ListIterator( gThirdBodyList, PrintMagicThirdBodyF77Vectorize, fp );*/
/*	fprintf( fp, "C\n" );*/
	if ( gCounter->steadyStates ) {
		fprintf( fp, "      CALL COMPSTEADYSTATES( KR, C, M, NGRIDPOINTS, NREACTIONS,\n" );
		fprintf( fp, "     $                             NSPECIES, NTHIRDBCONC )\nC\n" );
	}
/* gb 04.12.98 auskommentiert wegen neuer vektorisiert. Version */
/*	ListIterator( gReactionList, PrintMagicReacRateF77Vectorize, fp );*/
/*	fprintf( fp, "C\n" );*/
/*	ListIterator( gSpeciesList, PrintMagicProdRateF77Vectorize, fp );*/

	ListIterator( gReactionList, PrintMagicProdRateVectorize, fp );
	fprintf( fp, "      END\nC\nC\n" );
	


/*
	gLabel=10;
	
	fprintf( fp, "      SUBROUTINE GETLINDRATECOEFF( TEMP, PRESSURE,\n" );
	fprintf( fp, "     &    K0, KINF, FC, KR, NGRIDPOINTS, RGAS )\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      INTEGER NGRIDPOINTS\n" );
	fprintf( fp, "      DOUBLE PRECISION TEMP(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION K0(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION KINF(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION FC(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION KR(0:(NGRIDPOINTS+1))\n");
	fprintf( fp, "      DOUBLE PRECISION RGAS,PRESSURE\n");
	fprintf( fp, "C\n" );
	fprintf( fp, "      INTEGER K\n" );
	fprintf( fp, "      DOUBLE PRECISION NTMP\n" );
	fprintf( fp, "      DOUBLE PRECISION K0K,KLK\n" );
	fprintf( fp, "      DOUBLE PRECISION F\n" );
	fprintf( fp, "      DOUBLE PRECISION CONC\nC\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);
	fprintf( fp, "      CONC = PRESSURE / ( RGAS * TEMP(K) )\n" );
	fprintf( fp, "      NTMP = 0.75D+00 - 1.27D+00 * DLOG10( FC(K) )\n" );
	fprintf( fp, "      K0K = K0(K) * CONC / KINF(K)\n" );
	fprintf( fp, "      KLK = K0K / ( 1.0D+00 + K0K )\n" );
	fprintf( fp, "      F = DLOG10( K0K ) / NTMP\n" );
	fprintf( fp, "      F = FC(K) ** ( 1.0D+00 / ( F * F + 1.0D+00 ) )\n" );
	fprintf( fp, "      KR(K) = KINF(K) * F * KLK\n" );
	fprintf( fp, "%5d CONTINUE\n", gLabel);
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );
*/
/* CTCHZEROV */
	gLabel=10;
	
	fprintf( fp, "      SUBROUTINE CTCHZEROV( X, NGRIDPOINTS )\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INTEGER NGRIDPOINTS\n" );
	fprintf( fp, "      DOUBLE PRECISION X(0:(NGRIDPOINTS+1))\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INTEGER K\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);
	fprintf( fp, "      IF ( X(K).EQ.0.0D0 ) THEN\n" );
	fprintf( fp, "         X(K) = 1.0D-20\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "%5d CONTINUE\n", gLabel);
	fprintf( fp, "      RETURN\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* SOLVQDRTV */
	gLabel=10;
	
	fprintf( fp, "      SUBROUTINE SOLVQDRTV( A, B, C, Y, NGRIDPOINTS )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     SOLVES THE QUADRATIC EQUATION A x^2 + B x - C = 0\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INTEGER NGRIDPOINTS\n" );
	fprintf( fp, "      DOUBLE PRECISION A(0:(NGRIDPOINTS+1))\n" );
	fprintf( fp, "      DOUBLE PRECISION B(0:(NGRIDPOINTS+1))\n" );
	fprintf( fp, "      DOUBLE PRECISION C(0:(NGRIDPOINTS+1))\n" );
	fprintf( fp, "      DOUBLE PRECISION Y(0:(NGRIDPOINTS+1))\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INTEGER K\n" );
	fprintf( fp, "      DOUBLE PRECISION RAD\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      DO %d K=1,NGRIDPOINTS\n", gLabel);
	fprintf( fp, "      IF ( A(K).EQ.0.0D+00 ) THEN\n" );
	fprintf( fp, "        IF ( B(K).EQ.0.0D+00 ) THEN\n" );
	fprintf( fp, "          WRITE(*,*) '#WARNING: INVALID ARGUMENTS IN '\n" );
	fprintf( fp, "     $                ,'SUBROUTINE SOLVQDRTV '\n" );
	fprintf( fp, "          Y(K) = 0.0D+00\n" );
	fprintf( fp, "        ELSE\n" );
	fprintf( fp, "          Y(K) = C(K) / B(K)\n" );
	fprintf( fp, "        END IF\n" );
	fprintf( fp, "      ELSE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "        B(K) = B(K) / A(K)\n" );
	fprintf( fp, "        C(K) = C(K) / A(K)\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "        RAD = 0.25D+00 * B(K) * B(K) + C(K)\n" );
	fprintf( fp, "        IF ( RAD.GE.0.0D+00 ) THEN\n" );
	fprintf( fp, "          RAD = DSQRT( RAD ) + 0.5D+00 * B(K)\n" );
	fprintf( fp, "          IF ( RAD.NE.0.0D0 ) THEN\n" );
	fprintf( fp, "            Y(K) = C(K) / RAD\n" );
	fprintf( fp, "          ELSE\n" );
	fprintf( fp, "            Y(K)=-0.5D+00*B(K)+DSQRT(0.25D+00*B(K)*B(K)+C(K))\n" );
	fprintf( fp, "          ENDIF\n" );
	fprintf( fp, "        ELSE\n" );
	fprintf( fp, "          Y(K)=0.0D+00\n" );
	fprintf( fp, "        ENDIF\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "%5d CONTINUE\n", gLabel);
	fprintf( fp, "      RETURN\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* CTCHZERO */
	fprintf( fp, "      DOUBLE PRECISION FUNCTION CTCHZERO( X )\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "      DOUBLE PRECISION X\n" );
	fprintf( fp, "      IF ( X.EQ.0.0D0 ) THEN\n" );
	fprintf( fp, "         X = 1.0D-20\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "      CTCHZERO = X\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );


/* SOLVQDRT */
	fprintf( fp, "      DOUBLE PRECISION FUNCTION SOLVQDRT( A, B, C )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     SOLVES THE QUADRATIC EQUATION A x^2 + B x - C = 0\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      DOUBLE PRECISION RAD, A, B, C\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IF ( A.EQ.0.0D0 ) THEN\n" );
	fprintf( fp, "        IF ( B.EQ.0.0D0 ) THEN\n" );
	fprintf( fp, "          WRITE(*,*) '#WARNING: INVALID ARGUMENTS IN FUNCTION SOLVQDRT '\n" );
	fprintf( fp, "          SOLVQDRT = 0.0\n" );
	fprintf( fp, "        ELSE\n" );
	fprintf( fp, "          SOLVQDRT = C / B\n" );
	fprintf( fp, "        END IF\n" );
	fprintf( fp, "      ELSE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "        B = B / A\n" );
	fprintf( fp, "        C = C / A\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "        RAD = 0.25D0 * B * B + C\n" );
	fprintf( fp, "        IF ( RAD.GE.0.0D0 ) THEN\n" );
	fprintf( fp, "          RAD = DSQRT( RAD ) + 0.5 * B\n" );
	fprintf( fp, "          IF ( RAD.NE.0.0D0 ) THEN\n" );
	fprintf( fp, "            SOLVQDRT = C / RAD\n" );
	fprintf( fp, "          ELSE\n" );
	fprintf( fp, "            SOLVQDRT = -0.5D0 * B + DSQRT( 0.25 * B * B + C )\n" );
	fprintf( fp, "          ENDIF\n" );
	fprintf( fp, "        ELSE\n" );
	fprintf( fp, "          SOLVQDRT = 0.0D0\n" );
	fprintf( fp, "        ENDIF\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETMOLARMASS */
	fprintf( fp, "      SUBROUTINE GETMOLARMASS( MM )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     FILLS 'MM' WITH SPECIES MOLAR MASS IN KG/KMOLE\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      DOUBLE PRECISION MM(%d)\n"
				, gCounter->species );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	ListIterator( gSpeciesList, PrintMagicMolarMassF77, fp );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETSPECIESNAMES */
	fprintf( fp, "      SUBROUTINE GETSPECIESNAMES( NAMES )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	FILLS 'NAMES' WITH SPECIES IDENTIFIER/KG\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      CHARACTER *20 NAMES(%d)\n"
				, gCounter->species );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	ListIterator( gSpeciesList, PrintMagicNamesF77, fp );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );
/* GETNSPECIES */
	fprintf( fp, "      SUBROUTINE GETNSPECIES( NSPECIES )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	FILLS 'NSPECIES' WITH NUMBER OF SPECIES \n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      INTEGER NSPECIES\n");
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	fprintf( fp, "      NSPECIES = SEND - 1\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETNREACTIONS */
	fprintf( fp, "      SUBROUTINE GETNREACTIONS( NREACTIONS )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	FILLS 'NREACTIONS' WITH NUMBER OF REACTIONS \n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      INTEGER NREACTIONS\n");
	fprintf( fp, "C\n" );
	fprintf( fp, "      NREACTIONS = %d\n", gCounter->reactions );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETNSPECS */
	fprintf( fp, "      SUBROUTINE GETNSPECS(NSPECIES_NONS)\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     RETURNS THE NUMBER OF THE NON STEADY STATE SPECIES\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      INTEGER NSPECIES_NONS\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      NSPECIES_NONS = %d\n", gCounter->species - gCounter->steadyStates );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETMUCOEFF */
	fprintf( fp, "      SUBROUTINE GETMUCOEFF( MUCOEFF )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	FILLS 'MUCOEFF' WITH MUECOEFF IN KG/(M*S)\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      DOUBLE PRECISION MUCOEFF(%d)\n", gCounter->species );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	ListIterator( gSpeciesList, PrintMagicMuCoeffF77, fp );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" ); 

/* GETKOVEREPS */
	fprintf( fp, "      SUBROUTINE GETKOVEREPS( KOVEREPS )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	    FILLS 'KOVEREPS' WITH KOVEREPS\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      DOUBLE PRECISION KOVEREPS(%d)\n", gCounter->species );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	ListIterator( gSpeciesList, PrintMagicKOverEpsF77, fp );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* COMPSTEADYSTATES */
	fprintf( fp, "      SUBROUTINE COMPSTEADYSTATES( KR, C, M, NGRIDPOINTS, NREACTIONS,\n" );
	fprintf( fp, "     $                             NSPECIES, NTHIRDBCONC )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS FROM\n" );
	fprintf( fp, "C     THE CONCENTRATIONS OF COMPUTED SPECIES AND RATE COEFFICIENTS.\n" );
	fprintf( fp, "C     CONCENTRATIONS OF COMPUTED SPECIES MAY NOT BE ALTERED.\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INTEGER          NGRIDPOINTS, NREACTIONS, NSPECIES, NTHIRDBCONC\n" );
	fprintf( fp, "      DOUBLE PRECISION KR(0:(NGRIDPOINTS+1),NREACTIONS)\n" );
	fprintf( fp, "      DOUBLE PRECISION C(0:(NGRIDPOINTS+1),NSPECIES)\n" );
	fprintf( fp, "      DOUBLE PRECISION M(0:(NGRIDPOINTS+1),NTHIRDBCONC)\n" );
	fprintf( fp, "      INTEGER          K\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	if ( isReduced ) {
		fprintf( fp, "      INCLUDE '%sF.ss'\n", gOptions->base );
	}
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );
	PrintMagicThermoDataF77Vectorize( fp );

/* GETNTHIRDBCONC */
	fprintf( fp, "      INTEGER FUNCTION GETNTHIRDBCONC()\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     RETURNS THE NUMBER OF THIRD BODIES\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      GETNTHIRDBCONC = %d\n", gCounter->thirdBodies );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETNASAPOLYNOMIALS */
	PrintMagicNASAPolynomialsVectorize( fp );
	
	fclose( fp );
}

void WriteMagicF77( void )
{
	char	fName[128];
	FILE	*fp;
	Flag	isReduced = ( gCounter->steadyStates - gCounter->pahSpecies - gCounter->sootSpecies 
				|| gHeader->globalMechanism ) ? TRUE : FALSE;
	
	sprintf( fName, "%sF.f", gOptions->base );
	if ( !( fp = fopen( fName, "w" ) ) ) {
		fprintf( stderr, "#error: unable to open file %s\n", fName );
		return;
	}
	ListIterator( gReactionList, ReactionLabelsToUpper, NULL );
	SaveLabelsF77();
	
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C ======= %s =======\n", fName );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      SUBROUTINE PRODRATES( CDOT, W, K, C, M, TEMP,\n" );
	fprintf( fp, "     &    PRESSURE )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     THIS SUBROUTINE COMPUTES RATES OF PRODUCTION CDOT\n" );
	fprintf( fp, "C     IN [KMOLE/(M^3S)]. THE PARAMETERS W ( REACTION RATE ),\n" );
	fprintf( fp, "C     K ( RATE COEFFICIENT ) AND M ( THIRD BODY CONCENTRATIONS ) ARE\n" );
	fprintf( fp, "C     JUST WORK SPACE FOR THIS FUNCTION.\n" );
	fprintf( fp, "C     C CONTAINS THE CONCENTRATIONS OF NON STEADY STATE SPECIES IN\n" );
	fprintf( fp, "C     [KMOLE/M^3] AND IS WORKSPACE FOR THE STEADY STATE CONCENTRATIONS,\n" );
	fprintf( fp, "C     WHICH ARE COMPUTED IN THIS FUNCTION.\n" );
	fprintf( fp, "C     TEMP IS THE TEMPERATURE IN [K] AND\n" );
	fprintf( fp, "C     PRESSURE IS THE PRESSURE IN [PA].\n" );
	fprintf( fp, "C     CALLED FUNCTIONS ARE 'GETLINDRATECOEFF', 'COMPSTEADYSTATES',\n" );
	fprintf( fp, "C     'CTCHZERO'\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "      DOUBLE PRECISION CDOT(%d), W(%d), K(%d),\n",
		gCounter->species - gCounter->steadyStates, gCounter->reactions, gCounter->reactions );
	fprintf( fp, "     &    C(%d), M(%d), TEMP, PRESSURE\n",
		gCounter->species, gCounter->thirdBodies);
	fprintf( fp, "      INTEGER I\n" );
	fprintf( fp, "      DOUBLE PRECISION GETLINDRATECOEFF, LT\n" );
	fprintf( fp, "      DOUBLE PRECISION RGAS, RT, CONCDEFAULT\n" );
	fprintf( fp, "      PARAMETER ( RGAS = 8314.34, CONCDEFAULT = -1.0 )\n" );
	ListIterator( gReactionList, PrintDeclarationBroadeningFactorsF77, fp );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	fprintf( fp, "      LT = DLOG( TEMP )\n" );
	fprintf( fp, "      RT = RGAS * TEMP \n" );
	fprintf( fp, "C\n" );
	ListIterator( gThirdBodyList, PrintMagicThirdBodyF77, fp );
	fprintf( fp, "\nC\nC\n" );
#ifdef VECTOR
	fprintf( fp, "C\n" );
	fprintf( fp, "      DO 100 I=1,%d\n", gCounter->reactions );
	fprintf( fp, "        K(I)=A(I)*DEXP(N(I)*LT-E(I)/RT)\n" );
	fprintf( fp, " 100  CONTINUE\n" );
	fprintf( fp, "C\n" );

	ListIterator( gReactionList, PrintMagicRateCoeffF77New, fp );
#else
	ListIterator( gReactionList, PrintMagicRateCoeffF77, fp );
#endif
	fprintf( fp, "C\n" );
	if ( gCounter->steadyStates ) {
		fprintf( fp, "      CALL COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )\nC\n" );
	}
	ListIterator( gReactionList, PrintMagicReacRateF77, fp );
	fprintf( fp, "C\n" );
	ListIterator( gSpeciesList, PrintMagicProdRateF77, fp );
	fprintf( fp, "\n      END\nC\nC\n" );
	
	fprintf( fp, "      DOUBLE PRECISION FUNCTION GETLINDRATECOEFF( TEMP, PRESSURE,\n" );
	fprintf( fp, "     &    K0, KINF, FC, CONCIN )\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "      INTEGER ITROE\n" );
	fprintf( fp, "      DOUBLE PRECISION TEMP, PRESSURE, K0, KINF, FC\n" );
	fprintf( fp, "      DOUBLE PRECISION R\n" );
	fprintf( fp, "      PARAMETER (R = 8314.34)\n" );
	fprintf( fp, "C     [J / kmole K]\n" );
	fprintf( fp, "      DOUBLE PRECISION NTMP,CCOEFF,DCOEFF,LGKNULL\n" );
	fprintf( fp, "      DOUBLE PRECISION KL\n" );
	fprintf( fp, "      DOUBLE PRECISION F\n" );
	fprintf( fp, "      DOUBLE PRECISION CONC, CONCIN\nC\n" );
#ifdef TROE
	fprintf( fp, "      ITROE = 1\n" );
#else
	fprintf( fp, "      ITROE = 0\n" );
#endif
	fprintf( fp, "      IF (CONCIN.GT.0.0) THEN\n" );
	fprintf( fp, "        CONC = CONCIN\n" );
	fprintf( fp, "      ELSE\n" );
	fprintf( fp, "        CONC = PRESSURE / ( R * TEMP )\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "      NTMP = 0.75 - 1.27 * DLOG10( FC )\n" );

	fprintf( fp, "      IF (ITROE.EQ.1) THEN\n" );
	fprintf( fp, "        CCOEFF = - 0.4 - 0.67 * DLOG10( FC )\n" );
	fprintf( fp, "        DCOEFF = 0.14\n" );
	fprintf( fp, "        K0 = K0 * CONC / MAX(KINF, 1.0D-60)\n" );
	fprintf( fp, "        LGKNULL = DLOG10( K0 )\n" );
	fprintf( fp, "        F=(LGKNULL+CCOEFF)/(NTMP-DCOEFF*(LGKNULL+CCOEFF))\n" );
	fprintf( fp, "        F = FC**(1.0 / ( F * F + 1.0 ))\n" );
	fprintf( fp, "        GETLINDRATECOEFF = KINF * F * K0 / ( 1.0 + K0 )\n" );
	fprintf( fp, "      ELSE\n" );
	fprintf( fp, "        K0 = K0 * CONC / KINF\n" );
	fprintf( fp, "        KL = K0 / ( 1.0 + K0 )\n" );
	fprintf( fp, "        F = DLOG10( K0 ) / NTMP\n" );
	fprintf( fp, "        F = FC ** ( 1.0 / ( F * F + 1.0 ) )\n" );
	fprintf( fp, "        GETLINDRATECOEFF = KINF * F * KL\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

	fprintf( fp, "      DOUBLE PRECISION FUNCTION CTCHZERO( X )\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "      DOUBLE PRECISION X\n" );
	fprintf( fp, "      IF ( X.EQ.0.0D0 ) THEN\n" );
	fprintf( fp, "         X = 1.0D-20\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "      CTCHZERO = X\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* SOLVQDRT */
	fprintf( fp, "      DOUBLE PRECISION FUNCTION SOLVQDRT( A, B, C )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	SOLVES THE QUADRATIC EQUATION A x^2 + B x - C = 0\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      DOUBLE PRECISION RAD, A, B, C\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IF ( A.EQ.0.0D0 ) THEN\n" );
	fprintf( fp, "        IF ( B.EQ.0.0D0 ) THEN\n" );
	fprintf( fp, "          WRITE(*,*) '#WARNING: INVALID ARGUMENTS IN FUNCTION SOLVQDRT '\n" );
	fprintf( fp, "          SOLVQDRT = 0.0\n" );
	fprintf( fp, "        ELSE\n" );
	fprintf( fp, "          SOLVQDRT = C / B\n" );
	fprintf( fp, "        END IF\n" );
	fprintf( fp, "      ELSE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "        B = B / A\n" );
	fprintf( fp, "        C = C / A\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "        RAD = 0.25D0 * B * B + C\n" );
	fprintf( fp, "        IF ( RAD.GE.0.0D0 ) THEN\n" );
	fprintf( fp, "          RAD = DSQRT( RAD ) + 0.5 * B\n" );
	fprintf( fp, "          IF ( RAD.NE.0.0D0 ) THEN\n" );
	fprintf( fp, "            SOLVQDRT = C / RAD\n" );
	fprintf( fp, "          ELSE\n" );
	fprintf( fp, "            SOLVQDRT = -0.5D0 * B + DSQRT( 0.25 * B * B + C )\n" );
	fprintf( fp, "          ENDIF\n" );
	fprintf( fp, "        ELSE\n" );
	fprintf( fp, "          SOLVQDRT = 0.0D0\n" );
	fprintf( fp, "        ENDIF\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETMOLARMASS */
	fprintf( fp, "      SUBROUTINE GETMOLARMASS( MM )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	FILLS 'MM' WITH SPECIES MOLAR MASS IN KG/KMOLE\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      DOUBLE PRECISION MM(%d)\n"
				, gCounter->species );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	ListIterator( gSpeciesList, PrintMagicMolarMassF77, fp );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETSPECIESNAMES */
	fprintf( fp, "      SUBROUTINE GETSPECIESNAMES( NAMES )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	FILLS 'NAMES' WITH SPECIES IDENTIFIER/KG\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      CHARACTER *20 NAMES(%d)\n"
				, gCounter->species );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	ListIterator( gSpeciesList, PrintMagicNamesF77, fp );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETNSPECIES */
	fprintf( fp, "      SUBROUTINE GETNSPECIES( NSPECIES )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	FILLS 'NSPECIES' WITH NUMBER OF SPECIES \n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      INTEGER NSPECIES\n");
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	fprintf( fp, "      NSPECIES = SEND - 1\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETNREACTIONS */
	fprintf( fp, "      SUBROUTINE GETNREACTIONS( NREACTIONS )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	FILLS 'NREACTIONS' WITH NUMBER OF REACTIONS \n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      INTEGER NREACTIONS\n");
	fprintf( fp, "C\n" );
	fprintf( fp, "      NREACTIONS = %d\n", gCounter->reactions );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETNSPECS */
	fprintf( fp, "      SUBROUTINE GETNSPECS(NSPECIES_NONS)\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     RETURNS THE NUMBER OF THE NON STEADY STATE SPECIES\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      INTEGER NSPECIES_NONS\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      NSPECIES_NONS = %d\n", gCounter->species - gCounter->steadyStates );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );

/* GETMUCOEFF */
	fprintf( fp, "      SUBROUTINE GETMUCOEFF( MUCOEFF )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	FILLS 'MUCOEFF' WITH MUECOEFF IN KG/(M*S)\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      DOUBLE PRECISION MUCOEFF(%d)\n", gCounter->species );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	ListIterator( gSpeciesList, PrintMagicMuCoeffF77, fp );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" ); 

/* GETKOVEREPS */
	fprintf( fp, "      SUBROUTINE GETKOVEREPS( KOVEREPS )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C	    FILLS 'KOVEREPS' WITH KOVEREPS\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "C\n" );
	fprintf( fp, "      DOUBLE PRECISION KOVEREPS(%d)\n", gCounter->species );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	fprintf( fp, "C\n" );
	ListIterator( gSpeciesList, PrintMagicKOverEpsF77, fp );
	fprintf( fp, "C\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );


/* COMPSTEADYSTATES */
	fprintf( fp, "      SUBROUTINE COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "C     THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS FROM\n" );
	fprintf( fp, "C     THE CONCENTRATIONS OF COMPUTED SPECIES AND RATE COEFFICIENTS.\n" );
	fprintf( fp, "C     CONCENTRATIONS OF COMPUTED SPECIES MAY NOT BE ALTERED.\n" );
	fprintf( fp, "C-----------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	if ( isReduced ) {
		fprintf( fp, "      DOUBLE PRECISION CTCHZERO,SOLVQDRT\n" );
	}
	fprintf( fp, "      INTEGER NSPECIN, NSPEC, INOW\n" );
	fprintf( fp, "      DOUBLE PRECISION K(%d), C(%d), M(%d)\n",
		gCounter->reactions, gCounter->species, gCounter->thirdBodies );
	fprintf( fp, "      DOUBLE PRECISION TEMP, PRESSURE, CTOT\n");
	fprintf( fp, "      DOUBLE PRECISION R\n" );
	fprintf( fp, "      PARAMETER (R = 8314.34)\n" );
	fprintf( fp, "      INCLUDE '%sF.h'\n", gOptions->base );
	if ( isReduced ) {
		fprintf( fp, "      INCLUDE '%sF.ss'\n", gOptions->base );
		fprintf( fp, "      NSPECIN = %d\n",gCounter->species - gCounter->steadyStates );
		fprintf( fp, "      NSPEC = %d\n",gCounter->species );
		fprintf( fp, "      CTOT = PRESSURE / ( R * TEMP )\n");
		fprintf( fp, "      DO INOW=NSPECIN+1,NSPEC\n");
		fprintf( fp, "        IF ( C(INOW).gt.CTOT ) THEN\n");
		fprintf( fp, "          C(INOW) = CTOT\n");
		fprintf( fp, "        END IF\n");
		fprintf( fp, "      END DO\n");
	}
	fprintf( fp, "      END\n" );
	fprintf( fp, "C\nC\n" );
	PrintMagicThermoDataF77( fp );
	
	fclose( fp );
}

void WriteMagicF90( void )
{
	char	fName[128];
	FILE	*fp;
	Flag	isReduced = ( gCounter->steadyStates - gCounter->pahSpecies - gCounter->sootSpecies 
				|| gHeader->globalMechanism ) ? TRUE : FALSE;
	
	sprintf( fName, "%sF.f90", gOptions->base );
	if ( !( fp = fopen( fName, "w" ) ) ) {
		fprintf( stderr, "#error: unable to open file %s\n", fName );
		return;
	}
	ListIterator( gReactionList, ReactionLabelsToUpper, NULL );
	SaveLabelsF90();
	
	fprintf( fp, "!----------------------------------------------------------\n" );
	fprintf( fp, "! ======= %s =======\n", fName );
	fprintf( fp, "!----------------------------------------------------------\n" );
	fprintf( fp, "subroutine PRODRATES( CDOT, W, K, C, M, TEMP, PRESSURE)\n" );
	fprintf( fp, "!----------------------------------------------------------\n" );
	fprintf( fp, "!     THIS SUBROUTINE COMPUTES RATES OF PRODUCTION CDOT\n" );
	fprintf( fp, "!     IN [KMOLE/(M^3S)]. THE PARAMETERS W ( REACTION RATE ),\n" );
	fprintf( fp, "!     K ( RATE COEFFICIENT ) AND M ( THIRD BODY CONCENTRATIONS ) ARE\n" );
	fprintf( fp, "!     JUST WORK SPACE FOR THIS FUNCTION.\n" );
	fprintf( fp, "!     C CONTAINS THE CONCENTRATIONS OF NON STEADY STATE SPECIES IN\n" );
	fprintf( fp, "!     [KMOLE/M^3] AND IS WORKSPACE FOR THE STEADY STATE \n" );
	fprintf( fp, "!     CONCENTRATIONS, WHICH ARE COMPUTED IN THIS FUNCTION.\n" );
	fprintf( fp, "!     TEMP IS THE TEMPERATURE IN [K] AND\n" );
	fprintf( fp, "!     PRESSURE IS THE PRESSURE IN [PA].\n" );
	fprintf( fp, "!     CALLED FUNCTIONS ARE 'GETLINDRATECOEFF', 'COMPSTEADYSTATES',\n" );
	fprintf( fp, "!     'CTCHZERO'\n" );
	fprintf( fp, "!----------------------------------------------------------\n" );
	fprintf( fp, "      implicit none\n" );
	fprintf( fp, "      include '%sF90.h'\n", gOptions->base );
	fprintf( fp, "      real(DP) :: CDOT(%d), W(%d), K(%d), &\n",
		gCounter->species - gCounter->steadyStates, gCounter->reactions, gCounter->reactions );
	fprintf( fp, "      C(%d), M(%d), TEMP, PRESSURE\n",
		gCounter->species, gCounter->thirdBodies);
	fprintf( fp, "      integer ::  I\n" );
	fprintf( fp, "      real(DP) :: GETLINDRATECOEFF, LT, RT\n" );
	fprintf( fp, "      real(DP), parameter ::  RGAS = 8314.34, CONCDEFAULT = -1.0 \n" );
	ListIterator( gReactionList, PrintDeclarationBroadeningFactorsF90, fp );
	fprintf( fp, "\n" );
	fprintf( fp, "      LT = DLOG( TEMP )\n" );
	fprintf( fp, "      RT = RGAS * TEMP \n" );
	fprintf( fp, "\n" );
	ListIterator( gThirdBodyList, PrintMagicThirdBodyF90, fp );
	fprintf( fp, "\n\n\n" );
#ifdef VECTOR
	fprintf( fp, "\n" );
	fprintf( fp, "      DO 100 I=1,%d\n", gCounter->reactions );
	fprintf( fp, "        K(I)=A(I)*DEXP(N(I)*LT-E(I)/RT)\n" );
	fprintf( fp, " 100  CONTINUE\n" );
	fprintf( fp, "\n" );

	ListIterator( gReactionList, PrintMagicRateCoeffF77New, fp );
#else
	ListIterator( gReactionList, PrintMagicRateCoeffF90, fp );
#endif
	fprintf( fp, "\n" );
	if ( gCounter->steadyStates ) {
		fprintf( fp, "      CALL COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )\n\n" );
	}
	ListIterator( gReactionList, PrintMagicReacRateF77, fp );
	fprintf( fp, "\n" );
	ListIterator( gSpeciesList, PrintMagicProdRateF90, fp );
	fprintf( fp, "\n      END\n\n\n" );
	
	fprintf( fp, "      DOUBLE PRECISION FUNCTION GETLINDRATECOEFF( TEMP, PRESSURE, &\n" );
	fprintf( fp, "\t   K0, KINF, FC, CONCIN )\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "      INTEGER ITROE\n" );
	fprintf( fp, "      DOUBLE PRECISION TEMP, PRESSURE, K0, KINF, FC\n" );
	fprintf( fp, "      DOUBLE PRECISION R\n" );
	fprintf( fp, "      PARAMETER (R = 8314.34)\n" );
	fprintf( fp, "!     [J / kmole K]\n" );
	fprintf( fp, "      DOUBLE PRECISION NTMP,CCOEFF,DCOEFF,LGKNULL\n" );
	fprintf( fp, "      DOUBLE PRECISION KL\n" );
	fprintf( fp, "      DOUBLE PRECISION F\n" );
	fprintf( fp, "      DOUBLE PRECISION CONC, CONCIN\n\n" );
#ifdef TROE
	fprintf( fp, "      ITROE = 1\n" );
#else
	fprintf( fp, "      ITROE = 0\n" );
#endif
	fprintf( fp, "      IF (CONCIN.GT.0.0) THEN\n" );
	fprintf( fp, "        CONC = CONCIN\n" );
	fprintf( fp, "      ELSE\n" );
	fprintf( fp, "        CONC = PRESSURE / ( R * TEMP )\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "      NTMP = 0.75 - 1.27 * DLOG10( FC )\n" );

	fprintf( fp, "      IF (ITROE.EQ.1) THEN\n" );
	fprintf( fp, "        CCOEFF = - 0.4 - 0.67 * DLOG10( FC )\n" );
	fprintf( fp, "        DCOEFF = 0.14\n" );
	fprintf( fp, "        K0 = K0 * CONC / MAX(KINF, 1.0D-60)\n" );
	fprintf( fp, "        LGKNULL = DLOG10( K0 )\n" );
	fprintf( fp, "        F=(LGKNULL+CCOEFF)/(NTMP-DCOEFF*(LGKNULL+CCOEFF))\n" );
	fprintf( fp, "        F = FC**(1.0 / ( F * F + 1.0 ))\n" );
	fprintf( fp, "        GETLINDRATECOEFF = KINF * F * K0 / ( 1.0 + K0 )\n" );
	fprintf( fp, "      ELSE\n" );
	fprintf( fp, "        K0 = K0 * CONC / KINF\n" );
	fprintf( fp, "        KL = K0 / ( 1.0 + K0 )\n" );
	fprintf( fp, "        F = DLOG10( K0 ) / NTMP\n" );
	fprintf( fp, "        F = FC ** ( 1.0 / ( F * F + 1.0 ) )\n" );
	fprintf( fp, "        GETLINDRATECOEFF = KINF * F * KL\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "\n\n" );

	fprintf( fp, "      DOUBLE PRECISION FUNCTION CTCHZERO( X )\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "      DOUBLE PRECISION X\n" );
	fprintf( fp, "      IF ( X.EQ.0.0D0 ) THEN\n" );
	fprintf( fp, "         X = 1.0D-20\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "      CTCHZERO = X\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "\n\n" );

/* SOLVQDRT */
	fprintf( fp, "      DOUBLE PRECISION FUNCTION SOLVQDRT( A, B, C )\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "!	SOLVES THE QUADRATIC EQUATION A x^2 + B x - C = 0\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      DOUBLE PRECISION RAD, A, B, C\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      IF ( A.EQ.0.0D0 ) THEN\n" );
	fprintf( fp, "        IF ( B.EQ.0.0D0 ) THEN\n" );
	fprintf( fp, "          WRITE(*,*) '#WARNING: INVALID ARGUMENTS IN FUNCTION SOLVQDRT '\n" );
	fprintf( fp, "          SOLVQDRT = 0.0\n" );
	fprintf( fp, "        ELSE\n" );
	fprintf( fp, "          SOLVQDRT = C / B\n" );
	fprintf( fp, "        END IF\n" );
	fprintf( fp, "      ELSE\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "        B = B / A\n" );
	fprintf( fp, "        C = C / A\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "        RAD = 0.25D0 * B * B + C\n" );
	fprintf( fp, "        IF ( RAD.GE.0.0D0 ) THEN\n" );
	fprintf( fp, "          RAD = DSQRT( RAD ) + 0.5 * B\n" );
	fprintf( fp, "          IF ( RAD.NE.0.0D0 ) THEN\n" );
	fprintf( fp, "            SOLVQDRT = C / RAD\n" );
	fprintf( fp, "          ELSE\n" );
	fprintf( fp, "            SOLVQDRT = -0.5D0 * B + DSQRT( 0.25 * B * B + C )\n" );
	fprintf( fp, "          ENDIF\n" );
	fprintf( fp, "        ELSE\n" );
	fprintf( fp, "          SOLVQDRT = 0.0D0\n" );
	fprintf( fp, "        ENDIF\n" );
	fprintf( fp, "      END IF\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "\n\n" );

/* GETMOLARMASS */
	fprintf( fp, "      SUBROUTINE GETMOLARMASS( MM )\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "!	FILLS 'MM' WITH SPECIES MOLAR MASS IN KG/KMOLE\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      DOUBLE PRECISION MM(%d)\n"
				, gCounter->species );
	fprintf( fp, "      INCLUDE '%sF90.h'\n", gOptions->base );
	fprintf( fp, "\n" );
	ListIterator( gSpeciesList, PrintMagicMolarMassF77, fp );
	fprintf( fp, "\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "\n\n" );

/* GETSPECIESNAMES */
	fprintf( fp, "      SUBROUTINE GETSPECIESNAMES( NAMES )\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "!	FILLS 'NAMES' WITH SPECIES IDENTIFIER/KG\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      CHARACTER *20 NAMES(%d)\n"
				, gCounter->species );
	fprintf( fp, "      INCLUDE '%sF90.h'\n", gOptions->base );
	fprintf( fp, "\n" );
	ListIterator( gSpeciesList, PrintMagicNamesF77, fp );
	fprintf( fp, "\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "\n\n" );

/* GETNSPECIES */
	fprintf( fp, "      SUBROUTINE GETNSPECIES( NSPECIES )\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "!	FILLS 'NSPECIES' WITH NUMBER OF SPECIES \n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      INTEGER NSPECIES\n");
	fprintf( fp, "      INCLUDE '%sF90.h'\n", gOptions->base );
	fprintf( fp, "\n" );
	fprintf( fp, "      NSPECIES = SEND - 1\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "\n\n" );

/* GETNREACTIONS */
	fprintf( fp, "      SUBROUTINE GETNREACTIONS( NREACTIONS )\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "!	FILLS 'NREACTIONS' WITH NUMBER OF REACTIONS \n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      INTEGER NREACTIONS\n");
	fprintf( fp, "\n" );
	fprintf( fp, "      NREACTIONS = %d\n", gCounter->reactions );
	fprintf( fp, "\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "\n\n" );

/* GETNSPECS */
	fprintf( fp, "      SUBROUTINE GETNSPECS(NSPECIES_NONS)\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "!     RETURNS THE NUMBER OF THE NON STEADY STATE SPECIES\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "      implicit none\n" );
	fprintf( fp, "      integer ::  NSPECIES_NONS\n" );
	fprintf( fp, "      include '%sF90.h'\n", gOptions->base );
	fprintf( fp, "\n" );
	fprintf( fp, "      NSPECIES_NONS = %d\n", gCounter->species - gCounter->steadyStates );
	fprintf( fp, "      END\n" );
	fprintf( fp, "\n\n" );

/* GETMUCOEFF */
	fprintf( fp, "      SUBROUTINE GETMUCOEFF( MUCOEFF )\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "!	FILLS 'MUCOEFF' WITH MUECOEFF IN KG/(M*S)\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      implicit none\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      include '%sF90.h'\n", gOptions->base );
	fprintf( fp, "      real(DP) :: MUCOEFF(%d)\n", gCounter->species );
	fprintf( fp, "\n" );
	ListIterator( gSpeciesList, PrintMagicMuCoeffF77, fp );
	fprintf( fp, "\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "\n\n" ); 

/* GETKOVEREPS */
	fprintf( fp, "      SUBROUTINE GETKOVEREPS( KOVEREPS )\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "!	    FILLS 'KOVEREPS' WITH KOVEREPS\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      implicit none\n" );
	fprintf( fp, "\n" );
	fprintf( fp, "      include '%sF90.h'\n", gOptions->base );
	fprintf( fp, "      real(DP) :: KOVEREPS(%d)\n", gCounter->species );
	fprintf( fp, "\n" );
	ListIterator( gSpeciesList, PrintMagicKOverEpsF77, fp );
	fprintf( fp, "\n" );
	fprintf( fp, "      END\n" );
	fprintf( fp, "\n\n" );


/* COMPSTEADYSTATES */
	fprintf( fp, "      SUBROUTINE COMPSTEADYSTATES( K, C, M, TEMP, PRESSURE )\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "!     THIS SUBROUTINE COMPUTES THE STEADY STATE CONCENTRATIONS FROM\n" );
	fprintf( fp, "!     THE CONCENTRATIONS OF COMPUTED SPECIES AND RATE COEFFICIENTS.\n" );
	fprintf( fp, "!     CONCENTRATIONS OF COMPUTED SPECIES MAY NOT BE ALTERED.\n" );
	fprintf( fp, "!------------------------------------------------------------------\n" );
	fprintf( fp, "      IMPLICIT NONE\n" );
	if ( isReduced ) {
		fprintf( fp, "      DOUBLE PRECISION CTCHZERO,SOLVQDRT\n" );
	}
	fprintf( fp, "      INCLUDE '%sF90.h'\n", gOptions->base );
	fprintf( fp, "      integer  ::  NSPECIN, NSPEC, INOW\n" );
	fprintf( fp, "      real(DP) :: K(%d), C(%d), M(%d)\n",
		gCounter->reactions, gCounter->species, gCounter->thirdBodies );
	fprintf( fp, "      real(DP) :: TEMP, PRESSURE, CTOT\n");
	fprintf( fp, "      real(DP), parameter ::  R= 8314.34\n" );

	if ( isReduced ) {
		fprintf( fp, "      INCLUDE '%sF90.ss'\n", gOptions->base );
		fprintf( fp, "      NSPECIN = %d\n",gCounter->species - gCounter->steadyStates );
		fprintf( fp, "      NSPEC = %d\n",gCounter->species );
		fprintf( fp, "      CTOT = PRESSURE / ( R * TEMP )\n");
		fprintf( fp, "      DO INOW=NSPECIN+1,NSPEC\n");
		fprintf( fp, "        IF ( C(INOW).gt.CTOT ) THEN\n");
		fprintf( fp, "          C(INOW) = CTOT\n");
		fprintf( fp, "        END IF\n");
		fprintf( fp, "      END DO\n");
	}
	fprintf( fp, "      END\n" );
	fprintf( fp, "\n\n" );
	PrintMagicThermoDataF90( fp );
	
	fclose( fp );
}

void WriteRedMechReactions( void )
{
	char	fName[128];
	FILE	*fp;
	
	sprintf( fName, "%s.redmech", gOptions->base );
	if ( !( fp = fopen( fName, "w" ) ) ) {
		fprintf( stderr, "#error: unable to open file %s\n", fName );
		return;
	}
	
	ListIterator( gReactionList, PrintRedMechReaction, fp );
	fprintf( fp, "\n\n" );
	ListIterator( gSpeciesList, PrintRedMechProdRate, fp );
	
	fclose( fp );
}

/*void WriteMathematica( void )
{
	if ( gOptions->mathematica ) {
		ListIterator( gReactionList, PrintMathReactions, fp );	
	}	
}*/

static void Initialize( void )
{
	gAtomsList = NewList( "Atoms" );
	gSpeciesList = NewList( "Species" );
	gReactionList = NewList( "Reactions" );
	gPAHReactionList = NewList( "undefined" );
	gSootReactionList = NewList( "undefined" );
	gAllowAtomList = NewList( "Allowed Atoms" );
	gUsedThirdBodyList = NewList( "Used Third Bodies" );
	gThirdBodyList = NewList( "Third Body Coefficients" );
	gDimensionList = NewList( "Allowed Dimensions" );
	gAllowTempExpList = NewList( "Allowed Temperatur Parameter in exponents of units" );
	gAllowReacExpList = NewList( "Allowed order of reaction Parameter in exponents of units" );
	gUsedDimensionList = NewList( "Used Dimensions" );
	gCommentList = NewList( "Comments" );
	gLabelFollowingCommentList = NewList( "Label that follows the comment" );
	gFcNameList = NewList( "Names of Broadening" );
	gFcContentsList = NewList( "Contents of Broadening" );
	gFcLineNumberList = NewList( "Line Number" );
	gSteadyStateList = NewList( "Steady States" );

	gOptions = NewOptions();
	gLineContents[0] = '\0';
	gHeader = InitHeader( gHeader );
	gCounter = NewCounter();
	InitDimension( gDimensionList );
	InitAllowedExponents();
	gLabelIncrement	=	10;
}


static void Usage( int exit_code )
{
	FILE	*fp = (exit_code) ? stderr : stdout;
	
	fprintf( fp, "\nUsage: %s -i <mechanism> [-o outFile] [-t <thermo data>] [-R <steady states>] [abcdhprsuCLMS3VX] [-Z <num>]\n", gProgram );
	fprintf( fp, "         Options:\n" );
	fprintf( fp, "              -a:  print a list of allowed and used atoms\n" );
	fprintf( fp, "              -b:  print broadening factors\n" );
	fprintf( fp, "              -c:  do not check consistency of forward and backward reactions\n" );
	fprintf( fp, "                   and do not compute missing constants for Arrhenius functions\n" );
	fprintf( fp, "                   for pressure dependent reactions\n" );
	fprintf( fp, "                   (forward and backward reactions must have the same label,\n" );
	fprintf( fp, "                    except for the last char, which must be either 'f' or 'b')\n" );
	fprintf( fp, "              -h:  print some help information\n" );
	fprintf( fp, "              -W:  write FORTRAN 77 source term function\n" );
	fprintf( fp, "              -w:  write FORTRAN 90 source term function\n" );
	fprintf( fp, "              -p:  print progress information\n" );
	fprintf( fp, "              -r:  print all reactions\n" );
	fprintf( fp, "              -s:  print all species\n" );
	fprintf( fp, "              -u:  print a list of allowed a/o specified units\n" );
	fprintf( fp, "              -3:  print a list of used and defined third-body efficiencies\n" );
	fprintf( fp, "              -L:  generate a LaTeX file with reactions\n" );
	fprintf( fp, "              -M:  generate output for Mathematica\n" );
/* C originally for Josef */
	fprintf( fp, "              -C:  generate C source term function\n" );
	fprintf( fp, "              -S:  use string comparisons for matching species\n" );
	fprintf( fp, "              -V:  generate vectorizing FORTRAN code\n" );
	fprintf( fp, "              -X:  do not generate a binary output file (overrides -o option)\n" );
	fprintf( fp, "\n              -d [Scanner,Parser]: debug scanner a/o parser\n" );
	fflush( fp );

	CleanUp();
	exit( exit_code );
}


static void HandleCommandLine( int argc, char *argv[] )
{
	extern char 	*optarg;
	const char		*myOpts = "abchwWprsuCLMS3VXR:Z:i:o:t:d:";
	const char		*shit_man = "memory exhausted (HandleCommandLine)";
	int				anOpt;
	int				errors = 0;

	/*	Dispatch command line arguments.										*/
	while ( (anOpt = getopt( argc, argv, myOpts )) != EOF ) {
	
		switch ( anOpt ) {
		
			/*	Options with arguments										*/
			case 'i':
				gOptions->inName = optarg;
				break;
			case 'o':
				gOptions->outName = optarg;
				break;
			case 't':
				gOptions->thermoName = optarg;
				break;
			case 'd':
				if ( strstr( optarg, "Scanner" ) ) {
					gOptions->debugScanner = TRUE;
				}
				if ( strstr( optarg, "Parser" ) ) {
					yydebug = TRUE;
				}
				break;
			case 'R':
				gOptions->steadyStates = optarg;
			case 'Z':
				gOptions->Z = atoi(optarg);
				break;
			
			/*	Simple flags												*/
			case 'a':	gOptions->printAtoms = TRUE;			break;
			case 'b':	gOptions->printBroadening = TRUE;		break;
			case 'c':	gOptions->useForwardBackward = FALSE;	break;
			case 'h':	Usage( 0 );
			case 'C':	gOptions->outputWriter = TRUE;			break;
		    case 'w':	gOptions->outputWriterF90 = TRUE;		break;
			case 'W':	gOptions->outputWriterF77 = TRUE;		break;
			case 'p':	gOptions->progress = TRUE;				break;
			case 'r':	gOptions->printReactions = TRUE;		break;
			case 's':	gOptions->printSpecies = TRUE;			break;
			case 'u':	gOptions->printUnits = TRUE;			break;
			case 'L':	gOptions->latexTable = TRUE;			break;
			case 'M':	gOptions->mathematica = TRUE;			break;
			case 'S':	gSpeciesCompareFunc = FindSpeciesName;	break;
			case 'V':	gOptions->vectorizeF77 = TRUE;			break;
			case 'X':	gOptions->noBinaryOutput = TRUE;		break;
			case '3':	gOptions->printThird = TRUE;			break;
			
			/*	Handle errors												*/
			case '?':	++errors;							break;
			default:	FatalError( "something wrong in HandleCommandLine()" );
			
		}	/*	switch ( anOpt ) {	*/
		
	}	/* while ( (anOpt = getopt(...	*/
	
	if ( errors ) Usage( 2 );


	/*	Take care of additional initializations.							*/
	
	if ( gOptions->inName == NULL ) {
		fprintf( stderr, "\n### Required input file is missing.\n" );
		Usage( 2 );
	} else {
		int		base_length = strlen(gOptions->inName) - 5;	/*	without '\0'	*/
		char	*start = NULL;
		if ( strcmp( gOptions->inName + base_length, ".mech" ) != 0 ) {
			/*	Use full name, if the input file doesn't end in ".mech".	*/
			base_length += 5;
		}
		gOptions->base = (char *)malloc( base_length + 1 );
		if ( gOptions->base == NULL ) FatalError( shit_man );
		if ( start = strrchr( gOptions->inName, '/' ) ) {
			strncpy( gOptions->base, start+1, base_length );
			base_length -= start+1 - gOptions->inName;
		}
		else {
			strncpy( gOptions->base, gOptions->inName, base_length );
		}
		*(gOptions->base + base_length) = '\0';
	}

	/*	Try to open mechanism file for reading.								*/
	if ( ( gOptions->reactionInFile = fopen( gOptions->inName, "r" ) ) == NULL ) {
		fprintf( stderr, "### Unable to open mechanism file \"%s\".\n", gOptions->inName );
		CleanUp();
		exit( 2 );
	}
	yyin = gOptions->reactionInFile;
	
	/*	Try to open the thermo library for binary reading.					*/
	if ( ( gOptions->thermoInFile = fopen( gOptions->thermoName, "rb" ) ) == NULL ){
		if ( ( gOptions->thermoInFile = fopen( DefaultThermoFile( gOptions->thermoName ), "rb" ) ) == NULL ){
			fprintf( stderr, "### Unable to open thermo data file \"%s\".\n", gOptions->thermoName );
			CleanUp();
			exit( 2 );
		}
	}
	
	if ( gOptions->noBinaryOutput == FALSE ) {

		if ( !gOptions->outName ) {	/* We need to generate the file name.	*/
			const char *suffix = ".pre";
			gOptions->outName = (char *)malloc( strlen(gOptions->base)+strlen(suffix)+1 );
			if ( gOptions->outName == NULL ) FatalError( shit_man );
			sprintf( gOptions->outName, "%s%s", gOptions->base, suffix );
		}

		if ( ( gOptions->outFile = fopen( gOptions->outName, "wb" ) ) == NULL ){
			fprintf( stderr, "### Can't open output file \"%s\".\n", gOptions->outName );
			CleanUp();
			exit( 2 );
		}
#ifdef applec
		fsetfileinfo( gOptions->outName, 'MPS ', 'OBJ ' );
#endif
	}
}


static char *DefaultThermoFile( const char* name )
{
	char *path = NULL;
	char *fullName = NULL;
	
	path = getenv( "myData" );
	if ( path ) {
		if ( path[strlen( path )-1] == '/' ) {
			fullName = (char *)malloc( strlen(path) + strlen(name) + 1 );
			if ( fullName ) {
				strcpy( fullName, path );
				strcat( fullName, name );
			}
		}
		else {
			fullName = (char *)malloc( strlen(path) + strlen(name) + 2 );
			if ( fullName ) {
				strcpy( fullName, path );
				strcat( fullName, "/" );
				strcat( fullName, name );
			}
		}
	}
	else {
		fullName = (char *)malloc( strlen(name) + 1 );
		strcpy( fullName, name );
	}
#ifdef qDebug
	fprintf( stderr, "# Default thermo library: \"%s\"\n", fullName );
#endif
	
	return fullName;
}



OptionsPtr NewOptions( void )
{
	OptionsPtr	opts = NULL;

	if ( (opts = (OptionsPtr) malloc( sizeof( Options ) ) ) == NULL )
		FatalError( "memory exhausted (NewOptions)" );
	
	memset( opts, 0, sizeof(*opts) );

	opts->thermoName = DefaultThermoFile( "thermdat.bin" );	
	/*	Checking forward/backward stuff is now on default turned on.		*/
	opts->useForwardBackward = TRUE;
	
	return opts;
}


void FreeOptions( OptionsPtr opts )
{
	free( opts );
}


void PrintLists( void )
{
	if ( gOptions->printAtoms ) {
		ListIterator( gAllowAtomList, PrintString, gAllowAtomList );
		ListIterator( gAtomsList, PrintAtoms, gAtomsList );	
	}	
	if ( gOptions->printSpecies ) {
		ListIterator( gSpeciesList, PrintSpecies, gSpeciesList );	
	}	
	if ( gOptions->printThird ) {
		ListIterator( gUsedThirdBodyList, PrintString, gUsedThirdBodyList );	
		ListIterator( gThirdBodyList, PrintThirdBody, gThirdBodyList );	
	}	
	if ( gOptions->printUnits ) {
		ListIterator( gDimensionList, PrintDimensionTable, gDimensionList );
		ListIterator( gUsedDimensionList, PrintDimension, gUsedDimensionList );	
	}	
	if ( gOptions->printReactions ) {
		ListIterator( gReactionList, PrintReaction, gReactionList );	
	}	
	if ( gOptions->printReactions ) {
		ListIterator( gPAHReactionList, PrintReaction, gPAHReactionList );	
	}	
	if ( gOptions->printReactions ) {
		ListIterator( gSootReactionList, PrintReaction, gSootReactionList );	
	}	
	if ( gOptions->printBroadening ) {
		ListIterator( gFcNameList, PrintBroadening, gFcContentsList );	
	}	
	if ( gOptions->printSpecies ) {
		ListIterator( gSpeciesList, PrintSpeciesNames, gSpeciesList );	
	}	

}


void CleanUp( void )
{
	if ( gOptions->reactionInFile ) {
		fclose( gOptions->reactionInFile );
	}
	if ( gOptions->thermoInFile ) {
		fclose( gOptions->thermoInFile );
	}
	FreeList( gUsedDimensionList, FreeDimension );
	FreeList( gAllowReacExpList, FreeString );
	FreeList( gAllowTempExpList, FreeString );
	FreeList( gDimensionList, FreeDimension );
	FreeList( gThirdBodyList, FreeThirdBody );
	FreeList( gUsedThirdBodyList, FreeString );
	FreeList( gAllowAtomList, FreeString );
	FreeList( gReactionList, FreeReaction );
	FreeList( gSpeciesList, FreeSpecies );
	FreeList( gAtomsList, FreeAtoms );

	FreeOptions( gOptions );
}

void WriteScanManMech( void )
{
	LinkPtr			link = NULL;
	LinkPtr 		speciesLink = NULL;
	LinkPtr 		dimLinkA = NULL;
	LinkPtr 		dimLinkE = NULL;
	LinkPtr			thirdBodyLink = NULL;
	ReactionPtr 	reaction = NULL;
	DimensionPtr	dimA = NULL;
	DimensionPtr	dimE = NULL;
	SpeciesPtr		species = NULL;
	FILE			*fp = NULL;
	char			*fName = NULL;
	int				numberOfTables = 0;
	int				i, j, k;
	int				tableCounter;
	int				reactionCounter;
	int				reactionPerPage;
	Flag			first = TRUE;
	char			*thirdBody = NULL;
	int	base_length = strlen(gOptions->inName) - 5;	/*	without '\0'	*/

	if ( strcmp( gOptions->inName + base_length, ".mech" ) != 0 ) {
		/*	Use full name, if the input file doesn't end in ".mech".	*/
		base_length += 5;
	}
	fName = (char *)malloc( base_length + 5 ); /*	fName.mec	*/
	if ( fName == NULL ) FatalError( "shit_manmemory allocation failed" );
	fName[0] = '\0';
	strncpy( fName, gOptions->inName, base_length );
	fName[base_length] = '\0';
	strcat( fName, ".mec" );

	if ( !( fp = fopen( fName, "w" ) ) ) {
		fprintf( stderr, "warning: can't open file '%s'\n", fName );
		return;
	}

		for ( i = 0; i < gReactionList->items; ++i ) {
			link = Find_n_th_Item( i+1, gReactionList );
			reaction = ( ReactionPtr )link->item;
			
			if ( link = Find_n_th_Item( NumberOfListItem( gLabelFollowingCommentList, FindString, reaction->label ), gCommentList ) ) {
				thirdBody = ( char * )link->item;
				fprintf( fp, "\n#P%s\n", thirdBody );
			}
			fprintf( fp, "%s: ", reaction->label );
			
			first = TRUE;
			for ( k = 0; k < reaction->numberOfSpecies; ++k ) {
				speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[k] );
				if ( !speciesLink ) {
					fprintf( stderr, " #something wrong in WriteLatexTable\n" );
					goto Exit;
				}
				species = ( SpeciesPtr ) speciesLink->item;
		
				if ( reaction->speciesCoeff[k] > 0.0 ) {
					if ( !first ) {
						fprintf( fp, " + " );
					}
					if ( fabs( reaction->speciesCoeff[k] ) != 1.0 ) {
						fprintf( fp, "%g ", fabs( reaction->speciesCoeff[k] ) );
					}
					fprintf( fp, "%s", species->name );
					first = FALSE;
				}
			}
			if ( reaction->withThirdBody ) {
				if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
					thirdBody = ( char * ) thirdBodyLink->item;
					if ( atoi( &thirdBody[1] ) == 1 ) {
						fprintf( fp, " + M'" );
					}
					else if ( thirdBody[1] == '_' ) {
						fprintf( fp, " + %s", &thirdBody[2] );
					}
					else {
						fprintf( fp, " + M%d", atoi( &thirdBody[1] ) );
					}
				}
				else {
					fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
					goto Exit;
				}
			}
			
			if ( reaction->partialEquilibrium ) {
				fprintf( fp, "$==$" );
			}
			else {
				fprintf( fp, " -> " );
			}
			first = TRUE;
			for ( k = 0; k < reaction->numberOfSpecies; ++k ) {				speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[k] );
				if ( !speciesLink ) {
					fprintf( stderr, " #something wrong in WriteLatexTable\n" );
					goto Exit;
				}
				species = ( SpeciesPtr ) speciesLink->item;
		
				if ( reaction->speciesCoeff[k] < 0.0 ) {
					if ( !first ) {
						fprintf( fp, " + " );
					}
					if ( fabs( reaction->speciesCoeff[k] ) != 1.0 ) {
						fprintf( fp, "%g ", fabs( reaction->speciesCoeff[k] ) );
					}
					fprintf( fp, "%s", species->name );
					first = FALSE;
				}
			}
			if ( reaction->withThirdBody ) {
				if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
					thirdBody = ( char * ) thirdBodyLink->item;
					if ( atoi( &thirdBody[1] ) == 1 ) {
						fprintf( fp, " + M'" );
					}
					else if ( thirdBody[1] == '_' ) {
						fprintf( fp, " + %s", &thirdBody[2] );
					}
					else {
						fprintf( fp, " + M%d", atoi( &thirdBody[1] ) );
					}
				}
				else {
					fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
					goto Exit;
				}
			}
			fprintf( fp, "\n\t\t\t\t\t\t\t\t{" );

			dimLinkA = ListIterator( gUsedDimensionList, FindDimension, "A" );
			if ( dimLinkA ) {
				dimA = ( Dimension * ) dimLinkA->item;
			}
			else {
				fprintf( stderr, "#something's wrong with the unit of A\n" );
				goto Exit;
			}
			fprintf( fp, "\tA = %9.4E"
				, reaction->a * dimA->value * 1.0e-6 
				* pow( dimA->orderOfReacValue * 1.0e6, (Double)GetOrderOfReaction( reaction, FALSE ) ) );
			fprintf( fp, "\tn = %4.3f", reaction->n );
			dimLinkE = ListIterator( gUsedDimensionList, FindDimension, "E" );
			if ( dimLinkE ) {
				dimE = ( Dimension * ) dimLinkE->item;
			}
			else {
				fprintf( stderr, "#something's wrong with the unit of E\n" );
				goto Exit;
			}
			fprintf( fp, "\tE = %8.5G", reaction->e * dimE->value / 1000.0 );
			if ( reaction->withLindemann ) {
				fprintf( fp, "\tAi = %9.4E", reaction->aInf * dimA->value * 1.0e-6
				* pow( dimA->orderOfReacValue * 1.0e6, (Double)GetOrderOfReaction( reaction, TRUE ) ) );
			 	fprintf( fp, "\tni = %4.3f", reaction->nInf );
				fprintf( fp, "\tEi = %8.5G", reaction->eInf * dimE->value / 1000.0 );
				if ( reaction->fca ) {
					fprintf( fp, "\n\t\t\t\t\tfcA = %7.4f\tfcTA = %7.4f"
							, reaction->fca, reaction->fcTa );
				}
				if ( reaction->fcb ) {
					fprintf( fp, "\n\t\t\t\t\tfcB = %7.4f\tfcTB = %7.4f"
							, reaction->fcb, reaction->fcTb );
				}
				if ( reaction->fcc ) {
					fprintf( fp, "\n\t\t\t\t\tfcC = %7.4f\tfcTC = %7.4f"
							, reaction->fcc, reaction->fcTc );
				}
			}
			fprintf( fp, "\t}\n" );
		}

Exit:
	fclose( fp );
}

#ifdef OLDWRONGLATEXTABLE
void WriteLatexTable( void )
{
	const int		maxReacPerPage = 33; /* changed to 38 from 29 and back to 33 */
	LinkPtr			link = NULL;
	LinkPtr 		speciesLink = NULL;
	LinkPtr			thirdBodyLink = NULL;
	ReactionPtr 	reaction = NULL;
	SpeciesPtr		species = NULL;
	FILE			*fp = NULL;
	char			*fName = NULL;
	int				numberOfTables = 0;
	int				i, j, k;
	int				tableCounter;
	int				reactionCounter;
	int				reactionPerPage;
	Flag			first = TRUE;
	char			*thirdBody = NULL;
	int	base_length = strlen(gOptions->inName) - 5;	/*	without '\0'	*/

	if ( strcmp( gOptions->inName + base_length, ".mech" ) != 0 ) {
		/*	Use full name, if the input file doesn't end in ".mech".	*/
		base_length += 5;
	}
	fName = (char *)malloc( base_length + 5 ); /*	fName.tex	*/
	if ( fName == NULL ) FatalError( "shit_manmemory allocation failed" );
	fName[0] = '\0';
	strncpy( fName, gOptions->inName, base_length );
	fName[base_length] = '\0';
	strcat( fName, ".tex" );

	if ( !( fp = fopen( fName, "w" ) ) ) {
		fprintf( stderr, "warning: can't open file '%s'\n", fName );
		return;
	}

	fprintf( fp, "\\documentstyle[10pt, fleqn]{article}\n" );
	fprintf( fp, "\\sloppy\n" );
	fprintf( fp, "\\frenchspacing\n" );
	fprintf( fp, "\\renewcommand{\\baselinestretch}{1.2}\n" );
	fprintf( fp, "\\addtolength{\\topmargin}{-50pt}\n" );
	fprintf( fp, "\\addtolength{\\textheight}{140pt}\n" );
	fprintf( fp, "\\pagestyle{empty}\n" );
	fprintf( fp, "\\hfuzz=140truept\n\n" );
	fprintf( fp, "\\begin{document}\n" );
	fprintf( fp, "\\oddsidemargin 1cm\n" );

#ifdef SKIPOMPARSE

	numberOfTables = gReactionList->items / maxReacPerPage;
	if ( gReactionList->items % maxReacPerPage ) {
	   ++numberOfTables;
	}
	for ( tableCounter = 0; tableCounter < numberOfTables; ++ tableCounter ) {
		reactionPerPage = MIN( maxReacPerPage, gReactionList->items - tableCounter * maxReacPerPage );
		
		fprintf( fp, "\\begin{center}\n" );
		fprintf( fp, "\\begin{tabular}{|c|rcll|c|c|c|}\n" );
		fprintf( fp, "\\hline\n" );
		fprintf( fp, "Number&\\multicolumn{4}{|c|}{Reaction}&$A$&$n$&$E$\\\\\n" );
		fprintf( fp, "\\hline\n" );
		fprintf( fp, "\\hline\n" );

		for ( i = 0; i < reactionPerPage; ++i ) {
			reactionCounter = maxReacPerPage * tableCounter + i + 1;
			
			link = Find_n_th_Item( reactionCounter, gReactionList );
			reaction = ( ReactionPtr )link->item;
			
			if ( link = Find_n_th_Item( NumberOfListItem( gLabelFollowingCommentList, FindString, reaction->label ), gCommentList ) ) {
				thirdBody = ( char * )link->item;
				fprintf( fp, "\\multicolumn{8}{|c|}{%s}\\\\\n", thirdBody );
				fprintf( fp, "\\hline\n" );
			}
			fprintf( fp, "%s", reaction->label );
		  fprintf( fp, "&" );
			
			first = TRUE;
			for ( k = 0; k < reaction->numberOfSpecies; ++k ) {
				speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[k] );
				if ( !speciesLink ) {
					fprintf( stderr, " #something wrong in WriteLatexTable\n" );
					goto Exit;
				}
				species = ( SpeciesPtr ) speciesLink->item;
		
				if ( reaction->speciesCoeff[k] > 0.0 ) {
					if ( !first ) {
						fprintf( fp, " $+$ " );
					}
					if ( fabs( reaction->speciesCoeff[k] ) != 1.0 ) {
						fprintf( fp, "%g", fabs( reaction->speciesCoeff[k] ) );
					}
					for ( j = 0; j < strlen( species->name ); ++j ) {
						if ( isdigit( species->name[j] ) ) {
							fprintf( fp, "$_{%d}$", atoi( &species->name[j] ) );
							while ( isdigit( species->name[++j] ) );
							--j;
						}
						else {
							fprintf( fp, "%c", species->name[j] );
						}
					}
					first = FALSE;
				}
			}
			if ( reaction->withThirdBody ) {
			
				if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
					thirdBody = ( char * ) thirdBodyLink->item;
					fprintf( fp, " $+$ M$" );
					for ( j = 0; j < atoi( &thirdBody[1] ); ++j ) {
						fprintf( fp, "'" );
					}
					fprintf( fp, "$" );
				}
				else {
					fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
					goto Exit;
				}
			}
			
		  fprintf( fp, "&" );
			if ( reaction->partialEquilibrium ) {
				fprintf( fp, "$==$" );
			}
			else {
				fprintf( fp, "$\\rightarrow$" );
			}
		  fprintf( fp, "&" );
			first = TRUE;
		  	if ( !reaction->withLindemann ) {
				fprintf( fp, "\\multicolumn{2}{l}{" );
			}
			for ( k = 0; k < reaction->numberOfSpecies; ++k ) {
				speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[k] );
				if ( !speciesLink ) {
					fprintf( stderr, "#something wrong in WriteLatexTable\n" );
					goto Exit;
				}
				species = ( SpeciesPtr ) speciesLink->item;
		
				if ( reaction->speciesCoeff[k] < 0.0 ) {
					if ( !first ) {
						fprintf( fp, " $+$ " );
					}
					if ( fabs( reaction->speciesCoeff[k] ) != 1.0 ) {
						fprintf( fp, "%g", fabs( reaction->speciesCoeff[k] ) );
					}
					for ( j = 0; j < strlen( species->name ); ++j ) {
						if ( isdigit( species->name[j] ) ) {
							fprintf( fp, "$_{%d}$", atoi( &species->name[j] ) );
							while ( isdigit( species->name[++j] ) );
							--j;
						}
						else {
							fprintf( fp, "%c", species->name[j] );
						}
					}
					first = FALSE;
				}
			}
			if ( reaction->withThirdBody ) {
				if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
					thirdBody = ( char * ) thirdBodyLink->item;
					fprintf( fp, " $+$ M$" );
					for ( j = 0; j < atoi( &thirdBody[1] ); ++j ) {
						fprintf( fp, "'" );
					}
					fprintf( fp, "$" );
				}
				else {
					fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
					goto Exit;
				}
			}
		  	if ( !reaction->withLindemann ) {
				fprintf( fp, "}\\vline" );
			}
			else {
		  		fprintf( fp, "&$k_{0}$" );
			}
		  fprintf( fp, "&" );
		  	fprintf( fp, "%9.3E", reaction->a );
		  fprintf( fp, "&" );
		  	fprintf( fp, "%4.2f", reaction->n );
		  fprintf( fp, "&" );
		  	fprintf( fp, "%4.3G", reaction->e );
		  fprintf( fp, "\\\\\n" );
		  if ( reaction->withLindemann ) {
				/*fprintf( fp, "\\multicolumn{4}{c}{}&$k_{\\infty}$&" );*/
				fprintf( fp, "&&&&$k_{\\infty}$&" );
				fprintf( fp, "%9.3E", reaction->aInf );
			  fprintf( fp, "&" );
				fprintf( fp, "%4.2f", reaction->nInf );
			  fprintf( fp, "&" );
				fprintf( fp, "%4.3G", reaction->eInf );
			  fprintf( fp, "\\\\\n" );
		  }
		  fprintf( fp, "\\hline\n" );
		}

		  fprintf( fp, "\\end{tabular}\n" );
		  fprintf( fp, "\\end{center}\n" );
		  if ( tableCounter != numberOfTables - 1 ) {
			fprintf( fp, "\\newpage\n" );
		}
    }
	fprintf( fp, "\\end{document}\n" );
#endif /* SKIPOMPARSE */

Exit:
	fclose( fp );
#ifdef applec
	fsetfileinfo( fName, '*TEX', 'TEXT' );
#endif
}
#else
void WriteLatexTable( void )
{
/* Function written by Adelbert Grudno */

	const int		maxLinePerPage = 38;
	LinkPtr			link = NULL;
	LinkPtr 		speciesLink = NULL;
	LinkPtr			thirdBodyLink = NULL;
	ReactionPtr 	reaction = NULL;
	SpeciesPtr		species = NULL;
	FILE			*fp = NULL;
	char			*fName = NULL;
	int				i, j, k;
	int				reactionCounter = 1;
	Flag			first;
	char			*thirdBody = NULL;
	char			*commentLine = NULL;
	int	base_length = strlen(gOptions->inName) - 5;	/*	without '\0'	*/
	char			buffer[32];
	int				rcount = 1;

	
	if ( strcmp( gOptions->inName + base_length, ".mech" ) != 0 ) {
		/*	Use full name, if the input file doesn't end in ".mech".	*/
		base_length += 5;
	}
	fName = (char *)malloc( base_length + 5 ); /*	fName.tex	*/
	if ( fName == NULL ) FatalError( "shit_manmemory allocation failed" );
	fName[0] = '\0';
	strncpy( fName, gOptions->inName, base_length );
	fName[base_length] = '\0';
	strcat( fName, ".tex" );

	if ( !( fp = fopen( fName, "w" ) ) ) {
		fprintf( stderr, "warning: can't open file '%s'\n", fName );
		return;
	}
 
	fprintf( fp, "\\documentclass[12pt]{article}\n" );
	fprintf( fp, "\\sloppy\n" );
    fprintf( fp, "\\topmargin = 0mm\n" );
    fprintf( fp, "\\headheight = 0mm\n" );
    fprintf( fp, "\\headsep = 0mm\n" );
    fprintf( fp, "\\topskip = 0mm\n" );
    fprintf( fp, "\\oddsidemargin = 0mm\n" );
    fprintf( fp, "\\textwidth = 6.5in\n" );
    fprintf( fp, "\\textheight = 8.5in\n" );
    fprintf( fp, "%\\pagestyle{empty}\n" );
	fprintf( fp, "%\\oddsidemargin 1cm\n" );
    fprintf( fp, "\\begin{document}\n" );

	while (reactionCounter <= gReactionList->items) {
		fprintf( fp, "\\begin{center}\n" );
#ifdef REFERENCE
		fprintf( fp, "\\begin{tabular}{|c|rcll|c|c|c|c|}\n" );
		fprintf( fp, "\\hline\n" );
		fprintf( fp, "Number&\\multicolumn{4}{|c|}{Reaction}&$A$&$n$&$E$&Reference\\\\\n" );documentclass[12pt]{article} 
#else
		fprintf( fp, "\\begin{tabular}{|c|rcll|c|c|c|}\n" );
		fprintf( fp, "\\hline\n" );
		fprintf( fp, "Number&\\multicolumn{4}{|c|}{Reaction}&$A$&$n$&$E$\\\\\n" );
#endif		
		fprintf( fp, "\\hline\n" );
		fprintf( fp, "\\hline\n" );

		for ( i = 0; (i < maxLinePerPage) && (reactionCounter <= gReactionList->items); ++i, ++reactionCounter ) {
			link = Find_n_th_Item( reactionCounter, gReactionList );
			reaction = ( ReactionPtr )link->item;

			if ( reaction->hasForward ) {
			  /*				continue;*/
			}

			
#ifndef SUPPRESSCOMMENTS
			if ( link = Find_n_th_Item( NumberOfListItem( gLabelFollowingCommentList, FindString, reaction->label ), gCommentList ) ) {
				commentLine = ( char * )link->item;
#ifdef REFERENCE
				fprintf( fp, "\\multicolumn{9}{|c|}{%s}\\\\\n", commentLine );
#else
				fprintf( fp, "\\multicolumn{8}{|c|}{%s}\\\\\n", commentLine );
#endif		
				fprintf( fp, "\\hline\n" );
				++i;
			}
#endif
#ifdef ENUMLATEXTABLE
			fprintf( fp, "%d", rcount );
			if ( reaction->hasBackward ) {
				fprintf( fp, "f" );
				/*				rcount++;*/
			}
			else if ( reaction->hasForward ) {
				fprintf( fp, "b" );
								rcount++;
			}
			else {
				rcount++;
			}
			
#else
			fprintf( fp, "%s", reaction->label );
#endif
		  fprintf( fp, "&" );
			
			/** begin left hand side of reaction **/
			for ( k = 0, first = TRUE; k < reaction->numberOfSpecies; ++k) {
				speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[k] );
				if ( !speciesLink ) {
					fprintf( stderr, " #something wrong in WriteLatexTable\n" );
					goto Exit;
				}
				species = ( SpeciesPtr ) speciesLink->item;
		
				if ( reaction->speciesCoeff[k] > 0.0 ) {
					if ( !first ) {
						fprintf( fp, " $+$ " );
					}
					if ( fabs( reaction->speciesCoeff[k] ) != 1.0 ) {
						fprintf( fp, "%g\\,", fabs( reaction->speciesCoeff[k] ) );
					}
					if ( strlen( GetSpeciesName( buffer, species->name ) ) != strlen(species->name) ) {
					  fprintf( fp, "%s", GetSpeciesName( buffer, species->name ) );
					}
					else {
					for ( j = 0; j < strlen( GetSpeciesName( buffer, species->name ) ); ++j ) {
						if ( isdigit( GetSpeciesName( buffer, species->name )[j] ) ) {
							if ( strchr( GetSpeciesName( buffer, species->name ), '-' ) == 0 
								|| j > strchr( GetSpeciesName( buffer, species->name ), '-' ) - GetSpeciesName( buffer, species->name ) ) {
								fprintf( fp, "$_{%d}$", atoi( &GetSpeciesName( buffer, species->name )[j] ) );
								while ( isdigit( GetSpeciesName( buffer, species->name )[++j] ) );
								--j;
							} 
							else {
								fprintf( fp, "%c", GetSpeciesName( buffer, species->name )[j] );
							}
						}
						else {
							fprintf( fp, "%c", GetSpeciesName( buffer, species->name )[j] );
						}
					}
					}
					first = FALSE;
				}
			}
			if ( reaction->withThirdBody ) {
			
				if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
					thirdBody = ( char * ) thirdBodyLink->item;
					fprintf( fp, " $+$ M" );
					if ( atoi( &thirdBody[1] ) < 4 ) {
					  fprintf( fp, "$" );
					for ( j = 0; j < atoi( &thirdBody[1] ); ++j ) {
						fprintf( fp, "'" );
					}
					fprintf( fp, "$" );
				}
				else {
						fprintf( fp, "%d", atoi( &thirdBody[1] ) );
					}
				}
				else {
					fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
					goto Exit;
				}
			}
			/** end left hand side of reaction **/
			
		  fprintf( fp, "&" );
			if ( reaction->partialEquilibrium ) {
				fprintf( fp, "$==$" );
			}
			else {
				fprintf( fp, "$\\rightarrow$" );
			}
		  fprintf( fp, "&" );
			if ( !reaction->withLindemann ) {
				fprintf( fp, "\\multicolumn{2}{l}{" );
			}

			/** begin right hand side of reaction **/
			for ( k = 0, first = TRUE; k < reaction->numberOfSpecies; ++k ) {
				speciesLink = ListIterator( gSpeciesList, FindSpeciesNumber, &reaction->speciesNumber[k] );
				if ( !speciesLink ) {
					fprintf( stderr, "#something wrong in WriteLatexTable\n" );
					goto Exit;
				}
				species = ( SpeciesPtr ) speciesLink->item;
		
				if ( reaction->speciesCoeff[k] < 0.0 ) {
					if ( !first ) {
						fprintf( fp, " $+$ " );
					}
					if ( fabs( reaction->speciesCoeff[k] ) != 1.0 ) {
						fprintf( fp, "%g\\,", fabs( reaction->speciesCoeff[k] ) );
					}
                    if ( strlen( GetSpeciesName( buffer, species->name ) ) != strlen(species->name) ) {
                      fprintf( fp, "%s", GetSpeciesName( buffer, species->name ) );
                    }
                    else {
					for ( j = 0; j < strlen( GetSpeciesName( buffer, species->name ) ); ++j ) {
						if ( isdigit( GetSpeciesName( buffer, species->name )[j] ) ) {
							if ( strchr( GetSpeciesName( buffer, species->name ), '-' ) == 0 
								|| j > strchr( GetSpeciesName( buffer, species->name ), '-' ) - GetSpeciesName( buffer, species->name ) ) {
								fprintf( fp, "$_{%d}$", atoi( &GetSpeciesName( buffer, species->name )[j] ) );
								while ( isdigit( GetSpeciesName( buffer, species->name )[++j] ) );
								--j;
							} 
							else {
								fprintf( fp, "%c", GetSpeciesName( buffer, species->name )[j] );
							}
						}
						else {
							fprintf( fp, "%c", GetSpeciesName( buffer, species->name )[j] );
						}
					}
					}
					first = FALSE;
				}
			}
			if ( reaction->withThirdBody ) {
				if ( thirdBodyLink = Find_n_th_Item( reaction->thirdBodyNumber+1, gUsedThirdBodyList ) ) {
					thirdBody = ( char * ) thirdBodyLink->item;
					fprintf( fp, " $+$ M" );
					if ( atoi( &thirdBody[1] ) < 4 ) {
					  fprintf( fp, "$" );
					for ( j = 0; j < atoi( &thirdBody[1] ); ++j ) {
						fprintf( fp, "'" );
					}
					fprintf( fp, "$" );
				}
				else {
						fprintf( fp, "%d", atoi( &thirdBody[1] ) );
					}
				}
				else {
					fprintf( stderr, "#something's wrong with the third body no.%d\n", reaction->thirdBodyNumber );
					goto Exit;
				}
			}
			/** end right hand side of reaction **/

			/** begin rate constants **/
			if ( !reaction->withLindemann ) {
				fprintf( fp, "}\\vline" );
			}
			else {
				fprintf( fp, "&$k_{0}$" );
			}
		  fprintf( fp, "&" );
			fprintf( fp, "%9.3E", reaction->a );
		  fprintf( fp, "&" );
			fprintf( fp, "%4.2f", reaction->n );
		  fprintf( fp, "&" );
			fprintf( fp, "%4.3G", reaction->e );
#ifdef REFERENCE
		  fprintf( fp, "&\\cite{baulch}" );
#endif		
		  fprintf( fp, "\\\\\n" );
		  if ( reaction->withLindemann ) {
				/*fprintf( fp, "\\multicolumn{4}{c}{}&$k_{\\infty}$&" );*/
				fprintf( fp, "&&&&$k_{\\infty}$&" );
				fprintf( fp, "%9.3E", reaction->aInf );
			  fprintf( fp, "&" );
				fprintf( fp, "%4.2f", reaction->nInf );
			  fprintf( fp, "&" );
				fprintf( fp, "%4.3G", reaction->eInf );
#ifdef REFERENCE
			  fprintf( fp, "&\\cite{baulch}" );
#endif		
			  fprintf( fp, "\\\\\n" );
				++i;
		  }
		  fprintf( fp, "\\hline\n" );
			/** end rate constants **/
		}

		fprintf( fp, "\\end{tabular}\n" );
		fprintf( fp, "\\end{center}\n" );
		fprintf( fp, "\\newpage\n" );
	}

	fprintf( fp, "\\end{document}\n" );

 Exit:
	fclose( fp );
}
#endif

char *GetSpeciesName( char *buffer, char *name )
{
	char	buffer2[32];
	char	*isomere;

	strcpy( buffer2, name );
	if ( ( isomere = strrchr( buffer2, '-' ) ) ) {
		if (isomere - buffer2 > 3) *isomere = '\0';
/*		if (isomere - buffer2 > 3) *isomere = '\0';*/
	}
	strcpy( buffer, buffer2 );

	return buffer;
}

void SetCounter( void )
{
	int			i, nSteadyStates = 0, nSpecies = gSpeciesList->items;
	LinkPtr 	link = gSpeciesList->firstItem;
	SpeciesPtr	species;

	gCounter = NewCounter();
	
	gCounter->atoms = gAtomsList->items;
	gCounter->species = nSpecies;
	
	gCounter->pahSpecies = 0;
	gCounter->sootSpecies = 0;

	for ( i = 0; i < nSpecies; ++i ) {
		species = ( SpeciesPtr ) link->item;
		if ( species->isSteadyState ) ++nSteadyStates;

		if ( !FindListItem( gReactionList, species, FindSpeciesInReaction ) ) {
			if ( FindListItem( gPAHReactionList, species, FindSpeciesInReaction ) ) {
				/*fprintf( stderr, "species '%s' is a PAH species\n", species->name );*/
				++gCounter->pahSpecies;
			}
			else {
				if ( FindListItem( gSootReactionList, species, FindSpeciesInReaction ) ) {
					/*fprintf( stderr, "species '%s' is a soot species\n", species->name );*/
					++gCounter->sootSpecies;
				}
				/*else {
					fprintf( stderr, "species '%s' is inert\n", species->name );
				}*/
			}
		}
		link = link->next;
	}	
	gCounter->steadyStates = nSteadyStates;
	
	gCounter->reactions = gReactionList->items;
	gCounter->pahReactions = gPAHReactionList->items;
	gCounter->sootReactions = gSootReactionList->items;
	gCounter->thirdBodies = gThirdBodyList->items;
	gCounter->dimensions = gUsedDimensionList->items;
}


void WriteOutput( void )
{
	FILE	*fp = gOptions->outFile;
/*	int 	i;*/
	
	/*	Sollte die Ausgabedatei nicht erst hier gešffnet werden???	*/
#ifdef SKIPOMPARSE
	
	if ( !fwrite( gHeader, sizeof( Header ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping header\n" );
		exit( 2 );
	}

	SetCounter();
	if ( !fwrite( gCounter, sizeof( Counter ), 1, fp ) ) {
		fprintf( stderr, "# error while dumping counter\n" );
		exit( 2 );
	}

	if ( ListIterator( gAtomsList, DumpAtoms, gAtomsList ) ) {
		fprintf( stderr, "# error while dumping list of used atoms\n" );
		exit( 2 );
	}
	if ( ListIterator( gSpeciesList, DumpSpecies, gSpeciesList ) ) {
		fprintf( stderr, "# error while dumping list of used species\n" );
		exit( 2 );
	}
	if ( ListIterator( gReactionList, DumpReaction, gReactionList ) ) {
		fprintf( stderr, "# error while dumping list of reactions\n" );
		exit( 2 );
	}
	strcpy( gbuffer, gPAHReactionList->listName );
	if ( !fwrite( gbuffer, gHeader->maxLenOfString, 1, gOptions->outFile ) ) {
		fprintf( stderr, "# error while dumping label of reaction\n");
		exit( 2 );
	}
	if ( ListIterator( gPAHReactionList, DumpReaction, gPAHReactionList ) ) {
		fprintf( stderr, "# error while dumping list of pahreactions\n" );
		exit( 2 );
	}
	strcpy( gbuffer, gSootReactionList->listName );
	if ( !fwrite( gbuffer, gHeader->maxLenOfString, 1, gOptions->outFile ) ) {
		fprintf( stderr, "# error while dumping label of reaction\n");
		exit( 2 );
	}
	if ( ListIterator( gSootReactionList, DumpReaction, gPAHReactionList ) ) {
		fprintf( stderr, "# error while dumping list of sootreactions\n" );
		exit( 2 );
	}
	if ( ListIterator( gThirdBodyList, DumpThirdBody, gThirdBodyList ) ) {
		fprintf( stderr, "# error while dumping list of used thirdbody efficiencies\n");
		exit( 2 );
	}
	if ( ListIterator( gUsedDimensionList, DumpDimension, gUsedDimensionList ) ) {
		fprintf( stderr, "# error while dumping list of used units\n" );
		exit( 2 );
	}
#endif /* SKIPOMPARSE */
	fclose( fp );
}



CounterPtr NewCounter( void )
{
	CounterPtr	counter = NULL;
	
	if ( !( counter = ( CounterPtr )malloc( sizeof( Counter ) ) ) ) {
		fprintf( stderr, "# memory allocation of Counter failed\n");
		exit(2);
	}
	
	return counter;
}

HeaderPtr InitHeader( HeaderPtr header )
{
	if ( !( header = ( HeaderPtr )malloc( sizeof( Header ) ) ) ) {
		fprintf( stderr, "# memory allocation of Header failed\n");
		exit(2);
	}
	strcpy( header->version, VERSION );
	header->maxLenOfString = 128;
	gbuffer = ( char * ) malloc( header->maxLenOfString );
	
	return header;
}

long GetFileSize( FILE *fp )
{
	long	nBytes = -1;
	
#ifdef HP
	fpos_t	position;

	if ( fgetpos( fp, &position ) != 0 ) {
		return 0;
	}
#elif defined (applec) || defined (powerc)
	fpos_t	position;

	if ( fgetpos( fp, &position ) != 0 ) {
		return 0;
	}
#elif defined SUNGXX
	fpos_t	position;

	if ( fgetpos( fp, &position ) != 0 ) {
		return 0;
	}
#else
	fpos_t	position;

	if ( fgetpos( fp, &position ) != 0 ) {
		return 0;
	}
#endif
	if ( fseek( fp, 0L, SEEK_END ) != 0 ) {
		return 0;
	}

	if ( ( nBytes = ftell( fp ) ) == -1L ) {
		return 0;
	}
	
#ifdef HP
	if ( fsetpos( fp, &position ) != 0 ) {
		return 0;
	}
#elif defined (applec) || defined (powerc)
	if ( fsetpos( fp, &position ) != 0 ) {
		return 0;
	}
#elif defined SUNGXX
	if ( fsetpos( fp, &position ) != 0 ) {
		return 0;
	}
/*#elif defined SUNOLD*/
/*	if ( rewind( fp ) == -1L ) {*/
/*		return 0;*/
/*	}*/
#else
	if ( fsetpos( fp, &position ) != 0 ) {
		return 0;
	}
#endif

	return nBytes;
}

void ReadThermoProperties( void )
{
	Flag				error = FALSE;
	int					i;
	int					nSpecies = gSpeciesList->items;
	SpeciesRecordHigh	dummy;
	FILE				*fp = gOptions->thermoInFile;
	LinkPtr 			link = NULL;
	SpeciesPtr			species = NULL;
	int					nRecords;
	Flag				warningHeaderOut = FALSE;

	if ( ( nRecords = GetFileSize( fp ) / sizeof( SpeciesRecord ) ) == 0 ) {
		fprintf( stderr, "# error while reading size of thermoprops inputfile" );
		CleanUp();
		exit( 2 );
	}

	dummy.speciesRecordPtr = NewSpeciesRecordArray( nRecords );
	dummy.len = nRecords;
	
	if ( fread( dummy.speciesRecordPtr, sizeof(SpeciesRecord), nRecords, fp ) != nRecords ) {
		fprintf( stderr, "# error while reading SpeciesRecords\n" );
		CleanUp();
		exit( 2 );
	}

	/* convert species names to uppercase */
        for ( i = 0; i < dummy.len; ++i ) {
          UpperString( dummy.speciesRecordPtr[i].name );
        }
	
#ifdef SHOWTHERMOFILE
	for ( i = 0; i < dummy.len; ++i ) {
		fprintf( stderr, "%s\tM = %g\n"
					, dummy.speciesRecordPtr[i].name
					, dummy.speciesRecordPtr[i].M );
	}
#endif

	ListIterator( gSpeciesList, GetThermoProps, &dummy );

	for ( i = 0; i < nSpecies; ++i ) {
		link = Find_n_th_Item( i+1, gSpeciesList );
		species = ( SpeciesPtr ) link->item;
		if ( !species->molarMass ) {
			/* try to match isomere part of species name */
			if ( !MatchIsomereName( link, &dummy ) ) {
				SetZeroProp(species);			
				fprintf( stderr, "### Warning: set NASA polynomial coefficients of '%s' to zero\n", species->name );
				
			}
		}
	}

	if ( warningHeaderOut ) {
		FreeSpeciesRecordArray( dummy.speciesRecordPtr );
		exit(2);
		/*return;*/
	}

	FreeSpeciesRecordArray( dummy.speciesRecordPtr );

}


static void ConflictWarning( const char *conflict, const char *chosen, const char *reason )
{
	fprintf( stderr, "### Warning: Conflict matching species \"%s\":\n", conflict );
	fprintf( stderr, "#            resolved by choosing \"%s\", because %s.\n", chosen, reason );
}

int GetThermoProps( LinkPtr link, void *var )
{
	SpeciesPtr 			species = ( SpeciesPtr ) link->item;
	SpeciesRecordHigh	recordStruct = *( ( SpeciesRecordHigh * ) var );
	SpeciesRecordPtr	record = recordStruct.speciesRecordPtr;
	int					len = recordStruct.len;
	int					i;
	char				buffer[32];
	Flag				foundnotexact = FALSE;
	static IntVectorPtr	*composition = NULL;

	if ( link == gSpeciesList->firstItem ) {
		composition = NewIntVectorArray( len, gAllowAtomList->items );
		for ( i = 0; i < len; ++i ) {
			GetComposition( record[i].name, composition[i] );
		}
	}
	 
	buffer[0] = '\0';

	for ( i = 0; i < len; ++i ) {
		gComposition = composition[i];
		if ( FindSpecies( link, record[i].name ) ) {
			if ( buffer[0] == '\0' ) {
				/* matched and first */
				CopySpecies( species, &record[i] );
				if ( strcmp( species->name, record[i].name ) == 0 ) {
					/* exact match */
#				ifdef DEBUGCHOICE
					fprintf( stderr, "take %s for %s\n", record[i].name, species->name );
#				endif
					return 0;
				}
				else {
					foundnotexact = TRUE;
					strcpy( buffer, record[i].name );
				}

#				ifdef DEBUG
				PrintSpecies( link, gSpeciesList );
#				endif
			} 
			else {
				/* matched and not first */
				if ( strcmp( species->name, record[i].name ) == 0 ) {
					/* the current matches exactly */
#				ifdef DEBUGCHOICE
					fprintf( stderr, "take %s for %s\n", record[i].name, species->name );
#				endif
					CopySpecies( species, &record[i] );
					return 0;
				}
			}
		}
	}

	if ( foundnotexact ) {
		fprintf( stderr, "### Warning: No exact match for species '%s' in thermo data file\n", species->name );
		fprintf( stderr, "#            resolved by choosing first similar species '%s'\n", buffer );
		return 0;
	}
	
	if ( link == gSpeciesList->lastItem ) {
		FreeIntVectorArray( composition, gAllowAtomList->items );
	}

	return 0;
}

Flag MatchIsomereName( LinkPtr link, SpeciesRecordHigh	*recordStruct )
{
	SpeciesPtr 			species = ( SpeciesPtr ) link->item;
	SpeciesRecordPtr	record = recordStruct->speciesRecordPtr;
	int					len = recordStruct->len;
	int					i;
	char				buffer[32];
	char				buffer2[32];
	char				*isomere;

	strcpy( buffer2, species->name );
	if ( !( isomere = strrchr( buffer2, '-' ) ) ) {
		return FALSE;
	}
	*isomere = '\0';
	strcpy( buffer, buffer2 );
	UpperString( buffer );

	for ( i = 0; i < len; ++i ) {
		strcpy( buffer2, record[i].name );
		UpperString( buffer2 );
		if ( strcmp( buffer, buffer2 ) == 0 ) {
			CopySpecies( species, &record[i] );
			fprintf( stderr, "### Thermo properties of '%s' are used for species '%s'\n", record[i].name, species->name );
			return TRUE;
		}
	}
	
	return FALSE;
}

void SortSpecies( void )
{
	int			i, nSpecies = gSpeciesList->items;
	ListPtr 	steadyStateList = NULL;
	SpeciesPtr	species = NULL;
	LinkPtr		link = NULL;
	
	steadyStateList = NewList( "Species" );
	
	/* copy steady states to new list and remove them from specieslist */
	link = gSpeciesList->firstItem;
	for ( i = 0; i < nSpecies; ++i ) {
		species = ( SpeciesPtr ) link->item;
		if ( species->isSteadyState ) {
			AddItem( species, steadyStateList );
			RemoveLinkFromList( gSpeciesList, link );
		}
		link = link->next;
	}

	if ( steadyStateList->items ) {
		/* append steady states to gSpeciesList */
		if ( gSpeciesList->firstItem ) {
			gSpeciesList->lastItem->next = steadyStateList->firstItem;
		}
		else {
			gSpeciesList->firstItem = steadyStateList->firstItem;
		}
		steadyStateList->firstItem->last = gSpeciesList->lastItem;
		gSpeciesList->lastItem = steadyStateList->lastItem;
		gSpeciesList->items += steadyStateList->items;
	}	
	FreeList( steadyStateList, NULL );
}

void CopySpecies( SpeciesPtr species, SpeciesRecordPtr record )
{
	int i;
	
	for ( i = 0; i < 7; ++i ) {
		species->coeffHot[i] = record->hot[i];
		species->coeffCold[i] = record->cold[i];
	}

	species->molarMass = record->M;
	species->k_over_eps = record->k_over_eps;
	species->sigma = record->sigma;
	species->muCoeff = record->mu_coeff;
	species->dCoeff = NewVector( gSpeciesList->items );
	species->omegaCoeff = NewVector( gSpeciesList->items );

/*	fprintf( stderr, "%s: M = %g\tk_over_eps = %g\tsigma = %g\n",*/
/*			species->name, species->molarMass, species->k_over_eps, species->sigma );*/
}

int GetComposition( char *name, IntVectorPtr composition )
{
	int		i;
	int		len = strlen( name );
	int		currentLoc = len-1;
	ListPtr singleChar = NULL;
	ListPtr doubleChar = NULL;
	char	buffer[3];

	buffer[0] = '\0';
	buffer[1] = '\0';
	buffer[2] = '\0';
	
	for ( i = 0; i < len-1; ++i ) {
		if ( name[currentLoc-1] == '-' ) {
			break;
		}
		else {
			--currentLoc;
		}
	}

	singleChar = NewList( "allowed atoms consisting of one char" );
	doubleChar = NewList( "allowed atoms consisting of two chars" );

	CreateSortedList( singleChar, 1 );
	CreateSortedList( doubleChar, 2 );

	while ( currentLoc < len ) {
		/* set buffer */
		if ( currentLoc < len-1 && isalpha( name[currentLoc+1] ) ) {
			/* shift two */
			buffer[0] = name[currentLoc];
			buffer[1] = name[currentLoc+1];
		}
		else {
			/* shift one */
			buffer[1] = buffer[0];
			buffer[0] = name[currentLoc];
		}

		/* first try to match doubleChar */
		if ( buffer[1] != '\0' ) {
			MatchBuffer( doubleChar, buffer, &currentLoc, name, composition );
		}

		if ( buffer[0] != '\0' ) {
			/* then if not matched try to match singleChar */
			buffer[1] = '\0';
			MatchBuffer( singleChar, buffer, &currentLoc, name, composition );
		}
		if ( buffer[0] != '\0' ) {
			/* species couldn't be matched -> initialize composition */
			for ( i = 0; i < composition->len; ++i ) {
				composition->vec[i] = 0;
			}
			return 0;
		}
	}

#	ifdef DEBUGCOMPOSITION
	for ( i = 0; i < composition->len; ++i ) {
		fprintf( stdout, "%s contains %d of atom number %d\n", name, composition->vec[i], i );
	}
	fprintf( stdout, "\n" );
#	endif

	FreeList( doubleChar, FreeString );
	FreeList( singleChar, FreeString );

	return 1;
}

int MatchBuffer( ListPtr list, String buffer, int *currentLoc, String name,
					IntVectorPtr composition )
{
	int 	i;
	int		index;
	int		number;
	LinkPtr link = NULL;
	String	atom = NULL;
	
	for ( i = 0; i < list->items; ++i ) {
		link = Find_n_th_Item( i+1, list );
		atom = ( String ) link->item;
		if ( strcmp( buffer, atom ) == 0 ) {
			/* matched */
			/* set pointer */
			*currentLoc += strlen( buffer );
			/* init buffer */
			buffer[1] = '\0';
			buffer[0] = '\0';
			/* find index */
			index = 0;
			while ( isdigit( name[*currentLoc] ) ) {
				index = index * 10 + name[(*currentLoc)++] - '0'; 
			}
			index = MAX( 1, index );
			/* find number of atom */
			if ( !( link = ListIterator( gAtomsList, FindAtom, atom ) ) ) {
				return -1;
			}
			number = ( ( AtomsPtr )link->item )->number;
				
			composition->vec[number] += index;
			break;
		}
	}
	return 0;
}

void CreateSortedList( ListPtr charList, int number )
{
	LinkPtr link = NULL;
	String	atom = NULL;
	int		i;

	for ( i = 0; i < gAllowAtomList->items; ++i ) {
		link = Find_n_th_Item( i+1, gAllowAtomList );
		atom = ( String ) link->item;
		if ( strlen( atom ) == number && ListIterator( gAtomsList, FindAtom, atom ) ) {
			AddString( charList, atom);
		}
	}
}

void ComputeOmegaDCoeff( SpeciesPtr species1, int j )
{
	LinkPtr link = NULL;
	SpeciesPtr	species2 = NULL;
	
	link = Find_n_th_Item( j+1, gSpeciesList );
	species2 = ( SpeciesPtr )link->item;
	
	/*
	 *	Computes the constant part of the diffusion coefficient approximation
	 */

	species1->dCoeff->vec[j] = 0.075320612 * sqrt( 1.0 / species1->molarMass + 1.0 / species2->molarMass ) / 
			( species1->sigma + species2->sigma ) / ( species1->sigma + species2->sigma );

	/*
	 *	Computes constant part of Sto§integral approximation for the diffusion 
	 *	coefficient
	 */
	species1->omegaCoeff->vec[j] = sqrt( species1->k_over_eps * species2->k_over_eps );
}

int ComputeDiffusionConstants( LinkPtr link, void *var )
{
	int	j;
	SpeciesPtr species = ( SpeciesPtr )link->item;
	ListPtr list = ( ListPtr ) var;
	
	for ( j = 0; j < list->items; ++j ) {
		ComputeOmegaDCoeff( species, j );
	}
	
	return 0;
}

SpeciesRecordPtr NewSpeciesRecordArray( int n )
{
	SpeciesRecordPtr	p = NULL;
	Double				dPtr = 0.0;
	 
	if ( !(p = (SpeciesRecordPtr)calloc( n, sizeof(SpeciesRecord) ) ) ) {
		fprintf( stderr, "Couldn't allocate memory for SpeciesRecordArray\n" );
	}
	
	return p;
}

/*
	FreeSpeciesRecordArray releases the memory which was previously
	allocated with NewSpeciesRecordArray.
*/

void FreeSpeciesRecordArray( SpeciesRecordPtr ptr )
{
	free( ptr );
}

IntVectorPtr *NewIntVectorArray( int number, int length )
{
	IntVectorPtr	*ptr = NULL;
	int				i;
	
	/* create array of IntVectorPtr */
	ptr = ( IntVectorPtr * ) calloc( number, sizeof( IntVectorPtr ) );
	
	for ( i = 0; i < number; ++i ) {
		ptr[i] = NewIntVector( length );
	}
	
	return ptr;
}

void FreeIntVectorArray( IntVectorPtr *ptr, int number )
{
	int i;
	
	for ( i = 0; i < number; ++i ) {
		DisposeIntVector( ptr[i] );
	}
	free( ptr );
}


void FatalError( const char *s )
{
	CleanUp();
	fprintf( stderr, "\n### Fatal Error: %s.\n", s );
	exit( 2 );
}


void Warning( const char *s )
{
	fprintf( stderr, "### Warning: %s.\n", s );
}

void SortSteadyStates( void )
{
	int 				i;
	int					nItems = NumberOfItems( gSpeciesList );
	IntVectorPtr		indexPtr = NewIntVector( gSpeciesList->items );
	int					*index = indexPtr->vec;
	LinkPtr				link;
	LinkPtr				nextLink;
	SpeciesPtr			species;
	
	if ( gSpeciesList->items == 0 ) {
		fprintf( stderr, "# error in SortSteadyStates: species list is empty\n" );
		exit( 2 );
	}
	
	/* sort species	*/
	link = gSpeciesList->firstItem;
	
	for ( i = 0; i < nItems; ++i, link = nextLink ) {
		nextLink = link->next;
		species = ( SpeciesPtr ) link->item;
		if ( species->isSteadyState ) {
			RemoveLinkFromList( gSpeciesList, link );
			AppendLink( link, gSpeciesList );
		}
	}
	
	/* build index table	*/
	link = gSpeciesList->firstItem;
	for ( i = 0; i < nItems; ++i ) {
		species = ( SpeciesPtr ) link->item;
		index[species->number] = i;
		link = link->next;
	}
	
	/* update lists	*/
	ListIterator( gSpeciesList, UpdateSpeciesList, index );
	ListIterator( gReactionList, UpdateReactionList, index );
	ListIterator( gThirdBodyList, UpdateThirdBodyList, index );
	ListIterator( gPAHReactionList, UpdateReactionList, index );
	ListIterator( gSootReactionList, UpdateReactionList, index );
}


int UpdateReactionList( LinkPtr link, void *var )
{
	int				*index = ( int* ) var;
	ReactionPtr		reaction = ( ReactionPtr ) link->item;
	int				i;
	
	for ( i = 0; i < reaction->numberOfSpecies; ++i ) {
		reaction->speciesNumber[i] = index[reaction->speciesNumber[i]];
	}
	
	return 0;
}

int UpdateSpeciesList( LinkPtr link, void *var )
{
	int				*index = ( int* ) var;
	SpeciesPtr		species = ( SpeciesPtr ) link->item;
	
	species->number = index[species->number];
	
	return 0;
}

int UpdateThirdBodyList( LinkPtr link, void *var )
{
	int				*index = ( int* ) var;
	ThirdBodyPtr	thirdbody = ( ThirdBodyPtr ) link->item;
	int				*specNum = thirdbody->speciesNumber->vec;
	Double			*specCoeff = thirdbody->speciesCoeff->vec;
	int				len = thirdbody->speciesNumber->len;
	int				i;
	SpeciesPtr		spec;
	
	for ( i = 0; i < len; ++i ) {
	  	specNum[i] = index[specNum[i]]; 
		spec = ( SpeciesPtr )Find_n_th_Item( specNum[i]+1, gSpeciesList )->item;
/*		fprintf( stderr, "Species '%s' is%s steady state\n"*/
/*				, spec->name, (spec->isSteadyState)?" ":" not" );*/
		if ( spec->isSteadyState ) {
			specCoeff[i] = 0.0;
		}
	}
	
	return 0;
}

void SetZeroProp (SpeciesPtr species)
{
	int i;

	for ( i = 0; i < 7; ++i ) {
 	  species->coeffHot[i] = 0;
	  species->coeffCold[i] = 0;
	}
	species->molarMass = 0.;
	species->k_over_eps = 0.;
	species->sigma = 0.;
	species->muCoeff = 0.;
	species->dCoeff = NewVector(gSpeciesList->items);
	species->omegaCoeff = NewVector(gSpeciesList->items);
}
