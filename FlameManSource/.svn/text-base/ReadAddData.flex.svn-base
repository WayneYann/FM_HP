%{  

#include "FlameMaster.h"
#include <ctype.h>
#undef YY_DECL
#define YY_DECL int FirstInput::yylex( void )

#undef DEBUG
#ifdef DEBUG
#undef YY_USER_ACTION
#define YY_USER_ACTION 
fprintf( stderr, "%s\n", ( char * )yytext  );
#endif
%}

DOUBLE		[0-9\.Ee\-]+
INT			[0-9\-]+
SPECIES		[0-9A-Za-z\-]+
STRING		[^\n]+

%x scComment
%x scBC
%x scBCMode
%x scGetDouble
%x scGetInt
%x scGetFlag
%x scGetString
%x scSetSpeciesAddress
%x scReadSpecies
%x scReadFlameType
%x scGlobalReaction

%%
	BoundaryInputPtr	currentBoundary;
	int					bcCond;
	Double				*dAddress = NULL;
	int					*iAddress = NULL;
	char				**sAddress = NULL;
	Flag				*flagAddress = NULL;
	int					startCondition;
	Flag				skipWhites = TRUE;
	Flag				toUpper = TRUE;

	yy_did_buffer_switch_on_eof = yy_did_buffer_switch_on_eof; // suppress warning
%{
/* eat up whitespace  */
%}
<scBC,scBCMode,scGlobalReaction,scGetDouble,scGetInt,scGetFlag,scGetString,scSetSpeciesAddress,scReadSpecies,scReadFlameType,INITIAL>[ \t\n]+

%{
/* comments  */
%}
<scBC>#	{
													startCondition = scBC;
													BEGIN(scComment);
									}
<scBCMode>#	{
													startCondition = scBCMode;
													BEGIN(scComment);
									}
<scGetDouble>#	{
													startCondition = scGetDouble;
													BEGIN(scComment);
									}
<scGetInt>#	{
													startCondition = scGetInt;
													BEGIN(scComment);
									}
<scGetFlag>#	{
													startCondition = scGetFlag;
													BEGIN(scComment);
									}
<scGetString>#	{
													startCondition = scGetString;
													BEGIN(scComment);
									}
<scSetSpeciesAddress>#	{
													startCondition = scSetSpeciesAddress;
													BEGIN(scComment);
									}
<scReadSpecies>#	{
													startCondition = scReadSpecies;
													BEGIN(scComment);
									}
<scReadFlameType>#	{
													startCondition = scReadFlameType;
													BEGIN(scComment);
									}
<scGlobalReaction>#	{
													startCondition = scGlobalReaction;
													BEGIN(scComment);
									}
<INITIAL>#	{
													startCondition = INITIAL;
													BEGIN(scComment);
									}
<scComment>\n     	{ 
						BEGIN(startCondition);
                	}
<scComment>[^\n]  

%{
/* INITIAL  */
%}
<INITIAL>"/*"		{
		int c;
		
		for( ; ; ) {
			while ( ( c = yyinput() ) != '*' && c != EOF ) /* null statement */;
			if ( c == '*' ) {
				while ( ( c = yyinput() ) == '*' )  /* null statement */;
				if ( c == '/' ) {
					break;
				}
			}
		}
	}


<INITIAL>InitialCond(ition)?[ \t\n]*"{"	{
		currentBoundary = fInitialCond;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scBCMode)"  );
		}
		BEGIN(scBCMode);
	}

<INITIAL>(leftb?c?[ \t\n]*"{")|(Fuel[ \t\n]*Side[ \t\n]*"{")|(Far[ \t\n]*Field[ \t\n]*"{")|(Unburnt[ \t\n]*Side[ \t\n]*"{")		{
		currentBoundary = leftBoundary;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scBC)"  );
		}
		BEGIN(scBC);
	}

<INITIAL>(rightb?c?[ \t\n]*"{")|(Oxidizer[ \t\n]*Side[ \t\n]*"{")|(Particle[ \t\n]*"{")	{
		currentBoundary = rightBoundary;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scBC)"  );
		}
		BEGIN(scBC);
	}

%{
/* scBC  */
%}

<scBC>"}"	{
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}

<scBC>(gradient[ \t\n]*"{")|(neumann[ \t\n]*"{")	{
		bcCond = kNeumann;		
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scBCMode)"  );
		}
		BEGIN(scBCMode);
	}

<scBC>dirichlet[ \t\n]*"{"	{
		bcCond = kDirichlet;		
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scBCMode)"  );
		}
		BEGIN(scBCMode);
	}


%{
/* scBCMode  */
%}

<scBCMode>"}"	{
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scBC)"  );
		}
		BEGIN(scBC);
	}

<scBCMode>U|f"'"	{
		currentBoundary->fBcFlag[fUVelocityOffset] = bcCond;
		dAddress = &currentBoundary->fValue[fUVelocityOffset];
		startCondition = scBCMode;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}

<scBCMode>V|f	{
		currentBoundary->fBcFlag[fVVelocityOffset] = bcCond;
		dAddress = &currentBoundary->fValue[fVVelocityOffset];
		startCondition = scBCMode;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}

<scBCMode>T	{
		currentBoundary->fBcFlag[fTemperatureOffset] = bcCond;
		dAddress = &currentBoundary->fValue[fTemperatureOffset];
		startCondition = scBCMode;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}


<scBCMode>p	{
		currentBoundary->fBcFlag[fPressureOffset] = bcCond;
		dAddress = &currentBoundary->fValue[fPressureOffset];
		startCondition = scBCMode;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}

<scBCMode>y->	{
		if ( currentBoundary->fMixtureSpecification == 0 || currentBoundary->fMixtureSpecification == kMassFraction ) {
			currentBoundary->fMixtureSpecification = kMassFraction;
		}
		else {
			fprintf( stderr, "%s%d\n", "error: mixture at left boundary now specified as MassFraction has previously specified with value ", currentBoundary->fMixtureSpecification  );
			exit(2);
		}
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scSetSpeciesAddress)"  );
		}
		BEGIN(scSetSpeciesAddress);
	}

<scBCMode>(massflux->)|(epsilon->)	{
		if ( currentBoundary->fMixtureSpecification == 0 || currentBoundary->fMixtureSpecification == kMassFlux ) {
			currentBoundary->fMixtureSpecification = kMassFlux;
		}
		else {
			fprintf( stderr, "%s%d\n", "error: mixture at left boundary now specified as MassFlux has previously specified with value ", currentBoundary->fMixtureSpecification  );
			exit(2);
		}
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scSetSpeciesAddress)"  );
		}
		BEGIN(scSetSpeciesAddress);
	}

<scBCMode>x->	{
		if ( currentBoundary->fMixtureSpecification == 0 || currentBoundary->fMixtureSpecification == kMolarFraction ) {
			currentBoundary->fMixtureSpecification = kMolarFraction;
		}
		else {
			fprintf( stderr, "%s%d\n", "error: mixture at left boundary now specified as MolarFraction has previously specified with value ", currentBoundary->fMixtureSpecification  );
			exit(2);
		}
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scSetSpeciesAddress)"  );
		}
		BEGIN(scSetSpeciesAddress);
	}

%{
/* scGlobalReaction  */
%}

<scGlobalReaction>[^#]+;	{
   		if ( !( fGlobalReaction = new char[yyleng + 1] ) ) {
			FatalError( "memory allocation in yylex failed" );
		}
		strcpy( fGlobalReaction, ( char * ) yytext );
		fGlobalReaction[yyleng-1] = '\0';
		SkipWhites( fGlobalReaction );
		UpperString( fGlobalReaction );
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}

	
%{
/* scSetSpeciesAddress  */
%}

<scSetSpeciesAddress>{SPECIES}	{
		int number = currentBoundary->fSpecifiedSpeciesBCs;
		char *name = currentBoundary->speciesName[number];

		if ( currentBoundary->fBcFlagSpecies == kNone || currentBoundary->fBcFlagSpecies == bcCond ) {
			currentBoundary->fBcFlagSpecies = bcCond;
		}
		else {
			fprintf( stderr, "%s%d%s%d\n", "error: boundary now specified as ", bcCond, " has previously specified with value ", currentBoundary->fBcFlagSpecies  );
		}
		strcpy( name, ( char * ) yytext );
		SkipWhites( name );
		UpperString( name );
		dAddress = &currentBoundary->fValueSpecies[number];

		++currentBoundary->fSpecifiedSpeciesBCs;
		startCondition = scBCMode;

		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}


%{
/* scGetDouble  */
%}

<scGetDouble>=|is /* skip */

<scGetDouble>{DOUBLE}	{
		*dAddress = atof( ( char * ) yytext );
		if ( fScannerProgress ) {
			fprintf( stderr, "%s%d%s\n", "BEGIN(", startCondition, ")"  );
		}
		BEGIN(startCondition);
	}


%{
/* scGetInt  */
%}

<scGetInt>=|is /* skip */

<scGetInt>{INT}	{
		*iAddress = atoi( ( char * ) yytext );
		if ( fScannerProgress ) {
			fprintf( stderr, "%s%d%s\n", "BEGIN(", startCondition, ")"  );
		}
		BEGIN(startCondition);
	}

%{
/* scGetFlag  */
%}

<scGetFlag>= /* skip */
<scGetFlag>is /* skip */

<scGetFlag>TRUE	{
		*flagAddress = TRUE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s%d%s\n", "BEGIN(", startCondition, ")"  );
		}
		BEGIN(startCondition);
	}
<scGetFlag>FALSE	{
		*flagAddress = FALSE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s%d%s\n", "BEGIN(", startCondition, ")"  );
		}
		BEGIN(startCondition);
	}

%{
/* scGetString  */
%}

<scGetString>{STRING}	{
		if ( yyleng > 127 ) {
			fprintf( stderr, "%s\n", "allowed length of a string is 127 chars"  );
			exit( 2 );
		}
		strcpy( stringBuffer, ( char * ) yytext );
		if ( skipWhites ) {
			SkipWhites( stringBuffer );
		}
		else {
			SkipIntroWhites( stringBuffer );
		}
		if ( toUpper ) {
			UpperString( stringBuffer );
		}
		*sAddress = new char[ strlen( stringBuffer ) + 1 ];
		strcpy( *sAddress, stringBuffer );
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}

%{
/* scReadSpecies  */
%}

<scReadSpecies>{SPECIES}	{
		if ( yyleng > 127 ) {
			fprintf( stderr, "%s\n", "allowed length of a string is 127 chars"  );
			exit( 2 );
		}
		strcpy( stringBuffer, ( char * ) yytext );
		SkipWhites( stringBuffer );
		UpperString( stringBuffer );
		*sAddress = new char[ strlen( stringBuffer ) + 1 ];
		strcpy( *sAddress, stringBuffer );
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}

%{
/* scReadFlameType  */
%}

<scReadFlameType>(Counterflow)?[\n\t ]*Diffusion[\n\t ]*in[\n\t ]*Mixt?u?r?e?[\n\t ]*Fract?i?o?n?[\n\t ]*Space	{
		fFlameType = kCountDiffMix;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>Counterflow[\n\t ]*Premixed[\n\t ]*((with)?|(in)?|(on)?)[\n\t ]*Physical[\n\t ]*Coordinate?	{
		fFlameType = kCountPremPhys;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>Counterflow[\n\t ]*Premixed[\n\t ]*((with)?|(in)?|(on)?)[\n\t ]*Similarity[\n\t ]*Coordinate?	{
		fFlameType = kCountPremSim;
		fPremConfiguration = kFull;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>Counterflow[\n\t ]*Premixed[\n\t ]*Back[\n\t ]*to[\n\t ]*Back[\n\t ]*((with)?|(in)?|(on)?)[\n\t ]*Similarity[\n\t ]*Coordinate?	{
		fFlameType = kCountPremSim;
		fPremConfiguration = kBackToBack;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>Counterflow[\n\t ]*Premixed[\n\t ]*Fresh[\n\t ]*to[\n\t ]*Burned[\n\t ]*((with)?|(in)?|(on)?)[\n\t ]*Similarity[\n\t ]*Coordinate?	{
		fFlameType = kCountPremSim;
		fPremConfiguration = kFreshToBurned;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>Flamesheet	{
		fFlameType = kStartUpDiffusion;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>Diffusion[\n\t ]*((with)?|(in)?|(on)?)[\n\t ]*Physical[\n\t ]*Coordinate?	{
		fFlameType = kDiffusionPhys;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>EigenValueDiffusion[\n\t ]*((with)?|(in)?|(on)?)[\n\t ]*Physical[\n\t ]*Coordinate?	{
		fFlameType = kDiffPhysEigen;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>CounterFlowDiffusion	{
		fFlameType = kStartUpDiffusion2;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>CounterFlowPremixed	{
		fPremConfiguration = 1;
		fFlameType = kStartUpDiffusion2;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>CounterFlowDiffusion[\n\t ]*(with)?[\n\t ]*Continuation	{
		fFlameType = kCountDiffCont;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>UnstretchedPremixed	{
		fFlameType = kUnstrPremPhys;
		if ( fScannerProgress ) {
			cerr << "BEGIN(INITIAL)" << NEWL;
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>Isochor[\n\t ]*Homo(geneous)?[\n\t ]*Reactor	{
		fFlameType = kHomoIsoChor;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>Isobar[\n\t ]*Homo(geneous)?[\n\t ]*Reactor	{
		fFlameType = kHomoIsoBar;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>Transient[\t\n ]*Flamelet	{
		fFlameType = kTransFlamelet;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>Transient[\t\n ]*1DIsoChor	{
		fFlameType = kTrans1DIsoChor;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}
<scReadFlameType>Transient[\t\n ]*Coal[\t\n ]*Particle	{
		fFlameType = kTransCoalParticle;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(INITIAL)"  );
		}
		BEGIN(INITIAL);
	}

%{
/* additional input  */
%}

<INITIAL>ExactBackward	{
		flagAddress = &fExactBackward;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>AdditionalOutput	{
		flagAddress = &fAdditionalOutput;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Equidistant(Output)?	{
		flagAddress = &fEquidistant;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>ArclengthCont(in)?(uation)?	{
		flagAddress = &fArcLengthContin;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>StrainRateCont(in)?(uation)?	{
		flagAddress = &fStrainRateContin;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>N(umber)?(Of)?Outputs	{
		iAddress = &fNOutputs;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetInt)"  );
		}
		BEGIN(scGetInt);
	}
<INITIAL>Art(ificial)?Source	{
		dAddress = &fArtificialSource;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Art(ificial)?SourceTime	{
		dAddress = &fArtSourceTime;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Const(ant)?LewisNumber	{
		flagAddress = &fConstantLewisNumber;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>(Compute)?WithRadiation	{
		flagAddress = &fWithRadiation;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>ThermoDiffusion	{
		flagAddress = &fThermoDiffusion;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>(Compute)?WithSoot	{
		flagAddress = &fWithSoot;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Nucleation	{
		flagAddress = &fNucleation;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Condensation	{
		flagAddress = &fCondensation;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>SootRadiation	{
		flagAddress = &fSootRadiation;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>SootUp(date)?ProdRate(s)?	{
		flagAddress = &fSootUpdateProdRates;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>SizeDep(endent)?Diff(usion)?	{
		flagAddress = &fSizeDepDiff;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Surf(ace)?Dep(endent)?Coag(ulation)?	{
		flagAddress = &fSurfDepCoag;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Coagulation	{
		flagAddress = &fCoagulation;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>SurfaceGrowth	{
		flagAddress = &fSurfaceGrowth;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>SurfaceOxidation	{
		flagAddress = &fSurfaceOxidation;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>ThermoPhoresis	{
		flagAddress = &fThermoPhoresis;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>PAHOHOxidation	{
		flagAddress = &fOHPAHOxidation;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>PAHO2Oxidation	{
		flagAddress = &fO2PAHOxidation;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Coag(ulation)?Fact(or)?	{
		dAddress = &fCoagFact;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>NSootMoments	{
		iAddress = &fNSootMoments;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetInt)"  );
		}
		BEGIN(scGetInt);
	}
<INITIAL>FlameIsAxiSymmetric	{
		flagAddress = &fIsAxiSymmetric;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>NoDiff(usivity)?Corr(ection)?	{
		flagAddress = &fNoDiffCorr;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>ClipNegativeConc(entration)?s	{
		flagAddress = &fClipNegativeConcs;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}



<INITIAL>WriteBT	{
		flagAddress = &fWriteBT;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>WriteResi?d?u?u?m?	{
		flagAddress = &fWriteResiduum;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>WatchGridding	{
		flagAddress = &fWatchGridding;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>WriteEverySol(ution)?	{
		flagAddress = &fWriteEverySolution;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>OutputPath[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fOutputPath;
		skipWhites = TRUE;
		toUpper = FALSE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetString)"  );
		}
		BEGIN(scGetString);
	}



<INITIAL>WriteFullResi?d?u?a?l?s?	{
		flagAddress = &fWriteFullRes;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>ComputeUnphysicalChain|ConstMassFlux	{
		flagAddress = &fCompUnPhysChain;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>(NOfUnphys(ical)?Chain)|FirstSensRate	{
		iAddress = &fNUnPhysChain;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetInt)"  );
		}
		BEGIN(scGetInt);
	}
<INITIAL>(Use)?ModifiedNewton	{
		flagAddress = &fUseModifiedNewton;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>(Use)?NumericalJac	{
		flagAddress = &fUseNumericalJac;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			cerr << "BEGIN(scGetFlag)" << NEWL;
		}
		BEGIN(scGetFlag);
	}
<INITIAL>(Use)?Sec(ond)?Ord(er)?Jac	{
		flagAddress = &fUseSecOrdJac;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			cerr << "BEGIN(scGetFlag)" << NEWL;
		}
		BEGIN(scGetFlag);
	}
<INITIAL>(Use)?NumericalDM	{
		flagAddress = &fUseNumericalDM;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Sens(itivity)?Obj(ect)?[ \t\n]*(is|=)	{
		if ( fNSensObj >= 1000 ) {
			FATAL( "too many sens objs" );
		}
		sAddress = &fSensObj[fNSensObj++];
		skipWhites = TRUE;
		toUpper = TRUE;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetString)"  );
		}
		BEGIN(scGetString);
	}
<INITIAL>ReactionFluxAnal(ysis)?	{
		flagAddress = &fReactionFluxes;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>PrintRHSSpecies	{
		flagAddress = &fPrintRHSSpecies;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>PrintRHSTemp(erature)?	{
		flagAddress = &fPrintRHSTemp;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>KeepMassFrac(tion)?s	{
		flagAddress = &fKeepMassFracs;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>LiquidPool(BC)?	{
		flagAddress = &fLiquidPoolBC;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>StrainRate|MassFlowRate	{
		if ( nStrainRates >= 100 ) {
			FATAL( "too many strainrates" );
		}
		dAddress = &fStrainRate[nStrainRates++];
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Pressure	{
		if ( fNPressures >= 100 ) {
			FATAL( "too many strainrates" );
		}
		dAddress = &fPressure[fNPressures++];
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Phi	{
		if ( fNPhi >= 100 ) {
			FATAL( "too many equivalence ratios" );
		}
		dAddress = &fPhi[fNPhi++];
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>((Scalar)?[ \t]*Diss(ipation)?[ \t]*Rate)|Chi	{
		dAddress = &fDissRate[nDissRates++];
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>MinStrainRate|MassFluxLiquidPool	{
		dAddress = &fMinStrainRate;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>DeltaSCont(inuation)?	{
		dAddress = &fDeltaSCont;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>InitialEquation?s	{
		iAddress = &fInitialEquations;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetInt)"  );
		}
		BEGIN(scGetInt);
	}
<INITIAL>maxgridpoints	{
		iAddress = &fMaxGridPoints;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetInt)"  );
		}
		BEGIN(scGetInt);
	}
<INITIAL>(Initial)?GridPoints	{
		iAddress = &fInitialGridPoints;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetInt)"  );
		}
		BEGIN(scGetInt);
	}
<INITIAL>DampFlag	{
		flagAddress = &fDampFlag;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>ContinFlag	{
		flagAddress = &fContinFlag;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>TimeDepFlag	{
		flagAddress = &fTimeDepFlag;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>DeltaNewGrid	{
		iAddress = &fDeltaNewGrid;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetInt)"  );
		}
		BEGIN(scGetInt);
	}
<INITIAL>TolRes|RelTol	{
		dAddress = &fTolRes;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>TolDy|AbsTol	{
		dAddress = &fTolDy;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>MaxIter	{
		iAddress = &fMaxIter;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetInt)"  );
		}
		BEGIN(scGetInt);
	}
<INITIAL>TStart	{
		dAddress = &fStart;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>TEnd	{
		dAddress = &fEnd;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Left	{
		dAddress = &fLeft;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Right	{
		dAddress = &fRight;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Gamma {
		iAddress = &fGamma;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetInt)"  );
		}
		BEGIN(scGetInt);
	}
<INITIAL>Kappa	{
		dAddress = &fKappa;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Tau	{
		dAddress = &fTauGrid;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>CL	{
		dAddress = &fCl;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Hv	{
		dAddress = &fHv;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Tb	{
		dAddress = &fT_B;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>OneSolu?t?i?o?n?OneGrid	{
		flagAddress = &fOneSolOneGrid;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>AdjustComp(utational)?Domain	{
		flagAddress = &fAdjustComputationalDomain;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>PrintMol(e|(ar))Fractions?	{
		flagAddress = &fPrintMolarFractions;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Sens(itivity)?Anal(ysis)?	{
		flagAddress = &fSensAnal;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Sens(itivity)?Anal(ysis)?Spec(ies)?	{
		flagAddress = &fSensAnalSpec;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Sens(itivity)?Anal(ysis)?Reac(tion)?	{
		flagAddress = &fSensAnalReac;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Sens(itivity)?Anal(ysis)?Fac(tor)?	{
		dAddress = &fSensAnalFac;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Sens(itivity)?Obj(ect)?All	{
		flagAddress = &fSensObjAll;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Sens(itivity)?Max	{
		flagAddress = &fSensMax;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>Sens(itivity)?Final	{
		flagAddress = &fSensFinal;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>SteadyStatesNumerical	{
		flagAddress = &fSteadyStatesNumerical;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetFlag)"  );
		}
		BEGIN(scGetFlag);
	}
<INITIAL>R	{
		dAddress = &fR;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Q	{
		dAddress = &fQ;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>GridCorr(ection)?start	{
		dAddress = &fStart;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>GridCorr(ection)?End	{
		dAddress = &fEnd;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>GridCorr(ection)?Alpha	{
		dAddress = &fAlpha;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>LambdaMin	{
		dAddress = &fLambdaMin;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>DeltaTMax	{
		dAddress = &fDeltaTMax;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>DeltaTStart	{
		dAddress = &fDeltaTStart;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>ContSteps	{
		iAddress = &fContSteps;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetInt)"  );
		}
		BEGIN(scGetInt);
	}
<INITIAL>Date[ \t\n]*Created(is|=)	{
		startCondition = INITIAL;
		sAddress = &fDateCreated;
		skipWhites = FALSE;
		toUpper = FALSE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetString)"  );
		}
		BEGIN(scGetString);
	}
<INITIAL>author[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fAuthor;
		skipWhites = FALSE;
		toUpper = FALSE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetString)"  );
		}
		BEGIN(scGetString);
	}
<INITIAL>fuel[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fFuelName[fNFuels++];
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scReadSpecies)"  );
		}
		BEGIN(scReadSpecies);
	}
<INITIAL>(ToSpecies|GotoFuel)[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fToSpeciesName;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scReadSpecies)"  );
		}
		BEGIN(scReadSpecies);
	}
<INITIAL>(FromSpecies)[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fFromSpeciesName;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scReadSpecies)"  );
		}
		BEGIN(scReadSpecies);
	}
<INITIAL>Cont(inuation)?Bound	{
		dAddress = &fContBound;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>Cont(inuation)?Inc	{
		dAddress = &fContInc;
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetDouble)"  );
		}
		BEGIN(scGetDouble);
	}
<INITIAL>oxidizer[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fOxName;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scReadSpecies)"  );
		}
		BEGIN(scReadSpecies);
	}
<INITIAL>MechanismFile[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fReactionFile;
		skipWhites = TRUE;
		toUpper = FALSE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetString)"  );
		}
		BEGIN(scGetString);
	}
<INITIAL>StartProfilesFile[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fStartProfileFile;
		skipWhites = TRUE;
		toUpper = FALSE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetString)"  );
		}
		BEGIN(scGetString);
	}
<INITIAL>(CAInFile|AddInFile)[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fAddFileNo1;
		skipWhites = TRUE;
		toUpper = FALSE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetString)"  );
		}
		BEGIN(scGetString);
	}
<INITIAL>LewisNumberFile[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fLewisNumberFile;
		skipWhites = TRUE;
		toUpper = FALSE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetString)"  );
		}
		BEGIN(scGetString);
	}
<INITIAL>outputfile[ \t\n]*(is|=) {
		startCondition = INITIAL;
		sAddress = &fOutFileName;
		skipWhites = TRUE;
		toUpper = FALSE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetString)"  );
		}
		BEGIN(scGetString);
	}
<INITIAL>FlameT?y?p?e?[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scReadFlameType)"  );
		}
		BEGIN(scReadFlameType);
	}
<INITIAL>Cont(in)?(uation)?Type[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fContinType;
		skipWhites = TRUE;
		toUpper = TRUE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetString)"  );
		}
		BEGIN(scGetString);
	}
<INITIAL>Cont(in)?(uation)?Side[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		sAddress = &fContinSide;
		skipWhites = TRUE;
		toUpper = TRUE;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGetString)"  );
		}
		BEGIN(scGetString);
	}
<INITIAL>GlobalReaction[ \t\n]*(is|=)	{
		startCondition = INITIAL;
		if ( fScannerProgress ) {
			fprintf( stderr, "%s\n", "BEGIN(scGlobalReaction)"  );
		}
		BEGIN(scGlobalReaction);
	}

<INITIAL,scSetSpeciesAddress,scReadSpecies,scReadFlameType,scGlobalReaction>.   	{
					fprintf( stderr, "a '%s' is a complete surprise for me in this context\n", ( char * )yytext );
				}
<INITIAL,scGetString,scGetFlag,scGetInt,scGetDouble,scBCMode,scBC,scComment>.   	{
					fprintf( stderr, "a '%s' is a complete surprise for me in this context\n", ( char * )yytext );
				}
%{
/*  end of file  */
%}
<<EOF>>             { yyrestart( yyin ); return 0; }
%%

int TInputData::GetSpeciesNumber( void )
{
	int		addressNumber;
	char	species[32];
		
	if ( yyleng > 31 ) {
		fprintf( stderr, "%s%s\n", "error: length of species greater than 31 chars:", ( char * )yytext  );
		exit( 2 );
	}
	else {
		strcpy( species, ( char * ) yytext );
		SkipWhites( species );
		UpperString( species );
	}
	if ( ( addressNumber = FindSpecies( species ) ) == -1 ) {
		fprintf( stderr, "%s%s\n", "error: can't find specified species ", species  );
		exit( 2 );
	}
	return addressNumber;
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

void SkipIntroWhites( char * label )
{
	unsigned int	i, j = 0, length = 0;
	Flag			letterFound = FALSE;
	length = strlen( label );
	for ( i = 0; i < strlen( label ); ++i ) {
		if ( !letterFound && isspace( label[i] ) ) {
			--length;
		}
		else {
			letterFound = TRUE;
			label[j++] = label[i];
		}
	}
	label[length] = '\0';
}

void UpperString( char *string )
{
	while ( *string ) {
		*string = toupper( *string );
		++string;
	}
}

int yywrap()
	{
	return 1;
	}
