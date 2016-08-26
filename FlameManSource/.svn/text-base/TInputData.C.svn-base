#include "FlameMaster.h"

#undef DEBUG

#define DEFAULTINPUTFILE "FlameMaster.input"

void BoundaryInput::InitBoundaryInput( int len )
{
	int	i;
	fLen = len;
	fSpecifiedSpeciesBCs = 0;
	fMixtureSpecification = 0;
	fBcFlagSpecies = kNone;

   	if ( !( speciesName = new String[fLen] ) ) FatalError( "memory allocation of BoundaryInput failed" );
   	if ( !( fValueSpecies = new Double[fLen] ) ) FatalError( "memory allocation of BoundaryInput failed" );
   	if ( !( fBcFlag = new int[fNVars] ) ) FatalError( "memory allocation of BoundaryInput failed" );
   	if ( !( fValue = new Double[fNVars] ) ) FatalError( "memory allocation of BoundaryInput failed" );

	for ( i = 0; i < fLen; ++i ) {
	   	if ( !( speciesName[i] = new char[31] ) ) FatalError( "memory allocation of BoundaryInput failed" );
		fValueSpecies[i] = 0.0;
	}
	for ( i = 0; i < fNVars; ++i ) {
		fBcFlag[i] = kNone;
		fValue[i] = 0.0;
	}
}

BoundaryInput::~BoundaryInput( void )
{
	for ( int i = 0; i < fLen; ++i ) {
		delete speciesName[i];
	}

	delete fValue;
	delete fBcFlag;
	delete fValueSpecies;
	delete speciesName;
	//	delete fBcFlagSpecies;
}

BoundaryInput& BoundaryInput::operator=( const BoundaryInput& boundary )
{
// take fLen of first argument
	int	i;

	fSpecifiedSpeciesBCs = boundary.fSpecifiedSpeciesBCs;
	fMixtureSpecification = boundary.fMixtureSpecification;
	fBcFlagSpecies = boundary.fBcFlagSpecies;
	
	for ( i = 0; i < fLen; ++i ) {
		strcpy( speciesName[i], boundary.speciesName[i] );
		fValueSpecies[i] = boundary.fValueSpecies[i];
	}
	for ( i = 0; i < fNVars; ++i ) {
		fBcFlag[i] = boundary.fBcFlag[i];
		fValue[i] = boundary.fValue[i];
	}

	return *this;
}

void BoundaryInput::PrintBoundary( FILE *fp )
{
	int		i;
	char	spec[31];

	if ( fBcFlag[fUVelocityOffset] == 1 ) {
		fprintf( fp, "U = %g\n", fValue[fUVelocityOffset] );
	}
	else if ( fBcFlag[fUVelocityOffset] == 2 ) {
		fprintf( fp, "U' = %g\n", fValue[fUVelocityOffset] );
	}
	else {
		fprintf( fp, "U not specified\n" );
	}

	if ( fBcFlag[fVVelocityOffset] == 1 ) {
		fprintf( fp, "V = %g\n", fValue[fVVelocityOffset] );
	}
	else if ( fBcFlag[fVVelocityOffset] == 2 ) {
		fprintf( fp, "V' = %g\n", fValue[fVVelocityOffset] );
	}
	else {
		fprintf( fp, "V not specified\n" );
	}

	if ( fBcFlag[fTemperatureOffset] == 1 ) {
		fprintf( fp, "T = %g\n", fValue[fTemperatureOffset] );
	}
	else if ( fBcFlag[fTemperatureOffset] == 2 ) {
		fprintf( fp, "T' = %g\n", fValue[fTemperatureOffset] );
	}
	else {
		fprintf( fp, "T not specified\n" );
	}

	if ( fMixtureSpecification == kMassFlux ) {
		strcpy( spec, "epsilon" );
	}
	else if ( fMixtureSpecification == kMassFraction ) {
		strcpy( spec, "Y" );
	}
	else if ( fMixtureSpecification == kMolarFraction ) {
		strcpy( spec, "X" );
	}
	for ( i = 0; i < fSpecifiedSpeciesBCs; ++i ) {
		if ( fBcFlagSpecies == 1 ) {
			fprintf( fp, "%s->%s = %g\n", spec, speciesName[i], fValueSpecies[i] );
		}
		else if ( fBcFlagSpecies == 2 ) {
			fprintf( fp, "%s'->%s = %g\n", spec, speciesName[i], fValueSpecies[i] );
		}
		else {
			fprintf( fp, "%s not specified\n", speciesName[i] );
		}
	}
}

void FirstInput::InitFirstInput( int argc, char *argv[] )
{
	fprintf( stderr, "executing: '%s", argv[0] );
	for ( int jj = 1; jj < argc; ++jj ) {
		fprintf( stderr, " %s", argv[jj] );
	}
	fprintf( stderr, "'\n" );

	fprintf( stderr, "\nFlameMaster Version %s written by Heinz Pitsch\n\n", (char *)VERSION );

	// Exact Backward Reaction Constants
	fExactBackward = FALSE;

	// 0D-stuff
	fAdditionalOutput = FALSE;
	fEquidistant = FALSE;
	fNOutputs = -1;
	fArtificialSource = 0.0;
	fArtSourceTime = 0.0;
	fpA = NULL;
	fpS = NULL;
	fGlobalReaction = NULL;
	fScannerProgress = FALSE;
	fAddFileNo1 = NULL;
	fLewisNumberFile = NULL;
	fMechanismFileComm = NULL;
	fStartProfileFileComm = NULL;
	fPressureComm = -1.0e5;
	fParameterComm = -1;
	fPremConfiguration = -1;
	fCoagFact = 1.0;
	fOutputPathComm = NULL;
	
	if ( argc ) {
		ParseCommandLine( argc, argv );
	}
	else {
		if ( strlen( *argv ) > 0 ) {
			strcpy( fAdditionalInfoFile, *argv );
			if ( ( fpA = fopen( fAdditionalInfoFile, "r" ) ) == NULL ){
				fprintf( stderr, "#error: can't open specified input file '%s'\n", fAdditionalInfoFile);
				exit( 2 );
			}
			else {
				yyin = fpA;
			}
		}
		else {
			if ( ( fpA = fopen( DEFAULTINPUTFILE, "r" ) ) == 0 ) {
				fprintf( stderr, "#error: can't open specified input file '%s'\n", DEFAULTINPUTFILE);
				exit( 2 );
			}
			else {
				strcpy( fAdditionalInfoFile, DEFAULTINPUTFILE );
				cerr << "use input file '" << DEFAULTINPUTFILE << "'"  << NEWL;
			}
		}
	}
	
   	if ( !( fInitialCond = new BoundaryInput( this, 100 ) ) ) FatalError( "memory allocation of BoundaryInput failed" );
   	if ( !( leftBoundary = new BoundaryInput( this, 100 ) ) ) FatalError( "memory allocation of BoundaryInput failed" );
   	if ( !( rightBoundary = new BoundaryInput( this, 100 ) ) ) FatalError( "memory allocation of BoundaryInput failed" );
	
// set default values of
	// T1DFlame

	fDateCreated = "";
	fAuthor = "";
	fFlameType = kNotSpecified;
	fContinType = NULL;
	fContinSide = NULL;
	fContBound = 0.0;
	fConstantLewisNumber = FALSE;
	fWithRadiation = FALSE;
	fArcLengthContin = FALSE;
	fStrainRateContin = FALSE;
	fThermoDiffusion = FALSE;
	fWithSoot = FALSE;
	fNSootMoments = 0;
	fCompUnPhysChain = FALSE;
	fNUnPhysChain = -1;
	fIsAxiSymmetric = FALSE;
	fClipNegativeConcs = TRUE;
	fNoDiffCorr = FALSE;
	fWriteFullRes = FALSE;
	fUseModifiedNewton = FALSE;
	//	fStrainRate = 10.0;

	// Sensitivity Analysis
	fSensAnal = FALSE;
	fSensObjAll = FALSE;
	fSensAnalSpec = FALSE;
	fSensAnalReac = FALSE;
	fSensAnalFac = 0.0;  
	fSensMax = FALSE;
	fSensFinal = FALSE;
	fNSensObj = 0;
	fReactionFluxes = FALSE;

	fPrintRHSSpecies = FALSE;
	fPrintRHSTemp = FALSE;
	nStrainRates = 0;
	nDissRates = 0;
	fKeepMassFracs = FALSE;
	fLiquidPoolBC = FALSE;
	fDeltaSCont = 20.0;
	fMinStrainRate = -1.0;
	fNPressures = 0;
	fNPhi = 0;
	fFromSpeciesName = NULL;
	fToSpeciesName = NULL;
	fOutFileName = NULL;
	fReactionFile = NULL;
	fStartProfileFile = NULL;
	fContInc = 0.0;
	fPrintMolarFractions = FALSE;
	fSteadyStatesNumerical = FALSE;
	fUseNumericalJac = TRUE;
	fUseSecOrdJac = FALSE;
	fUseNumericalDM = FALSE;
	fNFuels = 0;
	
	fNucleation = TRUE;
	fCondensation = TRUE;
	fCoagulation = TRUE;
	fSurfaceGrowth = TRUE;
	fSurfaceOxidation = TRUE;
	fThermoPhoresis = TRUE;
	fOHPAHOxidation = TRUE;
	fO2PAHOxidation = TRUE;
	fSootRadiation = TRUE;
	fSootUpdateProdRates = TRUE;
	fSizeDepDiff = TRUE;
	fSurfDepCoag = FALSE;
	// TNewton
	fWriteBT = FALSE;
	fWriteResiduum = FALSE;
	fWatchGridding = FALSE;
	fWriteEverySolution = FALSE;
	fOutputPath = "./";

	fInitialEquations = 0;
	fMaxGridPoints = 301;
	fInitialGridPoints = 101;
	fDampFlag = FALSE;
	fTimeDepFlag = FALSE;
	fContinFlag = FALSE;
	fDeltaNewGrid = 3;
	fTolRes = -1.0;
    fTolDy = -1.0;
    fMaxIter = 100;
    fLeft = 0.0;
    fRight = 0.0;
    fGamma = 0;
    fKappa = 1.0;
    fTauGrid = 1.0e-4;

    fCl = 0.0;
    fHv = 0.0;
    fT_B = 0.0;

	// TAdaptiveGrid
	fOneSolOneGrid = FALSE;
    fR = 60.0;
    fQ = 0.35;
    fStart = 0.0;
    fEnd = 0.0;
    fAlpha = 0.0;
	fAdjustComputationalDomain = TRUE;

	// TDamp
    fLambdaMin = 0.01;

	// TTime
    fDeltaTStart = -1.0;
    fDeltaTMax = 0.0;

	// TContinuation
    fContSteps = 10;

	yylex();
	
	fclose (fpA);
/*	fprintf( stderr, "first sweep finished\n" );*/
/*	fprintf( stderr, "now use input file %s\n", fAdditionalInfoFile );*/
/*	fpA = NULL;*/
/*	if ( ( fpA = fopen( fAdditionalInfoFile, "r" ) ) == NULL ){*/
/*		fprintf( stderr, "#error: can't open specified input file '%s'\n", fAdditionalInfoFile);*/
/*		exit( 2 );*/
/*	}*/
/*	else {*/
/*		yyin = fpA;*/
/*	}*/
/*	yylex();*/

	if ( fPressureComm > 0.0 ) {
		fPressure[0] = fPressureComm;
		fNPressures = 1;
	}
	
/*	char	*fileName = GetFullPath( fReactionFile, kFileName );
	if ( !fReactionFile || ( fpR = fopen( fileName, "r" ) ) == NULL ) {
		fprintf( stderr, "#error: can't open input file for reactions '%s'\n", fReactionFile );
		exit( 2 );
	}
	delete fileName;
	
*/
	char	*fileName;
	
	if ( fFlameType != kStartUpDiffusion && fFlameType != kHomoIsoChor && fFlameType != kHomoIsoBar && fFlameType != kTransFlamelet && fFlameType != kTrans1DIsoChor && fFlameType != kTransCoalParticle ) {
//	if ( fFlameType != kStartUpDiffusion && fFlameType != kHomoIsoChor && fFlameType != kHomoIsoBar && fFlameType != kTrans1DIsoChor ) {
		if ( !fStartProfileFileComm ) {
			fileName = GetFullPath( fStartProfileFile, kFileName );
			if ( !fStartProfileFile || ( fpS = fopen( fStartProfileFile, "r" ) ) == NULL ) {
				fprintf( stderr, "#error: can't open input file for startprofiles '%s'\n", fStartProfileFile );
				exit( 2 );
			}
			else {
//				fprintf( stderr, "use startprofiles file '%s'\n", fStartProfileFile );
			}
			delete fileName;
		}
		else {
			if ( ( fpS = fopen( fStartProfileFileComm, "r" ) ) == NULL ) {
				fprintf( stderr, "#error: can't open input file for startprofiles '%s'\n", fStartProfileFileComm );
				exit( 2 );
			}
			else {
				fStartProfileFile = fStartProfileFileComm;
//				fprintf( stderr, "use startprofiles file '%s'\n", fStartProfileFileComm );
			}
		}
	}

	if ( fFlameType == kStartUpDiffusion ) {
		fVariablesWithoutSpecies = 4;
	}
	else if ( fFlameType == kDiffusionPhys ) {
		fVariablesWithoutSpecies = 4;
	}
	else if ( fFlameType == kDiffPhysEigen ) {
		fVariablesWithoutSpecies = 5;
	}
	else if ( fFlameType == kStartUpDiffusion2 ) {
		fVariablesWithoutSpecies = 4;
	}
	else if ( fFlameType == kCountPremPhys ) {
		fVariablesWithoutSpecies = 3;
	}
	else if ( fFlameType == kCountPremSim ) {
		fVariablesWithoutSpecies = 4;
	}
	else if ( fFlameType == kCountDiffMix ) {
		fVariablesWithoutSpecies = 2;
	}
	else if ( fFlameType == kCountDiffCont ) {
		fVariablesWithoutSpecies = 5;
	}
	else if ( fFlameType == kUnstrPremPhys ) {
		fVariablesWithoutSpecies = 2;
	}
	else if ( fFlameType == kHomoIsoChor ) {
		fVariablesWithoutSpecies = 1;
	}
	else if ( fFlameType == kHomoIsoBar ) {
		fVariablesWithoutSpecies = 1;
	}
	else if ( fFlameType == kTransFlamelet ) {
		fVariablesWithoutSpecies = 1;
	}
	else if ( fFlameType == kTrans1DIsoChor ) {
		fVariablesWithoutSpecies = 4;
	}
	else if ( fFlameType == kTransCoalParticle ) {
		fVariablesWithoutSpecies = 4;
	}
	else{
		cerr << "error: no flametype specified" << NEWL;
		exit(2);
	}

	if ( fWithSoot ) {
		fVariablesWithoutSpecies += fNSootMoments;
	}
	else {
		fNSootMoments = 0;
	}
}

FirstInput::~FirstInput( void )
{
	delete leftBoundary;
	delete rightBoundary;
}

void FirstInput::PrintAdditionalData( void )
{
	int	i;
	FILE	*fp = NULL;


	fp = fopen( "FirstInput.tout", "w" );
	if ( !( fp = fopen( "FirstInput.tout", "w" ) ) ) { 
		cerr << "#warning: unable to open file 'FirstInput.tout'" << NEWL;
		return;
	}
	
	fprintf( fp, "FlameType is : " );
	if ( fFlameType == kStartUpDiffusion ) {
		fprintf( fp, "StartUpDiffusion\n" );
	}
	else if ( fFlameType == kDiffusionPhys ){
		fprintf( fp, "Diffusion\n" );
	}
	else if ( fFlameType == kStartUpDiffusion2 ){
		fprintf( fp, "second stage of StartUpDiffusion\n" );
	}

	fprintf( fp, "Left Boundary:\n" );

	leftBoundary->PrintBoundary( fp );
	fprintf( fp, "Right Boundary:\n" );
	rightBoundary->PrintBoundary( fp );

	fprintf( fp, "fInitialEquations = %d\n", fInitialEquations );
	fprintf( fp, "fVariablesWithoutSpecies = %d\n", fVariablesWithoutSpecies );
	
	fprintf( fp, "%d different StrainRates defined\n", nStrainRates );
	for ( i = 0; i < nStrainRates; ++i ) {
		fprintf( fp, "\t%g", fStrainRate[i] );
	}
	fprintf( fp, "\n" );
	
	fprintf( fp, "%d different scalar dissipation rates defined\n", nDissRates );
	for ( i = 0; i < nDissRates; ++i ) {
		fprintf( fp, "\t%g", fDissRate[i] );
	}
	fprintf( fp, "\n" );
	
	fprintf( fp, "fPressure = %g\n", fPressure[0] );
	fprintf( fp, "fMaxGridPoints = %d\n", fMaxGridPoints );
	fprintf( fp, "fInitialGridPoints = %d\n", fInitialGridPoints );
	fprintf( fp, "fDampFlag = %d\n", fDampFlag );
	fprintf( fp, "fTimeDepFlag = %d\n", fTimeDepFlag );
	fprintf( fp, "fContinFlag = %d\n", fContinFlag );
	fprintf( fp, "fDeltaNewGrid = %d\n", fDeltaNewGrid );
	fprintf( fp, "fTolRes = %g\n", fTolRes );
	fprintf( fp, "fTolDy = %g\n", fTolDy );
	fprintf( fp, "fMaxIter = %d\n", fMaxIter );
	fprintf( fp, "fLeft = %g\n", fLeft );
	fprintf( fp, "fRight = %g\n", fRight );
	fprintf( fp, "fR = %g\n", fR );
	fprintf( fp, "fQ = %g\n", fQ );
	fprintf( fp, "fContSteps = %d\n", fContSteps );
	fprintf( fp, "ReactionInputFile is %s\n", (fReactionFile[0])?fReactionFile:"not specified" );
	fprintf( fp, "OutputFile is %s\n", (fOutFileName[0])?fOutFileName:"not specified" );
	fprintf( fp, "GlobalReaction is %s\n", fGlobalReaction );
	fclose( fp );
}

void FirstInput::ParseCommandLine( int argc, char *argv[] )
{
	int 	error = FALSE;
		
   	while ( --argc > 0 && (*++argv)[0] == '-' ) {

		if ( *(argv[0]+1) == 'i' ) {
		    strcpy( fAdditionalInfoFile, *++argv ); --argc;
		    //WAS:	    strcpy( fAdditionalInfoFile, *++argv ); --argc;
//		    fAdditionalInfoFile = *++argv; --argc;
			if ( ( fpA = fopen( *argv, "r" ) ) == NULL ){
				fprintf( stderr, "#error: can't open specified input file '%s'\n", *argv);
				exit( 2 );
			}
			else {
				yyin = fpA;
			}
		}
		else if ( *(argv[0]+1) == 's' ) {
		   	fStartProfileFileComm = *++argv; --argc;
		}
		else if ( *(argv[0]+1) == 'r' ) {
		   	fMechanismFileComm = *++argv; --argc;
		}
		else if ( *(argv[0]+1) == 'o' ) {
		   	fOutputPathComm = new char[strlen(*++argv)+1];
			strcpy( fOutputPathComm, *argv );
			--argc;
		}
		else if ( *(argv[0]+1) == 'h' ) {
			error = TRUE;
		}
		else if ( *(argv[0]+1) == 'P' ) {
			fParameterComm = atof( *++argv ); --argc;
		}
		else if ( *(argv[0]+1) == 'd' ) {
			fScannerProgress = TRUE;
		}
		else if ( *(argv[0]+1) == 'p' ) {
			fPressureComm = atof( *++argv ) * 1.0e5; --argc;
		}
		else if ( *(argv[0]+1) == 'C' ) {
			fCoagFact = atof( *++argv ); --argc;
		}
		else {
			error = TRUE;
		}
	}
		
	if ( argc != 0 ) {
		error = TRUE;
	}
	
	if ( error ) {
		fprintf( stderr, "try: FlameMaster -i <arg> -s <arg> -o <arg> -r <arg> -p <arg> -P <arg> -d -h\n" );
		fprintf( stderr, "            -i:	following input file, default is 'FlameMaster.input'\n" );
		fprintf( stderr, "            -s:	following startprofiles file\n" );
		fprintf( stderr, "            -o:	following output path\n" );
		fprintf( stderr, "            -r:	following mechanism file\n" );
		fprintf( stderr, "            -p:	following real value for pressure p\n" );
		fprintf( stderr, "            -P:	following real value for parameter p\n" );
		fprintf( stderr, "            -d:	debug scanner\n" );
		fprintf( stderr, "            -h:	print online help\n" );
		exit(2);
	}

	if ( !fpA ) {
		if ( ( fpA = fopen( DEFAULTINPUTFILE, "r" ) ) == 0 ) {
			fprintf( stderr, "#error: can't open specified input file '%s'\n", DEFAULTINPUTFILE);
			exit( 2 );
		}
		else {
//			char *ingile;
//			ingile = new char( strlen(DEFAULTINPUTFILE) + 1 );
//		   	strcpy( ingile, DEFAULTINPUTFILE );
//			fAdditionalInfoFile = ingile;
//			int lenn = strlen(ingile);
//			fAdditionalInfoFile = new char( lenn + 1 );
		   	strcpy( fAdditionalInfoFile, DEFAULTINPUTFILE );
			cerr << "use input file '" << DEFAULTINPUTFILE << "'"  << NEWL;
		}
	}
	yyin = fpA;
}
   
void TInputData::InitInputData( FirstInputPtr firstInp )
{
	int	i;
	// the storage of the following structs is allocated in the ReadFunctions
	fCounter = NULL;
	fAtoms = NULL;
	fSpecies = NULL;
	fReaction = NULL;
	fThirdBody = NULL;
	fDimension = NULL;
	fHeader = NULL;
	fBuffer = NULL;
	fpR = firstInp->fpR;
//	fpO = firstInp->fpO;

	if ( !firstInp->fMechanismFileComm ) {
		if ( firstInp->fReactionFile ) {
			fReactionFile = new char[strlen( firstInp->fReactionFile ) + 1 ];
			strcpy( fReactionFile, firstInp->fReactionFile );
		}
		else {
			fReactionFile = NULL;
		}
	}
	else {
		fReactionFile = new char[strlen( firstInp->fMechanismFileComm ) + 1 ];
		strcpy( fReactionFile, firstInp->fMechanismFileComm );
	}

	if ( firstInp->fOutFileName ) {
		fOutFileName = new char[strlen( firstInp->fOutFileName ) + 1 ];
		strcpy( fOutFileName, firstInp->fOutFileName );
	}
	else {
		fOutFileName = NULL;
	}

	ReadInReactions();

	if ( firstInp->fGlobalReaction ) {
		fGlobalReaction = new char[strlen( firstInp->fGlobalReaction ) + 1];
		if ( !fGlobalReaction ) FatalError( "memory allocation in InitInputData failed" );
		strcpy( fGlobalReaction, firstInp->fGlobalReaction );
	}
	else {
		fGlobalReaction = NULL;
	}
	fInitialCond = new BoundaryInput( firstInp, firstInp->fInitialCond->fSpecifiedSpeciesBCs );
   	if ( !fInitialCond ) FatalError( "memory allocation of BoundaryInput failed" );
	leftBoundary = new BoundaryInput( firstInp, firstInp->leftBoundary->fSpecifiedSpeciesBCs );
   	if ( !leftBoundary ) FatalError( "memory allocation of BoundaryInput failed" );
	rightBoundary = new BoundaryInput( firstInp, firstInp->rightBoundary->fSpecifiedSpeciesBCs );
   	if ( !rightBoundary ) FatalError( "memory allocation of BoundaryInput failed" );

	*fInitialCond = *firstInp->fInitialCond;
	*leftBoundary = *firstInp->leftBoundary;
	*rightBoundary = *firstInp->rightBoundary;

	// Exact Backward Reaction Rates
	fExactBackward = firstInp->fExactBackward;
	
	// 0D-stuff
	fAdditionalOutput = firstInp->fAdditionalOutput;
	fEquidistant = firstInp->fEquidistant;
	fNOutputs = firstInp->fNOutputs;
	fArtificialSource = firstInp->fArtificialSource;
	fArtSourceTime = firstInp->fArtSourceTime;

	// set values of
	// T1DFlame
	fDateCreated = firstInp->fDateCreated;
	fAuthor = firstInp->fAuthor;
	fFlameType = firstInp->fFlameType;

	if ( firstInp->fContinType ) {
		if ( strcmp( firstInp->fContinType, "TEMPERATURE" ) == 0 ) {
			fContinType = kTemperature;
		}
		else if ( strcmp( firstInp->fContinType, "VELOCITY" ) == 0 ){
			fContinType = kVelocity;
		}
		else if ( strcmp( firstInp->fContinType, "MOMENTUM" ) == 0 ){
			fContinType = kMomentum;
		}
		else if ( strcmp( firstInp->fContinType, "SESHATEMP" ) == 0 ){
			fContinType = kSeshaTemp;
		}
		else if ( strcmp( firstInp->fContinType, "POLYTROPE" ) == 0 ){
			fContinType = kPolytrope;
		}
		else if ( strcmp( firstInp->fContinType, "PRESSURE" ) == 0 ){
			fContinType = kPressure;
		}
	}
	else {
		fContinType = kNoCont;
	}

	if ( firstInp->fContinSide ) {
		if ( strcmp( firstInp->fContinSide, "LEFT" ) == 0 ) {
			fContinSide = kLeft;
		}
		else if ( strcmp( firstInp->fContinSide, "RIGHT" ) == 0 ){
			fContinSide = kRight;
		}
		else if ( strcmp( firstInp->fContinSide, "BOTHSIDES" ) == 0 ){
			fContinSide = kBothSides;
		}
		else if ( strcmp( firstInp->fContinSide, "" ) != 0 ) {
			cerr << "#error: " << firstInp->fContinSide 
					<< " is not of type ContinuationSide" << NEWL;
			exit( 2 );
		}
	}
	else {
		fContinSide = kNoSide;
	}

	fContBound = firstInp->fContBound;

	fConstantLewisNumber = firstInp->fConstantLewisNumber;
	fWithRadiation = firstInp->fWithRadiation;
	fArcLengthContin = firstInp->fArcLengthContin;
	fStrainRateContin = firstInp->fStrainRateContin;

	fThermoDiffusion = firstInp->fThermoDiffusion;
	fWithSoot = firstInp->fWithSoot;
	fNSootMoments = firstInp->fNSootMoments;
	if ( fWithSoot ) {
		fNucleation = firstInp->fNucleation;
		fCondensation = firstInp->fCondensation;
		fCoagulation = firstInp->fCoagulation;
		fSurfaceGrowth = firstInp->fSurfaceGrowth;
		fSurfaceOxidation = firstInp->fSurfaceOxidation;
		fThermoPhoresis = firstInp->fThermoPhoresis;
		fOHPAHOxidation = firstInp->fOHPAHOxidation;
		fO2PAHOxidation = firstInp->fO2PAHOxidation;
		fCoagFact = firstInp->fCoagFact;
		fSootRadiation = firstInp->fSootRadiation;
		fSootUpdateProdRates = firstInp->fSootUpdateProdRates;
		fSizeDepDiff = firstInp->fSizeDepDiff;
		fSurfDepCoag = firstInp->fSurfDepCoag;
	}

	fCompUnPhysChain = firstInp->fCompUnPhysChain;
	fNUnPhysChain = firstInp->fNUnPhysChain;
	fIsAxiSymmetric = firstInp->fIsAxiSymmetric;
	fClipNegativeConcs = firstInp->fClipNegativeConcs;
	fNoDiffCorr = firstInp->fNoDiffCorr;
	fWriteFullRes = firstInp->fWriteFullRes;
	fPrintMolarFractions = firstInp->fPrintMolarFractions;
	fSteadyStatesNumerical = firstInp->fSteadyStatesNumerical;
	fUseModifiedNewton = firstInp->fUseModifiedNewton;
	fUseNumericalJac = firstInp->fUseNumericalJac;
	fUseSecOrdJac = firstInp->fUseSecOrdJac;
	fUseNumericalDM = firstInp->fUseNumericalDM;
	fMinStrainRate = firstInp->fMinStrainRate;
	fVariablesWithoutSpecies = firstInp->fVariablesWithoutSpecies;

	if ( !firstInp->nStrainRates || ( fFlameType == kCountDiffCont && firstInp->nStrainRates != 1 ) ) {
		firstInp->nStrainRates = 0;
		fStrainRate = NULL;
	}
	else {
		fStrainRate = NewVector( firstInp->nStrainRates );
	}
	for ( i = 0; i < firstInp->nStrainRates; ++i ) {
		fStrainRate->vec[i] = firstInp->fStrainRate[i];
	}

	if ( !firstInp->fNPressures ) {
		fPressure = NULL;
	}
	else {
		fPressure = NewVector( firstInp->fNPressures );
	}
	for ( i = 0; i < firstInp->fNPressures; ++i ) {
		fPressure->vec[i] = firstInp->fPressure[i];
	}

	if ( firstInp->nDissRates ) {
		fDissRate = NewVector( firstInp->nDissRates );
	}
	else {
		fDissRate = NULL;
	}
	for ( i = 0; i < firstInp->nDissRates; ++i ) {
		fDissRate->vec[i] = firstInp->fDissRate[i];
	}

	fParameterComm = firstInp->fParameterComm;
	fPremConfiguration = firstInp->fPremConfiguration;
	
	if ( !firstInp->fNPhi ) {
		fPhi = NULL;
	}
	else {
		fPhi = NewVector( firstInp->fNPhi );
	}
	for ( i = 0; i < firstInp->fNPhi; ++i ) {
		fPhi->vec[i] = firstInp->fPhi[i];
	}

	fDeltaSCont = firstInp->fDeltaSCont;
//	fChi = firstInp->fChi;

	fFuelIndex = NewIntVector( firstInp->fNFuels );
	for ( i = 0; i < firstInp->fNFuels; ++i ) {
		if ( ( fFuelIndex->vec[i] = FindSpecies( firstInp->fFuelName[i] ) ) == -1 ) {
			cerr << "error: can't find specified fuel: " << firstInp->fFuelName[i] << NEWL;
			exit( 2 );
		}
	}
	if ( ( fOxIndex = FindSpecies( firstInp->fOxName ) ) == -1 ) {
		cerr << "error: can't find specified oxidizer: " << firstInp->fOxName << NEWL;
		exit( 2 );
	}
//	if ( fFlameType == kStartUpDiffusion || fFlameType == kStartUpDiffusion2 
//			|| fFlameType == kDiffusionPhys || fFlameType == kDiffPhysEigen || fFlameType == kCountDiffCont ) {
		if ( ( fH2OIndex = FindSpecies( "H2O" ) ) == -1 ) {
//			cerr << "error: can't find species 'H2O'" << NEWL;
//			exit( 2 );
		}
		if ( ( fCO2Index = FindSpecies( "CO2" ) ) == -1 ) {
//			cerr << "warning: can't find species 'CO2', check if the code is able to handle this case" << NEWL;
		}
//	}	
	
	if ( firstInp->fFromSpeciesName ) {
		if ( ( fFromSpeciesIndex = FindSpecies( firstInp->fFromSpeciesName ) ) == -1 ) {
			cerr << "warning: can't find specified 'FromSpecies': " << firstInp->fFromSpeciesName << NEWL;
			fFromSpeciesIndex = -1;
		}
	}
	else {
		fFromSpeciesIndex = -1;
	}
	
	if ( firstInp->fToSpeciesName ) {
		if ( ( fToSpeciesIndex = FindSpecies( firstInp->fToSpeciesName ) ) == -1 ) {
			cerr << "warning: can't find specified 'ToSpecies': " << firstInp->fToSpeciesName << NEWL;
			fToSpeciesIndex = -1;
		}
	}
	else {
		fToSpeciesIndex = -1;
	}
	

	// copy objects for sensitivity analysis
	fSensAnal = firstInp->fSensAnal;
	fSensAnalSpec = firstInp->fSensAnalSpec;
	fSensObjAll = firstInp->fSensObjAll;
	fSensAnalReac = firstInp->fSensAnalReac;
	fSensAnalFac = firstInp->fSensAnalFac;
	fSensMax = firstInp->fSensMax;
	fSensFinal = firstInp->fSensFinal;
	fNSensObj = firstInp->fNSensObj;
	fSensObj = new char* [fNSensObj];
	CopyStringArray( fSensObj, firstInp->fSensObj, fNSensObj );
	fReactionFluxes = firstInp->fReactionFluxes;
	fPrintRHSSpecies = firstInp->fPrintRHSSpecies;
	fPrintRHSTemp = firstInp->fPrintRHSTemp;
	fKeepMassFracs = firstInp->fKeepMassFracs;
	fLiquidPoolBC = firstInp->fLiquidPoolBC;
	
	fContInc = firstInp->fContInc;	

	// TNewton
	fWriteBT = firstInp->fWriteBT;
	fWriteResiduum = firstInp->fWriteResiduum;
	fWatchGridding = firstInp->fWatchGridding;
	fWriteEverySolution = firstInp->fWriteEverySolution;
	if ( firstInp->fOutputPathComm ) {
		fOutputPath = firstInp->fOutputPathComm;
	}
	else {
		fOutputPath = firstInp->fOutputPath;
	}

	fNVariables = fCounter->species - fCounter->steadyStates + fVariablesWithoutSpecies;

	if ( firstInp->fInitialEquations != 0 ) {
		if ( firstInp->fInitialEquations < 0 ) {
			fInitialEquations = fNVariables + firstInp->fInitialEquations;
		}
		else {
			fInitialEquations = ( int )MIN( fNVariables, firstInp->fInitialEquations );
		}
	}
	else {
		fInitialEquations = fNVariables;
	}
	fMaxGridPoints = ( firstInp->fMaxGridPoints % 2 ) ? firstInp->fMaxGridPoints : firstInp->fMaxGridPoints + 1;
	fInitialGridPoints = firstInp->fInitialGridPoints;
	fDampFlag = firstInp->fDampFlag;
	fTimeDepFlag = firstInp->fTimeDepFlag;
	fContinFlag = firstInp->fContinFlag;
	fDeltaNewGrid = firstInp->fDeltaNewGrid;
	fTolRes = firstInp->fTolRes;
    fTolDy = firstInp->fTolDy;
    fMaxIter = firstInp->fMaxIter;
    fLeft = firstInp->fLeft;
    fRight = firstInp->fRight;
    fGamma = firstInp->fGamma;
    fKappa = firstInp->fKappa;
    fTauGrid = firstInp->fTauGrid;

    fCl = firstInp->fCl;
    fHv = firstInp->fHv;
    fT_B = firstInp->fT_B;

	// TAdaptiveGrid
	fOneSolOneGrid = firstInp->fOneSolOneGrid;
    fR = firstInp->fR;
    fQ = firstInp->fQ;
    fStart = firstInp->fStart;
    fEnd = firstInp->fEnd;
    fAlpha = firstInp->fAlpha;
	fAdjustComputationalDomain = firstInp->fAdjustComputationalDomain;

	// TDamp
    fLambdaMin = firstInp->fLambdaMin;

	// TTime
    fDeltaTStart = firstInp->fDeltaTStart;
    fDeltaTMax = firstInp->fDeltaTMax;

	// TContinuation
    fContSteps = firstInp->fContSteps;
	
	
	// Files
	fAddFileNo1 = firstInp->fAddFileNo1;
	fLewisNumberFile = firstInp->fLewisNumberFile;
	if ( firstInp->fStartProfileFile && strlen( firstInp->fStartProfileFile ) ) {
		fStartProfileFile = new char[strlen( firstInp->fStartProfileFile )+1];
		strcpy( fStartProfileFile, firstInp->fStartProfileFile );
	}
	else {
		fStartProfileFile = NULL;
	}
	
	fpR = firstInp->fpR;		// filepointer to fReactionFile
	fpS = firstInp->fpS;		// filepointer to fStartProfileFile
//	fpO = firstInp->fpO;		// filepointer to fReactionFile

//	PrintAdditionalData();
}

TInputData::~TInputData( void )
{
	delete fStartProfileFile;
	delete fOutFileName;
	delete fReactionFile;

	if ( fPhi ) {
		DisposeVector( fPhi );
	}
	if ( fDissRate ) {
		DisposeVector( fDissRate );
	}
	if ( fPressure ) {
		DisposeVector( fPressure );
	}
	if ( fStrainRate ) {
		DisposeVector( fStrainRate );
	}

	delete rightBoundary;
	delete leftBoundary;
	delete fGlobalReaction;

	FreeDimensionArray( fDimension, fCounter->dimensions );
	FreeThirdBodyArray( fThirdBody, fCounter->thirdBodies );
	FreeReactionArray( fReaction, fCounter->reactions );
	FreeSpeciesArray( fSpecies, fCounter->species );
	FreeAtomsArray( fAtoms, fCounter->atoms );
	delete fCounter;
	delete fBuffer;
	
	// sensitivity stuff
	if ( !fSensObjAll ){
	  for ( int i = 0; i < fNSensObj; ++i ) {
	    delete fSensObj[i];
	  }
	  delete fSensObj;
	}
	delete fHeader;
}

void TInputData::PrintAdditionalData( void )
{
	FILE	*fp = NULL;
	char	*outPath = GetFullPath( fOutputPath, kPath );
	char	*fileName = new char[ strlen(outPath) + 15 ];

	sprintf( fileName, "%sinputData.tout", outPath );
	if ( !( fp = fopen( fileName, "w" ) ) ) {
		cerr << "#error: unable to open file '" << fileName << "'" << NEWL; 
		exit(2);
	}
	
	delete fileName;
	delete outPath;
	
	fprintf( fp, "FlameType is : " );
	if ( fFlameType == kStartUpDiffusion ) {
		fprintf( fp, "StartUpDiffusion\n" );
	}
	else if ( fFlameType == kDiffusionPhys ){
		fprintf( fp, "Diffusion in physical coordinates\n" );
	}

	fprintf( fp, "\n" );
	fprintf( fp, "Left Boundary:\n" );
	leftBoundary->PrintBoundary( fp );
	fprintf( fp, "\n" );
	fprintf( fp, "Right Boundary:\n" );
	rightBoundary->PrintBoundary( fp );
	
	
	if ( fStrainRate ) {
		fprintf( fp, "%d different StrainRates defined\n", fStrainRate->len );
		for ( int i = 0; i < fStrainRate->len; ++i ) {
			fprintf( fp, "\t%g", fStrainRate->vec[i] );
		}
	}
	fprintf( fp, "\n" );

	if ( fDissRate ) {
		fprintf( fp, "%d different scalar dissipation rates defined\n", fDissRate->len );
		for ( int i = 0; i < fDissRate->len; ++i ) {
			fprintf( fp, "\t%g", fDissRate->vec[i] );
		}
	}
	fprintf( fp, "\n" );

	fprintf( fp, "%d different pressures defined\n", fPressure->len );
	for ( int i = 0; i < fPressure->len; ++i ) {
		fprintf( fp, "\t%g", fPressure->vec[i] );
	}
	fprintf( fp, "\n" );
	
	if ( fFlameType == kStartUpDiffusion || fFlameType == kStartUpDiffusion2 ) {
		for ( int i = 0; i < fFuelIndex->len; ++i ) {
			fprintf( fp, "fFuelIndex[%d] = %d\n", i, fFuelIndex->vec[i] );
		}
		fprintf( fp, "fOxIndex = %d\n", fOxIndex );
		fprintf( fp, "fH2OIndex = %d\n", fH2OIndex );
		fprintf( fp, "fCO2Index = %d\n", fCO2Index );
	}
	
	fprintf( fp, "fVariablesWithoutSpecies = %d\n", fVariablesWithoutSpecies );
	fprintf( fp, "fInitialEquations = %d\n", fInitialEquations );
	fprintf( fp, "fMaxGridPoints = %d\n", fMaxGridPoints );
	fprintf( fp, "fInitialGridPoints = %d\n", fInitialGridPoints );
	fprintf( fp, "fDampFlag = %d\n", fDampFlag );
	fprintf( fp, "fTimeDepFlag = %d\n", fTimeDepFlag );
	fprintf( fp, "fContinFlag = %d\n", fContinFlag );
	fprintf( fp, "fDeltaNewGrid = %d\n", fDeltaNewGrid );
	fprintf( fp, "fTolRes = %g\n", fTolRes );
	fprintf( fp, "fTolDy = %g\n", fTolDy );
	fprintf( fp, "fMaxIter = %d\n", fMaxIter );
	fprintf( fp, "fLeft = %g\n", fLeft );
	fprintf( fp, "fRight = %g\n", fRight );
	fprintf( fp, "fR = %g\n", fR );
	fprintf( fp, "fQ = %g\n", fQ );
	fprintf( fp, "fLambdaMin = %g\n", fLambdaMin );
	fprintf( fp, "fContSteps = %d\n", fContSteps );
	fclose( fp );
}


//// ReadInAdditionalData
//

void TInputData::ReadInAdditionalData( void )
{
  // EMPTY //
}

void TInputData::ReadInReactions( void )
{
/* open input file */
	char	*fileName = GetFullPath( fReactionFile, kFileName );
	if ( !fReactionFile ) {
		fprintf( stderr, "#error: no reaction file specified\n" );
		exit( 2 );
	}
	if ( ( fpR = fopen( fileName, "rb" ) ) == NULL ) {
		if ( ( fpR = fopen( MyDataFile( fileName ), "r" ) ) == NULL ) {
		fprintf( stderr, "#error: can't open input file for reactions '%s'\n", fReactionFile );
		exit( 2 );
	}
	else {
		fprintf( stderr, "use mechanism file '%s'\n", fileName );
	}
	}
	else {
		fprintf( stderr, "use mechanism file '%s'\n", fileName );
	}
	delete fileName;
	
/* header */
	fHeader = ReadHeader( );

	if ( strcmp( fHeader->version, (char *)VERSION ) == 0 ) {
		fCounter = ReadCounter();
		fAtoms = ReadAtoms();
		fSpecies = ReadSpecies();
		fReaction = ReadReaction( fCounter->reactions );

		if ( fread( fBuffer, fHeader->maxLenOfString, 1, fpR ) != 1 ) {
			FatalError( "invalid PAHSymbol in mechanism input file\n" );
		}
		fPAHSymbol = ( String )malloc( strlen( fBuffer ) + 1 );
		if ( !fPAHSymbol ) FatalError( "malloc failed\n" );
		strcpy( fPAHSymbol, fBuffer );
		fPAHReaction = ReadReaction( fCounter->pahReactions );

		if ( fread( fBuffer, fHeader->maxLenOfString, 1, fpR ) != 1 ) {
			FatalError( "invalid SootSymbol in mechanism input file\n" );
		}
		fSootSymbol = ( String )malloc( strlen( fBuffer ) + 1 );
		if ( !fSootSymbol ) FatalError( "malloc failed\n" );
		strcpy( fSootSymbol, fBuffer );
		fSootReaction = ReadReaction( fCounter->sootReactions );
		
		fThirdBody = ReadThirdBody();
		fDimension = ReadDimension();
	}
	else {
		cerr << "##error: input file for the mechanism '" << fReactionFile 
							<< "'\n\thas been generated with ScanMan Version "
							<< fHeader->version << "\n\t-> synchronize versions\n";
		fclose( fpR );
		exit( 2 );
	}
/* close input file */
	fclose( fpR );
}


HeaderPtr TInputData::NewHeader( void )
{	
	HeaderPtr header = NULL;
	
	if ( !( header = ( HeaderPtr )calloc( 1, sizeof( Header ) ) ) ) {
		fprintf( stderr, "# memory allocation of Header_struct failed\n");
		exit(2);
	}
	
	return header;
}

HeaderPtr TInputData::ReadHeader( void )
{
	HeaderPtr	header = NULL;
	
	header = NewHeader();
	
	if ( fread( header, sizeof( Header ), 1, fpR ) != 1 ) {
		FatalError( "invalid struct header in mechanism input file\n" );
	}
	fBuffer = ( String ) calloc( 1, header->maxLenOfString );

	return header;
}

CounterPtr TInputData::ReadCounter( void )
{
	CounterPtr	counter = new Counter;
   	if ( !counter ) FatalError( "memory allocation of Counter failed" );
	
	if ( fread( counter, sizeof( Counter ), 1, fpR ) != 1 ) {
		FatalError( "invalid struct counter in mechanism input file\n" );
	}

	return counter;
}

AtomsPtr TInputData::NewAtomArray( int len )
{
	AtomsPtr atom = NULL;
	
	if ( !( atom = ( AtomsPtr )calloc( 1, len * sizeof( Atoms ) ) ) ) {
		fprintf( stderr, "# memory allocation of Atoms_struct_array failed\n");
		exit(2);
	}
		
	return atom;
}

AtomsPtr TInputData::ReadAtoms( void )
{
	AtomsPtr		atoms = NULL;
	int 			i;
	
	atoms = NewAtomArray( fCounter->atoms );
	for ( i = 0; i < fCounter->atoms; ++i ) {
		if ( fread( &atoms[i], sizeof( Atoms ), 1, fpR ) != 1 ) {
			FatalError( "invalid struct atoms in mechanism input file\n" );
		}
		if ( fread( fBuffer, fHeader->maxLenOfString, 1, fpR ) != 1 ) {
			FatalError( "invalid struct atoms in mechanism input file\n" );
		}
		atoms[i].name = ( String )malloc( strlen( fBuffer ) + 1 );
		if ( !atoms[i].name ) FatalError( "malloc failed\n" );
		strcpy( atoms[i].name, fBuffer );
	}
	
	return atoms;
}

void TInputData::FreeAtomsArray( AtomsPtr atoms, int len )
{
	for ( int i = 0; i < len; ++i ) {
		free( fAtoms[i].name );
	}
	delete atoms;
}

SpeciesPtr TInputData::NewSpeciesArray( int len )
{
	SpeciesPtr species = NULL;
	
	if ( !( species = ( SpeciesPtr )calloc( 1, len * sizeof( Species ) ) ) ) {
		fprintf( stderr, "# memory allocation of Species_struct_array failed\n");
		exit(2);
	}
	
	return species;
}

SpeciesPtr TInputData::ReadSpecies( void )
{
	SpeciesPtr		species = NULL;
	int 			vectorLen;
	int 			i;
	
	species = NewSpeciesArray( fCounter->species );
	for ( i = 0; i < fCounter->species; ++i ) {
		if ( fread( &species[i], sizeof( Species ), 1, fpR ) != 1 ) {
			FatalError( "invalid struct species in mechanism input file\n" );
		}
		if ( fread( fBuffer, fHeader->maxLenOfString, 1, fpR ) != 1 ) {
			FatalError( "invalid struct species in mechanism input file\n" );
		}
		species[i].name = ( String )malloc( strlen( fBuffer )+1 );
		if ( !species[i].name ) FatalError( "malloc failed\n" );
		strcpy( species[i].name, fBuffer );
		
		if ( fread( &vectorLen, sizeof( int ), 1, fpR ) != 1 ) {
			FatalError( "invalid struct species in mechanism input file\n" );
		}
		species[i].composition = NewIntVector( vectorLen );
		if ( fread( species[i].composition->vec, sizeof( int ), vectorLen, fpR ) != vectorLen ) {
			FatalError( "invalid struct species in mechanism input file\n" );
		}
		
		if ( fread( &vectorLen, sizeof( int ), 1, fpR ) != 1 ) {
			FatalError( "invalid struct species in mechanism input file\n" );
		}
		species[i].dCoeff = NewVector( vectorLen );
		if ( fread( species[i].dCoeff->vec, sizeof( Double ), vectorLen, fpR ) != vectorLen ) {
			FatalError( "invalid struct species in mechanism input file\n" );
		}
		
		if ( fread( &vectorLen, sizeof( int ), 1, fpR ) != 1 ) {
			FatalError( "invalid struct species in mechanism input file\n" );
		}
		species[i].omegaCoeff = NewVector( vectorLen );
		if ( fread( species[i].omegaCoeff->vec, sizeof( Double ), vectorLen, fpR ) != vectorLen ) {
			FatalError( "invalid struct species in mechanism input file\n" );
		}
	}
	
	return species;
}

void TInputData::FreeSpeciesArray( SpeciesPtr species, int len )
{
	for ( int i = 0; i < len; ++i ) {
		DisposeVector( species[i].omegaCoeff );
		DisposeVector( species[i].dCoeff );
		DisposeIntVector( species[i].composition );
		free( species[i].name );
	}
	delete species;
}

ReactionPtr TInputData::NewReactionArray( int len )
{
	ReactionPtr reaction = NULL;
	
	if ( !len ) {
		return NULL;
	}
	if ( !( reaction = ( ReactionPtr )calloc( 1, len * sizeof( Reaction ) ) ) ) {
		fprintf( stderr, "# memory allocation of Reaction_struct_array failed\n");
		fprintf( stderr, "# tried to allocate %d ReactionStructs of size %d\n", len, sizeof( Reaction ) );
		exit(2);
	}
	
	return reaction;
}

ReactionPtr TInputData::ReadReaction( int len )
{
	ReactionPtr			reaction = NULL;
	int 				i;
	
	reaction = NewReactionArray( len );
	for ( i = 0; i < len; ++i ) {
		if ( fread( &reaction[i], sizeof( Reaction ), 1, fpR ) != 1 ) {
			FatalError( "invalid struct reaction in mechanism input file\n" );
		}
		if ( fread( fBuffer, fHeader->maxLenOfString, 1, fpR ) != 1 ) {
			FatalError( "invalid struct reaction in mechanism input file\n" );
		}
		reaction[i].label = ( String )malloc( strlen( fBuffer ) + 1 );
		if ( !reaction[i].label ) FatalError( "malloc failed\n" );
		strcpy( reaction[i].label, fBuffer );
	}
	
	return reaction;
}

void TInputData::FreeReactionArray( ReactionPtr reaction, int len )
{
	for ( int i = 0; i < len; ++i ) {
		free( reaction[i].label );
	}
	delete reaction;
}

ThirdBodyPtr TInputData::NewThirdBodyArray( int len )
{
	ThirdBodyPtr thirdBody = NULL;
	
	if ( len == 0 ) {
		return NULL;
	}
	if ( !( thirdBody = ( ThirdBodyPtr )calloc( 1, len * sizeof( ThirdBody ) ) ) ) {
		fprintf( stderr, "# memory allocation of ThirdBody_struct_array failed\n");
		cerr << "# malloc tried to allocate " << len << " * " << sizeof( ThirdBody )
		   << " bytes in the free store" << NEWL;
		exit(2);
	}
	
	return thirdBody;
}

ThirdBodyPtr TInputData::ReadThirdBody( void )
{
	ThirdBodyPtr		thirdBody = NULL;
	int 				arrayLen;
	int 				i;
	
	thirdBody = NewThirdBodyArray( fCounter->thirdBodies );
	
	for ( i = 0; i < fCounter->thirdBodies; ++i ) {
		if ( fread( &thirdBody[i], sizeof( ThirdBody ), 1, fpR ) != 1 ) {
			FatalError( "invalid struct thirdBody in mechanism input file\n" );
		}
		
		if ( fread( fBuffer, fHeader->maxLenOfString, 1, fpR ) != 1 ) {
			FatalError( "invalid struct thirdBody in mechanism input file\n" );
		}
		thirdBody[i].name = ( String )malloc( strlen( fBuffer ) + 1 );
		if ( !thirdBody[i].name ) FatalError( "malloc failed\n" );
		strcpy( thirdBody[i].name, fBuffer );
		
		if ( fread( &arrayLen, sizeof( int ), 1, fpR ) != 1 ) {
			FatalError( "invalid struct thirdBody in mechanism input file\n" );
		}
		thirdBody[i].speciesNumber = NewIntVector( arrayLen );
		if ( fread( thirdBody[i].speciesNumber->vec, sizeof( int ), arrayLen, fpR ) != arrayLen ) {
			FatalError( "invalid struct thirdBody in mechanism input file\n" );
		}

		if ( fread( &arrayLen, sizeof( int ), 1, fpR ) != 1 ) {
			FatalError( "invalid struct thirdBody in mechanism input file\n" );
		}
		thirdBody[i].speciesCoeff = NewVector( arrayLen );
		if ( fread( thirdBody[i].speciesCoeff->vec, sizeof( Double ), arrayLen, fpR ) != arrayLen ) {
			FatalError( "invalid struct thirdBody in mechanism input file\n" );
		}
	}
	
	return thirdBody;
}

void TInputData::FreeThirdBodyArray( ThirdBodyPtr thirdBody, int len )
{
	for ( int i = 0; i < len; ++i ) {
		DisposeVector( thirdBody[i].speciesCoeff );
		DisposeIntVector( thirdBody[i].speciesNumber );
		free( thirdBody[i].name );
	}
	delete thirdBody;
}

DimensionPtr TInputData::NewDimensionArray( int len )
{
	DimensionPtr dimension = NULL;
	
	if ( len == 0 ) {
		return NULL;
	}
	if ( !( dimension = ( DimensionPtr )calloc( 1, len * sizeof( Dimension ) ) ) ) {
		fprintf( stderr, "# memory allocation of Dimension_struct_array failed\n");
		exit(2);
	}
	
	return dimension;
}

DimensionPtr TInputData::ReadDimension( void )
{
	DimensionPtr		dimension = NULL;
	int 				i;
	
	dimension = NewDimensionArray( fCounter->dimensions );
	for ( i = 0; i < fCounter->dimensions; ++i ) {
		if ( fread( &dimension[i], sizeof( Dimension ), 1, fpR ) != 1 ) {
			FatalError( "invalid struct dimension in mechanism input file\n" );
		}
		if ( fread( fBuffer, fHeader->maxLenOfString, 1, fpR ) != 1 ) {
			FatalError( "invalid struct dimension in mechanism input file\n" );
		}
		dimension[i].name = ( String )malloc( strlen( fBuffer ) + 1 );
		if ( !dimension[i].name ) FatalError( "malloc failed\n" );
		strcpy( dimension[i].name, fBuffer );
	
		dimension[i].parameterValues = NULL;
	}
	
	return dimension;
}

void TInputData::FreeDimensionArray( DimensionPtr dimension, int len )
{
	for ( int i = 0; i < len; ++i ) {
		free( dimension[i].name );
	}
	delete dimension;
}

int TInputData::FindAtomIndex( ConstString atom )
{
	int	nOfAtoms = fCounter->atoms;
	
	for ( int i = 0; i < nOfAtoms; ++ i ) {
		if ( strcmp( atom, fAtoms[i].name ) == 0 ) {
			return i;
		}
	}
	return -1;
}

int TInputData::FindSpeciesCoefficient( int speciesIndex, int reactionIndex )
{
	int	nOfSpecies = fReaction[reactionIndex].numberOfSpecies;
	
	for ( int i = 0; i < nOfSpecies; ++ i ) {
		if ( fReaction[reactionIndex].speciesNumber[i] == speciesIndex ) {
			return (int)fReaction[reactionIndex].speciesCoeff[i];
		}
	}
	return -1;
}

int TInputData::FindSpecies( ConstString species )
{
	int	nOfSpecies = fCounter->species;
	
	for ( int i = 0; i < nOfSpecies; ++ i ) {
		if ( strcmp( species, fSpecies[i].name ) == 0 ) {
			return i;
		}
	}
	return -1;
}
