#include "FlameMaster.h"
#include "TTrans1DIsoChorSolver.h"
#include "dassl.h"
#include "MapMan.h"
#include "Interrupt.h"

#undef DEBUGRES
#undef DEBUGGRID
#undef DEBUGINITIAL
#undef DEBUGWEIGHTS
#undef DEBUGGPWEIGHTS
#define DEBUGYPRIME

#undef VECTORTOL
#undef DPDT
#define DELTAPROG 1

void TTrans1DIsoChorSolver::InitTTrans1DIsoChorSolver( void )
{
	int	i;

// set output file
	if ( !fInputData->fOutFileName || strcmp( fInputData->fOutFileName, "sterr" ) == 0 ) {
		fOutFilePtr = stderr;
		fprintf( stderr, "write output to 'stderr'\n" );
	}
	else if ( strcmp( fInputData->fOutFileName, "stdout" ) == 0 ) {
		fOutFilePtr = stdout;
		fprintf( stderr, "write output to 'stdout'\n" );
	}
	else {
		fOutFilePtr = GetOutfile( fInputData->fOutFileName, kText );
		fprintf( stderr, "write output to file '%s'\n\n"
			, GetOutfileName( fInputData->fOutFileName, kText ) );
	}

	fFirstCall = TRUE;
	fMaxOrd = 5;
	fEquidistant = fInputData->fEquidistant;
	fGamma = fInputData->fGamma;
	fKappa = fInputData->fKappa;
	fArtificialSource = fInputData->fArtificialSource;
	fTauGrid = fInputData->fTauGrid;
	fPrintMolarFractions = fInputData->fPrintMolarFractions;
	if ( fInputData->fDeltaTMax > 0.0 ) {
		fDeltaTMax = fInputData->fDeltaTMax;
	}
	else {
		fDeltaTMax = 1.0e-4;
	}

	fNGridPoints = fInputData->fInitialGridPoints-2;

	fFirstTimeStep = fInputData->fDeltaTStart;

	int	nOfSpecies = fSpecies->GetNOfSpecies();
	int	nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();

#ifdef VECTORTOL
	fRTol = New1DArray( fNOfEquations );
	fATol = New1DArray( fNOfEquations );
	for ( i = 0; i < fNOfEquations; ++i ) {
		fRTol[i] = ( fInputData->fTolRes > 0.0 ) ? fInputData->fTolRes : 1.0e-4;
		fATol[i] = ( fInputData->fTolDy > 0.0 ) ? fInputData->fTolDy : 1.0e-12;
	}
//	fRTol[GetFuelIndex()+fFirstSpecies] = 1.0e-6;
//	fATol[GetFuelIndex()+fFirstSpecies] = 1.0e-13;
#else
	fRTol = New1DArray( 1 );
	fATol = New1DArray( 1 );
	fRTol[0] = ( fInputData->fTolRes > 0.0 ) ? fInputData->fTolRes : 1.0e-4;
	fATol[0] = ( fInputData->fTolDy > 0.0 ) ? fInputData->fTolDy : 1.0e-12;
#endif
	
// workspace for ddassl
	fNOfSolver = 1;
	int ML = 2 * fNOfEquations;
	int	NEQ = fNGridPoints * fNOfEquations;
	fLRW =  40 + ( fMaxOrd + 4 ) * NEQ 
				+ ( 2 * ML + ML + 1 ) * NEQ + 2 * ( NEQ / ( ML + ML + 1 ) + 1 );
	fLIW = 20 + NEQ;

	fDasslNEq = NEQ;
	fTime = NewVector( fNOfSolver );
	fInfo = new IntVectorPtr[fNOfSolver];
	fRWrk = new VectorPtr[fNOfSolver];
	fIWrk = new IntVectorPtr[fNOfSolver];
	fNActualStep = new int*[fNOfSolver];
	fNActualOrd = new int*[fNOfSolver];
	fMaxStepTaken = New1DArray( fNOfSolver );
   	if ( !fRWrk || !fIWrk || !fInfo || !fNActualStep || !fNActualOrd ) {
		FatalError( "memory allocation of TTrans1DIsoChorSolver failed" );
	}
	for ( i = 0; i < fNOfSolver; ++i ) {
		fInfo[i] = NewIntVector( 16 );
		fRWrk[i] = NewVector( fLRW );
		fIWrk[i] = NewIntVector( fLIW );
		fNActualStep[i] = &fIWrk[i]->vec[kF11];
		fNActualOrd[i] = &fIWrk[i]->vec[kF8];
	}
	InitDassl( fNOfSolver );
	if ( fInfo[0]->vec[kF10] == 1 ) {
		fprintf( fOutFilePtr, "###attention: clip negative concentrations\n" );
	}
	fSolution = NewMatrix( fNOfEquations, fNGridPoints+2, kColumnPointers );
	fSolution->mat = &fSolution->mat[kNext];

	fSolPrime = NewMatrix( fNOfEquations, fNGridPoints, kColumnPointers );
	
// other workspace
	
	fSolTime = NewVector( fNGridPoints+2 );
	fSolMassFracs = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
	fSolTemp = NewVector( fNGridPoints+2 );
	fSolPressure = NewVector( fNGridPoints+2 );
	fSolR = NewVector( fNGridPoints+2 );
	fSolPsi = NewVector( fNGridPoints+2 );
	fSolOldTime = NewVector( fNGridPoints+2 );
	fSolOldMassFracs = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
	fSolOldTemp = NewVector( fNGridPoints+2 );
	fSolOldPressure = NewVector( fNGridPoints+2 );
	fSolOldR = NewVector( fNGridPoints+2 );
	fSolOldPsi = NewVector( fNGridPoints+2 );
	fh = New1DArray( fNGridPoints );
	fhm = New1DArray( fNGridPoints );
	fhnenn = New1DArray( fNGridPoints );
	fFDWCurr = New1DArray( fNGridPoints );
	fFDWPlus = New1DArray( fNGridPoints );
	fFDWMinus = New1DArray( fNGridPoints );
	fWPlus = New1DArray( fNGridPoints );
	fWMinus = New1DArray( fNGridPoints );
	fGPDistWeight = New1DArray( fNGridPoints );
	fGPDens = NewVector( fNGridPoints+2 );
	fMassFracsWork = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
	fTempWork = NewVector( fNGridPoints+2 );
	fPressureWork = NewVector( fNGridPoints+2 );
	fRWork = NewVector( fNGridPoints+2 );
	fPsiWork = NewVector( fNGridPoints+2 );
	fOutSolWork = NewVector( fNGridPoints+2 );
	fRho = New1DArray( fNGridPoints );
	fMixHeatCap = New1DArray( fNGridPoints );

	fSolTime->vec = &fSolTime->vec[kNext];
	fSolMassFracs->mat = &fSolMassFracs->mat[kNext];
	fSolTemp->vec = &fSolTemp->vec[kNext];
	fSolPressure->vec = &fSolPressure->vec[kNext];
	fSolR->vec = &fSolR->vec[kNext];
	fSolPsi->vec = &fSolPsi->vec[kNext];
	fSolOldTime->vec = &fSolOldTime->vec[kNext];
	fSolOldMassFracs->mat = &fSolOldMassFracs->mat[kNext];
	fSolOldTemp->vec = &fSolOldTemp->vec[kNext];
	fSolOldPressure->vec = &fSolOldPressure->vec[kNext];
	fSolOldR->vec = &fSolOldR->vec[kNext];
	fSolOldPsi->vec = &fSolOldPsi->vec[kNext];
	fGPDens->vec = &fGPDens->vec[kNext];
	fMassFracsWork->mat = &fMassFracsWork->mat[kNext];
	fTempWork->vec = &fTempWork->vec[kNext];
	fPressureWork->vec = &fPressureWork->vec[kNext];
	fRWork->vec = &fRWork->vec[kNext];
	fPsiWork->vec = &fPsiWork->vec[kNext];

	// set variable names
	fVariableNames = new String[fNOfEquations];

	fVariableNames[fR] = new char[2];
	strcpy( fVariableNames[fR], "R" );
	fVariableNames[fPsi] = new char[4];
	strcpy( fVariableNames[fPsi], "Psi" );
	fVariableNames[fPressure] = new char[9];
	strcpy( fVariableNames[fPressure], "Pressure" );
	fVariableNames[fTemperature] = new char[2];
	strcpy( fVariableNames[fTemperature], "T" );
	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies->GetNames()[i] ) + 1];
		strcpy( fVariableNames[fFirstSpecies + i], fSpecies->GetNames()[i] );
	}
	
	FILE *fpNames = GetOutfile( "VarNames", kText );
	for ( i = 0; i < fNOfEquations; ++i ) {
		fprintf( fpNames, "%s\n", fVariableNames[i] );
	}
	fclose( fpNames );

	CompLewisNumbers( fInputData->fLewisNumberFile );
	fflush( fOutFilePtr );
}

TTrans1DIsoChorSolver::~TTrans1DIsoChorSolver( void )
{
	int nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();

	for ( int i = 0; i < nOfSpeciesIn+fVariablesWithoutSpecies; ++i ) {
		delete fVariableNames[i];
	}
	delete fVariableNames;
	
	fPsiWork->vec = &fPsiWork->vec[kPrev];
	fRWork->vec = &fRWork->vec[kPrev];
	fPressureWork->vec = &fPressureWork->vec[kPrev];
	fTempWork->vec = &fTempWork->vec[kPrev];
	fMassFracsWork->mat = &fMassFracsWork->mat[kPrev];
	fGPDens->vec = &fGPDens->vec[kPrev];
	fSolOldPsi->vec = &fSolOldPsi->vec[kPrev];
	fSolOldR->vec = &fSolOldR->vec[kPrev];
	fSolOldPressure->vec = &fSolOldPressure->vec[kPrev];
	fSolOldTemp->vec = &fSolOldTemp->vec[kPrev];
	fSolOldMassFracs->mat = &fSolOldMassFracs->mat[kPrev];
	fSolOldTime->vec = &fSolOldTime->vec[kPrev];
	fSolPsi->vec = &fSolPsi->vec[kPrev];
	fSolR->vec = &fSolR->vec[kPrev];
	fSolPressure->vec = &fSolPressure->vec[kPrev];
	fSolTemp->vec = &fSolTemp->vec[kPrev];
	fSolMassFracs->mat = &fSolMassFracs->mat[kPrev];
	fSolTime->vec = &fSolTime->vec[kPrev];

	Free1DArray( fMixHeatCap );
	Free1DArray( fRho );

	DisposeVector( fOutSolWork );
	DisposeVector( fPsiWork );
	DisposeVector( fRWork );
	DisposeVector( fPressureWork );
	DisposeVector( fTempWork );
	DisposeMatrix( fMassFracsWork );
	DisposeVector( fGPDens );
	Free1DArray( fGPDistWeight );
	Free1DArray( fWMinus );
	Free1DArray( fWPlus );
	Free1DArray( fFDWMinus );
	Free1DArray( fFDWPlus );
	Free1DArray( fFDWCurr );
	Free1DArray( fhnenn );
	Free1DArray( fhm );
	Free1DArray( fh );
	DisposeVector( fSolOldPsi );
	DisposeVector( fSolOldR );
	DisposeVector( fSolOldPressure );
	DisposeVector( fSolOldTemp );
	DisposeMatrix( fSolOldMassFracs );
	DisposeVector( fSolOldTime );
	DisposeVector( fSolPsi );
	DisposeVector( fSolR );
	DisposeVector( fSolPressure );
	DisposeVector( fSolTemp );
	DisposeMatrix( fSolMassFracs );
	DisposeVector( fSolTime );

	DisposeMatrix( fSolPrime );
	fSolution->mat = &fSolution->mat[kPrev];
	DisposeMatrix( fSolution );

	DisposeIntVector( fIWrk[0] );
	DisposeVector( fRWrk[0] );
	DisposeIntVector( fInfo[0] );

	Free1DArray( fMaxStepTaken );
	delete fNActualOrd;
	delete fNActualStep;

	delete fIWrk;
	delete fRWrk;
	delete fInfo;
	DisposeVector( fTime );

	Free1DArray( fATol );
	Free1DArray( fRTol );
}

void TTrans1DIsoChorSolver::SetInfo( int i )
{
	int		*info = fInfo[i]->vec;

	info[kF1] = 0;	// first call
#ifdef VECTORTOL
	info[kF2] = 1;	// RTOL, ATOL are vectors
#else
	info[kF2] = 0;	// RTOL, ATOL are scalars
#endif
	info[kF3] = 1;	// solution at TOUT, no intermediates
	info[kF4] = 1;	// integration can be carried out beyond TOUT ( and interpolated )
//	if ( fUseNumericalDM ) {
		info[kF5] = 0;	// numerical jacobian
//	}
//	else {
//		info[kF5] = 1;
//	}
	info[kF6] = 1;	// full ( not banded ) jacobian
	info[kF7] = 1;	// no choice for maximum stepsize
	info[kF8] = 1;	// no choice for first stepsize
	if ( fMaxOrd == 5 ) {
		info[kF9] = 0;	// choose default for MAXORD ( = 5 )
	}
	else {
		info[kF9] = 1;	// don't choose default for MAXORD ( = 5 )
	}
	info[kF10] = 0;	// solve without non negative constraint
	info[kF11] = 1;	// initial values are consistent
}

void TTrans1DIsoChorSolver::InitDassl( int nGridPoints )
{

	for ( int i = 0; i < nGridPoints; ++i ) {
		SetInfo( i );
		InitDasslWorkSpace( i );
	}
}

void TTrans1DIsoChorSolver::InitDasslWorkSpace( int i )
{
	int		*info = fInfo[i]->vec;
	int		*iWork = fIWrk[i]->vec;
	Double	*rWork = fRWrk[i]->vec;

	if ( info[kF6] == 1 ) {
		iWork[kF1] = 2 * fNOfEquations;
		iWork[kF2] = 2 * fNOfEquations;
	}
	if ( info[kF8] == 1 ) {
		rWork[kF3] = fFirstTimeStep;
	}
	if ( info[kF9] == 1 ) {
		iWork[kF3] = fMaxOrd;
	}
	if ( info[kF7] == 1 ) {
		rWork[kF2] = fDeltaTMax;
	}
}

void TTrans1DIsoChorSolver::SetMaxTimeStep( int kAct, Double /*t*/ )
{
	fRWrk[kAct]->vec[kF2] = fDeltaTMax;
}

void TTrans1DIsoChorSolver::Initialize( Double timeStart
					, ConstStringArray names
					, Double **startSolution, int vars
					, Double *grid, int gridPointsA
					, Double pressureStart, Double phiStart
					, Double firstTimeStep )
{
	int 	k;
	Double	*time = fTime->vec;
	
	if ( firstTimeStep > 1.0e-13 ) {
		fprintf( fOutFilePtr, "initial timestep is %g s\n", firstTimeStep );
		for ( k = 0; k < fNOfSolver; ++k ) {
			fRWrk[k]->vec[kF3] = firstTimeStep;
		}
	}
	fTEnd = timeStart;
	if ( timeStart > 0.0 ) {
		fprintf( fOutFilePtr, "start at time %g s\n", timeStart );
	}
	for ( k = 0; k < fNOfSolver; ++k ) {
		time[k] = timeStart;
	}
	fPressStart = fPressEnd = pressureStart;
	fPhiEnd = phiStart;
	
	fRLeft = startSolution[kCurr][fR];
	fRRight = startSolution[gridPointsA-1][fR];

//	MakeGrid( &fSolGrid->vec[kPrev], fSolGrid->len, fRLeft, fRRight, fEquidistant );

	SetInitial( startSolution, gridPointsA, vars );
	SetInitialPsi();
	fPsiLeft = fSolution->mat[kPrev][fPsi];
	fPsiRight = fSolution->mat[fNGridPoints][fPsi];
	SetFDWeights( fSolution->mat[0] );
	
#ifdef DEBUGINITIAL
	FILE	*fp = GetOutfile( "InitialSol", kData );
	PrintSolution( fp, gridPointsA-2, fNOfEquations, fSolution->mat, fVariableNames );
	fclose( fp );
#endif
	
	fTempLeftEnd = fSolution->mat[kPrev][fTemperature];
	fTempRightEnd = fSolution->mat[fNGridPoints][fTemperature];
	fRLeftEnd = fSolution->mat[kPrev][fR];
	fRRightEnd = fSolution->mat[fNGridPoints][fR];
}

Flag TTrans1DIsoChorSolver::Solve( Double timeEnd, Double pressureEnd, Double phiEnd
					, Double temLeftEnd, Double tempRightEnd, int deltaStepsOut )
{
	Flag	error = FALSE;
	Flag	leave = FALSE;
		
	SetEndValues( timeEnd, pressureEnd, phiEnd, temLeftEnd, tempRightEnd );

	if ( fFirstCall ) {
		leave = FirstImpliStep();
	}
	error = leave;
	if ( fTime->vec[0] + MIN( ( fTEnd - fTStart ) * 1.0e-7, 1.0e-10 ) >= fTEnd ) {
		leave = TRUE;
	}
#	ifdef HP
		if ( gExit ) {
			DoExit();
		}
#	endif
	if ( leave ) {
		return error;
	}
	
	do {
		leave = OneImpliStep();
		if ( deltaStepsOut > 0 && *fNActualStep[0] % deltaStepsOut == 0 ) {
			FILE	*fp = GetOutputFile( fTime->vec[0], NULL, NULL, TFlame::kText );
			WriteFlameletFile( fp, NULL, NULL );
			fclose( fp );
		}
		error = leave;
		if ( fTime->vec[0] + MIN( ( fTEnd - fTStart ) * 1.0e-7, 1.0e-12 ) >= fTEnd ) {
			leave = TRUE;
		}
#	ifdef HP
		if ( gExit ) {
			DoExit();
		}
#	endif
	} while ( !leave );

#ifdef DEBUGYPRIME
	FILE	*fp = GetOutfile( "YPrime", kData );
	PrintSolution( fp, fNGridPoints-2, fNOfEquations, &fSolPrime->mat[1], fVariableNames );
	fclose( fp );
#endif

#	ifdef DEBUGWEIGHTS
	fp = GetOutfile( "GridDistWeights", kData );
	fprintf( fp, "*\nx\tw\n" );
	for ( int k = 0; k < fNGridPoints; ++k ) {
		fprintf( fp, "%g\t%g\n", fSolution->mat[k][fPsi], fGPDistWeight[k] );
	}
	fclose( fp );
#	endif
	
	return error;
}

void TTrans1DIsoChorSolver::SetEndValues( Double timeEnd, Double pressureEnd, Double phiEnd
					, Double tempLeftEnd, Double tempRightEnd )
{
	fTStart = fTEnd;
	fPressStart = fPressEnd;
	fPhiStart = fPhiEnd;
	fTempLeftStart = fTempLeftEnd;
	fTempRightStart = fTempRightEnd;

	fTEnd = timeEnd;
	fPressEnd = pressureEnd;
	fPhiEnd = phiEnd;
	fTempLeftEnd = tempLeftEnd;
	fTempRightEnd = tempRightEnd;
	
	fprintf( fOutFilePtr, "solve system from t = %g with p = %g\tphi = %g\tTleft = %g\tTright = %g\n"
		, fTStart, fPressStart, fPhiStart, fTempLeftStart, fTempRightStart );
	fprintf( fOutFilePtr, "                 to t = %g with p = %g\tphi = %g\tTleft = %g\tTright = %g\n"
		, fTEnd, fPressEnd, fPhiEnd, fTempLeftEnd, fTempRightEnd );

	fSolOldTime->vec[kPrev] = fSolOldTime->vec[fNGridPoints] = fTStart;
	fSolTime->vec[kPrev] = fSolTime->vec[fNGridPoints] = fTEnd;

	fSolOldTemp->vec[kPrev] = fTempLeftStart;
	fSolOldTemp->vec[fNGridPoints] = fTempRightStart;

	fSolTemp->vec[kPrev] = fTempLeftEnd;
	fSolTemp->vec[fNGridPoints] = fTempRightEnd;

#ifdef DPDT
	fDPdt = ( fPressEnd - fPressStart ) / ( fTEnd - fTStart );
#else
	fDPdt = 0.0;
#endif

	for ( int k = 0; k < fNOfSolver; ++k ) {
		if ( fInfo[k]->vec[kF4] == 1 ) {
			fRWrk[k]->vec[kF1] = fTEnd;
		}
	}
}

void TTrans1DIsoChorSolver::SetInitial( Double **startSol, int gridPointsA, int vars )
{
	int i;
	Double **sol = &fSolution->mat[-1];

	for ( i = 0; i < (gridPointsA) * vars; ++i ) {
		sol[i] = startSol[i];
	}

// save initial solution to fSol
	for ( i = -1; i <= fNGridPoints; ++i ) {
		SaveSolution( i, fTime->vec[0], fSolution->mat[i] );
		// second call puts solution down to fSolOld
		SaveSolution( i, fTime->vec[0], fSolution->mat[i] ); 
	}
}

void TTrans1DIsoChorSolver::SetOutSolution( void )
{
	int			i, k;
	int			nOfSpecies = fSpecies->GetNOfSpecies();
	Double		**theY = fMassFracsWork->mat;
	Double		*theTemp = fTempWork->vec;
	Double		*thePressure = fPressureWork->vec;
	Double		*theR = fRWork->vec;
	Double		*thePsi = fPsiWork->vec;
	Double		*oTemp = fSolOldTemp->vec;
	Double		*nTemp = fSolTemp->vec;
	Double		*oPressure = fSolOldPressure->vec;
	Double		*nPressure = fSolPressure->vec;
	Double		*oR = fSolOldR->vec;
	Double		*nR = fSolR->vec;
	Double		*oPsi = fSolOldPsi->vec;
	Double		*nPsi = fSolPsi->vec;
	Double		**oY = fSolOldMassFracs->mat;
	Double		**nY = fSolMassFracs->mat;
	Double		*oTime = fSolOldTime->vec;
	Double		*nTime = fSolTime->vec;
	Double		currTime = GetCurrentTime();
	Double		oTimek, nTimek, *theYk, *oYk, *nYk;

	for ( k = -1; k <= fNGridPoints; ++k ) {
		oTimek = oTime[k];
		nTimek = nTime[k];
		theYk = theY[k];
		oYk = oY[k];
		nYk = nY[k];
		theTemp[k] = Interpol( currTime, oTemp[k], oTimek, nTemp[k], nTimek );
		thePressure[k] = Interpol( currTime, oPressure[k], oTimek, nPressure[k], nTimek );
		theR[k] = Interpol( currTime, oR[k], oTimek, nR[k], nTimek );
		thePsi[k] = Interpol( currTime, oPsi[k], oTimek, nPsi[k], nTimek );
		for ( i = 0; i < nOfSpecies; ++i ) {
			theYk[i] = Interpol( currTime, oYk[i], oTimek, nYk[i], nTimek );
		}
	}
}

void TTrans1DIsoChorSolver::PrintSolution( FILE *fp, int nGPoints, int nEq
									, Double **sol, ConstStringArray names )
{
	int		i, k;
	char	format[128];
	int 	*len = new int[nEq];
	
	fprintf( fp, "*\n" );
	for ( i = 0; i < nEq; ++i ) {
		len[i] = maxint( strlen( names[i] ), 12 );
		sprintf( format, "%s%%-%ds", "\t", len[i] );
		fprintf( fp, format, names[i] );
	}
	fprintf( fp, "\n" );
	
	for ( k = -1; k < nGPoints+1; ++k ) {
		for ( i = 0; i < nEq; ++i ) {
			sprintf( format, "%s%%-%dE", "\t", len[i] );
            fprintf( fp, format, sol[k][i] );
		}
		fprintf( fp, "\n" );
	}
	
	delete len;
}

void makegrid( int **object, Double *grid, int *gridPoints
							, Double *left, Double *right, int *equidistant )
{
	TTrans1DIsoChorSolverPtr flame = *( TTrans1DIsoChorSolverPtr *)object;

	flame->MakeGrid( grid, *gridPoints, *left, *right, ( *equidistant ) ? TRUE : FALSE );
}

void TTrans1DIsoChorSolver::MakeGrid( VectorPtr theGrid, Double left, Double right, Flag equidistant )
{
	MakeGrid( theGrid->vec, theGrid->len, left, right, equidistant );
}

void TTrans1DIsoChorSolver::MakeGrid( Double *grid, int len
							, Double left, Double right, Flag equidistant )
{
	int			 	gridPoints = len - 2;
	Double			delta;
	
	if ( right <= left ) {
		cerr << "#warning: MakeGrid called with invalid parameters, using defaults\n";
		left = 0.0;
		right = 1.0;
	}
	
	delta = ( right - left ) / ( ( Double ) ( gridPoints + 1 ) );

	grid[0] = left;
	grid[gridPoints+1] = right;

	for ( int k = 1; k <= gridPoints; ++k ) {
		grid[k] = grid[k-1] + delta;
	}
	
#	ifdef DEBUGGRID
	static int	gridCount = 0;
	char		fName[128];
	sprintf( fName, "grid%2d", gridCount++ );
	FILE		*fp = GetOutfile( fName, kData );
	
	fprintf( fp, "*\nx\td\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "%g\t1.0\n", grid[k] );
	}
	
	fclose( fp );
#	endif
}

void TTrans1DIsoChorSolver::SetFDWeights( Double *sol )
{
	int		k = 0;
	int		gpOff = k * fNOfEquations;
	fh[k] = sol[gpOff+fPsi+fNOfEquations] - sol[gpOff+fPsi];
	fhm[k] = sol[gpOff+fPsi] - fPsiLeft;

	for ( k = 1; k < fNGridPoints-1; ++k ) {
		gpOff = k * fNOfEquations;
		fh[k] = sol[gpOff+fPsi+fNOfEquations] - sol[gpOff+fPsi];
		fhm[k] = sol[gpOff+fPsi] - sol[gpOff+fPsi-fNOfEquations];
	}

	k = fNGridPoints-1;
	gpOff = k * fNOfEquations;
	fh[k] = fPsiRight - sol[gpOff+fPsi];
	fhm[k] = sol[gpOff+fPsi] - sol[gpOff+fPsi-fNOfEquations];

	for ( k = 0; k < fNGridPoints; ++k ) {
		fhnenn[k] = fh[k] * fhm[k] * ( fh[k] + fhm[k] );
		fFDWCurr[k] = ( fh[k] * fh[k] - fhm[k] * fhm[k] ) / fhnenn[k];
		fFDWMinus[k] = - fh[k] * fh[k] / fhnenn[k];
		fFDWPlus[k] = fhm[k] * fhm[k] / fhnenn[k];
		fWMinus[k] = - fh[k] / fhnenn[k];
		fWPlus[k] = fhm[k] / fhnenn[k];
		fGPDens->vec[k] = 1.0 / fh[k];
	}
	fGPDens->vec[kPrev] = 1.0 / ( sol[fPsi] - fPsiLeft );

#	ifdef DEBUGGPWEIGHTS
	FILE		*fp = GetOutfile( "GridPointWeights", kData );
	fprintf( fp, "*\nx\tw\n" );
	for ( k = -1; k < fNGridPoints+1; ++k ) {
		fprintf( fp, "%g\t%g\n", sol[k][fPsi], fGPDens->vec[k] );
	}
	fclose( fp );
#	endif
}

int TTrans1DIsoChorSolver::GetActualPoint( Double tEnd )
{
	int		minPoint = 0;
	Double	*time = fSolTime->vec;
	Double	minValue = time[0];
	
	for ( int k = 1; k < fNGridPoints; ++k ) {
		if ( time[k] < minValue ) {
			minPoint = k;
			minValue = time[k];
		}
	}
	
	if ( minValue < tEnd ) {
		return minPoint;
	}
	else {
		return -1;
	}
}

void TTrans1DIsoChorSolver::SaveSolution( int k, Double time, Double *y )
{
	int		nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double	*YNew = fSolMassFracs->mat[k];
	Double	*YOld = fSolOldMassFracs->mat[k];
	
	// save old solution
	fSolOldTime->vec[k] = fSolTime->vec[k];
	copy( nOfSpeciesIn, YNew, 1, YOld, 1 );
	fSolOldTemp->vec[k] = fSolTemp->vec[k];
	fSolOldPressure->vec[k] = fSolPressure->vec[k];
	fSolOldR->vec[k] = fSolR->vec[k];
	fSolOldPsi->vec[k] = fSolPsi->vec[k];

	// save new solution
	fSolTime->vec[k] = time;
	copy( nOfSpeciesIn, &y[fFirstSpecies], 1, YNew, 1 );

	fSolTemp->vec[k] = y[fTemperature];
	fSolPressure->vec[k] = y[fPressure];
	fSolR->vec[k] = y[fR];
	fSolPsi->vec[k] = y[fPsi];
}

Flag TTrans1DIsoChorSolver::FirstImpliStep( void )
{
	Flag	leave;
	
	leave = OneImpliStep();
	
	fInfo[0]->vec[kF11] = 0;

	fFirstCall = FALSE;
	return leave;
}

Flag TTrans1DIsoChorSolver::OneImpliStep( void )
{
#ifdef UNICOS
	SDASSL( ::ResTrans1DIsoChorImpliSolver, &fDasslNEq, &fTime->vec[0], fSolution->mat[0]
			, fSolPrime->mat[0], &fTEnd, fInfo[0]->vec, fRTol, fATol, &fIdid
			, fRWrk[0]->vec, &fLRW, fIWrk[0]->vec, &fLIW, NULL, ( int * ) this
			, NULL );
#else
	DDASSL( ::ResTrans1DIsoChorImpliSolver, &fDasslNEq, &fTime->vec[0], fSolution->mat[0]
			, fSolPrime->mat[0], &fTEnd, fInfo[0]->vec, fRTol, fATol, &fIdid
			, fRWrk[0]->vec, &fLRW, fIWrk[0]->vec, &fLIW, NULL, ( int * ) this
			, NULL );
#endif

	fMaxStepTaken[0] = MAX( fMaxStepTaken[0], fRWrk[0]->vec[kF7] );
	if ( *fNActualStep[0] % DELTAPROG == 0 ) {
		fprintf( fOutFilePtr, "stp = %5d    dt = %8.2e s    ord = %1d    time = %.4g ms\ttemp = %6.1f K\n"
			, *fNActualStep[0], fRWrk[0]->vec[kF7], *fNActualOrd[0], fTime->vec[0] * 1000.0
			, fSolution->mat[LocationOfMax( fNGridPoints, &fSolution->mat[0][fTemperature], fSolution->phys_rows )][fTemperature] );
		fflush( fOutFilePtr );
	}
	

	if ( fIdid < 1 ) {
		fprintf( fOutFilePtr, "#error: ddassl error no. %d occured at timestep no. %d\n"
			, fIdid, *fNActualStep[0] );
		fprintf( fOutFilePtr, "number of error test failures is %d\n", fIWrk[0]->vec[kF14] );
		fprintf( fOutFilePtr, "number of convergence test failures is %d\n", fIWrk[0]->vec[kF15] );
		fprintf( fOutFilePtr, "dassl suggested order %d and stepsize %g for the next step\n"
			,fIWrk[0]->vec[kF7], fRWrk[0]->vec[kF3] );
/*		cerr << "#error: ddassl error no. " << fIdid << " occured in timestep no. " 
				<< *fNActualStep[0] << NEWL;
		cerr << "number of error test failures is " << fIWrk[0]->vec[kF14] << NEWL;
		cerr << "number of convergence test failures is " << fIWrk[0]->vec[kF15] << NEWL;
		cerr << "dassl suggested order " << fIWrk[0]->vec[kF7] << " and stepsize " 
				<< fRWrk[0]->vec[kF3] << " for the next step" << NEWL;*/
		return TRUE;
	}
	else {
		if ( fRWrk[0]->vec[kF3] > fTime->vec[0] ) {
			fprintf( fOutFilePtr, "#warning: integration carried out beyond tEnd\n" );
			fprintf( fOutFilePtr, "info[4] = %d\n", fInfo[0]->vec[kF4] );
			fprintf( fOutFilePtr, "fIdid = %d\n", fIdid );
/*			cerr << "#warning: integration carried out beyond tEnd" << NEWL;
			cerr << "info[4] = " << fInfo[0]->vec[kF4] << NEWL;
			cerr << "fIdid = " << fIdid << NEWL;*/
		}

		for ( int k = 0; k < fNGridPoints; ++k ) {
			SaveSolution( k, fTime->vec[0], fSolution->mat[k] );
		}
		return FALSE;
	}
}

void ResTrans1DIsoChorImpliSolver( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar )
{
	( ( TTrans1DIsoChorSolverPtr ) iPar )->ResTrans1DIsoChorImpliSolver( T, y, yPrime, delta
			, iRes, rPar, iPar );
}

void TTrans1DIsoChorSolver::ResTrans1DIsoChorImpliSolver( Double *time, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double */*rPar*/, int */*iPar*/ )
{
	int						k;
	int						nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	int						i, gpOff, gpROff, gpPsiOff, gpPressureOff, gpTempOff, gpSpecOff;
	Double					sumEnthFlux, sumHeatSource;

	Double					*prodRate = fSpecies->GetProductionRate()->vec;
	Double					*enth = fSpecies->GetEnthalpy()->vec;
	Double					*molarMass = fSpecies->GetMolarMass()->vec;
	Double					*heatCap = fSpecies->GetHeatCapacity()->vec;
	Double					*Le = fSpecies->GetLewisNumber()->vec;

	Double					rCurr, rPrev;
	Double					tempCurr, tempPrev, tempNext;
	Double					pressureCurr, pressureNext;
	Double					diffCurr, diffPrev, diffNext;
	Double					*gpDens = fGPDens->vec;
	Double					**nY = fSolMassFracs->mat;

	Double					kappa_ = fKappa * ( 1.0 + fKappa );
	
	// copy YF to nY
	for ( k = 0; k < fNGridPoints; ++k ) {
		copy( nSpeciesIn, &y[k*fNOfEquations+fFirstSpecies], 1, nY[k], 1 );
		if ( y[k*fNOfEquations+fTemperature] < 100.0 || y[k*fNOfEquations+fTemperature] > 5000.0 ) {
			fprintf( fOutFilePtr, "###Warning: Temp = %g @ gp = %d \n", y[k*fNOfEquations+fTemperature], k );
			*iRes = -1;
			return;
		}
		ComputeProperties( y[k*fNOfEquations+fTemperature], nY[k], y[k*fNOfEquations+fPressure], 
						   fRho[k], kDensFromPress, NULL );
		fMixHeatCap[k] = fProperties->GetMixHeatCapacity();
	}

	SetFDWeights( y );
	ComputeGPDistributionWeights( y );
		
	// first point
	k = 0;
	gpOff = k * fNOfEquations;
	gpROff = gpOff + fR;
	gpPsiOff = gpOff + fPsi;
	gpPressureOff = gpOff + fPressure;
	gpSpecOff = gpOff + fFirstSpecies;
	gpTempOff = gpOff + fTemperature;
	rCurr = y[gpROff];
	rPrev = fRLeftEnd;
	pressureCurr = y[gpPressureOff];
	pressureNext = y[gpPressureOff+fNOfEquations];
	tempCurr = y[gpTempOff];
	tempNext = y[gpTempOff+fNOfEquations];
	diffCurr = fRho[k] * modPow( rCurr, 2 * fGamma ) * LambdaOverCp( tempCurr );
	diffNext = fRho[k+1] * modPow( y[gpROff+fNOfEquations], 2 * fGamma ) * LambdaOverCp( tempNext );

	UpdateThermoProps( nY[k], tempCurr, pressureCurr, fRho[k], kDensFromPress, NULL );

	// psi equation
	delta[gpPsiOff] = - 2.0 * yPrime[gpPsiOff] + yPrime[gpPsiOff+fNOfEquations]; 

	// r equation
	delta[gpROff] = ( rCurr - rPrev ) / fhm[k] - 1.0 / ( fRho[k] * modPow( rCurr, fGamma ) ); 

	// pressure equation
	delta[gpPressureOff] = ( pressureNext - pressureCurr ) / fh[k]; 

	// species equation
	sumEnthFlux = 0.0;
	sumHeatSource = 0.0;
	for ( i = 0; i < nSpeciesIn; ++i ) {
		delta[gpSpecOff+i] = yPrime[gpSpecOff+i]
//			- ( nY[k+1][i] - nY[k][i] ) * fFDWPlus[k] * yPrime[gpPsiOff]
//			- ( nY[k+1][i] - nY[k][i] ) / ( fh[k] + fhm[k] ) * yPrime[gpPsiOff]
			- ( diffCurr + diffNext ) / Le[i] * ( nY[k+1][i] - nY[k][i] ) * fWPlus[k]
			- prodRate[i] / fRho[k];
		if ( ( -yPrime[gpPsiOff] ) > 0.0 ) {
			delta[gpSpecOff+i] += 0.0;
		}
		else {
			delta[gpSpecOff+i] += - ( nY[k+1][i] - nY[k][i] ) / fh[k] * yPrime[gpPsiOff];
		}
		sumHeatSource += enth[i] * prodRate[i];
		sumEnthFlux += heatCap[i] / Le[i] * fFDWPlus[k] * ( nY[k+1][i] - nY[k][i] );
	}

	// temperature equation
	delta[gpTempOff] = yPrime[gpTempOff]
//			- ( tempNext - tempCurr ) * fFDWPlus[k] * yPrime[gpPsiOff]
//			- ( tempNext - tempCurr ) / ( fh[k] + fhm[k] ) * yPrime[gpPsiOff]
			- 1.0 / ( fRho[k] * fMixHeatCap[k] ) * yPrime[gpPressureOff]
//			- 1.0 / ( fRho[k] * fMixHeatCap[k] ) * ( yPrime[gpPressureOff]
//			  - ( pressureNext - pressureCurr ) * fFDWPlus[k] * yPrime[gpPsiOff] )
//			  - ( pressureNext - pressureCurr ) / ( fh[k] + fhm[k] ) * yPrime[gpPsiOff] )
			- ( diffCurr + diffNext ) * ( tempNext - tempCurr ) * fWPlus[k]
			- diffCurr / fMixHeatCap[k] * fFDWPlus[k] * ( fMixHeatCap[k+1] - fMixHeatCap[k] )
									   * fFDWPlus[k] * ( tempNext - tempCurr)
			- diffCurr / fMixHeatCap[k] * sumEnthFlux * fFDWPlus[k] * ( tempNext - tempCurr)
			+ 1.0 / ( fRho[k] * fMixHeatCap[k] ) * sumHeatSource; 
	if ( ( -yPrime[gpPsiOff] ) > 0.0 ) {
		delta[gpTempOff] += 0.0;
	}
	else {
		delta[gpTempOff] += - ( tempNext - tempCurr ) / fh[k] * yPrime[gpPsiOff];
		delta[gpTempOff] += 1.0 / ( fRho[k] * fMixHeatCap[k] ) * ( pressureNext - pressureCurr ) / fh[k] * yPrime[gpPsiOff];
	}
	if ( fProperties->GetRadiation() ) {
		delta[gpTempOff] -= ( fProperties->GetRadiation()->GetRadiation() ) / ( fRho[k] * fMixHeatCap[k] );
	}
	if ( fArtificialSource ) {
		delta[gpTempOff] -= GetTempSource( rCurr, time[kCurr] ) / ( fRho[k] * fMixHeatCap[k] );
	}

// for debug purposes
#ifdef DEBUGRES
		fprintf( stdout, "gp = %d\n", k );
		for ( int j = 0; j < nSpeciesIn+fVariablesWithoutSpecies; ++j ) {
			fprintf( stdout, "%s\t%g\t%g\n", fVariableNames[j]
									, y[gpOff + j]
									, delta[gpOff + j] );
		}
		fprintf( stdout, "\n" );
#endif

	
	// rest points
	for ( k = 1; k < fNGridPoints-1; ++k ) {
		gpOff = k * fNOfEquations;
		gpROff = gpOff + fR;
		gpPsiOff = gpOff + fPsi;
		gpPressureOff = gpOff + fPressure;
		gpSpecOff = gpOff + fFirstSpecies;
		gpTempOff = gpOff + fTemperature;
		rCurr = y[gpROff];
		rPrev = y[gpROff-fNOfEquations];
		pressureCurr = y[gpPressureOff];
		pressureNext = y[gpPressureOff+fNOfEquations];
		tempCurr = y[gpTempOff];
		tempPrev = y[gpTempOff-fNOfEquations];
		tempNext = y[gpTempOff+fNOfEquations];
		diffCurr = fRho[k] * modPow( rCurr, 2 * fGamma ) * LambdaOverCp( tempCurr );
		diffPrev = fRho[k-1] * modPow( rPrev, 2 * fGamma ) * LambdaOverCp( tempPrev );
		diffNext = fRho[k+1] * modPow( y[gpROff+fNOfEquations], 2 * fGamma ) * LambdaOverCp( tempNext );
	
		UpdateThermoProps( nY[k], tempCurr, pressureCurr, fRho[k], kDensFromPress, NULL );

		// psi equation
		delta[gpPsiOff] = fTauGrid * (
			( kappa_ * gpDens[k-2] * gpDens[k-2] / fGPDistWeight[k-1] ) * yPrime[gpPsiOff-2*fNOfEquations]
		  + ( - kappa_ * gpDens[k-1] * gpDens[k-1] / fGPDistWeight[k]
		  	+ ( - kappa_ * gpDens[k-2] * gpDens[k-2] 
			  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fGPDistWeight[k-1] ) * yPrime[gpPsiOff-fNOfEquations]
		  + ( ( kappa_ * gpDens[k-1] * gpDens[k-1]
		      + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fGPDistWeight[k] 
			+ ( kappa_ * gpDens[k] * gpDens[k] 
			  + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] * gpDens[k-1] ) / fGPDistWeight[k-1] ) * yPrime[gpPsiOff]
		  + ( - kappa_ * gpDens[k] * gpDens[k] / fGPDistWeight[k-1]
		  	+ ( - kappa_ * gpDens[k+1] * gpDens[k+1] 
			  - ( 1.0 + 2.0 * kappa_ ) * gpDens[k] * gpDens[k] ) / fGPDistWeight[k] ) * yPrime[gpPsiOff+fNOfEquations]
		  + ( kappa_ * gpDens[k+1] * gpDens[k+1] / fGPDistWeight[k] ) * yPrime[gpPsiOff+2*fNOfEquations]
		  )
		  + ( - kappa_ * gpDens[k-1] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k] - kappa_ * gpDens[k+1] ) / fGPDistWeight[k]
		  - ( - kappa_ * gpDens[k-2] + ( 1.0 + 2.0 * kappa_ ) * gpDens[k-1] - kappa_ * gpDens[k] ) / fGPDistWeight[k-1];
	
		// r equation
		delta[gpROff] = ( rCurr - rPrev ) / fhm[k] - 1.0 / ( fRho[k] * modPow( rCurr, fGamma ) ); 
	
		// pressure equation
		delta[gpPressureOff] = ( pressureNext - pressureCurr ) / fh[k]; 
	
		// species equation
		sumEnthFlux = 0.0;
		sumHeatSource = 0.0;
		for ( i = 0; i < nSpeciesIn; ++i ) {
			delta[gpSpecOff+i] = yPrime[gpSpecOff+i] 
//				- ( ( nY[k+1][i] - nY[k][i] ) * fFDWPlus[k] + 
//				    ( nY[k][i] - nY[k-1][i] ) * fFDWMinus[k] ) * yPrime[gpPsiOff]
//				- ( nY[k][i+1] - nY[k-1][i] ) / ( fh[k] + fhm[k] ) * yPrime[gpPsiOff]
				- ( ( diffCurr + diffNext ) / Le[i] * ( nY[k+1][i] - nY[k][i] ) * fWPlus[k]
				  + ( diffPrev + diffCurr ) / Le[i] * ( nY[k][i] - nY[k-1][i] ) * fWMinus[k] )
				- prodRate[i] / fRho[k];
			if ( ( -yPrime[gpPsiOff] ) > 0.0 ) {
				delta[gpSpecOff+i] += - ( nY[k][i] - nY[k-1][i] ) / fhm[k] * yPrime[gpPsiOff];
			}
			else {
				delta[gpSpecOff+i] += - ( nY[k+1][i] - nY[k][i] ) / fh[k] * yPrime[gpPsiOff];
			}
			sumHeatSource += enth[i] * prodRate[i];
			sumEnthFlux += heatCap[i] / Le[i] * ( fFDWPlus[k] * ( nY[k+1][i] - nY[k][i] )
												+ fFDWMinus[k] * ( nY[k][i] - nY[k-1][i] ) );
		}
	
		// temperature equation
		delta[gpTempOff] = yPrime[gpTempOff]
//				- ( ( tempNext - tempCurr ) * fFDWPlus[k]
//				  + ( tempCurr - tempPrev ) * fFDWMinus[k] ) * yPrime[gpPsiOff]
//				- ( tempNext - tempPrev ) / ( fh[k] + fhm[k] ) * yPrime[gpPsiOff]
				- 1.0 / ( fRho[k] * fMixHeatCap[k] ) * yPrime[gpPressureOff]
//				- 1.0 / ( fRho[k] * fMixHeatCap[k] ) * ( yPrime[gpPressureOff]
//				  - ( ( pressureNext - pressureCurr ) * fFDWPlus[k]
//				    + ( pressureCurr - y[gpPressureOff-fNOfEquations] ) * fFDWMinus[k] ) * yPrime[gpPsiOff] )
//				  - ( pressureNext - y[gpPressureOff-fNOfEquations] ) / ( fh[k] + fhm[k] ) * yPrime[gpPsiOff] )
				- ( ( diffCurr + diffNext ) * ( tempNext - tempCurr ) * fWPlus[k]
				  + ( diffPrev + diffCurr ) * ( tempCurr - tempPrev ) * fWMinus[k] )
				- diffCurr / fMixHeatCap[k] * ( fFDWPlus[k] * ( fMixHeatCap[k+1] - fMixHeatCap[k] )
										  + fFDWMinus[k] * ( fMixHeatCap[k] - fMixHeatCap[k-1] ) )
										* ( fFDWPlus[k] * ( tempNext - tempCurr)
										  + fFDWMinus[k] * ( tempCurr - tempPrev ) )
				- diffCurr / fMixHeatCap[k] * sumEnthFlux * ( fFDWPlus[k] * ( tempNext - tempCurr)
														+ fFDWMinus[k] * ( tempCurr - tempPrev ) )
				+ 1.0 / ( fRho[k] * fMixHeatCap[k] ) * sumHeatSource; 
		if ( ( -yPrime[gpPsiOff] ) > 0.0 ) {
			delta[gpTempOff] += - ( tempCurr - tempPrev ) / fhm[k] * yPrime[gpPsiOff];
			delta[gpTempOff] += 1.0 / ( fRho[k] * fMixHeatCap[k] ) * ( pressureCurr - y[gpPressureOff-fNOfEquations] ) / fhm[k] * yPrime[gpPsiOff];
		}
		else {
			delta[gpTempOff] += - ( tempNext - tempCurr ) / fh[k] * yPrime[gpPsiOff];
			delta[gpTempOff] += 1.0 / ( fRho[k] * fMixHeatCap[k] ) * ( pressureNext - pressureCurr ) / fh[k] * yPrime[gpPsiOff];
		}
		if ( fProperties->GetRadiation() ) {
			delta[gpTempOff] -= ( fProperties->GetRadiation()->GetRadiation() ) / ( fRho[k] * fMixHeatCap[k] );
		}
		if ( fArtificialSource ) {
			delta[gpTempOff] -= GetTempSource( rCurr, time[kCurr] ) / ( fRho[k] * fMixHeatCap[k] );
		}

// for debug purposes
#ifdef DEBUGRES
		fprintf( stdout, "gp = %d\n", k );
		for ( j = 0; j < nSpeciesIn+fVariablesWithoutSpecies; ++j ) {
			fprintf( stdout, "%s\t%g\t%g\n", fVariableNames[j]
									, y[gpOff + j]
									, delta[gpOff + j] );
		}
		fprintf( stdout, "\n" );
#endif

	}

	// last point
	k = fNGridPoints-1;
	gpOff = k * fNOfEquations;
	gpROff = gpOff + fR;
	gpPsiOff = gpOff + fPsi;
	gpPressureOff = gpOff + fPressure;
	gpSpecOff = gpOff + fFirstSpecies;
	gpTempOff = gpOff + fTemperature;
	rCurr = y[gpROff];
	rPrev = y[gpROff-fNOfEquations];
	pressureCurr = y[gpPressureOff];
	tempCurr = y[gpTempOff];
	tempPrev = y[gpTempOff-fNOfEquations];
	diffCurr = fRho[k] * modPow( rCurr, 2 * fGamma ) * LambdaOverCp( tempCurr );
	diffPrev = fRho[k-1] * modPow( rPrev, 2 * fGamma ) * LambdaOverCp( tempPrev );

	UpdateThermoProps( nY[k], tempCurr, pressureCurr, fRho[k], kDensFromPress, NULL );

	// psi equation
	delta[gpPsiOff] = yPrime[gpPsiOff-fNOfEquations] - 2.0 * yPrime[gpPsiOff]; 

	// r equation
	delta[gpROff] = ( rCurr - rPrev ) / fhm[k] - 1.0 / ( fRho[k] * modPow( rCurr, fGamma ) ); 

	// pressure equation
	delta[gpPressureOff] = ( fRRightEnd - rCurr ) / fh[k] - 1.0 / ( fRho[k] * modPow( fRRightEnd, fGamma ) ); 

	// species equation
	sumEnthFlux = 0.0;
	sumHeatSource = 0.0;
	for ( i = 0; i < nSpeciesIn; ++i ) {
		delta[gpSpecOff+i] = yPrime[gpSpecOff+i] 
//			- ( nY[k][i] - nY[k-1][i] ) * fFDWMinus[k] * yPrime[gpPsiOff]
//			- ( nY[k][i] - nY[k-1][i] ) / ( fh[k] + fhm[k] ) * yPrime[gpPsiOff]
			- ( diffPrev + diffCurr ) / Le[i] * ( nY[k][i] - nY[k-1][i] ) * fWMinus[k]
			- prodRate[i] / fRho[k];
		if ( ( -yPrime[gpPsiOff] ) > 0.0 ) {
			delta[gpSpecOff+i] += - ( nY[k][i] - nY[k-1][i] ) / fhm[k] * yPrime[gpPsiOff];
		}
		else {
			delta[gpSpecOff+i] += 0.0;
		}
		sumHeatSource += enth[i] * prodRate[i];
		sumEnthFlux += heatCap[i] / Le[i] * fFDWMinus[k] * ( nY[k][i] - nY[k-1][i] );
	}

	// temperature equation
	delta[gpTempOff] = yPrime[gpTempOff]
//			- ( tempCurr - tempPrev ) * fFDWMinus[k] * yPrime[gpPsiOff]
//			- ( tempCurr - tempPrev ) / ( fh[k] + fhm[k] ) * yPrime[gpPsiOff]
			- 1.0 / ( fRho[k] * fMixHeatCap[k] ) * yPrime[gpPressureOff]
//			- 1.0 / ( fRho[k] * fMixHeatCap[k] ) * ( yPrime[gpPressureOff]
//			  - ( pressureCurr - y[gpPressureOff-fNOfEquations] ) * fFDWMinus[k] * yPrime[gpPsiOff] )
//			  - ( pressureCurr - y[gpPressureOff-fNOfEquations] ) / ( fh[k] + fhm[k] ) * yPrime[gpPsiOff] )
			- ( diffPrev + diffCurr ) * ( tempCurr - tempPrev ) * fWMinus[k]
			- diffCurr / fMixHeatCap[k] * fFDWMinus[k] * ( fMixHeatCap[k] - fMixHeatCap[k-1] )
									* fFDWMinus[k] * ( tempCurr - tempPrev )
			- diffCurr / fMixHeatCap[k] * sumEnthFlux * fFDWMinus[k] * ( tempCurr - tempPrev )
			+ 1.0 / ( fRho[k] * fMixHeatCap[k] ) * sumHeatSource; 
	if ( ( -yPrime[gpPsiOff] ) > 0.0 ) {
		delta[gpTempOff] += - ( tempCurr - tempPrev ) / fhm[k] * yPrime[gpPsiOff];
		delta[gpTempOff] += 1.0 / ( fRho[k] * fMixHeatCap[k] ) * ( pressureCurr - y[gpPressureOff-fNOfEquations] ) / fhm[k] * yPrime[gpPsiOff];
	}
	else {
		delta[gpTempOff] += 0.0;
	}
	if ( fProperties->GetRadiation() ) {
		delta[gpTempOff] -= ( fProperties->GetRadiation()->GetRadiation() ) / ( fRho[k] * fMixHeatCap[k] );
	}
	if ( fArtificialSource ) {
		delta[gpTempOff] -= GetTempSource( rCurr, time[kCurr] ) / ( fRho[k] * fMixHeatCap[k] );
	}

// for debug purposes
#ifdef DEBUGRES
		fprintf( stdout, "gp = %d\n", k );
		for ( j = 0; j < nSpeciesIn+fVariablesWithoutSpecies; ++j ) {
			fprintf( stdout, "%s\t%g\t%g\n", fVariableNames[j]
									, y[gpOff + j]
									, delta[gpOff + j] );
		}
		fprintf( stdout, "\n" );
		exit( 2 );
#endif
}


Double TTrans1DIsoChorSolver::Interpol( Double t, Double valOld, Double tOld, Double valNew, Double tNew )
{
	if ( valNew == valOld ) {
		return valNew;
	}
	return valOld + ( valNew - valOld ) / ( tNew - tOld ) * ( t - tOld );
}

Double TTrans1DIsoChorSolver::GetPressure( void )
{
	return fSolPressure->vec[0];
}

int	TTrans1DIsoChorSolver::GetVariableIndex( const char *name )
{
	for ( int i = 0; i < fNOfEquations; ++ i ) {
		if ( strcmp( name, fVariableNames[i] ) == 0 ) {
			return i;
		}
	}

	return -1;
}

int	TTrans1DIsoChorSolver::GetVariableIndex( const char *name, ConstStringArray array, int len )
{
	for ( int i = 0; i < len; ++ i ) {
		if ( strcmp( name, array[i] ) == 0 ) {
			return i;
		}
	}

	return -1;
}

#undef WRITEGLOBALPRODRATE

void TTrans1DIsoChorSolver::WriteFlameletFile( FILE *fp, char *head, char *tail )
{
	T0DPropertiesPtr	props = GetProperties();
	T0DSpeciesPtr		species = GetSpecies();
	Double				*molarMass = species->GetMolarMass()->vec;
	int					i, k;
	int					nOfSpecies = species->GetNOfSpecies();
	int					nSpeciesIn = species->GetNSpeciesInSystem();
	time_t				theDate;
	char				buffer[80];
	char				**names = species->GetNames();
	Flag				fpOpen = FALSE;
	Double				**Y = fMassFracsWork->mat;
	Double				*temp = fTempWork->vec;
	Double				*pressure = fPressureWork->vec;
	Double				*psi = fPsiWork->vec;
	Double				*r = fRWork->vec;
	Double				theTime = GetCurrentTime();
	char				tmpName[127];

	SetOutSolution();

	if ( !fp ) {
		fpOpen = TRUE;
		fp = GetOutputFile( theTime, head, tail, TFlame::kText );
	}

	// write header
	fprintf( fp, "header\n\n" );

	fprintf( fp, "title = \"transient isochor premixed flame\"\n" );
	fprintf( fp, "mechanism = \"%s\"\n", fInputData->fReactionFile );
	fprintf( fp, "author = \"%s\"\n", GetAuthor() );
	time( &theDate );
	strcpy( buffer, ctime( &theDate ) );
	if ( buffer[strlen(buffer)-1] == '\n' )
		buffer[strlen(buffer)-1] = '\0';
	fprintf( fp, "date = \"%s\"\n\n", buffer );
	for ( i = 0; i < GetNFuels(); ++i ) {
		fprintf( fp, "fuel = \"%s\"\n", fVariableNames[fFirstSpecies+GetFuelIndex( i )] );
	}
	fprintf( fp, "time = %g [ms]\n", theTime * 1.0e3 );
	fprintf( fp, "pressure = %g [bar]\n", GetPressure() / 1.0e5 );
	fprintf( fp, "Phi = %g [1/s]\n", GetPhi() );


	fprintf( fp, "Tmax = %g [K]\n", temp[LocationOfMax( fNGridPoints+2
												, &temp[kPrev] )-1] );

	Double	EIFuel = 0.0;
	for ( i = 0; i < GetNFuels(); ++i ) {
		EIFuel += ComputeEmissionIndex( GetFuelIndex( i ) );
	}
	fprintf( fp, "EmissionIndexFuel = %g [kg/m^3s]\n", EIFuel );
	int	index = species->FindSpecies( "NO" );
	if ( index > -1 ) {
		fprintf( fp, "EmissionIndexNO = %g [kg/sm^3]\n"
				, ComputeEmissionIndex( index ) / EIFuel );
	}
	
	index = species->FindSpecies( "NO2" );
	if ( index > -1 ) {
		fprintf( fp, "EmissionIndexNO2 = %g [kg/sm^3]\n"
				, ComputeEmissionIndex( index ) / EIFuel );
	}
	
	index = species->FindSpecies( "N2O" );
	if ( index > -1 ) {
		fprintf( fp, "EmissionIndexN2O = %g [kg/sm^3]\n"
				, ComputeEmissionIndex( index ) );
	}

	fprintf( fp, "CurrTimeStep = %g [s]\n", fRWrk[0]->vec[kF7] );

	fprintf( fp, "\nFuelSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", fTempRightEnd );
	for ( i = 0; i < nOfSpecies; ++i ) {
		if ( fabs( Y[fNGridPoints][i] ) > 1.0e-10 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], Y[fNGridPoints][i] );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "OxidizerSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", fTempLeftEnd );
	for ( i = 0; i < nOfSpecies; ++i ) {
		if ( fabs( Y[kPrev][i] ) > 1.0e-10 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], Y[0][i] );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "numOfSpecies = %d\n", nOfSpecies );
	fprintf( fp, "gridPoints = %d\n\n", fNGridPoints+2 );

	fprintf( fp, "body\n" );

	
// write solution
	PrintFlameletVector( fNGridPoints+2, &r[kPrev], "r", fp );
	PrintFlameletVector( fNGridPoints+2, &psi[kPrev], "psi", fp );
	PrintFlameletVector( fNGridPoints+2, &pressure[kPrev], "pressure", fp );
	PrintFlameletVector( fNGridPoints+2, &temp[kPrev], "temperature [K]", fp );

	// write massfractions of species
	for ( i = 0; i < nOfSpecies; ++i ) {
		sprintf( tmpName, "massfraction-%s", names[i] );
		PrintFlameletVector( fNGridPoints+2, &Y[kPrev][i], tmpName
				, fp, fMassFracsWork->phys_rows );
	}

	if ( fPrintMolarFractions ) {

		MatrixPtr	moleMat = NewMatrix( nOfSpecies, fNGridPoints+2, kColumnPointers );
		Double		**molefracs = &moleMat->mat[kNext];
		for ( k = -1; k <= fNGridPoints; ++k ) {
			fProperties->ComputeMixtureMolarMass( props->GetMixMolarMassRef()
					, Y[k], molarMass, nSpeciesIn );
			for ( i = 0; i < nOfSpecies; ++i ) {
				molefracs[k][i] = Y[k][i] * props->GetMixMolarMass() / molarMass[i];
			}
		}
		for ( i = 0; i < nOfSpecies; ++i ) {
			sprintf( tmpName, "molarfraction-%s", names[i] );
			PrintFlameletVector( fNGridPoints+2, &molefracs[kPrev][i], tmpName
					, fp, moleMat->phys_rows );
		}
		DisposeMatrix( moleMat );
	}
	
// density
	Double	*density = New1DArray( fNGridPoints+2 );
	for ( k = -1; k <= fNGridPoints; ++k ) {
		fProperties->ComputeMixtureMolarMass( props->GetMixMolarMassRef()
				, Y[k], molarMass, nSpeciesIn );
		density[k+1] = GetPressure() * props->GetMixMolarMass() / ( RGAS * temp[k] );
	}
	PrintFlameletVector( fNGridPoints+2, density, "density [kg/m^3]", fp );

	Free1DArray( density );
	
	// Production Rate of NO
	index = species->FindSpecies( "NO" );
	if ( index > -1 ) {
		PrintProdRate( index, fp );
	}
	
#ifdef WRITEGLOBALPRODRATE
	PrintProdRateGlobalReac( theTime );	
#endif
	
	fprintf( fp, "trailer\n" );
	
	Double	*Le = species->GetLewisNumber()->vec;
	for ( i = 0; i < nOfSpecies; ++i ) {
		fprintf( fp, "%s\t%g\n", names[i], Le[i] );
	}
	if ( fpOpen ) {
		fclose( fp );
	}
}

void TTrans1DIsoChorSolver::PrintProdRate( int speciesIndex, FILE *fp )
{
	int			k;
	Double		**Y = fSolMassFracs->mat;
	Double		*T = fSolTemp->vec;
	Double		*prodRate = fSpecies->GetProductionRate()->vec;
	Double		pressure = GetPressure();
	Double		*prodOfR = New1DArray( fNGridPoints+2 );
	char		buffer[128];
	sprintf( buffer, "ProdRate-%s", fSpecies->GetNames()[speciesIndex] );

	for ( k = 0; k < fNGridPoints; ++k ) {
		UpdateThermoProps( Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, NULL );
		prodOfR[k] = prodRate[speciesIndex];
	}
	PrintFlameletVector( fNGridPoints+2, prodOfR, buffer, fp );
	
	Free1DArray( prodOfR );
}

Double TTrans1DIsoChorSolver::ComputeEmissionIndex( int speciesIndex )
{
	int			k;
	Double		**Y = fSolMassFracs->mat;
	Double		*T = fSolTemp->vec;
	Double		*r = fSolR->vec;
	Double		*prodRate = fSpecies->GetProductionRate()->vec;
	Double		pressure = GetPressure();
	char		**names = fSpecies->GetNames();
	Double		sum;

	sum = 0.0;  // m[0] = 0, m[L] = 0
	for ( k = 1; k < fNGridPoints; ++k ) {
		UpdateThermoProps( Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, NULL );
		sum += prodRate[speciesIndex] * ( r[k+1] - r[k-1] );
	}			
	return 0.5 * sum;
}

FILE *TTrans1DIsoChorSolver::GetOutputFile( Double theTime, const char *head, const char *tail, FileType type )
{
	int 		nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();
	char		*name = new char[64];
	FILE 		*fp;
	
	sprintf( name, "%s%s%.8s_p%.2dt%8.2ems%s"
					, ( head ) ? head : "", ( head ) ? "_" : ""
					, fSpecies->GetNames()[GetFuelIndex()]
					, ( int ) floor( GetPressure() * 1.0e-5 + 0.5 )	// in [bar]
					, theTime * 1000.0
					, ( tail ) ? tail : "" );
	
	fp = GetOutfile( name, type );
	delete name;
	
	return fp;
}

FILE *TTrans1DIsoChorSolver::GetOutputFile( char */*head*/, char */*tail*/, FileType /*type*/ )
{ 
	fprintf( fOutFilePtr, "wrong instance of function 'TTrans1DIsoChorSolver::GetOutputFile'\n" );
	exit( 2 );

	return NULL;
}


void TTrans1DIsoChorSolver::DoExit( void )
{
	FILE	*fp = GetOutfile( "interrupt", kText );
	WriteFlameletFile( fp, NULL, NULL );
	fprintf( fOutFilePtr, "\nprogram stopped by user\noutput written to %s\n"
				, GetOutfileName( "interrupt", kText ) );
	fclose( fp );
	exit( 2 );
}

Double TTrans1DIsoChorSolver::GetCurrentTime( void )
{
	int	theTime = GetActualPoint( fTEnd );

	return ( theTime >= 0.0 ) ? fSolTime->vec[theTime] : fTEnd;
}


#ifdef NOCCLINK
extern "C" int _main();
#endif

#ifdef WRITEGLOBALPRODRATE
void TTrans1DIsoChorSolver::PrintProdRateGlobalReac( Double time )
{
#include "iOctaneNOx.shortnoign.h"

	int			k;
	Double		**Y = fSolMassFracs->mat;
	Double		*T = fSolTemp->vec;
	Double		*r = fSolR->vec;
	Double		pressure = GetPressure();
	Double		*w = GetReaction()->GetReactionRate()->vec;
	Double		reacRate;
	FILE		*fp = GetOutputFile( time, "GlRa", NULL, TFlame::kData );
	
//	header
	fprintf( fp, "*\n" );
	fprintf( fp, "%s", "y" );
	fprintf( fp, "\t%s", "I.    N\\d2\\n + O\\d2\\n = 2 NO" );
	fprintf( fp, "\t%s", "IIa.  N\\d2\\n + C\\d2\\nH\\d4\\n + O\\d2\\n + OH = NO + HCN + H\\d2\\nO + CO H\\d2\\n" );
	fprintf( fp, "\t%s", "IIb.  HCN + O\\d2\\n + H\\d2\\nO = NO + CO + H\\d2\\n + OH" );
	fprintf( fp, "\t%s", "IIIa. N\\d2\\n + 2 OH + M = N\\d2\\nO + H\\d2\\nO + M" );
	fprintf( fp, "\t%s", "IIIb. N\\d2\\nO + 2OH = 2 NO + H\\d2\\nO" );
	fprintf( fp, "\t%s", "IV.   NO + H\\d2\\nO = NO\\d2\\n + H\\d2\\n" );
	fprintf( fp, "\t%s", "Va.   NO + C\\d2\\nH\\d4\\n + OH = HNCO + CO + 2H\\d2\\n" );
	fprintf( fp, "\t%s", "Vb.   HNCO + NO = N\\d2\\n + CO + OH" );
	fprintf( fp, "\t%s", "VI.   NH\\d3\\n + NO = N\\d2\\n + H\\d2\\n + OH" );

	for ( k = 0; k < fNGridPoints; ++k ) {
		UpdateThermoProps( Y[k], T[k], pressure, fProperties->GetDensityRef()
										, kDensFromPress, NULL );
		fprintf( fp, "\n%-12E", r[k] );
		// I
		reacRate = w[rN7f] - w[rN7b] - w[rN101f] + w[rN101b] - w[rN104f] + w[rN104b] 
				 - w[rN37f] + w[rN37b];
		fprintf( fp, "\t%-12E", reacRate );
		// IIa
		reacRate = w[rN110f] - w[rN110b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// IIb
		reacRate = - w[rN111] - w[rN113f] + w[rN113b] + w[rN120f] - w[rN120b]
			 	   - w[rN130f] + w[rN130b] - w[rN131f] + w[rN131b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// IIIa
		reacRate = w[rN6f] - w[rN6b] + w[rN9f] - w[rN9b] - w[rN25f] + w[rN25b]
				 + w[rN40f] - w[rN40b] - w[rN45f] + w[rN45b] - w[rN46yf] + w[rN46yb]
				 - w[rN49f] + w[rN49b] + w[rN60f] - w[rN60b] + w[rN63f] - w[rN63b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// IIIb
		reacRate = w[rN17f] - w[rN17b] + w[rN18f] - w[rN18b] + w[rN19f] - w[rN19b]
				 + w[rN20f] - w[rN20b] + w[rN47f] - w[rN47b] - w[rN71f] + w[rN71b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// IV
		reacRate = - w[rN84f] + w[rN84b] + w[rN85f] - w[rN85b] - w[rN87f] + w[rN87b]
				   - w[rN88f] + w[rN88b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// Va
		reacRate = - w[rN62f] + w[rN62b] - w[rN65f] + w[rN65b] + w[rN74f] - w[rN74b]
				   + w[rN112f] - w[rN112b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// Vb
		reacRate = w[rN60f] - w[rN60b] + w[rN61f] - w[rN61b] + w[rN63f] - w[rN63b];
		fprintf( fp, "\t%-12E", reacRate );
		
		// VI
		reacRate = w[rN1f] - w[rN1b] + w[rN2f] - w[rN2b] + w[rN3f] - w[rN3b] + 
				   w[rN4f] - w[rN4b];
		fprintf( fp, "\t%-12E", reacRate );
		
	}
	
	fclose( fp );
}
#endif

void TTrans1DIsoChorSolver::SetInitialPsi( void )
{
	int			nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double		**sol = fSolution->mat;
	Double		**nY = fSolMassFracs->mat;
	Double		rho;
	
	copy( nSpeciesIn, &sol[kPrev][fFirstSpecies], 1, nY[kPrev], 1 );
	ComputeProperties( sol[kPrev][fTemperature], nY[kPrev], sol[kPrev][fPressure], rho, kDensFromPress, NULL );
	sol[kPrev][fPsi] = rho * modPow( fRLeft, fGamma + 1 ) / ( fGamma + 1 );
	for ( int k = 0; k < fNGridPoints; ++k ) {
		copy( nSpeciesIn, &sol[k][fFirstSpecies], 1, nY[k], 1 );
		ComputeProperties( sol[k][fTemperature], nY[k], sol[k][fPressure], rho, kDensFromPress, NULL );
		sol[k][fPsi] = sol[k-1][fPsi] + rho * modPow( sol[k][fR], fGamma )
						   * ( sol[k][fR] - sol[k-1][fR] );
	}
	sol[fNGridPoints][fPsi] = sol[fNGridPoints-1][fPsi] + rho * modPow( sol[fNGridPoints][fR], fGamma )
						   * ( sol[fNGridPoints][fR] - sol[fNGridPoints-1][fR] );
}

void TTrans1DIsoChorSolver::ComputeGPDistributionWeights( Double *sol )
{
	int			nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double		tempMax = 3000.0; 					// estimated
	Double		cY = fPsiRight * fPsiRight;
	Double		cT = cY / ( tempMax * tempMax );
	Double		dTdPsi; 
	int			gpOff;
	
	for ( int k = 0; k < fNGridPoints-1; ++k ) {
		gpOff = k * fNOfEquations;
		fGPDistWeight[k] = 1.0;
		dTdPsi = ( sol[gpOff+fTemperature+fNOfEquations] - sol[gpOff+fTemperature] ) / fh[k];
		fGPDistWeight[k] += cT * dTdPsi * dTdPsi;
//		for ( int j = 0; j < nSpeciesIn; ++j ) {
//			dYdPsi = ( sol[gpOff+fNOfEquations+fFirstSpecies+j] - sol[gpOff+fFirstSpecies+j] ) / fh[k];
//			fGPDistWeight[k] += cY * dYdPsi * dYdPsi;
//		}
		fGPDistWeight[k] = sqrt( fGPDistWeight[k] );
	}
	fGPDistWeight[fNGridPoints-1] = 1.0;
}

Double TTrans1DIsoChorSolver::GetTempSource( Double r, Double time )
{
	Double rs = 1.0e-3;
	Double ti = 5.0e-5;
	Double D = 3.0e6;	// units are J/m
	Double dummy = pow( r / rs, 8.0 );
	
	dummy = ( dummy > 100.0 ) ? 100.0 : dummy;
	if ( time > 0.0 && time <= ti ) {
		return D / ti * exp( - dummy );
	}
	else {
		return 0.0;
	}
}
