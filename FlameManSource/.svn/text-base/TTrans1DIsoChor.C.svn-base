#include "FlameMaster.h"
#include "TTrans1DIsoChor.h"
//#include "ListTool.h"
#include "Spline.h"

#undef MAPPING

void TTrans1DIsoChor::InitTTrans1DIsoChor( void )
{
	int i;
	fNOutputs = fInputData->fNOutputs;
	fNGridPoints = fInputData->fInitialGridPoints-2;
	fEquidistant = fInputData->fEquidistant;
	fTStart = fInputData->fStart;
	fTEnd = fInputData->fEnd;
	fDeltaTStart = fInputData->fDeltaTStart;
	fPressStart = 0.0;
	fTimeCut = 1.3e-3;
//	fTimeCut = 1.0e-5;
	
	if ( fInputData->fPhi ) {
		SetPhi( fInputData->fPhi->vec[0] );
		if ( fInputData->fPhi->len > 1 ) {
			cerr << "#warning: more than one equivalence ratio specified, use" 
					<< GetPhi() << NEWL;
		}
		if ( GetPhi() < 0.0 ) {
			cerr << "initial equivalence ratio missing, use given massfractions " << NEWL;
		}
		else {
			cerr << "initial equivalence ratio is " << GetPhi() << NEWL;
		}
	}
	else {
		cerr << "#error: no equivalence ratio specified" << NEWL;
		exit( 2 );
	}
	
	fGridSol = NewVector( fNGridPoints+2 );
	fStartSol = NewMatrix( fVariables, fNGridPoints+2, kColumnPointers );
	
	fGridSol->vec = &fGridSol->vec[kNext];
	fStartSol->mat = &fStartSol->mat[kNext];

	fBCFlagLeft = New1DIntArray( fNOfEquations );
	fBCFlagRight = New1DIntArray( fNOfEquations );

	
	int	nOfSpeciesIn = fSpecies->GetNSpeciesInSystem();
	fVariableNames = new String[fVariablesWithoutSpecies + nOfSpeciesIn];
	if ( !fVariableNames ) {
		cerr << "#error memory allocation of TTrans1DIsoChor failed" << NEWL;
		exit( 2 );
	}
	fVariableNames[fR] = new char[2];
	strcpy( fVariableNames[fR], "R" );
	fVariableNames[fPsi] = new char[4];
	strcpy( fVariableNames[fPsi], "psi" );
	fVariableNames[fPressure] = new char[9];
	strcpy( fVariableNames[fPressure], "Pressure" );
	fVariableNames[fTemperature] = new char[2];
	strcpy( fVariableNames[fTemperature], "T" );
	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies->GetNames()[i] ) + 1];
		strcpy( fVariableNames[fFirstSpecies + i], fSpecies->GetNames()[i] );
	}

	MakeGrid( &fGridSol->vec[kPrev], fNGridPoints+2, fInputData->fLeft, fInputData->fRight, fEquidistant );
	SetInitialBC();

	ReadStartProfiles( fInputData );

	fTempLeftEnd = fStartSol->mat[kPrev][fTemperature];
	fTempRightEnd = fStartSol->mat[fNGridPoints][fTemperature];
	fRLeftEnd = fStartSol->mat[kPrev][fR];
	fRRightEnd = fStartSol->mat[fNGridPoints][fR];

	fTCurr = fTStart;
		
	fDeltaT = ( fTEnd - fTStart ) / ( ( Double ) fNOutputs );
	if ( fDeltaT <= 0.0 ) {
		cerr << "#error: t_start = " << fTStart << " greater as t_end = " 
					<< fTEnd << NEWL;
		exit( 2 );
	}

	ffpEI = NULL;

	char	buffer[32];
	sprintf( buffer, "Scalars_phi%d", ( int ) GetPhi() );
	ffpEI = GetOutfile( buffer, kData );
	fprintf( ffpEI, "*\ntime [s]" );
	for ( i = 0; i < GetNFuels(); ++i ) {
		fprintf( ffpEI, "\tEI_%s [kg/sm^3]", fSpecies->GetNames()[GetFuelIndex( i )] );
	}
	if ( fSpecies->FindSpecies( "NO" ) > -1 ) {
		fprintf( ffpEI, "\tEI_NO [kg/sm^3]" );
	}
	if ( fSpecies->FindSpecies( "NO2" ) > -1 ) {
		fprintf( ffpEI, "\tEI_NO2 [kg/sm^3]" );
	}
	if ( fSpecies->FindSpecies( "N2O" ) > -1 ) {
		fprintf( ffpEI, "\tEI_N2O [kg/sm^3]" );
	}
	fprintf( ffpEI, "\n" );
	fflush( ffpEI );

	FILE	*fp = GetOutfile( "InitialValues", kData );
	PrintSolution( fp, fNGridPoints, fNOfEquations
					, fStartSol->mat, fVariableNames );
	fclose( fp );

#ifdef HP
	InstallInterruptHP();
#endif
}

TTrans1DIsoChor::~TTrans1DIsoChor( void )
{
	Free1DIntArray( fBCFlagRight );
	Free1DIntArray( fBCFlagLeft );

	fStartSol->mat = &fStartSol->mat[kPrev];
	fGridSol->vec = &fGridSol->vec[kPrev];

	DisposeMatrix( fStartSol );
	DisposeVector( fGridSol );
	
	if ( ffpEI ) {
		fclose( ffpEI );
	}
}

void TTrans1DIsoChor::Solve( void )
{
	FILE	*fp;
	Flag 	leave = FALSE;
	Double	bound = fTEnd - 1.0e-20;
	Double	tempEnd;
	
	Initialize( fTStart, fVariableNames, &fStartSol->mat[kPrev], fNOfEquations
					, &fGridSol->vec[kPrev], fNGridPoints+2
					, TFlame::GetPressure(), GetPhi(), fDeltaTStart );

	cerr << "dump output" << NEWL;
	fp = GetOutputFile( fTCurr, NULL, NULL, TFlame::kText );
	WriteFlameletFile( fp, NULL, NULL );
	fclose( fp );

#ifdef MAPPING
	cerr << NEWL << "stop here for mapping" << NEWL << NEWL;
	exit(2);
#endif

	cerr << NEWL << "start computation" << NEWL << NEWL;

	int	jCountOutputs = 0;
	do {
		fTCurr += fDeltaT;
		++jCountOutputs;
		tempEnd = fTempRightEnd;
		leave = TTrans1DIsoChorSolver::Solve( fTCurr, TFlame::GetPressure(), GetPhi(), fTempLeftEnd, tempEnd, -1 );
		
		if ( jCountOutputs % 1 == 0 ) {
			cerr << "dump output" << NEWL;
			
			fp = GetOutputFile( fTCurr, ( leave ) ? "Sol" : NULL, NULL, TFlame::kText );
			WriteFlameletFile( fp, NULL, NULL );
			fclose( fp );
		}
		WriteScalars( fTCurr );

	} while ( fTCurr < bound && !leave );

	cerr << "done" << NEWL;

	return;
}

Double TTrans1DIsoChor::GetTempEnd( Double theTime, Double timeStart, Double timeCut
									, Double tempStart, Double tempEnd )
{
	if ( theTime < timeCut ) {
		return Interpol( theTime, tempStart, timeStart, tempEnd, timeCut );
	}
	else {
		return tempEnd;
	}
}

void TTrans1DIsoChor::SetInitialBC( void )
{
	int					i;
	Double				mixMolarMass;
	int					nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	SpeciesPtr			species = fInputData->GetSpecies();
	BoundaryInputPtr	right = fInputData->leftBoundary;
	BoundaryInputPtr	left = fInputData->rightBoundary;
	int					inpTOffset = fInputData->fTemperatureOffset;
	int					*speciesIndexLeft = NULL;
	int					*speciesIndexRight = NULL;
	int					leftSpecifiedSpecies = left->fSpecifiedSpeciesBCs;
	int					rightSpecifiedSpecies = right->fSpecifiedSpeciesBCs;
	Double				*yleft = fStartSol->mat[kPrev];
	Double				*yright = fStartSol->mat[fNGridPoints];
	
	speciesIndexLeft = new int[left->fSpecifiedSpeciesBCs];
	if ( !speciesIndexLeft ) FatalError( "memory allocation of T1DIsoChor failed" );
	speciesIndexRight = new int[right->fSpecifiedSpeciesBCs];
	if ( !speciesIndexRight ) FatalError( "memory allocation of T1DIsoChor failed" );

	//	set speciesIndex
	for ( i = 0; i < leftSpecifiedSpecies; ++i ) {
		speciesIndexLeft[i] = fInputData->FindSpecies( left->speciesName[i] );
	}
	for ( i = 0; i < rightSpecifiedSpecies; ++i ) {
		speciesIndexRight[i] = fInputData->FindSpecies( right->speciesName[i] );
	}
	
	// set BCFlags
	fBCFlagLeft[fTemperature] = left->fBcFlag[inpTOffset];
	fBCFlagRight[fTemperature] = right->fBcFlag[inpTOffset];
	
	if ( left->fBcFlagSpecies == kNeumann || right->fBcFlagSpecies == kNeumann ) {
		cerr << "#error: invalid boundary condition type "
			<< "'Neumann' for species concentration" << NEWL;			
	}
	for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
		fBCFlagLeft[i] = left->fBcFlagSpecies;
		fBCFlagRight[i] = right->fBcFlagSpecies;
	}

	// set value
	yleft[fTemperature] = left->fValue[inpTOffset];
	yright[fTemperature] = right->fValue[inpTOffset];
	yleft[fPressure] = fInputData->fPressure->vec[0];
	yright[fPressure] = fInputData->fPressure->vec[0];
	yleft[fR] = fInputData->fLeft;
	yright[fR] = fInputData->fRight;
	
	for ( i = 0; i < leftSpecifiedSpecies; ++i ) {
		yleft[speciesIndexLeft[i]+fFirstSpecies] = left->fValueSpecies[i];
	}
	if ( left->fMixtureSpecification == kMolarFraction ) {
		// first compute molar mass of mixture
		for ( i = 0, mixMolarMass = 0; i < nSpeciesInSystem; ++i ) {
			mixMolarMass += species[i].molarMass * yleft[i+fFirstSpecies];
		}
		// compute massfractions
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			yleft[i+fFirstSpecies] *= species[i].molarMass / mixMolarMass;
		}
		for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
			fBCFlagLeft[i] = kMassFraction;
		}
	}

	for ( i = 0; i < rightSpecifiedSpecies; ++i ) {
		yright[speciesIndexRight[i]+fFirstSpecies] = right->fValueSpecies[i];
	}
	if ( right->fMixtureSpecification == kMolarFraction ) {
		// first compute molar mass of mixture
		for ( i = 0, mixMolarMass = 0; i < nSpeciesInSystem; ++i ) {
			mixMolarMass += species[i].molarMass * yright[i+fFirstSpecies];
		}
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			yright[i+fFirstSpecies] *= species[i].molarMass / mixMolarMass;
		}
		for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
			fBCFlagRight[i] = kMassFraction;
		}
	}
	delete speciesIndexRight;
	delete speciesIndexLeft;
}

void TTrans1DIsoChor::SetInitialValues( void )
{
	Double	**y = fStartSol->mat;
	Double	*x = fGridSol->vec;
	Double	*yLeft = fStartSol->mat[kPrev];
	Double	*yRight = fStartSol->mat[fNGridPoints];
	Double	left = fGridSol->vec[kPrev];
	Double	right = fGridSol->vec[fNGridPoints];

	for ( int k = 0; k < fNGridPoints; ++k ) {
		for ( int j = 0; j < fNOfEquations; ++j ) {
			y[k][j] = Interpol( x[k], yLeft[j], left, yRight[j], right );
		}
	}
}

void TTrans1DIsoChor::ReadStartProfiles( TInputDataPtr inp )
{
	StartProfilePtr	sp = NULL;
	char 	*insp = inp->fStartProfileFile;
	FILE	*fpS = NULL;
	char	*fileName = GetFullPath( insp, kFileName );

	sp = new StartProfile;	

	if ( !insp || ( fpS = fopen( fileName, "r" ) ) == NULL ) {
		fprintf( stderr, "interpolate initial solution from boundary conditions\n" );
		SetInitialValues();
	} 
	else {
		if ( ( fpS = fopen( fileName, "r" ) ) == NULL ) {
			fprintf( stderr, "startprofiles file '%s' not found\n", fileName );
			exit( 2 );
		} 
		else {
			fprintf( stderr, "read initial solution from file '%s'\n", fileName );
			::ReadStartProfiles( sp, fpS );
			SetInitialValues( fInputData, sp );
			CleanReadStartProfile();
			delete sp;
			fclose( fpS );
		}
	}
	delete fileName;
	
}

void TTrans1DIsoChor::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp )
{
	int 				i, j, k;
	int					variable, speciesIndex;
	Flag				XSet = FALSE;
	Flag				chooseInputGrid = FALSE;
	int					gridPointsIn = sp->gridPoints;	// including boundaries
	int					nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	Double				*xIn = New1DArray( gridPointsIn );
	Double				*yIn =  New1DArray( gridPointsIn );
	Double				*yInFloat = sp->data;
	char				*string = sp->labels;
	Double				**y = fStartSol->mat;
	Double				*yWork = New1DArray( fNGridPoints+2 );
	yWork = &yWork[kNext];
	SplinePtr			theSpline = NULL;
	Double				leftSlope;
	Double				rightSlope;
	Double				*x = fGridSol->vec;
	struct _parameter	*param = GetParameter( "pressure" );
	
	SetInitialValues(); // default values

//	choose grid from input or equidistant
	if ( gridPointsIn >= fNGridPoints+2 ) {
		fNGridPoints = gridPointsIn - 2;
		chooseInputGrid = TRUE;
	}
	else {
		chooseInputGrid = FALSE;
	}
	
// find independent coordinate
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strncmp( string, "r", 1 ) == 0 && XSet == FALSE ) {
			if ( chooseInputGrid ) {
				cerr << "choose inputGrid" << NEWL;
				for ( j = -1; j <= gridPointsIn-2; ++j ) {
					x[j] = yInFloat[i*gridPointsIn + j+1];		// implicit cast from float to Double
				}
				XSet = TRUE;
			}
			else { // choose own Grid, but read x for interpolation
				cerr << "choose own Grid" << NEWL;
				for ( j = 0; j < gridPointsIn; ++j ) {
					xIn[j] = yInFloat[i*gridPointsIn + j];		// implicit cast from float to Double
				}
				XSet = TRUE;
			}
		}
		string += strlen(string) + 1;
	}

// error checking
	if ( !XSet ) {
		cerr << "error: can't find coordinate 'x'" << NEWL;
		exit(2);
	}

// reset string
	string = sp->labels;
	
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strncmp( string, "r", 2 ) == 0 ) {
			variable = fR;
		}
		if ( strncmp( string, "psi", 4 ) == 0 ) {
			variable = fPsi;
		}
		if ( strncmp( string, "pressure", 9 ) == 0 ) {
			variable = fPressure;
		}
		if ( strncmp( string, "temperature", 11 ) == 0 ) {
			variable = fTemperature;
		}
		else if ( strncmp( string, "massfraction-", 13 ) == 0 ){
			string += 13;
			UpperString( string );
			if ( ( speciesIndex = inp->FindSpecies( string ) ) >= 0 ) {
//				if ( speciesIndex < nSpeciesInSystem ) {
					variable = fFirstSpecies + speciesIndex;
//				}
//				else {
//					string += strlen(string) + 1;
//					continue;
//				}
			}
			else {
				cerr << "warning: no match for species " << string << NEWL;
				string += strlen(string) + 1;
				continue;
			}
		}
		else {
			string += strlen(string) + 1;
			continue;
		}

		string += strlen(string) + 1;
		if ( chooseInputGrid ) {
			for ( k = -1; k <= gridPointsIn-2; ++k ) {
				y[k][variable] = yInFloat[i*gridPointsIn + k + 1];	// copy workspace to vector of solution
			}
		}
		else {
			for ( k = 0; k < gridPointsIn; ++k ) {	// store vector in workspace
				yIn[k] = yInFloat[i * gridPointsIn + k];	// implicit cast from float to Double
			}
		

#ifdef MAPPING
			yIn[0] = yIn[1];
			yIn[gridPointsIn-1] = yIn[gridPointsIn-2];
#endif

			leftSlope = ( yIn[1] - yIn[0] ) / ( xIn[1] - xIn[0] );
			rightSlope = ( yIn[gridPointsIn-1] - yIn[gridPointsIn-2] ) / ( xIn[gridPointsIn-1] - xIn[gridPointsIn-2] );
			theSpline = ComputeSimpleSpline( xIn, yIn, gridPointsIn, FALSE, leftSlope, FALSE, rightSlope, NULL, TRUE );
			SplineInterpolate( theSpline, &x[kPrev], &yWork[kPrev], fNGridPoints+2 );
			for ( k = -1; k <= fNGridPoints; ++k ) {
				y[k][variable] = yWork[k];	// copy workspace to vector of solution
			}
		}
	}
	
	FreeSpline( theSpline );
	yWork = &yWork[kPrev];
	Free1DArray( yWork );
	Free1DArray( yIn );
	Free1DArray( xIn );
	
// set scalars
	if ( param ) {
		fPressStart = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "bar" ) == 0 ) {
			fPressStart *= 1.0e5;
		}
	}
	else { // choose default
		cerr << "#warning: no input value for 'pressure', use " << fPressStart << NEWL;
	}
	
	param = GetParameter( "time" );
	if ( param ) {
		fTStart = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "ms" ) == 0 ) {
			fTStart *= 1.0e-3;
		}
	}
	else { // choose default
		cerr << "#warning: no input value for 'time', use " << fTStart << NEWL;
	}

	param = GetParameter( "currtimestep" );
	if ( param ) {
		fDeltaTStart = (Double)param->what.quantity.value;
		if ( strcmp( param->what.quantity.unit, "ms" ) == 0 ) {
			fDeltaTStart *= 1.0e-3;
		}
	}
	else { // choose default
		cerr << "#warning: no input value for 'CurrTimeStep', use default" << NEWL;
	}
	
	fTempLeftStart = y[kPrev][fTemperature];
	fTempRightStart = y[fNGridPoints][fTemperature];

	FILE	*fp = GetOutfile( "initialguess", TFlame::kData );
	PrintSolution( fp, fNGridPoints, fNOfEquations, y, fVariableNames );
	fclose(fp);
}

void TTrans1DIsoChor::WriteScalars( Double time )
{
	int		index;

	fprintf( ffpEI, "%g", time );

	Double	EIFuel = 0.0;
	for ( int i = 0; i < GetNFuels(); ++i ) {
		EIFuel = ComputeEmissionIndex( GetFuelIndex( i ) );
		fprintf( ffpEI, "\t%g", EIFuel );
	}
	
	index = fSpecies->FindSpecies( "NO" );
	if ( index > -1 ) {
		fprintf( ffpEI, "\t%g", ComputeEmissionIndex( index ) );
	}
	
	index = fSpecies->FindSpecies( "NO2" );
	if ( index > -1 ) {
		fprintf( ffpEI, "\t%g", ComputeEmissionIndex( index ) );
	}
	
	index = fSpecies->FindSpecies( "N2O" );
	if ( index > -1 ) {
		fprintf( ffpEI, "\t%g", ComputeEmissionIndex( index ) );
	}

	fprintf( ffpEI, "\n" );
	fflush( ffpEI );
}
