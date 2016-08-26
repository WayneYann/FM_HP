#include "FlameMaster.h"
#include "ListTool.h"
#include "Spline.h"
#include "TCountPremFlameSim.h"


#define FCONVLEFTTORIGHT

#define UPWINDCONVECTION

#undef FULLDIFFUSION

#define MOLARDIFF

void TCountPremFlameSim::InitCountPremFlameSim( void )
{
	int i;
	TBVPSolverPtr	solver = GetSolver();
	TNewtonPtr		bt = solver->bt;
	TGridPtr		fine = bt->GetGrid()->GetFine();
	TGridPtr		coarse = bt->GetGrid()->GetCoarse();
	int				nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	int				maxGridPoints = bt->GetMaxGridPoints();

#ifdef FULLDIFFUSION
	flame->UnSetDiffCorr();
#endif

	fDeltaArcLength = 0.0;
	fdTds = 0.0;
	fdlnStrainrateds = 0.0;
	fDeltaTRef = 20.0;
	fNFlameletsCount = 0;
	fMaxFlamelets = 100;
	fDeltaStrainrateRef = log( 1.3 ); // will give a/a_old = 1.3
	fArcLengthContin = fInputData->fArcLengthContin;
	fStrainRateContin = fInputData->fStrainRateContin;

	if ( fSoot ) {
		fSoot->SetMomentsOffset( fSootMoments );
	}

	if ( fInputData->fParameterComm >= 0.0 ) {
		if ( fStrainRate ){
			DisposeVector( fStrainRate );
		}
		fStrainRate = NewVector( 1 );
		fStrainRate->vec[0] = fInputData->fParameterComm;
		fStrainRate->len = 0;
		cerr << "use initial strain rate from command line: a = " << fStrainRate->vec[0] << NEWL;
	}

	fPremConfiguration = fInputData->fPremConfiguration; // options are kBackToBack, kFreshToBurned, kFull

	fprintf( stderr, "Configuration is %s\n"
		, (fPremConfiguration == kBackToBack) ? "back to back" : ((fPremConfiguration == kBackToBack) ? "fresh to burned" : "full premixed twin flame" ));

	fVariableNames = new String[fVariablesWithoutSpecies + nSpeciesInSystem];

	fVariableNames[fVVelocity] = new char[2];
	strcpy( fVariableNames[fVVelocity], "f" );
	fVariableNames[fUVelocity] = new char[3];
	strcpy( fVariableNames[fUVelocity], "f'" );
	fVariableNames[fTemperature] = new char[2];
	strcpy( fVariableNames[fTemperature], "T" );
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies->GetNames()[i] ) + 1];
		strcpy( fVariableNames[fFirstSpecies + i], fSpecies->GetNames()[i] );
	}
	fVariableNames[fLnStrainrate] = new char[4];
	strcpy( fVariableNames[fLnStrainrate], "lnA" );
	if ( fSoot ) {
		int	offset = fSoot->GetOffsetSootMoments();
		for ( i = 0; i < fSoot->GetNSootMoments(); ++i ) {
			fVariableNames[offset + i] = new char[8];
			sprintf( fVariableNames[offset + i], "M%d/rho", i );
		}
	}
	fPrintMolarFractions = fInputData->fPrintMolarFractions;
	
//	always solve for some extra equation now
	fSolLnStrainrate = NewVector( bt->GetMaxGridPoints() + 2 );
	fSolLnStrainrate->vec = &fSolLnStrainrate->vec[kNext];
	fSolLnStrainrate->len -= 2;
	fSavedLnStrainrate = NewVector( bt->GetMaxGridPoints() + 2 );
	fSavedLnStrainrate->vec = &fSavedLnStrainrate->vec[kNext];
	fSavedLnStrainrate->len -= 2;

//	vectors of solution
	fSolV = NewVector( maxGridPoints + 2 );
	fSolU = NewVector( maxGridPoints + 2 );

	fSolV->vec = &fSolV->vec[kNext];
	fSolU->vec = &fSolU->vec[kNext];

	fSolV->len -= 2;
	fSolU->len -= 2;

//	saved solution
	fSavedV = NewVector( maxGridPoints + 2 );
	fSavedU = NewVector( maxGridPoints + 2 );

	fSavedV->vec = &fSavedV->vec[kNext];
	fSavedU->vec = &fSavedU->vec[kNext];

	fSavedV->len -= 2;
	fSavedU->len -= 2;

	bt->SetUtFuncs( CountPremSimCheckCompDomain, CountPremSimJacRest, CountPremSimJacRest
					, CountPremSimRHSRest, CountPremSimRHSRest, CountPremSimRHSRest
					, CountPremSimOutput, CountPremSimPostIter
					, SetCountPremSimNodeInfo, CountPremSimPostConv
					, GetCountPremSimVarNames
					, CountPremSimUpdateLeftBoundary, CountPremSimUpdateRightBoundary );
	T1DFlame::SetInitialBC( fine, fInputData );
	T1DFlame::SetInitialBC( coarse, fInputData );
	SetInitialBC( fine, fInputData );
	SetInitialBC( coarse, fInputData );
	cerr << "initial equivalence ratio is " << GetPhi() << NEWL;
	fMassFraction = new TMassFraction( this );

	ReadStartProfiles( fInputData );
	CheckBC();
	CheckInitialGuess();
	UpdateSolution( fine->GetY(), fine->GetYLeft(), fine->GetYRight() );


	if ( GetArcLengthCont() ) {
		if ( fStrainRate && fStrainRate->phys_len > 1 ) {
			if (fStrainRate->vec[1] >= fStrainRate->vec[0] ) {
				fArcUp = TRUE;
			}
			else {
				fArcUp = FALSE;
			}
		}
		else {
			fArcUp = TRUE;
		}
		fTempContStart = fSolTemp->vec[GetTMaxLoc()];
	}
	else {
		if ( fStrainRateContin ) {
			int	inflectionPoint = LocationOfMaxSlope( fSolTemp->vec, fine->GetX()->vec, fine->GetNGridPoints() );
			if ( bt->GetGrid()->IsFine() && inflectionPoint%2 == 0) ++inflectionPoint;	// make it odd for the fine grid ( even gridpoints are removed for the coarse grid! )
			SetTMaxLoc( inflectionPoint );
			fTempContStart = fSolTemp->vec[inflectionPoint];
			fLnStrainrateContStart = fine->GetX()->vec[inflectionPoint];
//			fprintf(stderr, "TPointstart at x = %g and ng = %d = %g\n", fine->GetX()->vec[inflectionPoint]
//											, inflectionPoint, fSolTemp->vec[inflectionPoint]);
		}
	}

	// fLnStrainrateContStart needs to be set even if ArcLengthCont is FALSE
	fLnStrainrateContStart = log( fStrainRate->vec[0] );

	// compute enthalpy at boundary
//	fEnthalpy = 0;
//	Double *hi = fFlameNode->enthalpy;
//	Double *Yi = fSolMassFracs->mat[kPrev];
//
//	SetFlameNode( kPrev );
//	ComputeProperties( fFlameNode, fSolTemp->vec[kPrev], Yi, GetPressure() );
//
//	for (int i = 0; i < nSpeciesInSystem; ++i) {
//		fEnthalpy += Yi[i] * hi[i];
//	}
//	fprintf(stderr, "Enthalpy = %g\n", fEnthalpy );
}

TCountPremFlameSim::~TCountPremFlameSim( void )
{
	int	nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();

	delete fMassFraction;

	fSolLnStrainrate->vec = &fSolLnStrainrate->vec[kPrev];
	DisposeVector( fSolLnStrainrate );
	fSavedLnStrainrate->vec = &fSavedLnStrainrate->vec[kPrev];
	DisposeVector( fSavedLnStrainrate );

	fSavedU->vec = &fSavedU->vec[kPrev];
	fSavedV->vec = &fSavedV->vec[kPrev];

	DisposeVector( fSavedU );
	DisposeVector( fSavedV );

	fSolV->vec = &fSolV->vec[kPrev];
	fSolU->vec = &fSolU->vec[kPrev];

	DisposeVector( fSolU );
	DisposeVector( fSolV );

	for ( int i = 0; i < nSpeciesInSystem+fVariablesWithoutSpecies; ++i ) {
		delete fVariableNames[i];
	}
	delete fVariableNames;
}

void CountPremSimJacRest( void */* object */, NodeInfoPtr /* nodeInfo */ )
{
}

void CountPremSimRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	if ( !flame->RHSAction( nodeInfo, rhsMode ) ) {
		return;
	}
	
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int 	fLnStrainrate = flame->GetOffsetLnStrainrate();
	int		eqLoop, speciesEq;
	int		M = nodeInfo->nOfEquations;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
    Double  hnenn = nodeInfo->hnenn;
	Double	*rhs = nodeInfo->rhs;
	Double	*YPrev = flameNode->Y[kPrev];
	Double	*Y = flameNode->Y[kCurr];
	Double	*YNext = flameNode->Y[kNext];
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
//	Double	strainRate = flame->GetStrainRate();
	Double	strainRate = exp( y[fLnStrainrate] );
	Double	*enthalpy = flameNode->enthalpy;
	Double	*mixDensity = flameNode->mixDensity;
	Double	mixDensityInf = flameNode->rhoInf;
	Double	mixViscosityInf = flameNode->viscosityInf;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*productionRate = flameNode->productionRate;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	oneOverATimesRho = 1.0 / ( strainRate * mixDensity[kCurr] );
	Double	idealGasCoeff = flame->GetPressure() * *flameNode->mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	mixHeatCapacity = *flameNode->mixHeatCapacity;
	Double	constThermDiffCoeff = constMassDiffCoeff / flameNode->mixHeatCapacity[kCurr];
	Double	sumCpDdYdx;
	Double	sumMH;

	if ( flame->GetSoot() ) {
		flame->GetSoot()->FillRHS( flame, nodeInfo, kSimilarity );
	}
	
// first fill all convection terms
	// first equation ( mass )
#ifdef FCONVLEFTTORIGHT
	rhs[fVVelocity] += FirstDerivUpwind( y[fVVelocity], yPrev[fVVelocity], hm );
#else
	rhs[fVVelocity] += FirstDerivUpwind( yNext[fVVelocity], y[fVVelocity], h );
#endif

#ifdef UPWINDCONVECTION
	// second equation ( momentum )
	rhs[fUVelocity] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fUVelocity], y[fUVelocity], yNext[fUVelocity], hm, h, FALSE );
	// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		rhs[speciesEq] += NonlinearConvectUpwind( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h, FALSE );
	}
	if ( fTemperature < M ) {
		rhs[fTemperature] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h, FALSE );
	}
#else
	// second equation ( momentum )
	rhs[fUVelocity] += NonlinearConvectCentral( y[fVVelocity], yPrev[fUVelocity], y[fUVelocity], yNext[fUVelocity], hm, h );
	// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		rhs[speciesEq] += NonlinearConvectCentral( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
	}
	if ( fTemperature < M ) {
		rhs[fTemperature] += NonlinearConvectCentral( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
	}
#endif

// mass equation
	rhs[fVVelocity] -= y[fUVelocity];

// momentum equation
	rhs[fUVelocity] += mixDensityInf / idealGasCoeff * y[fTemperature];
	rhs[fUVelocity] -= y[fUVelocity] * y[fUVelocity];
	rhs[fUVelocity] += constMassDiffCoeff * flame->StandardDiffusion( fUVelocity, mixViscosity, nodeInfo );
	
// fFirstSpecies to fFirstSpecies + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		eqLoop = speciesEq - fFirstSpecies;
#ifdef FULLDIFFUSION
		rhs[speciesEq] -= constMassDiffCoeff * flame->SecondDerivBinSpecDiff( speciesEq, nodeInfo );
#else
		rhs[speciesEq] += constMassDiffCoeff * flame->SecondDerivSpeciesDiffusion( speciesEq, nodeInfo );
#ifdef MOLARDIFF
		rhs[speciesEq] += constMassDiffCoeff * flame->SecondDerivXDiffusion( speciesEq, nodeInfo );
#endif
#endif

		if ( flame->UseDiffCorr() ) {
			rhs[speciesEq] -= constMassDiffCoeff * flame->NewDiffCorr( speciesEq, nodeInfo );
#ifdef MOLARDIFF
			rhs[speciesEq] -= flame->GetSolver()->bt->GetParameter() * constMassDiffCoeff * flame->NewDiffCorrX( speciesEq, nodeInfo );
#endif
		}
		if ( flame->fThermoDiffusion ) {
			rhs[speciesEq] += constMassDiffCoeff * flame->ThermoDiffusion( eqLoop, kSimilarity, nodeInfo );
		}
		rhs[speciesEq] += productionRate[eqLoop] * oneOverATimesRho;
	}
	
// energy equation
	if ( fTemperature < M ) {
		Double	diffCorrHeat = 0.0;
		Double	sumCpYD = 0.0;
		diffCorrHeat = ( flame->UseDiffCorr() ) ? flameNode->mixHeatCapacity[kCurr] : 0.0;
		Double	oneOverARhoCp = 1.0 / ( strainRate * mixDensity[kCurr] * mixHeatCapacity );
		sumCpDdYdx = 0.0;
		sumMH = 0.0;
		Double	sumYH = 0.0;

		rhs[fTemperature] += constThermDiffCoeff * flame->StandardDiffusion( fTemperature, flameNode->mixConductivity, nodeInfo );
//		rhs[fTemperature] += constThermDiffCoeff * (SecondDerivWeightedDiffusion( fTemperature, mixDensity[kCurr] * flameNode->mixConductivity[kCurr], nodeInfo ) 
//								+ FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) 
//								* FirstDeriv( mixDensity[kPrev] * flameNode->mixConductivity[kPrev], mixDensity[kCurr] * flameNode->mixConductivity[kCurr], mixDensity[kNext] * flameNode->mixConductivity[kNext], hm, h ) );

		for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
			sumCpDdYdx += ( heatCapacity[eqLoop] - diffCorrHeat ) * diffusivity[eqLoop] 
						* FirstDeriv( YPrev[eqLoop], Y[eqLoop], YNext[eqLoop], hm, h );
// check_ if next statement needs to be disabled for diffusivity correction switched off
			sumCpYD += ( heatCapacity[eqLoop] - diffCorrHeat ) * Y[eqLoop] * diffusivity[eqLoop];
			sumMH += productionRate[eqLoop] * enthalpy[eqLoop];
			sumYH += Y[eqLoop] * enthalpy[eqLoop];
		}
#	ifdef FULLDIFFUSION
		rhs[fTemperature] -=  1.0 * constMassDiffCoeff * flame->HeatFluxBinSpecDiff( nodeInfo )
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#	else
		rhs[fTemperature] +=  constThermDiffCoeff * sumCpDdYdx * mixDensity[kCurr] * mixDensity[kCurr]
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
		rhs[fTemperature] +=  constThermDiffCoeff * sumCpYD * mixDensity[kCurr] * mixDensity[kCurr]
						/ flameNode->mixMolarMass[kCurr]
						* FirstDeriv( flameNode->mixMolarMass[kPrev], flameNode->mixMolarMass[kCurr], flameNode->mixMolarMass[kNext], hm, h )
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#	endif
		rhs[fTemperature] -= sumMH * oneOverARhoCp;


//		rhs[fTemperature] = (flame->fEnthalpy - sumYH ) / hnenn;

		
		if ( flame->fProperties->GetRadiation() ) {
			rhs[fTemperature] += flameNode->radiation[kCurr] * oneOverARhoCp;
			if ( flame->GetSoot() ) {
				rhs[fTemperature] -= oneOverARhoCp * flame->GetSoot()->GetSootRadiation( y[fTemperature], flameNode->moments );
			}
		}
	}

	if ( flame->GetArcLengthCont() ) {
		if ( nodeInfo->gridPoint < flame->GetTMaxLoc() ) {
			rhs[fLnStrainrate] += FirstDerivUpwind( yNext[fLnStrainrate], y[fLnStrainrate], h );
		}
		else if ( nodeInfo->gridPoint == flame->GetTMaxLoc() ) {
			Double	dT = (y[fTemperature] - flame->GetTempContStart()) / flame->GetDeltaTref();
			Double	dlnStrainrate = (y[fLnStrainrate] - flame->GetStrainrateContStart()) / flame->GetDeltaStrainrateref();
			if ( flame->GetDeltaArcLength() == 0.0 ) {
				rhs[fLnStrainrate] += dlnStrainrate;
			}
			else {
				rhs[fLnStrainrate] += dT * flame->GetdTds() + dlnStrainrate * flame->GetdlnStrainrateds() 
								- flame->GetDeltaArcLength();
			}
		}
		else {
			rhs[fLnStrainrate] += FirstDerivUpwind( y[fLnStrainrate], yPrev[fLnStrainrate], hm );
		}
	}
	else {
		if ( flame->fStrainRateContin ) {
			if ( nodeInfo->gridPoint < flame->GetTMaxLoc() ) {
				rhs[fLnStrainrate] += FirstDerivUpwind( yNext[fLnStrainrate], y[fLnStrainrate], h );
			}
			else if ( nodeInfo->gridPoint == flame->GetTMaxLoc() ) {
				rhs[fLnStrainrate] += y[fTemperature] - flame->GetTempContStart();
			}
			else {
				rhs[fLnStrainrate] += FirstDerivUpwind( y[fLnStrainrate], yPrev[fLnStrainrate], hm );
			}
		}
		else {
			rhs[fLnStrainrate] += y[fLnStrainrate] - flame->GetStrainrateContStart();
		}
	}

	TTimePtr tim = flame->GetSolver()->time;
	for ( eqLoop = 0; eqLoop < M; ++eqLoop ) {
		if ( flame->GetSolver()->bt->GetTimedepFlag() && eqLoop != fVVelocity 
				&& !flame->GetSolver()->time->GetTimeConverged() ) {
			rhs[eqLoop] -= ( y[eqLoop] - tim->GetYOld()->mat[nodeInfo->gridPoint][eqLoop] ) 
								/ tim->GetDeltaT();
		}
		rhs[eqLoop] *= - hnenn;
	}
}

int CountPremSimPostIter( void *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	int			fUVelocity = flame->GetOffsetUVelocity();
	int			fVVelocity = flame->GetOffsetVVelocity();
	int			fTemperature = flame->GetOffsetTemperature();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();
	int			i;
	TNewtonPtr	bt = flame->GetSolver()->bt;
	NodeInfoPtr nodeInfo = bt->GetNodeInfo();
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = currGrid->GetNGridPoints();
	int			*bcFlagLeft = currGrid->GetBcFlagLeft();
	int			*bcFlagRight = currGrid->GetBcFlagRight();
	int			mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
	int			mixtureSpecificationRight = flame->GetMixtureSpecificationRight();
	MatrixPtr	yMat = currGrid->GetY();
	Double		**y = yMat->mat;
	VectorPtr	xVec = currGrid->GetX();
	Double		*x = xVec->vec;
	VectorPtr	yLeftVec = currGrid->GetYLeft();
	Double		*yLeft = yLeftVec->vec;
	VectorPtr	yRightVec = currGrid->GetYRight();
	Double		*yRight = yRightVec->vec;
	Double		*yLast = y[nGridPoints-1];
	Double		hLast = bt->GetRight() - x[nGridPoints-1];
	Double		*bcLeft = currGrid->GetBcLeft()->vec;
	Double		*bcRight = currGrid->GetBcRight()->vec;
	Double		hFirst = x[0] - bt->GetLeft();
	T1DSpeciesPtr	species = flame->GetSpecies();
	T1DPropertiesPtr	properties = flame->GetProperties();
	Double		*diffusivity = species->GetDiffusivity()->mat[0];
	Double		*density =  properties->GetDensity()->vec;
	Double		viscosityInf = properties->GetViscosity()->vec[nGridPoints];
	Double		pressure = flame->GetPressure();
	Double		&rhoLeft = density[-1];
	Double		&rhoRight = density[nGridPoints];
	Double		parameter = bt->GetParameter();
	Double		*temp = flame->GetTemperature()->vec;
	Double		**Y = flame->GetMassFracs()->mat;
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();
	int			variables = bt->GetNVariables();
	
// first set temperature and massfractions 
	for ( i = 0; i < nGridPoints; ++i ) {
		if (  fTemperature < bt->GetNEquations() ) {
			break;
		}
	}

	for ( i = 0; i < nGridPoints; ++i ) {
		if ( y[i][fTemperature] > 10000 ) {
			FILE *fp = flame->GetOutfile( "Error", TFlame::kData );
			flame->GetSolver()->bt->PrintSolution( x, y, flame->GetVariableNames(), fp );
			fclose( fp );
		}
		if ( flame->CheckSolution( y[i][fTemperature], &y[i][fFirstSpecies], nSpeciesInSystem ) ) {
			return 1;
		}
	}
	if ( flame->GetSoot() ) {
		flame->GetSoot()->PostIter( flame );
	}

//	fprintf(stderr, "strainrate = %g\n\n", flame->GetStrainRate() );

//	update properties
	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, temp[kPrev], Y[kPrev], pressure );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, temp[nGridPoints], Y[nGridPoints], pressure );

// left boundary
	CountPremSimUpdateLeftBoundary( object );

	CountPremSimUpdateRightBoundary( object );

	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	
	flame->UpdateThermoProps();


//	flame->SetTMaxLoc( flame->GetZRefLoc( currGrid->GetX() ) );
//	flame->SetTMaxLoc( LocationOfMax( nGridPoints, temp ) );
	if ( flame->GetArcLengthCont() || !flame->fStrainRateContin ) {
		flame->SetTMaxLoc( 3 );
	}

	return 0;
}

void TCountPremFlameSim::UpdateDimensions( int len )
{
	T1DFlame::UpdateDimensions( len );
	fSolV->len = len;
	fSolU->len = len;
	fSolLnStrainrate->len = len;
}

void TCountPremFlameSim::UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
{
	int		nGridPoints = yMat->cols;
	Double	*V = fSolV->vec;
	Double	*U = fSolU->vec;
	Double	*lnStrainrate = fSolLnStrainrate->vec;
	Double	**y = yMat->mat;
	Double	*yLeft = yLeftVec->vec;
	Double	*yRight = yRightVec->vec;

	UpdateDimensions( nGridPoints );

	T1DFlame::UpdateSolution( yMat, yLeftVec, yRightVec );
	
	V[kPrev] = yLeft[fVVelocity];
	U[kPrev] = yLeft[fUVelocity];
	lnStrainrate[kPrev] = yLeft[fLnStrainrate];
	for ( int k = 0; k < nGridPoints; ++k ) {
		V[k] = y[k][fVVelocity];
		U[k] = y[k][fUVelocity];
		lnStrainrate[k] = y[k][fLnStrainrate];
	}
	V[nGridPoints] = yRight[fVVelocity];
	U[nGridPoints] = yRight[fUVelocity];
	lnStrainrate[nGridPoints] = yRight[fLnStrainrate];
}

void TCountPremFlameSim::UpdateSolutionOnePoint( Double *y, int gridPoint )
{
	T1DFlame::UpdateSolution( y, gridPoint );
	
	fSolU->vec[gridPoint] = y[fUVelocity];
	fSolV->vec[gridPoint] = y[fVVelocity];
	fSolLnStrainrate->vec[gridPoint] = y[fLnStrainrate];
}

void TCountPremFlameSim::SolutionToSolver( void )
{
	TNewtonPtr	bt = fSolver->bt;
	TGridPtr	grid = bt->GetGrid()->GetFine();
	int		nGridPoints = fSolV->len;
	Double	*V = fSolV->vec;
	Double	*U = fSolU->vec;
	Double	*lnStrainrate = fSolLnStrainrate->vec;
	Double	**y = grid->GetY()->mat;

	T1DFlame::SolutionToSolver();
	
	for ( int k = 0; k < nGridPoints; ++k ) {
		y[k][fVVelocity] = V[k];
		y[k][fUVelocity] = U[k];
		y[k][fLnStrainrate] = lnStrainrate[k];
	}
	
	CountPremSimPostIter( this );
}

void TCountPremFlameSim::SaveSolution( void )
{
	int		k;
	int		len = fSolV->len;
	Double	*v = fSolV->vec;
	Double	*saveV = fSavedV->vec;
	Double	*u = fSolU->vec;
	Double	*saveU = fSavedU->vec;
	Double	*lnStrainrate = fSolLnStrainrate->vec;
	Double	*saveLnStrainrate = fSavedLnStrainrate->vec;

	T1DFlame::SaveSolution();
	fSavedV->len = fSolV->len;
	fSavedU->len = fSolU->len;
	fSavedLnStrainrate->len = fSolLnStrainrate->len;

	for ( k = -1; k <= len; ++k ) {
		saveV[k] = v[k];
		saveU[k] = u[k];
		saveLnStrainrate[k] = lnStrainrate[k];
	}

	fSavedTempContStart = fSolTemp->vec[GetTMaxLoc()];
	fSavedLnStrainrateContStart = fSolLnStrainrate->vec[1];
}

void TCountPremFlameSim::RestoreSolution( void )
{
	int		k;
	int		len = fSavedV->len;
	Double	*v = fSolV->vec;
	Double	*saveV = fSavedV->vec;
	Double	*u = fSolU->vec;
	Double	*saveU = fSavedU->vec;
	Double	*lnStrainrate = fSolLnStrainrate->vec;
	Double	*saveLnStrainrate = fSavedLnStrainrate->vec;

	UpdateDimensions( len );

	T1DFlame::RestoreSolution();

	for ( k = -1; k <= len; ++k ) {
		v[k] = saveV[k];
		u[k] = saveU[k];
		lnStrainrate[k] = saveLnStrainrate[k];
	}
	
	fTempContStart = fSavedTempContStart;
	fLnStrainrateContStart = fSavedLnStrainrateContStart;

	SolutionToSolver();
}

int	TCountPremFlameSim::GetOffsetVVelocity( void )
{
	return fVVelocity; 
}

int	TCountPremFlameSim::GetOffsetUVelocity( void )
{
	return fUVelocity; 
}

int	TCountPremFlameSim::GetOffsetLnStrainrate( void )
{
	return fLnStrainrate; 
}

int	TCountPremFlameSim::GetOffsetTemperature( void )
{
	return fTemperature; 
}

int	TCountPremFlameSim::GetOffsetFirstSpecies( void ) 
{
	return fFirstSpecies;
}

int	TCountPremFlameSim::GetOffsetMixFrac( void )
{
	cerr << "#error: class has no member fMixtureFraction" << NEWL;
	exit( 2 );
	return 0; 
}

ConstStringArray TCountPremFlameSim::GetVariableNames( void )
{
	return ( ConstStringArray ) fVariableNames;
}

int TCountPremFlameSim::GetVariablesWithoutSpecies( void )
{
	return fVariablesWithoutSpecies;
}

void TCountPremFlameSim::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp )
{
	int 				i, j, k;
	TBVPSolverPtr		solver = GetSolver();
	TNewtonPtr			bt = solver->bt;
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
	TAdaptiveGridPtr	adapGrid = bt->GetGrid();
	TGridPtr			grid = adapGrid->GetFine();
	int					variables = bt->GetNVariables();
	int					nGridPoints = grid->GetNGridPoints();
	int					maxGridPoints = bt->GetMaxGridPoints();
	int					initialGridPoints = bt->GetInitialGridpoints();
	int					nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	int					speciesIndex;
	MatrixPtr			yMat = grid->GetY();
	VectorPtr			yLeftVec = grid->GetYLeft();
	VectorPtr			yRightVec = grid->GetYRight();
	Double			 	*yLeft = yLeftVec->vec;
	Double 				*yRight = yRightVec->vec;
	Double				left = inp->fLeft;
	Double				right = inp->fRight;
	Double				*locX = grid->GetX()->vec;
	int					gridPointsIn = sp->gridPoints;	// including boundaries
	Double				*yWork = adapGrid->GetWorkVector()->vec;
	Flag				etaSet = FALSE;
	Flag				chooseInputGrid = FALSE;
	Double				*xIn = new Double[gridPointsIn];
	if ( !xIn ) FatalError( "memory allocation of TCountPremFlameSim failed" );
	Double				*yIn =  new Double[gridPointsIn];
	if ( !yIn ) FatalError( "memory allocation of TCountPremFlameSim failed" );
	Double				*yInFloat = sp->data;
	Double				**y = grid->GetY()->mat;
	int					variable;
	char				*string = sp->labels;
	SplinePtr			theSpline = NULL;
	Double				leftSlope;
	Double				rightSlope;
	FILE				*fp;
	Double				*temp = GetTemperature()->vec;
	Double				**Y = GetMassFracs()->mat;
	Double				strainRateIn;
	struct _parameter	*param = GetParameter( "strainrate" );
	

	if ( !fStrainRate || fStrainRate->vec[0] <= 0.0 ) {
	// get strainrate from input file
		if ( param ) {
			strainRateIn = (Double)param->what.quantity.value;
		}
		else { // choose default
			cerr << "#warning: no value for 'strainrate' in inputfile" << NEWL;
			strainRateIn = 10.0;
		}
		
		if ( !fStrainRate ) {
			fStrainRate = NewVector( 1 );
		}
		fStrainRate->vec[0] = strainRateIn;
		fStrainRate->len = 0;
		cerr << "initial strainrate is " << GetStrainRate() << NEWL;
	}
	
//	choose grid from input or equidistant
	if ( gridPointsIn <= maxGridPoints+2 && gridPointsIn > initialGridPoints+2 && ( gridPointsIn % 2 ) != 0 ) {
		grid->AdjustNGridPoints( gridPointsIn-2 );
		solver->UpdateAllDimensions( gridPointsIn-2 );	
		nGridPoints = grid->GetNGridPoints();
		chooseInputGrid = TRUE;
	}
	else {
		nGridPoints = initialGridPoints;
		grid->AdjustNGridPoints( nGridPoints );
		solver->UpdateAllDimensions( nGridPoints );	
		chooseInputGrid = FALSE;
	}
	
// find independent coordinate
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strcmp( string, "eta" ) == 0 ) {
			if ( chooseInputGrid ) {
				cerr << "choose inputGrid" << NEWL;
				for ( j = 0; j < gridPointsIn-2; ++j ) {
					locX[j] = yInFloat[i*gridPointsIn + j+1];		// implicit cast from float to Double
				}
				bt->SetLeft( yInFloat[i*gridPointsIn] );
				bt->SetRight( yInFloat[(i+1)*gridPointsIn-1] );
				left = bt->GetLeft();
				right = bt->GetRight();
				etaSet = TRUE;
			}
			else { // choose own Grid, but read x for interpolation
				cerr << "choose own Grid" << NEWL;
				for ( j = 0; j < gridPointsIn; ++j ) {
					xIn[j] = yInFloat[i*gridPointsIn + j];		// implicit cast from float to Double
				}
				bt->SetLeft( left );
				bt->SetRight( right );
				if ( xIn[0] > left ) {
					xIn[0] = left;
				}
				if ( xIn[gridPointsIn-1] < right ) {
					xIn[gridPointsIn-1] = right;
				}
				grid->Make_equi_Grid();
				etaSet = TRUE;
			}
		}
		string += strlen(string) + 1;
	}

// set default values
	for ( i = 0; i < nGridPoints; ++i ) {
		for ( j = 0; j < variables; ++j ) {
			y[i][j] = yLeft[j] + ( yRight[j] - yLeft[j] ) / ( right - left ) * locX[i];
		}
	}

// error checking
	if ( !etaSet ) {
		cerr << "error: can't find coordinate 'eta'" << NEWL;
		exit(2);
	}
// reset string
	string = sp->labels;
	
	for ( i = 0; i < sp->variables; ++i ) {
		if ( ( strcmp( string, "df/deta" ) == 0 ) || ( strcmp( string, "u" ) == 0 )
					|| ( strcmp( string, "f'" ) == 0 ) ) {
			variable = fUVelocity;
		}
		else if ( strcmp( string, "f" ) == 0 ) {
			variable = fVVelocity;
		}
		else if ( strncmp( string, "temperature", 11 ) == 0 ) {
			variable = fTemperature;
		}
		else if ( GetSoot() && strncmp( string, "conc-soot", 9 ) == 0 ) {
			string += 9;
//			cerr << "string = " << string << NEWL;
			int	num = atoi( string );
			if ( num < GetSoot()->GetNSootMoments() ) {
				variable = num + GetSoot()->GetOffsetSootMoments();
			}
			else {
				string += strlen(string) + 1;
				continue;
			}
		}
		else if ( strncmp( string, "massfraction-", 13 ) == 0 ){
			string += 13;
			UpperString( string );
			if ( ( speciesIndex = inp->FindSpecies( string ) ) >= 0 ) {
				if ( speciesIndex < nSpeciesInSystem ) {
					variable = fFirstSpecies + speciesIndex;
				}
				else {
					string += strlen(string) + 1;
					continue;
				}
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
			for ( k = 0; k < gridPointsIn-2; ++k ) {	// transform to current coordinate system
				if ( variable != fVVelocity || yInFloat[i*gridPointsIn] < yInFloat[i*gridPointsIn + 3] ) {
					y[k][variable] = yInFloat[i*gridPointsIn + k+1];	// copy workspace to vector of solution
				}
				else {
					y[k][variable] = -yInFloat[i*gridPointsIn + k+1];	// copy workspace to vector of solution and transform f
				}
			}
		}
		else {
			for ( k = 0; k < gridPointsIn; ++k ) {	// store vector in workspace
				yIn[k] = yInFloat[i * gridPointsIn + k];	// implicit cast from float to Double
			}
		
			if ( variable == fVVelocity ) {
				yIn[0] = yIn[1] - ( yIn[2] - yIn[1] ) / ( xIn[2] - xIn[1] ) * ( xIn[1] - xIn[0] );
				yIn[gridPointsIn-1] = yIn[gridPointsIn-2] + ( yIn[gridPointsIn-2] - yIn[gridPointsIn-3] ) 
										/ ( xIn[gridPointsIn-2] - xIn[gridPointsIn-3] ) 
										* ( xIn[gridPointsIn-1] - xIn[gridPointsIn-2] );
			}
					
			leftSlope = ( yIn[1] - yIn[0] ) / ( xIn[1] - xIn[0] );
			rightSlope = ( yIn[gridPointsIn-1] - yIn[gridPointsIn-2] ) / ( xIn[gridPointsIn-1] - xIn[gridPointsIn-2] );
			theSpline = ComputeSimpleSpline( xIn, yIn, gridPointsIn, FALSE, leftSlope, FALSE, rightSlope, NULL, TRUE );
			SplineInterpolate( theSpline, locX, yWork, nGridPoints );
			for ( k = 0; k < nGridPoints; ++k ) {	// transform to current coordinate system
				if ( variable != fVVelocity || yInFloat[i*gridPointsIn] < yInFloat[i*gridPointsIn + 3] ) {
					y[k][variable] = yWork[k];	// copy workspace to vector of solution
				}
				else {
					y[k][variable] = -yWork[k];	// copy workspace to vector of solution and transform f
				}
			}
		}
		//	set bc
		if ( variable == fTemperature ) {
			if ( yLeft[fTemperature] <= 0.0 ) {
				yLeft[fTemperature] = yInFloat[i*gridPointsIn];
			}
			if ( yRight[fTemperature] <= 0.0 ) {
				yRight[fTemperature] = yInFloat[(i+1)*gridPointsIn-1];
			}
		}
	}

// set left to zero
	for ( k = 0; k < nGridPoints; ++k ) {
		locX[k] -= bt->GetLeft();
	}
	bt->SetRight( bt->GetRight() - bt->GetLeft() );
	bt->SetLeft( 0.0 );
	
// set pressure
	if ( GetPressure() <= 0.0 ) {
		param = GetParameter( "pressure" );
		Double	thePressure;
		if ( param ) {
			thePressure = (Double)param->what.quantity.value;
			if ( strcmp( param->what.quantity.unit, "bar" ) == 0 ) {
				thePressure *= 1.0e5;
			}
			SetPressure( thePressure );
		}
		else { // exit
			cerr << "#error: no value for 'pressure' in inputfile" << NEWL;
			exit(2);
		}
	}

// set initial Boundary values
//	update properties
	UpdateSolution( yMat, yLeftVec, yRightVec );
	if ( GetSpecies()->IsConstantLewisNumber() ) {
		CompLewisNumbers( GetSpecies()->GetLewisNumberFile() );
	}
	SetFlameNode( kPrev );
	ComputeProperties( fFlameNode, temp[kPrev], Y[kPrev], GetPressure() );
	SetFlameNode( nGridPoints );
	ComputeProperties( fFlameNode, temp[nGridPoints], Y[nGridPoints], GetPressure() );
	
	if ( fPremConfiguration == kBackToBack ) { // symmetric around stagnation plane, solve only right part
	}
	else if ( fPremConfiguration == kFreshToBurned ) { // zero gradient for species and temp
	}
	else if ( fPremConfiguration == kFull ) {
		Double	*rho = GetProperties()->GetDensity()->vec;
#ifdef FCONVLEFTTORIGHT
		yLeft[fUVelocity] = sqrt( rho[nGridPoints] / rho[-1] );
		yLeft[fVVelocity] = y[0][fVVelocity] - ( locX[0] - bt->GetLeft() ) * yLeft[fUVelocity];
#else
		yRight[fUVelocity] = 1.0;
		yRight[fVVelocity] = y[nGridPoints-1][fVVelocity] + (bt->GetRight() - locX[nGridPoints-1]) * yRight[fUVelocity];
#endif
	}
	else {
		fprintf( stderr, "#error in premixed configuration\n" );
		exit(2);
	}

	for ( i = 0; i < nGridPoints; ++i ) {
		y[i][fLnStrainrate] = log( fStrainRate->vec[0] );
	}
	yLeft[fLnStrainrate] = yRight[fLnStrainrate] = log( fStrainRate->vec[0] );

	CountPremSimPostIter( this );

	FreeSpline( theSpline );
	delete yIn;
	delete xIn;

	adapGrid->SetSolutionScaler();
	
	fp = GetOutfile( "initialguess", TFlame::kData );
	bt->PrintSolution( locX, y, GetVariableNames(), fp );
	fclose(fp);
}


Double TCountPremFlameSim::HeatFluxBinSpecDiff( NodeInfoPtr nodeInfo )
{
	int		j, i;
	int		nVarj;
	int		nSpecIn = fSpecies->GetNSpeciesInSystem();
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	***Dij = fFlameNode->binDij;
	Double	*rho = fFlameNode->mixDensity;
	Double	*Wi = fSpecies->GetMolarMass()->vec;
	Double	W2C = fFlameNode->mixMolarMass[kCurr] * fFlameNode->mixMolarMass[kCurr];
	Double	WP = fFlameNode->mixMolarMass[kPrev];
	Double	WC = fFlameNode->mixMolarMass[kCurr];
	Double	WN = fFlameNode->mixMolarMass[kNext];
	Double	cpMix = fFlameNode->mixHeatCapacity[kCurr];
	Double	*cp = fFlameNode->heatCapacity;
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	sumk = 0.0, sum = 0.0;
	Double	rho2overcpM2 = rho[kCurr] * rho[kCurr] / W2C / cpMix;
	
	for ( i = 0; i < nSpecIn; ++i ) {
		sumk = 0.0;
		for ( j = 0; j < nSpecIn; ++j ) {
			if ( j == i ) continue;
			nVarj = j + fFirstSpecies;
			sumk += Dij[kCurr][i][j] * FirstDeriv( yPrev[nVarj] * WP
												, y[nVarj] * WC
												, yNext[nVarj] * WN, hm, h );
		}
		sum += rho2overcpM2 * Wi[i] * cp[i] * sumk;
	}
	
	return sum;
}

Double TCountPremFlameSim::SecondDerivBinSpecDiff( int nVariable, NodeInfoPtr nodeInfo )
{
	int		j, i = nVariable - fFirstSpecies;
	int		nVarj;
	int		nSpecIn = fSpecies->GetNSpeciesInSystem();
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	***Dij = fFlameNode->binDij;
	Double	*rho = fFlameNode->mixDensity;
	Double	*Wi = fSpecies->GetMolarMass()->vec;
	Double	W2P = fFlameNode->mixMolarMass[kPrev] * fFlameNode->mixMolarMass[kPrev];
	Double	W2C = fFlameNode->mixMolarMass[kCurr] * fFlameNode->mixMolarMass[kCurr];
	Double	W2N = fFlameNode->mixMolarMass[kNext] * fFlameNode->mixMolarMass[kNext];
	Double	WP = fFlameNode->mixMolarMass[kPrev];
	Double	WC = fFlameNode->mixMolarMass[kCurr];
	Double	WN = fFlameNode->mixMolarMass[kNext];
	Double	aNext, aCurr, aPrev;
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	sum = 0.0;
	
	for ( j = 0; j < nSpecIn; ++j ) {
		if ( j == i ) continue;
		nVarj = j + fFirstSpecies;
		aCurr = rho[kCurr] * rho[kCurr] * Wi[i] * Wi[j] / W2C * Dij[kCurr][i][j];
		aPrev = rho[kPrev] * rho[kPrev] * Wi[i] * Wi[j] / W2P * Dij[kPrev][i][j];
		aNext = rho[kNext] * rho[kNext] * Wi[i] * Wi[j] / W2N * Dij[kNext][i][j];
		sum += ( ( aCurr + aNext ) * hm * ( yNext[nVarj] * WN - y[nVarj] * WC )
				+ ( aPrev + aCurr ) * h * ( yPrev[nVarj] * WP - y[nVarj] * WC ) ) / Wi[j];
	}
	
	return sum / nodeInfo->hnenn;

}

Double TCountPremFlameSim::SecondDerivSpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo )
{
	int		speciesIndex = nVariable - fFirstSpecies;
	Double	yPrev = nodeInfo->yPrev[nVariable];
	Double	y = nodeInfo->y[nVariable];
	Double	yNext = nodeInfo->yNext[nVariable];
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	diffPlus = diffusivity[speciesIndex] * mixDensity[kCurr] * mixDensity[kCurr]
					+ diffusivityNext[speciesIndex] * mixDensity[kNext] * mixDensity[kNext];
	Double	diffMinus = diffusivityPrev[speciesIndex] * mixDensity[kPrev] * mixDensity[kPrev]
					+ diffusivity[speciesIndex] * mixDensity[kCurr] * mixDensity[kCurr];
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	

	return ( diffPlus * hm * ( yNext - y ) + diffMinus * h * ( yPrev - y ) ) 
				/ nodeInfo->hnenn;
}

/* Double TCountPremFlameSim::SecondDerivSpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo ) */
/* { */
/* 	int		speciesIndex = nVariable - fFirstSpecies; */
/* 	Double	yPrev = nodeInfo->yPrev[nVariable]; */
/* 	Double	y = nodeInfo->y[nVariable]; */
/* 	Double	yNext = nodeInfo->yNext[nVariable]; */
/* 	Double	*diffusivity = fFlameNode->diffusivity; */
/* 	Double	*diffusivityNext = fFlameNode->diffusivityNext; */
/* 	Double	*diffusivityPrev = fFlameNode->diffusivityPrev; */
/* 	Double	*mixDensity = fFlameNode->mixDensity; */
/* 	Double	diffPlus = diffusivityNext[speciesIndex] * mixDensity[kNext] * mixDensity[kNext]; */
/* 	Double	diffCurr = diffusivity[speciesIndex] * mixDensity[kCurr] * mixDensity[kCurr]; */
/* 	Double	diffMinus = diffusivityPrev[speciesIndex] * mixDensity[kPrev] * mixDensity[kPrev]; */
/* 	Double	hm = nodeInfo->hm; */
/* 	Double	h = nodeInfo->h; */
/* 	 */
/*  */
/* 	return SecondDerivWeightedDiffusion( nVariable, diffCurr, nodeInfo ) + FirstDeriv( yPrev, y, yNext, hm, h ) * FirstDeriv( diffMinus, diffCurr, diffPlus, hm, h ); */
/* } */

Double TCountPremFlameSim::SecondDerivXDiffusion( int nVariable, NodeInfoPtr nodeInfo )
{
	int		speciesIndex = nVariable - fFirstSpecies;
	Double	MPrev = fFlameNode->mixMolarMass[kPrev];
	Double	M = fFlameNode->mixMolarMass[kCurr];
	Double	MNext = fFlameNode->mixMolarMass[kNext];
	Double	yPrev = nodeInfo->yPrev[nVariable];
	Double	y = nodeInfo->y[nVariable];
	Double	yNext = nodeInfo->yNext[nVariable];
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	
	Double	diffPlusX = diffusivity[speciesIndex] * mixDensity[kCurr] * mixDensity[kCurr] * y / M 
					+ diffusivityNext[speciesIndex] * mixDensity[kNext] * mixDensity[kNext] * yNext / MNext;
	Double	diffMinusX = diffusivityPrev[speciesIndex] * mixDensity[kPrev] * mixDensity[kPrev] * yPrev / MPrev
					+ diffusivity[speciesIndex] * mixDensity[kCurr] * mixDensity[kCurr] * y / M;

	return ( diffPlusX * hm * ( MNext - M ) + diffMinusX * h * ( MPrev - M ) ) 
				/ nodeInfo->hnenn;
}

void TCountPremFlameSim::PrintRHSTemp( TNewtonPtr bt )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		k;
    int         		N = currentGrid->GetNGridPoints();
	VectorPtr 			physXVec = NewVector( N + 2 );
	Double				*physX = physXVec->vec;
	FILE				*fp = GetOutfile( "tempRHS", TFlame::kData );
	
	UpdateThermoProps();
		
	EtaToX( bt, physXVec );
		
	fprintf( fp, "*\n%-12s\t%-12s\t%-12s\t%-12s", "eta", "y", "convection", "diffusion" );
	fprintf( fp, "\t%-12s", "entflux" );
	if ( UseDiffCorr() ) {
		fprintf( fp, "\t%-12s", "diffCorr" );
	}
	fprintf( fp, "\t%-12s", "production" );
	if ( fProperties->GetRadiation() ) {
		fprintf( fp, "\t%-12s", "radiation" );
		if ( GetSoot() ) {
			fprintf( fp, "\t%-12s", "sootRad" );
		}
	}
	fprintf( fp, "\n" );
	
	for ( k = 0; k < N; ++k ){
		bt->SetNodeInfo( this, k );
		PrintRHSTemp( bt, nodeInfo, physX[k+1], fp );
	}
    fclose( fp );
}

void TCountPremFlameSim::PrintRHSTemp( TNewtonPtr bt, NodeInfoPtr nodeInfo, Double physX, FILE *fp )
{
	int		eqLoop, speciesEq;
    int		M = bt->GetNEquations();
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	strainRate = GetStrainRate();
	Double	*enthalpy = fFlameNode->enthalpy;
	Double	mixDensity = *fFlameNode->mixDensity;
	Double	mixDensityInf = fFlameNode->rhoInf;
	Double	mixViscosityInf = fFlameNode->viscosityInf;
	Double	*productionRate = fFlameNode->productionRate;
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*heatCapacity = fFlameNode->heatCapacity;
	Double	oneOverATimesRho = 1.0 / ( strainRate * mixDensity );
	Double	idealGasCoeff = GetPressure() * *fFlameNode->mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	mixHeatCapacity = *fFlameNode->mixHeatCapacity;
	Double	constThermDiffCoeff = constMassDiffCoeff / fFlameNode->mixHeatCapacity[kCurr];
	Double	sumCpDdYdx;
	Double	sumMH;
	Double	*diffCorr = fFlameNode->diffCorr;

	fprintf( fp, "%-.6e\t%-.6e", *nodeInfo->x, physX );
#ifdef UPWINDCONVECTION
	fprintf( fp, "\t%-.6e", NonlinearConvectUpwind( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) );
#else
	fprintf( fp, "\t%-.6e", NonlinearConvectCentral( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) );
#endif

// energy equation
	Double	sumCpY = 0.0;
	Double	corrCoeff = 0.0;
	Double	oneOverARhoCp = 1.0 / ( strainRate * mixDensity * mixHeatCapacity );
	sumCpDdYdx = 0.0;
	sumMH = 0.0;

	fprintf( fp, "\t%-.6e", constThermDiffCoeff * StandardDiffusion( fTemperature, fFlameNode->mixConductivity, nodeInfo ) );
	if ( UseDiffCorr() ) {
		corrCoeff = constMassDiffCoeff * mixDensity * mixDensity
								* diffCorr[kCurr] / mixHeatCapacity;
	}
	for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
		speciesEq = fFirstSpecies+eqLoop;
		sumCpDdYdx += heatCapacity[eqLoop] * diffusivity[eqLoop] 
					* FirstDeriv( yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
		if ( UseDiffCorr() ) {
			sumCpY += heatCapacity[eqLoop] * y[speciesEq];
		}
		sumMH += productionRate[eqLoop] * enthalpy[eqLoop];
	}
	fprintf( fp, "\t%-.6e", constThermDiffCoeff * sumCpDdYdx * mixDensity * mixDensity
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) );
	if ( UseDiffCorr() ) {
		fprintf( fp, "\t%-.6e", -corrCoeff * sumCpY
					* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) );
	}
	fprintf( fp, "\t%-.6e", -sumMH * oneOverARhoCp );
	if ( fProperties->GetRadiation() ) {
		fprintf( fp, "\t%-.6e", fFlameNode->radiation[kCurr] * oneOverATimesRho / mixHeatCapacity );
		if ( GetSoot() ) {
			fprintf( fp, "\t%-.6e", -oneOverATimesRho / mixHeatCapacity * GetSoot()->GetSootRadiation( y[fTemperature], fFlameNode->moments ) );
		}
	}
	fprintf( fp, "\n" );
}

void TCountPremFlameSim::PrintRHSSpecies( TNewtonPtr bt )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		i, k;
    int         		N = currentGrid->GetNGridPoints();
	VectorPtr 			physXVec = NewVector( N + 2 );
	Double				*physX = physXVec->vec;
	char				**names = GetSpecies()->GetNames();
	int					nSpecIn = fSpecies->GetNSpeciesInSystem();
	int					num = 7;
	int					start = 0, end;
	int					fileCount = 1;
	FILE				*fp;
	
	UpdateThermoProps();
	
	EtaToX( bt, physXVec );
		
	do {
		end = ( int )MIN( start + 255 / num, nSpecIn );
		sprintf( GetOutFileBuff(), "%sspeciesRHS%d.dout", GetOutputPath(), fileCount++ );
		if ( !( fp = fopen( GetOutFileBuff(), "w") ) ) { 
			cerr << "#warning: unable to open file " << GetOutFileBuff() << NEWL;
			exit(2);
		}
	
		fprintf( fp, "*\n%-12s\t%-12s", "eta", "y" );
		for ( i = start; i < end; ++i ) {
			fprintf( fp, "\tConv_%-7s\tDiff_%-7s", names[i], names[i] );
			if ( UseDiffCorr() ) {
				fprintf( fp, "\tDiCo_%-7s", names[i] );
			}
			if ( fThermoDiffusion ) {
				fprintf( fp, "\tThDi_%-7s", names[i] );
			}
			fprintf( fp, "\tProd_%-7s\tCons_%-7s\tSource_%-7s", names[i], names[i], names[i] );
		}
		fprintf( fp, "\n" );
			
		for ( k = 0; k < N; ++k ){
			bt->SetNodeInfo( this, k );
			PrintRHSSpecies( start, end, nodeInfo, physX[k+1], fp );
		}
		fclose( fp );
		start = end;
	} while ( end < nSpecIn );
}

void TCountPremFlameSim::PrintRHSSpecies( int start, int end, NodeInfoPtr nodeInfo, Double physX, FILE *fp )
{
	int		i, j;
	int		speciesEq;
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	strainRate = GetStrainRate();
	Double	mixDensityInf = fFlameNode->rhoInf;
	Double	mixViscosityInf = fFlameNode->viscosityInf;
	Double	idealGasCoeff = GetPressure() * *fFlameNode->mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = 1.0 / ( strainRate * idealGasCoeff );
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	diffTerm;
	Double	*reactionRate = fFlameNode->reactionRate;
	Double	productionRate;
	Double	source;
	Double	sink;
	int		*nOfUsedReactions = fSpecies->GetNOfUsedReactions()->vec;
	VectorPtr	*nu = fSpecies->GetNu();
	IntVectorPtr	*usedReactions = fSpecies->GetUsedReactions();
	Double	*molarMass =  fSpecies->GetMolarMass()->vec;

	fprintf( fp, "%-.6e\t%-.6e", *nodeInfo->x, physX );

	for ( i = start; i < end; ++i ) {
		speciesEq = fFirstSpecies+i;

#ifdef UPWINDCONVECTION
		fprintf( fp, "\t%-.6e", NonlinearConvectUpwind( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h ) );
#else
		fprintf( fp, "\t%-.6e", NonlinearConvectCentral( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h ) );
#endif

		diffTerm = constMassDiffCoeff * SecondDerivSpeciesDiffusion( speciesEq, nodeInfo );
		fprintf( fp, "\t%-.6e", diffTerm );
		if ( UseDiffCorr() ) {
			diffTerm = -constMassDiffCoeff * NewDiffCorr( speciesEq, nodeInfo );
			fprintf( fp, "\t%-.6e", diffTerm );
		}
		if ( fThermoDiffusion ) {
			diffTerm = constMassDiffCoeff * ThermoDiffusion( i, kSimilarity, nodeInfo );
			fprintf( fp, "\t%-.6e", diffTerm );
		}

		sink = source = 0.0;
		for ( j = 0; j < nOfUsedReactions[i]; ++j ) {
			productionRate = nu[i]->vec[j] * reactionRate[usedReactions[i]->vec[j]];
			if ( productionRate > 0.0 ) {
				sink += productionRate;
			}
			else {
				source += productionRate;
			}
		}
		sink *= - molarMass[i] * y[fTemperature] * ROverAPM;
		source *= - molarMass[i] * y[fTemperature] * ROverAPM;
		fprintf( fp, "\t%-.6e\t%-.6e\t%-.6e", source, -sink, source+sink );
	}
	fprintf( fp, "\n" );
}

void CountPremSimOutput( void *object, FILE *fp, char* tail )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	TNewtonPtr		bt = flame->GetSolver()->bt;
	T1DPropertiesPtr	props = flame->GetProperties();
	TSpeciesPtr		species = flame->GetSpecies();
	NodeInfoPtr		nodeInfo = bt->GetNodeInfo();
	Double			*rho = props->GetDensity()->vec;
	TGridPtr		currentGrid = bt->GetGrid()->GetCurrentGrid();
	Double			*x = currentGrid->GetX()->vec;
	Double			*mixMolarMass = props->GetMolarMass()->vec;
	Double			*molarMass = species->GetMolarMass()->vec;
	Double			**massFracs = flame->GetMassFracs()->mat;
	Double			*temp = flame->GetTemperature()->vec;
	Double			*V = flame->GetV()->vec;
	Double			*U = flame->GetU()->vec;
	Double			**y = currentGrid->GetY()->mat;
	Double			*yLeft = currentGrid->GetYLeft()->vec,
					*yRight = currentGrid->GetYRight()->vec;
	int				i, k;
	int				gridPoints = currentGrid->GetNGridPoints();
	int				nOfSpecies = flame->GetSpecies()->GetNOfSpecies();
	int				nOfVariables = bt->GetNVariables();
	int				nOfEquations = bt->GetNEquations();
	int				firstSpecies = flame->GetOffsetFirstSpecies();
	int				tempOffset = flame->GetOffsetTemperature();
	int				fUVelocity = flame->GetOffsetUVelocity();
	int				fVVelocity = flame->GetOffsetVVelocity();
	time_t			theDate;
	char			buffer[80];
	ConstStringArray	varNames = flame->GetVariableNames();
	char			**names = species->GetNames();
	VectorPtr 		physXVec = NewVector( gridPoints + 2 );
	Double			*physX = physXVec->vec;
	VectorPtr 		LewisFuelVec = NewVector( gridPoints + 2 );
	Double			*LewisFuel = &LewisFuelVec->vec[kNext];
	Flag			fpOpen = FALSE;
	
	if ( !fp ) {
		fpOpen = TRUE;
		fp = flame->GetOutputFile( NULL, tail, TFlame::kNone );
	}

	flame->EtaToX( bt, physXVec );

// write header
	fprintf( fp, "header\n\n" );

	fprintf( fp, "title = \"planar counterflow diffusion flame\"\n" );
	fprintf( fp, "author = \"%s\"\n", flame->GetAuthor() );
	time( &theDate );
	strcpy( buffer, ctime( &theDate ) );
	if ( buffer[strlen(buffer)-1] == '\n' )
		buffer[strlen(buffer)-1] = '\0';
	fprintf( fp, "date = \"%s\"\n\n", buffer );
	for ( i = 0; i < flame->GetNFuels(); ++i ) {
		fprintf( fp, "fuel = \"%s\"\n", varNames[firstSpecies+flame->GetFuelIndex( i )] );
	}
	fprintf( fp, "pressure = %g [bar]\n", flame->GetPressure() / 1.0e5 );
	fprintf( fp, "strainrate = %g [1/s]\n", flame->GetStrainRate() );

	if ( species->IsConstantLewisNumber() ) {
		fprintf( fp, "ConstantLewisNumbers = \"True\"\n" );
	}
	
	fprintf( fp, "burningVelocity = %g [cm/s]\n", flame->GetBurningVelocity() * 100.0 );
	fprintf( fp, "Tmax = %g [K]\n", temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
	
	Double	EIFuel = 0.0;
	for ( i = 0; i < flame->GetNFuels(); ++i ) {
		EIFuel += flame->ComputeEmissionIndex( flame->GetFuelIndex( i ), &physX[1] );
	}

	int	indNO = species->FindSpecies( "NO" );
	if ( indNO > -1 ) {
		fprintf( fp, "EmissionIndexNO = %g [g/kg]\n"
				, -1000.0 * flame->ComputeEmissionIndex( indNO, &physX[1] ) 
				/ EIFuel );
	}
	
	int	indNO2 = species->FindSpecies( "NO2" );
	if ( indNO2 > -1 ) {
		fprintf( fp, "EmissionIndexNO2 = %g [g/kg]\n"
				, -1000.0 * flame->ComputeEmissionIndex( indNO2, &physX[1] ) 
				/ EIFuel );
	}
	fprintf( fp, "dTds = %g\n", flame->GetdTds() );
	fprintf( fp, "dlnads = %g\n", flame->GetdlnStrainrateds() );
	fprintf( fp, "TContStart = %g\n", flame->GetTempContStart() );
	fprintf( fp, "dlnaStart = %g\n", flame->GetStrainrateContStart() );
	fprintf( fp, "ds = %g\n", flame->GetDeltaArcLength() );

	
//  write bc
	Double	locMoleMass;
	fprintf( fp, "Unburned\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", temp[kPrev] );
	for ( i = 0; i < nOfSpecies; ++i ) {
		if ( fabs( massFracs[kPrev][i] ) > 1.0e-5 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], massFracs[kPrev][i] );
		}
	}
	for ( i = 0; i < nOfSpecies; ++i ) { // write X_i
		locMoleMass = massFracs[kPrev][i] * mixMolarMass[kPrev] / molarMass[i];
		if ( fabs( locMoleMass ) > 1.0e-5 ) {
			fprintf( fp, "\tMolefraction-%s = %g\n", names[i], locMoleMass );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "numOfSpecies = %d\n", nOfSpecies );
	fprintf( fp, "gridPoints = %d\n\n", gridPoints+2 );

	fprintf( fp, "body\n" );

// write independent coordinate
	fprintf( fp, "eta\n" );
	fprintf( fp, "\t%-.6e", bt->GetLeft() );
	for ( k = 0; k < gridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", x[k] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", bt->GetRight() );
	
// write physical coordinate
	
	fprintf( fp, "y [m]\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", physX[k] );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}
		
// write solution
	// write f-Velocity, f'-Velocity, mixture fraction and temperature
	flame->PrintFlameletVector( gridPoints+2, &V[kPrev], "f", fp );
	flame->PrintFlameletVector( gridPoints+2, &U[kPrev], "df/deta", fp );
	flame->PrintFlameletVector( gridPoints+2, &temp[kPrev], "temperature [K]", fp );

	// write massfractions of species
	for ( i = 0; i < nOfSpecies; ++i ) {
		fprintf( fp, "massfraction-%s\n", names[i] );
		for ( k = 0; k < gridPoints+2; ++k ) {
			fprintf( fp, "\t%-.6e", massFracs[k-1][i] );
			if ( (k+1) % 5 == 0 ) {
				fprintf( fp, "\n" );
			}
		}
		if ( k % 5 ) {
			fprintf( fp, "\n" );
		}
	}

	if ( flame->fPrintMolarFractions ) {
		Double	locMolarMass;
		for ( i = 0; i < nOfSpecies; ++i ) {
			locMolarMass = molarMass[i];
			fprintf( fp, "molefraction-%s\n", names[i] );
			fprintf( fp, "\t%-.6e", massFracs[kPrev][i] * mixMolarMass[-1] / locMolarMass );
			for ( k = 0; k < gridPoints; ++k ) {
				fprintf( fp, "\t%-.6e", massFracs[k][i] * mixMolarMass[k] / locMolarMass );
				if ( (k+2) % 5 == 0 ) {
					fprintf( fp, "\n" );
				}
			}
			fprintf( fp, "\t%-.6e\n", massFracs[gridPoints][i] * mixMolarMass[gridPoints] / locMolarMass );
		}
	}

	if ( flame->fSoot ) {
		flame->GetSoot()->PrintFlameletFile( gridPoints, flame, fp );
	}
	
//	write V
	Double	coeff = -sqrt( flame->GetStrainRate() * flame->fFlameNode->rhoInf * flame->fFlameNode->viscosityInf );
	fprintf( fp, "V [kg/(s*m^2)]\n" );
	fprintf( fp, "\t%-.6e", coeff * V[kPrev] );
	for ( k = 0; k < gridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", coeff * V[k] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", coeff * V[gridPoints] );

//	write density
	flame->PrintFlameletVector( gridPoints+2, &rho[kPrev], "density [kg/m^3]", fp );

//	write heat capacity
	flame->PrintFlameletVector( gridPoints+2, &props->GetHeatCapacity()->vec[kPrev], "cp [J/m^3K]", fp );

//	write heat conductivity
	flame->PrintFlameletVector( gridPoints+2, &props->GetConductivity()->vec[kPrev], "lambda [W/m]", fp );

//	write viscosity
	flame->PrintFlameletVector( gridPoints+2, &props->GetViscosity()->vec[kPrev], "mu [kg/s^2m]", fp );
	
//	write mDot_F
	Double	**mDot = flame->GetSpecies()->GetProductionRate()->mat;
	fprintf( fp, "mDot_%s [kg/m^3s]\n", names[flame->GetFuelIndex()] );
	fprintf( fp, "\t%-.6e", 0.0 );
	for ( k = 0; k < gridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", mDot[k][flame->GetFuelIndex()] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );

//	write molar mass
	flame->PrintFlameletVector( gridPoints+2, &mixMolarMass[kPrev], "MolarMass [kg/kmole]", fp );

//	write enthalpy
	fprintf( fp, "TotalEnthalpy [J/kg]\n" );
	int		nSpeciesIn = flame->GetSpecies()->GetNSpeciesInSystem();
	Double	hTot = 0.0;
	Double	**ent = flame->GetSpecies()->GetEnthalpy()->mat;
	for ( k = 0; k < gridPoints+2; ++k ) {
		hTot = 0.0;
		for ( i = 0; i < nSpeciesIn; ++i ) {
			hTot += massFracs[k-1][i] * ent[k-1][i];
		}
		fprintf( fp, "\t%-.6e", hTot );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}

	int	nOfSpeciesIn = flame->GetSpecies()->GetNSpeciesInSystem();

	fprintf( fp, "trailer\n" );
// write Le numbers at maximum consumption
	fprintf( fp, "Lewis numbers evaluated at maximum consumption\n" );
	int	maxLoc;
	int	off = flame->GetSpecies()->GetProductionRate()->phys_rows;
	Double	**prodRates = flame->GetSpecies()->GetProductionRate()->mat;
	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		flame->CompLewis( i, LewisFuel );
		maxLoc = LocationOfMin( gridPoints, &prodRates[0][i], off );
		if ( prodRates[maxLoc][i] < 0.0 ) {

			fprintf( fp, "%s\t%g\n", names[i], LewisFuel[maxLoc] );
		}
		else{
			fprintf( fp, "%s\t-1\n", names[i] );
		}
	}
// write Le numbers at maximum formation
	fprintf( fp, "Lewis numbers evaluated at maximum formation\n" );
	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		flame->CompLewis( i, LewisFuel );
		maxLoc = LocationOfMax( gridPoints, &prodRates[0][i], off );
		if ( prodRates[maxLoc][i] > 0.0 ) {

			fprintf( fp, "%s\t%g\n", names[i], LewisFuel[maxLoc] );
		}
		else{
			fprintf( fp, "%s\t-1\n", names[i] );
		}
	}
// write Le numbers at maximum mass fraction
	fprintf( fp, "Lewis numbers evaluated at maximum mass fraction\n" );
	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		flame->CompLewis( i, LewisFuel );
		maxLoc = LocationOfMax( gridPoints, &massFracs[0][i], flame->GetMassFracs()->phys_rows );
		fprintf( fp, "%s\t%g\n", names[i], LewisFuel[maxLoc] );
	}

	if ( species->IsConstantLewisNumber() ) {
		fprintf( fp, "Lewis numbers used in current calculation\n" );
		Double	*Le = species->GetLewisNumber()->vec;
		for ( i = 0; i < nOfSpecies; ++i ) {
			fprintf( fp, "%s\t%g\n", names[i], Le[i] );
		}
	}
	if ( nOfEquations < nOfVariables) {
		fprintf( fp, "number of converged equations is %d\n", nOfEquations );
	}
	
	DisposeVector( physXVec );
	
	if ( fpOpen ) {
		fclose( fp );
	}
}

void TCountPremFlameSim::CompLewis( int which, Double *Le )
{
	int				k;
	int				nOfGridPoints = fSolver->bt->GetGrid()->GetCurrentGrid()->GetNGridPoints();
	Double			*lambda = fProperties->GetConductivity()->vec;
	Double			*density = fProperties->GetDensity()->vec;
	Double			*cp = fProperties->GetHeatCapacity()->vec;
#ifdef FULLDIFFUSION
	Double			***D = fSpecies->GetBinDij()->tensor;
	int				ind_N2 = fSpecies->FindSpecies( "N2" );
	if ( ind_N2 < 0 ) {
		fprintf( stderr, "#error: no Species N2 in function 'TCountPremFlameSim::CompLewis'\n" ); 
		return;
	}
#else
	Double			**D = fSpecies->GetDiffusivity()->mat;
#endif
	
	for ( k = -1; k <= nOfGridPoints; ++k ) {
#ifdef FULLDIFFUSION
		Le[k] = lambda[k] / ( density[k] * cp[k] * D[k][which][ind_N2] );
#else
		Le[k] = lambda[k] / ( density[k] * cp[k] * D[k][which] );
		if ( isnan( Le[k] ) ) {
			fprintf( stderr, "Le[k] = %g\n", Le[k] );
			fprintf( stderr, "lam = %g\n", lambda[k] );
			fprintf( stderr, "den = %g\n", density[k] );
			fprintf( stderr, "cp = %g\n", cp[k] );
			fprintf( stderr, "D = %g\n", D[k][which] );
		}
#endif
	}
}

Double TCountPremFlameSim::NewDiffCorr( int nVariable, NodeInfoPtr nodeInfo )
{
// returns     sum_j ( d/dy(rho^2 Y_k D_j dY_j/dy) )

	int		i, k;
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	*Y = fFlameNode->Y[kCurr];
	Double	*YPrev = fFlameNode->Y[kPrev];
	Double	*YNext = fFlameNode->Y[kNext];
	Double	*density = fFlameNode->mixDensity;
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	k = nVariable-fFirstSpecies;
	Double	coeffCurr = density[kCurr] * density[kCurr] * Y[k];
	Double	coeffPrev = density[kPrev] * density[kPrev] * YPrev[k];
	Double	coeffNext = density[kNext] * density[kNext] * YNext[k];
	Double	diffPlus, diffMinus;
	Double	value = 0.0;
	Double	sumYCurr = 0.0;
	Double	sumYNext = 0.0;
	Double	sumYPrev = 0.0;
	

	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		sumYCurr += Y[i];
		sumYPrev += YPrev[i];
		sumYNext += YNext[i];
	}

	coeffCurr /= sumYCurr;
	coeffPrev /= sumYPrev;
	coeffNext /= sumYNext;
	
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		diffPlus = coeffCurr * diffusivity[i] + coeffNext * diffusivityNext[i];
		diffMinus = coeffCurr * diffusivity[i] + coeffPrev * diffusivityPrev[i];
		value += ( diffPlus * hm * ( YNext[i] - Y[i] ) 
					+ diffMinus * h * ( YPrev[i] - Y[i] ) );
	}

	return value / nodeInfo->hnenn;
}

Double TCountPremFlameSim::NewDiffCorrX( int nVariable, NodeInfoPtr nodeInfo )
{
// returns     sum_j ( d/dy(rho^2 Y_k / M * D_j Y_j dM/dy) )

	int		i, k;
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	MPrev = fFlameNode->mixMolarMass[kPrev];
	Double	M = fFlameNode->mixMolarMass[kCurr];
	Double	MNext = fFlameNode->mixMolarMass[kNext];
	Double	*Y = fFlameNode->Y[kCurr];
	Double	*YPrev = fFlameNode->Y[kPrev];
	Double	*YNext = fFlameNode->Y[kNext];
	Double	*density = fFlameNode->mixDensity;
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	k = nVariable-fFirstSpecies;
	Double	coeffCurr = density[kCurr] * density[kCurr] * Y[k] / M;
	Double	coeffPrev = density[kPrev] * density[kPrev] * YPrev[k] / MPrev;
	Double	coeffNext = density[kNext] * density[kNext] * YNext[k] / MNext;
	Double	diffPlus, diffMinus;
	Double	value = 0.0;
	Double	sumYCurr = 0.0;
	Double	sumYNext = 0.0;
	Double	sumYPrev = 0.0;
	
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		sumYCurr += Y[i];
		sumYPrev += YPrev[i];
		sumYNext += YNext[i];
	}
	coeffCurr /= sumYCurr;
	coeffPrev /= sumYPrev;
	coeffNext /= sumYNext;
	
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		diffPlus = coeffCurr * diffusivity[i] * Y[i] + coeffNext * diffusivityNext[i] * YNext[i];
		diffMinus = coeffCurr * diffusivity[i] * Y[i] + coeffPrev * diffusivityPrev[i] * YPrev[i];
		value += ( diffPlus * hm * ( MNext - M ) 
					+ diffMinus * h * ( MPrev - M ) );
	}

	return value / nodeInfo->hnenn;
}

void SetCountPremSimNodeInfo( int k, void *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	
	flame->SetFlameNode( k );
}

void CountPremSimPostConv( void *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	TNewtonPtr 				bt = flame->GetSolver()->bt;
	Double					*temp = flame->GetTemperature()->vec;
	TGridPtr				currentGrid = bt->GetGrid()->GetCurrentGrid();
	Double					*x = currentGrid->GetX()->vec;
	int						gridPoints = currentGrid->GetNGridPoints();
	int						isConverged = bt->GetConvergeNewton();
	Double					coeffSpecies = 1.0 / flame->GetStrainRate();
	Double					coeffTemp = 1.0 / flame->GetStrainRate();
	int						nSpeciesIn = flame->GetSpecies()->GetNSpeciesInSystem();
	Double					*Le = flame->GetSpecies()->GetLewisNumber()->vec;

//	fprintf( stderr, "Tmax = %g\n"
//			, temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1]);
	fprintf( stderr, "Tmax = %g and a = %g\n\n\n"
			, temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1], flame->GetStrainRate() );

	if ( isConverged ) {
		flame->SaveSolution();
		if ( isConverged && flame->CheckComputationalDomain() ) {
			flame->fSolver->ReInit();
			bt->WriteOutput( object, NULL, "tmp" );
			return;
		}	
	
		if ( flame->fPrintRHSSpecies ) {
			flame->PrintRHSSpecies( flame->GetSolver()->bt );
		}
		if ( flame->fPrintRHSTemp ) {
			flame->PrintRHSTemp( flame->GetSolver()->bt );
		}

		if ( flame->fSensAnal ) {
			flame->SensitivityAnalysis( coeffSpecies, coeffTemp, kSimilarity );
		}
		if ( flame->fReactionFluxes ) {
			flame->ReactionFluxes( kSimilarity );
			flame->GetReaction()->PrintReactionRates( flame );
			flame->fReaction->PrintRateCoeffs( flame );
			flame->fReaction->PrintDetailedHeatRelease( flame );
		}
		if ( flame->GetArcLengthCont() ) {
			char	tail[28];
			
			sprintf( tail, "Tmax%04.0f", temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
			bt->WriteOutput( object, NULL, tail );
			
			flame->IncNFlameletCount();
			if ( flame->GetDeltaArcLength() == 0.0 ) {
				flame->SetdlnStrainrateds( 1.0 );
				flame->SetdTds( 0.0 );
				flame->SetDeltaArcLength( flame->GetArcUp() ? flame->GetdlnStrainrateds() : -flame->GetdlnStrainrateds() );
			}
			else {
				Double dlnStrainrateds = ( log( flame->GetStrainRate() ) - flame->GetStrainrateContStart() ) 
										/ ( flame->GetDeltaStrainrateref() * flame->GetDeltaArcLength() );
				Double dTds = ( temp[flame->GetTMaxLoc()] - flame->GetTempContStart() ) 
										/ ( flame->GetDeltaTref() * flame->GetDeltaArcLength() );
				if ( flame->fStrainRateContin ) {
					dTds = 0.0;
				}
				flame->SetdlnStrainrateds( dlnStrainrateds );
				flame->SetdTds( dTds );
				Double ds = 0;
				Double absdlnStrainrateds = abs( dlnStrainrateds );
				Double absdTds = abs( dTds );
				ds = 0.5 * ( absdTds + absdlnStrainrateds );
				flame->SetDeltaArcLength( (flame->GetArcUp()) ? ds : -ds );
			}
			flame->SetTempContStart( temp[flame->GetTMaxLoc()] );
			flame->SetStrainrateContStart( log( flame->GetStrainRate() ) );
//				fprintf( stderr, "GetdTds = %g GetTempContStart = %g\n", flame->GetdTds(), flame->GetTempContStart() );
//				fprintf( stderr, "GetdlnStrainrateds = %g GetStrainrateContStart %g\n", flame->GetdlnStrainrateds(), flame->GetStrainrateContStart() );
			if ( flame->GetNFlameletsCount() < flame->GetMaxFlamelets() ) {
				flame->fSolver->ReInit();
//				fprintf( stderr, "\n\nFlamelet #%d: Tmax is now %g a is now %g\n", flame->GetNFlameletsCount()+1, flame->GetTempContStart(), exp(flame->GetStrainrateContStart()) );
				fprintf( stderr, "\n\nStart flamelet #%d with arclength %g\n", flame->GetNFlameletsCount()+1, flame->GetDeltaArcLength() );
			}
			else{
				fprintf( stderr, "\n\n%d flamelets have been computed.\n", flame->GetMaxFlamelets() );
				fprintf( stderr, "The computation is stopped here. Increase the value of the\n" );
				fprintf( stderr, "variable fMaxFlamelets in file $FlameManSource/TCountPremFlameSim.C\n" );
				fprintf( stderr, "if you want to compute more flamelets. You can also\n" );
				fprintf( stderr, "just start from one of the last solutions if it is not at a turning point in the strain rate\n" );
			}
		}
		else {
			if ( flame->fStrainRateContin ) {
				char	tail[28];
				
				sprintf( tail, "Tmax%04.0f", temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
				bt->WriteOutput( object, NULL, tail );
				
				
				int	inflectionPoint = LocationOfMaxSlope( temp, x, gridPoints );
				Double slope = ( temp[inflectionPoint] - temp[inflectionPoint-1] ) / (x[inflectionPoint] - x[inflectionPoint-1]);
				Double dx = -( temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] - temp[gridPoints] ) / slope;
				int i = inflectionPoint;
				while (x[i] < x[inflectionPoint]+dx) ++i;
				Double dT = 0.5*(temp[inflectionPoint] - temp[i]);
				
				
				if ( bt->GetGrid()->IsFine() && inflectionPoint%2 == 0) ++inflectionPoint;	// make it odd for the fine grid ( even gridpoints are removed for the coarse grid! )
				flame->SetTMaxLoc( inflectionPoint );
				flame->SetStrainrateContStart( x[inflectionPoint] );
				flame->SetTempContStart( temp[inflectionPoint]+dT );
				flame->fSolver->ReInit();
	//			fprintf(stderr, "TPoint at x = %g and ng = %d T = %g\tdT = %g\tdx = %g\n", x[inflectionPoint]
	//										, inflectionPoint, temp[inflectionPoint], dT, dx);
			}
			else {
				flame->PostConvergence( object );
				CountPremSimPostIter( flame );
				flame->SetStrainrateContStart( log( flame->T1DFlame::GetStrainRate() ) );
			}
		}
	}
	else {
		flame->RestoreSolution();
		CountPremSimPostIter( flame );

		if ( flame->GetArcLengthCont() ) {
			flame->SetDeltaArcLength( 0.5 * flame->GetDeltaArcLength() );
			
 			if ( abs( flame->GetDeltaArcLength() ) > 0.05 * ( abs( flame->GetdlnStrainrateds() ) + abs(flame->GetdTds() ) ) ) {		 
				flame->fSolver->ReInit();
				fprintf( stderr, "\n\nStart from Tmax = %g and strainrate %g with ds = %g\n", temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1], flame->GetStrainRate(), flame->GetDeltaArcLength() );
//				fprintf( stderr, "GetdTds = %g GetTempContStart = %g\n", flame->GetdTds(), flame->GetTempContStart() );
//				fprintf( stderr, "GetdlnStrainrateds = %g GetStrainrateContStart %g\n", flame->GetdlnStrainrateds(), flame->GetStrainrateContStart() );
				FILE *fp = flame->GetOutfile( "recoveredsol", TFlame::kData );
				bt->PrintSolution( x, currentGrid->GetY()->mat, flame->GetVariableNames(), fp );
				fclose(fp);
			}
			else{
				flame->ReInitArcLengthCont();
				flame->PostConvergence( object );
				CountPremSimPostIter( flame );
			}
		}
		else {
			flame->PostConvergence( object );
			CountPremSimPostIter( flame );
			flame->SetStrainrateContStart( log( flame->T1DFlame::GetStrainRate() ) );
		}
	}

	if ( !flame->GetArcLengthCont() && bt->GetLeaveContin() ) {
		Flag leaveContinPrem = flame->PostConvTPremixed( isConverged );
		if ( !leaveContinPrem ) {
			flame->SetMassFracsOfPhi( flame, flame->GetPhi()
					, &currentGrid->GetBcRight()->vec[flame->GetOffsetFirstSpecies()], flame->fSpecies->GetNSpeciesInSystem()
					, flame->fSpecies->GetMolarMass()->vec, flame->fSpecies->GetNames() );
			flame->SetMassFracsOfPhi( flame, flame->GetPhi()
					, &currentGrid->GetYLeft()->vec[flame->GetOffsetFirstSpecies()], flame->fSpecies->GetNSpeciesInSystem()
					, flame->fSpecies->GetMolarMass()->vec, flame->fSpecies->GetNames() );
			flame->fSolver->ReInit();
		}
	}
}

ConstStringArray GetCountPremSimVarNames( void *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	
	return flame->GetVariableNames();
}

FILE *TCountPremFlameSim::GetOutputFile( char *head, char *tail, FileType type )
{
	int				fuelIndex = GetFuelIndex();
	char			*name = new char[64];
	FILE			*fp;
	char			**speciesNames = fSpecies->GetNames();
	int				tOxidizer = ( int ) fSolTemp->vec[fSolTemp->len];
	int				tFuel = ( int ) fSolTemp->vec[kPrev];
	Double			press = GetPressure() * 1.0e-5;
	Double			phi = GetPhi();
		
	if ( fPremConfiguration == kBackToBack ) { // symmetric around stagnation plane, solve only right part
		sprintf( name, "%s%s%.8s_p%.2d_%.1dphi%.1d_%.4da%.5d_%.1dtu%.4d%s"
						, ( head ) ? head : "", ( head ) ? "_" : ""
						, speciesNames[fuelIndex]
						, ( int ) floor( press )	// in [bar]
						, ( int ) ( ( press - ( floor( press ) ) ) * 10 + 0.5 )
						, ( int ) floor( phi ) , ( int ) ( ( phi - floor( phi ) ) * 1.e4 )
						, ( int ) floor( GetStrainRate() + 1e-10 ) 
						, ( int ) ( ( GetStrainRate() - floor(GetStrainRate()+ 1e-10) ) * 10 + 0.5 )
						, ( int )( tOxidizer ) 						// in [K]
						, ( tail ) ? tail : "" );
	}
	else if ( fPremConfiguration == kFreshToBurned ) { // zero gradient for species and temp
		sprintf( name, "%s%s%.8s_p%.2d_%.1da%.5d_%.1dtf%.4dto%.4d%s"
						, ( head ) ? head : "", ( head ) ? "_" : ""
						, speciesNames[fuelIndex]
						, ( int ) floor( press )	// in [bar]
						, ( int ) ( ( press - ( floor( press ) ) ) * 10 + 0.5 )
						, ( int ) floor( GetStrainRate() + 1e-10 ) 
						, ( int ) ( ( GetStrainRate() - floor(GetStrainRate()+ 1e-10) ) * 10 + 0.5 )
						, ( int )( tFuel )							// in [K]
						, ( int )( tOxidizer ) 						// in [K]
						, ( tail ) ? tail : "" );
	}
	else if ( fPremConfiguration == kFull ) {
		sprintf( name, "%s%s%.8s_p%.2d_%.1da%.5d_%.1dtf%.4dto%.4d%s"
						, ( head ) ? head : "", ( head ) ? "_" : ""
						, speciesNames[fuelIndex]
						, ( int ) floor( press )	// in [bar]
						, ( int ) ( ( press - ( floor( press ) ) ) * 10 + 0.5 )
						, ( int ) floor( GetStrainRate() + 1e-10 ) 
						, ( int ) ( ( GetStrainRate() - floor(GetStrainRate()+ 1e-10) ) * 10 + 0.5 )
						, ( int )( tFuel )							// in [K]
						, ( int )( tOxidizer ) 						// in [K]
						, ( tail ) ? tail : "" );
	}
	else {
		fprintf( stderr, "#error in premixed configuration\n" );
		exit(2);
	}
	fp = GetOutfile( name, type );
	delete name;

	return fp;
}

void TCountPremFlameSim::EtaToX( TNewtonPtr bt, VectorPtr xPhysVec )
{
	int			k;
	int			gridPoints = bt->GetCurrentGridPoints();
	int			kBefStoech = -1;
	NodeInfoPtr nodeInfo = bt->GetNodeInfo();
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();

	Double		factor = 0.5 * sqrt( fFlameNode->rhoInf * fFlameNode->viscosityInf / GetStrainRate() );
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		*rho = GetProperties()->GetDensity()->vec;
	Double		*x = grid->GetX()->vec;
	Double		**y = grid->GetY()->mat;
	Double		*xPhys = xPhysVec->vec;
	Double		zStoech = fMassFraction->GetZStoe();
	Double		xStoech;
	
	if ( xPhysVec->len < gridPoints + 2 ) {
		cerr << "#warning: Vector xPhys too short, values for physical grid are not computed" << NEWL;
		return;
	}
	
	xPhys[0] = 0.0;
	xPhys[1] = xPhys[0] + factor * ( x[0] - left ) * ( 1.0 / rho[-1] + 1.0 / rho[0] );
	for ( k = 1; k < gridPoints; ++k ) {
		xPhys[k+1] = xPhys[k] + factor * ( x[k] - x[k-1] ) * ( 1.0 / rho[k-1] + 1.0 / rho[k] );
	}
	xPhys[gridPoints+1] = xPhys[gridPoints] + factor * ( right - x[gridPoints-1] ) * ( 1.0 / rho[gridPoints-1] + 1.0 / rho[gridPoints] );
}


void CountPremSimUpdateLeftBoundary( void  *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	TNewtonPtr		bt = flame->GetSolver()->bt;
	TGridPtr 		currGrid = bt->GetGrid()->GetCurrentGrid();
	MatrixPtr		yMat = currGrid->GetY();
	Double			**y = yMat->mat;
	Double			*yLeft = currGrid->GetYLeft()->vec;
	Double			pressure = flame->GetPressure();
	Double			*yFirst = yMat->mat[0];
	Double			hFirst = currGrid->GetX()->vec[0] - bt->GetLeft();
	int				nGridPoints = currGrid->GetNGridPoints();
	int				fUVelocity = flame->GetOffsetUVelocity();
	int				fVVelocity = flame->GetOffsetVVelocity();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();
	int			variables = bt->GetNVariables();
	T1DSpeciesPtr	species = flame->GetSpecies();
	T1DPropertiesPtr	properties = flame->GetProperties();
	Double		*diffusivity = species->GetDiffusivity()->mat[0];
	Double		*density =  properties->GetDensity()->vec;
	Double		&rhoLeft = density[-1];
	Double		&rhoRight = density[nGridPoints];
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();
	int			*bcFlagLeft = currGrid->GetBcFlagLeft();
	int			mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
	Double		*bcLeft = currGrid->GetBcLeft()->vec;
	Double		viscosityInf = properties->GetViscosity()->vec[nGridPoints];

	if ( flame->fPremConfiguration == kBackToBack ) { // symmetric around stagnation plane, solve only right part
		for ( int j = 0; j < variables; ++j ) {
			yLeft[j] = y[0][j];
		}
		yLeft[fVVelocity] = 0;
	}
	else if ( flame->fPremConfiguration == kFreshToBurned ) { // zero gradient for species and temp
		for ( int j = 0; j < variables; ++j ) {
			if (j != fVVelocity ) yLeft[j] = y[0][j];
		}
		yLeft[fUVelocity] = sqrt( rhoRight / rhoLeft );
	}
	else if ( flame->fPremConfiguration == kFull ) {
		yLeft[fUVelocity] = sqrt( rhoRight / rhoLeft );
#ifndef FCONVLEFTTORIGHT
		yLeft[fVVelocity] = yFirst[fVVelocity] - hFirst * yLeft[fUVelocity];
#endif
/* 		if ( mixtureSpecificationLeft == kMassFlux ) { */
/* 			Double 		coeff; */
/* 			for ( int i = 0; i < nSpeciesInSystem; ++i ) { */
/* 				if ( bcFlagLeft[fFirstSpecies] == kDirichlet ) { */
/* 					coeff = rhoLeft * diffusivity[-1] / ( yLeft[fVVelocity] * rhoRight * viscosityInf * hFirst ); */
/* 					yLeft[fFirstSpecies+i] = ( bcLeft[fFirstSpecies+i] - y[0][fFirstSpecies+i] * coeff )  */
/* 												/ ( 1.0 - coeff ); */
/* 				} */
/* 				else if ( bcFlagLeft[fFirstSpecies] != kNone ) { */
/* 					cerr << "error: i can't handle boundary condition of kind " << bcFlagLeft[fFirstSpecies] << " for species no. " << i << NEWL; */
/* 				} */
/* 			} */
/* 		} */
	}
	else {
		fprintf( stderr, "#error in premixed configuration\n" );
		exit(2);
	}

	flame->UpdateSolutionOnePoint( yLeft, -1 );
	flame->SetFlameNode( -1 );
	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], pressure );
}

void CountPremSimUpdateRightBoundary( void *object )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	int				fUVelocity = flame->GetOffsetUVelocity();
	int				fVVelocity = flame->GetOffsetVVelocity();
	TNewtonPtr		bt = flame->GetSolver()->bt;
	TGridPtr 		currGrid = bt->GetGrid()->GetCurrentGrid();
	int				nGridPoints = currGrid->GetNGridPoints();
	MatrixPtr		yMat = currGrid->GetY();
	Double			**y = yMat->mat;
	Double			*yRight = currGrid->GetYRight()->vec;
	Double			*yLast = y[nGridPoints-1];
	Double			*x = currGrid->GetX()->vec;
	Double			hLast = bt->GetRight() - x[nGridPoints-1];
	Double			pressure = flame->GetPressure();
	int				variables = bt->GetNVariables();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();
	T1DSpeciesPtr	species = flame->GetSpecies();
	T1DPropertiesPtr	properties = flame->GetProperties();
	Double		*diffusivity = species->GetDiffusivity()->mat[0];
	Double		*density =  properties->GetDensity()->vec;
	Double		&rhoRight = density[nGridPoints];
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();
	int			*bcFlagRight = currGrid->GetBcFlagRight();
	int			mixtureSpecificationRight = flame->GetMixtureSpecificationRight();
	Double		*bcRight = currGrid->GetBcRight()->vec;
	Double		viscosityInf = properties->GetViscosity()->vec[nGridPoints];


// right boundary 
	if ( flame->fPremConfiguration == kBackToBack ) { // symmetric around stagnation plane, solve only right part
		yRight[fUVelocity] = 1.0;
		yRight[fVVelocity] = yLast[fVVelocity] + hLast * yRight[fUVelocity];
	}
	else if ( flame->fPremConfiguration == kFreshToBurned ) { // zero gradient for species and temp
		yRight[fUVelocity] = 1.0;
		yRight[fVVelocity] = yLast[fVVelocity] + hLast * yRight[fUVelocity];
	}
	else if ( flame->fPremConfiguration == kFull ) {
		yRight[fUVelocity] = 1.0;
#ifdef FCONVLEFTTORIGHT
		yRight[fVVelocity] = yLast[fVVelocity] + hLast * yRight[fUVelocity];
#endif
	}
	else {
		fprintf( stderr, "#error in premixed configuration\n" );
		exit(2);
	}

	if ( mixtureSpecificationRight == kMassFlux ) {
		Double 		coeff;
		for ( int i = 0; i < nSpeciesInSystem; ++i ) {
			if ( bcFlagRight[fFirstSpecies] == kDirichlet ) {
				coeff = rhoRight * diffusivity[nGridPoints] / ( yRight[fVVelocity] * rhoRight * viscosityInf * hLast );
				yRight[fFirstSpecies+i] = ( bcRight[fFirstSpecies+i] - y[nGridPoints-1][fFirstSpecies+i] * coeff ) 
											/ ( 1.0 - coeff );
				yRight[fFirstSpecies+i] = bcRight[fFirstSpecies+i]; 
			}
			else if ( bcFlagRight[fFirstSpecies] != kNone ) {
				cerr << "error: i can't handle boundary condition of kind " << bcFlagRight[fFirstSpecies] << " for species no. " << i << NEWL;
			}
		}
	}

	flame->UpdateSolutionOnePoint( yRight, nGridPoints );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], pressure );
}


Double TCountPremFlameSim::GetBurningVelocity( void )
{
	Double		*V = GetV()->vec;
	Double		*temp = GetTemperature()->vec;
	TGridPtr	currentGrid = GetSolver()->bt->GetGrid()->GetCurrentGrid();
	Double		*x = currentGrid->GetX()->vec;
	int			gridPoints = currentGrid->GetNGridPoints();
	int			inflectionPoint = LocationOfMaxSlope( temp, x, gridPoints );
	Double		rho_u = GetProperties()->GetDensity()->vec[gridPoints-1];
	Double		coeff = sqrt( GetStrainRate() * fFlameNode->rhoInf * fFlameNode->viscosityInf );

	Double		velocity = coeff * V[inflectionPoint] / rho_u;
	
//	fprintf( stderr, "T0 = %g\tsL = %g\n", temp[inflectionPoint], velocity );

	return velocity;
}

void CountPremSimCheckCompDomain( void *object, NodeInfoPtr /*nodeInfo*/ )
{
	TCountPremFlameSimPtr	flame = ( TCountPremFlameSimPtr )object;
	flame->CheckComputationalDomain();
}

int TCountPremFlameSim::CheckComputationalDomain( void )
{
	int			changed = 0;
	int				last;
	TBVPSolverPtr	solver = GetSolver();
	TNewtonPtr	bt = solver->bt;
	TGridPtr 	grid = bt->GetGrid()->GetFine();
	VectorPtr	xVec = grid->GetX();
	Double		*x = xVec->vec;
	Double		**y = grid->GetY()->mat;
	int			nGridPoints = grid->GetNGridPoints();
	Double		*temp = fSolTemp->vec;
	Double		right =  bt->GetRight();
	Double		left =  bt->GetLeft();
	Double		*yLeft = grid->GetYLeft()->vec;
	Double		*yRight = grid->GetYRight()->vec;
	Double		slope;
	Double		slope_min = 1.0e-6;	
	Double		slope_max = 1.0e-3;	
	Double		alpha = 0.1;
	Double		deltaX, deltaleft;
		
	// check right boundary
	slope = fabs( FirstDerivUpwind( temp[nGridPoints-1], temp[nGridPoints-2], x[nGridPoints-1] - x[nGridPoints-2]) / temp[nGridPoints-1] );
//	fprintf( stderr, "sloperight = %g\n", slope );
	if ( slope > slope_max ) {
		if ( !bt->GetGrid()->GetOneSolOneGrid() ) {
			cerr << "#warning: cannot enlarge computational domain for equidistant grid" << NEWL;
			changed = MAX(changed, 0);
		}
		else {
			deltaX = right - x[nGridPoints-1];
			right += ( right - left ) * alpha;;
			bt->SetRight( right );
			yRight[fVVelocity] = y[nGridPoints-1][fVVelocity] + 
								 ( yRight[fVVelocity] - y[nGridPoints-1][fVVelocity] ) / deltaX *
								 ( right - x[nGridPoints-1] );
			cerr << "enlarge the computational domain by " << alpha * 100.0 << " percent on the right" << NEWL << NEWL;
			changed = MAX(changed, 1);
		}
	}
	else if ( slope < slope_min ) {
		if ( !bt->GetGrid()->GetOneSolOneGrid() ) {
			cerr << "#warning: cannot cut computational domain for equidistant grid" << NEWL;
			changed = MAX(changed, 0);
		}
		else {
			right -= ( right - left ) * alpha;
			bt->SetRight( right );
			for ( last = nGridPoints - 1; x[last] > right; --last );
			grid->AdjustNGridPoints( last );
			solver->UpdateAllDimensions( last );
			CountPremSimPostIter( this );
			cerr << "cut the computational domain by " << alpha * 100.0 << " percent" << NEWL << NEWL;
			changed = MAX(changed, 1);
		}
	}

	if ( fPremConfiguration != kBackToBack ) {
		// check left boundary
		slope = fabs( FirstDerivUpwind( temp[1], temp[0], x[1] - x[0]) / temp[0] );
//		fprintf( stderr, "slopeleft = %g\n", slope );
		if ( slope > slope_max ) {
			if ( !bt->GetGrid()->GetOneSolOneGrid() ) {
				cerr << "#warning: cannot enlarge computational domain for equidistant grid" << NEWL;
				changed = MAX(changed, 0);
			}
			else {
				deltaX = x[0] - left;
				left -= ( right - left ) * alpha;
				bt->SetLeft( left );
				yLeft[fVVelocity] = y[0][fVVelocity] + 
									 ( yLeft[fVVelocity] - y[0][fVVelocity] ) / deltaX *
									 ( x[0] - left );
				cerr << "enlarge the computational domain by " << alpha * 100.0 << " percent on the left" << NEWL << NEWL;
				changed = MAX(changed, 1);
			}
		}
		else if ( slope < slope_min ) {
			if ( !bt->GetGrid()->GetOneSolOneGrid() ) {
				cerr << "#warning: cannot cut computational domain for equidistant grid" << NEWL;
				changed = MAX(changed, 0);
			}
			else {
				left += ( right - left ) * alpha;
				bt->SetLeft( left );
				for ( last = 0; x[last] < left; ++last );
				grid->AdjustNGridPoints( nGridPoints-last+1 );
				solver->UpdateAllDimensions( last );
				CountPremSimPostIter( this );
				cerr << "cut the computational domain by " << alpha * 100.0 << " percent" << NEWL << NEWL;
				changed = MAX(changed, 1);
			}
		}
	}
	
	return changed;
}

void TCountPremFlameSim::SetInitialBC( TGridPtr grid, TInputDataPtr inp )
{
	int					nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	BoundaryInputPtr	right = inp->rightBoundary;
	int					rightSpecifiedSpecies = right->fSpecifiedSpeciesBCs;
	Double				*yleft = grid->GetYLeft()->vec;
	Double				*bcleft = grid->GetBcLeft()->vec;
	Double				*bcRight = grid->GetBcRight()->vec;
	Double				*yright = grid->GetYRight()->vec;
	Double				*phi = GetPhiVector()->vec;

	if ( phi[0] > 0 ) {
// right
		SetMassFracsOfPhi( this, phi[0], &bcRight[fFirstSpecies]
			, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		SetMassFracsOfPhi( this, phi[0], &yright[fFirstSpecies]
			, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		SetMixtureSpecificationRight( kMassFlux );

// left only for full
		if ( fPremConfiguration == kFull ) {
			SetMassFracsOfPhi( this, phi[0], &bcleft[fFirstSpecies]
				, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
			SetMassFracsOfPhi( this, phi[0], &yleft[fFirstSpecies]
				, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
			SetMixtureSpecificationLeft( kMassFlux );
		}
	}
	else if ( rightSpecifiedSpecies > 0 ) {
		SetPhiOfMassFracs( this, phi, &bcRight[fFirstSpecies], fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		fprintf( stderr, "phi is %g\n", phi[0] );
	}
	else {
		cerr << "#error: boundary conditions for the unburnt missing!" << NEWL;
		exit( 2 );
	}
}

int TCountPremFlameSim::GetTMaxLoc( void ) 
{ 
	if ( fArcLengthContin || !fStrainRateContin ) {
		return fTmaxLoc; 
	}
	else {
		int i = 0;
		Double *x = fSolver->bt->GetGrid()->GetCurrentGrid()->GetX()->vec;
		while (x[i] < fLnStrainrateContStart) ++i;
	
		return i;
	}
}
