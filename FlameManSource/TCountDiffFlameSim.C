#include "FlameMaster.h"
#include "ListTool.h"
#include "Spline.h"
#include "TCountDiffFlameSim.h"


#define UPWINDCONVECTION

#undef FULLDIFFUSION

void TCountDiffFlameSim::InitCountDiffFlameSim( void )
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
		cerr << "use initial scalar dissipation rate from command line: a = " << GetStrainRate() << NEWL;
	}

	fVariableNames = new String[fVariablesWithoutSpecies + nSpeciesInSystem];

	fVariableNames[fVVelocity] = new char[2];
	strcpy( fVariableNames[fVVelocity], "f" );
	fVariableNames[fUVelocity] = new char[3];
	strcpy( fVariableNames[fUVelocity], "f'" );
	fVariableNames[fMixFrac] = new char[2];
	strcpy( fVariableNames[fMixFrac], "Z" );
	fVariableNames[fTemperature] = new char[2];
	strcpy( fVariableNames[fTemperature], "T" );
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies->GetNames()[i] ) + 1];
		strcpy( fVariableNames[fFirstSpecies + i], fSpecies->GetNames()[i] );
	}
	if ( fSoot ) {
		int	offset = fSoot->GetOffsetSootMoments();
		for ( i = 0; i < fSoot->GetNSootMoments(); ++i ) {
			fVariableNames[offset + i] = new char[8];
			sprintf( fVariableNames[offset + i], "M%d/rho", i );
		}
	}
	fPrintMolarFractions = fInputData->fPrintMolarFractions;
	
//	vectors of solution
	fSolV = NewVector( maxGridPoints + 2 );
	fSolU = NewVector( maxGridPoints + 2 );
	fSolMixFrac = NewVector( maxGridPoints + 2 );

	fSolV->vec = &fSolV->vec[kNext];
	fSolU->vec = &fSolU->vec[kNext];
	fSolMixFrac->vec = &fSolMixFrac->vec[kNext];

	fSolV->len -= 2;
	fSolU->len -= 2;
	fSolMixFrac->len -= 2;

//	saved solution
	fSavedV = NewVector( maxGridPoints + 2 );
	fSavedU = NewVector( maxGridPoints + 2 );
	fSavedMixFrac = NewVector( maxGridPoints + 2 );

	fSavedV->vec = &fSavedV->vec[kNext];
	fSavedU->vec = &fSavedU->vec[kNext];
	fSavedMixFrac->vec = &fSavedMixFrac->vec[kNext];

	fSavedV->len -= 2;
	fSavedU->len -= 2;
	fSavedMixFrac->len -= 2;

	bt->SetUtFuncs( CountDiffSimJacRest, CountDiffSimJacRest, CountDiffSimJacRest
					, CountDiffSimRHSRest, CountDiffSimRHSRest, CountDiffSimRHSRest
					, CountDiffSimOutput, CountDiffSimPostIter
					, SetCountDiffSimNodeInfo, CountDiffSimPostConv
					, GetCountDiffSimVarNames
					, CountDiffSimUpdateLeftBoundary, CountDiffSimUpdateRightBoundary );
	SetInitialBC( fine, fInputData );
	SetInitialBC( coarse, fInputData );
	fMassFraction = new TMassFraction( this );
	if ( !fMassFraction ) FatalError( "memory allocation of TMassFraction failed" );
	//	fMassFraction->PrintMassFraction();

	ReadStartProfiles( fInputData );
	CheckBC();
	CheckInitialGuess();
	UpdateSolution( fine->GetY(), fine->GetYLeft(), fine->GetYRight() );
}

TCountDiffFlameSim::~TCountDiffFlameSim( void )
{
	int	nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();

	delete fMassFraction;

	fSavedMixFrac->vec = &fSavedMixFrac->vec[kPrev];
	fSavedU->vec = &fSavedU->vec[kPrev];
	fSavedV->vec = &fSavedV->vec[kPrev];

	DisposeVector( fSavedMixFrac );
	DisposeVector( fSavedU );
	DisposeVector( fSavedV );

	fSolMixFrac->vec = &fSolMixFrac->vec[kPrev];
	fSolV->vec = &fSolV->vec[kPrev];
	fSolU->vec = &fSolU->vec[kPrev];

	DisposeVector( fSolMixFrac );
	DisposeVector( fSolU );
	DisposeVector( fSolV );

	for ( int i = 0; i < nSpeciesInSystem+fVariablesWithoutSpecies; ++i ) {
		delete fVariableNames[i];
	}
	delete fVariableNames;
}

void CountDiffSimJacRest( void *object, NodeInfoPtr nodeInfo )
{
	TCountDiffFlameSimPtr	flame = ( TCountDiffFlameSimPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		M = nodeInfo->nOfEquations;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
    Double  hnenn = nodeInfo->hnenn;
	Double	**a = nodeInfo->a;
	Double	*temp = flameNode->temp;
	Double	*y = nodeInfo->y;
	Double	strainRate = flame->GetStrainRate();
	Double	mixDensity = *flameNode->mixDensity;
	Double	mixHeatCapacity = *flameNode->mixHeatCapacity;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	*enthalpy = flameNode->enthalpy;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	*productionRate = flameNode->productionRate;
	Double	mixDensityInf = flameNode->rhoInf;
	Double	mixViscosityInf = flameNode->viscosityInf;
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	constThermDiffCoeff = constMassDiffCoeff / mixHeatCapacity;
	Double	oneOverATimesRho = 1.0 / ( strainRate * mixDensity );
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
//	Double	ROverAP = RGAS / ( strainRate * flame->GetPressure() );
//	Double	ROverAPM = ROverAP / mixMolarMass;
//	Double	energySourceCoeff = ROverAPM / mixHeatCapacity;
	Double	energySourceCoeff = oneOverATimesRho / mixHeatCapacity;
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

	if ( flame->GetSoot() ) {
		if ( !flame->fUseNumericalDM )  {
			flame->GetSoot()->UpdateJacobian( flame );
		}
		flame->GetSoot()->FillJacobi( flame, nodeInfo, kSimilarity );
	}

// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );

#ifdef UPWINDCONVECTION
	// second equation ( momentum )
	FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0, FALSE );
	// third equation ( mixturefraction )
	FillJacNonlinearConvectUpwind( fVVelocity, fMixFrac, nodeInfo, 1.0, FALSE );
	// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0, FALSE );
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0, FALSE );
	}
#else
	// second equation ( momentum )
	FillJacNonlinearConvectCentral( fVVelocity, fUVelocity, nodeInfo );
	// third equation ( mixturefraction )
	FillJacNonlinearConvectCentral( fVVelocity, fMixFrac, nodeInfo );
	// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectCentral( fVVelocity, speciesEq, nodeInfo );
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectCentral( fVVelocity, fTemperature, nodeInfo );
	}
#endif

	
// first equation ( mass )
	a[fUVelocity][fVVelocity] -= hnenn;
	
// second equation ( momentum )
	flame->FillJacDiffusion( fUVelocity, fUVelocity, constMassDiffCoeff, flameNode->mixViscosity, nodeInfo );
	a[fUVelocity][fUVelocity] -= 2.0 * y[fUVelocity] * hnenn;
	a[fTemperature][fUVelocity] += mixDensityInf / ( mixDensity * temp[kCurr] ) * hnenn;
	
	Double	theConst = mixDensityInf * mixMolarMass / mixDensity * hnenn;
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq; ++speciesEq ) {
		a[speciesEq][fUVelocity] += theConst / molarMass[speciesEq-fFirstSpecies];
	}
	
// third equation ( mixturefraction )
	flame->FillJacMixFracDiffusion( fMixFrac, fMixFrac, constMassDiffCoeff, nodeInfo );

// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, constMassDiffCoeff, nodeInfo );
		if ( flame->UseDiffCorr() ) {
			flame->FillJacDiffCorr( speciesEq, constMassDiffCoeff, nodeInfo, kNegative );
		}
		//	implicit source term
		for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq && speciesVar < M; ++ speciesVar ) {
			speciesIndexVar = speciesVar - fFirstSpecies;
			a[speciesVar][speciesEq] += ( dMdY[speciesIndexVar][speciesIndexEq]
					+ productionRate[speciesIndexEq] * mixMolarMass / molarMass[speciesIndexVar] )
					* oneOverATimesRho * hnenn;
		}
		a[fTemperature][speciesEq] += productionRate[speciesIndexEq] / y[fTemperature]
										* oneOverATimesRho * hnenn;
		a[fTemperature][speciesEq] += dMdT[speciesIndexEq]
										* oneOverATimesRho * hnenn;
	}

// last equation ( temperature )
	if ( fTemperature < M ) {
		int 	i, l;
		flame->FillJacDiffusion( fTemperature, fTemperature, constThermDiffCoeff, flameNode->mixConductivity, nodeInfo );

		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			//	d/dT sum (m_i * h_i ) 
			a[fTemperature][fTemperature] -= ( dMdT[i] * enthalpy[i] + 
						productionRate[i] * ( enthalpy[i] / temp[kCurr] + heatCapacity[i] ) ) 
						* energySourceCoeff * hnenn;
			for ( l = 0; l < nSpeciesInSystem; ++l ) {
				//	d/dY_l sum (m_i * h_i ) 
				a[l + fFirstSpecies][fTemperature] -= ( dMdY[l][i]
							+ productionRate[i] * mixMolarMass / molarMass[l] ) 
							* enthalpy[i] * energySourceCoeff * hnenn;
			}
		}
		
		if ( flame->fProperties->GetRadiation() ) {
			flame->fProperties->GetRadiation()->FillJacRadiation( energySourceCoeff, flame, nodeInfo );
		}
	}

	if ( flame->GetSolver()->bt->GetTimedepFlag() 
				&& !flame->GetSolver()->time->GetTimeConverged() ) {
		TTimePtr tim = flame->GetSolver()->time;
		for ( int eqLoop = 0; eqLoop < M; ++eqLoop ) {
			if ( eqLoop != fVVelocity ) {
				a[eqLoop][eqLoop] -= hnenn / tim->GetDeltaT();
			}
		}
	}
}

void CountDiffSimRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
	TCountDiffFlameSimPtr	flame = ( TCountDiffFlameSimPtr )object;
	if ( !flame->RHSAction( nodeInfo, rhsMode ) ) {
		return;
	}
	
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
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
	Double	strainRate = flame->GetStrainRate();
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
	rhs[fVVelocity] += FirstDerivUpwind( y[fVVelocity], yPrev[fVVelocity], hm );

#ifdef UPWINDCONVECTION
	// second equation ( momentum )
	rhs[fUVelocity] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fUVelocity], y[fUVelocity], yNext[fUVelocity], hm, h, FALSE );
	// third equation ( mixturefraction )
	rhs[fMixFrac] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fMixFrac], y[fMixFrac], yNext[fMixFrac], hm, h, FALSE );
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
	// third equation ( mixturefraction )
	rhs[fMixFrac] += NonlinearConvectCentral( y[fVVelocity], yPrev[fMixFrac], y[fMixFrac], yNext[fMixFrac], hm, h );
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
	
// mixture fraction equation
	rhs[fMixFrac] += constMassDiffCoeff * flame->SecondDerivMixFracDiffusion( fMixFrac, nodeInfo );

// fFirstSpecies to fFirstSpecies + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		eqLoop = speciesEq - fFirstSpecies;
#ifdef FULLDIFFUSION
		rhs[speciesEq] -= constMassDiffCoeff * flame->SecondDerivBinSpecDiff( speciesEq, nodeInfo );
#else
		rhs[speciesEq] += constMassDiffCoeff * flame->SecondDerivSpeciesDiffusion( speciesEq, nodeInfo );
		rhs[speciesEq] += constMassDiffCoeff * flame->SecondDerivXDiffusion( speciesEq, nodeInfo );
#endif

		if ( flame->UseDiffCorr() ) {
			rhs[speciesEq] -= constMassDiffCoeff * flame->NewDiffCorr( speciesEq, nodeInfo );
			rhs[speciesEq] -= flame->GetSolver()->bt->GetParameter() * constMassDiffCoeff * flame->NewDiffCorrX( speciesEq, nodeInfo );
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

		rhs[fTemperature] += constThermDiffCoeff * flame->StandardDiffusion( fTemperature, flameNode->mixConductivity, nodeInfo );

		for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
			sumCpDdYdx += ( heatCapacity[eqLoop] - diffCorrHeat ) * diffusivity[eqLoop] 
						* FirstDeriv( YPrev[eqLoop], Y[eqLoop], YNext[eqLoop], hm, h );
// check_ if next statement needs to be disabled for diffusivity correction switched off
			sumCpYD += ( heatCapacity[eqLoop] - diffCorrHeat ) * Y[eqLoop] * diffusivity[eqLoop];
			sumMH += productionRate[eqLoop] * enthalpy[eqLoop];
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
		
		if ( flame->fProperties->GetRadiation() ) {
			rhs[fTemperature] += flameNode->radiation[kCurr] * oneOverARhoCp;
			if ( flame->GetSoot() ) {
				rhs[fTemperature] -= oneOverARhoCp * flame->GetSoot()->GetSootRadiation( y[fTemperature], flameNode->moments );
			}
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

Double TCountDiffFlameSim::SecondDerivXMixFracDiff( int nVariable, NodeInfoPtr nodeInfo )
{
	Double	*W = fFlameNode->mixMolarMass;
	Double	*lambda = fFlameNode->mixConductivity;
	Double	*rho = fFlameNode->mixDensity;
	Double	*cp = fFlameNode->mixHeatCapacity;
	Double	diffPlus =  rho[kCurr] * lambda[kCurr] / cp[kCurr] / W[kCurr]
					+ rho[kNext] * lambda[kNext] / cp[kNext] / W[kNext];
	Double	diffMinus = rho[kPrev] * lambda[kPrev] / cp[kPrev] / W[kPrev]
					+ rho[kCurr] * lambda[kCurr] / cp[kCurr] / W[kCurr];
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	yPrev = nodeInfo->yPrev[nVariable];
	Double	y = nodeInfo->y[nVariable];
	Double	yNext = nodeInfo->yNext[nVariable];
	
	return ( diffPlus * hm * yNext * W[kNext] - ( diffPlus * hm + diffMinus * h ) * y * W[kCurr] 
				+ diffMinus * h * yPrev * W[kPrev] ) 
				/ ( h * hm * ( h + hm ) );
}

int CountDiffSimPostIter( void *object )
{
	TCountDiffFlameSimPtr	flame = ( TCountDiffFlameSimPtr )object;
	int			fUVelocity = flame->GetOffsetUVelocity();
	int			fVVelocity = flame->GetOffsetVVelocity();
	int			fTemperature = flame->GetOffsetTemperature();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();
	int			fMixFrac = flame->GetOffsetMixFrac();
	int			i;
	TNewtonPtr	bt = flame->GetSolver()->bt;
	NodeInfoPtr nodeInfo = bt->GetNodeInfo();
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = currGrid->GetNGridPoints();
	int			*bcFlagLeft = currGrid->GetBcFlagLeft();
	int			mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
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
	
// first set temperature and massfractions 
	for ( i = 0; i < nGridPoints; ++i ) {
		bt->SetNodeInfo( flame, i );		
		flame->fMassFraction->UpdateMixFracDependence( flame, nodeInfo, TRUE );
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

//	update properties
	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, temp[kPrev], Y[kPrev], pressure );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, temp[nGridPoints], Y[nGridPoints], pressure );

// left boundary
	if ( mixtureSpecificationLeft == kMassFlux ) {
		Double 		coeff;
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			if ( bcFlagLeft[fFirstSpecies] == kDirichlet ) {
				coeff = rhoLeft * diffusivity[-1] / ( yLeft[fVVelocity] * rhoRight * viscosityInf * hFirst );
				yLeft[fFirstSpecies+i] = ( bcLeft[fFirstSpecies+i] - y[0][fFirstSpecies+i] * coeff ) 
											/ ( 1.0 - coeff );
			}
			else {
				cerr << "error: i can't handle boundary condition of kind " << bcFlagLeft[fFirstSpecies] << " for species no. " << i << NEWL;
			}
		}
	}
	
	if ( mixtureSpecificationLeft == kMassFlux ) {	// means tsuji bc
		yLeft[fUVelocity] = sqrt( rhoRight / rhoLeft ) * ( 1.0 - parameter );
	}
	else {
		yLeft[fUVelocity] = sqrt( rhoRight / rhoLeft );
	}
	flame->fMassFraction->SetMassFractionsBC( fMixFrac, fFirstSpecies, yLeft, yRight );

	
// right boundary 
	yRight[fUVelocity] = 1.0;
	yRight[fVVelocity] = yLast[fVVelocity] + hLast * yRight[fUVelocity];

	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	
	flame->UpdateThermoProps();


	return 0;
}

#include "TofZ.h"

void TCountDiffFlameSim::UpdateDimensions( int len )
{
	T1DFlame::UpdateDimensions( len );
	fSolV->len = len;
	fSolU->len = len;
	fSolMixFrac->len = len;
}

void TCountDiffFlameSim::UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
{
	int		nGridPoints = yMat->cols;
	Double	*V = fSolV->vec;
	Double	*U = fSolU->vec;
	Double	*Z = fSolMixFrac->vec;
	Double	**y = yMat->mat;
	Double	*yLeft = yLeftVec->vec;
	Double	*yRight = yRightVec->vec;

	UpdateDimensions( nGridPoints );

	T1DFlame::UpdateSolution( yMat, yLeftVec, yRightVec );
	
	V[kPrev] = yLeft[fVVelocity];
	U[kPrev] = yLeft[fUVelocity];
	Z[kPrev] = yLeft[fMixFrac];
	for ( int k = 0; k < nGridPoints; ++k ) {
		V[k] = y[k][fVVelocity];
		U[k] = y[k][fUVelocity];
		Z[k] = y[k][fMixFrac];
	}
	V[nGridPoints] = yRight[fVVelocity];
	U[nGridPoints] = yRight[fUVelocity];
	Z[nGridPoints] = yRight[fMixFrac];
}

void TCountDiffFlameSim::UpdateSolutionOnePoint( Double *y, int gridPoint )
{
	T1DFlame::UpdateSolution( y, gridPoint );
	
	fSolU->vec[gridPoint] = y[fUVelocity];
	fSolV->vec[gridPoint] = y[fVVelocity];
	fSolMixFrac->vec[gridPoint] = y[fMixFrac];
}

void TCountDiffFlameSim::SolutionToSolver( void )
{
	TNewtonPtr	bt = fSolver->bt;
	TGridPtr	grid = bt->GetGrid()->GetFine();
	int		nGridPoints = fSolV->len;
	Double	*V = fSolV->vec;
	Double	*U = fSolU->vec;
	Double	*Z = fSolMixFrac->vec;
	Double	**y = grid->GetY()->mat;

	T1DFlame::SolutionToSolver();
	
	for ( int k = 0; k < nGridPoints; ++k ) {
		y[k][fVVelocity] = V[k];
		y[k][fUVelocity] = U[k];
		y[k][fMixFrac] = Z[k];
	}
	
	CountDiffSimPostIter( this );
}

void TCountDiffFlameSim::SaveSolution( void )
{
	int		k;
	int		len = fSolV->len;
	Double	*v = fSolV->vec;
	Double	*saveV = fSavedV->vec;
	Double	*u = fSolU->vec;
	Double	*saveU = fSavedU->vec;
	Double	*mixFrac = fSolMixFrac->vec;
	Double	*saveMixFrac = fSavedMixFrac->vec;

	T1DFlame::SaveSolution();
	fSavedV->len = fSolV->len;
	fSavedU->len = fSolU->len;
	fSavedMixFrac->len = fSolMixFrac->len;

	for ( k = -1; k <= len; ++k ) {
		saveV[k] = v[k];
		saveU[k] = u[k];
		saveMixFrac[k] = mixFrac[k];
	}
}

void TCountDiffFlameSim::RestoreSolution( void )
{
	int		k;
	int		len = fSavedV->len;
	Double	*v = fSolV->vec;
	Double	*saveV = fSavedV->vec;
	Double	*u = fSolU->vec;
	Double	*saveU = fSavedU->vec;
	Double	*mixFrac = fSolMixFrac->vec;
	Double	*saveMixFrac = fSavedMixFrac->vec;

	UpdateDimensions( len );

	T1DFlame::RestoreSolution();

	for ( k = -1; k <= len; ++k ) {
		v[k] = saveV[k];
		u[k] = saveU[k];
		mixFrac[k] = saveMixFrac[k];
	}
	
	SolutionToSolver();
}

int	TCountDiffFlameSim::GetOffsetVVelocity( void )
{
	return fVVelocity; 
}

int	TCountDiffFlameSim::GetOffsetUVelocity( void )
{
	return fUVelocity; 
}

int	TCountDiffFlameSim::GetOffsetTemperature( void )
{
	return fTemperature; 
}

int TCountDiffFlameSim::GetOffsetMixFrac( void )
{
	return fMixFrac;
}

int	TCountDiffFlameSim::GetOffsetFirstSpecies( void ) 
{
	return fFirstSpecies;
}

ConstStringArray TCountDiffFlameSim::GetVariableNames( void )
{
	return ( ConstStringArray ) fVariableNames;
}

int TCountDiffFlameSim::GetVariablesWithoutSpecies( void )
{
	return fVariablesWithoutSpecies;
}

void TCountDiffFlameSim::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp )
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
	Flag				oxidizerFound = FALSE;
	Flag				chooseInputGrid = FALSE;
	Double				*xIn = new Double[gridPointsIn];
	if ( !xIn ) FatalError( "memory allocation of TCountDiffFlameSim failed" );
	Double				*yIn =  new Double[gridPointsIn];
	if ( !yIn ) FatalError( "memory allocation of TCountDiffFlameSim failed" );
	Double				*yInFloat = sp->data;
	Double				**y = grid->GetY()->mat;
	int					variable;
	char				*string = sp->labels;
	SplinePtr			theSpline = NULL;
	Double				leftSlope;
	Double				rightSlope;
	int					oxidizerSide; // this program assumes that oxidizerSide = kRight
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
	
// find independent coordinate and oxidizerSide
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
		else if ( strcmp( string, "massfraction-o2" ) == 0 ) {
			if ( yInFloat[i * gridPointsIn] < yInFloat[(i+1) * gridPointsIn - 1] ) {
				oxidizerSide = kRight;
			}
			else {
				oxidizerSide = kLeft;
			}
			oxidizerFound = TRUE;
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
	if ( !oxidizerFound ) {
		cerr << "error: can't find massfraction of oxidizer 'massfraction-o2'" << NEWL;
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
		else if ( strcmp( string, "z" ) == 0 ) {
			variable = fMixFrac;
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
			if ( oxidizerSide == kRight ) {
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
				for ( k = 0; k < gridPointsIn-2; ++k ) { // turn vector and transform to current coordinate system
					if ( variable != fVVelocity || yInFloat[i*gridPointsIn] > yInFloat[i*gridPointsIn + 3] ) {
						y[k][variable] = yInFloat[(i+1)*gridPointsIn - k-2];	// copy workspace to vector of solution
					}
					else {
						y[k][variable] = -yInFloat[(i+1)*gridPointsIn - k-2];	// copy workspace to vector of solution
					}
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
			if ( oxidizerSide == kRight ) {
				for ( k = 0; k < nGridPoints; ++k ) {	// transform to current coordinate system
					if ( variable != fVVelocity || yInFloat[i*gridPointsIn] < yInFloat[i*gridPointsIn + 3] ) {
						y[k][variable] = yWork[k];	// copy workspace to vector of solution
					}
					else {
						y[k][variable] = -yWork[k];	// copy workspace to vector of solution and transform f
					}
				}
			}
			else {
				for ( k = 0; k < nGridPoints; ++k ) { // turn vector and transform to current coordinate system
					if ( variable != fVVelocity || yInFloat[i*gridPointsIn] > yInFloat[i*gridPointsIn + 3] ) {
						y[k][variable] = yWork[nGridPoints-k-1];	// copy workspace to vector of solution
					}
					else {
						y[k][variable] = -yWork[nGridPoints-k-1];	// copy workspace to vector of solution
					}
				}
			}
		}
		//	set bc
		if ( oxidizerSide == kRight ) {
			if ( variable == fTemperature ) {
				if ( yLeft[fTemperature] <= 0.0 ) {
					yLeft[fTemperature] = yInFloat[i*gridPointsIn];
				}
				if ( yRight[fTemperature] <= 0.0 ) {
					yRight[fTemperature] = yInFloat[(i+1)*gridPointsIn-1];
				}
			}
		}
		else {
			if ( variable == fTemperature ) {
				if ( yLeft[fTemperature] <= 0.0 ) {
					yLeft[fTemperature] = yInFloat[(i+1)*gridPointsIn-1];
				}
				if ( yRight[fTemperature] <= 0.0 ) {
					yRight[fTemperature] = yInFloat[i*gridPointsIn];
				}
			}
		}
	}

// set left to zero
	if (  oxidizerSide == kRight ) {
		for ( k = 0; k < nGridPoints; ++k ) {
			locX[k] -= bt->GetLeft();
		}
		bt->SetRight( bt->GetRight() - bt->GetLeft() );
		bt->SetLeft( 0.0 );
	}
	else {
		for ( k = 0; k < nGridPoints; ++k ) {
			yWork[k] = bt->GetRight() - locX[nGridPoints-k-1]; 	// turn sign and sort for ascending values
		}
		for ( k = 0; k < nGridPoints; ++k ) {
			locX[k] = yWork[k];					// copy workspace to x
		}
		bt->SetRight( bt->GetRight() - bt->GetLeft() );
		bt->SetLeft( 0.0 );
	}
	
// set initial Boundary values and pressure

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
//	update properties
	UpdateSolution( yMat, yLeftVec, yRightVec );
	if ( GetSpecies()->IsConstantLewisNumber() ) {
		CompLewisNumbers( GetSpecies()->GetLewisNumberFile() );
	}

	SetFlameNode( kPrev );
	ComputeProperties( fFlameNode, temp[kPrev], Y[kPrev], GetPressure() );
	SetFlameNode( nGridPoints );
	ComputeProperties( fFlameNode, temp[nGridPoints], Y[nGridPoints], GetPressure() );
	
	if ( inp->leftBoundary->fBcFlag[inp->fVVelocityOffset] == kNone ) {	// f' - g = 0
		int		mixtureSpecificationLeft = GetMixtureSpecificationLeft();
		Double	*rho = GetProperties()->GetDensity()->vec;
		// first set g
		if ( mixtureSpecificationLeft == kMassFlux ) {	// means tsuji bc
			yLeft[fUVelocity] = 0.0;
		}
		else {
			yLeft[fUVelocity] = sqrt( rho[nGridPoints] / rho[-1] );
		}
		yLeft[fVVelocity] = y[0][fVVelocity] - ( locX[0] - bt->GetLeft() ) 
							* yLeft[fUVelocity];
	}
	CountDiffSimPostIter( this );

	FreeSpline( theSpline );
	delete yIn;
	delete xIn;

	adapGrid->SetSolutionScaler();
	
	fp = GetOutfile( "initialguess", TFlame::kData );
	bt->PrintSolution( locX, y, GetVariableNames(), fp );
	fclose(fp);
}


void TCountDiffFlameSim::FillJacSpeciesDiffusion( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
{
// fills the jacobian with     constCoeff * d/dy ( rho * rho * diffusivity * df/dy )

	int		speciesIndex = nVariable - fFirstSpecies;
	Double	diff = fFlameNode->diffusivity[speciesIndex];
	Double	diffNext = fFlameNode->diffusivityNext[speciesIndex];
	Double	diffPrev = fFlameNode->diffusivityPrev[speciesIndex];
	Double	*rho = fFlameNode->mixDensity;

	if ( sign == kNegative ) {
		constCoeff *= -1.0;
	}
	Double	diffPlusHm = nodeInfo->hm * constCoeff * ( diff * rho[kCurr] * rho[kCurr]
										+ diffNext * rho[kNext] * rho[kNext] );
	Double	diffMinusH = nodeInfo->h * constCoeff * ( diffPrev * rho[kPrev] * rho[kPrev]
					+ diff * rho[kCurr] * rho[kCurr] );

	nodeInfo->a[nVariable][nVariable] -= ( diffPlusHm + diffMinusH );
	if ( !nodeInfo->lastPoint ) {
		nodeInfo->b[nVariable][nVariable] += diffPlusHm;
	}
	if ( !nodeInfo->firstPoint ) {
		nodeInfo->c[nVariable][nVariable] += diffMinusH;
	}
}

Double TCountDiffFlameSim::HeatFluxBinSpecDiff( NodeInfoPtr nodeInfo )
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

Double TCountDiffFlameSim::SecondDerivBinSpecDiff( int nVariable, NodeInfoPtr nodeInfo )
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

Double TCountDiffFlameSim::SecondDerivSpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo )
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

Double TCountDiffFlameSim::SecondDerivXDiffusion( int nVariable, NodeInfoPtr nodeInfo )
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

void TCountDiffFlameSim::PrintRHSTemp( TNewtonPtr bt )
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

void TCountDiffFlameSim::PrintRHSTemp( TNewtonPtr bt, NodeInfoPtr nodeInfo, Double physX, FILE *fp )
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

void TCountDiffFlameSim::PrintRHSSpecies( TNewtonPtr bt )
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

void TCountDiffFlameSim::PrintRHSSpecies( int start, int end, NodeInfoPtr nodeInfo, Double physX, FILE *fp )
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

void CountDiffSimOutput( void *object, FILE *fp, char* tail )
{
	TCountDiffFlameSimPtr	flame = ( TCountDiffFlameSimPtr )object;
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
	Double			*Z = flame->GetZ()->vec;
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
	int				fMixFrac = flame->GetOffsetMixFrac();
	time_t			theDate;
	char			buffer[80];
	ConstStringArray	varNames = flame->GetVariableNames();
	char			**names = species->GetNames();
	VectorPtr 		physXVec = NewVector( gridPoints + 2 );
	Double			*physX = physXVec->vec;
	VectorPtr 		scalarDissVec = NewVector( gridPoints + 2 );
	VectorPtr 		scalarDissVecB = NewVector( gridPoints + 2 );
	VectorPtr 		ZBilgerVec = NewVector( gridPoints + 2 );
	VectorPtr 		LewisFuelVec = NewVector( gridPoints + 2 );
	Double			*scalarDiss = scalarDissVec->vec;
	Double			*scalarDissB = scalarDissVecB->vec;
	Double			*ZBilger = &ZBilgerVec->vec[kNext];
	Double			*LewisFuel = &LewisFuelVec->vec[kNext];
	Double			stoechScalarDiss;
	Double			stoechScalarDissB;
	Flag			fpOpen = FALSE;
	Flag			isPremixed = (flame->GetInputData()->fPremConfiguration == 1) ? TRUE : FALSE;
	
	if ( !fp ) {
		fpOpen = TRUE;
		fp = flame->GetOutputFile( NULL, tail, TFlame::kNone );
	}

	if ( !isPremixed ) {
		ZBilger[kPrev] = Z[kPrev];
		ZBilger[gridPoints] = Z[gridPoints];
		for ( k = 0; k < gridPoints; ++k ) {
			ZBilger[k] = flame->ComputeZBilger( massFracs[k], massFracs[-1], massFracs[gridPoints] );
		}
		flame->ScalarDissipation( ZBilger, bt, &scalarDissVecB->vec[kNext], &stoechScalarDissB );
		flame->ScalarDissipation( Z, bt, &scalarDissVec->vec[kNext], &stoechScalarDiss );
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
	if ( !isPremixed ) {
		fprintf( fp, "chiStoichAnal = %g [1/s]\n", flame->GetChiStoichAnal() );
		Double	r1 = sqrt( rho[gridPoints] / InterpolOne( flame->GetZStoich(), &ZBilger[kPrev]
				, &rho[kPrev], gridPoints+2 ) );
		Double	r2 = 0.75 * ( r1 + 1.0 ) * ( r1 + 1.0 ) / ( 2.0 * r1 + 1.0 );
		fprintf( fp, "chiStoichAnalKW = %g [1/s]\n", flame->GetChiStoichAnal() * r2 );
	}
	if ( species->IsConstantLewisNumber() ) {
		fprintf( fp, "ConstantLewisNumbers = \"True\"\n" );
	}
	
	if ( !isPremixed ) {
		fprintf( fp, "FlameLocZ = %g [K]\n", Z[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
	}
	fprintf( fp, "FlameLocBilger = %g [K]\n", ZBilger[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
	fprintf( fp, "Tmax = %g [K]\n", temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
	if ( !isPremixed ) {
		fprintf( fp, "stoichScalarDissRate = %g [1/s]\n", stoechScalarDiss );
		fprintf( fp, "stoichScalarDissRateBilger = %g [1/s]\n", stoechScalarDissB );
	}
	
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
	
//  write bc
	Double	locMoleMass;
	fprintf( fp, "\nOxidizerSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", temp[gridPoints] );
	for ( i = 0; i < nOfSpecies; ++i ) { // write Y_i
		if ( fabs( massFracs[gridPoints][i] ) > 1.0e-5 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], massFracs[gridPoints][i] );
		}
	}
	for ( i = 0; i < nOfSpecies; ++i ) { // write X_i
		locMoleMass = massFracs[gridPoints][i] * mixMolarMass[gridPoints] / molarMass[i];
		if ( fabs( locMoleMass ) > 1.0e-5 ) {
			fprintf( fp, "\tMolefraction-%s = %g\n", names[i], locMoleMass );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "FuelSide\n" );
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
	flame->PrintFlameletVector( gridPoints+2, &Z[kPrev], "Z", fp );
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


	if ( !isPremixed ) {
//	write scalar dissipation rate
	fprintf( fp, "ScalarDissRate [1/s]\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", scalarDiss[k] );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( (k+1) % 5 ) {
		fprintf( fp, "\n" );
	}

//	write mixture fraction and scalar dissipation rate
	flame->PrintFlameletVector( gridPoints+2, &ZBilger[kPrev], "ZBilger", fp );
	flame->PrintFlameletVector( gridPoints+2, scalarDissB, "ScalarDissRateBilger [1/s]", fp );
	
//	write chi analyt
	fprintf( fp, "ChiAnalytOverChiStoich\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", flame->GetChiOverChiStoich( ZBilger[k-1] ) );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( (k+1) % 5 ) {
		fprintf( fp, "\n" );
	}
	
//	write chi K&W
	Double rhoStoich = InterpolOne( flame->GetZStoich(), &ZBilger[kPrev]
			, &rho[kPrev], gridPoints+2 );
	Double	rat, fact;
	Double	ratStoich = sqrt( rho[gridPoints] / rhoStoich );
	Double	factStoich = 1.5 * ( ratStoich + 1.0 ) * ( ratStoich + 1.0 ) / ( 2.0 * ratStoich + 1.0 );
	fprintf( fp, "ChiKWOverChiKWStoich\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		rat = sqrt( rho[gridPoints] / rho[k-1] );
		fact = 1.5 * ( rat + 1.0 ) * ( rat + 1.0 ) / ( 2.0 * rat + 1.0 );
		fprintf( fp, "\t%-.6e", fact / factStoich * flame->GetChiOverChiStoich( ZBilger[k-1] ) );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( (k+1) % 5 ) {
		fprintf( fp, "\n" );
	}
	
//	write chi K&W factor
	fprintf( fp, "ChiKWFact\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		rat = sqrt( rho[gridPoints] / rho[k-1] );
		fact = 1.5 * ( rat + 1.0 ) * ( rat + 1.0 ) / ( 2.0 * rat + 1.0 );
		fprintf( fp, "\t%-.6e", fact / factStoich );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( (k+1) % 5 ) {
		fprintf( fp, "\n" );
	}

	}
	
//	write V_DY fuel
	
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
			fprintf( fp, "@ Z = %g\n", Z[maxLoc] );
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
			fprintf( fp, "@ Z = %g\n", Z[maxLoc] );
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
		fprintf( fp, "@ Z = %g\n", Z[maxLoc] );
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
	
	DisposeVector( ZBilgerVec );
	DisposeVector( scalarDissVecB );
	DisposeVector( scalarDissVec );
	DisposeVector( physXVec );
	
	if ( fpOpen ) {
		fclose( fp );
	}
}

void TCountDiffFlameSim::CompLewis( int which, Double *Le )
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
		fprintf( stderr, "#error: no Species N2 in function 'TCountDiffFlameSim::CompLewis'\n" ); 
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

Double TCountDiffFlameSim::ComputeAbsStrainRate( int k1, Double f1, Double f2, Double fp1, Double fp2, Double h )
{
	Double	*rho = GetProperties()->GetDensity()->vec;
	Double	rho1 = rho[k1];
	Double	rho2 = rho[k1+1];
	Double	strainRate = GetStrainRate();
	Double	fpm = 0.5 * ( fp1 + fp2 );
	
	Double fact = 0.5 * ( rho1 + rho2 ) 
					* ( f2/rho2 - f1/rho1 ) / h;
	return strainRate * sqrt( fpm * fpm + fact * fact ); 
}

void TCountDiffFlameSim::FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
{
// fills the jacobian with     constCoeff * sum_j ( d/dy(rho^2 Y_k D_j dY_j/dy) )

	if ( sign == kNegative ) {
		constCoeff *= -1.0;
	}

	int		i, lVar;
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	**a = nodeInfo->a;
	Double	**b = nodeInfo->b;
	Double	**c = nodeInfo->c;
	Double	*Y = fFlameNode->Y[kCurr];
	Double	*YPrev = fFlameNode->Y[kPrev];
	Double	*YNext = fFlameNode->Y[kNext];
	Double	*density = fFlameNode->mixDensity;
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	lVar = nVariable-fFirstSpecies;
	Double	coeffCurr = constCoeff * density[kCurr] * density[kCurr] * Y[lVar];
	Double	coeffPrev = constCoeff * density[kPrev] * density[kPrev] * YPrev[lVar];
	Double	coeffNext = constCoeff * density[kNext] * density[kNext] * YNext[lVar];
	Double	diffPlusHm, diffMinusH;
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
		lVar = fFirstSpecies + i;
		// d/dY_l
		diffPlusHm = hm * ( coeffCurr * diffusivity[i] + coeffNext * diffusivityNext[i] );
		diffMinusH = h * ( coeffCurr * diffusivity[i] + coeffPrev * diffusivityPrev[i] );
		a[lVar][nVariable] -= diffPlusHm + diffMinusH;
		if ( !nodeInfo->lastPoint ) {
			b[lVar][nVariable] += diffPlusHm;
		}
		if ( !nodeInfo->firstPoint ) {
			c[lVar][nVariable] += diffMinusH;
		}
	}
}

Double TCountDiffFlameSim::NewDiffCorr( int nVariable, NodeInfoPtr nodeInfo )
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

Double TCountDiffFlameSim::NewDiffCorrX( int nVariable, NodeInfoPtr nodeInfo )
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

void TCountDiffFlameSim::ScalarDissipation( Double *Z, TNewtonPtr bt, Double *scalarDiss, Double *stoechScalarDiss )
{
	int				k;
	int				nOfGridPoints = bt->GetCurrentGridPoints();
	TGridPtr		grid = bt->GetGrid()->GetCurrentGrid();
	NodeInfoPtr		nodeInfo = bt->GetNodeInfo();
	Double			**yMat = grid->GetY()->mat;
	Double			left = bt->GetLeft();
	Double			right = bt->GetRight();
	Double			*x = grid->GetX()->vec;
	Double			*density = fProperties->GetDensity()->vec;
	Double			*lambda = fProperties->GetConductivity()->vec;
	Double			*cp = fProperties->GetHeatCapacity()->vec;
	Double			rhoInf = fFlameNode->rhoInf;
	Double			viscosityInf = fFlameNode->viscosityInf;
	Double			factor = 2.0 * GetStrainRate() / ( rhoInf * viscosityInf );
	Double			dZdx;
	Double			zStoe = fMassFraction->GetZStoe();
	
	dZdx = ( Z[0] - Z[kPrev] ) / ( x[0] - left );
	scalarDiss[kPrev] = factor * density[kPrev] * lambda[kPrev] / cp[kPrev]  * dZdx * dZdx;
	for ( k = 0; k < nOfGridPoints; ++k ) {
		bt->SetNodeInfo( this, k );
		dZdx = FirstDeriv( Z[k-1], Z[k], Z[k+1], nodeInfo->hm, nodeInfo->h );
		scalarDiss[k] = factor * density[k] * lambda[k] / cp[k]  * dZdx * dZdx;
		if ( ( Z[k] - zStoe ) * ( Z[k-1] - zStoe ) <= 0.0 ) {
			*stoechScalarDiss = scalarDiss[k-1] + ( scalarDiss[k] - scalarDiss[k-1] ) 
							/ ( Z[k] - Z[k-1] ) * ( zStoe - Z[k-1] );
		}
	}
	dZdx = ( Z[nOfGridPoints] - Z[nOfGridPoints-1] ) / ( right - x[nOfGridPoints-1] );
	scalarDiss[nOfGridPoints] = factor * density[nOfGridPoints] * lambda[nOfGridPoints] / cp[nOfGridPoints-1]  * dZdx * dZdx;
}

void SetCountDiffSimNodeInfo( int k, void *object )
{
	TCountDiffFlameSimPtr	flame = ( TCountDiffFlameSimPtr )object;
	
	flame->SetFlameNode( k );
}

void CountDiffSimPostConv( void *object )
{
	TCountDiffFlameSimPtr	flame = ( TCountDiffFlameSimPtr )object;
	TNewtonPtr 				bt = flame->GetSolver()->bt;
	Double					*temp = flame->GetTemperature()->vec;
	Double					*Z = flame->GetZ()->vec;
	int						gridPoints = bt->GetGrid()->GetCurrentGrid()->GetNGridPoints();
	int						isConverged = bt->GetConvergeNewton();
	Double					coeffSpecies = 1.0 / flame->GetStrainRate();
	Double					coeffTemp = 1.0 / flame->GetStrainRate();
	int						nSpeciesIn = flame->GetSpecies()->GetNSpeciesInSystem();
	Double					*Le = flame->GetSpecies()->GetLewisNumber()->vec;
	fprintf( stderr, "Tmax = %g @ Z = %g\n"
			, temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1]
			, Z[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
	if ( isConverged ) {
		flame->SaveSolution();
		if ( flame->GetInputData()->fAdjustComputationalDomain ) {
			if( flame->AdjustGrid( CountDiffSimPostIter ) ) {
				return;
			}
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
	}
	else {
		flame->RestoreSolution();
		CountDiffSimPostIter( flame );
	}

	flame->PostConvergence( object );
	CountDiffSimPostIter( flame );
}

ConstStringArray GetCountDiffSimVarNames( void *object )
{
	TCountDiffFlameSimPtr	flame = ( TCountDiffFlameSimPtr )object;
	
	return flame->GetVariableNames();
}

FILE *TCountDiffFlameSim::GetOutputFile( char *head, char *tail, FileType type )
{
	int				fuelIndex = GetFuelIndex();
	char			*name = new char[64];
	FILE			*fp;
	char			**speciesNames = fSpecies->GetNames();
	int				tOxidizer = ( int ) fSolTemp->vec[fSolTemp->len];
	int				tFuel = ( int ) fSolTemp->vec[kPrev];
	Double			press = GetPressure() * 1.0e-5;
		
	sprintf( name, "%s%s%.8s_p%.2d_%.1da%.5d_%.1dtf%.4dto%.4d%s"
					, ( head ) ? head : "", ( head ) ? "_" : ""
					, speciesNames[fuelIndex]
					, ( int ) floor( press )	// in [bar]
					, ( int ) ( ( press - ( floor( press ) ) ) * 10 + 0.5 )
					, ( int ) floor( GetStrainRate() ) 
					, ( int ) ( ( GetStrainRate() - floor(GetStrainRate()) ) * 10 + 0.5 )
					, ( int )( tFuel )							// in [K]
					, ( int )( tOxidizer ) 						// in [K]
					, ( tail ) ? tail : "" );
	fp = GetOutfile( name, type );
	delete name;

	return fp;
}

void TCountDiffFlameSim::EtaToX( TNewtonPtr bt, VectorPtr xPhysVec )
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
	
	// linear Interpolation
	
	for ( k = 0; k < gridPoints-1; ++k ) {
		if ( ( y[k][fMixFrac] - zStoech ) * ( y[k+1][fMixFrac] - zStoech ) <= 0.0 ) {
			kBefStoech = k;
			break;
		}
	}
	
	if ( kBefStoech < 0 ) {
		cerr << "##warning: can't find the point of stoichiometric mixture fraction" << NEWL;
		return;
	}
	
	xStoech = xPhys[kBefStoech+1] + ( xPhys[kBefStoech+2] - xPhys[kBefStoech+1] ) 
						/ ( y[kBefStoech+1][fMixFrac] - y[kBefStoech][fMixFrac] )
						* ( zStoech - y[kBefStoech][fMixFrac] );
	
	for ( k = 0; k < gridPoints+2; ++k ) {
		xPhys[k] -= xStoech;
	}
}


Double EntFluxFuncDiffSim( int /*equation*/, NodeInfoPtr nodeInfo, void *object, Flag /*theFlag*/  )
{
	Double	entFlux = 0.0;
	TCountDiffFlameSimPtr	flame = ( TCountDiffFlameSimPtr )object;
	TFlameNodePtr			flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int		eqLoop, speciesEq;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*mixDensity = flameNode->mixDensity;
	Double	mixDensityInf = flameNode->rhoInf;
	Double	mixViscosityInf = flameNode->viscosityInf;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	constThermDiffCoeff = 1.0 / ( *flameNode->mixHeatCapacity * mixDensityInf * mixViscosityInf );
	Double	*diffCorr = flameNode->diffCorr;
	Double	h = nodeInfo->h;
	Double	hm = nodeInfo->hm;
	Double	sumCpY = 0.0;
	Double	entCoeff = 0.0;
	Double	sumCpDdYdx = 0.0;
			
	entCoeff = constThermDiffCoeff * mixDensity[kCurr] * mixDensity[kCurr]
					* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );

	for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
		speciesEq = eqLoop + fFirstSpecies;
		sumCpDdYdx += heatCapacity[eqLoop] * diffusivity[eqLoop] 
					* FirstDeriv( yPrev[speciesEq], y[speciesEq], y[speciesEq], hm, h );
		if ( flame->UseDiffCorr() ) {
			sumCpY += heatCapacity[eqLoop] * y[speciesEq];
		}
	}
	entFlux = sumCpDdYdx;
	if ( flame->UseDiffCorr() ) {
		entFlux -= diffCorr[kCurr] * sumCpY;
	}
	entFlux *= entCoeff;
	
	return entFlux;
}

void CountDiffSimUpdateLeftBoundary( void  *object )
{
	TCountDiffFlameSimPtr	flame = ( TCountDiffFlameSimPtr )object;
	TNewtonPtr		bt = flame->GetSolver()->bt;
	TGridPtr 		currGrid = bt->GetGrid()->GetCurrentGrid();
	MatrixPtr		yMat = currGrid->GetY();
	Double			*yLeft = currGrid->GetYLeft()->vec;
	Double			pressure = flame->GetPressure();

}

void CountDiffSimUpdateRightBoundary( void *object )
{
	TCountDiffFlameSimPtr	flame = ( TCountDiffFlameSimPtr )object;
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

	yRight[fUVelocity] = 1.0;
	yRight[fVVelocity] = yLast[fVVelocity] + hLast * yRight[fUVelocity];

	flame->UpdateSolutionOnePoint( yRight, nGridPoints );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], pressure );
}

Double TCountDiffFlameSim::GetChiOverChiStoich( Double z )
{
	Double	zStoi = GetZStoich();
	Double	zStoiNew = zStoi - 0.5;
	Double	zStoiNew2 = zStoiNew * zStoiNew;
	Double	zStoiNew4 = zStoiNew2 * zStoiNew2;
	Double	zStoiNew6 = zStoiNew4 * zStoiNew2;
	Double	zStoiNew8 = zStoiNew6 * zStoiNew2;
	Double	zNew = z - 0.5;
	Double	zNew2 = zNew * zNew;
	Double	zNew4 = zNew2 * zNew2;
	Double	zNew6 = zNew4 * zNew2;
	Double	zNew8 = zNew6 * zNew2;
	
	Double	duedxst = 1.0 / ( 1.00539 * ( 12.9041 - 
					82.123 * zStoiNew2 +
                    115.29 * zStoiNew4 -
                    201.898 * zStoiNew6 +
                    912.136 * zStoiNew8 ) );
	Double	chiNow = 1.00539       *  duedxst *
                    (12.9041       - 
                    82.123  * zNew2 +
                    115.29  * zNew4 -
                    201.898 * zNew6 +
                    912.136 * zNew8 );

    if ( chiNow < 1.0e-10 ) {
		chiNow = 1.0e-10;
	}

	return chiNow;
}

Double TCountDiffFlameSim::GetChiStoichAnal( void )
{
	const Double pi = 3.1415927;
	Double	zStoi = GetZStoich();
	Double	zStoiNew = zStoi - 0.5;
	Double	zStoiNew2 = zStoiNew * zStoiNew;
	Double	zStoiNew4 = zStoiNew2 * zStoiNew2;
	Double	zStoiNew6 = zStoiNew4 * zStoiNew2;
	Double	zStoiNew8 = zStoiNew6 * zStoiNew2;
	
	Double	fOfZ = -6.4096 * zStoiNew2 
					+ 9.8518 * zStoiNew4 
					-21.72 * zStoiNew6 
					+ 83.162 * zStoiNew8 
					+ 1.0007;

	return GetStrainRate() * fOfZ / pi;
}

Double TCountDiffFlameSim::SecondDerivConstMixFracDiff( int nVariable, NodeInfoPtr nodeInfo )
{
	Double	lambdainf = fProperties->GetConductivity()->vec[fSolver->bt->GetCurrentGridPoints()];
	Double	cpinf = fProperties->GetHeatCapacity()->vec[fSolver->bt->GetCurrentGridPoints()];
	Double	lambdaOverCpInf = lambdainf / cpinf;
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	yPrev = nodeInfo->yPrev[nVariable];
	Double	y = nodeInfo->y[nVariable];
	Double	yNext = nodeInfo->yNext[nVariable];
	Double	*rho = fFlameNode->mixDensity;
	Double	diffPlus =  rho[kCurr]
					+ rho[kNext];
	Double	diffMinus = rho[kPrev]
					+ rho[kCurr];
	
	return lambdaOverCpInf * ( diffPlus * hm * yNext - ( diffPlus * hm + diffMinus * h ) * y + diffMinus * h * yPrev ) 
				/ ( h * hm * ( h + hm ) );
}
