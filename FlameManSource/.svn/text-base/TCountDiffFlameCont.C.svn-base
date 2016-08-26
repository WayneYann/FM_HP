//#define FOOL_SOFTBENCH( x ) 

#include "FlameMaster.h"
#include "ListTool.h"
#include "Spline.h"
#include "TCountDiffFlameCont.h"

#undef UPWINDCONVECTION

#define MOLARDIFFUSION
#define DIFFUSIVITYCORRECTION

#define ENTHALPYFLUX

#undef ENTFLUXINJAC

#undef DEBUGNEWGUESS

void CountDiffContPostConvPre( void *object )
{
	TCountDiffFlameContPrePtr	flame = ( TCountDiffFlameContPrePtr )object;
	TNewtonPtr	bt = flame->GetSolver()->bt;
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = currGrid->GetNGridPoints();
	int			count = flame->GetStrainRateVector()->len;
	int			tmaxLoc;
	int			isConverged = bt->GetConvergeNewton();
	Double		*temp = flame->GetTemperature()->vec;
	static Double	firstStrainRate = flame->fStrainRate->vec[0];

	tmaxLoc = LocationOfMax( nGridPoints+2, &temp[kPrev] ) - 1;
	flame->fTmax[count] = temp[tmaxLoc];
	flame->fAInv[count] = 1.0 / flame->GetStrainRate();

	CountDiffSimPostConv( object );
	
	if ( bt->GetLeaveContin() ) {
		if ( bt->GetConvergeNewton()
			&& flame->fStrainRate->vec[0] > firstStrainRate ) {
			flame->PostAll();
		}
		else {
			if ( isConverged ) {
				flame->fStrainRate->vec[0] *= 1.01;
				bt->UnSetLeaveContin();
				bt->InitNIter();
		//		bt->GetGrid()->UnSetSolHasNewGrid();
				cerr << "strainRate is now " << flame->GetStrainRate() << NEWL;
			}
			else {
				exit(2);
			}
		}
	}
}

void TCountDiffFlameContPre::PostAll( void )
{
	TNewtonPtr	bt = GetSolver()->bt;
	TGridPtr	grid = bt->GetGrid()->GetFine();
	MatrixPtr	yFine = grid->GetY();
	VectorPtr	xFine = grid->GetX();
	VectorPtr	yLeft = grid->GetYLeft();
	VectorPtr	yRight = grid->GetYRight();
	int			nSpeciesIn = fSpecies->GetNSpeciesInSystem();
	Double		dT = fTmax[1] - fTmax[0];
	Double		dAInv = fAInv[1] - fAInv[0];
	Double		ds = sqrt( dT * dT + dAInv * dAInv );

	fdTds = dT / ds;
	fdAInvds = dAInv / ds;
	
	fSolPack->fSol = fSolPack->CopyIt( yFine );
	fSolPack->fCoord = fSolPack->CopyIt( xFine );
	fSolPack->fSolLeft = fSolPack->CopyIt( yLeft );
	fSolPack->fSolRight = fSolPack->CopyIt( yRight );
	fSolPack->fSolNames = fSolPack->CopyIt( fVariableNames, fVariablesWithoutSpecies + nSpeciesIn );
	fSolPack->fLeft = bt->GetGrid()->GetLeft();
	fSolPack->fRight = bt->GetGrid()->GetRight();
	fSolPack->fPressure = GetPressure();
}

void TCountDiffFlameContPre::InitCountDiffFlameContPre( FirstInputPtr /*firstInp*/ )
{
	TNewtonPtr		bt = GetSolver()->bt;
	TGridPtr		fine = bt->GetGrid()->GetFine();

/*	if ( firstInp->nStrainRates != 2 ) {
		cerr << "#error: two different strainrates have to be specified" << NEWL;
		exit( 2 );
	}*/

	bt->SetUtFuncs( CountDiffSimJacRest, CountDiffSimJacRest, CountDiffSimJacRest
					, CountDiffSimRHSRest, CountDiffSimRHSRest, CountDiffSimRHSRest
					, CountDiffSimOutput, CountDiffSimPostIter
					, SetCountDiffSimNodeInfo, CountDiffContPostConvPre
					, GetCountDiffSimVarNames
					, CountDiffSimUpdateLeftBoundary, CountDiffSimUpdateRightBoundary);
	fSolPack = new TSolution;
}

TCountDiffFlameContPre::~TCountDiffFlameContPre( void )
{
//	delete fSolPack;
}

void TCountDiffFlameCont::InitCountDiffFlameCont( TSolutionPtr sol, Double strainRate )
{
	TBVPSolverPtr	solver = GetSolver();
	TNewtonPtr		bt = solver->bt;
	TGridPtr		fine = bt->GetGrid()->GetFine();
	TGridPtr		coarse = bt->GetGrid()->GetCoarse();
	int				nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	int				maxGridPoints = bt->GetMaxGridPoints();

	fDeltaS = fInputData->fDeltaSCont;
	fSConverged = 0.0;		//	value of input file
	fSNotConverged = -1.0;	//	means not set
	fCompUnPhysChain = fInputData->fCompUnPhysChain;
	fNUnPhysChain = fInputData->fNUnPhysChain;

	if ( fInputData->fMinStrainRate < 0.0 ) {
		fMinStrainRate = 10.0;
	}
	else {
		fMinStrainRate = fInputData->fMinStrainRate;
	}

	fVariableNames = new String[fVariablesWithoutSpecies + nSpeciesInSystem];
   	if ( !fVariableNames ) FatalError( "memory allocation of TCountDiffFlameCont failed" );

	fVariableNames[fVVelocity] = new char[2];
	strcpy( fVariableNames[fVVelocity], "f" );
	fVariableNames[fUVelocity] = new char[3];
	strcpy( fVariableNames[fUVelocity], "f'" );
	fVariableNames[fMixFrac] = new char[2];
	strcpy( fVariableNames[fMixFrac], "Z" );
	fVariableNames[fTemperature] = new char[2];
	strcpy( fVariableNames[fTemperature], "T" );
	fVariableNames[fInvStrainRate] = new char[4];
	strcpy( fVariableNames[fInvStrainRate], "1/a" );
	for ( int i = 0; i < nSpeciesInSystem; ++i ) {
		fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies->GetNames()[i] ) + 1];
		strcpy( fVariableNames[fFirstSpecies + i], fSpecies->GetNames()[i] );
	}
	
//	vectors of solution
	fSolV = NewVector( maxGridPoints + 2 );
	fSolU = NewVector( maxGridPoints + 2 );
	fSolMixFrac = NewVector( maxGridPoints + 2 );
	fSolInvStrainRate = NewVector( maxGridPoints + 2 );

	fSolV->vec = &fSolV->vec[kNext];
	fSolU->vec = &fSolU->vec[kNext];
	fSolMixFrac->vec = &fSolMixFrac->vec[kNext];
	fSolInvStrainRate->vec = &fSolInvStrainRate->vec[kNext];

	fSolV->len -= 2;
	fSolU->len -= 2;
	fSolMixFrac->len -= 2;
	fSolInvStrainRate->len -= 2;

//	saved solution
	fSavedV = NewVector( maxGridPoints + 2 );
	fSavedU = NewVector( maxGridPoints + 2 );
	fSavedMixFrac = NewVector( maxGridPoints + 2 );
	fSavedInvStrainRate = NewVector( maxGridPoints + 2 );

	fSavedV->vec = &fSavedV->vec[kNext];
	fSavedU->vec = &fSavedU->vec[kNext];
	fSavedMixFrac->vec = &fSavedMixFrac->vec[kNext];
	fSavedInvStrainRate->vec = &fSavedInvStrainRate->vec[kNext];

	fSavedV->len -= 2;
	fSavedU->len -= 2;
	fSavedMixFrac->len -= 2;
	fSavedInvStrainRate->len -= 2;

	bt->SetUtFuncs( CountDiffContJacRest, CountDiffContJacRest, CountDiffContJacRest
					, CountDiffContRHSRest, CountDiffContRHSRest, CountDiffContRHSRest
					, CountDiffContOutput, CountDiffContPostIter
					, SetCountDiffContNodeInfo, CountDiffContPostConv
					, GetCountDiffContVarNames
					, CountDiffContUpdateLeftBoundary, CountDiffContUpdateRightBoundary);
	SetInitialBC( fine, fInputData );
	SetInitialBC( coarse, fInputData );
	fMassFraction = new TMassFraction( this );
	if ( !fMassFraction ) FatalError( "memory allocation of TMassFraction failed" );
	if ( sol ) {
		SetInitialValues( fInputData, sol, strainRate );
	}
	else {
		ReadStartProfiles( fInputData );
	}
	CheckBC();
	CheckInitialGuess();
	UpdateSolution( fine->GetY(), fine->GetYLeft(), fine->GetYRight() );	
	fTmaxStart = fSolTemp->vec[fTmaxLoc];
	sprintf( GetOutFileBuff(), "%sExt%.8s_p%.2dTox%.4d.dout"
			, GetOutputPath()
			, fSpecies->GetNames()[GetFuelIndex()]
			, ( int )( GetPressure() * 1.0e-5 )
			, ( int )( fine->GetYRight()->vec[fTemperature] ) );
	if ( !( fpExtCurve = fopen( GetOutFileBuff(), "w") ) ) { 
		cerr << "#warning: unable to open file " << GetOutFileBuff() << NEWL;
		exit(2);
	}
	fprintf( fpExtCurve, "*\n%-12s\t%-12s\t%-12s\n", "a [s\\u-1\\n]", "1/a [s]", "T\\dmax\\n [K]" );
}

TCountDiffFlameCont::~TCountDiffFlameCont( void )
{
	int	nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();

	fclose ( fpExtCurve );

	delete fMassFraction;

	fSavedInvStrainRate->vec = &fSavedInvStrainRate->vec[kPrev];
	fSavedMixFrac->vec = &fSavedMixFrac->vec[kPrev];
	fSavedU->vec = &fSavedU->vec[kPrev];
	fSavedV->vec = &fSavedV->vec[kPrev];

	DisposeVector( fSavedInvStrainRate );
	DisposeVector( fSavedMixFrac );
	DisposeVector( fSavedU );
	DisposeVector( fSavedV );

	fSolInvStrainRate->vec = &fSolInvStrainRate->vec[kPrev];
	fSolMixFrac->vec = &fSolMixFrac->vec[kPrev];
	fSolV->vec = &fSolV->vec[kPrev];
	fSolU->vec = &fSolU->vec[kPrev];

	DisposeVector( fSolInvStrainRate );
	DisposeVector( fSolMixFrac );
	DisposeVector( fSolU );
	DisposeVector( fSolV );

	for ( int i = 0; i < nSpeciesInSystem+fVariablesWithoutSpecies; ++i ) {
		delete fVariableNames[i];
	}
	delete fVariableNames;
}

void CountDiffContJacRest( void *object, NodeInfoPtr nodeInfo )
{
	TCountDiffFlameContPtr	flame = ( TCountDiffFlameContPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fInvStrainRate = flame->GetOffsetInvStrainRate();
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
	Double	mixDensity = *flameNode->mixDensity;
	Double	mixHeatCapacity = *flameNode->mixHeatCapacity;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	*enthalpy = flameNode->enthalpy;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	*productionRate = flameNode->productionRate;
	Double	mixDensityInf = flameNode->rhoInf;
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * flameNode->viscosityInf );
	Double	constThermDiffCoeff = constMassDiffCoeff / mixHeatCapacity;
	Double	oneOverATimesRho = y[fInvStrainRate] / mixDensity;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	energySourceCoeff = oneOverATimesRho / mixHeatCapacity;
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

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
	
//	Double	theConst = mixDensityInf * mixMolarMass / mixDensity * hnenn;
//	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq; ++speciesEq ) {
//		a[speciesEq][fUVelocity] += theConst / molarMass[speciesEq-fFirstSpecies];
//	}
	
// third equation ( mixturefraction )
	flame->FillJacMixFracDiffusion( fMixFrac, fMixFrac, constMassDiffCoeff, nodeInfo );

// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, constMassDiffCoeff, nodeInfo );
#ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			flame->FillJacDiffCorr( speciesEq, constMassDiffCoeff, nodeInfo, kNegative );
		}
#endif
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
		a[fInvStrainRate][speciesEq] += productionRate[speciesIndexEq] / mixDensity * hnenn;
	}

// inverse strainrate equation
	if ( nodeInfo->gridPoint < flame->fTmaxLoc ) {
		FillJacFirstDerivDown( fInvStrainRate, fInvStrainRate, nodeInfo );
	}
	else if ( nodeInfo->gridPoint == flame->fTmaxLoc ) {
//		cerr << yPrev[fTemperature] << TAB << y[fTemperature] << TAB << yNext[fTemperature] << NEWL;
		a[fTemperature][fInvStrainRate] += flame->fdTds * hnenn;
		a[fInvStrainRate][fInvStrainRate] += flame->fdAInvds * hnenn;
	}
	else {
		FillJacFirstDerivUp( fInvStrainRate, fInvStrainRate, nodeInfo );
	}

// last equation ( temperature )
	if ( fTemperature < M ) {
		int 	i, l;
		Double	sumMH = 0.0;
		flame->FillJacDiffusion( fTemperature, fTemperature, constThermDiffCoeff, flameNode->mixConductivity, nodeInfo );

		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			sumMH += productionRate[i] * enthalpy[i];
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
		a[fInvStrainRate][fTemperature] -= sumMH / ( mixDensity * mixHeatCapacity ) * hnenn;
		if ( flame->fProperties->GetRadiation() ) {
//			a[fInvStrainRate][fTemperature] += flameNode->radiation[kCurr] / ( mixDensity * mixHeatCapacity ) * hnenn;
//			a[fTemperature][fTemperature] += flameNode->radiation[kCurr] * oneOverATimesRho / ( mixHeatCapacity * temp[kCurr] ) * hnenn;
			flame->fProperties->GetRadiation()->FillJacRadiation( energySourceCoeff, flame, nodeInfo );
		}
#ifdef ENTFLUXINJAC
#ifdef ENTHALPYFLUX
		Double	*diffusivity = flameNode->diffusivity;
		Double	sumCpY = 0.0;
		Double	sumCpDdYdx = 0.0;
		Double	entFluxCoeff = constThermDiffCoeff * mixDensity * mixDensity
					* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
		Double	entTerm;
		Double	**Y = flameNode->Y;
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			sumCpDdYdx += heatCapacity[i] * diffusivity[i] 
						* FirstDeriv( Y[kPrev][i], Y[kCurr][i], Y[kNext][i], hm, h );
#	ifdef DIFFUSIVITYCORRECTION
			if ( flame->UseDiffCorr() ) {
				sumCpY += heatCapacity[i] * Y[kCurr][i];
			}
#	endif
		}
		entTerm = sumCpDdYdx;
#	ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			entTerm -= flameNode->diffCorr[kCurr] * sumCpY;
		}
#	endif
		entTerm *= entFluxCoeff;
		entTerm *= 2.0 * hnenn;
		for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq; ++speciesEq ) {
			i = speciesEq - fFirstSpecies;
			a[speciesEq][fTemperature] += dfdyUpwind( speciesEq, fTemperature, EntFluxFuncDiffSim, nodeInfo, flame ) * hnenn;
			a[speciesEq][fTemperature] -= entTerm 
					* mixMolarMass / molarMass[speciesEq-fFirstSpecies];
		}
//		for ( int eqLoop = 0; eqLoop < M; ++eqLoop ){
//			a[eqLoop][fTemperature] += dfdyUpwind( eqLoop, fTemperature, EntFluxFuncDiffSim, nodeInfo, flame ) * hnenn;
//		}
#endif
#endif
	}
}

Double EntFluxFuncDiffCont( int /*equation*/, NodeInfoPtr nodeInfo, void *object, Flag /*theFlag*/  )
{
	Double	entFlux = 0.0;
#ifdef ENTHALPYFLUX
	TCountDiffFlameContPtr	flame = ( TCountDiffFlameContPtr )object;
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
#	ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			sumCpY += heatCapacity[eqLoop] * y[speciesEq];
		}
#	endif
	}
	entFlux = sumCpDdYdx;
#	ifdef DIFFUSIVITYCORRECTION
	if ( flame->UseDiffCorr() ) {
		entFlux -= diffCorr[kCurr] * sumCpY;
	}
#	endif
	entFlux *= entCoeff;
#endif
	
	return entFlux;
}

void CountDiffContRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
	TCountDiffFlameContPtr	flame = ( TCountDiffFlameContPtr )object;
	if ( !flame->RHSAction( nodeInfo, rhsMode ) ) {
		return;
	}
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fInvStrainRate = flame->GetOffsetInvStrainRate();
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
	Double	*temp = flameNode->temp;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*enthalpy = flameNode->enthalpy;
	Double	*mixDensity = flameNode->mixDensity;
	Double	mixDensityInf = flameNode->rhoInf;
	Double	mixViscosityInf = flameNode->viscosityInf;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*productionRate = flameNode->productionRate;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	oneOverATimesRho = y[fInvStrainRate] / mixDensity[kCurr];
	Double	idealGasCoeff = flame->GetPressure() * *flameNode->mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = y[fInvStrainRate] / idealGasCoeff;
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	mixHeatCapacity = *flameNode->mixHeatCapacity;
	Double	constThermDiffCoeff = constMassDiffCoeff / flameNode->mixHeatCapacity[kCurr];
	Double	sumCpDdYdx;
	Double	sumMH;
	
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
		rhs[speciesEq] += constMassDiffCoeff * flame->SpeciesDiffusion( speciesEq, nodeInfo );
#ifdef MOLARDIFFUSION
		rhs[speciesEq] += constMassDiffCoeff * flame->SecondDerivXDiffusion( speciesEq, nodeInfo );
#endif
#ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			rhs[speciesEq] -= constMassDiffCoeff * flame->NewDiffCorr( speciesEq, nodeInfo );
#	ifdef MOLARDIFFUSION
			rhs[speciesEq] -= flame->GetSolver()->bt->GetParameter() * constMassDiffCoeff * flame->NewDiffCorrX( speciesEq, nodeInfo );
#	endif
		}
#endif
		if ( flame->fThermoDiffusion ) {
			rhs[speciesEq] += constMassDiffCoeff * flame->ThermoDiffusion( eqLoop, kSimilarity, nodeInfo );
		}
		rhs[speciesEq] += productionRate[eqLoop] * y[fTemperature] * ROverAPM;
	}
	
// inverse strainrate equation
	if ( nodeInfo->gridPoint < flame->fTmaxLoc ) {
		rhs[fInvStrainRate] += FirstDerivUpwind( yNext[fInvStrainRate], y[fInvStrainRate], h );
	}
	else if ( nodeInfo->gridPoint == flame->fTmaxLoc ) {
		rhs[fInvStrainRate] += flame->fdTds * ( temp[kCurr] - flame->fTmaxStart )
								+ flame->fdAInvds * ( y[fInvStrainRate] - flame->fAInvStart )
								- ( flame->fArcLength - flame->fArcLengthStart );
	}
	else {
		rhs[fInvStrainRate] += FirstDerivUpwind( y[fInvStrainRate], yPrev[fInvStrainRate], hm );
	}

// energy equation
	if ( fTemperature < M ) {
#ifdef ENTHALPYFLUX
#	ifdef MOLARDIFFUSION
		Double	sumCpYD = 0.0;
#	endif
#endif
		Double	diffCorrHeat = 0.0;
#	ifdef DIFFUSIVITYCORRECTION
		diffCorrHeat = ( flame->UseDiffCorr() ) ? flameNode->mixHeatCapacity[kCurr] : 0.0;
#	endif
		Double	oneOverARhoCp = oneOverATimesRho / mixHeatCapacity;
		sumCpDdYdx = 0.0;
		sumMH = 0.0;

		rhs[fTemperature] += constThermDiffCoeff * flame->StandardDiffusion( fTemperature, flameNode->mixConductivity, nodeInfo );
#ifdef DIFFUSIVITYCORRECTION
/*		if ( flame->UseDiffCorr() ) {*/
/*			corrCoeff = constMassDiffCoeff * mixDensity[kCurr] * mixDensity[kCurr]*/
/*									* diffCorr[kCurr] / mixHeatCapacity;*/
/*		}*/
#endif

		for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
			sumCpDdYdx += ( heatCapacity[eqLoop] - diffCorrHeat ) * diffusivity[eqLoop] 
						* FirstDeriv( YPrev[eqLoop], Y[eqLoop], YNext[eqLoop], hm, h );
#ifdef ENTHALPYFLUX
#	ifdef MOLARDIFFUSION
			sumCpYD += ( heatCapacity[eqLoop] - diffCorrHeat ) * Y[eqLoop] * diffusivity[eqLoop];
#	endif
#endif
			sumMH += productionRate[eqLoop] * enthalpy[eqLoop];
		}
#ifdef ENTHALPYFLUX
		rhs[fTemperature] +=  constThermDiffCoeff * sumCpDdYdx * mixDensity[kCurr] * mixDensity[kCurr]
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#		ifdef MOLARDIFFUSION
		rhs[fTemperature] +=  constThermDiffCoeff * sumCpYD * mixDensity[kCurr] * mixDensity[kCurr]
						/ flameNode->mixMolarMass[kCurr]
						* FirstDeriv( flameNode->mixMolarMass[kPrev], flameNode->mixMolarMass[kCurr], flameNode->mixMolarMass[kNext], hm, h )
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#		endif
#endif
		rhs[fTemperature] -= sumMH * oneOverARhoCp;
		
		if ( flame->fProperties->GetRadiation() ) {
			rhs[fTemperature] += flameNode->radiation[kCurr] * oneOverARhoCp;
		}
		
	}

	for ( eqLoop = 0; eqLoop < M; ++eqLoop ){
		rhs[eqLoop] *= - hnenn;
	}
}

Double TCountDiffFlameCont::SecondDerivXDiffusion( int nVariable, NodeInfoPtr nodeInfo )
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

Double TCountDiffFlameCont::NewDiffCorrX( int nVariable, NodeInfoPtr nodeInfo )
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

Double EntFluxFunc( int /*j*/, NodeInfoPtr nodeInfo, void *object, Flag /*theFlag*/ )
{
	TCountDiffFlameContPtr	flame = ( TCountDiffFlameContPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		eqLoop;
	int 	fTemperature = flame->GetOffsetTemperature();
	Double	*diffusivity = flameNode->diffusivity;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	mixDensity = flameNode->mixDensity[kCurr];
	Double	*YPrev = flameNode->Y[kPrev];
	Double	*Y = flameNode->Y[kCurr];
	Double	*YNext = flameNode->Y[kNext];
	Double	*temp = flameNode->temp;
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
	Double	val;
	Double	sumCpDdYdx = 0.0;
	Double	sumCpY = 0.0;
	Double	coeff = mixDensity * mixDensity / ( flameNode->mixHeatCapacity[kCurr] 
					* flameNode->rhoInf * flameNode->viscosityInf );
	
	flame->ComputeDiffusivityCorrection( flameNode->Y, nodeInfo );
	for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
		sumCpDdYdx += heatCapacity[eqLoop] * diffusivity[eqLoop] 
					* FirstDeriv( YPrev[eqLoop], Y[eqLoop], YNext[eqLoop], hm, h );
#ifdef DIFFUSIVITYCORRECTION
		if ( flame->UseDiffCorr() ) {
			sumCpY += heatCapacity[eqLoop] * Y[eqLoop];
		}
#endif
	}

	val =  coeff * sumCpDdYdx
						* FirstDeriv( temp[kPrev], temp[kCurr], temp[kNext], hm, h );

#ifdef DIFFUSIVITYCORRECTION
	if ( flame->UseDiffCorr() ) {
		val -= coeff * flameNode->diffCorr[kCurr] * sumCpY
						* FirstDeriv( temp[kPrev], temp[kCurr], temp[kNext], hm, h );
	}
#endif
	return val;
}

int CountDiffContPostIter( void *object )
{
	TCountDiffFlameContPtr	flame = ( TCountDiffFlameContPtr )object;
	int			fUVelocity = flame->GetOffsetUVelocity();
	int			fVVelocity = flame->GetOffsetVVelocity();
	int			fTemperature = flame->GetOffsetTemperature();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();
	int			fMixFrac = flame->GetOffsetMixFrac();
	int			fInvStrainRate = flame->GetOffsetInvStrainRate();
	int			i;
	TNewtonPtr	bt = flame->GetSolver()->bt;
	NodeInfoPtr nodeInfo = bt->GetNodeInfo();
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = currGrid->GetNGridPoints();
	int			*bcFlagLeft = currGrid->GetBcFlagLeft();
	int			mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
	MatrixPtr	yMat = currGrid->GetY();
	Double		**y = yMat->mat;
	Double		*x = currGrid->GetX()->vec;
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
		if ( flame->CheckSolution( y[i][fTemperature], &y[i][fFirstSpecies], nSpeciesInSystem ) ) {
			return 1;
		}
		if ( y[i][fInvStrainRate] <= 0.0 ) {
			y[i][fInvStrainRate] = flame->fMinStrainRate;
		}
	}


	if ( !bt->GetGrid()->IsFine() ) {
		// don't use fInvStrainRate for grid generation
		TGridPtr	fine = bt->GetGrid()->GetFine();
		Double		**yFine = fine->GetY()->mat;
		for ( i = 0; i < nGridPoints; ++i ) {
			y[i][fInvStrainRate] = yFine[2*i+1][fInvStrainRate];
		}
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
	
	yLeft[fInvStrainRate] = y[0][fInvStrainRate];

// right boundary 
	yRight[fUVelocity] = 1.0;
	yRight[fVVelocity] = yLast[fVVelocity] + hLast * yRight[fUVelocity];
	yRight[fInvStrainRate] = yLast[fInvStrainRate];

	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	
	flame->fTmaxLoc = LocationOfMax( nGridPoints+2, &temp[kPrev] ) - 1;

	flame->UpdateThermoProps();
	cerr << "a = " << 1.0 / yLeft[fInvStrainRate] << TAB;

	return 0;
}

#include "TofZ.h"

void TCountDiffFlameCont::UpdateDimensions( int len )
{
	T1DFlame::UpdateDimensions( len );
	fSolV->len = len;
	fSolU->len = len;
	fSolMixFrac->len = len;
	fSolInvStrainRate->len = len;
}

void TCountDiffFlameCont::UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
{
	int		nGridPoints = yMat->cols;
	Double	*V = fSolV->vec;
	Double	*U = fSolU->vec;
	Double	*Z = fSolMixFrac->vec;
	Double	*oneOverA = fSolInvStrainRate->vec;
	Double	**y = yMat->mat;
	Double	*yLeft = yLeftVec->vec;
	Double	*yRight = yRightVec->vec;

	UpdateDimensions( nGridPoints );

	T1DFlame::UpdateSolution( yMat, yLeftVec, yRightVec );
	
	V[kPrev] = yLeft[fVVelocity];
	U[kPrev] = yLeft[fUVelocity];
	Z[kPrev] = yLeft[fMixFrac];
	oneOverA[kPrev] = yLeft[fInvStrainRate];
	for ( int k = 0; k < nGridPoints; ++k ) {
		V[k] = y[k][fVVelocity];
		U[k] = y[k][fUVelocity];
		Z[k] = y[k][fMixFrac];
		oneOverA[k] = y[k][fInvStrainRate];
	}
	V[nGridPoints] = yRight[fVVelocity];
	U[nGridPoints] = yRight[fUVelocity];
	Z[nGridPoints] = yRight[fMixFrac];
	oneOverA[nGridPoints] = yRight[fInvStrainRate];
}

void TCountDiffFlameCont::UpdateSolutionOnePoint( Double *y, int gridPoint )
{
	T1DFlame::UpdateSolution( y, gridPoint );
	
	fSolU->vec[gridPoint] = y[fUVelocity];
	fSolV->vec[gridPoint] = y[fVVelocity];
	fSolMixFrac->vec[gridPoint] = y[fMixFrac];
	fSolInvStrainRate->vec[gridPoint] = y[fInvStrainRate];
}

void TCountDiffFlameCont::SolutionToSolver( void )
{
	TNewtonPtr	bt = fSolver->bt;
	TGridPtr	grid = bt->GetGrid()->GetFine();
	int		nGridPoints = fSolV->len;
	Double	*V = fSolV->vec;
	Double	*U = fSolU->vec;
	Double	*Z = fSolMixFrac->vec;
	Double	*oneOverA = fSolInvStrainRate->vec;
	Double	**y = grid->GetY()->mat;

	T1DFlame::SolutionToSolver();
	
	for ( int k = 0; k < nGridPoints; ++k ) {
		y[k][fVVelocity] = V[k];
		y[k][fUVelocity] = U[k];
		y[k][fMixFrac] = Z[k];
		y[k][fInvStrainRate] = oneOverA[k];
	}
	
	CountDiffContPostIter( this );
}

void TCountDiffFlameCont::SaveSolution( void )
{
	int		k;
	int		len = fSolV->len;
	Double	*v = fSolV->vec;
	Double	*saveV = fSavedV->vec;
	Double	*u = fSolU->vec;
	Double	*saveU = fSavedU->vec;
	Double	*mixFrac = fSolMixFrac->vec;
	Double	*saveMixFrac = fSavedMixFrac->vec;
	Double	*oneOverA = fSolInvStrainRate->vec;
	Double	*saveOneOverA = fSavedInvStrainRate->vec;

	T1DFlame::SaveSolution();
	fSavedV->len = fSolV->len;
	fSavedU->len = fSolU->len;
	fSavedMixFrac->len = fSolMixFrac->len;
	fSavedInvStrainRate->len = fSolInvStrainRate->len;

	for ( k = -1; k <= len; ++k ) {
		saveV[k] = v[k];
		saveU[k] = u[k];
		saveMixFrac[k] = mixFrac[k];
		saveOneOverA[k] = oneOverA[k];
	}
	fSaveddTds = fdTds;
	fSaveddAInvds = fdAInvds;
}

void TCountDiffFlameCont::RestoreSolution( void )
{
	int		k;
	int		len = fSolV->len;
	Double	*v = fSolV->vec;
	Double	*saveV = fSavedV->vec;
	Double	*u = fSolU->vec;
	Double	*saveU = fSavedU->vec;
	Double	*mixFrac = fSolMixFrac->vec;
	Double	*saveMixFrac = fSavedMixFrac->vec;
	Double	*oneOverA = fSolInvStrainRate->vec;
	Double	*saveOneOverA = fSavedInvStrainRate->vec;

	UpdateDimensions( len );

	T1DFlame::RestoreSolution();

	for ( k = -1; k <= len; ++k ) {
		v[k] = saveV[k];
		u[k] = saveU[k];
		mixFrac[k] = saveMixFrac[k];
		oneOverA[k] = saveOneOverA[k];
	}
	fdTds = fSaveddTds;
	fdAInvds = fSaveddAInvds;

	SolutionToSolver();
}

int	TCountDiffFlameCont::GetOffsetVVelocity( void )
{
	return fVVelocity; 
}

int	TCountDiffFlameCont::GetOffsetUVelocity( void )
{
	return fUVelocity; 
}

int	TCountDiffFlameCont::GetOffsetTemperature( void )
{
	return fTemperature; 
}

int TCountDiffFlameCont::GetOffsetMixFrac( void )
{
	return fMixFrac;
}

int TCountDiffFlameCont::GetOffsetInvStrainRate( void )
{
	return fInvStrainRate;
}

int	TCountDiffFlameCont::GetOffsetFirstSpecies( void ) 
{
	return fFirstSpecies;
}

ConstStringArray TCountDiffFlameCont::GetVariableNames( void )
{
	return (ConstStringArray)fVariableNames;
}

int TCountDiffFlameCont::GetVariablesWithoutSpecies( void )
{
	return fVariablesWithoutSpecies;
}

void TCountDiffFlameCont::SetInitialValues( TInputDataPtr inp, TSolutionPtr sol, Double strainRate )
{
	int 				i, j, k;
	TBVPSolverPtr		solver = GetSolver();
	TNewtonPtr			bt = solver->bt;
	TAdaptiveGridPtr	adapGrid = bt->GetGrid();
	TGridPtr			grid = adapGrid->GetFine();
	int					variables = bt->GetNVariables();
	int					nGridPoints = sol->fCoord->len;
	MatrixPtr			yMat = grid->GetY();
	Double				**y = yMat->mat;
	VectorPtr			yLeftVec = grid->GetYLeft();
	VectorPtr			yRightVec = grid->GetYRight();
	Double			 	*yLeft = yLeftVec->vec;
	Double 				*yRight = yRightVec->vec;
	VectorPtr			xVec = grid->GetX();
	Double				*x = xVec->vec;
	FILE				*fp;
	Flag 				found;
	
	SetPressure( sol->fPressure );

	int		varsIn = sol->fNOfSolNames;
	char	**namesIn = sol->fSolNames;
	Double	*xIn = sol->fCoord->vec;
	Double	*yLeftIn = sol->fSolLeft->vec;
	Double	*yRightIn = sol->fSolRight->vec;
	Double	**yIn = sol->fSol->mat;

	Double	oneOverA = 1.0 / strainRate;
	for ( k = 0; k < nGridPoints; ++k ) {
		y[k][fInvStrainRate] = oneOverA; 
	}
	fAInvOld = fAInvStart = oneOverA;
//	fArcLength = fDeltaS;
	fArcLength = 0.0;
	fArcLengthStart = 0.0;

	grid->AdjustNGridPoints( nGridPoints );
	solver->UpdateAllDimensions( nGridPoints );	
	for ( k = 0; k < nGridPoints; ++k ) {
		x[k] = xIn[k];
	}

	bt->SetLeft( sol->fLeft );
	bt->SetRight( sol->fRight );

	for ( i = 0; i < varsIn; ++i ) { // loop input
		// find variable
		found = FALSE;
		for ( j = 0; j < variables; ++j ) {
			if ( strcmp( namesIn[i], fVariableNames[j] ) == 0 ) {
				found = TRUE;
				break;
			}
		}
		if ( !found ) {
			cerr << "#error: variable " << namesIn[i] << " couldn't be found" << NEWL;
			exit( 2 );
		}
		// else found, now copy
		yLeft[j] = yLeftIn[i];
		yRight[j] = yRightIn[i];
		yMat->cols = sol->fSol->cols;
		for ( k = 0; k < nGridPoints; ++k ) {
			y[k][j] = yIn[k][i];
		}
	}
	
//	update properties
	Double	*temp = GetTemperature()->vec;
	Double	**Y = GetMassFracs()->mat;

	UpdateSolution( yMat, yLeftVec, yRightVec );
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
		yLeft[fVVelocity] = y[0][fVVelocity] - ( x[0] - bt->GetLeft() ) 
							* yLeft[fUVelocity];
	}

	if ( GetSpecies()->IsConstantLewisNumber() ) {
		CompLewisNumbers( GetSpecies()->GetLewisNumberFile() );
	}
	CountDiffContPostIter( this );
	bt->GetGrid()->SetSolutionScaler();

	fp = GetOutfile( "initialguess2", TFlame::kData );
	bt->PrintSolution( x, y, GetVariableNames(), fp );
	fclose(fp);

	SaveSolution();
}

void TCountDiffFlameCont::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp )
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
	if ( !xIn ) FatalError( "memory allocation of TCountDiffFlameCont failed" );
	Double				*yIn =  new Double[gridPointsIn];
	if ( !yIn ) FatalError( "memory allocation of TCountDiffFlameCont failed" );
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
	struct _parameter	*param = GetParameter( "strainrate" );
//	Flag				invStrainRateSet = FALSE;

// get strainrate
	if ( param ) {
		fAInvStart = 1.0 / (Double)param->what.quantity.value;
	}
	else { // choose default
		cerr << "#warning: no value for 'strainrate' in inputfile" << NEWL;
		fAInvStart = 1.0 / 100.0;
	}
// get Tmax
	param = GetParameter( "tmax" );
	if ( param ) {
		fTmaxStart = (Double)param->what.quantity.value;
	}
	else { // choose default
		cerr << "#warning: no value for 'Tmax' in inputfile" << NEWL;
		fTmaxStart = 2000.0;
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
//		else if ( strncmp( string, "InvStrainRate", 11 ) == 0 ) {
//			variable = fInvStrainRate;
//			invStrainRateSet = TRUE;
//		}
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
	
	Double	oneOverA = fAInvStart;
	for ( k = 0; k < nGridPoints; ++k ) {
		y[k][fInvStrainRate] = oneOverA; 
	}
//	fArcLength = fDeltaS;
	fArcLength = 0.0;
	fAInvOld = fAInvStart;
	fArcLengthStart = 0.0;
/*	fdTds = -1.0;
	fdAInvds = -1.0e-5;*/
/*	fdTds = -0.5;
	fdAInvds = -0.5;*/
/*	fdTds = -1.0;
	fdAInvds = -0.00227;*/

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
	if ( GetSpecies()->IsConstantLewisNumber() ) {
		CompLewisNumbers( GetSpecies()->GetLewisNumberFile() );
	}
	CountDiffContPostIter( this );

	FreeSpline( theSpline );
	delete yIn;
	delete xIn;

	adapGrid->SetSolutionScaler();
	
	fp = GetOutfile( "initialguess", TFlame::kData );
	bt->PrintSolution( locX, y, GetVariableNames(), fp );
	fclose(fp);
	
	SaveSolution();
}

void TCountDiffFlameCont::FillJacSpeciesDiffusion( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
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

Double TCountDiffFlameCont::SpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo )
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

void TCountDiffFlameCont::PrintRHSTemp( TNewtonPtr bt )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		k;
    int         		N = currentGrid->GetNGridPoints();
	char				fName[32];
	FILE				*fp = NULL;
	
#if defined (applec) || defined (powerc)
    RotateCursor( 32 * bt->GetNIter() );
#endif
	
	UpdateThermoProps();
	
	sprintf( fName, "tempeq%.2d", bt->GetNIter() );
	if ( !( fp = fopen( fName, "w" ) ) ) { 
		cerr << "#warning: unable to open file " << fName << NEWL;
		exit(2);
	}
	
	fprintf( fp, "*\n%-12s\t%-12s\t%-12s\t%-12s\t%-12s\n", "x", "convection", "diffusion", "enthalpyFlux", "production" );
	
	for ( k = 0; k < N; ++k ){
#if defined (applec) || defined (powerc)
	RotateCursor( 32 );
#endif
		bt->SetNodeInfo( this, k );
		PrintRHSTemp( bt, nodeInfo, fp );
	}
    fclose( fp );
}

void TCountDiffFlameCont::PrintRHSTemp( TNewtonPtr bt, NodeInfoPtr nodeInfo, FILE *fp )
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
	Double	mixDensityInf = fFlameNode->rhoInf;
	Double	mixViscosityInf = fFlameNode->viscosityInf;
	Double	*productionRate = fFlameNode->productionRate;
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*heatCapacity = fFlameNode->heatCapacity;
	Double	idealGasCoeff = GetPressure() * *fFlameNode->mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = 1.0 / ( strainRate * idealGasCoeff );
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	mixHeatCapacity = *fFlameNode->mixHeatCapacity;
	Double	constThermDiffCoeff = constMassDiffCoeff / fFlameNode->mixHeatCapacity[kCurr];
	Double	sumCpDdYdx;
	Double	sumMH;

	fprintf( fp, "%-.6e", *nodeInfo->x );
#ifdef UPWINDCONVECTION
	fprintf( fp, "\t%-.6e", NonlinearConvectUpwind( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) );
#else
	fprintf( fp, "\t%-.6e", NonlinearConvectCentral( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h ) );
#endif

// energy equation
	sumCpDdYdx = 0.0;
	sumMH = 0.0;

	fprintf( fp, "\t%-.6e", constThermDiffCoeff * StandardDiffusion( fTemperature, fFlameNode->mixConductivity, nodeInfo ) );
	for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
		speciesEq = fFirstSpecies+eqLoop;
		sumCpDdYdx += heatCapacity[eqLoop] * diffusivity[eqLoop] 
					* FirstDeriv( yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
		sumMH += productionRate[eqLoop] * enthalpy[eqLoop];
	}
	fprintf( fp, "\t%-.6e", constThermDiffCoeff * sumCpDdYdx / ( y[fTemperature] * y[fTemperature] )
					* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h )
					* idealGasCoeff * idealGasCoeff );
	fprintf( fp, "\t%-.6e", sumMH * y[fTemperature] * ROverAPM / mixHeatCapacity );
	fprintf( fp, "\n" );
}

void TCountDiffFlameCont::PrintRHSSpecies( TNewtonPtr bt )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		i, k;
    int         		N = currentGrid->GetNGridPoints();
	VectorPtr 			physXVec = NewVector( N + 2 );
	Double				*physX = physXVec->vec;
	char				fName[32];
	char				**names = GetSpecies()->GetNames();
	int					nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	FILE				*fp = NULL;
	
#if defined (applec) || defined (powerc)
    RotateCursor( 32 * bt->GetNIter() );
#endif
	
	UpdateThermoProps();
	
	EtaToX( bt, physXVec );
	
	sprintf( fName, "specieseq%.2d", bt->GetNIter() );
	if ( !( fp = fopen( fName, "w" ) ) ) { 
		cerr << "#warning: unable to open file " << fName << NEWL;
		exit(2);
	}
	
	fprintf( fp, "*\n%-12s\t%-12s", "eta", "y" );
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fprintf( fp, "\tConv_%-7s\tDiff_%-7s\tProd_%-7s\tCons_%-7s", names[i], names[i], names[i], names[i] );
	}
	fprintf( fp, "\n" );
		
	for ( k = 0; k < N; ++k ){
#if defined (applec) || defined (powerc)
	RotateCursor( 32 );
#endif
		bt->SetNodeInfo( this, k );
		PrintRHSSpecies( bt, nodeInfo, physX[k+1], fp );
	}
    fclose( fp );
}

void TCountDiffFlameCont::PrintRHSSpecies( TNewtonPtr bt, NodeInfoPtr nodeInfo, Double physX, FILE *fp )
{
	int		i, j;
	int		speciesEq;
    int		M = bt->GetNEquations();
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	strainRate = GetStrainRate();
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	mixDensityInf = fFlameNode->rhoInf;
	Double	mixViscosityInf = fFlameNode->viscosityInf;
	Double	idealGasCoeff = GetPressure() * *fFlameNode->mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = 1.0 / ( strainRate * idealGasCoeff );
	Double	constMassDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );
	Double	diffTerm;
	Double	*diffCorr = fFlameNode->diffCorr;
	Double	*reactionRate = fFlameNode->reactionRate;
	Double	productionRate;
	Double	source;
	Double	sink;
	int		*nOfUsedReactions = fSpecies->GetNOfUsedReactions()->vec;
	VectorPtr	*nu = fSpecies->GetNu();
	IntVectorPtr	*usedReactions = fSpecies->GetUsedReactions();
	Double	*molarMass =  fSpecies->GetMolarMass()->vec;

	fprintf( fp, "%-.6e\t%-.6e", *nodeInfo->x, physX );

	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		speciesEq = fFirstSpecies+i;

#ifdef UPWINDCONVECTION
		fprintf( fp, "\t%-.6e", NonlinearConvectUpwind( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h ) );
#else
		fprintf( fp, "\t%-.6e", NonlinearConvectCentral( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h ) );
#endif

		diffTerm = constMassDiffCoeff * SpeciesDiffusion( speciesEq, nodeInfo );
#ifdef DIFFUSIVITYCORRECTION
		diffTerm -= constMassDiffCoeff * FirstDeriv( mixDensity[kPrev] * mixDensity[kPrev] * yPrev[speciesEq] * diffCorr[kPrev], 
														  mixDensity[kCurr] * mixDensity[kCurr] * y[speciesEq] * diffCorr[kCurr],
														  mixDensity[kNext] * mixDensity[kNext] * yNext[speciesEq] * diffCorr[kNext], hm, h );
#endif
		fprintf( fp, "\t%-.6e", diffTerm );

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
		fprintf( fp, "\t%-.6e\t%-.6e", source, -sink );
	}
	fprintf( fp, "\n" );
}

void CountDiffContOutput( void *object, FILE *fp, char* tail )
{
	TCountDiffFlameContPtr	flame = ( TCountDiffFlameContPtr )object;
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
	Double			*oneOverA = flame->GetOneOverA()->vec;
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
	int				fInvStrainRate = flame->GetOffsetInvStrainRate();
	time_t			theDate;
	char			buffer[80];
	ConstStringArray	varNames = flame->GetVariableNames();
	char			**names = species->GetNames();
	VectorPtr 		physXVec = NewVector( gridPoints + 2 );
	Double			*physX = physXVec->vec;
	VectorPtr 		scalarDissVec = NewVector( gridPoints + 2 );
	Double			*scalarDiss = scalarDissVec->vec;
	Double			stoechScalarDiss;
	Flag			fpOpen = FALSE;
	
	if ( !fp ) {
		fpOpen = TRUE;
		fp = flame->GetOutputFile( NULL, tail, TFlame::kNone );
	}

	flame->ScalarDissipation( bt, scalarDissVec, &stoechScalarDiss );
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
		fprintf( fp, "\tfuel = \"%s\"\n", varNames[firstSpecies+flame->GetFuelIndex( i )] );
	}
	fprintf( fp, "pressure = %g [bar]\n", flame->GetPressure() / 1.0e5 );
	fprintf( fp, "strainrate = %g [1/s]\n", flame->GetStrainRate() );

	if ( species->IsConstantLewisNumber() ) {
		fprintf( fp, "ConstantLewisNumbers = \"True\"\n" );
	}
	
	fprintf( fp, "Tmax = %g [K]\n", temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
	fprintf( fp, "stoichScalarDissRate = %g [1/s]\n", stoechScalarDiss );
	
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

//	write Z Bilger
	fprintf( fp, "ZBilger\n" );
	fprintf( fp, "\t%-.6e", Z[kPrev] );
	for ( k = 0; k < gridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", flame->ComputeZBilger( massFracs[k], massFracs[-1], massFracs[gridPoints] ) );
//		fprintf( fp, "\t%-.6e", flame->ComputeZBilger( &y[k][firstSpecies], &yLeft[firstSpecies], &yRight[firstSpecies] ) );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", Z[gridPoints] );
		
	fprintf( fp, "trailer\n" );
	if ( nOfEquations < nOfVariables) {
		fprintf( fp, "number of converged equations is %d\n", nOfEquations );
	}
	
	DisposeVector( scalarDissVec );
	DisposeVector( physXVec );
	
	if ( fpOpen ) {
		fclose( fp );
	}
}

void TCountDiffFlameCont::FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
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

Double TCountDiffFlameCont::NewDiffCorr( int nVariable, NodeInfoPtr nodeInfo )
{
// returns     sum_j ( d/dy(rho^2 Y_k D_j dY_j/dy) )

	int		i;
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
	i = nVariable-fFirstSpecies;
	Double	coeffCurr = density[kCurr] * density[kCurr] * Y[i];
	Double	coeffPrev = density[kPrev] * density[kPrev] * YPrev[i];
	Double	coeffNext = density[kNext] * density[kNext] * YNext[i];
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

void TCountDiffFlameCont::ScalarDissipation( TNewtonPtr bt, VectorPtr scalarDissVec, Double *stoechScalarDiss )
{
	int				k;
	int				nOfGridPoints = bt->GetCurrentGridPoints();
	TGridPtr		grid = bt->GetGrid()->GetCurrentGrid();
	NodeInfoPtr		nodeInfo = bt->GetNodeInfo();
	Double			**yMat = grid->GetY()->mat;
	Double			left = bt->GetLeft();
	Double			right = bt->GetRight();
	Double			*x = grid->GetX()->vec;
	Double			zLeft = grid->GetYLeft()->vec[fMixFrac];
	Double			zRight = grid->GetYRight()->vec[fMixFrac];
	Double			*scalarDiss = scalarDissVec->vec;
	Double			*lambda = fProperties->GetConductivity()->vec;
	Double			*density = fProperties->GetDensity()->vec;
	Double			*cp = fProperties->GetHeatCapacity()->vec;
	Double			rhoInf = fFlameNode->rhoInf;
	Double			viscosityInf = fFlameNode->viscosityInf;
	Double			factor = 2.0 * GetStrainRate() / ( rhoInf * viscosityInf );
	Double			dZdx;
	Double			zStoe = fMassFraction->GetZStoe();
	
	dZdx = ( yMat[0][fMixFrac] - zLeft ) / ( x[0] - left );
	scalarDiss[0] = factor * density[0]  * lambda[0] / cp[0]  * dZdx * dZdx;
	for ( k = 0; k < nOfGridPoints; ++k ) {
		bt->SetNodeInfo( this, k );
		dZdx = FirstDeriv( nodeInfo->yPrev[fMixFrac], nodeInfo->y[fMixFrac], nodeInfo->yNext[fMixFrac], nodeInfo->hm, nodeInfo->h );
		scalarDiss[k+1] = factor * density[k]  * lambda[k] / cp[k]  * dZdx * dZdx;
		if ( ( nodeInfo->y[fMixFrac] - zStoe ) * ( nodeInfo->yPrev[fMixFrac] - zStoe ) <= 0.0 ) {
			*stoechScalarDiss = scalarDiss[k] + ( scalarDiss[k+1] - scalarDiss[k] ) 
							/ ( nodeInfo->y[fMixFrac] - nodeInfo->yPrev[fMixFrac] ) * ( zStoe - nodeInfo->yPrev[fMixFrac] );
		}
	}
	dZdx = ( zRight - yMat[nOfGridPoints-1][fMixFrac] ) / ( right - x[nOfGridPoints-1] );
	scalarDiss[nOfGridPoints+1] = factor * density[nOfGridPoints-1]  * lambda[nOfGridPoints-1] / cp[nOfGridPoints-1]  * dZdx * dZdx;
}

void SetCountDiffContNodeInfo( int k, void *object )
{
	TCountDiffFlameContPtr	flame = ( TCountDiffFlameContPtr )object;
	
	flame->SetFlameNode( k );
}

void CountDiffContPostConv( void *object )
{
	TCountDiffFlameContPtr	flame = ( TCountDiffFlameContPtr )object;
	TNewtonPtr 				bt = flame->GetSolver()->bt;
	int						isConverged = bt->GetConvergeNewton();
	
	if ( isConverged ) {
		flame->SaveSolution();

		if( flame->AdjustGrid( CountDiffContPostIter ) ) {
//			bt->WriteOutput( object, NULL, "Map" );
			return;
		}

		int i, k;
		TGridPtr 	currentGrid = bt->GetGrid()->GetCurrentGrid();
		int			currentGridPoints = currentGrid->GetNGridPoints();
		int			equations = bt->GetNEquations();
		Double		**y = currentGrid->GetY()->mat;
		MatrixPtr	dy = bt->GetDy();
		Double		**locDy = dy->mat;
		Double		ds = flame->GetDeltaS();
		Double		hnenn;
		Double		aInv;
		Flag		unPhysical;
		static int	unPhysCounter = 0;
				
		//	init
		flame->GetSolver()->ReInit();
/*		bt->InitNIter();
		bt->UnSetLeaveContin();
		bt->GetGrid()->UnSetSolHasNewGrid();*/

		//	output
		aInv = y[0][flame->GetOffsetInvStrainRate()];
		if ( flame->fArcLength - flame->fArcLengthStart ) {
			unPhysical = ( ( aInv - flame->fAInvOld ) / flame->fDeltaS > 0.0 ) ? TRUE : FALSE;
		}
		else {
			unPhysical = FALSE;
		}
		flame->fAInvOld = aInv;
		
		Double	tMax = y[flame->fTmaxLoc][flame->GetOffsetTemperature()];
		fprintf( flame->fpExtCurve, "%-12E\t%-12E\t%-12E\n", 1.0 / aInv, aInv, tMax );
		fflush( flame->fpExtCurve );
		if ( unPhysical ) {
			bt->WriteOutput( object, NULL, "noPhys" );
			unPhysCounter++;
			if ( !flame->ComputeUnphysicalChain() || unPhysCounter >= flame->fNUnPhysChain ) {
				bt->SetLeaveContin();
			}
		}
		else {
			bt->WriteOutput( object, NULL, "" );
		}

		if ( !unPhysical || ( unPhysical && flame->ComputeUnphysicalChain() 
				&& flame->GetStrainRate() > flame->GetMinStrainRate() 
				&& unPhysCounter < flame->fNUnPhysChain ) ) {

			//	fill rhs
			bt->SetNodeInfo( flame, flame->fTmaxLoc );
			hnenn = bt->GetNodeInfo()->hnenn;
			ClearMatrix( dy );
			dy->mat[flame->fTmaxLoc][flame->GetOffsetInvStrainRate()] = 1.0 * hnenn;
					
			//	BackSolve
			bt->BackSolve( dy, NULL, FALSE );
			
			//	set new arclength
			flame->fSConverged = flame->fArcLength;
			if ( flame->fSNotConverged < 0.0 ) {
/*				if ( flame->fUnPhysFound ) {
					flame->SetDeltaS( 0.5 * ds );
					ds = flame->GetDeltaS();
				}*/
				flame->fArcLength += ds;
			}
			else {
				flame->fArcLength = flame->fSNotConverged;
				flame->fSNotConverged = -1.0;
			}
			
			flame->fdTds = locDy[flame->fTmaxLoc][flame->GetOffsetTemperature()];
			flame->fdAInvds = locDy[flame->fTmaxLoc][flame->GetOffsetInvStrainRate()];
	
			//	update new guess
			for ( k = 0; k < currentGridPoints; ++k ){
				for ( i = 0; i < equations; ++i ){
					y[k][i] += locDy[k][i] * ds;
				}
			}		
			CountDiffContPostIter( flame );
	
	#ifdef DEBUGNEWGUESS
			//	write new guess
			FILE	*fpa;
			AMPrintOptions  prnt;
			DefaultAMPOpts( &prnt );
			
			fpa = GetOutfile( "dyds", TFlame::kText );
			PrintMatrix( dy, &prnt, fpa );
			fclose(fpa);
	
			FILE	*fp = NULL;
			fp = GetOutfile( "newguess", TFlame::kData );
			bt->PrintSolution( currentGrid->GetX()->vec, y, flame->GetVariableNames(), fp );
			fclose(fp);
	 #endif
	
			cerr << "ArcLength is now " << flame->fArcLength << TAB
					<< "d/ds(1/a) = " << flame->fdAInvds << NEWL;
		}
	}
	else {
		Double	interS = flame->fSConverged + ( flame->fArcLength - flame->fSConverged ) * 0.5;
		
		flame->RestoreSolution();
		CountDiffContPostIter( flame );
		flame->GetSolver()->ReInit();
/*		bt->InitNIter();
		bt->UnSetLeaveContin();
		bt->GetGrid()->UnSetSolHasNewGrid();*/

		if ( ( interS - flame->fSConverged ) < 1.0 ) {
			bt->SetLeaveContin();
		}
		else {
			flame->fSNotConverged = flame->fArcLength;
			flame->fArcLength = interS;
			cerr << "ArcLength is now " << flame->fArcLength << NEWL;
			bt->WriteOutput( object, NULL, "noC" );
		}
	}
}

ConstStringArray GetCountDiffContVarNames( void *object )
{
	TCountDiffFlameContPtr	flame = ( TCountDiffFlameContPtr )object;
	
	return flame->GetVariableNames();
}

Double TCountDiffFlameCont::GetStrainRate( void )
{
	return 1.0 / fSolInvStrainRate->vec[0];
}

FILE *TCountDiffFlameCont::GetOutputFile( char *head, char *tail, FileType type )
{
	TNewtonPtr		bt = fSolver->bt;
	TGridPtr		currentGrid = bt->GetGrid()->GetCurrentGrid();
	int				fuelIndex = GetFuelIndex();
	int				tempOffset = GetOffsetTemperature();
	char			*name = new char[64];
	FILE			*fp;
	char			**speciesNames = fSpecies->GetNames();
	int				tOxidizer = ( int ) currentGrid->GetYRight()->vec[tempOffset];
	int				tFuel = ( int ) currentGrid->GetYLeft()->vec[tempOffset];
		
	sprintf( name, "%s%s%.8s_p%.2da%.5dtf%.4dto%.4d%s"
					, ( head ) ? head : "", ( head ) ? "_" : ""
					, speciesNames[fuelIndex]
					, ( int ) floor( GetPressure() * 1.0e-5 + 0.5 )	// in [bar]
					, ( int ) floor( GetStrainRate() + 0.5 )			// in [1/s]
					, ( int )( tFuel )							// in [K]
					, ( int )( tOxidizer ) 						// in [K]
					, ( tail ) ? tail : "" );

	fp = GetOutfile( name, type );
	delete name;

	return fp;
}

char **TSolution::CopyIt( char **source, int len )
{
	char **dest = new String[len];

	fNOfSolNames = len;
   	if ( !dest ) FatalError( "memory allocation of char ** in function TSolution::SaveIt failed" );

	for ( int i = 0; i < len; ++i ) {
		dest[i] = new char[strlen( source[i] ) + 1];
		strcpy( dest[i], source[i] );
	}

	return dest;
}

VectorPtr TSolution::CopyIt( VectorPtr source )
{
	VectorPtr	dest = NewVector( source->phys_len );
	Double		*tmp = dest->vec;
	memcpy( dest, source, sizeof( Vector ) );
	dest->vec = tmp;
	copy_vec( dest, source );

	return dest;
}

MatrixPtr TSolution::CopyIt( MatrixPtr source )
{
	MatrixPtr	dest = NewMatrix( source->phys_rows, source->phys_cols, kColumnPointers );
	Double		**tmp = dest->mat;
	memcpy( dest, source, sizeof( Matrix ) );
	dest->mat = tmp;
	copy_mat( dest, source );

	return dest;
}

void TCountDiffFlameCont::EtaToX( TNewtonPtr bt, VectorPtr xPhysVec )
{
	int			k;
	int			gridPoints = bt->GetCurrentGridPoints();
	int			kBefStoech = -1;
	NodeInfoPtr nodeInfo = bt->GetNodeInfo();
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
#undef OLD

#ifdef OLD
	Double		factor = 2.0 * sqrt( fFlameNode->rhoInf * fFlameNode->viscosityInf / GetStrainRate() );
#else
	Double		factor = 0.5 * sqrt( fFlameNode->rhoInf * fFlameNode->viscosityInf / GetStrainRate() );
#endif
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
#ifdef OLD
	xPhys[1] = xPhys[0] + factor * ( x[0] - left ) / ( rho[-1] + rho[0] );
#else
	xPhys[1] = xPhys[0] + factor * ( x[0] - left ) * ( 1.0 / rho[-1] + 1.0 / rho[0] );
#endif
	for ( k = 1; k < gridPoints; ++k ) {
#ifdef OLD
		xPhys[k+1] = xPhys[k] + factor * ( x[k] - x[k-1] ) / ( rho[k-1] + rho[k] );
#else
		xPhys[k+1] = xPhys[k] + factor * ( x[k] - x[k-1] ) * ( 1.0 / rho[k-1] + 1.0 / rho[k] );
#endif
	}
#ifdef OLD
	xPhys[gridPoints+1] = xPhys[gridPoints] + factor * ( right - x[gridPoints-1] ) / ( rho[gridPoints-1] + rho[gridPoints] );
#else
	xPhys[gridPoints+1] = xPhys[gridPoints] + factor * ( right - x[gridPoints-1] ) * ( 1.0 / rho[gridPoints-1] + 1.0 / rho[gridPoints] );
#endif
	
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

void CountDiffContUpdateLeftBoundary( void  *object )
{
	TCountDiffFlameContPtr	flame = ( TCountDiffFlameContPtr )object;
	TNewtonPtr		bt = flame->GetSolver()->bt;
	TGridPtr 		currGrid = bt->GetGrid()->GetCurrentGrid();
	MatrixPtr		yMat = currGrid->GetY();
	Double			*yLeft = currGrid->GetYLeft()->vec;
	Double			pressure = flame->GetPressure();

//	flame->UpdateSolutionOnePoint( yLeft, kPrev );
//	flame->SetFlameNode( kPrev );
//	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
//							, flame->fFlameNode->Y[kCurr], pressure );
}

void CountDiffContUpdateRightBoundary( void *object )
{
	TCountDiffFlameContPtr	flame = ( TCountDiffFlameContPtr )object;
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
