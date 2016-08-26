#include "FlameMaster.h"
#include "ListTool.h"
#include "Spline.h"
#include "MapMan.h"
#include "TUnstrPremFlamePhys.h"

// following allows for full description of diffusion with D_ij
#undef FULLDIFFUSION

#define MOLARDIFFUSION

// following switches between upwind and central discretization of convection terms
#define UPWINDCONVECTION

// following works only if UPWINDCONVECTION is defined
#define MIXEDUPWINDCENTRAL

void TUnstrPremFlamePhys::InitTUnstrPremFlamePhys( void )
{
	int i;
	TBVPSolverPtr	solver = GetSolver();
	TNewtonPtr		bt = solver->bt;
	TGridPtr		fine = bt->GetGrid()->GetFine();
	TGridPtr		coarse = bt->GetGrid()->GetCoarse();
	int				nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	int				maxGridPoints = bt->GetMaxGridPoints();

	fRadFact = 1.0;

	if ( fSoot ) {
		fSoot->SetMomentsOffset( fSootMoments );
	}

//	names of variables
	fVariableNames = new String[fVariablesWithoutSpecies + nSpeciesInSystem];

	fVariableNames[fMassFlowRate] = new char[2];
	strcpy( fVariableNames[fMassFlowRate], "M" );
	fVariableNames[fTemperature] = new char[2];
	strcpy( fVariableNames[fTemperature], "T" );
	fPrintMolarFractions = fInputData->fPrintMolarFractions;
	
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies->GetNames()[i] ) + 1];
		strcpy( fVariableNames[fFirstSpecies + i], fSpecies->GetNames()[i] );
	}

	if ( fSoot ) {
		int	offset = fSoot->GetOffsetSootMoments();
		for ( i = 0; i < fSoot->GetNSootMoments(); ++i ) {
			fVariableNames[offset + i] = new char[8];
			if ( !fVariableNames[offset + i] ) FatalError( "new failed" );
			sprintf( fVariableNames[offset + i], "M%d/rho", i );
		}
	}

//	vectors of solution
	fSolM = NewVector( maxGridPoints + 2 );
	fSolM->vec = &fSolM->vec[kNext];
	fSolM->len -= 2;

//	vectors of saved solution
	fSavedM = NewVector( maxGridPoints + 2 );
	fSavedM->vec = &fSavedM->vec[kNext];
	fSavedM->len -= 2;
	
	fYLeftVec = NewVector( fSpecies->GetNOfSpecies() );
	fRhoY_iV_iPlus = NewVector( fSpecies->GetNSpeciesInSystem() );
	fYLeftVec->len = nSpeciesInSystem;

	fBC_Left = &fine->GetBcLeft()->vec[ fFirstSpecies ];	

	fNewtonInfoL = NewNewtonInfo( fYLeftVec->len, fYLeftVec );
	fNewtonInfoL->modified = FALSE;
	fNewtonInfoL->maxSteps = 100;
	SetNewtonFuncs( fNewtonInfoL, BCLeftNewtonFuncsTUnstr, NULL, NULL );

	fConstMassFlux = fInputData->fCompUnPhysChain;

	fNCutGrid = fNEnlargeGrid = 0;

	if ( fUseNumericalJac ) {
		bt->SetUtFuncs( NULL, NULL, NULL
					, UnstrPremPhysRHSRest, UnstrPremPhysRHSRest, UnstrPremPhysRHSRest 
					, UnstrPremPhysOutput, UnstrPremPhysPostIter
					, SetUnstrPremPhysNodeInfo, UnstrPremPhysPostConv
					, GetUnstrPremPhysVarNames
					, UnstrPremPhysUpdateLeftBoundary, UnstrPremPhysUpdateRightBoundary );
	}
	else {
		bt->SetUtFuncs( UnstrPremPhysJacFirst, UnstrPremPhysJacRest, UnstrPremPhysJacLast
					, UnstrPremPhysRHSRest, UnstrPremPhysRHSRest, UnstrPremPhysRHSRest 
					, UnstrPremPhysOutput, UnstrPremPhysPostIter
					, SetUnstrPremPhysNodeInfo, UnstrPremPhysPostConv
					, GetUnstrPremPhysVarNames
					, UnstrPremPhysUpdateLeftBoundary, UnstrPremPhysUpdateRightBoundary );
	}
	SetInitialBC( fine, fInputData );
	SetInitialBC( coarse, fInputData );
	ReadStartProfiles( fInputData );
	cerr << "initial equivalence ratio is " << GetPhi() << NEWL;
	CheckBC();
	CheckInitialGuess();
	UpdateSolution( fine->GetY(), fine->GetYLeft(), fine->GetYRight() );	
    SaveSolution();
}

void TUnstrPremFlamePhys::SetMassFlux( Double massFlux, MatrixPtr yMat
									, Double *yLeft, Double *yRight )
{
	int		k;
	int		nGridPoints = yMat->cols;
	Double	**y = yMat->mat;

	yLeft[fMassFlowRate] = massFlux;
	for ( k = 0; k < nGridPoints; ++k ) {
		y[k][fMassFlowRate] = massFlux;
	}
	yRight[fMassFlowRate] = massFlux;
}

TUnstrPremFlamePhys::~TUnstrPremFlamePhys( void )
{
	int	nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();

//	delete fMassFraction;
	FreeNewtonInfo( fNewtonInfoL );
	DisposeVector( fYLeftVec );
	DisposeVector( fRhoY_iV_iPlus );

	fSolM->vec = &fSolM->vec[kPrev];
	fSavedM->vec = &fSavedM->vec[kPrev];

	DisposeVector( fSavedM );
	DisposeVector( fSolM );

	for ( int i = 0; i < nSpeciesInSystem+fVariablesWithoutSpecies; ++i ) {
		delete fVariableNames[i];
	}
	delete fVariableNames;
}

void UnstrPremPhysJacFirst( void *object, NodeInfoPtr nodeInfo )
{
	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fMassFlowRate = flame->GetOffsetMassFlowRate();
	int		nEq = nodeInfo->nOfEquations;
	int		nOfSpecies = flame->GetSpecies()->GetNOfSpecies();
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
	int		mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
    Double  hnenn = nodeInfo->hnenn;
	Double	**a = nodeInfo->a;
	Double	*y = nodeInfo->y;
	Double	*mixDensity = flameNode->mixDensity;
	Double	*mixHeatCapacity = flameNode->mixHeatCapacity;
	Double	oneOverCp = 1.0 / mixHeatCapacity[kCurr];
	Double	*enthalpy = flameNode->enthalpy;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*diffusivityPrev = flameNode->diffusivityPrev;
	Double	*productionRate = flameNode->productionRate;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	*mixConductivity = flameNode->mixConductivity;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];
	Double	diffMinusH;
	Double	dYmdY;	// dY(i)/dY(i-1) calculated from massflux boundary (left) condition

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

	if ( flame->GetSoot() ) {
		if ( !flame->fUseNumericalDM )  {
			flame->GetSoot()->UpdateJacobian( flame );
		}
		else {
			flame->FilldMomdMomOnePoint( flameNode );
		}
		flame->GetSoot()->FillJacobi( flame, nodeInfo, kPhysical );
	}

// first fill all convection terms
// first equation ( mass )
	FillJacFirstDerivDown( fMassFlowRate, fMassFlowRate, nodeInfo, kNegative );

#ifdef UPWINDCONVECTION
	Double	hCoeff = nodeInfo->h * ( nodeInfo->h + nodeInfo->hm );

// second to one + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fMassFlowRate, speciesEq, nodeInfo, 1.0 );
		if ( mixtureSpecificationLeft == kMassFlux && y[fMassFlowRate] > 0.0 ) {
			speciesIndexEq = speciesEq - fFirstSpecies;
			dYmdY = flame->GetdYPrevdY( speciesIndexEq, nodeInfo );
			a[speciesEq][speciesEq] -= hCoeff * y[fMassFlowRate] * dYmdY;
		}
	}
	if ( fTemperature < nEq ) {
		FillJacNonlinearConvectUpwind( fMassFlowRate, fTemperature, nodeInfo, 1.0 );
	}
#else
	// second to one + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		FillJacNonlinearConvectCentral( fMassFlowRate, speciesEq, nodeInfo );
		if ( mixtureSpecificationLeft == kMassFlux && y[fMassFlowRate] > 0.0 ) {
			speciesIndexEq = speciesEq - fFirstSpecies;
			dYmdY = flame->GetdYPrevdY( speciesIndexEq, nodeInfo );
			a[speciesEq][speciesEq] -= ( nodeInfo->h * nodeInfo->h ) * y[fMassFlowRate] * dYmdY;
		}
	}
	if ( fTemperature < nEq ) {
		FillJacNonlinearConvectCentral( fMassFlowRate, fTemperature, nodeInfo );
	}
#endif
	
// fill all diffusion and source terms
// first equation ( mass )
	
// second to one + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, speciesIndexEq, 1.0, nodeInfo, kNegative );
		if ( mixtureSpecificationLeft == kMassFlux ) {
			dYmdY = flame->GetdYPrevdY( speciesIndexEq, nodeInfo );
			diffMinusH = ( diffusivityPrev[speciesIndexEq] * mixDensity[kPrev]
						 + diffusivity[speciesIndexEq] * mixDensity[kCurr] ) * nodeInfo->h;
			a[speciesEq][speciesEq] -= diffMinusH * dYmdY;
		}

		if ( flame->UseDiffCorr() ) {
			flame->FillJacDiffCorr( speciesEq, 1.0, nodeInfo );
		}
		//	implicit source term
		a[fTemperature][speciesEq] -= dMdT[speciesIndexEq] * hnenn;
		for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq && speciesVar < nEq; ++ speciesVar ) {
			speciesIndexVar = speciesVar - fFirstSpecies;
			a[speciesVar][speciesEq] -= dMdY[speciesIndexVar][speciesIndexEq] * hnenn;
		}
	}

// last equation ( temperature )
	if ( fTemperature < nEq ) {
		FillJacWithDiffusion( fTemperature, fTemperature, oneOverCp, mixConductivity, nodeInfo, kNegative );
		for ( speciesIndexEq = 0; speciesIndexEq < nSpeciesInSystem; ++speciesIndexEq ) {
			a[fTemperature][fTemperature] += oneOverCp * ( dMdT[speciesIndexEq] * enthalpy[speciesIndexEq] 
														 + productionRate[speciesIndexEq] * heatCapacity[speciesIndexEq] ) * hnenn;
			for ( speciesIndexVar = 0; speciesIndexVar < nSpeciesInSystem; ++speciesIndexVar ) {
				a[speciesIndexVar + fFirstSpecies][fTemperature] += oneOverCp * dMdY[speciesIndexVar][speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
			}
		}
		if ( flame->fProperties->GetRadiation() ) {
			flame->fProperties->GetRadiation()->FillJacRadiation( -flame->fRadFact * oneOverCp, flame, nodeInfo );
			if ( flame->GetSoot() && flame->GetSoot()->WithSootRadiation() ) {
				flame->GetSoot()->FillJacSootRadiation( oneOverCp, flame, nodeInfo );
			}
		}
	}
	if ( flame->GetSolver()->bt->GetTimedepFlag() 
				&& !flame->GetSolver()->time->GetTimeConverged() ) {
		TTimePtr tim = flame->GetSolver()->time;
		for ( int eqLoop = 0; eqLoop < nEq; ++eqLoop ) {
			if ( eqLoop != fMassFlowRate ) {
				a[eqLoop][eqLoop] += mixDensity[kCurr] * hnenn / tim->GetDeltaT();
			}
		}
	}
}

Double TUnstrPremFlamePhys::GetdYPrevdY( int speciesIndex, NodeInfoPtr nodeInfo )
{
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	M = nodeInfo->yPrev[fMassFlowRate];
	Double	coeff = mixDensity[kPrev] * diffusivityPrev[speciesIndex]
								/ ( nodeInfo->hm * M );
								
	Double	coeffCorr = 0.0;
	if ( UseDiffCorr() ) {
		coeffCorr = mixDensity[kPrev] / M * fFlameNode->diffCorr[kPrev];
	}

	return coeff / ( 1.0 + coeff + coeffCorr );
}

void UnstrPremPhysJacRest( void *object, NodeInfoPtr nodeInfo )
{
	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fMassFlowRate = flame->GetOffsetMassFlowRate();
	int		nEq = nodeInfo->nOfEquations;
	int		nOfSpecies = flame->GetSpecies()->GetNOfSpecies();
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
    Double  hnenn = nodeInfo->hnenn;
	Double	**a = nodeInfo->a;
	Double	*mixDensity = flameNode->mixDensity;
	Double	*mixHeatCapacity = flameNode->mixHeatCapacity;
	Double	oneOverCp = 1.0 / mixHeatCapacity[kCurr];
	Double	*enthalpy = flameNode->enthalpy;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	*productionRate = flameNode->productionRate;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	*mixConductivity = flameNode->mixConductivity;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

	if ( flame->GetSoot() ) {
		if ( !flame->fUseNumericalDM )  {
			flame->GetSoot()->UpdateJacobian( flame );
		}
		else {
			flame->FilldMomdMomOnePoint( flameNode );
		}
		flame->GetSoot()->FillJacobi( flame, nodeInfo, kPhysical );
	}

// first fill all convection terms
// first equation ( mass )
	if ( nodeInfo->gridPoint == flame->fTfixLoc ) {
		if ( !flame->fConstMassFlux && fTemperature < nEq ) {
			a[fTemperature][fMassFlowRate] = 1.0;	// if the energy equation is solved, give a BC here
		}
		else {
			a[fMassFlowRate][fMassFlowRate] = 1.0;	// if no energy eqation is solved, don't change M
		}
	}
	else if ( nodeInfo->gridPoint < flame->fTfixLoc ) {
		FillJacFirstDerivDown( fMassFlowRate, fMassFlowRate, nodeInfo, kNegative );
	}
	else if ( nodeInfo->gridPoint > flame->fTfixLoc ) {
		FillJacFirstDerivUp( fMassFlowRate, fMassFlowRate, nodeInfo );
	}
	

#ifdef UPWINDCONVECTION
// second to one + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fMassFlowRate, speciesEq, nodeInfo, 1.0 );
	}
	if ( fTemperature < nEq ) {
		FillJacNonlinearConvectUpwind( fMassFlowRate, fTemperature, nodeInfo, 1.0 );
	}
#else
	// second to one + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		FillJacNonlinearConvectCentral( fMassFlowRate, speciesEq, nodeInfo );
	}
	if ( fTemperature < nEq ) {
		FillJacNonlinearConvectCentral( fMassFlowRate, fTemperature, nodeInfo );
	}
#endif
	
// fill all diffusion and source terms
// first equation ( mass )
	
// second to one + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, speciesIndexEq, 1.0, nodeInfo, kNegative );

		if ( flame->UseDiffCorr() ) {
			flame->FillJacDiffCorr( speciesEq, 1.0, nodeInfo );
		}
		//	implicit source term
		a[fTemperature][speciesEq] -= dMdT[speciesIndexEq] * hnenn;
		for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq && speciesVar < nEq; ++ speciesVar ) {
			speciesIndexVar = speciesVar - fFirstSpecies;
			a[speciesVar][speciesEq] -= dMdY[speciesIndexVar][speciesIndexEq] * hnenn;
		}
	}

// last equation ( temperature )
	if ( fTemperature < nEq ) {
		FillJacWithDiffusion( fTemperature, fTemperature, oneOverCp, mixConductivity, nodeInfo, kNegative );
		for ( speciesIndexEq = 0; speciesIndexEq < nSpeciesInSystem; ++speciesIndexEq ) {
			a[fTemperature][fTemperature] += oneOverCp * ( dMdT[speciesIndexEq] * enthalpy[speciesIndexEq] 
														 + productionRate[speciesIndexEq] * heatCapacity[speciesIndexEq] ) * hnenn;
			for ( speciesIndexVar = 0; speciesIndexVar < nSpeciesInSystem; ++speciesIndexVar ) {
				a[speciesIndexVar + fFirstSpecies][fTemperature] += oneOverCp * dMdY[speciesIndexVar][speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
			}
		}
		if ( flame->fProperties->GetRadiation() ) {
			flame->fProperties->GetRadiation()->FillJacRadiation( -flame->fRadFact * oneOverCp, flame, nodeInfo );
			if ( flame->GetSoot() && flame->GetSoot()->WithSootRadiation() ) {
				flame->GetSoot()->FillJacSootRadiation( oneOverCp, flame, nodeInfo );
			}
		}
	}
	if ( flame->GetSolver()->bt->GetTimedepFlag() 
				&& !flame->GetSolver()->time->GetTimeConverged() ) {
		TTimePtr tim = flame->GetSolver()->time;
		for ( int eqLoop = 0; eqLoop < nEq; ++eqLoop ) {
			if ( eqLoop != fMassFlowRate ) {
				a[eqLoop][eqLoop] += mixDensity[kCurr] * hnenn / tim->GetDeltaT();
			}
		}
	}
}

void UnstrPremPhysJacLast( void *object, NodeInfoPtr nodeInfo )
{
	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fMassFlowRate = flame->GetOffsetMassFlowRate();
	int		nEq = nodeInfo->nOfEquations;
	int		nOfSpecies = flame->GetSpecies()->GetNOfSpecies();
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
	int		mixtureSpecificationRight = flame->GetMixtureSpecificationRight();
    Double  hnenn = nodeInfo->hnenn;
	Double	**a = nodeInfo->a;
	Double	**c = nodeInfo->c;
	Double	*y = nodeInfo->y;
	Double	*mixDensity = flameNode->mixDensity;
	Double	*mixHeatCapacity = flameNode->mixHeatCapacity;
	Double	oneOverCp = 1.0 / mixHeatCapacity[kCurr];
	Double	*enthalpy = flameNode->enthalpy;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*diffusivityNext = flameNode->diffusivityNext;
	Double	*productionRate = flameNode->productionRate;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	*mixConductivity = flameNode->mixConductivity;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

	if ( flame->GetSoot() ) {
		if ( !flame->fUseNumericalDM )  {
			flame->GetSoot()->UpdateJacobian( flame );
		}
		else {
			flame->FilldMomdMomOnePoint( flameNode );
		}
		flame->GetSoot()->FillJacobi( flame, nodeInfo, kPhysical );

		// update right boundary condition
		int	sootoff = flame->GetSoot()->GetOffsetSootMoments();
		int	nSootMom = flame->GetSoot()->GetNSootMoments();
		Double	diff = flameNode->diffSoot[kCurr];
		Double	diffNext = flameNode->diffSoot[kNext];
		Double	fact = nodeInfo->hm * ( diff * mixDensity[kCurr]
												+ diffNext * mixDensity[kNext] );
		if ( flame->GetSoot()->WithSizeDepDiff() ) {
		}
		else {
			for ( int i = 0; i < nSootMom; ++i ) {
				a[i+sootoff][i+sootoff] -= fact;
			}
		}
	}

// first fill all convection terms
// first equation ( mass )
	FillJacFirstDerivUp( fMassFlowRate, fMassFlowRate, nodeInfo );

#ifdef UPWINDCONVECTION
// second to one + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fMassFlowRate, speciesEq, nodeInfo, 1.0 );
	}
	if ( fTemperature < nEq ) {
		FillJacNonlinearConvectUpwind( fMassFlowRate, fTemperature, nodeInfo, 1.0 );
	}
#else
	// second to one + nOfSpecies equation ( species )
	Double	h = nodeInfo->h;
	Double	hm = nodeInfo->hm;
	Double	hh = h * h;
	Double	hcoeff = ( h + hm ) * ( h + hm );
	Double	dfdfm = - hcoeff / ( hh - hcoeff );
	Double	dfdfmm = hh / ( hh - hcoeff );
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		FillJacNonlinearConvectCentral( fMassFlowRate, speciesEq, nodeInfo );
		a[speciesEq][speciesEq] += hm * hm * dfdfm * y[fMassFlowRate];
		c[speciesEq][speciesEq] -= hh * dfdfmm * y[fMassFlowRate];
	}
	if ( fTemperature < nEq ) {
		FillJacNonlinearConvectCentral( fMassFlowRate, fTemperature, nodeInfo );
		a[fTemperature][fTemperature] += hm * hm * dfdfm * y[fMassFlowRate];
		a[fTemperature][fTemperature] -= hh * dfdfmm * y[fMassFlowRate];
	}
#endif
	
// fill all diffusion and source terms
// first equation ( mass )
	
// second to one + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, speciesIndexEq, 1.0, nodeInfo, kNegative );
		a[speciesEq][speciesEq] -= ( diffusivity[speciesIndexEq] * mixDensity[kCurr] 
								   + diffusivityNext[speciesIndexEq] * mixDensity[kNext] ) * nodeInfo->hm;

		if ( flame->UseDiffCorr() ) {
			flame->FillJacDiffCorr( speciesEq, 1.0, nodeInfo );
		}
		//	implicit source term
		a[fTemperature][speciesEq] -= dMdT[speciesIndexEq] * hnenn;
		for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq && speciesVar < nEq; ++ speciesVar ) {
			speciesIndexVar = speciesVar - fFirstSpecies;
			a[speciesVar][speciesEq] -= dMdY[speciesIndexVar][speciesIndexEq] * hnenn;
		}
	}

// last equation ( temperature )
	if ( fTemperature < nEq ) {
		FillJacWithDiffusion( fTemperature, fTemperature, oneOverCp, mixConductivity, nodeInfo, kNegative );
		a[fTemperature][fTemperature] -= oneOverCp * ( mixConductivity[kCurr] + mixConductivity[kNext] ) *
										 nodeInfo->hm;
		for ( speciesIndexEq = 0; speciesIndexEq < nSpeciesInSystem; ++speciesIndexEq ) {
			a[fTemperature][fTemperature] += oneOverCp * ( dMdT[speciesIndexEq] * enthalpy[speciesIndexEq] 
														 + productionRate[speciesIndexEq] * heatCapacity[speciesIndexEq] ) * hnenn;
			for ( speciesIndexVar = 0; speciesIndexVar < nSpeciesInSystem; ++speciesIndexVar ) {
				a[speciesIndexVar + fFirstSpecies][fTemperature] += oneOverCp * dMdY[speciesIndexVar][speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
			}
		}
		if ( flame->fProperties->GetRadiation() ) {
			flame->fProperties->GetRadiation()->FillJacRadiation( -flame->fRadFact * oneOverCp, flame, nodeInfo );
			if ( flame->GetSoot() && flame->GetSoot()->WithSootRadiation() ) {
				flame->GetSoot()->FillJacSootRadiation( oneOverCp, flame, nodeInfo );
			}
		}
	}
	if ( flame->GetSolver()->bt->GetTimedepFlag() 
				&& !flame->GetSolver()->time->GetTimeConverged() ) {
		TTimePtr tim = flame->GetSolver()->time;
		for ( int eqLoop = 0; eqLoop < nEq; ++eqLoop ) {
			if ( eqLoop != fMassFlowRate ) {
				a[eqLoop][eqLoop] += mixDensity[kCurr] * hnenn / tim->GetDeltaT();
			}
		}
	}
}

void UnstrPremPhysRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	if ( !flame->RHSAction( nodeInfo, rhsMode ) ) {
		return;
	}
	
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fMassFlowRate = flame->GetOffsetMassFlowRate();
	int		eqLoop, speciesEq;
	int		nEq = nodeInfo->nOfEquations;
	int		nOfSpecies = flame->GetSpecies()->GetNOfSpecies();
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
	Double	*enthalpy = flameNode->enthalpy;
	Double	*mixDensity = flameNode->mixDensity;
	Double	*productionRate = flameNode->productionRate;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	oneOverCp = 1.0 / flameNode->mixHeatCapacity[kCurr];
	Double	sumCpDdYdx, diffCorrHeat = 0.0, sumCpYD;
	Double	sumMH;

	if ( flame->GetSoot() ) {
		flame->GetSoot()->FillRHS( flame, nodeInfo, kPhysical );
	}

// first fill all convection terms
	// first equation ( mass )
	if ( nodeInfo->gridPoint == flame->fTfixLoc ) {
		if ( !flame->fConstMassFlux && fTemperature < nEq ) {
			rhs[fMassFlowRate] += ( y[fTemperature] - flame->fTfix ) / hnenn;	// if the energy equation is solved, give a BC here
		}
		else {
			rhs[fMassFlowRate] += ( y[fMassFlowRate] - flame->fConstMassFlowRate ) / hnenn;		// if no energy equation is solved, don't change M
		}
	}
	else if ( nodeInfo->gridPoint < flame->fTfixLoc ) {
		rhs[fMassFlowRate] -= FirstDerivUpwind( yNext[fMassFlowRate], y[fMassFlowRate], h );
	}
	else if ( nodeInfo->gridPoint > flame->fTfixLoc ) {
		rhs[fMassFlowRate] += FirstDerivUpwind( y[fMassFlowRate], yPrev[fMassFlowRate], hm );
	}


#ifdef UPWINDCONVECTION
	// second to one + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		rhs[speciesEq] += NonlinearConvectUpwind( y[fMassFlowRate], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
	}
#	ifdef MIXEDUPWINDCENTRAL	
	if ( fTemperature < nEq ) {
		rhs[fTemperature] += 0.5 * NonlinearConvectUpwind( y[fMassFlowRate], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h )
							+ 0.5 * NonlinearConvectCentral( y[fMassFlowRate], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
	}
#	else
	if ( fTemperature < nEq ) {
		rhs[fTemperature] += NonlinearConvectUpwind( y[fMassFlowRate], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
	}
#	endif

#else
	// second to one + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		rhs[speciesEq] += NonlinearConvectCentral( y[fMassFlowRate], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
	}
	if ( fTemperature < nEq ) {
		rhs[fTemperature] += NonlinearConvectCentral( y[fMassFlowRate], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
	}
#endif

// mass equation
	
// fFirstSpecies to fFirstSpecies + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		eqLoop = speciesEq - fFirstSpecies;
		rhs[speciesEq] -= productionRate[eqLoop];
	}
#ifdef FULLDIFFUSION
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < nEq; ++speciesEq ) {
		eqLoop = speciesEq - fFirstSpecies;
		rhs[speciesEq] += flame->SecondDerivPolySpecDiff( speciesEq, nodeInfo );
	}
#else
#endif
	flame->CompleteSpeciesDiffusion( nodeInfo );
	
// energy equation
	if ( fTemperature < nEq ) {
		sumCpDdYdx = 0.0;
		sumMH = 0.0;
		sumCpYD = 0.0;
		diffCorrHeat = ( flame->UseDiffCorr() ) ? flameNode->mixHeatCapacity[kCurr] : 0.0;

// diffusion
		rhs[fTemperature] -= oneOverCp * SecondDerivDiffusion( fTemperature, flameNode->mixConductivity, nodeInfo );
		
// compute all sums
		for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
			sumCpDdYdx += ( heatCapacity[eqLoop] - diffCorrHeat ) * diffusivity[eqLoop] 
						* FirstDeriv( YPrev[eqLoop], Y[eqLoop], YNext[eqLoop], hm, h );
#	ifdef MOLARDIFFUSION
			sumCpYD += ( heatCapacity[eqLoop] - diffCorrHeat ) * Y[eqLoop] * diffusivity[eqLoop];
#	endif
		}

		for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
			sumMH += productionRate[eqLoop] * enthalpy[eqLoop];
		}

//	enthalpy flux
#ifdef FULLDIFFUSION
		rhs[fTemperature] += flame->HeatFluxPolySpecDiff( nodeInfo )
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#else
		rhs[fTemperature] -= oneOverCp * mixDensity[kCurr] * sumCpDdYdx 
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#		ifdef MOLARDIFFUSION
		rhs[fTemperature] -= oneOverCp * sumCpYD * mixDensity[kCurr]
						/ flameNode->mixMolarMass[kCurr]
						* FirstDeriv( flameNode->mixMolarMass[kPrev], flameNode->mixMolarMass[kCurr], flameNode->mixMolarMass[kNext], hm, h )
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#		endif
#endif
		
		rhs[fTemperature] += oneOverCp * sumMH;
		
		if ( flame->fProperties->GetRadiation() ) {
			rhs[fTemperature] -= flame->fRadFact * oneOverCp * flameNode->radiation[kCurr];
			if ( flame->GetSoot() && flame->GetSoot()->WithSootRadiation() ) {
				rhs[fTemperature] += oneOverCp * flame->GetSoot()->GetSootRadiation( y[fTemperature], flameNode->moments );
			}
		}
		
	}

	TTimePtr tim = flame->GetSolver()->time;
	for ( eqLoop = 0; eqLoop < nEq; ++eqLoop ){
		if ( flame->GetSolver()->bt->GetTimedepFlag() && eqLoop != fMassFlowRate
				&& !flame->GetSolver()->time->GetTimeConverged() ) {
			rhs[eqLoop] += mixDensity[kCurr] * ( y[eqLoop] - tim->GetYOld()->mat[nodeInfo->gridPoint][eqLoop] ) 
								/ tim->GetDeltaT();
		}
		rhs[eqLoop] *= - hnenn;
	}
}

void TUnstrPremFlamePhys::FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
{
// fills the jacobian with     constCoeff * sum_j ( d/dy(rho Y_k D_j / sum(Y_i) * dY_j/dy) )

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
	Double	coeffCurr = constCoeff * density[kCurr] * Y[lVar];
	Double	coeffPrev = constCoeff * density[kPrev] * YPrev[lVar];
	Double	coeffNext = constCoeff * density[kNext] * YNext[lVar];
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
		if ( nodeInfo->firstPoint ) {
			if ( GetMixtureSpecificationLeft() == kMassFlux ) {
				a[lVar][nVariable] += diffMinusH * GetdYPrevdY( i, nodeInfo );
			}
		}
		else {
			c[lVar][nVariable] += diffMinusH;
		}
		if ( nodeInfo->lastPoint ) {
			a[lVar][nVariable] += diffPlusHm;
		}
		else {
			b[lVar][nVariable] += diffPlusHm;
		}
	}
}

Double TUnstrPremFlamePhys::DiffCorr( int nVariable, NodeInfoPtr nodeInfo )
{
// returns     sum_j ( d/dy(rho Y_k D_j / sum(Y_i) dY_j/dy) )

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
	Double	coeffCurr = density[kCurr] * Y[i];
	Double	coeffPrev = density[kPrev] * YPrev[i];
	Double	coeffNext = density[kNext] * YNext[i];
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

Double TUnstrPremFlamePhys::DiffCorrX( int nVariable, NodeInfoPtr nodeInfo )
{
// returns     sum_j ( d/dy(rho Y_k / M * D_j Y_j dM/dy) )

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
	Double	coeffCurr = density[kCurr] * Y[k] / M;
	Double	coeffPrev = density[kPrev] * YPrev[k] / MPrev;
	Double	coeffNext = density[kNext] * YNext[k] / MNext;
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
#ifdef DIFFCORRCORR
	coeffCurr /= sumYCurr;
	coeffPrev /= sumYPrev;
	coeffNext /= sumYNext;
#endif	
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		diffPlus = coeffCurr * diffusivity[i] * Y[i] + coeffNext * diffusivityNext[i] * YNext[i];
		diffMinus = coeffCurr * diffusivity[i] * Y[i] + coeffPrev * diffusivityPrev[i] * YPrev[i];
		value += ( diffPlus * hm * ( MNext - M ) 
					+ diffMinus * h * ( MPrev - M ) );
	}

	return value / nodeInfo->hnenn;
}

int UnstrPremPhysPostIter( void *object )
{
	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	int			fMassFlowRate = flame->GetOffsetMassFlowRate();
	int			fTemperature = flame->GetOffsetTemperature();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();
	int			i;
	TNewtonPtr	bt = flame->GetSolver()->bt;
	NodeInfoPtr nodeInfo = bt->GetNodeInfo();
    TAdaptiveGridPtr	grid = bt->GetGrid();
	TGridPtr 	currGrid = grid->GetCurrentGrid();
	int			nGridPoints = currGrid->GetNGridPoints();
	int			*bcFlagLeft = currGrid->GetBcFlagLeft();
	int			*bcFlagRight = currGrid->GetBcFlagRight();
	int			mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
	int			mixtureSpecificationRight = flame->GetMixtureSpecificationRight();
	MatrixPtr	yMat = currGrid->GetY();
	VectorPtr	yLeftVec = currGrid->GetYLeft();
	VectorPtr	yRightVec = currGrid->GetYRight();
	Double		**y = yMat->mat;
	Double		*x = currGrid->GetX()->vec;
	Double		*yLeft = yLeftVec->vec;
	Double		*yRight = yRightVec->vec;
	Double		hLast = bt->GetRight() - x[nGridPoints-1];
	Double		*bcLeft = currGrid->GetBcLeft()->vec;
	Double		*bcRight = currGrid->GetBcRight()->vec;
	Double		hFirst = x[0] - bt->GetLeft();
	T1DSpeciesPtr	species = flame->GetSpecies();
	T1DPropertiesPtr	properties = flame->GetProperties();
	Double		*diffusivityLeft = species->GetDiffusivity()->mat[-1];
	Double		*diffusivityRight = species->GetDiffusivity()->mat[nGridPoints];
	Double		lambdaLeft = properties->GetConductivity()->vec[-1];
	Double		lambdaRight = properties->GetConductivity()->vec[nGridPoints];
	Double		cpLeft = properties->GetHeatCapacity()->vec[-1];
	Double		cpRight = properties->GetHeatCapacity()->vec[nGridPoints];
	Double		&rhoLeft = properties->GetDensity()->vec[-1];
	Double		&rhoRight = properties->GetDensity()->vec[nGridPoints];
	Double		*temp = flame->GetTemperature()->vec;
	Double		**Y = flame->GetMassFracs()->mat;
	int			variables = bt->GetNVariables();
	int			nSpecies = species->GetNOfSpecies();
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();
	static int      iDidit = 0;

	for ( i = 0; i < nGridPoints; ++i ) {
		if ( flame->CheckSolution( y[i][fTemperature], &y[i][fFirstSpecies], nSpeciesInSystem ) ) {
			return 1;
		}
	}

//	update properties
	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, temp[kPrev], Y[kPrev], flame->GetPressure() );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, temp[nGridPoints], Y[nGridPoints], flame->GetPressure() );

	UnstrPremPhysUpdateLeftBoundary( flame );
	UnstrPremPhysUpdateRightBoundary( flame );
	
	flame->UpdateSolutionOnePoint( yLeft, kPrev );
	flame->UpdateSolutionOnePoint( yRight, nGridPoints );
	
	flame->UpdateThermoProps();
	
	if ( flame->fConstMassFlux ) {
		Double	mixMolarMass = 0.0;
		flame->fProperties->ComputeMixtureMolarMass( 
				mixMolarMass, &bcLeft[fFirstSpecies], flame->GetSpecies()->GetMolarMass()->vec, nSpeciesInSystem );
		Double	rho298 = flame->GetPressure() * mixMolarMass / ( RGAS * 298.0 );
		flame->fConstMassFlowRate = flame->GetStrainRate() * rho298;
		flame->SetMassFlux( flame->fConstMassFlowRate, currGrid->GetY()
				, currGrid->GetYLeft()->vec, currGrid->GetYLeft()->vec );
	}

	// set point of fixed temperature
	Double Tmin = temp[LocationOfMin( nGridPoints+2, &temp[kPrev] ) - 1];//printf("tmin = %g\t",Tmin);
	Double Tmax = temp[LocationOfMax( nGridPoints+2, &temp[kPrev] ) - 1];//printf("tmax = %g\t",Tmax);
	Double deltaTFix = (Tmax - Tmin)/2;
	if (deltaTFix > flame->fDeltaTFix) deltaTFix = flame->fDeltaTFix;
		for ( i = 0; i < nGridPoints; ++i ) {
		  if ( temp[i] > ( yLeft[fTemperature] +  deltaTFix /*flame->fDeltaTFix*/ ) ) {
				flame->fTfixLoc = i;
				break;
			}
		}

	if ( grid->IsFine() ) flame->fTfixLoc += flame->fTfixLoc % 2 - 1;	// make it odd for the fine grid ( even gridpoints are removed for the coarse grid! )
	flame->fTfix = temp[flame->fTfixLoc];

	return 0;
}

#include "TofZ.h"

void TUnstrPremFlamePhys::UpdateDimensions( int len )
{
	T1DFlame::UpdateDimensions( len );
	fSolM->len = len;
}

void TUnstrPremFlamePhys::UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
{
	int		nGridPoints = yMat->cols;
	Double	*M = fSolM->vec;
	Double	**y = yMat->mat;
	Double	*yLeft = yLeftVec->vec;
	Double	*yRight = yRightVec->vec;

	UpdateDimensions( nGridPoints );

	T1DFlame::UpdateSolution( yMat, yLeftVec, yRightVec );
	
	M[kPrev] = yLeft[fMassFlowRate];
	for ( int k = 0; k < nGridPoints; ++k ) {
		M[k] = y[k][fMassFlowRate];
	}
	M[nGridPoints] = yRight[fMassFlowRate];
}

void TUnstrPremFlamePhys::UpdateSolutionOnePoint( Double *y, int gridPoint )
{
	T1DFlame::UpdateSolution( y, gridPoint );

	fSolM->vec[gridPoint] = y[fMassFlowRate];
}

int	TUnstrPremFlamePhys::GetOffsetMassFlowRate( void )
{
	return fMassFlowRate; 
}

int	TUnstrPremFlamePhys::GetOffsetTemperature( void )
{
	return fTemperature; 
}

int	TUnstrPremFlamePhys::GetOffsetFirstSpecies( void ) 
{
	return fFirstSpecies;
}

ConstStringArray TUnstrPremFlamePhys::GetVariableNames( void )
{
	return ( ConstStringArray )fVariableNames;
}

int TUnstrPremFlamePhys::GetVariablesWithoutSpecies( void )
{
	return fVariablesWithoutSpecies;
}

void TUnstrPremFlamePhys::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp )
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
	int					nSpecies = fSpecies->GetNOfSpecies();
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
	Flag				ySet = FALSE;
	Flag				oxidizerFound = FALSE;
	Flag				chooseInputGrid = FALSE;
	Double				*xIn = new Double[gridPointsIn];
	if ( !xIn ) FatalError( "memory allocation of TUnstrPremFlamePhys failed" );
	Double				*yIn =  new Double[gridPointsIn];
	if ( !yIn ) FatalError( "memory allocation of TUnstrPremFlamePhys failed" );
	Double				*yInFloat = sp->data;
	Double				**y = yMat->mat;
	int					variable;
	char				*string = sp->labels;
	SplinePtr			theSpline = NULL;
	enum InputSoot		{ kNoIn = 0, kMassFracs, kMoleFracs };
	InputSoot			sootin = kNoIn;
	Double				leftSlope;
	Double				rightSlope;
	Double				*temp = GetTemperature()->vec;
	Double				**Y = GetMassFracs()->mat;
	FILE				*fp;
	TGridPtr			fine = adapGrid->GetFine();
	TGridPtr			coarse = adapGrid->GetCoarse();
	
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
		if ( strncmp( string, "y", 1 ) == 0 ) {
			if ( chooseInputGrid ) {
				cerr << "choose inputGrid" << NEWL;
				for ( j = 0; j < gridPointsIn-2; ++j ) {
					locX[j] = yInFloat[i*gridPointsIn + j+1];		// implicit cast from float to Double
				}
				bt->SetLeft( yInFloat[i*gridPointsIn] );
				bt->SetRight( yInFloat[(i+1)*gridPointsIn-1] );
				left = bt->GetLeft();
				right = bt->GetRight();
				ySet = TRUE;
			}
			else { // choose own Grid, but read x for interpolation
				cerr << "choose own Grid" << NEWL;
				for ( j = 0; j < gridPointsIn; ++j ) {
					xIn[j] = yInFloat[i*gridPointsIn + j];		// implicit cast from float to Double
				}
				if ( right - left == 0 ) {
					left = yInFloat[i*gridPointsIn];
					right = yInFloat[(i+1)*gridPointsIn-1];
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
				ySet = TRUE;
			}
		}
		else if ( strcmp( string, "massfraction-o2" ) == 0 ) {
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
	if ( !ySet ) {
		cerr << "error: can't find coordinate 'y'" << NEWL;
		exit(2);
	}
	if ( !oxidizerFound ) {
		cerr << "error: can't find massfraction of oxidizer 'massfraction-o2'" << NEWL;
		exit(2);
	}
	
// reset string
	string = sp->labels;
	
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strncmp( string, "massflowrate", 12 ) == 0 ) {
			variable = fMassFlowRate;
		}
		else if ( strncmp( string, "temperature", 11 ) == 0 ) {
			variable = fTemperature;
		}
		else if ( strncmp( string, "conc-soot", 9 ) == 0 && GetSoot() && sootin != kMassFracs ) {
			string += 9;
			char name[3];
			sootin = kMoleFracs;
			for ( j = 0; j < GetSoot()->GetNSootMoments(); ++j ) {
				sprintf( name, "%d", j );
				if ( strncmp( name, string, 1 ) == 0 ) {
					variable = j + GetSoot()->GetOffsetSootMoments();
					string += 1;
				}
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
				if ( strncmp( string, "SOOT", 4 ) == 0 && GetSoot() ) {
					string += 4;
					char *name = "0";
					for ( j = 0; j < GetSoot()->GetNSootMoments(); ++j ) {
						sprintf( name, "%d", j );
						if ( strncmp( name, string, 1 ) == 0 ) {
							sootin = kMassFracs;
							variable = j + GetSoot()->GetOffsetSootMoments();
							string += 1;
						}
					}
				}
				else {
					cerr << "warning: no match for species " << string << NEWL;
					string += strlen(string) + 1;
					continue;
				}
			}
		}
		else {
			string += strlen(string) + 1;
			continue;
		}

		string += strlen(string) + 1;
		if ( chooseInputGrid ) {
			for ( k = 0; k < gridPointsIn-2; ++k ) {
				y[k][variable] = yInFloat[i*gridPointsIn + k+1];	// copy workspace to vector of solution
			}
		}
		else {
			for ( k = 0; k < gridPointsIn; ++k ) {	// store vector in workspace
				yIn[k] = yInFloat[i * gridPointsIn + k];	// implicit cast from float to Double
			}
		
			leftSlope = ( yIn[1] - yIn[0] ) / ( xIn[1] - xIn[0] );
			rightSlope = ( yIn[gridPointsIn-1] - yIn[gridPointsIn-2] ) / ( xIn[gridPointsIn-1] - xIn[gridPointsIn-2] );
			theSpline = ComputeSimpleSpline( xIn, yIn, gridPointsIn, FALSE, leftSlope, FALSE, rightSlope, NULL, TRUE );
			SplineInterpolate( theSpline, locX, yWork, nGridPoints );
			for ( k = 0; k < nGridPoints; ++k ) {
				y[k][variable] = yWork[k];	// copy workspace to vector of solution
			}
		}
		if ( variable == fTemperature ) {
			if ( yLeft[fTemperature] <= 0.0 ) {
				grid->GetBcLeft()->vec[fTemperature] = 
					yLeft[fTemperature] = yInFloat[i*gridPointsIn];
			}
			if ( yRight[fTemperature] <= 0.0 ) {
				grid->GetBcRight()->vec[fTemperature] = 
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

// set initial Boundary values
	for ( j = 0; j < variables; ++j ) {
		yRight[j] = y[nGridPoints-1][j];
	}
	yLeft[fMassFlowRate] = y[0][fMassFlowRate];
	
// set initial Boundary values and pressure

	if ( GetPressure() <= 0.0 ) {
		struct _parameter	*param = GetParameter( "pressure" );
		Double	thePressure;
		if ( param ) {
			thePressure = (Double)param->what.quantity.value;
			if ( strcmp( param->what.quantity.unit, "bar" ) == 0 ) {
				thePressure *= 1.0e5;
			}
			SetPressure( thePressure );
			fprintf( stderr, "%s%g%s\n", "initial pressure is ", GetPressure()/1.0e5, " bar"  );
		}
		else { // exit
			cerr << "#error: no value for 'pressure' in inputfile" << NEWL;
			exit(2);
		}
	}

	
	if ( GetSoot() ) {
		if ( sootin == kMassFracs ) {
			int		sootOff = GetSoot()->GetOffsetSootMoments();
			Double	MSoot;
			for ( k = 0; k < nGridPoints; ++k ) {
				MSoot = 24;
				for ( j = 0; j < GetSoot()->GetNSootMoments(); ++j ) {
					y[k][j + GetSoot()->GetOffsetSootMoments()] /= MSoot;
				}
			}
		}
		else if ( sootin == kMoleFracs ) {
			Double molarMass;
			Double rho;
			for ( k = 0; k < nGridPoints; ++k ) {
				fProperties->ComputeMixtureMolarMass( molarMass, &y[k][fFirstSpecies]
								, fSpecies->GetMolarMass()->vec, nSpeciesInSystem );
				rho = GetPressure() * molarMass / ( RGAS * y[k][fTemperature] );
				for ( j = 0; j < GetSoot()->GetNSootMoments(); ++j ) {
					y[k][j + GetSoot()->GetOffsetSootMoments()] /= rho;
				}
			}
		}
	}
	if ( GetInputData()->fParameterComm >= 0.0 ) {
		VectorPtr	phi = GetPhiVector();
		
		DisposeVector( phi );
		phi = NewVector( 1 );
		phi->vec[0] = GetInputData()->fParameterComm;
		SetPhiVec( phi );
		InitPhi();
		SetMassFracsOfPhi( this, GetPhi(), &fine->GetBcLeft()->vec[fFirstSpecies]
			, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		SetMassFracsOfPhi( this, GetPhi(), &coarse->GetBcLeft()->vec[fFirstSpecies]
			, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		SetMassFracsOfPhi( this, GetPhi(), &fine->GetYLeft()->vec[fFirstSpecies]
			, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
	}
	
	if ( GetPhi() <= 0.0 ) {
		struct _parameter	*param = GetParameter( "fuel-air-equivalence-ratio" );
		Double	equivalenceRatio;
		if ( param ) {
			equivalenceRatio = (Double)param->what.quantity.value;
			SetPhi( equivalenceRatio );
			SetMassFracsOfPhi( this, GetPhi(), &fine->GetBcLeft()->vec[fFirstSpecies]
				, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
			SetMassFracsOfPhi( this, GetPhi(), &coarse->GetBcLeft()->vec[fFirstSpecies]
				, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
			SetMassFracsOfPhi( this, GetPhi(), &fine->GetYLeft()->vec[fFirstSpecies]
				, nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		}
		else { // exit
			cerr << "#error: no value for 'fuel-air-equivalence-ratio' in inputfile" << NEWL;
			exit(2);
		}
	}
//	update properties
	UpdateSolution( yMat, yLeftVec, yRightVec );
	if ( GetSpecies()->IsConstantLewisNumber() ) {
		CompLewisNumbers( GetSpecies()->GetLewisNumberFile() );
	}

	UnstrPremPhysPostIter( this );
	
	if ( fConstMassFlux && fStrainRate && fStrainRate->vec[0] > 0.0 ) { // start value for massflux speciefiled
		Double	mixMolarMass = 0.0;
		fProperties->ComputeMixtureMolarMass( 
				mixMolarMass, &fine->GetBcLeft()->vec[fFirstSpecies]
				, fSpecies->GetMolarMass()->vec, nSpeciesInSystem );
		Double	rho298 = GetPressure() * mixMolarMass / ( RGAS * 298.0 );
		fConstMassFlowRate = GetStrainRate() * rho298;
		SetMassFlux( fConstMassFlowRate, fine->GetY()
				, fine->GetYLeft()->vec, fine->GetYLeft()->vec );
		fprintf( stderr, "mixMolarMass(left) = %g\trho298 = %g\tv298 = %g m/s\tmassflux = %g\n"
						, mixMolarMass, rho298, GetStrainRate(), fConstMassFlowRate );
	}
	else if ( fConstMassFlux ) {
		fprintf( stderr, "error: no mass flow rate specified\n" );
		exit( 2 );
	}

	FreeSpline( theSpline );
	delete yIn;
	delete xIn;

	adapGrid->SetSolutionScaler();

	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fYLeftVec->vec[i] = yLeft[fFirstSpecies+i];
	}
	fp = GetOutfile( "initialguess", TFlame::kData );
	bt->PrintSolution( locX, y, GetVariableNames(), fp );
	fclose(fp);
}

Double TUnstrPremFlamePhys::SpeciesDiffusion( int nVariable, int speciesIndex, NodeInfoPtr nodeInfo )
{
// returns    d/dy( rho D_k dY_k/dy) )

	Double	yPrev = nodeInfo->yPrev[nVariable];
	Double	y = nodeInfo->y[nVariable];
	Double	yNext = nodeInfo->yNext[nVariable];
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	diffPlusHm = nodeInfo->hm * ( diffusivity[speciesIndex] * mixDensity[kCurr]
					+ diffusivityNext[speciesIndex] * mixDensity[kNext] );
	Double	diffMinusH = nodeInfo->h * ( diffusivityPrev[speciesIndex] * mixDensity[kPrev]
					+ diffusivity[speciesIndex] * mixDensity[kCurr] );
	
	return ( ( diffPlusHm * ( yNext - y ) + diffMinusH * ( yPrev - y ) ) 
				/ nodeInfo->hnenn );
}

Double TUnstrPremFlamePhys::SecondDerivXDiffusion( int nVariable, NodeInfoPtr nodeInfo )
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
	
	Double	diffPlusX = diffusivity[speciesIndex] * mixDensity[kCurr] * y / M 
					+ diffusivityNext[speciesIndex] * mixDensity[kNext] * yNext / MNext;
	Double	diffMinusX = diffusivityPrev[speciesIndex] * mixDensity[kPrev] * yPrev / MPrev
					+ diffusivity[speciesIndex] * mixDensity[kCurr] * y / M;

	return ( diffPlusX * hm * ( MNext - M ) + diffMinusX * h * ( MPrev - M ) ) 
				/ nodeInfo->hnenn;
}

void TUnstrPremFlamePhys::FillJacSpeciesDiffusion( int nVariable, int speciesIndex, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
{
// fills the jacobian with     constCoeff * d/dy ( rho * diffusivity * df/dy)

	if ( sign == kNegative ) {
		constCoeff *= -1.0;
	}

	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	diffPlusHm = constCoeff * ( diffusivity[speciesIndex] * mixDensity[kCurr]
					+ diffusivityNext[speciesIndex] * mixDensity[kNext] ) * nodeInfo->hm;
	Double	diffMinusH = constCoeff * ( diffusivityPrev[speciesIndex] * mixDensity[kPrev]
					+ diffusivity[speciesIndex] * mixDensity[kCurr] ) * nodeInfo->h;

	nodeInfo->a[nVariable][nVariable] -= ( diffPlusHm + diffMinusH );
	if ( !nodeInfo->lastPoint ) {
		nodeInfo->b[nVariable][nVariable] += diffPlusHm;
	}
	if ( !nodeInfo->firstPoint ) {
		nodeInfo->c[nVariable][nVariable] += diffMinusH;
	}
}

void UnstrPremPhysOutput( void *object, FILE *fp, char* tail )
{
	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	TNewtonPtr		bt = flame->GetSolver()->bt;
	T1DPropertiesPtr	props = flame->GetProperties();
	TSpeciesPtr		species = flame->GetSpecies();
	NodeInfoPtr		nodeInfo = bt->GetNodeInfo();
	Double			*rho = props->GetDensity()->vec;
	Double			*mixMolarMass = props->GetMolarMass()->vec;
	Double			*molarMass = species->GetMolarMass()->vec;
	TGridPtr		currentGrid = bt->GetGrid()->GetCurrentGrid();
	Double			*x = currentGrid->GetX()->vec;
	Double			**massFracs = flame->GetMassFracs()->mat;
	Double			*temp = flame->GetTemperature()->vec;
	Double			*M = flame->GetM()->vec;
	Double			sL = flame->GetBurningVelocity();
	Double			**y = currentGrid->GetY()->mat;
	Double			*yLeft = currentGrid->GetYLeft()->vec,
					*yRight = currentGrid->GetYRight()->vec;
	int				i, k;
	int				gridPoints = currentGrid->GetNGridPoints();
	int				nOfSpecies = species->GetNOfSpecies();
	int				nOfVariables = bt->GetNVariables();
	int				nOfEquations = bt->GetNEquations();
	int				firstSpecies = flame->GetOffsetFirstSpecies();
	int				tempOffset = flame->GetOffsetTemperature();
	int				fMassFlowRate = flame->GetOffsetMassFlowRate();
	time_t			theDate;
	char			buffer[80];
	ConstStringArray	varNames = flame->GetVariableNames();
	char			**names = species->GetNames();
	Flag			fpOpen = FALSE;
	
	if ( !fp ) {
		fpOpen = TRUE;
		fp = flame->GetOutputFile( NULL, tail, TFlame::kNone );
	}
	
// write header
	fprintf( fp, "header\n\n" );

	fprintf( fp, "title = \"unstretched freely propagating premixed flame\"\n" );
	fprintf( fp, "mechanism = \"%s\"\n", flame->GetInputData()->fReactionFile );
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
	fprintf( fp, "fuel-air-equivalence-ratio = %g\n", flame->GetPhi() );
	if ( flame->fConstMassFlux ) {
		fprintf( fp, "v0 = %g [m/s]\n", flame->GetStrainRate() );
	}

	if ( flame->GetSoot() ) {
		fprintf( fp, "Nucleation = \"%s\"\n", ( flame->GetSoot()->WithNucleation() ) ? "TRUE" : "FALSE" );
		fprintf( fp, "Condensation = \"%s\"\n", ( flame->GetSoot()->WithCondensation() ) ? "TRUE" : "FALSE" );
		fprintf( fp, "Coagulation = \"%s\"\n", ( flame->GetSoot()->WithCoagulation() ) ? "TRUE" : "FALSE" );
		fprintf( fp, "SurfaceGrowth = \"%s\"\n", ( flame->GetSoot()->WithSurfaceGrowth() ) ? "TRUE" : "FALSE" );
		fprintf( fp, "SurfaceOxidation = \"%s\"\n", ( flame->GetSoot()->WithSurfaceOxidation() ) ? "TRUE" : "FALSE" );
		fprintf( fp, "Thermophoresis = \"%s\"\n", ( flame->GetSoot()->WithThermoPhoresis() ) ? "TRUE" : "FALSE" );
		fprintf( fp, "SootRadiation = \"%s\"\n", ( flame->GetSoot()->WithSootRadiation() ) ? "TRUE" : "FALSE" );
		fprintf( fp, "SootUpdateProdRate = \"%s\"\n", ( flame->GetSoot()->WithUpdateProdRates() ) ? "TRUE" : "FALSE" );
		fprintf( fp, "SizeDepDiffusion = \"%s\"\n", ( flame->GetSoot()->WithSizeDepDiff() ) ? "TRUE" : "FALSE" );
	}
	
	
	if ( species->IsConstantLewisNumber() ) {
		fprintf( fp, "ConstantLewisNumbers = \"True\"\n" );
	}
	
	fprintf( fp, "Tmax = %g [K]\n\n", temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
	int	nTemp10 = 0;
	while ( x[nTemp10] < 10.0e-3 && nTemp10 < gridPoints ) ++nTemp10;
	fprintf( fp, "T10mm = %g [K]\n\n", temp[nTemp10] );
	if ( x[nTemp10] < 10.0e-3 ) {
		fprintf( fp, "xOfT10mm = %g [mm]\n\n", x[nTemp10]*1.0e-3 );
	} 
	
	Double	moleFrac;
	fprintf( fp, "unburnt\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", yLeft[tempOffset] );
	for ( i = 0; i < nOfSpecies; ++i ) {
		if ( fabs( massFracs[0][i] ) > 1.0e-4 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], massFracs[0][i] );
		}
	}
	for ( i = 0; i < nOfSpecies; ++i ) { // write X_i
		moleFrac = massFracs[kPrev][i] * mixMolarMass[kPrev] / molarMass[i];
		if ( fabs( moleFrac ) > 1.0e-5 ) {
			fprintf( fp, "\tMolefraction-%s = %g\n", names[i], moleFrac );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "burningVelocity = %g [cm/sec]\n", sL );
	fprintf( fp, "FlameThickness = %g [m]\n"
		, flame->ComputeFlameThickness( temp, x, gridPoints ) );
	
	fprintf( fp, "numOfSpecies = %d\n", nOfSpecies );
	fprintf( fp, "gridPoints = %d\n\n", gridPoints+2 );

	fprintf( fp, "body\n" );

// write independent coordinate
	fprintf( fp, "y [m]\n" );
	fprintf( fp, "\t%-.6e", bt->GetLeft() );
	for ( k = 0; k < gridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", x[k] );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", bt->GetRight() );
			
	if ( flame->fConstMassFlux ) {
	// write independent coordinate in [mm]
		fprintf( fp, "h [mm]\n" );
		fprintf( fp, "\t%-.6e", bt->GetLeft()*1000.0 );
		for ( k = 0; k < gridPoints; ++k ) {
			fprintf( fp, "\t%-.6e", x[k]*1000.0 );
			if ( (k+2) % 5 == 0 ) {
				fprintf( fp, "\n" );
			}
		}
		fprintf( fp, "\t%-.6e\n", bt->GetRight()*1000.0 );
				
		Double currTime = 0.0;
		fprintf( fp, "t [ms]\n" );
		fprintf( fp, "\t%-.6e", 0.0 );
		currTime += 0.5 * ( rho[kPrev]/M[kPrev] + rho[kCurr]/M[kCurr] ) * ( x[kCurr] - bt->GetLeft() );
		fprintf( fp, "\t%-.6e", currTime );
		for ( k = 1; k < gridPoints; ++k ) {
			currTime += 0.5 * ( rho[k-1]/M[k-1] + rho[k]/M[k] ) * ( x[k] - x[k-1] );
			fprintf( fp, "\t%-.6e", currTime*1000.0 );
			if ( (k+2) % 5 == 0 ) {
				fprintf( fp, "\n" );
			}
		}
		currTime += 0.5 * ( rho[gridPoints-1]/M[gridPoints-1] + rho[gridPoints]/M[gridPoints] ) 
								* ( bt->GetRight() - x[gridPoints-1] );
		fprintf( fp, "\t%-.6e\n", currTime*1000.0 );
	}

// write solution
	// write massflowrate and temperature
	flame->PrintFlameletVector( gridPoints+2, &M[kPrev], "massflowrate [kg/m^2s]", fp );
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
		// print mole fractions
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

		// print concentrations
	}
	
	if ( flame->fSoot ) {
		flame->GetSoot()->PrintFlameletFile( gridPoints, flame, fp );
	}

//	write enthalpy
	fprintf( fp, "TotalEnthalpy [J/kmole]\n" );
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

//	write density
	fprintf( fp, "density\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", rho[k-1] );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( (k+1) % 5 ) {
		fprintf( fp, "\n" );
	}

//	write eta
	Double *eta = New1DArray(gridPoints+2);
	eta = &eta[1];
	Double etacoeff = sqrt( 100.0/props->GetViscosity()->vec[kPrev]*rho[kPrev] );

	eta[-1] = 0.0;
	eta[0] = etacoeff * rho[0] * ( x[0] - bt->GetLeft() );
	for ( k = 2; k < gridPoints+1; ++k ) {
		eta[k-1] = eta[k-2] + etacoeff * rho[k-1] * ( x[k-1] - x[k-2] );
	}
	eta[gridPoints] = eta[gridPoints-1] + etacoeff * rho[gridPoints] * ( bt->GetRight() - x[gridPoints-1] );

	fprintf( fp, "eta\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", eta[k-1] );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( (k+1) % 5 ) {
		fprintf( fp, "\n" );
	}

//	write f
	Double *f = New1DArray(gridPoints+2);
	etacoeff = 1.0 / sqrt( 100.0 * props->GetViscosity()->vec[kPrev] * rho[kPrev] );
	f = &f[1];
	f[-1] = -M[kPrev] * etacoeff;
	f[gridPoints] = 0.0;
	for ( k = 0; k < gridPoints; ++k ) {
		f[k] = f[-1] + ( f[gridPoints] - f[-1] ) / ( eta[gridPoints] - eta[-1] ) * ( eta[k] - eta[-1] );
//		f[k] = ( f[gridPoints] - f[-1] ) / ( bt->GetRight() - bt->GetLeft() ) * ( x[k] - bt->GetLeft() );
	}
	fprintf( fp, "f\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", f[k-1] );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( (k+1) % 5 ) {
		fprintf( fp, "\n" );
	}
	
// write df/deta	
	Double *dfdeta = New1DArray(gridPoints+2);
	dfdeta = &dfdeta[1];
	dfdeta[-1] = 1.0;
	for ( k = 0; k < gridPoints; ++k ) {
		dfdeta[k] = ( f[k+1] - f[-1] ) / ( eta[k+1] - eta[-1] );
	}
	dfdeta[gridPoints] = dfdeta[gridPoints-1];
	fprintf( fp, "df/deta\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", dfdeta[k-1] );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( (k+1) % 5 ) {
		fprintf( fp, "\n" );
	}
	
	dfdeta = &dfdeta[-1];
	eta = &eta[-1];
	f = &f[-1];
	Free1DArray(dfdeta);
	Free1DArray(f);
	Free1DArray(eta);

//	write heat capacity
	flame->PrintFlameletVector( gridPoints+2, &props->GetHeatCapacity()->vec[kPrev], "cp [J/m^3K]", fp );

//	write heat conductivity
	flame->PrintFlameletVector( gridPoints+2, &props->GetConductivity()->vec[kPrev], "lambda [W/m]", fp );

//	write viscosity
	flame->PrintFlameletVector( gridPoints+2, &props->GetViscosity()->vec[kPrev], "mu [kg/s^2m]", fp );
	
	fprintf( fp, "trailer\n" );
	
	if ( flame->fRadFact != 1.0 ) {
		fprintf( fp, "RadFact = %g\n", flame->fRadFact );
	}

	if ( flame->fSoot ) {
		if ( flame->fSoot->GetCoagFact() < 1.0 ) {
			fprintf( fp, "CoagFact =  %g\n", flame->GetSoot()->GetCoagFact() );
		}
	}

	if ( species->IsConstantLewisNumber() ) {
		Double	*Le = species->GetLewisNumber()->vec;
		for ( i = 0; i < nOfSpecies; ++i ) {
			fprintf( fp, "%s\t%g\n", names[i], Le[i] );
		}
	}
	if ( nOfEquations < nOfVariables) {
		fprintf( fp, "number of converged equations is %d\n", nOfEquations );
	}
	
	if ( fpOpen ) {
		fclose( fp );
	}
}

void SetUnstrPremPhysNodeInfo( int k, void *object )
{
	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	
	flame->SetFlameNode( k );
}

void UnstrPremPhysPostConv( void *object )
{
	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	TNewtonPtr			bt = flame->GetSolver()->bt;
	TAdaptiveGridPtr	grid = bt->GetGrid();
	TGridPtr			fine = grid->GetFine();
	int					fFirstSpecies = flame->GetOffsetFirstSpecies();
	int					nSpeciesInSystem = flame->fSpecies->GetNSpeciesInSystem();
	int					isConverged = bt->GetConvergeNewton();
	Flag				leaveContinPrem;
	Double				*yleftFine = fine->GetYLeft()->vec;
	Double				*bcLeftFine = fine->GetBcLeft()->vec;
	
	// check temperature gradient at right boundary
	if ( isConverged ) {
		flame->SaveSolution();
		bt->WriteOutput( object, NULL, "" );
		fprintf( stderr, "burning velocity is %g cm/s\n", flame->GetBurningVelocity() );
		if ( flame->AdjustComputationalDomain() ) {
			if ( flame->CheckComputationalDomain() ) {
				flame->fSolver->ReInit();
				return;
			}
		}
		if ( flame->fPrintRHSSpecies ) {
			flame->PrintRHSSpecies( flame->GetSolver()->bt );
		}
        flame->PostConvergence( object );
		UnstrPremPhysPostIter( flame );
	} else {
		flame->RestoreSolution();
        flame->fSolver->ReInit();
		flame->PostConvergence( object );
	}
	flame->ResetNGridModifications();
	
	if ( bt->GetLeaveContin() ) {
		if ( isConverged ) {
			for ( int i = 0; i < flame->fNSensObj; ++i ) {
				flame->GetSpecies()->PrintProdRateTerms( flame->fSensObj[i], flame );
			}
			bt->WriteOutput( object, NULL, "" );
			if ( flame->fReactionFluxes ) {
				flame->GetReaction()->PrintReactionRates( flame );
				flame->fReaction->PrintRateCoeffs( flame );
				flame->ReactionFluxes( kPhysical );
			}
			if ( flame->fSensAnal ) {
				flame->SensitivityAnalysis( 1.0, 1.0, kPhysical );
			}
		}
		leaveContinPrem = flame->PostConvTPremixed( isConverged );
		if ( !leaveContinPrem ) {
			flame->SetMassFracsOfPhi( flame, flame->GetPhi()
					, &fine->GetBcLeft()->vec[fFirstSpecies], nSpeciesInSystem
					, flame->fSpecies->GetMolarMass()->vec, flame->fSpecies->GetNames() );
			flame->fSolver->ReInit();
		}
	}
}

ConstStringArray GetUnstrPremPhysVarNames( void *object )
{
	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	
	return flame->GetVariableNames();
}

FILE *TUnstrPremFlamePhys::GetOutputFile( char *head, char *tail, FileType type )
{
	int				fuelIndex = GetFuelIndex();
	char			*name = new char[64];
	if ( !name ) FatalError( "new failed" );
	FILE			*fp;
	char			**speciesNames = fSpecies->GetNames();
	int				tu = ( int ) fSolTemp->vec[kPrev];
	Double			phi = GetPhi();
	Double			press = GetPressure() * 1.0e-5;
		
	if ( fConstMassFlux ) {
		sprintf( name, "%s%s%.8s_p%.2d_%.1dphi%.1d_%.4dtu%.4dv%.4d%s"
						, ( head ) ? head : "", ( head ) ? "_" : ""
						, speciesNames[fuelIndex]
						, ( int ) floor( press )	// in [bar]
						, ( int ) ( ( press - ( floor( press ) ) ) * 10 + 0.5 )
						, ( int ) floor( phi ) , ( int ) ( ( phi - floor(phi) ) * 10000 + 0.5 )
						, ( int )( tu )							// in [K]
						, ( int ) floor( GetStrainRate() * 10000 + 0.5 )	// in [bar]
						, ( tail ) ? tail : "" );
	}
	else {
		sprintf( name, "%s%s%.8s_p%.2d_%.1dphi%.1d_%.4dtu%.4d%s"
						, ( head ) ? head : "", ( head ) ? "_" : ""
						, speciesNames[fuelIndex]
						, ( int ) floor( press )	// in [bar]
						, ( int ) ( ( press - ( floor( press ) ) ) * 10 + 0.5 )
						, ( int ) floor( phi ) , ( int ) ( ( phi - floor(phi) ) * 10000 + 0.5 )
						, ( int )( tu )							// in [K]
						, ( tail ) ? tail : "" );
	}

	fp = GetOutfile( name, type );
	delete name;

	return fp;
}

void TUnstrPremFlamePhys::SetInitialBC( TGridPtr grid, TInputDataPtr inp )
{
	int					i;
	Double				mixMolarMass;
	int					nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	SpeciesPtr			species = inp->GetSpecies();
	BoundaryInputPtr	right = inp->rightBoundary;
	BoundaryInputPtr	left = inp->leftBoundary;
	int					inpTOffset = inp->fTemperatureOffset;
	int					*speciesIndexLeft = NULL;
	int					*speciesIndexRight = NULL;
	int					leftSpecifiedSpecies = left->fSpecifiedSpeciesBCs;
	int					rightSpecifiedSpecies = right->fSpecifiedSpeciesBCs;
	int					*bcFlagLeft = grid->GetBcFlagLeft();
	int					*bcFlagRight = grid->GetBcFlagRight();
	Double				*yleft = grid->GetYLeft()->vec;
	Double				*yright = grid->GetYRight()->vec;
	Double				*bcLeft = grid->GetBcLeft()->vec;
	Double				*bcRight = grid->GetBcRight()->vec;
	Double				*phi = GetPhiVector()->vec;
	
	speciesIndexLeft = new int[left->fSpecifiedSpeciesBCs];
	if ( !speciesIndexLeft ) FatalError( "memory allocation of TUnstrPremFlamePhys failed" );
	speciesIndexRight = new int[right->fSpecifiedSpeciesBCs];
	if ( !speciesIndexRight ) FatalError( "memory allocation of TUnstrPremFlamePhys failed" );

	//	set speciesIndex
	for ( i = 0; i < leftSpecifiedSpecies; ++i ) {
		speciesIndexLeft[i] = inp->FindSpecies( left->speciesName[i] );
	}
	for ( i = 0; i < rightSpecifiedSpecies; ++i ) {
		speciesIndexRight[i] = inp->FindSpecies( right->speciesName[i] );
	}
	
	// set fMixtureSpecification
	SetMixtureSpecificationLeft( left->fMixtureSpecification );
	SetMixtureSpecificationRight( right->fMixtureSpecification );
	
	// set BCFlags
	bcFlagLeft[fTemperature] = left->fBcFlag[inpTOffset];
	bcFlagRight[fTemperature] = right->fBcFlag[inpTOffset];
	for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
		bcFlagLeft[i] = left->fBcFlagSpecies;
		bcFlagRight[i] = right->fBcFlagSpecies;
	}

	// set value
	yleft[fTemperature] = left->fValue[inpTOffset];
	yright[fTemperature] = right->fValue[inpTOffset];

	bcLeft[fTemperature] = left->fValue[inpTOffset];
	bcRight[fTemperature] = right->fValue[inpTOffset];

	for ( i = 0; i < leftSpecifiedSpecies; ++i ) {
		yleft[speciesIndexLeft[i]+fFirstSpecies] = left->fValueSpecies[i];
		bcLeft[speciesIndexLeft[i]+fFirstSpecies] = left->fValueSpecies[i];
	}
	if ( left->fMixtureSpecification == kMolarFraction ) {
		// first compute molar mass of mixture
		for ( i = 0, mixMolarMass = 0; i < nSpeciesInSystem; ++i ) {
			mixMolarMass += species[i].molarMass * yleft[i+fFirstSpecies];
		}
		// compute massfractions
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			yleft[i+fFirstSpecies] *= species[i].molarMass / mixMolarMass;
			bcLeft[i+fFirstSpecies] = yleft[i+fFirstSpecies];
		}
		for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
			bcFlagLeft[i] = kMassFraction;
		}
	}

	for ( i = 0; i < rightSpecifiedSpecies; ++i ) {
		yright[speciesIndexRight[i]+fFirstSpecies] = right->fValueSpecies[i];
		bcRight[speciesIndexRight[i]+fFirstSpecies] = right->fValueSpecies[i];
	}
	if ( right->fMixtureSpecification == kMolarFraction ) {
		// first compute molar mass of mixture
		for ( i = 0, mixMolarMass = 0; i < nSpeciesInSystem; ++i ) {
			mixMolarMass += species[i].molarMass * yright[i+fFirstSpecies];
		}
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			yright[i+fFirstSpecies] *= species[i].molarMass / mixMolarMass;
			bcRight[i+fFirstSpecies] = yright[i+fFirstSpecies];
		}
		for ( i = fFirstSpecies; i < nSpeciesInSystem+fFirstSpecies; ++i ) {
			bcFlagRight[i] = kMassFraction;
		}
	}

	if ( phi[0] > 0 ) {
		SetMassFracsOfPhi( this, phi[0], &bcLeft[fFirstSpecies], nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		SetMassFracsOfPhi( this, phi[0], &yleft[fFirstSpecies], nSpeciesInSystem, fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
		SetMixtureSpecificationLeft( kMassFlux );
	}
	else if ( leftSpecifiedSpecies > 0 ) {
		SetPhiOfMassFracs( this, phi, &bcLeft[fFirstSpecies], fSpecies->GetMolarMass()->vec, fSpecies->GetNames() );
	}
	else {
		cerr << "read boundary conditions for the unburnt from flamelet file" << NEWL;
		SetMixtureSpecificationLeft( kMassFlux );
	}

	delete speciesIndexRight;
	delete speciesIndexLeft;
}

Double TUnstrPremFlamePhys::GetBurningVelocity( void )
{
	// [cm/s]
	return GetM()->vec[kCurr] / GetProperties()->GetDensity()->vec[0] * 100.0;
}

int	TUnstrPremFlamePhys::GetOffsetVVelocity( void )
{
	return fMassFlowRate; 
}

int	TUnstrPremFlamePhys::GetOffsetUVelocity( void )
{
	cerr << "#error: class has no member fUVelocity" << NEWL;
	exit( 2 );
	return 0; 
}

int	TUnstrPremFlamePhys::GetOffsetMixFrac( void )
{
	cerr << "#error: class has no member fMixtureFraction" << NEWL;
	exit( 2 );
	return 0; 
}

void TUnstrPremFlamePhys::PrintRHSSpecies( TNewtonPtr bt )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
    int         		N = currentGrid->GetNGridPoints();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		i, k;
	char				**names = GetSpecies()->GetNames();
	int					nSpecIn = fSpecies->GetNSpeciesInSystem();
	int					num = 7;
	int					start = 0, end;
	int					fileCount = 1;
	FILE				*fp;
	
	UpdateThermoProps();
	
	do {
		end = ( int )MIN( start + 100 / num, nSpecIn );
		sprintf( GetOutFileBuff(), "%sspeciesRHS%d.dout", GetOutputPath(), fileCount++ );
		if ( !( fp = fopen( GetOutFileBuff(), "w") ) ) { 
			cerr << "#warning: unable to open file " << GetOutFileBuff() << NEWL;
			exit(2);
		}
#if defined (applec) || defined (powerc)
  	  RotateCursor( 32 * bt->GetNIter() );
#endif
		
		fprintf( fp, "*\n%-12s", "y" );
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
#if defined (applec) || defined (powerc)
		RotateCursor( 32 );
#endif
			bt->SetNodeInfo( this, k );
			PrintRHSSpecies( start, end, nodeInfo, fp );
		}
		fclose( fp );
		start = end;
	} while ( end < nSpecIn );
}

void TUnstrPremFlamePhys::PrintRHSSpecies( int start, int end, NodeInfoPtr nodeInfo,  FILE *fp )
{
	int		eqLoop, j;
	int		speciesEq;
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	diffTerm;
	Double	*reactionRate = fFlameNode->reactionRate;
	Double	productionRate;
	Double	source;
	Double	sink;
	int		*nOfUsedReactions = fSpecies->GetNOfUsedReactions()->vec;
	VectorPtr	*nu = fSpecies->GetNu();
	IntVectorPtr	*usedReactions = fSpecies->GetUsedReactions();
	Double	*molarMass =  fSpecies->GetMolarMass()->vec;

	fprintf( fp, "%-.6e", *nodeInfo->x );

	for ( eqLoop = start; eqLoop < end; ++eqLoop ) {
		speciesEq = fFirstSpecies+eqLoop;

#ifdef UPWINDCONVECTION
		fprintf( fp, "\t%-.6e", NonlinearConvectUpwind( y[fMassFlowRate], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h ) );
#else
		fprintf( fp, "\t%-.6e", NonlinearConvectCentral( y[fMassFlowRate], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h ) );
#endif

		diffTerm = -SpeciesDiffusion( speciesEq, eqLoop, nodeInfo );
		fprintf( fp, "\t%-.6e", diffTerm );
		if ( UseDiffCorr() ) {
			diffTerm = DiffCorr( speciesEq, nodeInfo );
			fprintf( fp, "\t%-.6e", diffTerm );
		}
		if ( fThermoDiffusion ) {
			diffTerm = ThermoDiffusion( eqLoop, kPhysical, nodeInfo );
			fprintf( fp, "\t%-.6e", diffTerm );
		}

		sink = source = 0.0;
		for ( j = 0; j < nOfUsedReactions[eqLoop]; ++j ) {
			productionRate = nu[eqLoop]->vec[j] * reactionRate[usedReactions[eqLoop]->vec[j]];
			if ( productionRate > 0.0 ) {
				sink += productionRate;
			}
			else {
				source += productionRate;
			}
		}
		sink *= - molarMass[eqLoop];
		source *= - molarMass[eqLoop];
		fprintf( fp, "\t%-.6e\t%-.6e\t%-.6e", source, -sink, source+sink );
	}
	fprintf( fp, "\n" );
}

void UnstrPremPhysUpdateLeftBoundary( void  *object )
{
	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	int				fMassFlowRate = flame->GetOffsetMassFlowRate();
	int				speciesOff = flame->GetOffsetFirstSpecies();
	TNewtonPtr		bt = flame->GetSolver()->bt;
	TGridPtr 		currGrid = bt->GetGrid()->GetCurrentGrid();
	MatrixPtr		yMat = currGrid->GetY();
	Double			**y = yMat->mat;
	Double			*yLeft = currGrid->GetYLeft()->vec;
	Double			*bcLeft = currGrid->GetBcLeft()->vec;
	int				nSpeciesInSystem = flame->fSpecies->GetNSpeciesInSystem();
	Double			pressure = flame->GetPressure();
	int				mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
	Double			mixDensity = flame->fProperties->GetDensity()->vec[kPrev];
	Double			*diffusivity = flame->fSpecies->GetDiffusivity()->mat[kPrev];
	Double			hm = currGrid->GetX()->vec[0] - bt->GetLeft();

	yLeft[fMassFlowRate] = y[0][fMassFlowRate];
	flame->UpdateSolutionOnePoint( y[0], 0 );
	if ( mixtureSpecificationLeft == kMassFlux && y[0][fMassFlowRate] > 0.0 ) {
		flame->CalcYLeft();
	}
	else {
		for ( int i = 0; i < nSpeciesInSystem; ++i ) {
			yLeft[speciesOff+i] = bcLeft[speciesOff+i];
		}
	}

	if ( flame->GetSoot() ) {
		int	nSootMoments = flame->GetSoot()->GetNSootMoments();
		int	sootOff = flame->GetSoot()->GetOffsetSootMoments();
		for ( int i = 0; i < nSootMoments; ++i ) {
			switch ( i ) {
				case 0:
					yLeft[i+sootOff] = 1.0e-20;
					break;
				default:
					yLeft[i+sootOff] = 9 * yLeft[i+sootOff-1];
					break;
			}
		}
	}

	flame->UpdateSolutionOnePoint( yLeft, kPrev );
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], pressure );
}

void UnstrPremPhysUpdateRightBoundary( void *object )
{
	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	TNewtonPtr		bt = flame->GetSolver()->bt;
	TGridPtr 		currGrid = bt->GetGrid()->GetCurrentGrid();
	int				nGridPoints = currGrid->GetNGridPoints();
	MatrixPtr		yMat = currGrid->GetY();
	Double			**y = yMat->mat;
	Double			*yLast = y[nGridPoints-1];
	Double			*yRight = currGrid->GetYRight()->vec;
	Double			pressure = flame->GetPressure();
	int				variables = bt->GetNVariables();
	Double			*x = currGrid->GetX()->vec;
	Double			h = bt->GetRight() - x[nGridPoints-1];

#ifdef UPWINDCONVECTION
	for ( int j = 0; j < variables; ++j ) {
		yRight[j] = yLast[j];
	}
#else
	Double			hm = x[nGridPoints-1] - x[nGridPoints-2];
	Double			hh = h * h;
	Double			hcoeff = ( h + hm ) * ( h + hm );
	Double			*yLastLast = y[nGridPoints-2];
	for ( int j = 0; j < variables; ++j ) {
		yRight[j] = ( hh * yLastLast[j] - hcoeff * yLast[j] ) / ( hh - hcoeff );
	}
#endif

	flame->UpdateSolutionOnePoint( yRight, nGridPoints );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], pressure );
}

int TUnstrPremFlamePhys::CheckComputationalDomain( void )
{
	int				last;
	TBVPSolverPtr	solver = GetSolver();
	TNewtonPtr		bt = solver->bt;
	TGridPtr 		grid = bt->GetGrid()->GetFine();
	VectorPtr		xVec = grid->GetX();
	Double			*x = xVec->vec;
	Double			xLeft = bt->GetLeft();
	int				nGridPoints = grid->GetNGridPoints();
	Double			*temp = fSolTemp->vec;
		
	// check left boundary
	const Double	critDeltaT = 50.0;
	Double			critLength;
	int				nThickness_opt = 4;	
	int				nThickness_min = 3;	
	int				nThickness_max = 5;	// nThickness_min * nThickness <= length of the zone >= flameThickness * nThickness_max
	Double			flameThickness = ComputeFlameThickness( temp, x, nGridPoints );
	int				crit = 0;
	int				nPointsToMove = 0;
	const int		nPointsToMove_max = 10;
	
	while ( fabs( temp[++crit] - temp[kPrev] ) < critDeltaT ) ;
	critLength = x[crit] - xLeft;
	if ( critLength < nThickness_min * flameThickness ) {
		// MapMan uses last dx for extension
		nPointsToMove = ( int ) floor( ( nThickness_opt * flameThickness - critLength ) / ( x[1] - x[0] ) ) + 1;
		if ( nPointsToMove <= 0 ) fprintf( stderr, "#warning: adjustment of the preheat zone failed, i have problems moving the grid from the right to the left\n" );
		nPointsToMove *= -1;	// from the right to the left
	}
	else if ( critLength > nThickness_max * flameThickness ) {
		nPointsToMove = crit;
		while ( x[crit] - x[--nPointsToMove] + x[kCurr] - xLeft < nThickness_opt * flameThickness );
		if ( nPointsToMove <= 0 ) fprintf( stderr, "#warning: adjustment of the preheat zone failed, i have problems moving the grid from the left to the right\n" );
	}
	if ( nPointsToMove ) {
		if ( !bt->GetGrid()->GetOneSolOneGrid() ) {
			cerr << "#warning: cannot adjust the preheat zone if regridding not allowed" << NEWL;
			return 0;
		}
		if ( abs( nPointsToMove ) > nPointsToMove_max ) nPointsToMove = nPointsToMove_max * SIGN( nPointsToMove );
		NodeMover::fromType	fromTo = NodeMover::fromLeftSide;
		if ( nPointsToMove < 0 ) {
			nPointsToMove = -nPointsToMove;
			fromTo = NodeMover::fromRightSide;
		}
		fprintf( stderr, "%s%d%s%s%s\n", "move grid by ", nPointsToMove, " points to the " 
			, ( ( fromTo == NodeMover::fromLeftSide ) ? "right" : "left" ) 
			, NEWL  );

		int			i, k;
		int			nOfVars = bt->GetNVariables();
		Double		**y = grid->GetY()->mat;
		Double		*yLeft = grid->GetYLeft()->vec;
		Double		*yRight = grid->GetYRight()->vec;
		MMDataBag	bag( nOfVars );
		
		//	init
		bag.Initialize();
		bag.SetOldInpedVar( xVec, "x" );
		for ( i = 0; i < nOfVars; ++i ) {
			bag.Insert( &y[0][i], nGridPoints, nOfVars, GetVariableNames()[i] );
		}

		//	set new grid
		VectorPtr	newXVec = NewVector( nGridPoints );
		NodeMover nm( xVec, newXVec, nPointsToMove, fromTo );
		nm.MoveIt();
		bag.SetNewInpedVar( newXVec, "xNew" );

		//	map
		VectorPtr	newYVec = NewVector( nGridPoints );
		Double		*newY = newYVec->vec;
		for ( i = 0; i < bag.NumElems(); ++i ) {
			bag[i].Map( newYVec );
			for ( k = 0; k < nGridPoints; ++k ) {
				y[k][i] = newY[k];
			}
		}
		
		//	shift new grid
		Double	shift;
		Double	right;
		Double	*xNew = newXVec->vec;
		if ( fromTo == NodeMover::fromLeftSide ) {
			shift = x[nPointsToMove] - x[0];
			right = bt->GetRight() + nPointsToMove * ( x[nGridPoints-1] - x[nGridPoints-2] ) - shift;
		}
		else {
			shift = - nPointsToMove * ( x[1] - x[0] );
			right = x[nGridPoints-nPointsToMove] - shift;
		}
		for ( k = 0; k < nGridPoints; ++k ) {
			x[k] = xNew[k] - shift;
		}
		bt->SetRight( right );
		
		UnstrPremPhysPostIter( this );
		
		//	clean up
		DisposeVector( newYVec );
		DisposeVector( newXVec );
		return 1;
	}
	
	// check right boundary
	int				flameLocation = LocationOfMaxSlope( temp, x, nGridPoints );
	Double			phi = GetPhi();
	Double			alpha;
	const Double	alpha_min = 0.02;				// minimum enlargement
	const Double	alpha_max = 1.0;				// maximum enlargement
	Double			left =  bt->GetLeft();
	Double			right =  bt->GetRight();
	Double			new_right;
	
	critLength = right - x[flameLocation];
// settings by mbo
	nThickness_opt = 50;	
	nThickness_min = 46;	
	nThickness_max = 54;	// nThickness_min * nThickness <= length of the zone >= flameThickness * nThickness_max
// settings by hp to just see what happens for methane flames
//	nThickness_opt = 100;	
//	nThickness_min = 92;	
//	nThickness_max = 108;	// nThickness_min * nThickness <= length of the zone >= flameThickness * nThickness_max
// settings by hp to check chlorohydrocarbon flames
//	nThickness_opt = 12;	
//	nThickness_min = 16;	
//	nThickness_max = 20;	// nThickness_min * nThickness <= length of the zone >= flameThickness * nThickness_max
// settings by hp to check rich flames
/*	nThickness_opt = 350;	*/
/*	nThickness_min = 320;	*/
/*	nThickness_max = 380;	// nThickness_min * nThickness <= length of the zone >= flameThickness * nThickness_max*/
// settings by hp to check rich flames
/*	nThickness_opt = 150;	*/
/*	nThickness_min = 120;	*/
/*	nThickness_max = 180;	// nThickness_min * nThickness <= length of the zone >= flameThickness * nThickness_max*/

fprintf( stderr, "flameThickness is %g\n", flameThickness );
fprintf( stderr, "critLength is %g\n", critLength );
	if ( critLength < nThickness_min * flameThickness ) {
		if ( !bt->GetGrid()->GetOneSolOneGrid() ) {
			cerr << "#warning: cannot enlarge computational domain if regridding not allowed" << NEWL;
			return 0;
		}
		else {
			new_right = x[flameLocation] + nThickness_opt * flameThickness;
			alpha = ( new_right - right ) / right;
			if ( fabs( alpha ) < alpha_min ) {
				alpha = alpha_min * SIGN( alpha );
				new_right = ( 1.0 + alpha ) * right;
			}
			else if ( fabs( alpha ) > alpha_max ) {
				alpha = alpha_max * SIGN( alpha );
				new_right = ( 1.0 + alpha ) * right;
			}
			
			bt->SetRight( new_right );
			cerr << "enlarge the computational domain by " << alpha * 100.0 << " percent" << NEWL << NEWL;
			if ( fNEnlargeGrid >= fNCutGrid + 2 ) {
				fprintf( stderr, "#warning: number of grid modifications exceeded\n" );
				return 0;
			}
			else {
				++fNEnlargeGrid;
				return 1;
			}
		}
	}
	else if ( critLength > nThickness_max * flameThickness ) {
		if ( !bt->GetGrid()->GetOneSolOneGrid() ) {
			cerr << "#warning: cannot cut computational domain if regridding not allowed" << NEWL;
			return 0;
		}
		else {
			// look for length close to optimum and cut the computational domain
			for ( last = nGridPoints - 2; last > 0; --last ) {
				if ( x[last] - x[flameLocation] < nThickness_opt * flameThickness )
				break;
			}
			alpha = ( right - x[last+1] ) / right;
			right = x[last+1];
			bt->SetRight( right );
			grid->AdjustNGridPoints( last+1 );
			solver->UpdateAllDimensions( last+1 );
			UnstrPremPhysPostIter( this );
			cerr << "cut the computational domain by " << alpha * 100.0 << " percent" << NEWL << NEWL;
			if ( fNCutGrid >= 2 ) {
				fprintf( stderr, "#warning: number of grid modifications exceeded\n" );
				return 0;
			}
			else {
				++fNCutGrid;
				return 1;
			}
		}
	}
	return 0;
}

void TUnstrPremFlamePhys::SaveSolution( void )
{
	int		k;
	int		len = fSolM->len;
	Double	*M = fSolM->vec;
	Double	*saveM = fSavedM->vec;

	T1DFlame::SaveSolution();
	fSavedM->len = fSolM->len;

	for ( k = -1; k <= len; ++k ) {
		saveM[k] = M[k];
	}
}

void TUnstrPremFlamePhys::RestoreSolution( void )
{
	int		k;
	int		len = fSavedM->len;
	Double	*M = fSolM->vec;
	Double	*saveM = fSavedM->vec;

	UpdateDimensions( len );

	T1DFlame::RestoreSolution();

	for ( k = -1; k <= len; ++k ) {
		M[k] = saveM[k];
	}
	
	SolutionToSolver();
}

void TUnstrPremFlamePhys::SolutionToSolver( void )
{
	TNewtonPtr	bt = fSolver->bt;
	TGridPtr	grid = bt->GetGrid()->GetFine();
	int		nGridPoints = fSolM->len;
	Double	*M = fSolM->vec;
	Double	**y = grid->GetY()->mat;

	T1DFlame::SolutionToSolver();
	
	for ( int k = 0; k < nGridPoints; ++k ) {
		y[k][fMassFlowRate] = M[k];
	}
	
	UnstrPremPhysPostIter( this );
}

void TUnstrPremFlamePhys::CalcYLeft( void )
{
	int			i;
	int			nSpecIn = fSpecies->GetNSpeciesInSystem();
	Double		*YLeftNew;
	Double		*yLeft = fSolver->bt->GetGrid()->GetCurrentGrid()->GetYLeft()->vec;

	TGridPtr 		currGrid = fSolver->bt->GetGrid()->GetCurrentGrid();
	MatrixPtr		yMat = currGrid->GetY();

	SetFlameNode( 0 );
	ComputeProperties( fFlameNode, fFlameNode->temp[kCurr]
							, fFlameNode->Y[kCurr], GetPressure() );

	YLeftNew = GetYLeft( &yLeft[fFirstSpecies] );

	for ( i = 0; i < nSpecIn; ++i ) {
		yLeft[ fFirstSpecies + i ] = YLeftNew[i];
	}

	SetFlameNode( 0 );
	ComputeProperties( fFlameNode, fFlameNode->temp[kCurr]
							, fFlameNode->Y[kCurr], GetPressure() );
}

Double *TUnstrPremFlamePhys::GetYLeft( Double *Yguess )
{
	int i;
	int	speciesIn = fSpecies->GetNSpeciesInSystem();
	
	for ( i = 0; i < speciesIn; ++i ) {
		fYLeftVec->vec[i] = MAX(1e-8, Yguess[i] );
	}
	Double 		coeff;
	TGridPtr 	currGrid = fSolver->bt->GetGrid()->GetCurrentGrid();
	Double		*bcLeft = currGrid->GetBcLeft()->vec;
	Double		*x = currGrid->GetX()->vec;
	Double		hFirst = x[0] - fSolver->bt->GetLeft();
	Double		*yFirst = currGrid->GetY()->mat[0];

	NewtonSolve( fNewtonInfoL, FALSE, FALSE, this );

	if ( !fNewtonInfoL->converged ) {
		for ( i = 0; i < speciesIn; ++i ) {
			coeff = fProperties->GetDensity()->vec[-1] * fSpecies->GetDiffusivity()->mat[-1][i] 
					/ ( currGrid->GetYLeft()->vec[fMassFlowRate] * hFirst );
			fYLeftVec->vec[i] = ( bcLeft[fFirstSpecies+i] + coeff * yFirst[fFirstSpecies+i] ) 
				/ ( 1.0 + coeff );
			fprintf( stderr, "not converged\n" );
		}
		NewtonSolve( fNewtonInfoL, FALSE, FALSE, this );
		if ( !fNewtonInfoL->converged ) {
			fprintf( stderr, "newton finally not converged, use simple mixing\n" );
			for ( i = 0; i < speciesIn; ++i ) {
				coeff = fProperties->GetDensity()->vec[-1] * fSpecies->GetDiffusivity()->mat[-1][i] 
						/ ( currGrid->GetYLeft()->vec[fMassFlowRate] * hFirst );
				fYLeftVec->vec[i] = ( bcLeft[fFirstSpecies+i] + coeff * yFirst[fFirstSpecies+i] ) 
					/ ( 1.0 + coeff );
				fprintf( stderr, "YL[%s] = %g\tY0[%s] = %g\n"
					, fSpecies->GetNames()[i], fYLeftVec->vec[i], fSpecies->GetNames()[i], yFirst[fFirstSpecies+i]);
			}
		}
	}

	SetFlameNode( kPrev );

	ComputeProperties( fFlameNode, fSolver->bt->GetGrid()->GetCurrentGrid()->GetYLeft()->vec[fTemperature]
							, fYLeftVec->vec, GetPressure() );

	return fYLeftVec->vec;
}

int BCLeftNewtonFuncsTUnstr( const VectorPtr /*x*/, VectorPtr fVec, void *object )
{
// called by NewtonSolve

// F = rho v ( Y_i - Y_i^Liquid ) + rho Y_i V_i 

	TUnstrPremFlamePhysPtr	flame = ( TUnstrPremFlamePhysPtr )object;
	TGridPtr 	currGrid = flame->GetSolver()->bt->GetGrid()->GetCurrentGrid();
	Double		*x = currGrid->GetX()->vec;
	Double		hFirst = x[0] - flame->GetSolver()->bt->GetLeft();
	int			speciesIn = flame->GetSpecies()->GetNSpeciesInSystem();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();

	Double		*yLeft = flame->GetYLeftVec()->vec;
	Double		*bcLeft = currGrid->GetBcLeft()->vec;
	Double		*f = fVec->vec;
	
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, currGrid->GetYLeft()->vec[flame->GetOffsetTemperature()]
							, yLeft, flame->GetPressure() );

#ifdef FULLDIFFUSION
		fprintf( stderr, "not implemented\n" ); exit(2);
#else
	Double	*rhoY_iV_iPlus = flame->fRhoY_iV_iPlus->vec;
	flame->CalcAllDiffVeloRhoYiNext( yLeft, hFirst );
#endif
	for ( int i = 0; i < speciesIn; ++i ) {
		f[i] = ( yLeft[i] - MAX( bcLeft[fFirstSpecies+i], 1.0e-15 ) )
			+ rhoY_iV_iPlus[i]
			/ ( currGrid->GetYLeft()->vec[flame->fMassFlowRate] );
	}

	return 0;
}

void TUnstrPremFlamePhys::CalcAllDiffVeloRhoYiNext( Double *YCurr, Double hNext )
{
// returns rho Y_i V_i between boundary and first gridpoint
	int		i, nSpecIn = fSpecies->GetNSpeciesInSystem();

	Double	*YNext = fFlameNode->Y[kNext];
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	*rhoY_iV_iPlus = fRhoY_iV_iPlus->vec;
	Double	rhoCurr = fFlameNode->mixDensity[kCurr];
	Double	rhoNext = fFlameNode->mixDensity[kNext];
	Double	diffPlus, rhoV_cPlus = 0.0;
	
	diffPlus = - 0.5 * ( rhoCurr /*/ WC*/ + rhoNext /*/ WN*/ );
	for ( i = 0; i < nSpecIn; ++i ) {
		rhoY_iV_iPlus[i] = diffPlus * 0.5 * ( diffusivity[i] + diffusivityNext[i] ) * ( YNext[i] /** WN*/ - YCurr[i] /** WC*/ ) / hNext;
		rhoV_cPlus -= rhoY_iV_iPlus[i];
	}
	for ( i = 0; i < nSpecIn; ++i ) {
		rhoY_iV_iPlus[i] += 0.5 * ( YCurr[i] + YNext[i] ) * rhoV_cPlus;
	}
}

void TUnstrPremFlamePhys::CompleteSpeciesDiffusion( NodeInfoPtr nodeInfo )
{
// returns    d/dy( rho D_k dY_k/dy)

	int		speciesEq, i;
	int 	fFirstSpecies = GetOffsetFirstSpecies();
	int		tempOff = GetOffsetTemperature();
	int		nSpeciesInSystem = GetSpecies()->GetNSpeciesInSystem();
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*YP = &nodeInfo->yPrev[fFirstSpecies];
	Double	*YC = &nodeInfo->y[fFirstSpecies];
	Double	*YN = &nodeInfo->yNext[fFirstSpecies];
	Double	*rhsOff = &nodeInfo->rhs[fFirstSpecies];
	Double	*rhs = nodeInfo->rhs;
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*W = fFlameNode->mixMolarMass;
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	diffMinus, diffPlus;
	Double	rhoV_cMin = 0.0, rhoV_cPlus = 0.0;
	Double	constDiffMinus, constDiffPlus;	
	Double	constDiffThermMinus, constDiffThermPlus;	
	Double	*thermDiffPrev, *thermDiffCurr, *thermDiffNext;
	Double	hPlushm = nodeInfo->h + nodeInfo->hm;

	if ( fThermoDiffusion ) {
		thermDiffPrev = fFlameNode->diffTherm[kPrev];
		thermDiffCurr = fFlameNode->diffTherm[kCurr];
		thermDiffNext = fFlameNode->diffTherm[kNext];
		constDiffThermMinus = -0.5 * ( 1.0 / y[tempOff] + 1.0 / yPrev[tempOff] ) * ( y[tempOff] - yPrev[tempOff] ) /  (nodeInfo->hm * hPlushm);
		constDiffThermPlus = -0.5 * ( 1.0 / y[tempOff] + 1.0 / yNext[tempOff] ) * ( yNext[tempOff] - y[tempOff] ) / (nodeInfo->h * hPlushm);
	}
	
	constDiffMinus = - 0.5 * ( mixDensity[kCurr] / W[kCurr] + mixDensity[kPrev] / W[kPrev] ) / (nodeInfo->hm * hPlushm);
	constDiffPlus = - 0.5 * ( mixDensity[kCurr] / W[kCurr] + mixDensity[kNext] / W[kNext] ) / (nodeInfo->h * hPlushm);
	Double	rhoY_iV_iMin, rhoY_iV_iPlus, YCW;
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		diffMinus = constDiffMinus * ( diffusivityPrev[i] + diffusivity[i] );
		diffPlus = constDiffPlus * ( diffusivity[i] + diffusivityNext[i] );
		YCW = YC[i] * W[kCurr];
	
		rhoY_iV_iMin = /*0.5 * */diffMinus * ( YCW - YP[i] * W[kPrev] );
		rhoY_iV_iPlus = /*0.5 * */diffPlus * ( YN[i] * W[kNext] - YCW );

		if ( fThermoDiffusion ) {
			rhoY_iV_iMin += constDiffThermMinus * /*0.5 **/ ( thermDiffPrev[i]/* / yPrev[tempOff]*/ + thermDiffCurr[i]/* / y[tempOff]*/ );
			rhoY_iV_iPlus += constDiffThermPlus * /*0.5 **/ ( thermDiffCurr[i]/* / y[tempOff]*/ + thermDiffNext[i]/* / yNext[tempOff]*/ );
		}

		rhsOff[i] += /*2.0 **/ ( rhoY_iV_iPlus - rhoY_iV_iMin );
		
		rhoV_cMin -= rhoY_iV_iMin;
		rhoV_cPlus -= rhoY_iV_iPlus;
	}
	rhoV_cMin *= 0.5;
	rhoV_cPlus *= 0.5;

	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq; ++speciesEq ) {
		rhs[speciesEq] += ( ( y[speciesEq] + yNext[speciesEq] ) * rhoV_cPlus
									- ( y[speciesEq] + yPrev[speciesEq] ) * rhoV_cMin )
									;
	}
}

