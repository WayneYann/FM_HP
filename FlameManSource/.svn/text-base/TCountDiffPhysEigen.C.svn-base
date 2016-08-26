#include "FlameMaster.h"
#include "ListTool.h"
#include "Spline.h"
#include "TCountDiffPhysEigen.h"

// following defines inflow velocity boundary conditions to be mass fluxes in unites [kg/s m^2] versus velocities
#define BCINFLUX

// simple mixing is Fick's law with individual diffusivities for all species.
#undef SIMPLEMIXING

// convection term upwind or central
#define UPWINDCONV

// don't touch this
#define CHECKIT

// following allows for recombination bc at a porous plug burner. Sets O, OH, and H to zero.
#undef RECOMBINATION

// following allows for diffusion description with D_ij
#undef FULLDIFFUSION




#ifdef SIMPLEMIXING
#	undef MOLARDIFFUSION
#	undef DIFFUSIVITYCORRECTION
#else
#	define MOLARDIFFUSION
#	define DIFFUSIVITYCORRECTION
#endif

void TCountDiffPhysEigen::InitTCountDiffPhysEigen( void )
{
	int i;
	TBVPSolverPtr	solver = GetSolver();
	TNewtonPtr		bt = solver->bt;
	TGridPtr		fine = bt->GetGrid()->GetFine();
	TGridPtr		coarse = bt->GetGrid()->GetCoarse();
	int				nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	int				maxGridPoints = bt->GetMaxGridPoints();

	if ( fSoot ) {
		fSoot->SetMomentsOffset( fSootMoments );
	}

//	names of variables
	fVariableNames = new String[fVariablesWithoutSpecies + nSpeciesInSystem];

	fVariableNames[fVVelocity] = new char[2];
	strcpy( fVariableNames[fVVelocity], "V" );
	fVariableNames[fUVelocity] = new char[2];
	strcpy( fVariableNames[fUVelocity], "G" );
	fVariableNames[fMixFrac] = new char[2];
	strcpy( fVariableNames[fMixFrac], "Z" );
	fVariableNames[fPStrain] = new char[2];
	strcpy( fVariableNames[fPStrain], "P" );
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
	fLiquidPoolBC = fInputData->fLiquidPoolBC;
	fKappa = fInputData->fKappa;

//	vectors of solution
	fSolV = NewVector( maxGridPoints + 2 );
	fSolG = NewVector( maxGridPoints + 2 );
	fSolMixFrac = NewVector( maxGridPoints + 2 );
	fSolP = NewVector( maxGridPoints + 2 );

	fSolV->vec = &fSolV->vec[kNext];
	fSolG->vec = &fSolG->vec[kNext];
	fSolMixFrac->vec = &fSolMixFrac->vec[kNext];
	fSolP->vec = &fSolP->vec[kNext];

	fSolV->len -= 2;
	fSolG->len -= 2;
	fSolMixFrac->len -= 2;
	fSolP->len -= 2;

//	saved solution
	fSavedV = NewVector( maxGridPoints + 2 );
	fSavedG = NewVector( maxGridPoints + 2 );
	fSavedMixFrac = NewVector( maxGridPoints + 2 );
	fSavedP = NewVector( maxGridPoints + 2 );

	fSavedV->vec = &fSavedV->vec[kNext];
	fSavedG->vec = &fSavedG->vec[kNext];
	fSavedMixFrac->vec = &fSavedMixFrac->vec[kNext];
	fSavedP->vec = &fSavedP->vec[kNext];

	fSavedV->len -= 2;
	fSavedG->len -= 2;
	fSavedMixFrac->len -= 2;
	fSavedP->len -= 2;

	fYLeftVec = NewVector( fSpecies->GetNOfSpecies() );
	fYRightVec = NewVector( fSpecies->GetNOfSpecies() );
	fRhoY_iV_iPlus = NewVector( fSpecies->GetNSpeciesInSystem() );
	fYLeftVec->len = nSpeciesInSystem;
	fYRightVec->len = nSpeciesInSystem;

	fBC_Left = &fine->GetBcLeft()->vec[ fFirstSpecies ];	

	fNewtonInfoL = NewNewtonInfo( fYLeftVec->len, fYLeftVec );
	fNewtonInfoL->modified = FALSE;
	fNewtonInfoL->maxSteps = 100;
	SetNewtonFuncs( fNewtonInfoL, BCLeftNewtonFuncs, NULL, NULL );

	fNewtonInfoR = NewNewtonInfo( fYRightVec->len, fYRightVec );
	fNewtonInfoR->modified = FALSE;
	fNewtonInfoR->maxSteps = 100;
	SetNewtonFuncs( fNewtonInfoR, BCRightNewtonFuncs, NULL, NULL );
	
	if ( fUseNumericalJac ) {
		bt->SetUtFuncs( NULL, NULL, NULL
					, DiffPhysEigenRHSRest, DiffPhysEigenRHSRest, DiffPhysEigenRHSRest 
					, DiffPhysEigenOutput, DiffPhysEigenPostIter
					, SetDiffPhysEigenNodeInfo, DiffPhysEigenPostConv
					, GetDiffPhysEigenVarNames
					, DiffPhysEigenUpdateLeftBoundary, DiffPhysEigenUpdateRightBoundary );
	}
	else {
		bt->SetUtFuncs( DiffPhysEigenJacFirst, DiffPhysEigenJacRest, DiffPhysEigenJacLast
					, DiffPhysEigenRHSRest, DiffPhysEigenRHSRest, DiffPhysEigenRHSRest 
					, DiffPhysEigenOutput, DiffPhysEigenPostIter
					, SetDiffPhysEigenNodeInfo, DiffPhysEigenPostConv
					, GetDiffPhysEigenVarNames );
	}
	SetInitialBC( fine, fInputData );
	SetInitialBC( coarse, fInputData );
	fMassFraction = new TMassFraction( this );
	if ( !fMassFraction ) FatalError( "memory allocation of TMassFraction failed" );
	ReadStartProfiles( fInputData );
	CheckBC();
	CheckInitialGuess();
	UpdateSolution( fine->GetY(), fine->GetYLeft(), fine->GetYRight() );	
}

TCountDiffPhysEigen::~TCountDiffPhysEigen( void )
{
	int	nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();

	delete fMassFraction;

	fSolP->vec = &fSolP->vec[kPrev];
	fSolMixFrac->vec = &fSolMixFrac->vec[kPrev];
	fSolV->vec = &fSolV->vec[kPrev];
	fSolG->vec = &fSolG->vec[kPrev];

	fSavedP->vec = &fSavedP->vec[kPrev];
	fSavedMixFrac->vec = &fSavedMixFrac->vec[kPrev];
	fSavedV->vec = &fSavedV->vec[kPrev];
	fSavedG->vec = &fSavedG->vec[kPrev];

	FreeNewtonInfo( fNewtonInfoR );
	FreeNewtonInfo( fNewtonInfoL );
	DisposeVector( fYLeftVec );
	DisposeVector( fRhoY_iV_iPlus );
	DisposeVector( fSavedP );
	DisposeVector( fSavedMixFrac );
	DisposeVector( fSavedG );
	DisposeVector( fSavedV );

	DisposeVector( fSolP );
	DisposeVector( fSolMixFrac );
	DisposeVector( fSolG );
	DisposeVector( fSolV );

	for ( int i = 0; i < nSpeciesInSystem+fVariablesWithoutSpecies; ++i ) {
		delete fVariableNames[i];
	}
	delete fVariableNames;
}

void DiffPhysEigenJacFirst( void *object, NodeInfoPtr nodeInfo )
{
	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fPStrain = flame->GetOffsetPStrain();
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		M = nodeInfo->nOfEquations;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
	int		mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
    Double  hnenn = nodeInfo->hnenn;
	Double	**a = nodeInfo->a;
	Double	**massFracs = flameNode->Y;
	Double	*Y = massFracs[kCurr];
	Double	*YPrev = massFracs[kPrev];
	Double	*temp = flameNode->temp;
	Double	*y = nodeInfo->y;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*mixDensity = flameNode->mixDensity;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*mixHeatCapacity = flameNode->mixHeatCapacity;
	Double	oneOverCp = 1.0 / mixHeatCapacity[kCurr];
	Double	*enthalpy = flameNode->enthalpy;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	*mixConductivity = flameNode->mixConductivity;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];
	Double	geometry = 1.0 + flame->GetGeometry();
	Double	diffMinusH;

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

	if ( flame->GetSoot() ) {
		if ( !flame->fUseNumericalDM )  {
			flame->GetSoot()->UpdateJacobian( flame );
		}
		flame->GetSoot()->FillJacobi( flame, nodeInfo, kPhysical );
	}

// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );
	FillJacFirstDerivDown( fPStrain, fPStrain, nodeInfo );
	if ( flame->fLiquidPoolBC ) {
		Double	hL = 360064.8;		// latent heat of vaporization of the fuel at T = 373 K ( boiling point )
		a[fTemperature][fVVelocity] -= h / hm * ( h + hm )
				* 0.5 * ( flameNode->mixConductivity[kPrev] + flameNode->mixConductivity[kCurr] ) / hL;
	}

	Double	hCoeff = nodeInfo->h * ( nodeInfo->h + nodeInfo->hm );

// second equation ( momentum )
	flame->FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0 );
// third equation ( mixturefraction )
	flame->FillJacNonlinearConvectUpwind( fVVelocity, fMixFrac, nodeInfo, 1.0 );
	if ( mixtureSpecificationLeft == kMassFlux && y[fVVelocity] > 0.0 ) {
		Double	zCoeff = mixConductivity[kPrev] / ( mixHeatCapacity[kPrev] * yPrev[fVVelocity] * nodeInfo->hm );
		a[fMixFrac][fMixFrac] -= hCoeff * y[fVVelocity] * zCoeff / ( 1.0 + zCoeff );
	}
// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		flame->FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0 );
		if ( mixtureSpecificationLeft == kMassFlux && y[fVVelocity] > 0.0 ) {
			a[speciesEq][speciesEq] -= hCoeff * y[fVVelocity] 
									* flame->GetdYPrevdY( speciesEq-fFirstSpecies, nodeInfo );
		}
	}
	if ( fTemperature < M ) {
		flame->FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0 );
	}

// first equation ( mass )
	a[fUVelocity][fVVelocity] += geometry * mixDensity[kCurr] * hnenn;
	a[fTemperature][fVVelocity] -= geometry * mixDensity[kCurr] * y[fUVelocity] / temp[kCurr] * hnenn;
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		a[speciesEq][fVVelocity] -= geometry * mixDensity[kCurr] * y[fUVelocity] 
									* mixMolarMass / molarMass[speciesEq-fFirstSpecies] * hnenn;
	}
	
// second equation ( momentum )
	FillJacWithDiffusion( fUVelocity, fUVelocity, 1.0, mixViscosity, nodeInfo, kNegative );
	a[fUVelocity][fUVelocity] += 2.0 * mixDensity[kCurr] * y[fUVelocity] * hnenn;
	a[fTemperature][fUVelocity] -= mixDensity[kCurr] * y[fUVelocity] * y[fUVelocity] / temp[kCurr] * hnenn;
	a[fPStrain][fUVelocity] -= hnenn;
	
// third equation ( mixturefraction )
	flame->FillJacMixFracDiffusion( fMixFrac, fMixFrac, nodeInfo, kNegative );
	if ( mixtureSpecificationLeft == kMassFlux ) {
		diffMinusH = ( mixConductivity[kPrev] / mixHeatCapacity[kPrev]
						+ mixConductivity[kCurr] / mixHeatCapacity[kCurr] ) * nodeInfo->h;
		Double	zCoeff = mixConductivity[kPrev] / ( mixHeatCapacity[kPrev] * yPrev[fVVelocity] * nodeInfo->hm );
		a[fMixFrac][fMixFrac] -= diffMinusH * zCoeff / ( 1.0 + zCoeff );
	}

// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, speciesIndexEq, 1.0, nodeInfo, kNegative );

		if ( mixtureSpecificationLeft == kMassFlux ) {
			Double	*diffusivity = flameNode->diffusivity;
			Double	*diffusivityPrev = flameNode->diffusivityPrev;
			diffMinusH = ( diffusivityPrev[speciesIndexEq] * mixDensity[kPrev]
							+ diffusivity[speciesIndexEq] * mixDensity[kCurr] ) * nodeInfo->h;
			a[speciesEq][speciesEq] -= diffMinusH * flame->GetdYPrevdY( speciesIndexEq, nodeInfo );
		}

#ifdef DIFFUSIVITYCORRECTION
		flame->FillJacDiffCorr( speciesEq, 1.0, nodeInfo );

		if ( mixtureSpecificationLeft == kMassFlux ) {
			Double	*diffusivity = flameNode->diffusivity;
			Double	*diffusivityPrev = flameNode->diffusivityPrev;
			Double	coeffCurr = mixDensity[kCurr] * Y[speciesIndexEq];
			Double	coeffPrev = mixDensity[kPrev] * YPrev[speciesIndexEq];
			Double	sumYCurr = 0.0;
			Double	sumYPrev = 0.0;
			
			for ( speciesVar = 0; speciesVar < nSpeciesInSystem; ++speciesVar ) {
				sumYCurr += Y[speciesVar];
				sumYPrev += YPrev[speciesVar];
			}
			coeffCurr /= sumYCurr;
			coeffPrev /= sumYPrev;
			
			for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq && speciesVar < M; ++ speciesVar ) {
				speciesIndexVar = speciesVar - fFirstSpecies;
				diffMinusH = nodeInfo->h * ( coeffCurr * diffusivity[speciesIndexVar] 
											+ coeffPrev * diffusivityPrev[speciesIndexVar] );
				a[speciesVar][speciesEq] -= diffMinusH * flame->GetdYPrevdY( speciesIndexVar, nodeInfo );
			}
		}
#endif
		//	implicit source term
		a[fTemperature][speciesEq] -= dMdT[speciesIndexEq] * hnenn;
		for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq && speciesVar < M; ++ speciesVar ) {
			speciesIndexVar = speciesVar - fFirstSpecies;
			a[speciesVar][speciesEq] -= dMdY[speciesIndexVar][speciesIndexEq] * hnenn;
		}
	}

// last equation ( temperature )
	if ( fTemperature < M ) {
		FillJacWithDiffusion( fTemperature, fTemperature, oneOverCp, mixConductivity, nodeInfo, kNegative );
		for ( speciesIndexEq = 0; speciesIndexEq < nSpeciesInSystem; ++speciesIndexEq ) {
			a[fTemperature][fTemperature] += oneOverCp * ( dMdT[speciesIndexEq] * enthalpy[speciesIndexEq] ) * hnenn;
			for ( speciesIndexVar = 0; speciesIndexVar < nSpeciesInSystem; ++speciesIndexVar ) {
				a[speciesIndexVar + fFirstSpecies][fTemperature] += oneOverCp * dMdY[speciesIndexVar][speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
			}
		}
		if ( flame->fProperties->GetRadiation() ) {
			flame->fProperties->GetRadiation()->FillJacRadiation( -oneOverCp, flame, nodeInfo );
			if ( flame->GetSoot() && flame->GetSoot()->WithSootRadiation() ) {
				flame->GetSoot()->FillJacSootRadiation( oneOverCp, flame, nodeInfo );
			}
		}
	}
	if ( flame->GetSolver()->bt->GetTimedepFlag() 
				&& !flame->GetSolver()->time->GetTimeConverged() ) {
		TTimePtr tim = flame->GetSolver()->time;
		for ( int eqLoop = 0; eqLoop < M; ++eqLoop ) {
			if ( eqLoop != fPStrain ) {
				a[eqLoop][eqLoop] += hnenn / tim->GetDeltaT();
			}
		}
	}
}

Double TCountDiffPhysEigen::GetdYPrevdY( int speciesIndex, NodeInfoPtr nodeInfo )
{
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	VWall = nodeInfo->yPrev[fVVelocity];
	Double	coeff = mixDensity[kPrev] * diffusivityPrev[speciesIndex]
								/ ( nodeInfo->hm * VWall );
								
#ifdef DIFFUSIVITYCORRECTION
	Double	coeffCorr = mixDensity[kPrev] / VWall * fFlameNode->diffCorr[kPrev];
#else
	Double	coeffCorr = 0.0;
#endif

	return coeff / ( 1.0 + coeff + coeffCorr );
}

void DiffPhysEigenJacRest( void *object, NodeInfoPtr nodeInfo )
{
	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fPStrain = flame->GetOffsetPStrain();
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
	Double	mixDensity = *flameNode->mixDensity;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	oneOverCp = 1.0 / flameNode->mixHeatCapacity[kCurr];
	Double	*enthalpy = flameNode->enthalpy;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];
	Double	geometry = 1.0 + flame->GetGeometry();

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

	if ( flame->GetSoot() ) {
		if ( !flame->fUseNumericalDM )  {
			flame->GetSoot()->UpdateJacobian( flame );
		}
		flame->GetSoot()->FillJacobi( flame, nodeInfo, kPhysical );
	}

// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );
	FillJacFirstDerivDown( fPStrain, fPStrain, nodeInfo );

	// second equation ( momentum )
	flame->FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0 );
	// third equation ( mixturefraction )
	flame->FillJacNonlinearConvectUpwind( fVVelocity, fMixFrac, nodeInfo, 1.0 );
	// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		flame->FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0 );
	}
	if ( fTemperature < M ) {
		flame->FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0 );
	}

	
// first equation ( mass )
	a[fUVelocity][fVVelocity] += geometry * mixDensity * hnenn;
	a[fTemperature][fVVelocity] -= geometry * mixDensity * y[fUVelocity] / temp[kCurr] * hnenn;
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		a[speciesEq][fVVelocity] -= geometry * mixDensity * y[fUVelocity] 
									* mixMolarMass / molarMass[speciesEq-fFirstSpecies] * hnenn;
	}
	
// second equation ( momentum )
	FillJacWithDiffusion( fUVelocity, fUVelocity, 1.0, mixViscosity, nodeInfo, kNegative );
	a[fUVelocity][fUVelocity] += 2.0 * mixDensity * y[fUVelocity] * hnenn;
	a[fTemperature][fUVelocity] -= mixDensity * y[fUVelocity] * y[fUVelocity] / temp[kCurr] * hnenn;
	a[fPStrain][fUVelocity] -= hnenn;
	
// third equation ( mixturefraction )
	flame->FillJacMixFracDiffusion( fMixFrac, fMixFrac, nodeInfo, kNegative );

// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, speciesIndexEq, 1.0, nodeInfo, kNegative );
#ifdef DIFFUSIVITYCORRECTION
		flame->FillJacDiffCorr( speciesEq, 1.0, nodeInfo );
#endif
		if ( flame->fThermoDiffusion ) {
			flame->FillJacThermoDiffusion( speciesEq, -1.0, kPhysical, nodeInfo );
		}
		//	implicit source term
		a[fTemperature][speciesEq] -= dMdT[speciesIndexEq] * hnenn;
		for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq && speciesVar < M; ++ speciesVar ) {
			speciesIndexVar = speciesVar - fFirstSpecies;
			a[speciesVar][speciesEq] -= dMdY[speciesIndexVar][speciesIndexEq] * hnenn;
		}
	}

// last equation ( temperature )
	if ( fTemperature < M ) {
		FillJacWithDiffusion( fTemperature, fTemperature, oneOverCp, flameNode->mixConductivity, nodeInfo, kNegative );
		for ( speciesIndexEq = 0; speciesIndexEq < nSpeciesInSystem; ++speciesIndexEq ) {
			a[fTemperature][fTemperature] += oneOverCp * ( dMdT[speciesIndexEq] * enthalpy[speciesIndexEq] ) * hnenn;
			for ( speciesIndexVar = 0; speciesIndexVar < nSpeciesInSystem; ++speciesIndexVar ) {
				a[speciesIndexVar + fFirstSpecies][fTemperature] += oneOverCp * dMdY[speciesIndexVar][speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
			}
		}
		if ( flame->fProperties->GetRadiation() ) {
			flame->fProperties->GetRadiation()->FillJacRadiation( -oneOverCp, flame, nodeInfo );
			if ( flame->GetSoot() && flame->GetSoot()->WithSootRadiation() ) {
				flame->GetSoot()->FillJacSootRadiation( oneOverCp, flame, nodeInfo );
			}
		}
	}
	if ( flame->GetSolver()->bt->GetTimedepFlag() 
				&& !flame->GetSolver()->time->GetTimeConverged() ) {
		TTimePtr tim = flame->GetSolver()->time;
		for ( int eqLoop = 0; eqLoop < M; ++eqLoop ) {
			if ( eqLoop != fPStrain ) {
				a[eqLoop][eqLoop] += hnenn / tim->GetDeltaT();
			}
		}
	}
}

void DiffPhysEigenJacLast( void *object, NodeInfoPtr nodeInfo )
{
	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fPStrain = flame->GetOffsetPStrain();
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		M = nodeInfo->nOfEquations;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
	int		mixtureSpecificationRight = flame->GetMixtureSpecificationRight();
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
    Double  hnenn = nodeInfo->hnenn;
	Double	**a = nodeInfo->a;
	Double	**massFracs = flameNode->Y;
	Double	*Y = massFracs[kCurr];
	Double	*YNext = massFracs[kNext];
	Double	*temp = flameNode->temp;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*mixDensity = flameNode->mixDensity;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*mixHeatCapacity = flameNode->mixHeatCapacity;
	Double	oneOverCp = 1.0 / mixHeatCapacity[kCurr];
	Double	*enthalpy = flameNode->enthalpy;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	*mixConductivity = flameNode->mixConductivity;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];
	Double	geometry = 1.0 + flame->GetGeometry();
	Double	diffPlusHm;
#ifdef RECOMBINATION
	char	**names = flame->GetSpecies()->GetNames();
#endif

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

	if ( flame->GetSoot() ) {
		if ( !flame->fUseNumericalDM )  {
			flame->GetSoot()->UpdateJacobian( flame );
		}
		flame->GetSoot()->FillJacobi( flame, nodeInfo, kPhysical );
	}

// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );

	FillJacFirstDerivDown( fVVelocity, fPStrain, nodeInfo );
	a[fUVelocity][fPStrain] += 
				( 1.0 + flame->GetGeometry() ) * mixDensity[kCurr] 
								* hm * h * ( h + hm );

	Double	hCoeff = nodeInfo->hm * ( nodeInfo->h + nodeInfo->hm );

// second equation ( momentum )
	flame->FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0 );
// third equation ( mixturefraction )
	flame->FillJacNonlinearConvectUpwind( fVVelocity, fMixFrac, nodeInfo, 1.0 );
	if ( mixtureSpecificationRight == kMassFlux && y[fVVelocity] < 0.0 ) {
		Double	zCoeff = mixConductivity[kNext] / ( mixHeatCapacity[kNext] * yNext[fVVelocity] * nodeInfo->h );
		a[fMixFrac][fMixFrac] -= hCoeff * y[fVVelocity] * zCoeff / ( 1.0 - zCoeff );
	}
// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		flame->FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0 );
		if ( mixtureSpecificationRight == kMassFlux && y[fVVelocity] < 0.0 ) {
#ifdef RECOMBINATION
			int ind = speciesEq - fFirstSpecies;
			if ( strcmp( names[ind], "H" ) != 0 
						&& strcmp( names[ind], "O" ) != 0 
						&& strcmp( names[ind], "OH" ) != 0 ) {
#endif
				a[speciesEq][speciesEq] += hCoeff * y[fVVelocity] 
										* flame->GetdYNextdY( speciesEq-fFirstSpecies, nodeInfo );
#ifdef RECOMBINATION
			}
#endif
		}
	}
	if ( fTemperature < M ) {
		flame->FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0 );
	}

// first equation ( mass )
	a[fUVelocity][fVVelocity] += geometry * mixDensity[kCurr] * hnenn;
	a[fTemperature][fVVelocity] -= geometry * mixDensity[kCurr] * y[fUVelocity] / temp[kCurr] * hnenn;
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		a[speciesEq][fVVelocity] -= geometry * mixDensity[kCurr] * y[fUVelocity] 
									* mixMolarMass / molarMass[speciesEq-fFirstSpecies] * hnenn;
	}
	
// second equation ( momentum )
	FillJacWithDiffusion( fUVelocity, fUVelocity, 1.0, mixViscosity, nodeInfo, kNegative );
	a[fUVelocity][fUVelocity] += 2.0 * mixDensity[kCurr] * y[fUVelocity] * hnenn;
	a[fTemperature][fUVelocity] -= mixDensity[kCurr] * y[fUVelocity] * y[fUVelocity] / temp[kCurr] * hnenn;
	a[fPStrain][fUVelocity] -= hnenn;
	
// third equation ( mixturefraction )
	flame->FillJacMixFracDiffusion( fMixFrac, fMixFrac, nodeInfo, kNegative );
	if ( mixtureSpecificationRight == kMassFlux ) {
		diffPlusHm = ( mixConductivity[kNext] / mixHeatCapacity[kNext]
						+ mixConductivity[kCurr] / mixHeatCapacity[kCurr] ) * nodeInfo->hm;
		Double	zCoeff = mixConductivity[kNext] / ( mixHeatCapacity[kNext] * yNext[fVVelocity] * nodeInfo->h );
		a[fMixFrac][fMixFrac] += diffPlusHm * zCoeff / ( 1.0 - zCoeff );
	}

// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, speciesIndexEq, 1.0, nodeInfo, kNegative );

		if ( mixtureSpecificationRight == kMassFlux ) {
#ifdef RECOMBINATION
			int ind = speciesEq - fFirstSpecies;
			if ( strcmp( names[ind], "H" ) != 0 
						&& strcmp( names[ind], "O" ) != 0 
						&& strcmp( names[ind], "OH" ) != 0 ) {
#endif
				Double	*diffusivity = flameNode->diffusivity;
				Double	*diffusivityNext = flameNode->diffusivityNext;
				Double	diffPlusHm = ( diffusivityNext[speciesIndexEq] * mixDensity[kNext]
								+ diffusivity[speciesIndexEq] * mixDensity[kCurr] ) * nodeInfo->hm;
				a[speciesEq][speciesEq] -= diffPlusHm * flame->GetdYNextdY( speciesIndexEq, nodeInfo );
#ifdef RECOMBINATION
			}
#endif
		}

#ifdef DIFFUSIVITYCORRECTION
		flame->FillJacDiffCorr( speciesEq, 1.0, nodeInfo );

		if ( mixtureSpecificationRight == kMassFlux ) {
			Double	*diffusivity = flameNode->diffusivity;
			Double	*diffusivityNext = flameNode->diffusivityNext;
			Double	coeffCurr = mixDensity[kCurr] * Y[speciesIndexEq];
			Double	coeffNext = mixDensity[kNext] * YNext[speciesIndexEq];
			Double	sumYCurr = 0.0;
			Double	sumYNext = 0.0;
			
			for ( speciesVar = 0; speciesVar < nSpeciesInSystem; ++speciesVar ) {
				sumYCurr += Y[speciesVar];
				sumYNext += YNext[speciesVar];
			}
			coeffCurr /= sumYCurr;
			coeffNext /= sumYNext;
			
			for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq && speciesVar < M; ++ speciesVar ) {
				speciesIndexVar = speciesVar - fFirstSpecies;
				diffPlusHm = nodeInfo->hm * ( coeffCurr * diffusivity[speciesIndexVar] 
											+ coeffNext * diffusivityNext[speciesIndexVar] );
				a[speciesVar][speciesEq] += diffPlusHm * flame->GetdYNextdY( speciesIndexVar, nodeInfo );
			}
		}
#endif
		//	implicit source term
		a[fTemperature][speciesEq] -= dMdT[speciesIndexEq] * hnenn;
		for ( speciesVar = fFirstSpecies; speciesVar < lastSpeciesEq && speciesVar < M; ++ speciesVar ) {
			speciesIndexVar = speciesVar - fFirstSpecies;
			a[speciesVar][speciesEq] -= dMdY[speciesIndexVar][speciesIndexEq] * hnenn;
		}
	}

// last equation ( temperature )
	if ( fTemperature < M ) {
		FillJacWithDiffusion( fTemperature, fTemperature, oneOverCp, flameNode->mixConductivity, nodeInfo, kNegative );
		for ( speciesIndexEq = 0; speciesIndexEq < nSpeciesInSystem; ++speciesIndexEq ) {
			a[fTemperature][fTemperature] += oneOverCp * ( dMdT[speciesIndexEq] * enthalpy[speciesIndexEq] ) * hnenn;
			for ( speciesIndexVar = 0; speciesIndexVar < nSpeciesInSystem; ++speciesIndexVar ) {
				a[speciesIndexVar + fFirstSpecies][fTemperature] += oneOverCp * dMdY[speciesIndexVar][speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
			}
		}
		if ( flame->fProperties->GetRadiation() ) {
			flame->fProperties->GetRadiation()->FillJacRadiation( -oneOverCp, flame, nodeInfo );
			if ( flame->GetSoot() && flame->GetSoot()->WithSootRadiation() ) {
				flame->GetSoot()->FillJacSootRadiation( oneOverCp, flame, nodeInfo );
			}
		}
	}
	if ( flame->GetSolver()->bt->GetTimedepFlag() 
				&& !flame->GetSolver()->time->GetTimeConverged() ) {
		TTimePtr tim = flame->GetSolver()->time;
		for ( int eqLoop = 0; eqLoop < M; ++eqLoop ) {
			if ( eqLoop != fPStrain ) {
				a[eqLoop][eqLoop] += hnenn / tim->GetDeltaT();
			}
		}
	}
}

Double TCountDiffPhysEigen::GetdYNextdY( int speciesIndex, NodeInfoPtr nodeInfo )
{
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	VWall = nodeInfo->yNext[fVVelocity];
	Double	coeff = mixDensity[kNext] * diffusivityNext[speciesIndex]
								/ ( nodeInfo->h * VWall );
								
#ifdef DIFFUSIVITYCORRECTION
	Double	coeffCorr = mixDensity[kNext] / VWall * fFlameNode->diffCorr[kNext];
#else
	Double	coeffCorr = 0.0;
#endif

	return - coeff / ( 1.0 - coeff + coeffCorr );
}

void DiffPhysEigenRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	if ( !flame->RHSAction( nodeInfo, rhsMode ) ) {
		return;
	}
	
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fPStrain = flame->GetOffsetPStrain();
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		eqLoop, speciesEq;
	int		M = nodeInfo->nOfEquations;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
    Double  h2 = h * h;
    Double  hm2 = hm * hm;
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
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*productionRate = flameNode->productionRate;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	oneOverCp = 1.0 / flameNode->mixHeatCapacity[kCurr];
	Double	sumCpDdYdx;
	Double	sumMH;
	
	if ( nodeInfo->lastPoint == TRUE ) {
//		DiffPhysEigenUpdateRightBoundary( flame );
	}


	if ( flame->GetSoot() ) {
		flame->GetSoot()->FillRHS( flame, nodeInfo, kPhysical );
	}


// first fill all convection terms
	// first equation ( mass )


	if ( nodeInfo->firstPoint == TRUE ) {
		if ( flame->fLiquidPoolBC ) {
			Double	hL = 360064.8;		// latent heat of vaporization of the fuel at T = 373 K ( boiling point )
			yPrev[fVVelocity] = FirstDerivUpwind( y[fTemperature], yPrev[fTemperature], hm ) 
					* 0.5 * ( flameNode->mixConductivity[kPrev] + flameNode->mixConductivity[kCurr] ) / hL;
		}
//		flame->CalcYLeft();
		Double	coeff;

		coeff = flameNode->mixConductivity[kPrev] / ( flameNode->mixHeatCapacity[kPrev] * yPrev[fVVelocity] * hm );
		yPrev[fMixFrac] = ( 1.0 + coeff * y[fMixFrac] )
						/ ( 1.0 + coeff );

		yPrev[fPStrain] = y[fPStrain];
	}

	if ( nodeInfo->firstPoint == TRUE ) {
		//flame->CalcYLeft();
		Double	VPrev;
		if ( flame->fLiquidPoolBC ) {
			Double	hL = 360064.8;		// latent heat of vaporization of the fuel at T = 373 K ( boiling point )
			VPrev = FirstDerivUpwind( y[fTemperature], yPrev[fTemperature], hm ) 
					* 0.5 * ( flameNode->mixConductivity[kPrev] + flameNode->mixConductivity[kCurr] ) / hL;
		}
		else {
			VPrev = flame->GetSolver()->bt->GetGrid()->GetCurrentGrid()->GetBcLeft()->vec[fVVelocity];
		}
		rhs[fVVelocity] += FirstDerivUpwind( y[fVVelocity], VPrev, hm );
		rhs[fPStrain] += FirstDerivUpwind( yNext[fPStrain], y[fPStrain], h );
	}
	else {
		rhs[fVVelocity] += FirstDerivUpwind( y[fVVelocity], yPrev[fVVelocity], hm );
	}
	if ( nodeInfo->lastPoint == TRUE ) {
		rhs[fPStrain] += FirstDerivUpwind( yNext[fVVelocity], y[fVVelocity], h ) 
/* ATTENTION HP */
		+ 
						( 1.0 + flame->GetGeometry() ) * mixDensity[kCurr] * y[fUVelocity]
						;
	}
	else {
		rhs[fPStrain] += FirstDerivUpwind( yNext[fPStrain], y[fPStrain], h );
	}

#ifdef UPWINDCONV
	// second equation ( momentum )
	rhs[fUVelocity] += flame->NonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo );
	// third equation ( mixturefraction )
	rhs[fMixFrac] += flame->NonlinearConvectUpwind( fVVelocity, fMixFrac, nodeInfo );
	// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		rhs[speciesEq] += flame->NonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo );
	}
	if ( fTemperature < M ) {
		rhs[fTemperature] += flame->NonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo );
	}
#else
	// second equation ( momentum )
	rhs[fUVelocity] += NonlinearConvectCentral( fVVelocity, yPrev[fUVelocity], y[fUVelocity], yNext[fUVelocity], hm, h );
	// third equation ( mixturefraction )
	rhs[fMixFrac] += NonlinearConvectCentral( fVVelocity, yPrev[fMixFrac], y[fMixFrac], yNext[fMixFrac], hm, h );
	// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		rhs[speciesEq] += NonlinearConvectCentral( fVVelocity, yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
	}
	if ( fTemperature < M ) {
		rhs[fTemperature] += NonlinearConvectCentral( fVVelocity, yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
	}
#endif

// mass equation
	rhs[fVVelocity] += ( 1.0 + flame->GetGeometry() ) * mixDensity[kCurr] * y[fUVelocity];

// momentum equation
	rhs[fUVelocity] -= SecondDerivDiffusion( fUVelocity, mixViscosity, nodeInfo );
	rhs[fUVelocity] -= y[fPStrain];
	rhs[fUVelocity] += mixDensity[kCurr] * y[fUVelocity] * y[fUVelocity];
	
// mixture fraction equation
	rhs[fMixFrac] -= flame->MixFracDiffusion( fMixFrac, nodeInfo );
	
// fFirstSpecies to fFirstSpecies + nOfSpecies equation ( species )
	Double	*locprodRate = &productionRate[-fFirstSpecies];
	int		locLoopEnd = minint( lastSpeciesEq, M );
	for ( speciesEq = fFirstSpecies; speciesEq < locLoopEnd; ++speciesEq ) {
		rhs[speciesEq] -= locprodRate[speciesEq];
	}
#ifdef FULLDIFFUSION
	for ( speciesEq = fFirstSpecies; speciesEq < locLoopEnd; ++speciesEq ) {
		rhs[speciesEq] += flame->SecondDerivBinSpecDiff( speciesEq, nodeInfo );
	}
#else
	flame->CompleteSpeciesDiffusion( nodeInfo );
#endif
	
// energy equation
	if ( fTemperature < M ) {
		sumCpDdYdx = 0.0;
		sumMH = 0.0;
		Double	diffCorrHeat = 0.0;
		Double	sumCpYD = 0.0;
#	ifdef DIFFUSIVITYCORRECTION
		diffCorrHeat = ( flame->UseDiffCorr() ) ? flameNode->mixHeatCapacity[kCurr] : 0.0;
#	endif

// diffusion
		rhs[fTemperature] -= oneOverCp * SecondDerivDiffusion( fTemperature, flameNode->mixConductivity, nodeInfo );
		
// compute all sums
		for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
			sumCpDdYdx += ( heatCapacity[eqLoop] - diffCorrHeat ) * diffusivity[eqLoop] 
						* FirstDerivSq( YPrev[eqLoop], Y[eqLoop], YNext[eqLoop], hm2, h2, hnenn );
#	ifdef MOLARDIFFUSION
			sumCpYD += ( heatCapacity[eqLoop] - diffCorrHeat ) * Y[eqLoop] * diffusivity[eqLoop];
#	endif
			sumMH += productionRate[eqLoop] * enthalpy[eqLoop];
		}

//	enthalpy flux
		rhs[fTemperature] -= oneOverCp * sumCpDdYdx * mixDensity[kCurr]
						* FirstDerivSq( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm2, h2, hnenn );
#	ifdef MOLARDIFFUSION
		rhs[fTemperature] -= oneOverCp * sumCpYD * mixDensity[kCurr]
						/ flameNode->mixMolarMass[kCurr]
						* FirstDerivSq( flameNode->mixMolarMass[kPrev], flameNode->mixMolarMass[kCurr], flameNode->mixMolarMass[kNext], hm2, h2, hnenn )
						* FirstDerivSq( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm2, h2, hnenn );
#	endif
		
		rhs[fTemperature] += oneOverCp * sumMH;
		
		if ( flame->fProperties->GetRadiation() ) {
			rhs[fTemperature] -= oneOverCp * flameNode->radiation[kCurr];
			if ( flame->GetSoot() && flame->GetSoot()->WithSootRadiation() ) {
				rhs[fTemperature] += oneOverCp * flame->GetSoot()->GetSootRadiation( y[fTemperature], flameNode->moments );
			}
		}
	}

	TTimePtr tim = flame->GetSolver()->time;
	Double	mHnenn = -hnenn;


	if ( flame->GetSolver()->bt->GetTimedepFlag() 
			&& !flame->GetSolver()->time->GetTimeConverged() ) {
		for ( eqLoop = 0; eqLoop < M; ++eqLoop ) {
			if ( eqLoop != fPStrain && eqLoop != fVVelocity ) {
				rhs[eqLoop] += ( y[eqLoop] - tim->GetYOld()->mat[nodeInfo->gridPoint][eqLoop] ) 
								/ tim->GetDeltaT();
			} 
		}
	}
		for ( eqLoop = 0; eqLoop < M; ++eqLoop ) {
			rhs[eqLoop] *= mHnenn;
		}
}

void TCountDiffPhysEigen::FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
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
		if ( !nodeInfo->lastPoint ) {
			b[lVar][nVariable] += diffPlusHm;
		}
		if ( !nodeInfo->firstPoint ) {
			c[lVar][nVariable] += diffMinusH;
		}
	}
}

Double TCountDiffPhysEigen::DiffCorr( int nVariable, NodeInfoPtr nodeInfo )
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
	

	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		diffPlus = coeffCurr * diffusivity[i] + coeffNext * diffusivityNext[i];
		diffMinus = coeffCurr * diffusivity[i] + coeffPrev * diffusivityPrev[i];
		value += ( diffPlus * hm * ( YNext[i] - Y[i] ) 
					+ diffMinus * h * ( YPrev[i] - Y[i] ) );
	}

	return value / nodeInfo->hnenn;
}

void TCountDiffPhysEigen::UpdateDimensions( int len )
{
	T1DFlame::UpdateDimensions( len );
	fSolV->len = len;
	fSolG->len = len;
	fSolMixFrac->len = len;
	fSolP->len = len;
}

void TCountDiffPhysEigen::UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
{
	int		nGridPoints = yMat->cols;
	Double	*V = fSolV->vec;
	Double	*G = fSolG->vec;
	Double	*Z = fSolMixFrac->vec;
	Double	*P = fSolP->vec;
	Double	**y = yMat->mat;
	Double	*yLeft = yLeftVec->vec;
	Double	*yRight = yRightVec->vec;

	UpdateDimensions( nGridPoints );

	T1DFlame::UpdateSolution( yMat, yLeftVec, yRightVec );
	
	V[kPrev] = yLeft[fVVelocity];
	G[kPrev] = yLeft[fUVelocity];
	Z[kPrev] = yLeft[fMixFrac];
	P[kPrev] = yLeft[fPStrain];
	for ( int k = 0; k < nGridPoints; ++k ) {
		V[k] = y[k][fVVelocity];
		G[k] = y[k][fUVelocity];
		Z[k] = y[k][fMixFrac];
		P[k] = y[k][fPStrain];
	}
	V[nGridPoints] = yRight[fVVelocity];
	G[nGridPoints] = yRight[fUVelocity];
	Z[nGridPoints] = yRight[fMixFrac];
	P[nGridPoints] = yRight[fPStrain];
}

int	TCountDiffPhysEigen::GetOffsetVVelocity( void )
{
	return fVVelocity; 
}

int	TCountDiffPhysEigen::GetOffsetUVelocity( void )
{
	return fUVelocity; 
}

int	TCountDiffPhysEigen::GetOffsetTemperature( void )
{
	return fTemperature; 
}

int TCountDiffPhysEigen::GetOffsetMixFrac( void )
{
	return fMixFrac; 
}

int TCountDiffPhysEigen::GetOffsetPStrain( void )
{
	return fPStrain; 
}

int	TCountDiffPhysEigen::GetOffsetFirstSpecies( void ) 
{
	return fFirstSpecies;
}

ConstStringArray TCountDiffPhysEigen::GetVariableNames( void )
{
	return ( ConstStringArray ) fVariableNames;
}

int TCountDiffPhysEigen::GetVariablesWithoutSpecies( void )
{
	return fVariablesWithoutSpecies;
}

#include "TofZ.h"

void TCountDiffPhysEigen::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp )
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
	Flag				USet = FALSE;
//	Flag				sootconcSet = FALSE;
	Flag				oxidizerFound = FALSE;
	Flag				chooseInputGrid = FALSE;
	Double				*xIn = new Double[gridPointsIn];
	if ( !xIn ) FatalError( "memory allocation of TCountDiffPhysEigen failed" );
	Double				*yIn =  new Double[gridPointsIn];
	if ( !yIn ) FatalError( "memory allocation of TCountDiffPhysEigen failed" );
	Double				*yInFloat = sp->data;
	Double				**y = yMat->mat;
	int					variable;
	char				*string = sp->labels;
	SplinePtr			theSpline = NULL;
	Double				leftSlope;
	Double				rightSlope;
	int					oxidizerSide; // this program assumes that oxidizerSide = kRight
	Double				*temp = GetTemperature()->vec;
	Double				**Y = GetMassFracs()->mat;
	FILE				*fp;
	parameter			*param = GetParameter( "strainrate" );
	Double				strainRateIn;

#ifdef BCINFLUX
		if ( yLeft[fVVelocity] == 0.0 ) {
			struct _parameter	*param = GetParameter( "vleft" );
			if ( param ) {
				yLeft[fVVelocity] = (Double)param->what.quantity.value;
				grid->GetBcLeft()->vec[fVVelocity] = yLeft[fVVelocity];
			}
			else {
				fprintf( stderr, "###error: no value for VLeft\n" );
				exit( 2 );
			}
		}
		if ( yRight[fVVelocity] == 0.0 ) {
			struct _parameter	*param = GetParameter( "vright" );
			if ( param ) {
				yRight[fVVelocity] = (Double)param->what.quantity.value;
				grid->GetBcRight()->vec[fVVelocity] = yRight[fVVelocity];
			}
			else {
				fprintf( stderr, "###error: no value for VRight\n" );
				exit( 2 );
			}
		}
#endif

// get strainrate
	if ( param ) {
		strainRateIn = (Double)param->what.quantity.value;
	}
	else {
		cerr << "#warning: no strainrate found in startprofilesfile" << NEWL
			 << "          choose default a = 100 [1/s]" << NEWL;
		
		strainRateIn = 100.0;
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
		cerr << "error: can't find coordinate 'y'" << NEWL;
		exit(2);
	}
	if ( !oxidizerFound ) {
		cerr << "error: can't find massfraction of oxidizer 'massfraction-o2'" << NEWL
					<< "choose default: oxidizerside is right" << NEWL;
		oxidizerSide = kRight;
//		exit(2);
	}
	
	
// set oxidizer side always to right side and thereby disable flipping	
	oxidizerSide = kRight;

// reset string
	string = sp->labels;
	
	for ( i = 0; i < sp->variables; ++i ) {
		if ( strcmp( string, "u" ) == 0 && USet == FALSE ) {
			variable = fUVelocity;
			USet = TRUE;
		}
		else if ( strcmp( string, "df/deta" ) == 0 && USet == FALSE ) {
			variable = fUVelocity;
			USet = TRUE;
		}
		else if ( strcmp( string, "f'" ) == 0 && USet == FALSE ) {
			variable = fUVelocity;
			USet = TRUE;
		}
		else if ( strcmp( string, "g" ) == 0 && USet == FALSE ) {
			variable = fUVelocity;
			USet = TRUE;
		}
		else if ( strncmp( string, "v ", 2 ) == 0 ) {
//			fprintf( stderr, "variable %s found\n", string );
			variable = fVVelocity;
		}
		else if ( strcmp( string, "z" ) == 0 ) {
			variable = fMixFrac;
		}
		else if ( strncmp( string, "temperature", 11 ) == 0 ) {
			variable = fTemperature;
		}
		else if ( strncmp( string, "p ", 2 ) == 0 ) {
			variable = fPStrain;
		}
		else if ( strncmp( string, "conc-soot", 9 ) == 0 && GetSoot() ) {
			string += 9;
			char name[2];
			for ( j = 0; j < GetSoot()->GetNSootMoments(); ++j ) {
				sprintf( name, "%.1d", j );
				if ( strncmp( name, string, 1 ) == 0 ) {
					variable = j + GetSoot()->GetOffsetSootMoments();
					string += 1;
//					sootconcSet = TRUE;
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
				for ( k = 0; k < gridPointsIn-2; ++k ) {
					y[k][variable] = yInFloat[i*gridPointsIn + k+1];	// copy workspace to vector of solution
				}
			}
			else {
				for ( k = 0; k < gridPointsIn-2; ++k ) { // turn vector
					y[k][variable] = yInFloat[(i+1)*gridPointsIn - k-2];	// copy workspace to vector of solution
				}
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
			if ( oxidizerSide == kRight ) {
				for ( k = 0; k < nGridPoints; ++k ) {
					y[k][variable] = yWork[k];	// copy workspace to vector of solution
				}
			}
			else {
				for ( k = 0; k < nGridPoints; ++k ) { // turn vector
					y[k][variable] = yWork[nGridPoints-k-1];	// copy workspace to vector of solution
				}
			}
		}
	}

// G = U * a
	for ( k = 0; k < nGridPoints; ++k ) {
		y[k][fUVelocity] *= strainRateIn;
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
	
//	update properties
	UpdateSolution( yMat, yLeftVec, yRightVec );
	SetFlameNode( kPrev );
	ComputeProperties( fFlameNode, temp[kPrev], Y[kPrev], GetPressure() );
	SetFlameNode( nGridPoints );
	ComputeProperties( fFlameNode, temp[nGridPoints], Y[nGridPoints], GetPressure() );

	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fYLeftVec->vec[i] = yLeft[fFirstSpecies+i];
		fYRightVec->vec[i] = yRight[fFirstSpecies+i];
	}

	
	if ( fSoot ) {
		Double	W;
		for ( k = 0; k < nGridPoints; ++k ) {
			fProperties->ComputeMixtureMolarMass( W, Y[k]
					, fSpecies->GetMolarMass()->vec, nSpeciesInSystem );
			for ( j = 0; j < GetSoot()->GetNSootMoments(); ++j ) {
				y[k][GetSoot()->GetOffsetSootMoments()+j] *= RGAS * y[k][fTemperature]
						/ ( GetPressure() * W );
			}
		}
	}

	if ( inp->leftBoundary->fBcFlag[inp->fVVelocityOffset] == kNone ) {	// V + a rho U = 0
		int		mixtureSpecificationLeft = GetMixtureSpecificationLeft();
		Double	*rho = GetProperties()->GetDensity()->vec;
		// first set g
		if ( mixtureSpecificationLeft == kMassFlux ) {	// means tsuji bc
			yLeft[fUVelocity] = 0.0;
		}
		else {
			yLeft[fUVelocity] = sqrt( rho[nGridPoints] / rho[-1] );
		}
		
		// V + a rho U = 0
		yLeft[fVVelocity] = y[0][fVVelocity] + ( locX[0] - bt->GetLeft() ) 
							* GetStrainRate() * rho[-1] * yLeft[fUVelocity];
	}
	else {
		Double	*rho = GetProperties()->GetDensity()->vec;
#ifdef BCINFLUX
		fprintf( stderr, "input values for V are assumed to have units [kg/s m^2]\n" );
		if ( !fLiquidPoolBC ) {
			fprintf( stderr, "input values for VLeft  = %g [kg/s m^2]\n", yLeft[fVVelocity] );
		}
		fprintf( stderr, "input values for VRight = %g [kg/s m^2]\n", yRight[fVVelocity] );
#else
		fprintf( stderr, "input values for V are assumed to have units [m/s]\n" );
		yLeft[fVVelocity] = yLeft[fVVelocity] * rho[-1];
		yRight[fVVelocity] = yRight[fVVelocity] * rho[nGridPoints];
		Double	*bcLeft = grid->GetBcLeft()->vec;
		Double	*bcRight = grid->GetBcRight()->vec;
		bcLeft[fVVelocity] = bcLeft[fVVelocity] * rho[-1];
		bcRight[fVVelocity] = bcRight[fVVelocity] * rho[nGridPoints];
#endif
	}

	

	if ( GetSpecies()->IsConstantLewisNumber() ) {
		CompLewisNumbers( GetSpecies()->GetLewisNumberFile() );
	}
	DiffPhysEigenPostIter( this );

	FreeSpline( theSpline );
	delete yIn;
	delete xIn;

	adapGrid->SetSolutionScaler();
	
	fp = GetOutfile( "initialguess", TFlame::kData );
	bt->PrintSolution( locX, y, GetVariableNames(), fp );
	fclose(fp);
}

Double TCountDiffPhysEigen::MixFracDiffusion( int nVariable, NodeInfoPtr nodeInfo )
{
	Double	*lambda = fFlameNode->mixConductivity;
	Double	*cp = fFlameNode->mixHeatCapacity;
	Double	diffPlusHm = ( lambda[kCurr] / cp[kCurr]
					+ lambda[kNext] / cp[kNext] ) * nodeInfo->hm;
	Double	diffMinusH = ( lambda[kPrev] / cp[kPrev]
					+ lambda[kCurr] / cp[kCurr] ) * nodeInfo->h;
	Double	yPrev = nodeInfo->yPrev[nVariable];
	Double	y = nodeInfo->y[nVariable];
	Double	yNext = nodeInfo->yNext[nVariable];
	
	return ( diffPlusHm * ( yNext - y ) + diffMinusH * ( yPrev - y ) ) 
				/ nodeInfo->hnenn;
}

void TCountDiffPhysEigen::CompleteSpeciesDiffusion( NodeInfoPtr nodeInfo )
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
									/*/ hPlushm*/;
	}
}

Double TCountDiffPhysEigen::SpeciesDiffusion( int nVariable, int speciesIndex, NodeInfoPtr nodeInfo )
{
// returns    d/dy( rho D_k dY_k/dy)

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

void TCountDiffPhysEigen::FillJacMixFracDiffusion( int nVariable, int nEquation, NodeInfoPtr nodeInfo, Flag sign )
{
//	fill jacobian with  d/dy( lambda / cp * dZ/dy )

	Double	*lambda = fFlameNode->mixConductivity;
	Double	*cp = fFlameNode->mixHeatCapacity;
	Double	diffPlusHm = ( lambda[kCurr] / cp[kCurr]
					+ lambda[kNext] / cp[kNext] ) * nodeInfo->hm;
	Double	diffMinusH = ( lambda[kPrev] / cp[kPrev]
					+ lambda[kCurr] / cp[kCurr] ) * nodeInfo->h;

	if ( sign == kNegative ) {
		diffPlusHm *= -1.0;
		diffMinusH *= -1.0;
	}

	nodeInfo->a[nVariable][nEquation] -= ( diffPlusHm + diffMinusH );
	if ( !nodeInfo->lastPoint ) {
		nodeInfo->b[nVariable][nEquation] += diffPlusHm;
	}
	if ( !nodeInfo->firstPoint ) {
		nodeInfo->c[nVariable][nEquation] += diffMinusH;
	}
}

void TCountDiffPhysEigen::FillJacSpeciesDiffusion( int nVariable, int speciesIndex, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
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

void DiffPhysEigenOutput( void *object, FILE *fp, char* tail )
{
	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
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
	Double			*V = flame->GetV()->vec;
	Double			*G = flame->GetG()->vec;
	Double			*Z = flame->GetZ()->vec;
	Double			*P = flame->GetP()->vec;
	Double			**y = currentGrid->GetY()->mat;
	Double			*yLeft = currentGrid->GetYLeft()->vec,
					*yRight = currentGrid->GetYRight()->vec;
	int				i, k;
	int				gridPoints = currentGrid->GetNGridPoints();
	int				nOfSpecies = species->GetNOfSpecies();
	int				nOfSpeciesIn = species->GetNSpeciesInSystem();
	int				nOfVariables = bt->GetNVariables();
	int				nOfEquations = bt->GetNEquations();
	int				firstSpecies = flame->GetOffsetFirstSpecies();
	int				tempOffset = flame->GetOffsetTemperature();
	int				fUVelocity = flame->GetOffsetUVelocity();
	int				fVVelocity = flame->GetOffsetVVelocity();
	int				fMixFrac = flame->GetOffsetMixFrac();
	int				fPStrain = flame->GetOffsetPStrain();
	time_t			theDate;
	char			buffer[80];
	ConstStringArray	varNames = flame->GetVariableNames();
	char			**names = species->GetNames();
	VectorPtr 		etaVec = NewVector( gridPoints + 2 );
	Double			*eta = etaVec->vec;
	VectorPtr 		scalarDissVec = NewVector( gridPoints + 2 );
	Double			*scalarDiss = scalarDissVec->vec;
	Double			stoechScalarDiss;
	Double			strainRate = flame->GetStrainRate();
	Flag			fpOpen = FALSE;
	if ( !fp ) {
		fpOpen = TRUE;
		fp = flame->GetOutputFile( NULL, tail, TFlame::kNone );
	}

	flame->ScalarDissipation( bt, scalarDissVec, &stoechScalarDiss );
// write header
	fprintf( fp, "header\n\n" );

	fprintf( fp, "title = \"planar counterflow diffusion flame\"\n" );
	fprintf( fp, "mechanism = \"%s\"\n", flame->GetInputData()->fReactionFile );
	fprintf( fp, "author = \"%s\"\n", flame->GetAuthor() );
	if ( flame->GetGeometry() == 1 ) {
		fprintf( fp, "geometry = \"axisymmetric\"\n" );
	}
	time( &theDate );
	strcpy( buffer, ctime( &theDate ) );
	if ( buffer[strlen(buffer)-1] == '\n' )
		buffer[strlen(buffer)-1] = '\0';
	fprintf( fp, "date = \"%s\"\n\n", buffer );
	for ( i = 0; i < flame->GetNFuels(); ++i ) {
		fprintf( fp, "fuel = \"%s\"\n", varNames[firstSpecies+flame->GetFuelIndex( i )] );
	}

	fprintf( fp, "pressure = %g [bar]\n", flame->GetPressure() / 1.0e5 );
	fprintf( fp, "strainrate = %g [1/s]\n", strainRate );
	fprintf( fp, "strainrateSesha = %g [1/s]\n", flame->GetSeshaStrainRate() );
	fprintf( fp, "strainrateLiquidPool = %g [1/s]\n", flame->GetLiquidPoolStrainRate() );
	fprintf( fp, "VLeft = %g [kg/m^2s]\n", yLeft[fVVelocity] );
	fprintf( fp, "VRight = %g [kg/m^2s]\n", yRight[fVVelocity] ); 
	fprintf( fp, "vLeft = %g [m/s]\n", yLeft[fVVelocity]/rho[-1] );
	fprintf( fp, "vRight = %g [m/s]\n", yRight[fVVelocity]/rho[gridPoints] ); 

	if ( species->IsConstantLewisNumber() ) {
		fprintf( fp, "ConstantLewisNumbers = \"True\"\n" );
	}
	if ( flame->fLiquidPoolBC ) {
		fprintf( fp, "LiquidPool = \"True\"\n" );
	}
	
	fprintf( fp, "Tmax = %g [K]\n", temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
	fprintf( fp, "StScalarDissRate = %g [1/s]\n", stoechScalarDiss );
	
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

//  write bc
	Double	locMoleMass;
	fprintf( fp, "\nOxidizerSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", yRight[tempOffset] );
	for ( i = 0; i < nOfSpeciesIn; ++i ) { // write Y_i
		if ( fabs( massFracs[gridPoints][i] ) > 1.0e-5 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], massFracs[gridPoints][i] );
		}
	}
	for ( i = 0; i < nOfSpeciesIn; ++i ) { // write X_i
		locMoleMass = massFracs[gridPoints][i] * mixMolarMass[gridPoints] / molarMass[i];
		if ( fabs( locMoleMass ) > 1.0e-5 ) {
			fprintf( fp, "\tMolefraction-%s = %g\n", names[i], locMoleMass );
		}
	}

	if ( flame->GetMixtureSpecificationRight() == kMassFlux ) {
		Double	*bcRight = currentGrid->GetBcRight()->vec;
		int		*bcFlagRight = currentGrid->GetBcFlagRight();
		for ( i = 0; i < nOfSpeciesIn; ++i ) { // write epsilon_i
			if ( bcFlagRight[i+firstSpecies] == kMassFlux && fabs( bcRight[i+firstSpecies] ) > 1.0e-5 ) {
				fprintf( fp, "\tMassflux-%s = %g\n", names[i], bcRight[i+firstSpecies] );
			}
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "FuelSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", yLeft[tempOffset] );
	for ( i = 0; i < nOfSpeciesIn; ++i ) {
		if ( fabs( massFracs[kPrev][i] ) > 1.0e-5 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], massFracs[kPrev][i] );
		}
	}
	for ( i = 0; i < nOfSpeciesIn; ++i ) { // write X_i
		locMoleMass = massFracs[kPrev][i] * mixMolarMass[kPrev] / molarMass[i];
		if ( fabs( locMoleMass ) > 1.0e-5 ) {
			fprintf( fp, "\tMolefraction-%s = %g\n", names[i], locMoleMass );
		}
	}
	if ( flame->GetMixtureSpecificationLeft() == kMassFlux ) {
		Double	*bcLeft = currentGrid->GetBcLeft()->vec;
		int		*bcFlagLeft = currentGrid->GetBcFlagLeft();
		for ( i = 0; i < nOfSpeciesIn; ++i ) { // write epsilon_i
			if ( bcFlagLeft[i+firstSpecies] == kMassFlux && fabs( bcLeft[i+firstSpecies] ) > 1.0e-5 ) {
				fprintf( fp, "\tMassflux-%s = %g\n", names[i], bcLeft[i+firstSpecies] );
			}
		}
	}
	fprintf( fp, "end\n\n" );

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
	
// write eta
	
	flame->XToEta( bt, etaVec );
	flame->OriginToZstoich( etaVec, flame->GetZ(), flame->fMassFraction->GetZStoe() );
//	flame->OriginToZstoich( etaVec, flame->GetZ(), flame->GetZStoich() );
	flame->PrintFlameletVector( gridPoints+2, eta, "eta", fp );
		
			
// write solution
	// write V-Velocity, U-Velocity, mixture fraction and temperature
	flame->PrintFlameletVector( gridPoints+2, &V[kPrev], "V [kg/(s*m*m)]", fp );
	flame->PrintFlameletVector( gridPoints+2, &G[kPrev], "G [1/s]", fp );
	flame->PrintFlameletVector( gridPoints+2, &Z[kPrev], "Z", fp );
	flame->PrintFlameletVector( gridPoints+2, &P[kPrev], "P [kg/(m^3*s^2)]", fp );
	flame->PrintFlameletVector( gridPoints+2, &temp[kPrev], "temperature [K]", fp );

	// write massfractions of species
	char	speciesName[127];
	for ( i = 0; i < nOfSpecies; ++i ) {
		sprintf( speciesName, "massfraction-%s", names[i] );
		flame->PrintFlameletVector( gridPoints+2, &massFracs[kPrev][i], speciesName
				, fp, flame->GetMassFracs()->phys_rows );
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

//	write f
	Double	coeff = -sqrt( flame->GetStrainRate() * flame->fFlameNode->rhoInf * flame->fFlameNode->viscosityInf );
	fprintf( fp, "f\n" );
		fprintf( fp, "\t%-.6e", V[kPrev] / MAX(1e-30,coeff) );
		for ( k = 0; k < gridPoints; ++k ) {
			fprintf( fp, "\t%-.6e", V[k] / MAX(1e-30,coeff) );
			if ( (k+2) % 5 == 0 ) {
				fprintf( fp, "\n" );
			}
		}
		fprintf( fp, "\t%-.6e\n", V[gridPoints] / MAX(1e-30,coeff) );
	
//	write U
	fprintf( fp, "U\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", G[k-1] / MAX(1e-30,strainRate) );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( (k+1) % 5 ) {
		fprintf( fp, "\n" );
	}
	
//	write v
	fprintf( fp, "velocity [m/s]\n" );
		fprintf( fp, "\t%-.6e", V[kPrev] / rho[kPrev] );
		for ( k = 0; k < gridPoints; ++k ) {
			fprintf( fp, "\t%-.6e", V[k] / rho[k] );
			if ( (k+2) % 5 == 0 ) {
				fprintf( fp, "\n" );
			}
		}
	fprintf( fp, "\t%-.6e\n", V[gridPoints] / rho[gridPoints] );

//	write heat release
	fprintf( fp, "HeatRelease [J/m^3 s]\n" );
	Double	heatrel = 0.0;
	Double	**prod = flame->GetSpecies()->GetProductionRate()->mat;
	Double	**ent = flame->GetSpecies()->GetEnthalpy()->mat;

	fprintf( fp, "\t%-.6e", 0.0 );

	for ( k = 0; k < gridPoints; ++k ) {
		heatrel = 0.0;
		for ( i = 0; i < nOfSpeciesIn; ++i ) {
			heatrel += ent[k][i] * prod[k][i];
		}
		fprintf( fp, "\t%-.6e", -heatrel );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", 0.0 );

//	write total enthalpy
	fprintf( fp, "TotEnt [J/m^3]\n" );

	for ( k = -1; k < gridPoints+1; ++k ) {
		heatrel = 0.0;
		for ( i = 0; i < nOfSpeciesIn; ++i ) {
			heatrel += ent[k][i] * massFracs[k][i];
		}
		fprintf( fp, "\t%-.6e", heatrel );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( (k+1) % 5 ) {
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

//	write scalar dissipation rate
	fprintf( fp, "chi\n" );
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
	fprintf( fp, "\t%-.6e", flame->ComputeZBilger( massFracs[kPrev], massFracs[kPrev], massFracs[gridPoints] ) );
	for ( k = 0; k < gridPoints; ++k ) {
		fprintf( fp, "\t%-.6e", flame->ComputeZBilger( massFracs[k], massFracs[kPrev], massFracs[gridPoints] ) );
		if ( (k+2) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	fprintf( fp, "\t%-.6e\n", flame->ComputeZBilger( massFracs[gridPoints], massFracs[kPrev], massFracs[gridPoints] ) );
		
	
	fprintf( fp, "trailer\n" );
	if ( flame->fSoot ) {
		if ( flame->fSoot->GetCoagFact() < 1.0 ) {
			fprintf( fp, "CoagFact =  %g\n", flame->GetSoot()->GetCoagFact() );
		}
	}
	if ( nOfEquations < nOfVariables) {
		fprintf( fp, "number of converged equations is %d\n", nOfEquations );
	}
	
	DisposeVector( scalarDissVec );
	DisposeVector( etaVec );
	if ( fpOpen ) {
		fclose( fp );
	}
}

int TCountDiffPhysEigen::ScalarDissipation( TNewtonPtr bt, VectorPtr scalarDissVec, Double *stoechScalarDiss )
{
// return first gridpoint larger than x_stoich 
	int				k;
	int				k_Stoich;
	int				nOfGridPoints = bt->GetCurrentGridPoints();
	TGridPtr		grid = bt->GetGrid()->GetCurrentGrid();
	NodeInfoPtr		nodeInfo = bt->GetNodeInfo();
	Double			left = bt->GetLeft();
	Double			right = bt->GetRight();
	Double			*x = grid->GetX()->vec;
	Double			*Z = fSolMixFrac->vec;
	Double			*scalarDiss = scalarDissVec->vec;
	Double			*lambda = fProperties->GetConductivity()->vec;
	Double			*density = fProperties->GetDensity()->vec;
	Double			*cp = fProperties->GetHeatCapacity()->vec;
	Double			dZdx;
	Double			zStoe = fMassFraction->GetZStoe();
	
	dZdx          = ( Z[0] - Z[kPrev] ) / ( x[0] - left );
	scalarDiss[0] = 2 * lambda[0] / (density[0]  * cp[0]) * dZdx * dZdx;
	for ( k = 0; k < nOfGridPoints; ++k ) {
		bt->SetNodeInfo( this, k );
		dZdx = FirstDeriv( Z[k-1], Z[k], Z[k+1], nodeInfo->hm, nodeInfo->h );
		scalarDiss[k+1] = 2 * lambda[k] / (density[k]  * cp[k]) * dZdx * dZdx;
		if ( ( Z[k] - zStoe ) * ( Z[k-1] - zStoe ) <= 0.0 ) {
			*stoechScalarDiss = scalarDiss[k] + ( scalarDiss[k+1] - scalarDiss[k] ) 
							/ ( Z[k] - Z[k-1] ) * ( zStoe - Z[k-1] );
			k_Stoich = k;
		}
	}
	dZdx = ( Z[nOfGridPoints] - Z[nOfGridPoints-1] ) / ( right - x[nOfGridPoints-1] );
	scalarDiss[nOfGridPoints+1] = 2 * lambda[nOfGridPoints-1] / (density[nOfGridPoints-1] * cp[nOfGridPoints-1]) * dZdx * dZdx;

	return k_Stoich;
}

void SetDiffPhysEigenNodeInfo( int k, void *object )
{
	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	
	flame->SetFlameNode( k );
}

void DiffPhysEigenPostConv( void *object )
{
	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	TNewtonPtr				bt = flame->GetSolver()->bt;
	TAdaptiveGridPtr		grid = bt->GetGrid();
	int						isConverged = bt->GetConvergeNewton();
	Flag					BCConverged = TRUE;
	Double					*temp = flame->GetTemperature()->vec;
	Double					VNow = grid->GetCurrentGrid()->GetBcLeft()->vec[flame->GetOffsetVVelocity()];
	
	if ( isConverged ) {
		flame->SaveSolution();
		if ( flame->fPrintRHSSpecies ) {
			flame->PrintRHSSpecies( flame->GetSolver()->bt );
		}
		if ( flame->fPrintRHSTemp ) {
			flame->PrintRHSTemp( flame->GetSolver()->bt );
		}
		if ( flame->fSensAnal ) {
			flame->SensitivityAnalysis( 1.0, 1.0, kPhysical );
		}
		if ( flame->fReactionFluxes ) {
			flame->ReactionFluxes( kPhysical );
			flame->GetReaction()->PrintReactionRates( flame );
			flame->fReaction->PrintRateCoeffs( flame );
			flame->fReaction->PrintDetailedHeatRelease( flame );
		}
		for ( int i = 0; i < flame->fNSensObj; ++i ) {
			if ( flame->fSpecies->FindSpecies( flame->fSensObj[i] ) >= 0 ) {
				flame->GetSpecies()->PrintProdRateTerms( flame->fSensObj[i], flame );
			}
		}
		if ( flame->fLiquidPoolBC ) {
			BCConverged = flame->UpdateBCLiquidPool( object );
		}

#ifdef CHECKIT
		BCConverged = BCConverged;
		fprintf( stderr, "Tmax = %g @ a_%s = %g 1/s\n"
			, temp[LocationOfMax( grid->GetCurrentGrid()->GetNGridPoints()+2, &temp[kPrev] ) - 1]
			, ( flame->fLiquidPoolBC ) ? "Liq" : "Sesha"
			, ( flame->fLiquidPoolBC ) ? flame->GetLiquidPoolStrainRate() : flame->GetSeshaStrainRate() );

		if ( !flame->fLiquidPoolBC ) {
			fprintf( stderr, "VLiq would be %g\n", flame->GetVLiquidPool() );
		}
		flame->PostConvergence( object );
		DiffPhysEigenPostIter( flame );
#else
		static int				liquidTries = 0;
		if ( BCConverged ) {
			liquidTries = 0;
			if ( flame->fLiquidPoolBC ) {
				fprintf( stderr, "Liquid pool boundary condition converged\n" );
			}
			else {
				fprintf( stderr, "Tmax = %g @ a_Liq = %g 1/s and a_Sesh = %g 1/s\n"
					, temp[LocationOfMax( grid->GetCurrentGrid()->GetNGridPoints()+2, &temp[kPrev] ) - 1]
					, flame->GetLiquidPoolStrainRate(), flame->GetSeshaStrainRate() );
				fprintf( stderr, "VLiq would be %g\n", flame->GetVLiquidPool() );
			}
			flame->PostConvergence( object );
			DiffPhysEigenPostIter( flame );
		}
		else {
			++liquidTries;
			fprintf( stderr, "liquidTries now %d\n", liquidTries );
			if ( liquidTries == 5 ) {
				fprintf( stderr, "LiquidPool condition does not converge\n" );
			}
			else {
				flame->fSolver->ReInit();
			}
		}
#endif
	}
	else {
		flame->RestoreSolution();
		flame->PostConvergence( object );
		DiffPhysEigenPostIter( flame );
	}
}

ConstStringArray GetDiffPhysEigenVarNames( void *object )
{
	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	
	return flame->GetVariableNames();
}

FILE *TCountDiffPhysEigen::GetOutputFile( char *head, char *tail, FileType type )
{
	int				fuelIndex = GetFuelIndex();
	char			*name = new char[64];
	FILE			*fp;
	char			**speciesNames = fSpecies->GetNames();
	int				tOxidizer = ( int ) fSolTemp->vec[fSolTemp->len];
	int				tFuel = ( int ) fSolTemp->vec[kPrev];
		
	sprintf( name, "%s%s%.8s_p%.2da%.5dtf%.4dto%.4d%s"
					, ( head ) ? head : "", ( head ) ? "_" : ""
					, speciesNames[fuelIndex]
					, ( int ) floor( GetPressure() * 1.0e-5 + 0.5 )	// in [bar]
					, ( int ) floor( ( ( fLiquidPoolBC ) 
								? GetLiquidPoolStrainRate() : GetSeshaStrainRate() ) + 0.5 )			// in [1/s]
					, ( int )( tFuel )							// in [K]
					, ( int )( tOxidizer ) 						// in [K]
					, ( tail ) ? tail : "" );

	fp = GetOutfile( name, type );
	delete name;

	return fp;
}

Double TCountDiffPhysEigen::GetStrainRate( void )
{
	if ( fSolP->vec[0] < 0.0 ) {
		fprintf( stderr, "warning: eigenvalue negative\n" );
		return 0.0;
	}
	return sqrt( fSolP->vec[0] / fFlameNode->rhoInf );
}

Double TCountDiffPhysEigen::GetSeshaStrainRate( void )
{
	TNewtonPtr		bt = GetSolver()->bt;
	TGridPtr		currentGrid = bt->GetGrid()->GetCurrentGrid();
	Double			*yLeft = currentGrid->GetYLeft()->vec,
					*yRight = currentGrid->GetYRight()->vec;
	Double			*rho = fProperties->GetDensity()->vec;
	int				gridPoints = currentGrid->GetNGridPoints();
	Double			rhoOx = rho[gridPoints];
	Double			rhoFu = rho[-1];
	Double			vOx = fabs( yRight[fVVelocity] / rhoOx );
	Double			vFu = fabs( yLeft[fVVelocity] / rhoFu );

	return 2.0 * vOx * ( 1.0 + vFu / vOx * sqrt( rhoFu/rhoOx ) ) 
						/ ( bt->GetRight() - bt->GetLeft() );
}

Double TCountDiffPhysEigen::GetLiquidPoolStrainRate( void )
{
	TNewtonPtr		bt = GetSolver()->bt;
	TGridPtr		currentGrid = bt->GetGrid()->GetCurrentGrid();
	Double			*yLeft = currentGrid->GetYLeft()->vec,
					*yRight = currentGrid->GetYRight()->vec;
	Double			*rho = fProperties->GetDensity()->vec;
	int				gridPoints = currentGrid->GetNGridPoints();
	Double			rhoOx = rho[gridPoints];
	Double			rhoFu = rho[-1];
	Double			vOx = fabs( yRight[fVVelocity] / rhoOx );
	Double			vFu = fabs( yLeft[fVVelocity] / rhoFu );

	return 2.0 * vOx / ( bt->GetRight() - bt->GetLeft() );
}

void TCountDiffPhysEigen::FillJacNonlinearConvectUpwind( int nVariable1, int nVariable2, NodeInfoPtr nodeInfo, Double coeff, Flag velocityPositive )
{
//	velocityPositive should have the value FALSE, if 'nVariable1' has the negative 
//	direction of the physical velocity

// fills the jacobian with     coeff * Y1 * dY2/dy
// it is assumed that 'nVariable1' is the index of variable Y1 
// and that Y2 and the current equation have the index 'nVariable2'

	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*yPrev = nodeInfo->yPrev;

	if ( y[nVariable1] * yNext[nVariable1] <= 0.0 || y[nVariable1] * yPrev[nVariable1] <= 0.0 ) {
		FillJacNonlinearConvectCentral( nVariable1, nVariable2, nodeInfo );
	}
	else {
		if ( ( y[nVariable1] > 0.0 && velocityPositive ) || ( y[nVariable1] < 0.0 && !velocityPositive ) ) {
			coeff *= h * ( h + hm );
			nodeInfo->a[nVariable1][nVariable2] += coeff * ( y[nVariable2] - yPrev[nVariable2] );
			nodeInfo->a[nVariable2][nVariable2] += coeff * y[nVariable1];
			if ( !nodeInfo->firstPoint ) {
				nodeInfo->c[nVariable2][nVariable2] -= coeff * y[nVariable1];
			}
		}
		else {
			coeff *= hm * ( h + hm );
			nodeInfo->a[nVariable1][nVariable2] += coeff * ( yNext[nVariable2] - y[nVariable2] );
			nodeInfo->a[nVariable2][nVariable2] -= coeff * y[nVariable1];
			if ( !nodeInfo->lastPoint ) {
				nodeInfo->b[nVariable2][nVariable2] += coeff * y[nVariable1];
			}
		}
	}
}

Double TCountDiffPhysEigen::NonlinearConvectUpwind( int nVariable1, int nVariable2, NodeInfoPtr nodeInfo, Flag velocityPositive )
{
//	velocityPositive should have the value FALSE, if 'y1' has the negative 
//	direction of the physical velocity

// returns     Y1 * dY2/dy

	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*yPrev = nodeInfo->yPrev;


	if ( ( y[nVariable1] > 0.0 && velocityPositive ) || ( y[nVariable1] < 0.0 && !velocityPositive ) ) {
		return ( y[nVariable1] * ( y[nVariable2] - yPrev[nVariable2] ) / hm );
	}
	else {
		return ( y[nVariable1] * ( yNext[nVariable2] - y[nVariable2] ) / h );
	}
}

void TCountDiffPhysEigen::PrintRHSTemp( TNewtonPtr bt )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		k;
    int         		N = currentGrid->GetNGridPoints();
	FILE				*fp = NULL;
	
	UpdateThermoProps();
	
	fp = GetOutfile( "tempeq", TFlame::kData );

	fprintf( fp, "*\n%-12s\t%-12s\t%-12s\t%-12s\t%-12s\t%-12s", "x", "convection", "diffusion", "entFlux", "entFluxDiCo", "production" );
	if ( fProperties->GetRadiation() ) {
		fprintf( fp, "\t%-12s", "radiation" );
		if ( GetSoot() ) {
			fprintf( fp, "\t%-12s", "sootrad" );
		}
	}
	fprintf( fp, "\n" );

	for ( k = 0; k < N; ++k ){
		bt->SetNodeInfo( this, k );
		PrintRHSTemp( nodeInfo, fp );
	}
    fclose( fp );
}

void TCountDiffPhysEigen::PrintRHSTemp( NodeInfoPtr nodeInfo, FILE *fp )
{
	int		eqLoop;
	Double	sumCpDdYdx;
	Double	sumMH;
	Double	sumCpY = 0.0;
	Double	oneOverCp = 1.0 / fFlameNode->mixHeatCapacity[kCurr];
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*heatCapacity = fFlameNode->heatCapacity;
	Double	*enthalpy = fFlameNode->enthalpy;
	Double	*productionRate = fFlameNode->productionRate;
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	*YPrev = fFlameNode->Y[kPrev];
	Double	*Y = fFlameNode->Y[kCurr];
	Double	*YNext = fFlameNode->Y[kNext];
	Double	*temp = fFlameNode->temp;
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;

	fprintf( fp, "%-.6e", *nodeInfo->x );
	fprintf( fp, "\t%-.6e", NonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo ) );

	sumCpDdYdx = 0.0;
	sumMH = 0.0;

	fprintf( fp, "\t%-.6e", oneOverCp * SecondDerivDiffusion( fTemperature, fFlameNode->mixConductivity, nodeInfo ) );

// compute all sums
	for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
		sumCpDdYdx += heatCapacity[eqLoop] * diffusivity[eqLoop] 
					* FirstDeriv( YPrev[eqLoop], Y[eqLoop], YNext[eqLoop], hm, h );
#ifdef DIFFUSIVITYCORRECTION
		sumCpY += heatCapacity[eqLoop] * Y[eqLoop];
#endif
	}

	for ( eqLoop = 0; eqLoop < nSpeciesInSystem; ++eqLoop ) {
		sumMH += productionRate[eqLoop] * enthalpy[eqLoop];
	}

	fprintf( fp, "\t%-.6e", -oneOverCp * mixDensity[kCurr] * sumCpDdYdx 
						* FirstDeriv( temp[kPrev], temp[kCurr], temp[kNext], hm, h ) );
	fprintf( fp, "\t%-.6e", oneOverCp * mixDensity[kCurr] * sumCpY
						* FirstDeriv( temp[kPrev], temp[kCurr], temp[kNext], hm, h )
						* fFlameNode->diffCorr[kCurr] );
	fprintf( fp, "\t%-.6e", oneOverCp * sumMH );
	
	if ( fProperties->GetRadiation() ) {
		fprintf( fp, "\t%-.6e", -oneOverCp * fFlameNode->radiation[kCurr] );
		if ( GetSoot() && GetSoot()->WithSootRadiation() ) {
			fprintf( fp, "\t%-.6e", oneOverCp * fSoot->GetSootRadiation( temp[kCurr], fFlameNode->moments ) );
		}
	}
	
	fprintf( fp, "\n" );
}

void TCountDiffPhysEigen::PrintRHSSpecies( TNewtonPtr bt )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		i, k;
    int         		N = currentGrid->GetNGridPoints();
	char				**names = GetSpecies()->GetNames();
	int					nSpecIn = fSpecies->GetNSpeciesInSystem();
	FILE				*fp = NULL;
	int					num = 5;
	int					start = 0, end;
	int					fileCount = 1;
		
	do {
		end = ( int )MIN( start + 255 / num, nSpecIn );
		sprintf( GetOutFileBuff(), "%sspeciesRHS%d.dout", GetOutputPath(), fileCount++ );
		if ( !( fp = fopen( GetOutFileBuff(), "w") ) ) { 
			cerr << "#warning: unable to open file " << GetOutFileBuff() << NEWL;
			exit(2);
		}
		
		fprintf( fp, "*\n%-12s", "y" );
		for ( i = start; i < end; ++i ) {
#ifdef DIFFUSIVITYCORRECTION
			fprintf( fp, "\tConv_%-7s\tDiff_%-7s\tDiCo_%-7s\tProd_%-7s\tCons_%-7s", names[i], names[i], names[i], names[i], names[i] );
#else
			if ( fThermoDiffusion ) {
				fprintf( fp, "\tConv_%-7s\tDiff_%-7s\tDifT_%-7s\tProd_%-7s\tCons_%-7s", names[i], names[i], names[i], names[i], names[i] );
			}
			else {
				fprintf( fp, "\tConv_%-7s\tDiff_%-7s\tProd_%-7s\tCons_%-7s", names[i], names[i], names[i], names[i] );
			}
#endif
		}
		fprintf( fp, "\n" );
			
		for ( k = 0; k < N; ++k ){
			bt->SetNodeInfo( this, k );
			PrintRHSSpecies( start, end, nodeInfo, fp );
		}
		fclose( fp );
		start = end;
	} while ( end < nSpecIn );
}

void TCountDiffPhysEigen::PrintRHSSpecies( int start, int end, NodeInfoPtr nodeInfo, FILE *fp )
{
	int		eqLoop, speciesEq;
	int		nSpeciesInSystem = GetSpecies()->GetNSpeciesInSystem();
#ifdef DIFFUSIVITYCORRECTION
	Double	diffCorr;
#endif
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

		fprintf( fp, "\t%-.6e", NonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo ) );

		diffTerm = -SpeciesDiffusion( speciesEq, eqLoop, nodeInfo );
#ifdef DIFFUSIVITYCORRECTION
		diffCorr = DiffCorr( speciesEq, nodeInfo );
		fprintf( fp, "\t%-.6e\t%-.6e", diffTerm, diffCorr );
#else
		fprintf( fp, "\t%-.6e", diffTerm );
		if ( fThermoDiffusion ) {
			fprintf( fp, "\t%-.6e", -ThermoDiffusion( eqLoop, kPhysical, nodeInfo ) );
		}
#endif

		int j;
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
		fprintf( fp, "\t%-.6e\t%-.6e", source, -sink );
	}
	fprintf( fp, "\n" );
}

void DiffPhysEigenUpdateLeftNewton( void  */*object*/ )
{}

void DiffPhysEigenUpdateLeftBoundary( void  *object )
{
	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	int			fUVelocity = flame->GetOffsetUVelocity();
	int			fVVelocity = flame->GetOffsetVVelocity();
	int			fTemperature = flame->GetOffsetTemperature();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();
	int			fMixFrac = flame->GetOffsetMixFrac();
	int			fPStrain = flame->GetOffsetPStrain();
	TNewtonPtr	bt = flame->GetSolver()->bt;
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = currGrid->GetNGridPoints();
	int			*bcFlagLeft = currGrid->GetBcFlagLeft();
	int			mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
	MatrixPtr	yMat = currGrid->GetY();
	VectorPtr	yLeftVec = currGrid->GetYLeft();
	VectorPtr	yRightVec = currGrid->GetYRight();
	Double		**y = yMat->mat;
	Double		*x = currGrid->GetX()->vec;
	Double		*yLeft = yLeftVec->vec;
	Double		*yFirst = y[0];
	Double		*bcLeft = currGrid->GetBcLeft()->vec;
	Double		hFirst = x[0] - bt->GetLeft();
	T1DSpeciesPtr	species = flame->GetSpecies();
	T1DPropertiesPtr	properties = flame->GetProperties();
	Double		*diffusivityLeft = species->GetDiffusivity()->mat[-1];
	Double		*thermDiffLeft = species->GetDiffTherm()->mat[kPrev];
	Double		*mixViscosity = properties->GetViscosity()->vec;
	Double		&lambdaLeft = properties->GetConductivity()->vec[-1];
	Double		&cpLeft = properties->GetHeatCapacity()->vec[-1];
	Double		&rhoLeft = properties->GetDensity()->vec[-1];
	Double		pressure = flame->GetPressure();
	Double		*temp = flame->GetTemperature()->vec;
	Double		**Y = flame->GetMassFracs()->mat;
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();

	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, temp[kPrev], Y[kPrev], pressure );
	flame->SetFlameNode( 0 );
	flame->ComputeProperties( flame->fFlameNode, temp[0], Y[0], pressure );

	if ( mixtureSpecificationLeft == kMassFlux ) {	// means tsuji bc
		Double 		coeff;
		yLeft[fUVelocity] = 0.0;

#ifdef CHECKIT
		if ( flame->fLiquidPoolBC ) {
			flame->UpdateBCLiquidPool( object );
		}
#endif

		flame->CalcYLeft();

		coeff = lambdaLeft / ( cpLeft * yLeft[fVVelocity] * hFirst );
		yLeft[fMixFrac] = ( 1.0 + coeff * yFirst[fMixFrac] )
						/ ( 1.0 + coeff );

		yLeft[fPStrain] = yFirst[fPStrain];
	}
	else {
		cerr << "#error: boundary conditions for species at fuel side have" 
			<< " to be of kind massflux; " << NEWL
			<< "actually it is of kind " << mixtureSpecificationLeft << NEWL;
		exit( 2 );
	}
	
	flame->SetStrainRate( flame->GetStrainRate() ); // set strainrate in T1DFlame to actual value
	
	
	if ( flame->GetSoot() ) {
		flame->GetSoot()->PostIter( flame );
	}

	flame->UpdateSolutionOnePoint( yLeft, kPrev );
	flame->SetFlameNode( kPrev );

	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], pressure );
}

void DiffPhysEigenUpdateRightBoundary( void  *object )
{
	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	int			fUVelocity = flame->GetOffsetUVelocity();
	int			fVVelocity = flame->GetOffsetVVelocity();
	int			fTemperature = flame->GetOffsetTemperature();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();
	int			fMixFrac = flame->GetOffsetMixFrac();
	int			fPStrain = flame->GetOffsetPStrain();
	TNewtonPtr	bt = flame->GetSolver()->bt;
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = currGrid->GetNGridPoints();
	int			*bcFlagRight = currGrid->GetBcFlagRight();
	int			mixtureSpecificationRight = flame->GetMixtureSpecificationRight();
	MatrixPtr	yMat = currGrid->GetY();
	VectorPtr	yLeftVec = currGrid->GetYLeft();
	VectorPtr	yRightVec = currGrid->GetYRight();
	Double		**y = yMat->mat;
	Double		*x = currGrid->GetX()->vec;
	Double		*yRight = yRightVec->vec;
	Double		*yLast = y[nGridPoints-1];
	Double		hLast = bt->GetRight() - x[nGridPoints-1];
	Double		*bcRight = currGrid->GetBcRight()->vec;
	Double		hFirst = x[0] - bt->GetLeft();
	T1DSpeciesPtr	species = flame->GetSpecies();
	T1DPropertiesPtr	properties = flame->GetProperties();
	Double		*diffusivityRight = species->GetDiffusivity()->mat[nGridPoints];
	Double		*thermDiffRight = species->GetDiffTherm()->mat[nGridPoints];
	Double		*mixViscosity = properties->GetViscosity()->vec;
	Double		&lambdaRight = properties->GetConductivity()->vec[nGridPoints];
	Double		&cpRight = properties->GetHeatCapacity()->vec[nGridPoints];
	Double		&rhoRight = properties->GetDensity()->vec[nGridPoints];
	Double		pressure = flame->GetPressure();
	Double		*temp = flame->GetTemperature()->vec;
	Double		**Y = flame->GetMassFracs()->mat;
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();

	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	flame->SetFlameNode( nGridPoints-1 );
	flame->ComputeProperties( flame->fFlameNode, temp[nGridPoints-1], Y[nGridPoints-1], pressure );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, temp[nGridPoints], Y[nGridPoints], pressure );

	if ( mixtureSpecificationRight == kMassFlux ) {	// means tsuji bc
		Double 		coeff;
		yRight[fUVelocity] = 0.0;

		flame->CalcYRight();
		
/* ATTENTION HP */
		coeff = lambdaRight / ( cpRight * yRight[fVVelocity] * hLast );
		yRight[fMixFrac] = - coeff * yLast[fMixFrac]
						/ ( 1.0 - coeff );

		yRight[fPStrain] = yLast[fPStrain];
	}
	else {
		cerr << "#error: boundary conditions for species at oxidizer side have" 
			<< " to be of kind massflux" << NEWL;
		exit( 2 );
	}
	
	flame->SetStrainRate( flame->GetStrainRate() ); // set strainrate in T1DFlame to actual value
	

	if ( flame->GetSoot() ) {
		flame->GetSoot()->PostIter( flame );
	}

	flame->UpdateSolutionOnePoint( yRight, nGridPoints );

	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], pressure );
}

int DiffPhysEigenPostIter( void *object )
{
	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	int			fTemperature = flame->GetOffsetTemperature();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();
	int			i;
	TNewtonPtr	bt = flame->GetSolver()->bt;
	NodeInfoPtr nodeInfo = bt->GetNodeInfo();
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = currGrid->GetNGridPoints();
	MatrixPtr	yMat = currGrid->GetY();
	VectorPtr	yLeftVec = currGrid->GetYLeft();
	VectorPtr	yRightVec = currGrid->GetYRight();
	Double		**y = yMat->mat;

	Double		pressure = flame->GetPressure();
	Double		*temp = flame->GetTemperature()->vec;
	Double		**Y = flame->GetMassFracs()->mat;
	int			nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	

	for ( i = 0; i < nGridPoints; ++i ) {
		if ( flame->CheckSolution( y[i][fTemperature], &y[i][fFirstSpecies], nSpeciesInSystem ) ) {
			return 1;
		}
	}

//	update properties
	flame->UpdateSolution( yMat, yLeftVec, yRightVec );

	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, temp[kPrev], Y[kPrev], pressure );
	flame->SetFlameNode( 0 );
	flame->ComputeProperties( flame->fFlameNode, temp[0], Y[0], pressure );
	flame->SetFlameNode( nGridPoints-1 );
	flame->ComputeProperties( flame->fFlameNode, temp[nGridPoints-1], Y[nGridPoints-1], pressure );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, temp[nGridPoints], Y[nGridPoints], pressure );

	DiffPhysEigenUpdateLeftBoundary( flame );

	DiffPhysEigenUpdateRightBoundary( flame );
	
	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	flame->UpdateThermoProps();

	// find fNLeftofStag
	int	&point = flame->fNLeftofStag;
	Double	*V = flame->fSolV->vec;
	for ( point = 1; point < nGridPoints-1; ++point ) { // start at 1 because stag point should not be at first gridpoint
		if ( V[point+1] * V[point] <= 0.0 ) {
			break;
		}
	}
	
	return 0;
}

void TCountDiffPhysEigen::UpdateSolutionOnePoint( Double *y, int gridPoint )
{
	T1DFlame::UpdateSolution( y, gridPoint );

	fSolV->vec[gridPoint] = y[fVVelocity];
	fSolG->vec[gridPoint] = y[fUVelocity];
	fSolMixFrac->vec[gridPoint] = y[fMixFrac];
	fSolP->vec[gridPoint] = y[fPStrain];
}

Double TCountDiffPhysEigen::GetVLiquidPool( void )
{
	Double				lambdaLeft =  0.5 * ( fProperties->GetConductivity()->vec[-1] 
									+ fProperties->GetConductivity()->vec[0] );
	TGridPtr 			currGrid = fSolver->bt->GetGrid()->GetCurrentGrid();
	Double				tempFirst = currGrid->GetY()->mat[0][fTemperature];
	VectorPtr			yLeftVec = currGrid->GetYLeft();
	Double				*x = currGrid->GetX()->vec;
	Double				*yLeft = yLeftVec->vec;
	Double				hFirst = x[0] - fSolver->bt->GetLeft();
	Double				hL = 360064.8;		// latent heat of vaporization of the fuel at T = 373 K ( boiling point )

	return lambdaLeft * ( tempFirst - yLeft[fTemperature] ) / ( hFirst * hL );
}

Flag TCountDiffPhysEigen::UpdateBCLiquidPool( void *object )
{

	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	TNewtonPtr	bt = flame->GetSolver()->bt;
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	VectorPtr	yLeftVec = currGrid->GetYLeft();
	Double		*x = currGrid->GetX()->vec;
	Double		*yLeft = yLeftVec->vec;
	Double		*bcLeft = currGrid->GetBcLeft()->vec;
	Double		hFirst = x[0] - bt->GetLeft();
	T1DPropertiesPtr	properties = flame->GetProperties();
	Double		lambdaLeft =  0.5 * ( properties->GetConductivity()->vec[-1] 
									+ properties->GetConductivity()->vec[0] );
	Double		pressure = flame->GetPressure();
	Double		*temp = flame->GetTemperature()->vec;
	Double		tempFirst = currGrid->GetY()->mat[0][fTemperature];
	Double		hL = 360064.8;		// latent heat of vaporization of the fuel at T = 373 K ( boiling point )
	Double		A = 9.01875;		// coefficients of the Antoine equation
	Double		B = 1264.37;
	Double		C = -56.514;
#ifndef CHECKIT
	Double		vRatio;
#endif
	
	yLeft[fTemperature] = bcLeft[fTemperature] = B / ( A - log10( pressure ) ) - C;
	yLeft[fVVelocity] = bcLeft[fVVelocity] = lambdaLeft * ( tempFirst - yLeft[fTemperature] ) / ( hFirst * hL );

#ifndef CHECKIT
	vRatio = yLeft[fVVelocity] / VOld;

	fprintf( stderr, "TempOld = %g\tTempNew = %g\n", TempOld, bcLeft[fTemperature] );
	fprintf( stderr, "VOld = %g\tVNew = %g\n", VOld, bcLeft[fVVelocity] );
	
	yLeft[fVVelocity] = bcLeft[fVVelocity] = ( vRatio < 1.0 ) ? ( MAX( 0.9 , vRatio ) * VOld ) 
																: ( MIN( 1.1 , vRatio ) * VOld );

	yLeft[fVVelocity] = bcLeft[fVVelocity] = 0.5 * ( yLeft[fVVelocity] + VOld );

	fprintf( stderr, "Tmax = %g @ a = %g 1/s\n"
		, temp[LocationOfMax( currGrid->GetNGridPoints()+2, &temp[kPrev] ) - 1]
		, flame->GetLiquidPoolStrainRate() );

	fprintf( stderr, "Set V = %g\n", bcLeft[fVVelocity] );

	if ( fabs( VOld - bcLeft[fVVelocity] ) / VOld > eps ||
		 fabs( TempOld - bcLeft[fTemperature] ) / TempOld > eps ) {
		return FALSE;
	}
	else {
		return TRUE;
	}
#endif
	return TRUE;
}

Double TCountDiffPhysEigen::SecondDerivBinSpecDiff( int nVariable, NodeInfoPtr nodeInfo )
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
		aCurr = rho[kCurr] / W2C * Dij[kCurr][i][j];
		aPrev = rho[kPrev] / W2P * Dij[kPrev][i][j];
		aNext = rho[kNext] / W2N * Dij[kNext][i][j];
		sum += ( ( aCurr + aNext ) * hm * ( yNext[nVarj] * WN - y[nVarj] * WC )
				+ ( aPrev + aCurr ) * h * ( yPrev[nVarj] * WP - y[nVarj] * WC ) );
	}
	
	return sum * Wi[i] / nodeInfo->hnenn;

}

void TCountDiffPhysEigen::CalcYLeft( void )
{
	int			i;
	int			nSpecIn = fSpecies->GetNSpeciesInSystem();
	Double		*YLeftNew;
	Double		*yLeft = fSolver->bt->GetGrid()->GetCurrentGrid()->GetYLeft()->vec;

	YLeftNew = GetYLeft( &yLeft[fFirstSpecies] );

	for ( i = 0; i < nSpecIn; ++i ) {
		yLeft[ fFirstSpecies + i ] = YLeftNew[i];
	}

	SetFlameNode( 0 );
	ComputeProperties( fFlameNode, fFlameNode->temp[kCurr]
							, fFlameNode->Y[kCurr], GetPressure() );
}

Double *TCountDiffPhysEigen::GetYLeft( Double *Yguess )
{
	int i;
	int	speciesIn = fSpecies->GetNSpeciesInSystem();
	
	for ( i = 0; i < speciesIn; ++i ) {
		fYLeftVec->vec[i] = MAX(1e-8, Yguess[i]);
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
					/ ( currGrid->GetYLeft()->vec[fVVelocity] * hFirst );
			fYLeftVec->vec[i] = ( bcLeft[fFirstSpecies+i] + coeff * yFirst[fFirstSpecies+i] ) 
				/ ( 1.0 + coeff );
		}
		NewtonSolve( fNewtonInfoL, FALSE, FALSE, this );
		if ( !fNewtonInfoL->converged ) {
		}
		else {
		}
	}

	SetFlameNode( kPrev );

	ComputeProperties( fFlameNode, fSolver->bt->GetGrid()->GetCurrentGrid()->GetYLeft()->vec[fTemperature]
							, fYLeftVec->vec, GetPressure() );

	return fYLeftVec->vec;
}

/////gb//
int BCLeftNewtonFuncs( const VectorPtr /*x*/, VectorPtr fVec, void *object )
{
// called by NewtonSolve

// F = rho v ( Y_i - Y_i^Liquid ) + rho Y_i V_i 
// or
// F = Y_i - Y_i^Liquid + ( rho Y_i V_i ) / v 

	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	TGridPtr 	currGrid = flame->GetSolver()->bt->GetGrid()->GetCurrentGrid();
	Double		*x = currGrid->GetX()->vec;
	Double		hFirst = x[0] - flame->GetSolver()->bt->GetLeft();
	int			speciesIn = flame->GetSpecies()->GetNSpeciesInSystem();
	int			fVVelocity = flame->GetOffsetVVelocity();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();

	Double		*yLeft = flame->GetYLeftVec()->vec;
	Double		*bcLeft = currGrid->GetBcLeft()->vec;
	Double		*f = fVec->vec;
	
	flame->SetFlameNode( 0 );
	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], flame->GetPressure() );

	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, currGrid->GetYLeft()->vec[flame->GetOffsetTemperature()]
							, yLeft, flame->GetPressure() );

#ifdef CHECKIT
	if ( flame->fLiquidPoolBC ) {
		flame->UpdateBCLiquidPool( flame );
	}
#endif

#ifndef FULLDIFFUSION
	Double	*rhoY_iV_iPlus = flame->fRhoY_iV_iPlus->vec;
	flame->CalcAllDiffVeloRhoYiNext( yLeft, hFirst );
#endif
	for ( int i = 0; i < speciesIn; ++i ) {
		f[i] = ( yLeft[i] - MAX( bcLeft[fFirstSpecies+i], 1.0e-15 ) )
#ifdef FULLDIFFUSION
			+ flame->GetPolyDiffVeloRhoYiNext( yLeft, fFirstSpecies+i, hFirst ) / currGrid->GetYLeft()->vec[fVVelocity];
#else
			+ rhoY_iV_iPlus[i]
				/ ( 0.5 * ( currGrid->GetYLeft()->vec[fVVelocity] 
							+ currGrid->GetYLeft()->vec[fVVelocity] ) );
#endif
	}

	return 0;
}

Double TCountDiffPhysEigen::GetDiffVeloRhoYiNext( Double *YCurr, int nVariable, Double hNext )
{
// returns rho Y_i V_i between boundary and first gridpoint
	int		j, i = nVariable - fFirstSpecies;

	int		nSpecIn = fSpecies->GetNSpeciesInSystem();

	Double	*YNext = fFlameNode->Y[kNext];
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	rhoCurr = fFlameNode->mixDensity[kCurr];
	Double	rhoNext = fFlameNode->mixDensity[kNext];
	Double	WC = fFlameNode->mixMolarMass[kCurr];
	Double	WN = fFlameNode->mixMolarMass[kNext];
	Double	sum = 0.0;
	Double	coeffCurr;
	Double	coeffNext;
	Double	sumYCurr = 0.0;
	Double	sumYNext = 0.0;
	
	sumYCurr = 1.0;
	sumYNext = 1.0;
	
// ordinary diffusion Y	
	coeffCurr = rhoCurr * diffusivity[i];
	coeffNext = rhoNext * diffusivityNext[i];
	sum -= ( coeffCurr + coeffNext ) * ( YNext[i] - YCurr[i] );

#ifdef MOLARDIFFUSION
// ordinary diffusion W	
	coeffCurr *= YCurr[i] / WC;
	coeffNext *= YNext[i] / WN;
	sum -= ( coeffCurr + coeffNext ) * ( WN - WC );
#endif

#ifdef DIFFUSIVITYCORRECTION
// diffusion correction Y	
	for ( j = 0; j < nSpecIn; ++j ) {
//		ordinary diffusion Y	
		coeffCurr = rhoCurr * YCurr[i] * diffusivity[j] / sumYCurr;
		coeffNext = rhoNext * YNext[i] * diffusivityNext[j] / sumYNext;
		sum += ( coeffCurr + coeffNext ) * ( YNext[j] - YCurr[j] );
		
#	ifdef MOLARDIFFUSION
//		ordinary diffusion W	
		coeffCurr *= YCurr[j] / WC;
		coeffNext *= YNext[j] / WN;
		sum += ( coeffCurr + coeffNext ) * ( WN - WC );
#	endif
	}
#endif
	
	return 0.5 * sum / hNext;
}

void TCountDiffPhysEigen::CalcAllDiffVeloRhoYiNext( Double *YCurr, Double hNext )
{
// returns rho Y_i V_i between boundary and first gridpoint
	int		i, nSpecIn = fSpecies->GetNSpeciesInSystem();

	Double	*YNext = fFlameNode->Y[kNext];
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	*rhoY_iV_iPlus = fRhoY_iV_iPlus->vec;
	Double	rhoCurr = fFlameNode->mixDensity[kCurr];
	Double	rhoNext = fFlameNode->mixDensity[kNext];
	Double	WC = fFlameNode->mixMolarMass[kCurr];
	Double	WN = fFlameNode->mixMolarMass[kNext];
	Double	diffPlus, rhoV_cPlus = 0.0;
	
	diffPlus = - 0.5 * ( rhoCurr / WC + rhoNext / WN );
	for ( i = 0; i < nSpecIn; ++i ) {
		rhoY_iV_iPlus[i] = diffPlus * 0.5 * ( diffusivity[i] + diffusivityNext[i] ) * ( YNext[i] * WN - YCurr[i] * WC ) / hNext;
		rhoV_cPlus -= rhoY_iV_iPlus[i];
	}
	for ( i = 0; i < nSpecIn; ++i ) {
		rhoY_iV_iPlus[i] += 0.5 * ( YCurr[i] + YNext[i] ) * rhoV_cPlus;
	}
}

Double TCountDiffPhysEigen::GetPolyDiffVeloRhoYiNext( Double *YCurr, int nVariable, Double hNext )
{
// returns rho Y_i V_i between boundary and first gridpoint
	int		j, i = nVariable - fFirstSpecies;
	int		nSpecIn = fSpecies->GetNSpeciesInSystem();
	Double	*YNext = fFlameNode->Y[kNext];
	Double	***Dij = fFlameNode->binDij;
	Double	*Wi = fSpecies->GetMolarMass()->vec;
	Double	WC = fFlameNode->mixMolarMass[kCurr];
	Double	WN = fFlameNode->mixMolarMass[kNext];
	Double	*rho = fFlameNode->mixDensity;
	Double	sum = 0.0;
	Double	coeffCurr = rho[kCurr] / ( WC * WC ) * Wi[i];
	Double	coeffNext = rho[kNext] / ( WN * WN ) * Wi[i];
	
	for ( j = 0; j < nSpecIn; ++j ) {
		if ( j == i ) continue;
		sum += 0.5 * ( coeffCurr * Dij[kCurr][i][j] + coeffNext * Dij[kNext][i][j] ) 
							* ( YNext[j] * WN - YCurr[j] * WC ) / hNext;
	}
	
	return sum;
}

Double TCountDiffPhysEigen::SecondDerivXDiffusion( int nVariable, NodeInfoPtr nodeInfo )
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

Double TCountDiffPhysEigen::DiffCorrX( int nVariable, NodeInfoPtr nodeInfo )
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
	
	Double	auxPlus = hm * ( MNext - M );
	Double	auxMinus = h * ( MPrev - M );
	Double	coeffCurrNow;
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		coeffCurrNow = coeffCurr * diffusivity[i] * Y[i];
		diffPlus = coeffCurrNow + coeffNext * diffusivityNext[i] * YNext[i];
		diffMinus = coeffCurrNow + coeffPrev * diffusivityPrev[i] * YPrev[i];
		value += ( diffPlus * auxPlus + diffMinus * auxMinus );
	}

	return value / nodeInfo->hnenn;
}

void TCountDiffPhysEigen::SaveSolution( void )
{
	int		k;
	int		len = fSolV->len;
	Double	*v = fSolV->vec;
	Double	*saveV = fSavedV->vec;
	Double	*G = fSolG->vec;
	Double	*saveG = fSavedG->vec;
	Double	*mixFrac = fSolMixFrac->vec;
	Double	*saveMixFrac = fSavedMixFrac->vec;
	Double	*P = fSolP->vec;
	Double	*saveP = fSavedP->vec;

	T1DFlame::SaveSolution();
	fSavedV->len = fSolV->len;
	fSavedG->len = fSolG->len;
	fSavedMixFrac->len = fSolMixFrac->len;
	fSavedP->len = fSolP->len;

	for ( k = -1; k <= len; ++k ) {
		saveV[k] = v[k];
		saveG[k] = G[k];
		saveMixFrac[k] = mixFrac[k];
		saveP[k] = P[k];
	}
}

void TCountDiffPhysEigen::SolutionToSolver( void )
{
	TNewtonPtr	bt = fSolver->bt;
	TGridPtr	grid = bt->GetGrid()->GetFine();
	int		nGridPoints = fSolV->len;
	Double	*V = fSolV->vec;
	Double	*G = fSolG->vec;
	Double	*Z = fSolMixFrac->vec;
	Double	*P = fSolP->vec;
	Double	**y = grid->GetY()->mat;

	T1DFlame::SolutionToSolver();
	
	for ( int k = 0; k < nGridPoints; ++k ) {
		y[k][fVVelocity] = V[k];
		y[k][fUVelocity] = G[k];
		y[k][fPStrain] = P[k];
		y[k][fMixFrac] = Z[k];
	}
	
	DiffPhysEigenPostIter( this );
}

void TCountDiffPhysEigen::RestoreSolution( void )
{
	int		k;
	int		len = fSavedV->len;
	Double	*v = fSolV->vec;
	Double	*saveV = fSavedV->vec;
	Double	*G = fSolG->vec;
	Double	*saveG = fSavedG->vec;
	Double	*mixFrac = fSolMixFrac->vec;
	Double	*saveMixFrac = fSavedMixFrac->vec;
	Double	*P = fSolP->vec;
	Double	*saveP = fSavedP->vec;

	UpdateDimensions( len );

	T1DFlame::RestoreSolution();

	for ( k = -1; k <= len; ++k ) {
		v[k] = saveV[k];
		G[k] = saveG[k];
		mixFrac[k] = saveMixFrac[k];
		P[k] = saveP[k];
	}
	
	SolutionToSolver();
}

void TCountDiffPhysEigen::CalcYRight( void )
{
	int			i;
	int			gridPoints = fSolver->bt->GetGrid()->GetCurrentGrid()->GetNGridPoints();
	int			nSpecIn = fSpecies->GetNSpeciesInSystem();
	Double		*YRightNew;
	Double		*yRight = fSolver->bt->GetGrid()->GetCurrentGrid()->GetYRight()->vec;

	YRightNew = GetYRight( &yRight[fFirstSpecies] );

	for ( i = 0; i < nSpecIn; ++i ) {
		yRight[ fFirstSpecies + i ] = YRightNew[i];
	}

	SetFlameNode( gridPoints );
	ComputeProperties( fFlameNode, fFlameNode->temp[kCurr]
							, fFlameNode->Y[kCurr], GetPressure() );
}

Double *TCountDiffPhysEigen::GetYRight( Double *Yguess )
{
	int i;
	int	speciesIn = fSpecies->GetNSpeciesInSystem();
	
	for ( i = 0; i < speciesIn; ++i ) {
		fYRightVec->vec[i] = MAX(1e-8, Yguess[i]);
	}
	Double 		coeff;
	TGridPtr 	currGrid = fSolver->bt->GetGrid()->GetCurrentGrid();
	Double		*bcRight = currGrid->GetBcRight()->vec;
	Double		*x = currGrid->GetX()->vec;
	int			gridPoints = currGrid->GetNGridPoints();
	Double		hLast = fSolver->bt->GetRight() - x[gridPoints-1];
	Double		*yLast = currGrid->GetY()->mat[gridPoints-1];

	NewtonSolve( fNewtonInfoR, FALSE, FALSE, this );
	if ( !fNewtonInfoR->converged ) {

		for ( i = 0; i < speciesIn; ++i ) {
			coeff = fProperties->GetDensity()->vec[gridPoints] * fSpecies->GetDiffusivity()->mat[gridPoints][i] 
					/ ( currGrid->GetYRight()->vec[fVVelocity] * hLast );
			fYRightVec->vec[i] = ( bcRight[fFirstSpecies+i] + coeff * yLast[fFirstSpecies+i] ) 
				/ ( 1.0 + coeff );
		}
		NewtonSolve( fNewtonInfoR, FALSE, FALSE, this );
	}

	SetFlameNode( gridPoints-1 );

	ComputeProperties( fFlameNode, fSolver->bt->GetGrid()->GetCurrentGrid()->GetYRight()->vec[fTemperature]
							, fYRightVec->vec, GetPressure() );

	return fYRightVec->vec;
}

/////gb//
int BCRightNewtonFuncs( const VectorPtr /*x*/, VectorPtr fVec, void *object )
{
// called by NewtonSolve

// F = rho v ( Y_i - Y_i^Liquid ) + rho Y_i V_i 
// or
// F = Y_i - Y_i^Liquid + ( rho Y_i V_i ) / v 

	TCountDiffPhysEigenPtr	flame = ( TCountDiffPhysEigenPtr )object;
	TGridPtr 	currGrid = flame->GetSolver()->bt->GetGrid()->GetCurrentGrid();
	Double		*x = currGrid->GetX()->vec;
	int			gridPoints = currGrid->GetNGridPoints();
	Double		hLast = flame->GetSolver()->bt->GetRight() - x[gridPoints-1];
	int			speciesIn = flame->GetSpecies()->GetNSpeciesInSystem();
	int			fVVelocity = flame->GetOffsetVVelocity();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();

	Double		*yRight = flame->GetYRightVec()->vec;
	Double		*bcRight = currGrid->GetBcRight()->vec;
	Double		*f = fVec->vec;
	
	flame->SetFlameNode( gridPoints );

	flame->ComputeProperties( flame->fFlameNode, currGrid->GetYRight()->vec[flame->GetOffsetTemperature()]
							, yRight, flame->GetPressure() );

#ifndef FULLDIFFUSION
	Double	*rhoY_iV_iPrev = flame->fRhoY_iV_iPlus->vec;
	flame->CalcAllDiffVeloRhoYiPrev( yRight, hLast );
#endif
	for ( int i = 0; i < speciesIn; ++i ) {
		f[i] = ( MAX( bcRight[fFirstSpecies+i], 1.0e-15 ) - yRight[i] )
#ifdef FULLDIFFUSION
			;fprintf( stderr, "FULLDIFFUSION at right boundary not yet implemented\n" );exit(2);
#else
			- rhoY_iV_iPrev[i]
				/ ( 0.5 * ( currGrid->GetYRight()->vec[fVVelocity] 
							+ currGrid->GetYRight()->vec[fVVelocity] ) );
#endif
	}

	return 0;
}

void TCountDiffPhysEigen::CalcAllDiffVeloRhoYiPrev( Double *YCurr, Double hPrev )
{
// returns rho Y_i V_i between boundary and first gridpoint
	int		i, nSpecIn = fSpecies->GetNSpeciesInSystem();

	Double	*YPrev = fFlameNode->Y[kPrev];
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*rhoY_iV_iMinus = fRhoY_iV_iPlus->vec;
	Double	rhoCurr = fFlameNode->mixDensity[kCurr];
	Double	rhoPrev = fFlameNode->mixDensity[kPrev];
	Double	WC = fFlameNode->mixMolarMass[kCurr];
	Double	WP = fFlameNode->mixMolarMass[kPrev];
	Double	diffMinus, rhoV_cMinus = 0.0;
	
	diffMinus = - 0.5 * ( rhoCurr / WC + rhoPrev / WP );
	for ( i = 0; i < nSpecIn; ++i ) {
		rhoY_iV_iMinus[i] = diffMinus * 0.5 * ( diffusivity[i] + diffusivityPrev[i] ) 
					* ( YCurr[i] * WC - YPrev[i] * WP ) / hPrev;
		rhoV_cMinus -= rhoY_iV_iMinus[i];
	}
	for ( i = 0; i < nSpecIn; ++i ) {
		rhoY_iV_iMinus[i] += 0.5 * ( YCurr[i] + YPrev[i] ) * rhoV_cMinus;
	}
}
