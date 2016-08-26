//#define FOOL_SOFTBENCH( x ) 

#include "FlameMaster.h"
#include "ListTool.h"
#include "Spline.h"
#include "TCountDiffFlamePhys.h"

#define UPWINDCONVECTION

#define DIFFUSIVITYCORRECTION

#define ENTHALPYFLUX

void TCountDiffFlamePhys::InitTCountDiffFlamePhys( void )
{
	TBVPSolverPtr	solver = GetSolver();
	TNewtonPtr		bt = solver->bt;
	TGridPtr		fine = bt->GetGrid()->GetFine();
	TGridPtr		coarse = bt->GetGrid()->GetCoarse();
	int				nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	int				maxGridPoints = bt->GetMaxGridPoints();

//	names of variables
	fVariableNames = new String[fVariablesWithoutSpecies + nSpeciesInSystem];

	fVariableNames[fVVelocity] = new char[2];
	strcpy( fVariableNames[fVVelocity], "V" );
	fVariableNames[fUVelocity] = new char[3];
	strcpy( fVariableNames[fUVelocity], "U" );
	fVariableNames[fMixFrac] = new char[2];
	strcpy( fVariableNames[fMixFrac], "Z" );
	fVariableNames[fTemperature] = new char[2];
	strcpy( fVariableNames[fTemperature], "T" );
	fPrintMolarFractions = fInputData->fPrintMolarFractions;
	
	for ( int i = 0; i < nSpeciesInSystem; ++i ) {
		fVariableNames[fFirstSpecies + i] = new char[strlen( fSpecies->GetNames()[i] ) + 1];
		strcpy( fVariableNames[fFirstSpecies + i], fSpecies->GetNames()[i] );
	}

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

	bt->SetUtFuncs( CountDiffPhysJacFirst, CountDiffPhysJacRest, CountDiffPhysJacLast
				, CountDiffPhysRHSRest, CountDiffPhysRHSRest, CountDiffPhysRHSRest 
				, CountDiffPhysOutput, CountDiffPhysPostIter
				, SetCountDiffPhysNodeInfo, CountDiffPhysPostConv
				, GetCountDiffPhysVarNames );
	SetInitialBC( fine, fInputData );
	SetInitialBC( coarse, fInputData );
	fMassFraction = new TMassFraction( this );
	if ( !fMassFraction ) FatalError( "memory allocation of TMassFraction failed" );
	ReadStartProfiles( fInputData );
	CheckBC();
	CheckInitialGuess();
	UpdateSolution( fine->GetY(), fine->GetYLeft(), fine->GetYRight() );	
}

TCountDiffFlamePhys::~TCountDiffFlamePhys( void )
{
	int	nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();

	delete fMassFraction;

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

void CountDiffPhysJacFirst( void *object, NodeInfoPtr nodeInfo )
{
	TCountDiffFlamePhysPtr	flame = ( TCountDiffFlamePhysPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		M = nodeInfo->nOfEquations;
	int		nOfSpecies = flame->GetSpecies()->GetNOfSpecies();
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
	int		mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
    Double  hnenn = nodeInfo->hnenn;
	Double	**a = nodeInfo->a;
	Double	**massFracs = flameNode->Y;
	Double	*Y = massFracs[kCurr];
	Double	*YPrev = massFracs[kPrev];
	Double	*temp = flameNode->temp;
	Double	*y = nodeInfo->y;
	Double	*yPrev = nodeInfo->yPrev;
	Double	strainRate = flame->GetStrainRate();
	Double	*mixDensity = flameNode->mixDensity;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*mixHeatCapacity = flameNode->mixHeatCapacity;
	Double	oneOverCp = 1.0 / mixHeatCapacity[kCurr];
	Double	*enthalpy = flameNode->enthalpy;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	*mixConductivity = flameNode->mixConductivity;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = 1.0 / ( strainRate * idealGasCoeff );
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];
	Double	geometry = 1.0 + flame->GetGeometry();
	Double	diffMinusH;

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );

#ifdef UPWINDCONVECTION
	Double	hCoeff = nodeInfo->h * ( nodeInfo->h + nodeInfo->hm );

// second equation ( momentum )
	FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0 );
// third equation ( mixturefraction )
	FillJacNonlinearConvectUpwind( fVVelocity, fMixFrac, nodeInfo, 1.0 );
	if ( mixtureSpecificationLeft == kMassFlux && y[fVVelocity] > 0.0 ) {
		Double	zCoeff = mixConductivity[kPrev] / ( mixHeatCapacity[kPrev] * yPrev[fVVelocity] * nodeInfo->hm );
		a[fMixFrac][fMixFrac] -= hCoeff * y[fVVelocity] * zCoeff / ( 1.0 + zCoeff );
	}
// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0 );
		if ( mixtureSpecificationLeft == kMassFlux && y[fVVelocity] > 0.0 ) {
			a[speciesEq][speciesEq] -= hCoeff * y[fVVelocity] 
									* flame->GetdYPrevdY( speciesEq-fFirstSpecies, nodeInfo );
		}
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0 );
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
	a[fUVelocity][fVVelocity] += geometry * strainRate * mixDensity[kCurr] * hnenn;
	a[fTemperature][fVVelocity] -= geometry * strainRate * mixDensity[kCurr] * y[fUVelocity] / temp[kCurr] * hnenn;
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		a[speciesEq][fVVelocity] -= geometry * strainRate * mixDensity[kCurr] * y[fUVelocity] 
									* mixMolarMass / molarMass[speciesEq-fFirstSpecies] * hnenn;
	}
	
// second equation ( momentum )
	FillJacWithDiffusion( fUVelocity, fUVelocity, 1.0, mixViscosity, nodeInfo, kNegative );
	a[fUVelocity][fUVelocity] += 2.0 * strainRate * mixDensity[kCurr] * y[fUVelocity] * hnenn;
	a[fTemperature][fUVelocity] -= y[fUVelocity] * y[fUVelocity] / ( temp[kCurr] * temp[kCurr] * ROverAPM ) * hnenn;
	
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
			
			for ( speciesVar = 0; speciesVar < nOfSpecies; ++speciesVar ) {
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
	}
}

Double TCountDiffFlamePhys::GetdYPrevdY( int speciesIndex, NodeInfoPtr nodeInfo )
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

void CountDiffPhysJacRest( void *object, NodeInfoPtr nodeInfo )
{
	TCountDiffFlamePhysPtr	flame = ( TCountDiffFlamePhysPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		M = nodeInfo->nOfEquations;
	int		nOfSpecies = flame->GetSpecies()->GetNOfSpecies();
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
    Double  hnenn = nodeInfo->hnenn;
	Double	**a = nodeInfo->a;
	Double	*temp = flameNode->temp;
	Double	*y = nodeInfo->y;
	Double	strainRate = flame->GetStrainRate();
	Double	mixDensity = *flameNode->mixDensity;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	oneOverCp = 1.0 / flameNode->mixHeatCapacity[kCurr];
	Double	*enthalpy = flameNode->enthalpy;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = 1.0 / ( strainRate * idealGasCoeff );
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];
	Double	geometry = 1.0 + flame->GetGeometry();

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

/*	FILE *fp = fopen( "dmdy.tout", "w" );
	Print2DArray( dMdY, nSpeciesInSystem, nSpeciesInSystem, kColumnPointers, gPrnt, fp );
	fclose( fp );
	fprintf( stderr, "dMdY dumped\n" );
	fp = fopen( "dmdt.tout", "w" );
	Print1DArray( dMdT, nSpeciesInSystem, gPrnt, fp );
	fclose( fp );
	fprintf( stderr, "dMdT dumped\n" );
	exit(2);*/
// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );

#ifdef UPWINDCONVECTION
	// second equation ( momentum )
	FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0 );
	// third equation ( mixturefraction )
	FillJacNonlinearConvectUpwind( fVVelocity, fMixFrac, nodeInfo, 1.0 );
	// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0 );
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0 );
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
	a[fUVelocity][fVVelocity] += geometry * strainRate * mixDensity * hnenn;
	a[fTemperature][fVVelocity] -= geometry * strainRate * mixDensity * y[fUVelocity] / temp[kCurr] * hnenn;
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		a[speciesEq][fVVelocity] -= geometry * strainRate * mixDensity * y[fUVelocity] 
									* mixMolarMass / molarMass[speciesEq-fFirstSpecies] * hnenn;
	}
	
// second equation ( momentum )
	FillJacWithDiffusion( fUVelocity, fUVelocity, 1.0, mixViscosity, nodeInfo, kNegative );
	a[fUVelocity][fUVelocity] += 2.0 * strainRate * mixDensity * y[fUVelocity] * hnenn;
	a[fTemperature][fUVelocity] -= y[fUVelocity] * y[fUVelocity] / ( temp[kCurr] * temp[kCurr] * ROverAPM ) * hnenn;
	
// third equation ( mixturefraction )
	flame->FillJacMixFracDiffusion( fMixFrac, fMixFrac, nodeInfo, kNegative );

// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, speciesIndexEq, 1.0, nodeInfo, kNegative );
#ifdef DIFFUSIVITYCORRECTION
		flame->FillJacDiffCorr( speciesEq, 1.0, nodeInfo );
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
		}
	}
}

void CountDiffPhysJacLast( void *object, NodeInfoPtr nodeInfo )
{
	TCountDiffFlamePhysPtr	flame = ( TCountDiffFlamePhysPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		M = nodeInfo->nOfEquations;
	int		nOfSpecies = flame->GetSpecies()->GetNOfSpecies();
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		speciesEq, speciesVar, speciesIndexEq, speciesIndexVar;
	int		lastSpeciesEq = nSpeciesInSystem + fFirstSpecies;
	int		mixtureSpecificationRight = flame->GetMixtureSpecificationRight();
    Double  hnenn = nodeInfo->hnenn;
	Double	**a = nodeInfo->a;
	Double	**massFracs = flameNode->Y;
	Double	*Y = massFracs[kCurr];
	Double	*YNext = massFracs[kNext];
	Double	*temp = flameNode->temp;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	strainRate = flame->GetStrainRate();
	Double	*mixDensity = flameNode->mixDensity;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*mixHeatCapacity = flameNode->mixHeatCapacity;
	Double	oneOverCp = 1.0 / mixHeatCapacity[kCurr];
	Double	*enthalpy = flameNode->enthalpy;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	*mixConductivity = flameNode->mixConductivity;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = 1.0 / ( strainRate * idealGasCoeff );
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];
	Double	geometry = 1.0 + flame->GetGeometry();
	Double	diffPlusHm;

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );

#ifdef UPWINDCONVECTION
	Double	hCoeff = nodeInfo->hm * ( nodeInfo->h + nodeInfo->hm );

// second equation ( momentum )
	FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0 );
// third equation ( mixturefraction )
	FillJacNonlinearConvectUpwind( fVVelocity, fMixFrac, nodeInfo, 1.0 );
	if ( mixtureSpecificationRight == kMassFlux && y[fVVelocity] < 0.0 ) {
		Double	zCoeff = mixConductivity[kNext] / ( mixHeatCapacity[kNext] * yNext[fVVelocity] * nodeInfo->h );
		a[fMixFrac][fMixFrac] -= hCoeff * y[fVVelocity] * zCoeff / ( 1.0 - zCoeff );
	}
// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0 );
		if ( mixtureSpecificationRight == kMassFlux && y[fVVelocity] < 0.0 ) {
			a[speciesEq][speciesEq] += hCoeff * y[fVVelocity] 
									* flame->GetdYNextdY( speciesEq-fFirstSpecies, nodeInfo );
		}
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0 );
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
	a[fUVelocity][fVVelocity] += geometry * strainRate * mixDensity[kCurr] * hnenn;
	a[fTemperature][fVVelocity] -= geometry * strainRate * mixDensity[kCurr] * y[fUVelocity] / temp[kCurr] * hnenn;
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		a[speciesEq][fVVelocity] -= geometry * strainRate * mixDensity[kCurr] * y[fUVelocity] 
									* mixMolarMass / molarMass[speciesEq-fFirstSpecies] * hnenn;
	}
	
// second equation ( momentum )
	FillJacWithDiffusion( fUVelocity, fUVelocity, 1.0, mixViscosity, nodeInfo, kNegative );
	a[fUVelocity][fUVelocity] += 2.0 * strainRate * mixDensity[kCurr] * y[fUVelocity] * hnenn;
	a[fTemperature][fUVelocity] -= y[fUVelocity] * y[fUVelocity] / ( temp[kCurr] * temp[kCurr] * ROverAPM ) * hnenn;
	
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
			Double	*diffusivity = flameNode->diffusivity;
			Double	*diffusivityNext = flameNode->diffusivityNext;
			Double	diffPlusHm = ( diffusivityNext[speciesIndexEq] * mixDensity[kNext]
							+ diffusivity[speciesIndexEq] * mixDensity[kCurr] ) * nodeInfo->hm;
			a[speciesEq][speciesEq] -= diffPlusHm * flame->GetdYNextdY( speciesIndexEq, nodeInfo );
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
			
			for ( speciesVar = 0; speciesVar < nOfSpecies; ++speciesVar ) {
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
		FillJacWithDiffusion( fTemperature, fTemperature, oneOverCp, mixConductivity, nodeInfo, kNegative );
		for ( speciesIndexEq = 0; speciesIndexEq < nSpeciesInSystem; ++speciesIndexEq ) {
			a[fTemperature][fTemperature] += oneOverCp * ( dMdT[speciesIndexEq] * enthalpy[speciesIndexEq] ) * hnenn;
			for ( speciesIndexVar = 0; speciesIndexVar < nSpeciesInSystem; ++speciesIndexVar ) {
				a[speciesIndexVar + fFirstSpecies][fTemperature] += oneOverCp * dMdY[speciesIndexVar][speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
			}
		}
	}
}

Double TCountDiffFlamePhys::GetdYNextdY( int speciesIndex, NodeInfoPtr nodeInfo )
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

	return -coeff / ( 1.0 - coeff + coeffCorr );
}

void CountDiffPhysRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
	TCountDiffFlamePhysPtr	flame = ( TCountDiffFlamePhysPtr )object;
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
	Double	strainRate = flame->GetStrainRate();
	Double	*enthalpy = flameNode->enthalpy;
	Double	*mixDensity = flameNode->mixDensity;
	Double	mixDensityInf = flameNode->rhoInf;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*productionRate = flameNode->productionRate;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	oneOverCp = 1.0 / flameNode->mixHeatCapacity[kCurr];
	Double	sumCpDdYdx;
	Double	sumMH;
	
// first fill all convection terms
	// first equation ( mass )
	rhs[fVVelocity] += FirstDerivUpwind( y[fVVelocity], yPrev[fVVelocity], hm );

#ifdef UPWINDCONVECTION
	// second equation ( momentum )
	rhs[fUVelocity] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fUVelocity], y[fUVelocity], yNext[fUVelocity], hm, h );
	// third equation ( mixturefraction )
	rhs[fMixFrac] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fMixFrac], y[fMixFrac], yNext[fMixFrac], hm, h );
	// fourth to four + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		rhs[speciesEq] += NonlinearConvectUpwind( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
	}
	if ( fTemperature < M ) {
		rhs[fTemperature] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
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
	rhs[fVVelocity] += ( 1.0 + flame->GetGeometry() ) * mixDensity[kCurr] * strainRate * y[fUVelocity];

// momentum equation
	rhs[fUVelocity] -= SecondDerivDiffusion( fUVelocity, mixViscosity, nodeInfo );
	rhs[fUVelocity] -= strainRate * mixDensityInf;
	rhs[fUVelocity] += strainRate * mixDensity[kCurr] * y[fUVelocity] * y[fUVelocity];
	
// mixture fraction equation
	rhs[fMixFrac] -= flame->MixFracDiffusion( fMixFrac, nodeInfo );
	
// fFirstSpecies to fFirstSpecies + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		eqLoop = speciesEq - fFirstSpecies;
		rhs[speciesEq] -= flame->SpeciesDiffusion( speciesEq, eqLoop, nodeInfo );
#ifdef DIFFUSIVITYCORRECTION
		rhs[speciesEq] += flame->DiffCorr( speciesEq, nodeInfo );
#endif
		rhs[speciesEq] -= productionRate[eqLoop];
	}
	
// energy equation
	if ( fTemperature < M ) {
		sumCpDdYdx = 0.0;
		sumMH = 0.0;
#ifdef DIFFUSIVITYCORRECTION
		Double sumCpY = 0.0;
#endif

// diffusion
		rhs[fTemperature] -= oneOverCp * SecondDerivDiffusion( fTemperature, flameNode->mixConductivity, nodeInfo );
		
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

#ifdef ENTHALPYFLUX
//	enthalpy flux
		rhs[fTemperature] -= oneOverCp * mixDensity[kCurr] * sumCpDdYdx 
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
#	ifdef DIFFUSIVITYCORRECTION
		rhs[fTemperature] += oneOverCp * mixDensity[kCurr] * sumCpY
						* FirstDeriv( yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h )
						* flameNode->diffCorr[kCurr];  
#	endif
#endif
		
		rhs[fTemperature] += oneOverCp * sumMH;
		
		if ( flame->fProperties->GetRadiation() ) {
			rhs[fTemperature] -= oneOverCp * flameNode->radiation[kCurr];
		}
		
	}

	//	if ( mode == kRHS ) {
		for ( eqLoop = 0; eqLoop < M; ++eqLoop ){
			rhs[eqLoop] *= - hnenn;
		}
	//	}
}

/*void TCountDiffFlamePhys::FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
{
// fills the jacobian with     constCoeff * sum_j ( d/dy(rho Y_k D_j/sum(Y_l) * dY_j/dy) )

	if ( sign == kNegative ) {
		constCoeff *= -1.0;
	}

	int		i, lVar, nSpecies = fSpecies->GetNOfSpecies();
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	**a = nodeInfo->a;
	Double	**b = nodeInfo->b;
	Double	**c = nodeInfo->c;
	Double	*y = nodeInfo->y;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*yNext = nodeInfo->yNext;
	Double	*density = fFlameNode->mixDensity;
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	coeffCurr = constCoeff  * density[kCurr] * y[nVariable];
	Double	coeffPrev = constCoeff  * density[kPrev] * yPrev[nVariable];
	Double	coeffNext = constCoeff  * density[kNext] * yNext[nVariable];
	Double	diffPlusHm, diffMinusH;
	Double	sumYCurr = 0.0;
	Double	sumYNext = 0.0;
	Double	sumYPrev = 0.0;
	
	for ( i = 0; i < nSpecies; ++i ) {
		lVar = i + fFirstSpecies;
		sumYCurr += y[lVar];
		sumYPrev += yPrev[lVar];
		sumYNext += yNext[lVar];
	}
	coeffCurr /= sumYCurr;
	coeffPrev /= sumYPrev;
	coeffNext /= sumYNext;
	
	for ( i = 0; i < nSpecies; ++i ) {
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
*/

/*Double TCountDiffFlamePhys::DiffCorr( int nVariable, NodeInfoPtr nodeInfo )
{
// returns     sum_j ( d/dy(rho Y_k D_j dY_j/dy) )

	int		i, lVar, nSpecies = fSpecies->GetNOfSpecies();
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	*y = nodeInfo->y;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*yNext = nodeInfo->yNext;
	Double	*density = fFlameNode->mixDensity;
	Double	*diffusivity = fFlameNode->diffusivity;
	Double	*diffusivityPrev = fFlameNode->diffusivityPrev;
	Double	*diffusivityNext = fFlameNode->diffusivityNext;
	Double	coeffCurr = density[kCurr] * y[nVariable];
	Double	coeffPrev = density[kPrev] * yPrev[nVariable];
	Double	coeffNext = density[kNext] * yNext[nVariable];
	Double	diffPlus, diffMinus;
	Double	value = 0.0;
	Double	sumYCurr = 0.0;
	Double	sumYNext = 0.0;
	Double	sumYPrev = 0.0;
	
	for ( i = 0; i < nSpecies; ++i ) {
		lVar = i + fFirstSpecies;
		sumYCurr += y[lVar];
		sumYPrev += yPrev[lVar];
		sumYNext += yNext[lVar];
	}
	coeffCurr /= sumYCurr;
	coeffPrev /= sumYPrev;
	coeffNext /= sumYNext;
	
	for ( i = 0; i < nSpecies; ++i ) {
		lVar = i + fFirstSpecies;
		diffPlus = coeffCurr * diffusivity[i] + coeffNext * diffusivityNext[i];
		diffMinus = coeffCurr * diffusivity[i] + coeffPrev * diffusivityPrev[i];
		value += ( diffPlus * hm * ( yNext[lVar] - y[lVar] ) 
					+ diffMinus * h * ( yPrev[lVar] - y[lVar] ) );
	}

	return value / nodeInfo->hnenn;
}
*/

void TCountDiffFlamePhys::FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
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

Double TCountDiffFlamePhys::DiffCorr( int nVariable, NodeInfoPtr nodeInfo )
{
// returns     sum_j ( d/dy(rho Y_k D_j / sum(Y_i) dY_j/dy) )

	int		i, nSpecies = fSpecies->GetNOfSpecies();
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

int CountDiffPhysPostIter( void *object )
{
	TCountDiffFlamePhysPtr	flame = ( TCountDiffFlamePhysPtr )object;
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
	Double		*yLast = y[nGridPoints-1];
	Double		*yFirst = y[0];
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
	Double		pressure = flame->GetPressure();
	Double		strainRate = flame->GetStrainRate();
	Double		*temp = flame->GetTemperature()->vec;
	Double		**Y = flame->GetMassFracs()->mat;
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();
	
// first set temperature and massfractions 
	for ( i = 0; i < nGridPoints; ++i ) {
		bt->SetNodeInfo( flame, i );		
		flame->fMassFraction->UpdateMixFracDependence( flame, nodeInfo, TRUE );
		if ( fTemperature < bt->GetNEquations() ) {
			break;
		}
	}

	for ( i = 0; i < nGridPoints; ++i ) {
		if ( flame->CheckSolution( y[i][fTemperature], &y[i][fFirstSpecies], nSpeciesInSystem ) ) {
			return 1;
		}
	}

//	update properties
	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, temp[kPrev], Y[kPrev], pressure );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, temp[nGridPoints], Y[nGridPoints], pressure );

// left boundary
	if ( mixtureSpecificationLeft == kMassFlux ) {	// means tsuji bc
		Double 		coeff;
		Double		corrCoeff;
		yLeft[fUVelocity] = 0.0;

		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			if ( bcFlagLeft[fFirstSpecies] == kDirichlet ) {
#ifdef DIFFUSIVITYCORRECTION
				corrCoeff = rhoLeft / yLeft[fVVelocity] * flame->GetDiffCorr()->vec[kPrev];
#else
				corrCoeff = 0.0;
#endif
				coeff = rhoLeft * diffusivityLeft[i] / ( yLeft[fVVelocity] * hFirst );
				yLeft[fFirstSpecies+i] = ( bcLeft[fFirstSpecies+i] + coeff * yFirst[fFirstSpecies+i] ) 
						/ ( 1.0 + coeff + corrCoeff );
			}
			else {
				cerr << "error: i can't handle boundary condition of kind " << bcFlagLeft[fFirstSpecies] 
					<< "for species no. " << i << "at the left boundary" << NEWL;
			}
		}

		coeff = lambdaLeft / ( cpLeft * yLeft[fVVelocity] * hFirst );
		yLeft[fMixFrac] = ( 1.0 + coeff * yFirst[fMixFrac] )
						/ ( 1.0 + coeff );
	}
	else {
		yLeft[fUVelocity] = sqrt( rhoRight / rhoLeft );
	}


// right boundary
	yRight[fVVelocity] = yLast[fVVelocity] - hLast * strainRate * rhoRight * yRight[fUVelocity];
	if ( mixtureSpecificationRight == kMassFlux ) {	// means tsuji bc
		Double 		coeff;
		Double		corrCoeff;
		yRight[fUVelocity] = 0.0;

		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			if ( bcFlagRight[fFirstSpecies] == kDirichlet ) {
#ifdef DIFFUSIVITYCORRECTION
				corrCoeff = rhoRight / yRight[fVVelocity] * flame->GetDiffCorr()->vec[nGridPoints];
#else
				corrCoeff = 0.0;
#endif
				coeff = rhoRight * diffusivityRight[i] / ( yRight[fVVelocity] * hLast );
				yRight[fFirstSpecies+i] = ( bcRight[fFirstSpecies+i] 
											- coeff * yLast[fFirstSpecies+i] ) 
												/ ( 1.0 - coeff + corrCoeff );
			}
			else {
				cerr << "error: i can't handle boundary condition of kind " << bcFlagRight[fFirstSpecies] 
					<< "for species no. " << i << "at the right boundary" << NEWL;
			}
		}
		
		coeff = lambdaRight / ( cpRight * yRight[fVVelocity] * hLast );
		yRight[fMixFrac] = - coeff * yRight[fMixFrac]
						/ ( 1.0 - coeff );
	}
	else {
		yRight[fUVelocity] = 1.0;
	}
	
	flame->fMassFraction->SetMassFractionsBC( fMixFrac, fFirstSpecies
				, ( mixtureSpecificationLeft == kMassFlux ) ? NULL : yLeft
				, ( mixtureSpecificationRight == kMassFlux ) ? NULL : yRight );

	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	
	flame->UpdateThermoProps();

	return 0;
}

#include "TofZ.h"

void TCountDiffFlamePhys::UpdateDimensions( int len )
{
	T1DFlame::UpdateDimensions( len );
	fSolV->len = len;
	fSolU->len = len;
	fSolMixFrac->len = len;
}

void TCountDiffFlamePhys::UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
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

int	TCountDiffFlamePhys::GetOffsetVVelocity( void )
{
	return fVVelocity; 
}

int	TCountDiffFlamePhys::GetOffsetUVelocity( void )
{
	return fUVelocity; 
}

int	TCountDiffFlamePhys::GetOffsetTemperature( void )
{
	return fTemperature; 
}

int TCountDiffFlamePhys::GetOffsetMixFrac( void )
{
	return fMixFrac; 
}

int	TCountDiffFlamePhys::GetOffsetFirstSpecies( void ) 
{
	return fFirstSpecies;
}

ConstStringArray TCountDiffFlamePhys::GetVariableNames( void )
{
	return ( ConstStringArray )fVariableNames;
}

int TCountDiffFlamePhys::GetVariablesWithoutSpecies( void )
{
	return fVariablesWithoutSpecies;
}

void TCountDiffFlamePhys::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp )
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
	Flag				oxidizerFound = FALSE;
	Flag				chooseInputGrid = FALSE;
	Double				*xIn = new Double[gridPointsIn];
	if ( !xIn ) FatalError( "memory allocation of TCountDiffFlamePhys failed" );
	Double				*yIn =  new Double[gridPointsIn];
	if ( !yIn ) FatalError( "memory allocation of TCountDiffFlamePhys failed" );
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
	Double				strainRateIn;
	struct _parameter	*param = GetParameter( "strainrate" );

	if ( GetStrainRate() < 0.0 ) {
	// get scalar strain rate from input file
		if ( param ) {
			strainRateIn = (Double)param->what.quantity.value;
		}
		else { // choose default
			cerr << "#warning: no value for 'strain dissipation rate' in inputfile" << NEWL;
			strainRateIn = 100.0;
		}
		fStrainRate->vec[0] = strainRateIn;
		fStrainRate->len = 0;
		cerr << "initial strain dissipation rate is " << GetStrainRate() << NEWL;
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
		cerr << "error: can't find massfraction of oxidizer 'massfraction-o2'" << NEWL;
		exit(2);
	}
	
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
		else if ( strncmp( string, "v", 1 ) == 0 ) {
			variable = fVVelocity;
		}
		else if ( strcmp( string, "z" ) == 0 ) {
			variable = fMixFrac;
		}
		else if ( strncmp( string, "temperature", 11 ) == 0 ) {
			variable = fTemperature;
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
	
// set initial Boundary values
//	UpdateSolution( yMat, yLeftVec, yRightVec );

//	update properties
	UpdateSolution( yMat, yLeftVec, yRightVec );
	SetFlameNode( kPrev );
	ComputeProperties( fFlameNode, temp[kPrev], Y[kPrev], GetPressure() );
	SetFlameNode( nGridPoints );
	ComputeProperties( fFlameNode, temp[nGridPoints], Y[nGridPoints], GetPressure() );

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
	if ( GetSpecies()->IsConstantLewisNumber() ) {
		CompLewisNumbers( GetSpecies()->GetLewisNumberFile() );
	}
	CountDiffPhysPostIter( this );

	FreeSpline( theSpline );
	delete yIn;
	delete xIn;

	adapGrid->SetSolutionScaler();
	
	fp = GetOutfile( "initialguess", TFlame::kData );
	bt->PrintSolution( locX, y, GetVariableNames(), fp );
	fclose(fp);
}

Double TCountDiffFlamePhys::MixFracDiffusion( int nVariable, NodeInfoPtr nodeInfo )
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

Double TCountDiffFlamePhys::SpeciesDiffusion( int nVariable, int speciesIndex, NodeInfoPtr nodeInfo )
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

void TCountDiffFlamePhys::FillJacMixFracDiffusion( int nVariable, int nEquation, NodeInfoPtr nodeInfo, Flag sign )
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

void TCountDiffFlamePhys::FillJacSpeciesDiffusion( int nVariable, int speciesIndex, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
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

void CountDiffPhysOutput( void *object, FILE *fp, char* tail )
{
	TCountDiffFlamePhysPtr	flame = ( TCountDiffFlamePhysPtr )object;
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
	Double			*U = flame->GetU()->vec;
	Double			*Z = flame->GetZ()->vec;
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
	int				fUVelocity = flame->GetOffsetUVelocity();
	int				fVVelocity = flame->GetOffsetVVelocity();
	int				fMixFrac = flame->GetOffsetMixFrac();
	time_t			theDate;
	char			buffer[80];
	ConstStringArray varNames = flame->GetVariableNames();
	char			**names = species->GetNames();
	VectorPtr 		etaVec = NewVector( gridPoints + 2 );
	Double			*eta = etaVec->vec;
	VectorPtr 		scalarDissVec = NewVector( gridPoints + 2 );
	Double			*scalarDiss = scalarDissVec->vec;
	Double			stoechScalarDiss;
	Flag			fpOpen = FALSE;
	
	if ( !fp ) {
		fpOpen = TRUE;
		fp = flame->GetOutputFile( NULL, tail, TFlame::kNone );
	}

	flame->ScalarDissipation( bt, scalarDissVec, &stoechScalarDiss );
	
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
	fprintf( fp, "strainratev1v2 = %g [1/s]\n", ( yLeft[fVVelocity]/rho[-1] - yRight[fVVelocity]/rho[gridPoints] ) 
													/ bt->GetRight() - bt->GetLeft() );
	fprintf( fp, "vLeft = %g [m/s]\n", yLeft[fVVelocity]/rho[-1] );
	fprintf( fp, "vRight = %g [m/s]\n", yRight[fVVelocity]/rho[gridPoints] ); 

	if ( species->IsConstantLewisNumber() ) {
		fprintf( fp, "ConstantLewisNumbers = \"True\"\n" );
	}
	
	fprintf( fp, "Tmax = %g [K]\n", temp[LocationOfMax( gridPoints+2, &temp[kPrev] ) - 1] );
	fprintf( fp, "StScalarDissRate = %g [1/s]\n", stoechScalarDiss );
	
	Double	EIFuel = 0.0;
	for ( i = 0; i < flame->GetNFuels(); ++i ) {
		EIFuel += flame->ComputeEmissionIndex( flame->GetFuelIndex( i ), x );
	}

	int	indNO = species->FindSpecies( "NO" );
	if ( indNO > -1 ) {
		fprintf( fp, "EmissionIndexNO = %g [g/kg]\n"
				, -1000.0 * flame->ComputeEmissionIndex( indNO, x ) 
				/ EIFuel );
	}
	
	int	indNO2 = species->FindSpecies( "NO2" );
	if ( indNO2 > -1 ) {
		fprintf( fp, "EmissionIndexNO2 = %g [g/kg]\n"
				, -1000.0 * flame->ComputeEmissionIndex( indNO2, x ) 
				/ EIFuel );
	}
	
#ifdef NOXPRODRATES
	if ( strstr( MECHANISM, "NO.red4" ) ) {
		fprintf( fp, "EmissionIndexThermalNO = %g [g/kg]\n"
				, -1000.0 * flame->fReaction->ComputeEIThermalNO( flame ) / EIFuel );
		fprintf( fp, "EmissionIndexPromptNO = %g [g/kg]\n"
				, -1000.0 * flame->fReaction->ComputeEIPromptNO( flame ) / EIFuel );
		fprintf( fp, "EmissionIndexNitrousNO = %g [g/kg]\n"
				, -1000.0 * flame->fReaction->ComputeEINitrousNO( flame ) / EIFuel );
		fprintf( fp, "EmissionIndexReburnNO = %g [g/kg]\n"
				, -1000.0 * flame->fReaction->ComputeEIReburnNO( flame ) / EIFuel );
	}
#endif
	
	fprintf( fp, "\nOxidizerSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", yRight[tempOffset] );
	for ( i = 0; i < nOfSpecies; ++i ) {
		if ( fabs( massFracs[gridPoints][i] ) > 1.0e-20 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], massFracs[gridPoints][i] );
		}
	}
	fprintf( fp, "end\n\n" );

	fprintf( fp, "FuelSide\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", yLeft[tempOffset] );
	for ( i = 0; i < nOfSpecies; ++i ) {
		if ( fabs( massFracs[kPrev][i] ) > 1.0e-20 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], massFracs[kPrev][i] );
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
	fprintf( fp, "eta\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", eta[k] );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}
		
			
// write solution
	// write V-Velocity, U-Velocity, mixture fraction and temperature
	flame->PrintFlameletVector( gridPoints+2, &V[kPrev], "V [kg/(s*m*m)]", fp );
	flame->PrintFlameletVector( gridPoints+2, &U[kPrev], "U", fp );
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
	
//	write f
	Double	coeff = -sqrt( flame->GetStrainRate() * flame->fFlameNode->rhoInf * flame->fFlameNode->viscosityInf );
	fprintf( fp, "f\n" );
		fprintf( fp, "\t%-.6e", V[kPrev] / coeff );
		for ( k = 0; k < gridPoints; ++k ) {
			fprintf( fp, "\t%-.6e", V[k] / coeff );
			if ( (k+2) % 5 == 0 ) {
				fprintf( fp, "\n" );
			}
		}
		fprintf( fp, "\t%-.6e\n", V[gridPoints] / coeff );
	
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
		
	
	fprintf( fp, "trailer\n" );
	if ( nOfEquations < nOfVariables) {
		fprintf( fp, "number of converged equations is %d\n", nOfEquations );
	}
	
	DisposeVector( scalarDissVec );
	DisposeVector( etaVec );
	if ( fpOpen ) {
		fclose( fp );
	}
}

void TCountDiffFlamePhys::ScalarDissipation( TNewtonPtr bt, VectorPtr scalarDissVec, Double *stoechScalarDiss )
{
	int				k;
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
	
	dZdx = ( Z[0] - Z[kPrev] ) / ( x[0] - left );
	scalarDiss[0] = density[0]  * lambda[0] / cp[0]  * dZdx * dZdx;
	for ( k = 0; k < nOfGridPoints; ++k ) {
		bt->SetNodeInfo( this, k );
		dZdx = FirstDeriv( Z[k-1], Z[k], Z[k+1], nodeInfo->hm, nodeInfo->h );
		scalarDiss[k+1] = density[k]  * lambda[k] / cp[k]  * dZdx * dZdx;
		if ( ( Z[k] - zStoe ) * ( Z[k-1] - zStoe ) <= 0.0 ) {
			*stoechScalarDiss = scalarDiss[k] + ( scalarDiss[k+1] - scalarDiss[k] ) 
							/ ( Z[k] - Z[k-1] ) * ( zStoe - Z[k-1] );
		}
	}
	dZdx = ( Z[nOfGridPoints] - Z[nOfGridPoints-1] ) / ( right - x[nOfGridPoints-1] );
	scalarDiss[nOfGridPoints+1] = density[nOfGridPoints-1]  * lambda[nOfGridPoints-1] / cp[nOfGridPoints-1]  * dZdx * dZdx;
}

/*void TCountDiffFlamePhys::XToEta( TNewtonPtr bt, VectorPtr etaVec )
{
	int			k;
	int			gridPoints = bt->GetCurrentGridPoints();
	int			kBefStoech = -1;
	NodeInfoPtr nodeInfo = bt->GetNodeInfo();
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		factor = 0.5 * sqrt( GetStrainRate() / ( fFlameNode->rhoInf * fFlameNode->viscosityInf ) );
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		*Z = fSolMixFrac->vec;
	Double		*rho = GetProperties()->GetDensity()->vec;
	Double		*x = grid->GetX()->vec;
	Double		*eta = etaVec->vec;
	Double		zStoech = fMassFraction->GetZStoe();
	Double		xStoech;
	
	if ( etaVec->len < gridPoints + 2 ) {
		cerr << "#warning: Vector etaVec too short, values for physical grid are not computed" << NEWL;
		return;
	}
	
	eta[0] = 0.0;
	eta[1] = eta[0] + factor * ( x[0] - left ) * ( rho[-1] + rho[0] );
	for ( k = 1; k < gridPoints; ++k ) {
		eta[k+1] = eta[k] + factor * ( x[k] - x[k-1] ) * ( rho[k-1] + rho[k] );
	}
	eta[gridPoints+1] = eta[gridPoints] + factor * ( right - x[gridPoints-1] ) * ( rho[gridPoints-1] + rho[gridPoints] );
	
	// linear Interpolation
	
	for ( k = 0; k < gridPoints-1; ++k ) {
		if ( ( Z[k] - zStoech ) * ( Z[k+1] - zStoech ) <= 0.0 ) {
			kBefStoech = k;
			break;
		}
	}
	
	if ( kBefStoech < 0 ) {
		cerr << "##warning: can't find the point of stoichiometric mixture fraction" << NEWL;
		return;
	}
	
	xStoech = eta[kBefStoech+1] + ( eta[kBefStoech+2] - eta[kBefStoech+1] ) 
						/ ( Z[kBefStoech+1] - Z[kBefStoech] )
						* ( zStoech - Z[kBefStoech] );
	
	for ( k = 0; k < gridPoints+2; ++k ) {
		eta[k] -= xStoech;
	}
}
*/
void SetCountDiffPhysNodeInfo( int k, void *object )
{
	TCountDiffFlamePhysPtr	flame = ( TCountDiffFlamePhysPtr )object;
	
	flame->SetFlameNode( k );
}

void CountDiffPhysPostConv( void *object )
{
	TCountDiffFlamePhysPtr	flame = ( TCountDiffFlamePhysPtr )object;
	TNewtonPtr 				bt = flame->GetSolver()->bt;
	int						isConverged = bt->GetConvergeNewton();
	Double					coeffSpecies = 1.0 / flame->GetStrainRate();
	Double					coeffTemp = 1.0 / flame->GetStrainRate();
	
	if ( isConverged ) {
		flame->SaveSolution();
		if( flame->AdjustGrid( CountDiffPhysPostIter ) ) {
//			bt->WriteOutput( object, NULL, "Map" );
			return;
		}

		if ( flame->fSensAnal ) {
			flame->SensitivityAnalysis( coeffSpecies, coeffTemp, kPhysical );
		}
//		flame->GetSpecies()->PrintProdRateTerms( "O2", flame );
//		flame->fSpecies->PrintDetailedProdRate( bt, flame->fReaction );
		if ( flame->fReactionFluxes ) {
			flame->ReactionFluxes( kPhysical );
			flame->GetReaction()->PrintReactionRates( flame );
			flame->fReaction->PrintRateCoeffs( flame );
			flame->fReaction->PrintDetailedHeatRelease( flame );
		}
	}
	else {
		flame->RestoreSolution();
		CountDiffPhysPostIter( flame );
	}

	flame->PostConvergence( object );
	CountDiffPhysPostIter( flame );
//	flame->fSpecies->PrintDetailedProdRate( bt, flame->fReaction );
}

ConstStringArray GetCountDiffPhysVarNames( void *object )
{
	TCountDiffFlamePhysPtr	flame = ( TCountDiffFlamePhysPtr )object;
	
	return flame->GetVariableNames();
}

FILE *TCountDiffFlamePhys::GetOutputFile( char *head, char *tail, FileType type )
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
					, ( int ) floor( GetStrainRate() + 0.5 )			// in [1/s]
					, ( int )( tFuel )							// in [K]
					, ( int )( tOxidizer ) 						// in [K]
					, ( tail ) ? tail : "" );

	fp = GetOutfile( name, type );
	delete name;

	return fp;
}

/*void CountDiffPhysUpdateLeftBoundary( void  *object )
{
	TCountDiffFlamePhysPtr	flame = ( TCountDiffFlamePhysPtr )object;
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
	Double		*yLast = y[nGridPoints-1];
	Double		*yFirst = y[0];
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
	Double		pressure = flame->GetPressure();
	Double		strainRate = flame->GetStrainRate();
	Double		*temp = flame->GetTemperature()->vec;
	Double		**Y = flame->GetMassFracs()->mat;
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();
	
//	update properties
	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, temp[kPrev], Y[kPrev], pressure );

// left boundary
	if ( mixtureSpecificationLeft == kMassFlux ) {	// means tsuji bc
		Double 		coeff;
		Double		corrCoeff;
		yLeft[fUVelocity] = 0.0;

		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			if ( bcFlagLeft[fFirstSpecies] == kDirichlet ) {
#ifdef DIFFUSIVITYCORRECTION
				corrCoeff = rhoLeft / yLeft[fVVelocity] * flame->GetDiffCorr()->vec[kPrev];
#else
				corrCoeff = 0.0;
#endif
				coeff = rhoLeft * diffusivityLeft[i] / ( yLeft[fVVelocity] * hFirst );
				yLeft[fFirstSpecies+i] = ( bcLeft[fFirstSpecies+i] + coeff * yFirst[fFirstSpecies+i] ) 
						/ ( 1.0 + coeff + corrCoeff );
			}
			else {
				cerr << "error: i can't handle boundary condition of kind " << bcFlagLeft[fFirstSpecies] 
					<< "for species no. " << i << "at the left boundary" << NEWL;
			}
		}

		coeff = lambdaLeft / ( cpLeft * yLeft[fVVelocity] * hFirst );
		yLeft[fMixFrac] = ( 1.0 + coeff * yFirst[fMixFrac] )
						/ ( 1.0 + coeff );
	}
	else {
		yLeft[fUVelocity] = sqrt( rhoRight / rhoLeft );
	}
	
	flame->fMassFraction->SetMassFractionsBC( fMixFrac, fFirstSpecies
				, ( mixtureSpecificationLeft == kMassFlux ) ? NULL : yLeft
				, NULL );

	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	
	flame->UpdateThermoProps();
}

void CountDiffPhysUpdateRightBoundary( void *object )
{
	TCountDiffFlamePhysPtr	flame = ( TCountDiffFlamePhysPtr )object;
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
	Double		*yLast = y[nGridPoints-1];
	Double		*yFirst = y[0];
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
	Double		pressure = flame->GetPressure();
	Double		strainRate = flame->GetStrainRate();
	Double		*temp = flame->GetTemperature()->vec;
	Double		**Y = flame->GetMassFracs()->mat;
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();
	
//	update properties
	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, temp[nGridPoints], Y[nGridPoints], pressure );


// right boundary
	yRight[fVVelocity] = yLast[fVVelocity] - hLast * strainRate * rhoRight * yRight[fUVelocity];
	if ( mixtureSpecificationRight == kMassFlux ) {	// means tsuji bc
		Double 		coeff;
		Double		corrCoeff;
		yRight[fUVelocity] = 0.0;

		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			if ( bcFlagRight[fFirstSpecies] == kDirichlet ) {
#ifdef DIFFUSIVITYCORRECTION
				corrCoeff = rhoRight / yRight[fVVelocity] * flame->GetDiffCorr()->vec[nGridPoints];
#else
				corrCoeff = 0.0;
#endif
				coeff = rhoRight * diffusivityRight[i] / ( yRight[fVVelocity] * hLast );
				yRight[fFirstSpecies+i] = ( bcRight[fFirstSpecies+i] 
											- coeff * yLast[fFirstSpecies+i] ) 
												/ ( 1.0 - coeff + corrCoeff );
			}
			else {
				cerr << "error: i can't handle boundary condition of kind " << bcFlagRight[fFirstSpecies] 
					<< "for species no. " << i << "at the right boundary" << NEWL;
			}
		}
		
		coeff = lambdaRight / ( cpRight * yRight[fVVelocity] * hLast );
		yRight[fMixFrac] = - coeff * yRight[fMixFrac]
						/ ( 1.0 - coeff );
	}
	else {
		yRight[fUVelocity] = 1.0;
	}
	
	flame->fMassFraction->SetMassFractionsBC( fMixFrac, fFirstSpecies
				, NULL
				, ( mixtureSpecificationRight == kMassFlux ) ? NULL : yRight );

	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	
	flame->UpdateThermoProps();
}
*/
