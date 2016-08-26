//#define FOOL_SOFTBENCH( x ) 

#include "FlameMaster.h"
#include "ListTool.h"
#include "Spline.h"
#include "TCountPremFlamePhys.h"

#define UPWINDCONVECTION

#define DIFFUSIVITYCORRECTION

#define ENTHALPYFLUX

void TCountPremFlamePhys::InitTCountPremFlamePhys( void )
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

	fSolV->vec = &fSolV->vec[kNext];
	fSolU->vec = &fSolU->vec[kNext];

	fSolV->len -= 2;
	fSolU->len -= 2;

	if ( fUseNumericalJac ) {
		bt->SetUtFuncs( NULL, NULL, NULL
					, CountPremPhysRHSRest, CountPremPhysRHSRest, CountPremPhysRHSRest 
					, CountPremPhysOutput, CountPremPhysPostIter
					, SetCountPremPhysNodeInfo, CountPremPhysPostConv
					, GetCountPremPhysVarNames
					, CountPremPhysUpdateLeftBoundary, CountPremPhysUpdateRightBoundary );
	}
	else {
		bt->SetUtFuncs( CountPremPhysJacFirst, CountPremPhysJacRest, CountPremPhysJacLast
					, CountPremPhysRHSRest, CountPremPhysRHSRest, CountPremPhysRHSRest 
					, CountPremPhysOutput, CountPremPhysPostIter
					, SetCountPremPhysNodeInfo, CountPremPhysPostConv
					, GetCountPremPhysVarNames
					, CountPremPhysUpdateLeftBoundary, CountPremPhysUpdateRightBoundary );
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

TCountPremFlamePhys::~TCountPremFlamePhys( void )
{
	int	nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();

	delete fMassFraction;

	fSolV->vec = &fSolV->vec[kPrev];
	fSolU->vec = &fSolU->vec[kPrev];

	DisposeVector( fSolU );
	DisposeVector( fSolV );

	for ( int i = 0; i < nSpeciesInSystem+fVariablesWithoutSpecies; ++i ) {
		delete fVariableNames[i];
	}
	delete fVariableNames;
}

void CountPremPhysJacFirst( void *object, NodeInfoPtr nodeInfo )
{
	TCountPremFlamePhysPtr	flame = ( TCountPremFlamePhysPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
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
	Double	*temp = flameNode->temp;
	Double	*y = nodeInfo->y;
	Double	strainRate = flame->GetStrainRate();
	Double	*mixDensity = flameNode->mixDensity;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*mixHeatCapacity = flameNode->mixHeatCapacity;
	Double	oneOverCp = 1.0 / mixHeatCapacity[kCurr];
	Double	*enthalpy = flameNode->enthalpy;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	*productionRate = flameNode->productionRate;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	*mixConductivity = flameNode->mixConductivity;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = 1.0 / ( strainRate * idealGasCoeff );
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];
	Double	geometry = 1.0 + flame->GetGeometry();
	Double	h = nodeInfo->h;
	Double	*diffusivity = flameNode->diffusivity;
	Double	*diffusivityPrev = flameNode->diffusivityPrev;

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );

#ifdef UPWINDCONVECTION

// second equation ( momentum )
	FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0 );
// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0 );
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0 );
	}
#else
	// second equation ( momentum )
	FillJacNonlinearConvectCentral( fVVelocity, fUVelocity, nodeInfo );
	a[fUVelocity][fUVelocity] -= h * h * y[fVVelocity];
	// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectCentral( fVVelocity, speciesEq, nodeInfo );
		a[speciesEq][speciesEq] -= h * h * y[fVVelocity];
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectCentral( fVVelocity, fTemperature, nodeInfo );
		a[fTemperature][fTemperature] -= h * h * y[fVVelocity];
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
	a[fUVelocity][fUVelocity] -= ( mixViscosity[kPrev] + mixViscosity[kCurr] ) * h;
	a[fUVelocity][fUVelocity] += 2.0 * strainRate * mixDensity[kCurr] * y[fUVelocity] * hnenn;
	a[fTemperature][fUVelocity] -= y[fUVelocity] * y[fUVelocity] / ( temp[kCurr] * temp[kCurr] * ROverAPM ) * hnenn;
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		a[speciesEq][fUVelocity] -= strainRate * mixDensity[kCurr] * y[fUVelocity] * y[fUVelocity] 
									* mixMolarMass / molarMass[speciesEq-fFirstSpecies] * hnenn;
	}
	
// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		speciesIndexEq = speciesEq - fFirstSpecies;
		flame->FillJacSpeciesDiffusion( speciesEq, speciesIndexEq, 1.0, nodeInfo, kNegative );
		a[speciesEq][speciesEq] -= ( diffusivityPrev[speciesIndexEq] * mixDensity[kPrev]
					+ diffusivity[speciesIndexEq] * mixDensity[kCurr] ) * h;

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
		FillJacWithDiffusion( fTemperature, fTemperature, oneOverCp, mixConductivity, nodeInfo, kNegative );
		a[fTemperature][fTemperature] -= oneOverCp * ( mixConductivity[kPrev] + mixConductivity[kCurr] ) * h;
		for ( speciesIndexEq = 0; speciesIndexEq < nSpeciesInSystem; ++speciesIndexEq ) {
			a[fTemperature][fTemperature] += oneOverCp * ( dMdT[speciesIndexEq] * enthalpy[speciesIndexEq] 
														 + productionRate[speciesIndexEq] * heatCapacity[speciesIndexEq] ) * hnenn;
			for ( speciesIndexVar = 0; speciesIndexVar < nSpeciesInSystem; ++speciesIndexVar ) {
				a[speciesIndexVar + fFirstSpecies][fTemperature] += oneOverCp * dMdY[speciesIndexVar][speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
			}
		}
	}
}

Double TCountPremFlamePhys::GetdYPrevdY( int speciesIndex, NodeInfoPtr nodeInfo )
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

void CountPremPhysJacRest( void *object, NodeInfoPtr nodeInfo )
{
	TCountPremFlamePhysPtr	flame = ( TCountPremFlamePhysPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
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
	Double	*productionRate = flameNode->productionRate;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = 1.0 / ( strainRate * idealGasCoeff );
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];
	Double	geometry = 1.0 + flame->GetGeometry();

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );

#ifdef UPWINDCONVECTION
	// second equation ( momentum )
	FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0 );
	// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0 );
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0 );
	}
#else
	// second equation ( momentum )
	FillJacNonlinearConvectCentral( fVVelocity, fUVelocity, nodeInfo );
	// third to three + nOfSpecies equation ( species )
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
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		a[speciesEq][fUVelocity] -= strainRate * mixDensity * y[fUVelocity] * y[fUVelocity] 
									* mixMolarMass / molarMass[speciesEq-fFirstSpecies] * hnenn;
	}

// third to three + nOfSpecies equation ( species )
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
			a[fTemperature][fTemperature] += oneOverCp * ( dMdT[speciesIndexEq] * enthalpy[speciesIndexEq] 
														 + productionRate[speciesIndexEq] * heatCapacity[speciesIndexEq] ) * hnenn;
			for ( speciesIndexVar = 0; speciesIndexVar < nSpeciesInSystem; ++speciesIndexVar ) {
				a[speciesIndexVar + fFirstSpecies][fTemperature] += oneOverCp * dMdY[speciesIndexVar][speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
			}
		}
	}
}

void CountPremPhysJacLast( void *object, NodeInfoPtr nodeInfo )
{
	TCountPremFlamePhysPtr	flame = ( TCountPremFlamePhysPtr )object;
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
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
	Double	*temp = flameNode->temp;
	Double	*y = nodeInfo->y;
	Double	strainRate = flame->GetStrainRate();
	Double	*mixDensity = flameNode->mixDensity;
	Double	*mixViscosity = flameNode->mixViscosity;
	Double	*mixHeatCapacity = flameNode->mixHeatCapacity;
	Double	oneOverCp = 1.0 / mixHeatCapacity[kCurr];
	Double	*enthalpy = flameNode->enthalpy;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	*productionRate = flameNode->productionRate;
	Double	*heatCapacity = flameNode->heatCapacity;
	Double	mixMolarMass = *flameNode->mixMolarMass;
	Double	*mixConductivity = flameNode->mixConductivity;
	Double	idealGasCoeff = flame->GetPressure() * mixMolarMass / RGAS; // rho = idealGasCoeff / T
	Double	ROverAPM = 1.0 / ( strainRate * idealGasCoeff );
	Double	**dMdY = flameNode->dMdY;
	Double	*dMdT = flameNode->dMdY[nSpeciesInSystem];
	Double	geometry = 1.0 + flame->GetGeometry();

	flame->FilldMdYOnePoint( flameNode );
	flame->FilldMdTOnePoint( flameNode );

// first fill all convection terms
	// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );

#ifdef UPWINDCONVECTION
// second equation ( momentum )
	FillJacNonlinearConvectUpwind( fVVelocity, fUVelocity, nodeInfo, 1.0 );
// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		FillJacNonlinearConvectUpwind( fVVelocity, speciesEq, nodeInfo, 1.0 );
	}
	if ( fTemperature < M ) {
		FillJacNonlinearConvectUpwind( fVVelocity, fTemperature, nodeInfo, 1.0 );
	}
#else
	// second equation ( momentum )
	FillJacNonlinearConvectCentral( fVVelocity, fUVelocity, nodeInfo );
	// third to three + nOfSpecies equation ( species )
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
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		a[speciesEq][fUVelocity] -= strainRate * mixDensity[kCurr] * y[fUVelocity] * y[fUVelocity] 
									* mixMolarMass / molarMass[speciesEq-fFirstSpecies] * hnenn;
	}
	
// third to three + nOfSpecies equation ( species )
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
		FillJacWithDiffusion( fTemperature, fTemperature, oneOverCp, mixConductivity, nodeInfo, kNegative );
		for ( speciesIndexEq = 0; speciesIndexEq < nSpeciesInSystem; ++speciesIndexEq ) {
			a[fTemperature][fTemperature] += oneOverCp * ( dMdT[speciesIndexEq] * enthalpy[speciesIndexEq] 
														 + productionRate[speciesIndexEq] * heatCapacity[speciesIndexEq] ) * hnenn;
			for ( speciesIndexVar = 0; speciesIndexVar < nSpeciesInSystem; ++speciesIndexVar ) {
				a[speciesIndexVar + fFirstSpecies][fTemperature] += oneOverCp * dMdY[speciesIndexVar][speciesIndexEq] * enthalpy[speciesIndexEq] * hnenn;
			}
		}
	}
}

Double TCountPremFlamePhys::GetdYNextdY( int speciesIndex, NodeInfoPtr nodeInfo )
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

void CountPremPhysRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode )
{
	TCountPremFlamePhysPtr	flame = ( TCountPremFlamePhysPtr )object;
	if ( !flame->RHSAction( nodeInfo, rhsMode ) ) {
		return;
	}
	
	TFlameNodePtr	flameNode = flame->fFlameNode;
	int 	fFirstSpecies = flame->GetOffsetFirstSpecies();
	int 	fTemperature = flame->GetOffsetTemperature();
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
	// third to three + nOfSpecies equation ( species )
	for ( speciesEq = fFirstSpecies; speciesEq < lastSpeciesEq && speciesEq < M; ++speciesEq ) {
		rhs[speciesEq] += NonlinearConvectUpwind( y[fVVelocity], yPrev[speciesEq], y[speciesEq], yNext[speciesEq], hm, h );
	}
	if ( fTemperature < M ) {
		rhs[fTemperature] += NonlinearConvectUpwind( y[fVVelocity], yPrev[fTemperature], y[fTemperature], yNext[fTemperature], hm, h );
	}
#else
	// second equation ( momentum )
	rhs[fUVelocity] += NonlinearConvectCentral( y[fVVelocity], yPrev[fUVelocity], y[fUVelocity], yNext[fUVelocity], hm, h );

	// third to three + nOfSpecies equation ( species )
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

		for ( eqLoop = 0; eqLoop < M; ++eqLoop ){
			rhs[eqLoop] *= - hnenn;
		}
}

void TCountPremFlamePhys::FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
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
		if ( !nodeInfo->lastPoint ) {
			a[lVar][nVariable] -= diffPlusHm;
			b[lVar][nVariable] += diffPlusHm;
		}
		if ( !nodeInfo->firstPoint ) {
			a[lVar][nVariable] -= diffMinusH;
			c[lVar][nVariable] += diffMinusH;
		}

/*		a[lVar][nVariable] -= diffPlusHm + diffMinusH;
		if ( !nodeInfo->lastPoint ) {
			b[lVar][nVariable] += diffPlusHm;
		}
		if ( !nodeInfo->firstPoint ) {
			c[lVar][nVariable] += diffMinusH;
		}
*/		
	}
}

Double TCountPremFlamePhys::DiffCorr( int nVariable, NodeInfoPtr nodeInfo )
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
	
	for ( i = 0; i < nSpecies; ++i ) {
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

int CountPremPhysPostIter( void *object )
{
	TCountPremFlamePhysPtr	flame = ( TCountPremFlamePhysPtr )object;
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
	VectorPtr	yLeftVec = currGrid->GetYLeft();
	VectorPtr	yRightVec = currGrid->GetYRight();
	Double		**y = yMat->mat;
	Double		*x = currGrid->GetX()->vec;
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
	int			variables = bt->GetNVariables();
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();
	
	for ( i = 0; i < nGridPoints; ++i ) {
		if ( flame->CheckSolution( y[i][fTemperature], &y[i][fFirstSpecies], nSpeciesInSystem ) ) {
			return 1;
		}
	}

	CountPremPhysUpdateLeftBoundary( flame );
	CountPremPhysUpdateRightBoundary( flame );

	if ( flame->GetSoot() ) {
		flame->GetSoot()->PostIter( flame );
	}

/*
// left boundary
	for ( int j = 0; j < variables; ++j ) {
		yLeft[j] = y[0][j];
	}
	yLeft[fVVelocity] = 0.0;

//	update properties
	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, temp[kPrev], Y[kPrev], pressure );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, temp[nGridPoints], Y[nGridPoints], pressure );
	
// right boundary 
	yRight[fUVelocity] = 1.0;
	yRight[fVVelocity] = yLast[fVVelocity] - hLast * strainRate * rhoRight * yRight[fUVelocity];
*/	

	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	
	flame->UpdateThermoProps();

	return 0;
}

#include "TofZ.h"

void TCountPremFlamePhys::UpdateDimensions( int len )
{
	T1DFlame::UpdateDimensions( len );
	fSolV->len = len;
	fSolU->len = len;
}

void TCountPremFlamePhys::UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
{
	int		nGridPoints = yMat->cols;
	Double	*V = fSolV->vec;
	Double	*U = fSolU->vec;
	Double	**y = yMat->mat;
	Double	*yLeft = yLeftVec->vec;
	Double	*yRight = yRightVec->vec;

	UpdateDimensions( nGridPoints );

	T1DFlame::UpdateSolution( yMat, yLeftVec, yRightVec );
	
	V[kPrev] = yLeft[fVVelocity];
	U[kPrev] = yLeft[fUVelocity];
	for ( int k = 0; k < nGridPoints; ++k ) {
		V[k] = y[k][fVVelocity];
		U[k] = y[k][fUVelocity];
	}
	V[nGridPoints] = yRight[fVVelocity];
	U[nGridPoints] = yRight[fUVelocity];
}

void TCountPremFlamePhys::UpdateSolutionOnePoint( Double *y, int gridPoint )
{
	int		i;
	int		nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();

	fSolU->vec[gridPoint] = y[fUVelocity];
	fSolV->vec[gridPoint] = y[fVVelocity];
	fSolTemp->vec[gridPoint] = y[fTemperature];
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fSolMassFracs->mat[gridPoint][i] = y[i+fFirstSpecies];
	}
}

int	TCountPremFlamePhys::GetOffsetVVelocity( void )
{
	return fVVelocity; 
}

int	TCountPremFlamePhys::GetOffsetUVelocity( void )
{
	return fUVelocity; 
}

int	TCountPremFlamePhys::GetOffsetTemperature( void )
{
	return fTemperature; 
}

int	TCountPremFlamePhys::GetOffsetFirstSpecies( void ) 
{
	return fFirstSpecies;
}

ConstStringArray TCountPremFlamePhys::GetVariableNames( void )
{
	return ( ConstStringArray )fVariableNames;
}

int TCountPremFlamePhys::GetVariablesWithoutSpecies( void )
{
	return fVariablesWithoutSpecies;
}

void TCountPremFlamePhys::SetInitialValues( TInputDataPtr inp, StartProfilePtr sp )
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
	if ( !xIn ) FatalError( "memory allocation of TCountPremFlamePhys failed" );
	Double				*yIn =  new Double[gridPointsIn];
	if ( !yIn ) FatalError( "memory allocation of TCountPremFlamePhys failed" );
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
	Double				*rho = GetProperties()->GetDensity()->vec;
	FILE				*fp;
	
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

// set left hand side of the grid to zero
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
	for ( j = 0; j < variables; ++j ) {
		yLeft[j] = y[0][j];
	}
	yLeft[fVVelocity] = 0.0;

//	update properties
	UpdateSolution( yMat, yLeftVec, yRightVec );
	SetFlameNode( kPrev );
	ComputeProperties( fFlameNode, temp[kPrev], Y[kPrev], GetPressure() );
	SetFlameNode( nGridPoints );
	ComputeProperties( fFlameNode, temp[nGridPoints], Y[nGridPoints], GetPressure() );

	yRight[fUVelocity] = 1.0;
	yRight[fVVelocity] = y[nGridPoints-1][fVVelocity] - ( bt->GetRight() - locX[nGridPoints-1] ) 
						* GetStrainRate() * rho[nGridPoints-1] * yRight[fUVelocity];
						
	if ( GetSpecies()->IsConstantLewisNumber() ) {
		CompLewisNumbers( GetSpecies()->GetLewisNumberFile() );
	}
	CountPremPhysPostIter( this );

	FreeSpline( theSpline );
	delete yIn;
	delete xIn;

	adapGrid->SetSolutionScaler();
	
	fp = GetOutfile( "initialguess", TFlame::kData );
	bt->PrintSolution( locX, y, GetVariableNames(), fp );
	fclose(fp);
}

Double TCountPremFlamePhys::SpeciesDiffusion( int nVariable, int speciesIndex, NodeInfoPtr nodeInfo )
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

void TCountPremFlamePhys::FillJacSpeciesDiffusion( int nVariable, int speciesIndex, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign )
{
// fills the jacobian with     constCoeff * d/dy ( rho * diffusivity * df/dy)

	if ( sign == kNegative ) {
		constCoeff *= -1.0;
	}

	Double	diff = fFlameNode->diffusivity[speciesIndex];
	Double	diffNext = fFlameNode->diffusivityNext[speciesIndex];
	Double	diffPrev = fFlameNode->diffusivityPrev[speciesIndex];
	Double	*mixDensity = fFlameNode->mixDensity;
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	diffPlusHm = constCoeff * ( diff * mixDensity[kCurr] + diffNext * mixDensity[kNext] ) * hm;
	Double	diffMinusH = constCoeff * ( diffPrev * mixDensity[kPrev] + diff * mixDensity[kCurr] ) * h;
	Double	*molarMass = GetSpecies()->GetMolarMass()->vec;
	int		nSpeciesInSystem = GetSpecies()->GetNSpeciesInSystem();

	nodeInfo->a[nVariable][nVariable] -= ( diffPlusHm + diffMinusH );
/*	nodeInfo->a[fTemperature][nVariable] -= constCoeff * mixDensity[kCurr] / temp[kCurr] 
		* diff * ( hm * ( yNext - y ) - h * ( y - yPrev ) );
	for ( int j = 0; j < nSpeciesInSystem; ++j ) {
		nodeInfo->a[fFirstSpecies + j][nVariable] -= constCoeff * mixDensity[kCurr] 
			* mixMolarMass[kCurr] / molarMass[j] * diff 
			* ( hm * ( yNext - y ) - h * ( y - yPrev ) ); 
	} */

	if ( !nodeInfo->lastPoint ) {
		nodeInfo->b[nVariable][nVariable] += diffPlusHm;
/*		nodeInfo->b[fTemperature][nVariable] -= constCoeff * mixDensity[kNext] / temp[kNext] 
				* diffNext * hm * ( yNext - y );
		for ( j = 0; j < nSpeciesInSystem; ++j ) {
			nodeInfo->b[fFirstSpecies + j][nVariable] -= constCoeff * mixDensity[kNext] 
					* mixMolarMass[kNext] / molarMass[j] * diffNext 
					* hm * ( yNext - y ); 
		} */
	}
	if ( !nodeInfo->firstPoint ) {
		nodeInfo->c[nVariable][nVariable] += diffMinusH;
/*		nodeInfo->c[fTemperature][nVariable] += constCoeff * mixDensity[kPrev] / temp[kPrev] 
				* diffPrev * h * ( y - yPrev );
		for ( j = 0; j < nSpeciesInSystem; ++j ) {
			nodeInfo->c[fFirstSpecies + j][nVariable] += constCoeff * mixDensity[kPrev] 
					* mixMolarMass[kPrev] / molarMass[j] * diffPrev 
					* h * ( y - yPrev ); 
		} */
	}
}

void CountPremPhysOutput( void *object, FILE *fp, char* tail )
{
	TCountPremFlamePhysPtr	flame = ( TCountPremFlamePhysPtr )object;
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
	time_t			theDate;
	char			buffer[80];
	ConstStringArray	varNames = flame->GetVariableNames();
	char			**names = species->GetNames();
	VectorPtr 		etaVec = NewVector( gridPoints + 2 );
	Double			*eta = etaVec->vec;
	Flag			fpOpen = FALSE;
	
	if ( !fp ) {
		fpOpen = TRUE;
		fp = flame->GetOutputFile( NULL, tail, TFlame::kNone );
	}
	
// write header
	fprintf( fp, "header\n\n" );

	fprintf( fp, "title = \"planar counterflow premixed flame\"\n" );
	fprintf( fp, "author = \"Marcus Bollig, ITM Aachen, Germany\"\n" );
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
	
	fprintf( fp, "\nunburnt\n" );
	fprintf( fp, "begin\n" );
	fprintf( fp, "\tTemperature = %g [K]\n", yRight[tempOffset] );
	for ( i = 0; i < nOfSpecies; ++i ) {
		if ( fabs( massFracs[gridPoints][i] ) > 1.0e-10 ) {
			fprintf( fp, "\tMassfraction-%s = %g\n", names[i], massFracs[gridPoints][i] );
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
	
	fprintf( fp, "trailer\n" );
	if ( nOfEquations < nOfVariables) {
		fprintf( fp, "number of converged equations is %d\n", nOfEquations );
	}
	
	DisposeVector( etaVec );
	if ( fpOpen ) {
		fclose( fp );
	}
}

void SetCountPremPhysNodeInfo( int k, void *object )
{
	TCountPremFlamePhysPtr	flame = ( TCountPremFlamePhysPtr )object;
	
	flame->SetFlameNode( k );
}

void CountPremPhysPostConv( void *object )
{
	TCountPremFlamePhysPtr	flame = ( TCountPremFlamePhysPtr )object;
	
	flame->PostConvergence( object );
}

ConstStringArray GetCountPremPhysVarNames( void *object )
{
	TCountPremFlamePhysPtr	flame = ( TCountPremFlamePhysPtr )object;
	
	return flame->GetVariableNames();
}

FILE *TCountPremFlamePhys::GetOutputFile( char *head, char *tail, FileType type )
{
	int				fuelIndex = GetFuelIndex();
	char			*name = new char[64];
	FILE			*fp;
	char			**speciesNames = fSpecies->GetNames();
	int				tOxidizer = ( int ) fSolTemp->vec[fSolTemp->len];
	Double			phi = GetPhi();
		
	sprintf( name, "%s%s%.8s_p%.2dphi%.1d-%.4da%.5dt%.4d%s"
					, ( head ) ? head : "", ( head ) ? "_" : ""
					, speciesNames[fuelIndex]
					, ( int ) floor( GetPressure() * 1.0e-5 + 0.5 )	// in [bar]
					, ( int ) floor( phi ) , ( int ) floor(( phi - floor(phi)) * 1.e4 )
					, ( int ) floor( GetStrainRate() + 0.5 )			// in [1/s]
					, ( int )( tOxidizer ) 						// in [K]
					, ( tail ) ? tail : "" );

	fp = GetOutfile( name, type );
	delete name;

	return fp;
}

int	TCountPremFlamePhys::GetOffsetMixFrac( void )
{
	cerr << "#error: class has no member fMixtureFraction" << NEWL;
	exit( 2 );
	return 0; 
}

void CountPremPhysUpdateLeftBoundary( void  *object )
{
	TCountPremFlamePhysPtr	flame = ( TCountPremFlamePhysPtr )object;
	int				fVVelocity = flame->GetOffsetVVelocity();
	TNewtonPtr		bt = flame->GetSolver()->bt;
	int				nVariables = bt->GetNVariables();
	TGridPtr 		currGrid = bt->GetGrid()->GetCurrentGrid();
	MatrixPtr		yMat = currGrid->GetY();
	Double			**y = yMat->mat;
	VectorPtr		yLeftVec = currGrid->GetYLeft();
	Double			*yLeft = yLeftVec->vec;
	Double			pressure = flame->GetPressure();

	for ( int i = 0; i < nVariables; ++i ) {
		yLeft[i] = y[0][i];
	}
	yLeft[fVVelocity] = 0.0;

	flame->UpdateSolutionOnePoint( yLeft, kPrev );
	flame->SetFlameNode( kPrev );
	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], pressure );
}

void CountPremPhysUpdateRightBoundary( void *object )
{
	TCountPremFlamePhysPtr	flame = ( TCountPremFlamePhysPtr )object;
	int				fVVelocity = flame->GetOffsetVVelocity();
	int				fUVelocity = flame->GetOffsetUVelocity();
	TNewtonPtr		bt = flame->GetSolver()->bt;
	TGridPtr 		currGrid = bt->GetGrid()->GetCurrentGrid();
	int				nGridPoints = currGrid->GetNGridPoints();
	Double			*x = currGrid->GetX()->vec;
	MatrixPtr		yMat = currGrid->GetY();
	Double			**y = yMat->mat;
	Double			*yLast = y[nGridPoints-1];
	VectorPtr		yRightVec = currGrid->GetYRight();
	Double			*yRight = yRightVec->vec;
	Double			hLast = bt->GetRight() - x[nGridPoints-1];
	T1DPropertiesPtr	properties = flame->GetProperties();
	Double			&rhoRight = properties->GetDensity()->vec[nGridPoints];
	Double			pressure = flame->GetPressure();
	Double			strainRate = flame->GetStrainRate();

	yRight[fUVelocity] = 1.0;
	yRight[fVVelocity] = yLast[fVVelocity] - hLast * strainRate * rhoRight * yRight[fUVelocity];

	flame->UpdateSolutionOnePoint( yRight, nGridPoints );
	flame->SetFlameNode( nGridPoints );
	flame->ComputeProperties( flame->fFlameNode, flame->fFlameNode->temp[kCurr]
							, flame->fFlameNode->Y[kCurr], pressure );
}
