//#define FOOL_SOFTBENCH( x ) 

#include "FlameMaster.h"
#include "TFlameSheet.h"

void TFlameSheet::InitTFlameSheet( void )
{
	TBVPSolverPtr	solver = GetSolver();
	TNewtonPtr		bt = solver->bt;
	TGridPtr		fine = bt->GetGrid()->GetFine();
	TGridPtr		coarse = bt->GetGrid()->GetCoarse();
	int				nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();
	StartProfilePtr sp = NULL;
	int				maxGridPoints = bt->GetMaxGridPoints();

	fVariableNames = new String[fVariablesWithoutSpecies + nSpeciesInSystem];

	fVariableNames[0] = new char[2];
	strcpy( fVariableNames[0], "f" );
	fVariableNames[1] = new char[3];
	strcpy( fVariableNames[1], "f'" );
	fVariableNames[2] = new char[2];
	strcpy( fVariableNames[2], "Z" );
	fVariableNames[3] = new char[2];
	strcpy( fVariableNames[3], "T" );
	
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

	bt->SetUtFuncs( FlameSheetJacRest, FlameSheetJacRest, FlameSheetJacRest
					, FlameSheetRHSRest, FlameSheetRHSRest, FlameSheetRHSRest 
					, FlameSheetOutput, FlameSheetPostIter, SetFlameSheetNodeInfo
					, FlameSheetPostConv, GetFlameSheetVarNames );
	SetInitialBC( fine, fInputData );
	SetInitialBC( coarse, fInputData );
	SetInitialValues( fInputData, sp ); // initial values of coarse are set during gridgeneration
	bt->SetNOfEquations( 3 );
	fMassFraction = new TMassFraction( this );
	if ( !fMassFraction ) FatalError( "memory allocation of TMassFraction failed" );
	UpdateSolution( fine->GetY(), fine->GetYLeft(), fine->GetYRight() );	
}

TFlameSheet::~TFlameSheet( void )
{ 
	int	nSpeciesInSystem = fSpecies->GetNSpeciesInSystem();

	delete fMassFraction;

	for ( int i = 0; i < nSpeciesInSystem+fVariablesWithoutSpecies; ++i ) {
		delete fVariableNames[i];
	}
	delete fVariableNames;
}

void FlameSheetJacRest( void *object, NodeInfoPtr nodeInfo )
{
	TFlameSheetPtr	flame = ( TFlameSheetPtr )object;
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
    Double  hnenn = h * hm * ( hm + h );
	Double	**a = nodeInfo->a;
	Double	*y = nodeInfo->y;
	Double	strainRate = flame->GetStrainRate();
	Double	mixDensityInf = flame->fFlameNode->rhoInf;
	Double	mixViscosityInf = flame->fFlameNode->viscosityInf;
	Double	fMixFracDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );

// first fill blocks with linear terms
// first equation ( mass )
	FillJacFirstDerivUp( fVVelocity, fVVelocity, nodeInfo );

	a[fUVelocity][fVVelocity] -= hnenn;

// second equation ( momentum )
	FillJacNonlinearConvectCentral( fVVelocity, fUVelocity, nodeInfo );
	flame->FillJacMassDiffusion( fUVelocity, fUVelocity, nodeInfo );
	a[fUVelocity][fUVelocity] -= 2.0 * y[fUVelocity] * hnenn;
	
// third equation ( energy )
	FillJacNonlinearConvectCentral( fVVelocity, fMixFrac, nodeInfo );
	flame->FillJacMixFracDiffusion( fMixFrac, fMixFrac, fMixFracDiffCoeff, nodeInfo );
}

void FlameSheetRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode /*rhsMode*/ )
{
	TFlameSheetPtr	flame = ( TFlameSheetPtr )object;
	int 	fMixFrac = flame->GetOffsetMixFrac();
	int 	fUVelocity = flame->GetOffsetUVelocity();
	int 	fVVelocity = flame->GetOffsetVVelocity();
	int		eqLoop;
	int		M = nodeInfo->nOfEquations;
    Double  h = nodeInfo->h;
    Double  hm = nodeInfo->hm;
    Double  hnenn = hm * h * ( hm + h );
	Double	*rhs = nodeInfo->rhs;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	strainRate = flame->GetStrainRate();
	Double	mixDensity = *flame->fFlameNode->mixDensity;
	Double	mixDensityInf = flame->fFlameNode->rhoInf;
	Double	mixViscosityInf = flame->fFlameNode->viscosityInf;
	Double	mixFracDiffCoeff = 1.0 / ( mixDensityInf * mixViscosityInf );

	// fill RHS with linear terms
	rhs[fVVelocity] += FirstDerivUpwind( y[fVVelocity], yPrev[fVVelocity], hm );
	rhs[fVVelocity] -= y[fUVelocity];

// second equation ( momentum )
	rhs[fUVelocity] += NonlinearConvectCentral( y[fVVelocity], yPrev[fUVelocity], y[fUVelocity], yNext[fUVelocity], hm, h );
	rhs[fUVelocity] += flame->SecondDerivMassDiffusion( fUVelocity, nodeInfo );
	rhs[fUVelocity] -= y[fUVelocity] * y[fUVelocity];
	rhs[fUVelocity] += mixDensityInf / mixDensity;
	
// third equation ( energy )
	rhs[fMixFrac] += NonlinearConvectCentral( y[fVVelocity], yPrev[fMixFrac], y[fMixFrac], yNext[fMixFrac], hm, h );
	rhs[fMixFrac] += mixFracDiffCoeff * flame->SecondDerivMixFracDiffusion( fMixFrac, nodeInfo );
		
	// fill RHS with nonlinear terms

	for ( eqLoop = 0; eqLoop < M; ++eqLoop ){
		rhs[eqLoop] *= - hnenn;
	}
}

void TMassFraction::InitTMassFraction( T1DFlamePtr flame )
{ 
	TSpeciesPtr 	species = flame->GetSpecies();
	TInputDataPtr	input = flame->GetInputData();
	TBVPSolverPtr	solver = flame->GetSolver();
	SpeciesPtr		speciesStruct = input->GetSpecies();
	ReactionPtr		reactionStruct = input->GetReactions();
	int				firstSpecies = flame->GetOffsetFirstSpecies();
	int				nSpeciesInSystem = species->GetNSpeciesInSystem();
//	MatrixPtr		massFracs = flame->GetMassFracs();
	TGridPtr		fineGrid = solver->bt->GetGrid()->GetFine();
	Double			*yLeft = fineGrid->GetYLeft()->vec;
	Double			*yRight = fineGrid->GetYRight()->vec;
	Double			*molarMass = NULL;
	Double			molarMassFuel = 0.0;
	Double			molarMassOx = 0.0;
	Double			molarMassH2O = 0.0;
	Double			molarMassCO2 = 0.0;
	
	
	if ( flame->GetNFuels() > 1 ) { // ( hp )
		fprintf( stderr, "###error: function 'TMassFraction::InitTMassFraction' is not valid for more than one fuel\n" );
		exit( 2 );
	}
 
 	indexFuel = input->fFuelIndex->vec[0];
	indexOx = input->fOxIndex;

	indexH2O = input->fH2OIndex;
	indexCO2 = input->fCO2Index;
	
	if ( indexFuel >= nSpeciesInSystem || indexOx >= nSpeciesInSystem 
			|| indexH2O >= nSpeciesInSystem || indexCO2 >= nSpeciesInSystem ) {
		cerr << "error: participants of the global reaction can't be steady state" << NEWL;
		cerr << "indexfuel = " << indexFuel << NEWL;
		cerr << "indexOx = " << indexOx << NEWL;
		cerr << "indexH2O = " << indexH2O << NEWL;
		cerr << "indexCO2 = " << indexCO2 << NEWL;
		cerr << "nSpeciesInSystem = " << nSpeciesInSystem << NEWL;
		exit(2);
	}
	
	// errorchecking
	if ( indexH2O == -1 && indexCO2 == -1 ) {
		fprintf( stderr, "###warning: no CO2 and no H2O\n" );
	}
	
	if ( indexH2O == -1 ) {
		//cerr << "error: there is no species with the name 'H2O'" << NEWL;
		//		exit(2);
	}
	if ( indexCO2 == -1 ) {
		//cerr << "error: there is no species with the name 'CO2'" << NEWL;
		//		exit(2);
	}
	
	molarMass = species->GetMolarMass()->vec;
	molarMassFuel = molarMass[input->fFuelIndex->vec[0]];
	molarMassOx = molarMass[input->fOxIndex];
	if ( indexH2O == -1 ) {
		molarMassH2O = 0.0;
	}
	else {
		molarMassH2O = molarMass[input->fH2OIndex];
	}
	
	if ( indexCO2 == -1 ) {
		molarMassCO2 = 0.0;
	}
	else {
		molarMassCO2 = molarMass[input->fCO2Index];
	}

	aFH = input->FindAtomIndex( "H" );
	aFH = speciesStruct[indexFuel].composition->vec[aFH];

	aFC = input->FindAtomIndex( "C" );
	aFC = speciesStruct[indexFuel].composition->vec[aFC];
	
	nuOx = flame->GetNu( input->fGlobalReaction, "O2" );
	nuFuel = flame->GetNu( input->fGlobalReaction, speciesStruct[indexFuel].name );
	//	nuOx = input->FindSpeciesCoefficient( indexOx, 0 );
	if ( nuOx == -1 ) {
		cerr << "error: there is no oxidizer in reaction no. '0'" << NEWL;
		exit(2);
	}
	//	nuFuel = input->FindSpeciesCoefficient( indexFuel, 0 );
	if ( nuFuel == -1 ) {
		cerr << "error: there is no fuel in reaction no. '0'" << NEWL;
		exit(2);
	}
	
	fNu = nuOx * molarMassOx / ( nuFuel * molarMassFuel );
	if ( yRight[indexOx+firstSpecies] == 0 ) {
		fYF1 = &yRight[indexFuel+firstSpecies];
		fYOx2 = &yLeft[indexOx+firstSpecies];
		fYF2 = &yLeft[indexFuel+firstSpecies];
		fYOx1 = &yRight[indexOx+firstSpecies];
	}
	else {
		fYF1 = &yLeft[indexFuel+firstSpecies];
		fYOx2 = &yRight[indexOx+firstSpecies];
		fYF2 = &yRight[indexFuel+firstSpecies];
		fYOx1 = &yLeft[indexOx+firstSpecies];
	}
	
//	fprintf( stderr, "fYF1 = %g\tfYOx2 = %g\tfYF2 = %g\tfYOx1 = %g\t\n", fYF1, fYOx2, fYF2, fYOx1 );
	
//	fZStoe = 1.0 / ( 1.0 + fNu * *fYF1 / *fYOx2 );
	fZStoe = ( *fYOx2 - fNu * *fYF2 ) / ( fNu * ( *fYF1 - *fYF2 ) + *fYOx2 - *fYOx1 );
	
	fYH2OSt = *fYF1 * fZStoe * 0.5 * aFH * molarMassH2O / molarMassFuel;
	fYCO2St = *fYF1 * fZStoe * aFC * molarMassCO2 / molarMassFuel;
}

void TMassFraction::PrintMassFraction( void )
{
	FILE *fp = NULL;
	
	if ( !( fp = fopen( "massfraction", "w" ) ) ) { 
		cerr << "#warning: unable to open file 'massfraction'" << NEWL;
		return;
	}
	fprintf( fp, "indexFuel = %d\n", indexFuel );
	fprintf( fp, "indexOx = %d\n", indexOx );
	fprintf( fp, "indexH2O = %d\n", indexH2O );
	fprintf( fp, "indexCO2 = %d\n", indexCO2 );
	
	fprintf( fp, "coefficient of fuel is %g\n", nuFuel );
	fprintf( fp, "coefficient of oxidizer is %g\n", nuOx );

	fprintf( fp, "fuel contains %d atoms of type C\n", aFC );
	fprintf( fp, "fuel contains %d atoms of type H\n", aFH );

	fprintf( fp, "massfraction of fuel at boundary is %g\n", *fYF1 );
	fprintf( fp, "massfraction of oxidizer at boundary is %g\n", *fYOx2 );
	fprintf( fp, "nu = %g\n", fNu );
	fprintf( fp, "Zstoi = %g\n", fZStoe );
	fprintf( fp, "YH2OSt = %g\n", fYH2OSt );
	fprintf( fp, "YCO2St = %g\n", fYCO2St );
	fclose( fp );
}

void TMassFraction::SetMassFractionsBC( int mixFracIndex, int firstSpeciesIndex, Double *yLeft, Double *yRight )
{
	Double	denom = fNu * fYF1[kCurr] + fYOx2[kCurr] - fNu * fYF2[kCurr] - fYOx1[kCurr];

	if ( yLeft ) {
		yLeft[mixFracIndex] = ( fNu * yLeft[indexFuel+firstSpeciesIndex] - yLeft[indexOx+firstSpeciesIndex] + fYOx2[kCurr] - fNu * fYF2[kCurr] )
								/ denom;
	}
	if ( yRight ) {
		yRight[mixFracIndex] = ( fNu * yRight[indexFuel+firstSpeciesIndex] - yRight[indexOx+firstSpeciesIndex] + fYOx2[kCurr] - fNu * fYF2[kCurr] )
								/ denom;
	}
}

int FlameSheetPostIter( void *object )
{
	TFlameSheetPtr	flame = ( TFlameSheetPtr ) object;
	int			i;
	int			fUVelocity = flame->GetOffsetUVelocity();
	int			fVVelocity = flame->GetOffsetVVelocity();
	int			fFirstSpecies = flame->GetOffsetFirstSpecies();
	int			fTemperature = flame->GetOffsetTemperature();
	TNewtonPtr	bt = flame->GetSolver()->bt;
	NodeInfoPtr nodeInfo = bt->GetNodeInfo();
	TGridPtr 	currGrid = bt->GetGrid()->GetFine();
	int			nGridPoints = currGrid->GetNGridPoints();
	int			*bcFlagLeft = currGrid->GetBcFlagLeft();
	int			mixtureSpecificationLeft = flame->GetMixtureSpecificationLeft();
	MatrixPtr	yMat = currGrid->GetY();
	VectorPtr	yLeftVec = currGrid->GetYLeft();
	VectorPtr	yRightVec = currGrid->GetYRight();
	Double		**y = yMat->mat;
	Double		*x = currGrid->GetX()->vec;
	Double		*yLeft = yLeftVec->vec;
	Double		*yRight = yRightVec->vec;
	Double		*yLast = y[nGridPoints-1];
	Double		hLast = bt->GetRight() - x[nGridPoints-1];
	Double		*bcLeft = currGrid->GetBcRight()->vec;
	Double		hFirst = x[0] - bt->GetLeft();
	T1DSpeciesPtr	species = flame->GetSpecies();
	T1DPropertiesPtr	properties = flame->GetProperties();
	Double		*diffusivity = species->GetDiffusivity()->mat[0];
	Double		rhoLeft = properties->GetDensity()->vec[-1];
	Double		rhoRight = properties->GetDensity()->vec[nGridPoints];
	Double		*temp = flame->GetTemperature()->vec;
	Double		**Y = flame->GetMassFracs()->mat;
	Double		pressure = flame->GetPressure();
	int			nSpeciesInSystem = species->GetNSpeciesInSystem();
	
// first set temperature and massfractions 
	for ( i = 0; i < nGridPoints; ++i ) {
		bt->SetNodeInfo( flame, i );		
		flame->fMassFraction->UpdateMixFracDependence( flame, nodeInfo );
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
	if ( mixtureSpecificationLeft == kMassFlux ) {
		Double 		coeff;
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			if ( bcFlagLeft[fFirstSpecies] == kDirichlet ) {
				coeff = diffusivity[i] / yLeft[fVVelocity] / hFirst;
				yLeft[fFirstSpecies+i] = ( bcLeft[fFirstSpecies+i] + coeff ) / ( 1.0 + coeff );
			}
			else {
				cerr << "error: i can't handle boundary condition of kind " << bcFlagLeft[fFirstSpecies] << "for species no. " << i << NEWL;
			}
		}
	}
	
	rhoLeft = properties->GetDensity()->vec[-1];
	rhoRight = properties->GetDensity()->vec[nGridPoints];
	yLeft[fUVelocity] = sqrt( rhoRight / rhoLeft );
	//	yLeft[fVVelocity] = yFirst[fVVelocity] - hFirst * yLeft[fUVelocity];
	
// right boundary 
	yRight[fUVelocity] = 1.0;
	yRight[fVVelocity] = yLast[fVVelocity] + hLast * yRight[fUVelocity];

	flame->UpdateSolution( yMat, yLeftVec, yRightVec );
	
	flame->UpdateThermoProps();

	return 0;
}

#include "TofZ.h"

void TMassFraction::UpdateMixFracDependence( T1DFlamePtr flame, NodeInfoPtr nodeInfo, Flag justTemperature )
{
	static int	minNOfSpecies = 4;
	int			firstSpeciesOffset = flame->GetOffsetFirstSpecies();
	int			temperatureOffset = flame->GetOffsetTemperature();
	int			mixFracOffset = flame->GetOffsetMixFrac();
//	Double		nuCO2 = nuFuel * aFC;
//	Double		nuH2O = 0.5 * nuFuel * aFH;
	Double		*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double		*enthalpy = flame->GetFlameNode()->enthalpy;
	Double		cpMix = *flame->GetFlameNode()->mixHeatCapacity;
	Double		Tu;
	Double		*yLeft = nodeInfo->bcLeft; 
	Double		*yRight = nodeInfo->bcRight;
	Double		*y = nodeInfo->y;
	TofZParams	tOfZ;
		
	fZStoe = 1.0 / ( 1.0 + fNu * *fYF1 / *fYOx2 );
	if ( indexH2O >= 0 ) fYH2OSt = *fYF1 * fZStoe * 0.5 * aFH * molarMass[indexH2O] / molarMass[indexFuel];
	if ( indexCO2 >= 0 ) fYCO2St = *fYF1 * fZStoe * aFC * molarMass[indexCO2] / molarMass[indexFuel];

	if ( temperatureOffset < nodeInfo->nOfEquations ) {
		return;
	}
	
	get_coefficients( &tOfZ );
	
/*	Q = nuFuel * molarMass[indexFuel] * enthalpy[indexFuel] */
/*			+ nuOx * molarMass[indexOx] * enthalpy[indexOx]*/
/*			- nuCO2 * molarMass[indexCO2] * enthalpy[indexCO2]*/
/*			- nuH2O * molarMass[indexH2O] * enthalpy[indexH2O];*/

	Tu = yRight[temperatureOffset] - y[mixFracOffset] * ( yRight[temperatureOffset] - yLeft[temperatureOffset] );

	if ( y[mixFracOffset] < fZStoe ) {
		set_coefficients( &tOfZ );
		y[temperatureOffset] = temperature( y[mixFracOffset] );
		y[temperatureOffset] += Tu - 300.0;
		//	y[temperatureOffset] = Tu + Q * *fYF1 * y[mixFracOffset] / cpMix / nuFuel / molarMass[indexFuel];
		if ( !justTemperature ) {
			y[firstSpeciesOffset+indexFuel] = 0.0;
			y[firstSpeciesOffset+indexOx] = *fYOx2 * ( 1.0 - y[mixFracOffset] / fZStoe );
			y[firstSpeciesOffset+indexCO2] = fYCO2St * y[mixFracOffset] / fZStoe;
			y[firstSpeciesOffset+indexH2O] = fYH2OSt * y[mixFracOffset] / fZStoe;
		}
	}
	else {
		set_coefficients( &tOfZ );
		y[temperatureOffset] = temperature( y[mixFracOffset] );
		y[temperatureOffset] += Tu - 300.0;
		//	y[temperatureOffset] = Tu + Q * *fYOx2 * ( 1.0 - y[mixFracOffset] ) / cpMix / nuOx / molarMass[indexOx];
		if ( !justTemperature ) {
			y[firstSpeciesOffset+indexFuel] = *fYF1 * ( y[mixFracOffset] - fZStoe ) / ( 1.0 - fZStoe );
			y[firstSpeciesOffset+indexOx] = 0.0;
			y[firstSpeciesOffset+indexCO2] = fYCO2St * ( 1.0 - y[mixFracOffset] ) / ( 1.0 - fZStoe );
			y[firstSpeciesOffset+indexH2O] = fYH2OSt * ( 1.0 - y[mixFracOffset] ) / ( 1.0 - fZStoe );
		}
	}
	y[temperatureOffset] = MAX( 10.0, y[temperatureOffset] );
	
	int nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	if ( !justTemperature ) {
		if ( nSpeciesInSystem > minNOfSpecies ) {
			Double 	massFractRest = ( 1.0 - y[firstSpeciesOffset+indexFuel] - y[firstSpeciesOffset+indexOx]
									 	- y[firstSpeciesOffset+indexCO2] - y[firstSpeciesOffset+indexH2O] ) 
										/ ( nSpeciesInSystem - minNOfSpecies );
			for ( int i = 0; i < nSpeciesInSystem; ++i ) {
				if ( i != indexFuel && i != indexOx
						&& i != indexCO2 && i != indexH2O ) {
					y[firstSpeciesOffset+i] = MAX( 0.0, massFractRest );
				}
			}
		}
	}
}

void TFlameSheet::UpdateDimensions( int len )
{
	T1DFlame::UpdateDimensions( len );
	fSolV->len = len;
	fSolU->len = len;
	fSolMixFrac->len = len;
}

void TFlameSheet::UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
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

int	TFlameSheet::GetOffsetVVelocity( void )
{
	return fVVelocity; 
}

int	TFlameSheet::GetOffsetUVelocity( void )
{
	return fUVelocity; 
}

int	TFlameSheet::GetOffsetTemperature( void )
{
	return fTemperature; 
}

int	TFlameSheet::GetOffsetMixFrac( void )
{ 
	return fMixFrac;
}

int	TFlameSheet::GetOffsetFirstSpecies( void ) 
{
	return fFirstSpecies;
}

ConstStringArray TFlameSheet::GetVariableNames( void )
{
	return ( ConstStringArray ) fVariableNames;
}

int TFlameSheet::GetVariablesWithoutSpecies( void )
{
	return fVariablesWithoutSpecies;
}

void TFlameSheet::SetInitialValues( TInputDataPtr inp, StartProfilePtr /*sp*/ )
{
	int 		i, j;
	TGridPtr	grid = fSolver->bt->GetGrid()->GetFine();
	int			variables = grid->GetY()->rows;
	int			nGridPoints = grid->GetNGridPoints();
	Double 		*yLeft = grid->GetYLeft()->vec;
	Double 		*yRight = grid->GetYRight()->vec;
	Double		left = inp->fLeft;
	Double		right = inp->fRight;
	Double		*locX = grid->GetX()->vec;
	Double		**y = grid->GetY()->mat;
	
	for ( i = 0; i < nGridPoints; ++i ) {
		for ( j = 0; j < variables; ++j ) {
			y[i][j] = yLeft[j] + ( yRight[j] - yLeft[j] ) / ( right - left ) * locX[i];
		}
	}
}

void FlameSheetOutput( void */*object*/, FILE */*fp*/, char* /*tail*/ )
{
	cerr << "#warning: try to write flameletfile, but nothing happens" << NEWL;
}

void SetFlameSheetNodeInfo( int k, void *object )
{
	TFlameSheetPtr	flame = ( TFlameSheetPtr )object;
	
	flame->SetFlameNode( k );
}

void FlameSheetPostConv( void *object )
{
	TFlameSheetPtr	flame = ( TFlameSheetPtr )object;
	
	flame->PostConvergence( object );
}

ConstStringArray GetFlameSheetVarNames( void *object )
{
	TFlameSheetPtr	flame = ( TFlameSheetPtr )object;
	
	return flame->GetVariableNames();
}

void TFlameSheet::FillJacMassDiffusion( int nVariable, int nEquation, NodeInfoPtr nodeInfo, Flag sign )
{
// fills the jacobian with     1/( rho * mu )ref * d/dy ( rho * mu * df/dy)

	Double	*density = fFlameNode->mixDensity;
	Double	*viscosity = fFlameNode->mixViscosity;
	Double	constCoeff = 1.0 / ( fFlameNode->rhoInf * fFlameNode->viscosityInf );
	Double	diffPlus = constCoeff * ( density[kCurr] * viscosity[kCurr] + density[kNext] * viscosity[kNext] );
	Double	diffMinus = constCoeff * ( density[kPrev] * viscosity[kPrev] + density[kCurr] * viscosity[kCurr] );
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;

	if ( sign == kPositive ) {
		nodeInfo->a[nVariable][nEquation] -= ( hm * diffPlus + h * diffMinus );
		if ( !nodeInfo->lastPoint ) {
			nodeInfo->b[nVariable][nEquation] += hm * diffPlus;
		}
		if ( !nodeInfo->firstPoint ) {
			nodeInfo->c[nVariable][nEquation] += h * diffMinus;
		}
	}
	else {
		nodeInfo->a[nVariable][nEquation] += ( hm * diffPlus + h * diffMinus );
		if ( !nodeInfo->lastPoint ) {
			nodeInfo->b[nVariable][nEquation] -= hm * diffPlus;
		}
		if ( !nodeInfo->firstPoint ) {
			nodeInfo->c[nVariable][nEquation] -= h * diffMinus;
		}
	}
}

Double TFlameSheet::SecondDerivMassDiffusion( int nVariable, NodeInfoPtr nodeInfo )
{
	Double	yPrev = nodeInfo->yPrev[nVariable];
	Double	y = nodeInfo->y[nVariable];
	Double	yNext = nodeInfo->yNext[nVariable];
	Double	*density = fFlameNode->mixDensity;
	Double	*viscosity = fFlameNode->mixViscosity;
	Double	constCoeff = 1.0 / ( fFlameNode->rhoInf * fFlameNode->viscosityInf );
	Double	diffPlus = constCoeff * ( density[kCurr] * viscosity[kCurr] + density[kNext] * viscosity[kNext] );
	Double	diffMinus = constCoeff * ( density[kPrev] * viscosity[kPrev] + density[kCurr] * viscosity[kCurr] );
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	
	return ( diffPlus * hm * yNext - ( diffPlus * hm + diffMinus * h ) * y + diffMinus * h * yPrev ) 
				/ nodeInfo->hnenn;
}

FILE *TFlameSheet::GetOutputFile( char *head, char *tail, FileType type )
{
	int				fuelIndex = GetFuelIndex();
	char			*name = new char[64];
	FILE			*fp;
	char			**speciesNames = fSpecies->GetNames();
	int				tFuel = ( int ) fSolTemp->vec[fSolTemp->len];
	int				tOxidizer = ( int ) fSolTemp->vec[kPrev];
		
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
