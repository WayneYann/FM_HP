#include "FlameMaster.h"

#undef FULLDIFFUSION

// following makes use of liquid fuel easier, where the assumption of adiabatic droplet
// evaporation can lead to negative temperatures. Keeping liquid cp constant leads to smaller
// temperature change during evaporation
#define CONSTLOWTCP


// if source code for production rates is generated with ScanMan, it can be tested
// here and compared with the solution of the mechanism itself
// For using production rate files produced by ScanMan, define PRODRATEFILE
// If Fortran files should be used, define PRODRATEFILEF77, for C files undefine PRODRATEFILEF77
// Needs to be set also in TFlame.C
#undef PRODRATEFILE
#undef PRODRATEFILEF77

// optimization
#define OPTPROPS
#define OPTIMIZEDELTAI
#define OPTOMEGAD
#define OPTDTHERM


#ifdef PRODRATEFILE
void ComputeThermoData( Double *h, Double *cp, Double T );
#endif

#ifdef PRODRATEFILEF77
#	define PRODRATES prodrates_
#	define COMPTHERMODATA compthermodata_
#	define GETMOLARMASS getmolarmass_
extern "C" {
	void PRODRATES( Double *cdot, Double *w, Double *k
								, Double *c, Double *M, Double *temp, Double *pressure );
	void COMPTHERMODATA( Double *h, Double *cp, Double *T );
	void GETMOLARMASS( Double *MM );
}

#endif

void TSpecies::InitSpecies( TInputDataPtr input, TReactionPtr reaction )
{
	int 		i;	
   	int			nOfSpecies = input->GetCounter()->species;
	int			nOfReactions = input->GetCounter()->reactions;
	SpeciesPtr	species = input->GetSpecies();
   
	fNSteadyStates = input->GetCounter()->steadyStates;
	
	fTThermoLimit = 273.15;

	fNOfUsedReactions = NewIntVector( nOfSpecies );
	CalcUsedReactions( reaction );
	fUsedReactions = new IntVectorPtr[nOfSpecies];
	if ( !fUsedReactions ) FatalError( "memory allocation of TSpecies failed" );
	fNu = new VectorPtr[nOfSpecies];
	if ( !fNu ) FatalError( "memory allocation of TSpecies failed" );
	for ( i = 0; i < nOfSpecies; ++i ) {
		if ( fNOfUsedReactions->vec[i] ) {
			fUsedReactions[i] = NewIntVector( fNOfUsedReactions->vec[i] ); // last element is fUsedReactions[nOfSpecies]->vec[fNOfUsedReactions->vec[i]]
			fNu[i] = NewVector( fNOfUsedReactions->vec[i] ); // last element is fNu[nOfSpecies]->vec[fNOfUsedReactions->vec[i]]
		}
		else {
			fUsedReactions[i] = NewIntVector( 1 );
			fNu[i] = NewVector( 1 );
			fUsedReactions[i]->len = 0;
			fNu[i]->len = 0;
		}
	}

	fLewisNumberFile = input->fLewisNumberFile;
	fLewisNumber = NewVector( nOfSpecies );
	
	fNames = new String[nOfSpecies];
	if ( !fNames ) FatalError( "memory allocation of TSpecies failed" );
	fCoeffLow = NewMatrix( NOFCOEFFS, nOfSpecies, kColumnPointers ); // last element is fCoeffLow[nOfSpecies][7]
	fCoeffHigh = NewMatrix( NOFCOEFFS, nOfSpecies, kColumnPointers ); // last element is fCoeffHigh[nOfSpecies][7]
	fHCoeffLow = NewMatrix( NOFCOEFFS-1, nOfSpecies, kColumnPointers ); // last element is fCoeffLow[nOfSpecies][7]
	fHCoeffHigh = NewMatrix( NOFCOEFFS-1, nOfSpecies, kColumnPointers ); // last element is fCoeffHigh[nOfSpecies][7]
	fMolarMass = NewVector( nOfSpecies );
// don't use GetNSpeciesInSystem() before fKOverEps is defined
	fKOverEps = NewVector( nOfSpecies );
	fSigma = NewVector( nOfSpecies );
	fMuCoeff = NewVector( nOfSpecies );
	fDCoeff = NewMatrix( nOfSpecies, nOfSpecies, kColumnPointers ); // fDCoeff[i][j] is the same as fDCoeff[j][i]
	fOmegaCoeff = NewMatrix( nOfSpecies, nOfSpecies, kColumnPointers ); // fOmegaCoeff[i][j] is the same as fOmegaCoeff[j][i]
    fInvDij = NewMatrix( nOfSpecies, nOfSpecies, kColumnPointers );
	
	fIsSteadyState = new Flag[nOfSpecies];

	int nSpeciesIn = GetNSpeciesInSystem();

	sqrtMjOverMi = New2DArray( nSpeciesIn, nSpeciesIn );
	sqrtMiOverMjPlusOne = New2DArray( nSpeciesIn, nSpeciesIn );

	FillSpecies( input, reaction );
}

TSpecies::~TSpecies( void )
{
	int j;
	int nOfSpecies = GetNOfSpecies();
	
	Free2DArray( sqrtMiOverMjPlusOne );
	Free2DArray( sqrtMjOverMi );

	delete fIsSteadyState;

	for ( j = 0; j < nOfSpecies; ++j ) {
		delete fNames[j];
	}

	DisposeMatrix( fInvDij );
   	DisposeMatrix( fOmegaCoeff );
	DisposeMatrix( fDCoeff );
	DisposeVector( fMuCoeff );
	DisposeVector( fSigma );
	DisposeVector( fKOverEps );
	DisposeVector( fMolarMass );
	DisposeMatrix( fHCoeffHigh );
	DisposeMatrix( fHCoeffLow );
	DisposeMatrix( fCoeffHigh );
	DisposeMatrix( fCoeffLow );
	delete fNames;
	DisposeVector( fLewisNumber );

	for ( j = 0; j < nOfSpecies; ++j ) {
		DisposeVector( fNu[j] );
		DisposeIntVector( fUsedReactions[j] );
	}
	delete fNu;
	delete fUsedReactions;

	DisposeIntVector( fNOfUsedReactions );
}

void TSpecies::FillSpecies( TInputDataPtr input, TReactionPtr reaction )
{
	int 			i, j, k, l;
   	int				nOfSpecies = GetNOfSpecies();
	int				nOfReactions = reaction->GetNOfReactions();
	VectorPtr		*nuReaction = reaction->GetNu();
	IntVectorPtr	*speciesNumber = reaction->GetSpeciesNumber();
	SpeciesPtr		species = input->GetSpecies(); 
		
	//	fill fUsedReactions, fNu
	for ( i = 0; i < nOfSpecies; ++i ) { // i loops species
		l = 0;
		for ( j = 0; j < nOfReactions; ++j ) { // j loops reactions
			for ( k = 0; k < speciesNumber[j]->len; ++k ) { // k loops reactions speciesvector
			   	if ( i == speciesNumber[j]->vec[k] ) { // then species i appears in reaction j
					fUsedReactions[i]->vec[l] = j;
					fNu[i]->vec[l] = nuReaction[j]->vec[k];
					++l;
				}
			}
		}
	}
	
	Double	*coeffLow;
	Double	*coeffHigh;
	Double	*coeffCold;
	Double	*coeffHot;
	for ( i = 0; i < nOfSpecies; ++i ) {
		// names
		fNames[i] = new char[strlen( species[i].name )+1];
		if ( !fNames[i] ) FatalError( "memory allocation of TSpecies failed" );
		strcpy( fNames[i], species[i].name );
		// constants
		
		coeffLow = fCoeffLow->mat[i];
		coeffHigh = fCoeffHigh->mat[i];
		coeffCold = species[i].coeffCold;
		coeffHot = species[i].coeffHot;
		for ( j = 0; j < NOFCOEFFS; ++ j ) {
	   		coeffLow[j] = coeffCold[j];
	   		coeffHigh[j] = coeffHot[j];
		}
		
		coeffLow = fHCoeffLow->mat[i];
		coeffHigh = fHCoeffHigh->mat[i];
	   	coeffLow[0] = coeffCold[0];
	   	coeffHigh[0] = coeffHot[0];
	   	coeffLow[1] = coeffCold[1] * 0.5;
	   	coeffHigh[1] = coeffHot[1] * 0.5;
	   	coeffLow[2] = coeffCold[2] / 3.0;
	   	coeffHigh[2] = coeffHot[2] / 3.0;
	   	coeffLow[3] = coeffCold[3] * 0.25;
	   	coeffHigh[3] = coeffHot[3] * 0.25;
	   	coeffLow[4] = coeffCold[4] * 0.2;
	   	coeffHigh[4] = coeffHot[4] * 0.2;
	   	coeffLow[5] = coeffCold[5];
	   	coeffHigh[5] = coeffHot[5];
		
#ifndef PRODRATEFILEF77
	   	fMolarMass->vec[i] = species[i].molarMass;
#endif
	   	fKOverEps->vec[i] = species[i].k_over_eps;
	   	fSigma->vec[i] = species[i].sigma;
	   	fMuCoeff->vec[i] = species[i].muCoeff;
		for ( j = 0; j < nOfSpecies; ++ j ) {
	   		fDCoeff->mat[i][j] = species[i].dCoeff->vec[j];
	   		fOmegaCoeff->mat[i][j] = species[i].omegaCoeff->vec[j];
		}
		if ( species[i].isSteadyState ) {
			fIsSteadyState[i] = TRUE;
		}
		else {
			fIsSteadyState[i] = FALSE;
		}
	}


#ifdef PRODRATEFILEF77
	GETMOLARMASS( fMolarMass->vec );
#endif

	Double	*M = fMolarMass->vec;
	int nSpeciesIn = GetNSpeciesInSystem();

	for ( i = 0; i < nSpeciesIn; ++i ) {
		for ( j = 0; j < nSpeciesIn; ++j ) {
			sqrtMjOverMi[j][i] = sqrt( M[j] / M[i] );
			sqrtMiOverMjPlusOne[i][j] = 0.3535534 / sqrt( M[i] / M[j] + 1.0 );
		}
	}

	fConstLewisNumber = input->fConstantLewisNumber;
}

void TSpecies::CalcUsedReactions( TReactionPtr reaction )
{
	int	i, j, k;
	
	for ( i = 0; i < fNOfUsedReactions->len; ++i ) {
		for ( j = 0; j < reaction->GetNOfReactions(); ++j ) {
			for ( k = 0; k < reaction->GetSpeciesNumber()[j]->len; ++k ) {
				if ( i == reaction->GetSpeciesNumber()[j]->vec[k] ) {
					++fNOfUsedReactions->vec[i];
					continue;
				}
			}
		}
	}
}

/*
	omega_mu() returns the omega_mu for a given dimensionless temperature t/(eps/k).
*/
Double TSpecies::omega_mu( Double t )
{
	static Double m1 = 3.3530622607;
	static Double m2 = 2.53272006;
	static Double m3 = 2.9024238575;
	static Double m4 = 0.11186138893;
	static Double m5 = 0.8662326188;			/* = -0.1337673812 + 1.0 */
	static Double m6 = 1.3913958626;
	static Double m7 = 3.158490576;
	static Double m8 = 0.18973411754;
	static Double m9 = 0.00018682962894;
	Double num, den;

	num = m1 + t*(m2 + t*(m3 + t*m4));
	den = m5 + t*(m6 + t*(m7 + t*(m8 + t*m9)));
	return num / den;
}

/*
	omega_D() returns the Stossintegral for a given dimensionless 
	Temperature t/(eps/k)
*/
Double TSpecies::omega_D( Double t )
{					
	static Double m1 = 6.8728271691;
	static Double m2 = 9.4122316321;
	static Double m3 = 7.7442359037;
	static Double m4 = 0.23424661229;
	static Double m5 = 1.45337701568;			/* = 1.0 + 0.45337701568 */
	static Double m6 = 5.2269794238;
	static Double m7 = 9.7108519575;
	static Double m8 = 0.46539437353;
	static Double m9 = 0.00041908394781;

	return (m1 + t * (m2 + t * (m3 + t * m4))) / (m5 + t * (m6 + t * (m7 + t * (m8 + t * m9))));
}	 

Double TSpecies::GetPosNegProductionRate( int which, Double *reactionRate, Flag pos )
{
	int		j;
	int		nSpeciesInSystem = GetNSpeciesInSystem();

	Double	productionRate = 0.0;
	Double	*nu;
	Double	*W = fMolarMass->vec;
	int		*usedReactions;
	int		nOfUsedReactions;

	nu = fNu[which]->vec;
	usedReactions = fUsedReactions[which]->vec;
	nOfUsedReactions = fNOfUsedReactions->vec[which];
	for ( j = 0; j < nOfUsedReactions; ++j ) {
		if ( pos == TRUE && nu[j] < 0.0 ) {
			productionRate += nu[j] * reactionRate[usedReactions[j]];
		}
		else if ( pos == FALSE && nu[j] > 0.0 ) {
			productionRate += nu[j] * reactionRate[usedReactions[j]];
		}
	}
	productionRate *= - W[which];

	return productionRate;
}



void TSpecies::ComputeProductionRates( Double *productionRate, Double *reactionRate )
{
	int		i, j;
	int		nSpeciesInSystem = GetNSpeciesInSystem();

	Double	*nu;
	Double	*W = fMolarMass->vec;
	int		*usedReactions;
	int		nOfUsedReactions;

	for ( i = 0; i < nSpeciesInSystem; ++i ) {
//		ComputeProductionRate( i, productionRate[i], reactionRate );
		nu = fNu[i]->vec;
		usedReactions = fUsedReactions[i]->vec;
		nOfUsedReactions = fNOfUsedReactions->vec[i];
		productionRate[i] = 0.0;
		for ( j = 0; j < nOfUsedReactions; ++j ) {
			productionRate[i] += nu[j] * reactionRate[usedReactions[j]];
		}
		productionRate[i] *= - W[i];
	}
}


#ifdef PRODRATEFILE
void ComputeProductionRates( Double *cdot, Double *w, Double *k
								, Double *c, Double *M, Double temp, Double pressure );

void TSpecies::ComputeTheProductionRates( Double *productionRate, Double *reactionRate
						, Double temp, Double pressure, Double *c
						, Double *k, Double *M )
{
	::ComputeProductionRates( productionRate, reactionRate, k
			, c, M, temp, pressure );
}
#endif

#ifdef PRODRATEFILEF77
void TSpecies::ComputeTheProductionRates( Double *productionRate, Double *reactionRate
						, Double temp, Double pressure, Double *c
						, Double *k, Double *M )
{
	::PRODRATES( productionRate, reactionRate, k
			, c, M, &temp, &pressure );
}
#endif

void TSpecies::ComputeNoProductionRates( Double *productionRate )
{
	int		nSpeciesInSystem = GetNSpeciesInSystem();

	for ( int i = 0; i < nSpeciesInSystem; ++i ) {
		productionRate[i] = 0.0;
	}
}

void TSpecies::ComputeProductionRate( int i, Double &productionRate, Double *reactionRate )
{
	productionRate = 0.0;
	Double	*nu = fNu[i]->vec;
	int		*usedReactions = fUsedReactions[i]->vec;
	int		nOfUsedReactions = fNOfUsedReactions->vec[i];
	
	
	for ( int j = 0; j < nOfUsedReactions; ++j ) {
		productionRate += nu[j] * reactionRate[usedReactions[j]];
	}
	productionRate *= - fMolarMass->vec[i];
}

void TSpecies::CompHeatRelease( Double *heatRel, Double *prodRate, Double *enthalpy )
{	
	*heatRel = 0.0;

	for ( int i = 0; i < GetNSpeciesInSystem(); ++i ) {
		*heatRel += prodRate[i] * enthalpy[i];
	}
}

void TSpecies::CompProdCons( Double *source, Double *sink, Double *reacRate )
{
	Double	productionRate;
	int		*nUsedReacs = fNOfUsedReactions->vec;
	int		*usedReacs;
	Double	*nu;
//	Double	*molarMass = fMolarMass->vec;
	
	for ( int i = 0; i < GetNSpeciesInSystem(); ++i ) {
		usedReacs = fUsedReactions[i]->vec;
		nu = fNu[i]->vec;
		for ( int j = 0; j < nUsedReacs[i]; ++j ) {
			productionRate = nu[j] * reacRate[usedReacs[j]];
			if ( productionRate > 0.0 ) {
				sink[i] -= productionRate;
			}
			else {
				source[i] -= productionRate;
			}
		}
//		sink[i] *= molarMass[i];
//		source[i] *= molarMass[i];
	}
}

int TSpecies::CompDetailedProdRate( int i, Double *prodRate, Double *reacRate, TReactionPtr reaction )
{
	int				k;
	int				entry = 0;
	int				usedReac;
	int				*nUsedReacs = fNOfUsedReactions->vec;
	int				*forwReacs = reaction->GetForwardReacs()->vec;
	int				*backReacs = reaction->GetBackwardReacs()->vec;
	Double			stoiCoeffBack;
	char			**labels = reaction->GetLabels();
	
	for ( int j = 0; j < nUsedReacs[i]; ++j ) {
		usedReac = fUsedReactions[i]->vec[j];
		if ( forwReacs[usedReac] < 0 ) { // means not a backward reaction
			prodRate[entry] = fNu[i]->vec[j] * reacRate[usedReac];
			if ( backReacs[usedReac] >= 0 ) { // means there is a corresponding backward reaction
				for ( k = 0; k < nUsedReacs[i]; ++k ) {
					if ( strcmp( labels[backReacs[usedReac]], labels[fUsedReactions[i]->vec[k]] ) == 0 ) {
						stoiCoeffBack = fNu[i]->vec[k];
						break;
					}
				}
				prodRate[entry] += stoiCoeffBack * reacRate[backReacs[usedReac]];
			}
			prodRate[entry] *= - 1.0;//fMolarMass->vec[i];
			++entry;
		}
	}
	
	return entry;
}

Double TSpecies::ComputeDCpiDT( int index, Double temp )
{
	Double	*coeff = ( temp > 1000.0 ) ? fCoeffHigh->mat[index] : fCoeffLow->mat[index];
	
	return ( RGAS / fMolarMass->vec[index] * ( coeff[1] + temp * ( 2.0 * coeff[2] 
											+ temp * ( 3.0 * coeff[3] + temp * 4.0 * coeff[4] ) ) ) );
}

int TSpecies::FindSpecies( const char *species )
{
	int nOfSpecies = GetNOfSpecies();

	for ( int i = 0; i < nOfSpecies; ++ i ) {
		if ( strcmp( species, fNames[i] ) == 0 ) {
			return i;
		}
	}
	return -1;
}

void TSpecies::WriteRoggsSymbolsData( void )
{
// inerte Spezies muessen am Ende der Liste stehen

	FILE	*fp = fopen( "roggssymbols.data", "w" );
	if ( !fp ) {
		cerr << "#warning: unable to open file 'roggssymbols.data'" << NEWL;
		return;
	}	
	fprintf( fp, "LIST OF SYM / CHNO..........................................\n" );
	for ( int i = 0; i < GetNOfSpecies(); ++i ) {
		fprintf( fp, "%-5s%-8s 0000\n", "Y", fNames[i] );
	}
	fprintf( fp, "-END OF SYM   0000\n" );
	fprintf( fp, "BOUNDS ON TEMPERATURE.............\n" );
	fprintf( fp, "0298.00    2500.00\n" );
	fprintf( fp, "-END OF BOUNDS\n" );
	fprintf( fp, "CHEMTP: TYPE OF CHEMISTRY.........\n" );
	fprintf( fp, "DETAILED CHEMISTRY\n" );
	fprintf( fp, "-END OF CHEMTP\n" );
	fprintf( fp, "NSTMEC: PARAMETER.................\n" );
	fprintf( fp, " -1\n" );
	fprintf( fp, "-END OF NSTMEC\n" );

	fclose( fp );
}

#ifdef OPTIMIZEDELTAI
void T1DSpecies::ComputeDeltaIOpt( TFlameNodePtr flameNode, Double *Y, Double **GijOverWj, Flag newTemp )
{
	int 	i;
	int		nSpeciesIn = GetNSpeciesInSystem();
	Double	*deltaI = flameNode->deltaI;
	Double	*M = fMolarMass->vec;
	
#define NEWOPTDELTAI
#ifdef NEWOPTDELTAI
	Double			*savedDeltaiY = flameNode->savedDeltaiY;
	int				changed = -1;	// -1 means nothing changed, -2 means temperature
									// or more than two Y_i changed, >=0 is changed species, if only
									// one or two changed
	int				changed2 = -1;	// -1 means nothing changed, -2 means temperature

	if ( newTemp ) {
		changed = -2;
	}
	else {
		for ( i = 0; i < nSpeciesIn; ++i ) {
			if ( Y[i] != savedDeltaiY[i] ) {
				if ( changed2 < 0 ) {
					if ( changed == -1 ) {
						changed = i;
					}
					else {
						changed2 = i;
					}
				}
				else {
					changed = -2;
					break;
				}
			}
		}
	}

	if ( changed >= 0 /*&& changed2 == -1*/ ) {
		for ( i = 0; i < nSpeciesIn; ++i ) {
			deltaI[i] = deltaI[i] + M[i] * GijOverWj[i][changed] * ( Y[changed] - savedDeltaiY[changed] );
		}
		if ( changed2 >= 0 ) {
			for ( i = 0; i < nSpeciesIn; ++i ) {
				deltaI[i] = deltaI[i] + M[i] * GijOverWj[i][changed2] * ( Y[changed2] - savedDeltaiY[changed2] );
			}
		}
		savedDeltaiY[changed] = Y[changed];
		if ( changed2 >= 0 ) {
			savedDeltaiY[changed2] = Y[changed2];
		}
	}
	else if ( changed == -2 ) {
		for ( i = 0; i < nSpeciesIn; ++i ) {
			savedDeltaiY[i] = Y[i];
			deltaI[i] = CompDeltaIOpt( i, nSpeciesIn, GijOverWj, M, Y );
		}
	}
	else {
/*		fprintf( stderr, "changed = %d\tchanged2 = %d\n", changed, changed2 );*/
/*		Double	check = 0.0;*/
/*		for ( i = 0; i < nSpeciesIn; ++i ) {*/
/*			check = CompDeltaIOpt( i, nSpeciesIn, GijOverWj, M, Y );*/
/*		//	fprintf( stderr, "deltaI[%d] = %g\tcheck = %g\n", i, deltaI[i], check );*/
/*		}*/
	}
#else
	for ( i = 0; i < nSpeciesIn; ++i ) {
		deltaI[i] = CompDeltaIOpt( i, nSpeciesIn, GijOverWj, M, Y );
	}
#endif
}

Double TSpecies::CompDeltaIOpt( int i, int nSpeciesInSystem, Double **GijOverWj, Double *M, Double *Y )
{
	Double	Delta_i = 0;

/*	Double	*loci = GijOverWj[i];*/
/*	for ( int j = 0; j < nSpeciesInSystem; ++j ){*/
/*		Delta_i += loci[j] * Y[j];*/
/*	}*/
/*	*/
/*	return Delta_i * M[i];*/

	Delta_i = DotProd( nSpeciesInSystem, GijOverWj[i], 1, Y, 1 );

	return Delta_i * M[i];
}
#endif

void TSpecies::ComputeDeltaI( Double *deltaI, Double *Y, Double *viscosity )
{
	int 	i;
	int		nSpeciesIn = GetNSpeciesInSystem();
	Double	*M = fMolarMass->vec;
	
	for ( i = 0; i < nSpeciesIn; ++i ) {
		deltaI[i] = CompDeltaI( i, nSpeciesIn, M, viscosity, Y );
	}
}

Double TSpecies::CompDeltaI( int i, int nSpeciesInSystem, Double *M, Double *mu, Double *Y )
{
	Double	G_ij, Delta_i = 0;

	for ( int j = 0; j < nSpeciesInSystem; ++j ){
		G_ij = sqrt( sqrtMjOverMi[j][i] * mu[i] / mu[j] ) + 1.0;
//	(old form before bug was fixed)	G_ij = sqrt( sqrtMjOverMi[j][i] * mu[j] / mu[i] ) + 1.0;
		G_ij *= G_ij *  sqrtMiOverMjPlusOne[i][j];
		Delta_i += G_ij / M[j] * Y[j];
	}
	
	return Delta_i * M[i];
}

void TSpecies::ReadLewisNumbers( const char *file, VectorPtr Le )
{
	char				inSpecies[32];
	Double				inLe;	
	char				buffer[80];
	int					i, empty = 0;
	FILE 				*fp = NULL;
	Double				*lewis = Le->vec;
	
	for ( i = 0; i < Le->len; ++i ) {
		lewis[i] = 1.0;
	}
	
	if ( !(fp = fopen( file, "r" )) ) {
		fprintf( stderr, "Couldn't open \"%s\"---all Lewis numbers will be set to 1.0\n", file );
		return;
	}
		
	while ( fgets( buffer, 79, fp ) != NULL ) {
		if ( EmptyLine( buffer ) ) {
			++empty;								/*	Skip empty lines.				*/
		}
		else if ( buffer[0] == '#' ) {				/*	Skip line comments.				*/
		}
		else if ( sscanf( buffer, "%s%lf", inSpecies, &inLe ) == 2 ) {
			if ( ( i = FindSpecies( inSpecies ) ) < 0 ) {}
			else
				lewis[i] = inLe;
			}
		}

}

Flag TSpecies::EmptyLine( char *s )
{
	while ( *s ) {
		if ( !isspace(*s) ) return FALSE;
		++s;
	}
	
	return TRUE;
}

void T0DSpecies::InitT0DSpecies( TInputDataPtr input )
{
   	int			nOfSpecies = input->GetCounter()->species;
	int			nSpeciesInSystem = nOfSpecies - input->GetCounter()->steadyStates;
   
   	fHeatCapacity = NewVector( nSpeciesInSystem );
   	fEnthalpy = NewVector( nSpeciesInSystem );
   	fViscosity = NewVector( nSpeciesInSystem );
   	fConductivity = NewVector( nSpeciesInSystem );
   	fProductionRate = NewVector( nSpeciesInSystem );
	fDeltaI = NewVector( nOfSpecies );
}

T0DSpecies::~T0DSpecies( void )
{
	DisposeVector( fDeltaI );
	DisposeVector( fProductionRate );
	DisposeVector( fConductivity );
	DisposeVector( fViscosity );
	DisposeVector( fEnthalpy );
	DisposeVector( fHeatCapacity );
}


Flag T0DSpecies::ComputeSpeciesProperties( Double temp )
{
	int i;
	int	nSpeciesInSystem = GetNSpeciesInSystem();
	
	if ( fabs( ( fCurrTemp - temp ) / temp ) < 1.0e-5 ) {
//	if ( fCurrTemp == temp ) {\\}
		return FALSE;
	}
	else {
		fCurrTemp = temp;
//		cerr << "new thermo props" << NEWL;
	}
	
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		ComputeSpeciesProperties( i, temp );
	}

#ifdef PRODRATEFILE
	ComputeThermoData( fEnthalpy->vec, fHeatCapacity->vec, temp );
#elif defined PRODRATEFILEF77
	COMPTHERMODATA( fEnthalpy->vec, fHeatCapacity->vec, &temp );
#endif
	return TRUE;
}

void T0DSpecies::ComputeSpeciesProperties( int number, Double temp )
{
//	computes heatCapacity, enthalpy for the
//	species with index 'number' for the temperature temp
	
   	Double		*coeff, *hCoeff;			/* pointer to nasa coefficients */
	Double		r_over_m = RGAS / fMolarMass->vec[number];		/* in J / (kg K) */
	Double		*heatCapacity = fHeatCapacity->vec;
	Double		*viscosity = fViscosity->vec;
	
	if ( temp > 1000.0 ) {
		coeff = fCoeffHigh->mat[number];
		hCoeff = fHCoeffHigh->mat[number];
	}
	else {
		coeff = fCoeffLow->mat[number];
		hCoeff = fHCoeffLow->mat[number];
	}
	if ( temp < fTThermoLimit ) {
		heatCapacity[number] = r_over_m * ( coeff[0]+fTThermoLimit*(coeff[1]+fTThermoLimit*(coeff[2]+fTThermoLimit*(coeff[3]+fTThermoLimit*coeff[4]))));
		
		fEnthalpy->vec[number] = r_over_m * ( fTThermoLimit*(hCoeff[0]+fTThermoLimit*(hCoeff[1]+fTThermoLimit*(hCoeff[2]+
			fTThermoLimit*(hCoeff[3]+fTThermoLimit*hCoeff[4]))))+hCoeff[5]);
#ifdef CONSTLOWTCP
		fEnthalpy->vec[number] += (temp-fTThermoLimit)*heatCapacity[number];
#else
		Double B = r_over_m * ( coeff[1]+fTThermoLimit*(2.0*coeff[2]+fTThermoLimit*(3.0*coeff[3]+4.0*fTThermoLimit*coeff[4])));
		Double cpRef = heatCapacity[number];
		Double A = cpRef - B * fTThermoLimit;
		heatCapacity[number] = A + B * temp;
		fEnthalpy->vec[number] += 0.5 * (heatCapacity[number]+cpRef) * (temp - fTThermoLimit);
#endif
	}
	else {
		heatCapacity[number] = r_over_m * ( coeff[0]+temp*(coeff[1]+temp*(coeff[2]+temp*(coeff[3]+temp*coeff[4]))));
		
		fEnthalpy->vec[number] = r_over_m * ( temp*(hCoeff[0]+temp*(hCoeff[1]+temp*(hCoeff[2]+
			temp*(hCoeff[3]+temp*hCoeff[4]))))+hCoeff[5]);
	}					

	viscosity[number] = fMuCoeff->vec[number] * sqrt( temp ) / omega_mu( temp * fKOverEps->vec[number] );

	fConductivity->vec[number] = viscosity[number] * ( heatCapacity[number] + 1.2 * r_over_m );
}

void T0DSpecies::PrintSpecies( T0DFlamePtr flame )
{
	FILE	*fp = flame->GetOutfile( "species", TFlame::kText );
	
   	for ( int i = 0; i < GetNOfSpecies(); ++i ) {
		PrintSpecies( i, flame, fp );
	}

	fclose( fp );
}

void T0DSpecies::PrintSpecies( int number, T0DFlamePtr flame, FILE *fp )
{
	int i;	
	
	fprintf( fp, "%s has number %d and appears in %d reactions, which are\n", fNames[number], number, fNOfUsedReactions->vec[number] );
	for ( i = 0; i < fNOfUsedReactions->vec[number]; ++i ) {
		flame->GetReaction()->PrintReactionEquation( fUsedReactions[number]->vec[i], this, fp );
		fprintf( fp, "\t%g\n", flame->GetReaction()->
				GetReactionRate()->vec[fUsedReactions[number]->vec[i]] );
//		fprintf( fp, "%d\n", fUsedReactions[number]->vec[i] );
	}

	fprintf( fp, "HeatCapacity = %g\n", fHeatCapacity->vec[number] );
	fprintf( fp, "Enthalpy = %g\n", fEnthalpy->vec[number] );
	
	fprintf( fp, "fCoeffLow =");
	for ( i = 0; i < 7; ++i ) {
		fprintf( fp, "\t%g", fCoeffLow->mat[number][i] );
	}
	fprintf( fp, "\n");

	fprintf( fp, "fCoeffHigh =");
	for ( i = 0; i < 7; ++i ) {
		fprintf( fp, "\t%g", fCoeffHigh->mat[number][i] );
	}
	fprintf( fp, "\n");

	fprintf( fp, "molarMass = %g\n", fMolarMass->vec[number] );

	if ( !fIsSteadyState[number] ) {
		fprintf( fp, "productionRate = %g\n", fProductionRate->vec[number] );
	}
	fprintf( fp, "\n\n\n");
	
}

#ifndef ZEROD
void T1DSpecies::InitT1DSpecies( TInputDataPtr input )
{
   	int			nOfSpecies = input->GetCounter()->species;
	int			nSpeciesInSystem = nOfSpecies - input->GetCounter()->steadyStates;
   	int			maxGridPoints = input->fMaxGridPoints;
   
   	fViscosity = NewMatrix( nOfSpecies, maxGridPoints+2, kColumnPointers );
	fViscosity->mat = &fViscosity->mat[kNext];
   	fHeatCapacity = NewMatrix( nOfSpecies, maxGridPoints+2, kColumnPointers );
	fHeatCapacity->mat = &fHeatCapacity->mat[kNext];
	fConductivity = NewMatrix( nOfSpecies, maxGridPoints+2, kColumnPointers );
	fConductivity->mat = &fConductivity->mat[kNext];
	fEnthalpy = NewMatrix( nOfSpecies, maxGridPoints+2, kColumnPointers );
	fEnthalpy->mat = &fEnthalpy->mat[kNext];
	fDiffusivity = NewMatrix( nOfSpecies, maxGridPoints+2, kColumnPointers );
	fDiffusivity->mat = &fDiffusivity->mat[kNext];
	fThermoDiffusivity = NewMatrix( nOfSpecies, maxGridPoints+2, kColumnPointers );
	fThermoDiffusivity->mat = &fThermoDiffusivity->mat[kNext];
    fBinDij = NewTensor( maxGridPoints+2, nSpeciesInSystem, nSpeciesInSystem, kColumnPointers );
	fBinDij->tensor = &fBinDij->tensor[kNext];
    fGijOverWj = NewTensor( maxGridPoints+2, nSpeciesInSystem, nSpeciesInSystem, kColumnPointers );
	fGijOverWj->tensor = &fGijOverWj->tensor[kNext];
    fOmegaDOverDCoeff = NewTensor( maxGridPoints+2, nSpeciesInSystem, nSpeciesInSystem, kColumnPointers );
	fOmegaDOverDCoeff->tensor = &fOmegaDOverDCoeff->tensor[kNext];

	if ( input->fThermoDiffusion ) {
		fDThermConst = NewTensor( maxGridPoints+2, nSpeciesInSystem, nSpeciesInSystem, kColumnPointers );
		fDThermConst->tensor = &fDThermConst->tensor[kNext];
	}
	else {
		fDThermConst = NULL;
	}
	
	fProductionRate = NewMatrix( nSpeciesInSystem, maxGridPoints, kColumnPointers );

	fSavedDeltaiY = NewMatrix( nOfSpecies, maxGridPoints+2, kColumnPointers );
	fSavedDeltaiY->mat = &fSavedDeltaiY->mat[kNext];
	fSavedY = NewMatrix( nOfSpecies, maxGridPoints+2, kColumnPointers );
	fSavedY->mat = &fSavedY->mat[kNext];
	fSumDiff = NewMatrix( nOfSpecies, maxGridPoints+2, kColumnPointers );
	fSumDiff->mat = &fSumDiff->mat[kNext];
	fDeltaI = NewMatrix( nOfSpecies, maxGridPoints+2, kColumnPointers );
	fDeltaI->mat = &fDeltaI->mat[kNext];
	fTempProp = NewVector( maxGridPoints+2 );
	fTempProp->vec = &fTempProp->vec[kNext];
	fPressureProp = NewVector( maxGridPoints+2 );
	fPressureProp->vec = &fPressureProp->vec[kNext];
}

T1DSpecies::~T1DSpecies( void )
{
	fTempProp->vec = &fTempProp->vec[kPrev];
	DisposeVector( fTempProp );
	fPressureProp->vec = &fPressureProp->vec[kPrev];
	DisposeVector( fPressureProp );
	fDeltaI->mat = &fDeltaI->mat[kPrev];
	DisposeMatrix( fDeltaI );
	fSumDiff->mat = &fSumDiff->mat[kPrev];
	DisposeMatrix( fSumDiff );
	fSavedDeltaiY->mat = &fSavedDeltaiY->mat[kPrev];
	DisposeMatrix( fSavedDeltaiY );
	fSavedY->mat = &fSavedY->mat[kPrev];
	DisposeMatrix( fSavedY );
	
	DisposeMatrix( fProductionRate );

	if ( fDThermConst ) {
		fDThermConst->tensor = &fDThermConst->tensor[kPrev];
		DisposeTensor( fDThermConst );
	}
	fOmegaDOverDCoeff->tensor = &fOmegaDOverDCoeff->tensor[kPrev];
	DisposeTensor( fOmegaDOverDCoeff );
	fGijOverWj->tensor = &fGijOverWj->tensor[kPrev];
	DisposeTensor( fGijOverWj );
	fBinDij->tensor = &fBinDij->tensor[kPrev];
	DisposeTensor( fBinDij );
	fThermoDiffusivity->mat = &fThermoDiffusivity->mat[kPrev];
	DisposeMatrix( fThermoDiffusivity );
	fDiffusivity->mat = &fDiffusivity->mat[kPrev];
	DisposeMatrix( fDiffusivity );
	fEnthalpy->mat = &fEnthalpy->mat[kPrev];
	DisposeMatrix( fEnthalpy );
	fConductivity->mat = &fConductivity->mat[kPrev];
	DisposeMatrix( fConductivity );
	fHeatCapacity->mat = &fHeatCapacity->mat[kPrev];
	DisposeMatrix( fHeatCapacity );
	fViscosity->mat = &fViscosity->mat[kPrev];
	DisposeMatrix( fViscosity );
}

void T1DSpecies::PrintSpecies( int k )
{
	FILE	*fp;
	
	if ( !( fp = fopen( "species.tout", "w" ) ) ) { 
		cerr << "#warning: unable to open file 'species.tout'" << NEWL;
		return;
	}
	fprintf( fp, "GridPoint no. %d\n", k );
   	for ( int i = 0; i < GetNOfSpecies(); ++i ) {
		PrintSpecies( i, k, fp );
	}
	fclose( fp );
}

void T1DSpecies::PrintProductionRate( TNewtonPtr bt )
{
// print productionrate in [mole/cm^3] 	

	int 		i, k;
	char 		fileName[32];
	FILE		*fp = NULL;
	Double		**mat = fProductionRate->mat;
	Double		*molarMass = fMolarMass->vec;
	char		**names = GetNames();
	TGridPtr	currentGrid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = currentGrid->GetX()->vec;
	int			gridPoints = currentGrid->GetNGridPoints();
	static int	counter = 0;

	if ( gridPoints > fProductionRate->cols ) {
		fprintf( stderr, "#warning: productionrate not printed, because it doesn't correspond to current grid\n" );
		return;
	}
	
	sprintf( fileName, "ProdRate%d.dout", ++counter );
	if ( !( fp = fopen( fileName, "w" ) ) ) { 
		cerr << "#warning: unable to open file " << fileName << NEWL;
		return;
	}
	if ( counter >= 10 ) counter = 0;

	fprintf( fp, "*\n" );
	fprintf( fp, "%-9s\t%-12s", "gridPoint", "eta" );
	for ( i = 0; i < GetNSpeciesInSystem(); ++i ) {
    	fprintf( fp, "\t%-12s", names[i] );
	}

	fprintf( fp, "\n" );
	
	for ( k = 0; k < gridPoints; ++k ) {
		fprintf( fp, "%9d\t%12.4E", k, x[k] );
		for ( i = 0; i < GetNSpeciesInSystem(); ++i ) {
			fprintf( fp, "\t%12.4E", mat[k][i] / molarMass[i] * 1.0e-3 );
		}
		fprintf( fp, "\n" );
	}
	
	fclose( fp );
}

void T1DSpecies::PrintSpecies( int number, int gridPoint, FILE *fp )
{
	int i;	
	static Flag init = FALSE;
	
	if ( !fp ) {
		if ( !init ) {
			fp = fopen( "species.tout", "w" );
			fprintf( fp, "GridPoint no. %d\n", gridPoint );
		}
		else {
			fp = fopen( "species.tout", "a" );
		}
	}
	if ( !fp ) {
		cerr << "#warning: unable to open file 'species.tout'" << NEWL;
		return;
	}
	init = TRUE;
	if ( !fIsSteadyState[number] ) {
		fprintf( fp, "%s has number %d and appears in %d reactions, which are\n", fNames[number], number, fNOfUsedReactions->vec[number] );
		for ( i = 0; i < fNOfUsedReactions->vec[number]; ++i ) {
			fprintf( fp, "%d\n", fUsedReactions[number]->vec[i] );
		}
	}
	else {
		fprintf( fp, "%s has number %d and is steady state\n", fNames[number], number );
	}
	fprintf( fp, "Viscosity = %g\n", fViscosity->mat[gridPoint][number] );
	fprintf( fp, "HeatCapacity = %g\n", fHeatCapacity->mat[gridPoint][number] );
	fprintf( fp, "Conductivity = %g\n", fConductivity->mat[gridPoint][number] );
	fprintf( fp, "Enthalpy = %g\n", fEnthalpy->mat[gridPoint][number] );
	fprintf( fp, "Diffusivity = %g\n", fDiffusivity->mat[gridPoint][number] );
	if ( IsConstantLewisNumber() ) {
		fprintf( fp, "LewisNumber = %g\n", fLewisNumber->vec[number] );
	}
	
	fprintf( fp, "fCoeffLow =");
	for ( i = 0; i < 7; ++i ) {
		fprintf( fp, "\t%g", fCoeffLow->mat[number][i] );
	}
	fprintf( fp, "\n");

	fprintf( fp, "fCoeffHigh =");
	for ( i = 0; i < 7; ++i ) {
		fprintf( fp, "\t%g", fCoeffHigh->mat[number][i] );
	}
	fprintf( fp, "\n");

	fprintf( fp, "molarMass = %g\n", fMolarMass->vec[number] );
	fprintf( fp, "kOverEps = %g\n", fKOverEps->vec[number] );
	fprintf( fp, "sigma = %g\n", fSigma->vec[number] );
	fprintf( fp, "muCoeff = %g\n", fMuCoeff->vec[number] );

	fprintf( fp, "fDCoeff =");
	for ( i = 0; i < GetNOfSpecies(); ++i ) {
		fprintf( fp, "\t%g", fDCoeff->mat[number][i] );
	}
	fprintf( fp, "\n");

	fprintf( fp, "fOmegaCoeff =");
	for ( i = 0; i < GetNOfSpecies(); ++i ) {
		fprintf( fp, "\t%g", fOmegaCoeff->mat[number][i] );
	}
	fprintf( fp, "\n");
	if ( !fIsSteadyState[number] ) {
		fprintf( fp, "productionRate = %g\n", fProductionRate->mat[gridPoint][number] );
	}
	fprintf( fp, "\n\n\n");
	
}

Flag T1DSpecies::ComputeSpeciesProperties( TFlameNodePtr flameNode, Double temp, Double pressure )
{
	int i;
	int	nSpeciesInSystem = GetNSpeciesInSystem();
	
#ifdef OPTPROPS
	if ( flameNode->tempProp[kCurr] == temp ) {
		return FALSE;
	}
	else {
		flameNode->tempProp[kCurr] = temp;
	}
#else
	flameNode->tempProp[kCurr] = temp;
#endif
	
#define INLINESPECIESPROPS
#ifdef INLINESPECIESPROPS
	Double		**coMat = fCoeffHigh->mat, **hCoMat = fCoeffLow->mat;	/* pointer to nasa coefficients */
	Double		*coeff, *hCoeff;			/* pointer to nasa coefficients */
	Double		*heatCapacity = flameNode->heatCapacity;
	Double		*vis = flameNode->viscosity;
	Double		*ent = flameNode->enthalpy;
	Double		*conduc = flameNode->conductivity;
	Double		*kOverEps = fKOverEps->vec;
	Double		*W = fMolarMass->vec;
	Double		*muCoeff = fMuCoeff->vec;
	Double		sqTemp = sqrt( temp );;
	Double		R_gas = 8314.4;
	Double		r_over_m;

	if ( temp < 10.0 ) {
		temp = 10.0;
	}

	if ( temp > 1000.0 ) {
		coMat = fCoeffHigh->mat;
		hCoMat = fHCoeffHigh->mat;
	}
	else {
		coMat = fCoeffLow->mat;
		hCoMat = fHCoeffLow->mat;
	}

/*	Double	t2 = temp * temp;*/
/*	Double	t3 = t2 * temp;*/
/*	Double	t4 = t3 * temp;*/
/*	Double	t5 = t4 * temp;*/
	for( i = 0; i < nSpeciesInSystem; ++i ) {
		coeff = coMat[i];
		hCoeff = hCoMat[i];
		r_over_m = R_gas / W[i];
		if ( temp < fTThermoLimit ) {
			heatCapacity[i] = r_over_m * ( coeff[0]+fTThermoLimit*(coeff[1]+fTThermoLimit*(coeff[2]+fTThermoLimit*(coeff[3]+fTThermoLimit*coeff[4]))));
			
			ent[i] = r_over_m * ( fTThermoLimit*(hCoeff[0]+fTThermoLimit*(hCoeff[1]+fTThermoLimit*(hCoeff[2]+
					fTThermoLimit*(hCoeff[3]+fTThermoLimit*hCoeff[4]))))+hCoeff[5]);
#ifdef CONSTLOWTCP
			ent[i] += (temp-fTThermoLimit)*heatCapacity[i];
#else
			Double B = r_over_m * ( coeff[1]+fTThermoLimit*(2.0*coeff[2]+fTThermoLimit*(3.0*coeff[3]+4.0*fTThermoLimit*coeff[4])));
			Double cpRef = heatCapacity[i];
			Double A = cpRef - B * fTThermoLimit;
			heatCapacity[i] = A + B * temp;
			ent[i] += 0.5 * (heatCapacity[i]+cpRef) * (temp - fTThermoLimit);
#endif
		}
		else {
			heatCapacity[i] = r_over_m * ( coeff[0]+temp*(coeff[1]+temp*(coeff[2]+temp*(coeff[3]+temp*coeff[4]))));
			
			ent[i] = r_over_m * ( temp*(hCoeff[0]+temp*(hCoeff[1]+temp*(hCoeff[2]+
					temp*(hCoeff[3]+temp*hCoeff[4]))))+hCoeff[5]);
		}
/*		heatCapacity[i] = (coeff[0] + coeff[1]*temp + coeff[2]*t2 + coeff[3]*t3 + coeff[4]*t4 ) * r_over_m;*/
/*		ent[i] = (hCoeff[5] + hCoeff[0]*temp + hCoeff[1]*t2 + hCoeff[2]*t3 + hCoeff[3]*t4 + hCoeff[4]*t4 ) * r_over_m;*/
								
		vis[i] = muCoeff[i] * sqTemp / omega_mu( temp * kOverEps[i] );
	
		conduc[i] = vis[i] * ( heatCapacity[i] + 1.2 * r_over_m );
	}
#else
	for( i = 0; i < nSpeciesInSystem; ++i ) {
		ComputeSpeciesProperties( flameNode, temp, pressure, i );
	}
#endif
#if defined (OPTIMIZEDELTAI) || defined (OPTOMEGAD) || defined (OPTDTHERM)
	int		j;
	int		nSpeciesIn = GetNSpeciesInSystem();
	Double	**GijOverWj = flameNode->GijOverWj;
	Double	*viscosity = flameNode->viscosity;
	Double	**omegaCoeff = fOmegaCoeff->mat;
	Double	**dCoeff = fDCoeff->mat;
	Double	*W_i = fMolarMass->vec;
	Double	pressureInvT3 = pressure / ( temp * sqrt( temp ) );
	flameNode->pressureProp[kCurr] = pressure;
#	ifdef OPTIMIZEDELTAI
	Double	*GijOverWj_i, *sqrtMiOverMjPlusOne_i;
	for( i = 0; i < nSpeciesInSystem; ++i ) {
		GijOverWj_i = GijOverWj[i];
		for ( j = 0; j <= i; ++j ) {
			GijOverWj_i[j] = sqrt( sqrtMjOverMi[j][i] * viscosity[i] / viscosity[j] );
//	(old form before bug was fixed)	GijOverWj_i[j] = sqrt( sqrtMjOverMi[j][i] * viscosity[j] / viscosity[i] );
		}
	}
	for( i = 0; i < nSpeciesInSystem; ++i ) {
		GijOverWj_i = GijOverWj[i];
		for ( j = i+1; j < nSpeciesInSystem; ++j ) {
			GijOverWj_i[j] = 1.0 / GijOverWj[j][i];
		}
	}
	for( i = 0; i < nSpeciesInSystem; ++i ) {
		GijOverWj_i = GijOverWj[i];
		sqrtMiOverMjPlusOne_i = sqrtMiOverMjPlusOne[i];
		for ( j = 0; j < nSpeciesInSystem; ++j ) {
			GijOverWj_i[j] += 1.0;
			GijOverWj_i[j] *= GijOverWj_i[j] * sqrtMiOverMjPlusOne_i[j]
								/ W_i[j];
/*			GijOverWj[i][j] = sqrt( sqrtMjOverMi[j][i] * viscosity[j] / viscosity[i] ) + 1.0;*/
/*			GijOverWj[i][j] *= GijOverWj[i][j] * sqrtMiOverMjPlusOne[i][j] */
/*								/ W_i[j];*/
		}
	}
#	endif
#	ifdef OPTOMEGAD
/*	GetFuncOne( flameNode->OmegaDOverDCoeff, nSpeciesInSystem, pressureInvT3, omegaCoeff,*/
/*				dCoeff );*/
	Double	**omegaDOverDCoeff = flameNode->OmegaDOverDCoeff;
	for( i = 0; i < nSpeciesInSystem; ++i ) {
		for ( j = 0; j < i; ++j ) {
			omegaDOverDCoeff[i][j] = pressureInvT3 * omega_D( temp * omegaCoeff[i][j] ) / dCoeff[i][j];
		}
		omegaDOverDCoeff[i][i] = 0.0;
	}
	for( i = 0; i < nSpeciesInSystem; ++i ) {
		for ( j = i+1; j < nSpeciesInSystem; ++j ) {
			omegaDOverDCoeff[i][j] = omegaDOverDCoeff[j][i];
		}
	}
#	endif
#endif
	return TRUE;
}

/*void TSpecies::GetFuncOne( Double **omegaDOverDCoeff, int nSpeciesInSystem, Double pressureInvT3, Double **omegaCoeff,*/
/*				Double **dCoeff )*/
/*{*/
/*	int i, j;*/
/*	for( i = 0; i < nSpeciesInSystem; ++i ) {*/
/*		for ( j = 0; j < i; ++j ) {*/
/*			omegaDOverDCoeff[i][j] = pressureInvT3 * omega_D( temp * omegaCoeff[i][j] ) / dCoeff[i][j];*/
/*		}*/
/*		omegaDOverDCoeff[i][i] = 0.0;*/
/*	}*/
/*	for( i = 0; i < nSpeciesInSystem; ++i ) {*/
/*		for ( j = i+1; j < nSpeciesInSystem; ++j ) {*/
/*			omegaDOverDCoeff[i][j] = omegaDOverDCoeff[j][i];*/
/*		}*/
/*	}*/
/*}*/

void T1DSpecies::ComputeSpeciesProperties( TFlameNodePtr flameNode, Double temp, Double /*pressure*/, int number )
{
//	computes heatCapacity, viscosity, enthalpy and thermal conductivity for the
//	species with index 'number' for the temperature temp
	
   	Double		*coeff, *hCoeff;			/* pointer to nasa coefficients */
	Double		R_gas = 8314.4;
	Double		r_over_m = R_gas / fMolarMass->vec[number];		/* in J / (kg K) */
	Double		*heatCapacity = flameNode->heatCapacity;
	Double		*viscosity = flameNode->viscosity;
	Double		*kOverEps = fKOverEps->vec;
	
	if ( temp < 10.0 ) {
		temp = 10.0;
	}
	if ( temp > 1000.0 ) {
		coeff = fCoeffHigh->mat[number];
		hCoeff = fHCoeffHigh->mat[number];
	}
	else {
		coeff = fCoeffLow->mat[number];
		hCoeff = fHCoeffLow->mat[number];
	}
	
	if ( temp < fTThermoLimit ) {
		heatCapacity[number] = r_over_m * ( coeff[0]+fTThermoLimit*(coeff[1]+fTThermoLimit*(coeff[2]+fTThermoLimit*(coeff[3]+fTThermoLimit*coeff[4]))));
	
		flameNode->enthalpy[number] = r_over_m * ( fTThermoLimit*(hCoeff[0]+fTThermoLimit*(hCoeff[1]+fTThermoLimit*(hCoeff[2]+
			fTThermoLimit*(hCoeff[3]+fTThermoLimit*hCoeff[4]))))+hCoeff[5]);
#ifdef CONSTLOWTCP
		flameNode->enthalpy[number] += (temp-fTThermoLimit)*heatCapacity[number];
#else
		Double B = r_over_m * ( coeff[1]+fTThermoLimit*(2.0*coeff[2]+fTThermoLimit*(3.0*coeff[3]+4.0*fTThermoLimit*coeff[4])));
		Double A = heatCapacity[number] - B * fTThermoLimit;
		Double cpRef = heatCapacity[number];
		heatCapacity[number] = A + B * temp;
		flameNode->enthalpy[number] += 0.5 * (heatCapacity[number]+cpRef) * (temp - fTThermoLimit);
#endif
	}
	else {
		heatCapacity[number] = r_over_m * ( coeff[0]+temp*(coeff[1]+temp*(coeff[2]+temp*(coeff[3]+temp*coeff[4]))));
	
		flameNode->enthalpy[number] = r_over_m * ( temp*(hCoeff[0]+temp*(hCoeff[1]+temp*(hCoeff[2]+
			temp*(hCoeff[3]+temp*hCoeff[4]))))+hCoeff[5]);
	}				
	viscosity[number] = fMuCoeff->vec[number] * sqrt( temp ) / omega_mu( temp * kOverEps[number] );

	flameNode->conductivity[number] = viscosity[number] * ( heatCapacity[number] + 1.2 * r_over_m );


}

void T1DSpecies::Compute_D( TFlameNodePtr flameNode )
{
	int		i;
	int		nSpeciesInSystem = GetNSpeciesInSystem();
	Double	*diffusivity = flameNode->diffusivity;
	Double	num = *flameNode->mixConductivity / ( *flameNode->mixDensity * *flameNode->mixHeatCapacity );
	Double	*lewis = fLewisNumber->vec;
	
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		diffusivity[i] = num / lewis[i];
	}
	
}

#ifdef FULLDIFFUSION
void T1DSpecies::Compute_D( TFlameNodePtr flameNode, Double temp, Double *Y, Double press, Flag newTemp )
{
	int		i, j, k;
	int		nSpeciesIn = GetNSpeciesInSystem();
	Double	**mat = NULL;
	Double	**inv = NULL;
	Double	*col = NULL;
	int		*index = NULL;
	Double	**dCoeff = fDCoeff->mat;
	Double	**omegaCoeff = fOmegaCoeff->mat;
	Double	**invDij = fInvDij->mat;
	Double	**binDij = flameNode->binDij[kCurr];
	Double	*mw = fMolarMass->vec;
	Double	pressureInvT3 = press / ( temp * sqrt( temp ) );
	const Double	zero = 0.0;
	
	mat = New2DArray( nSpeciesIn, nSpeciesIn );
	inv = New2DArray( nSpeciesIn, nSpeciesIn );
	index = New1DIntArray( nSpeciesIn );
	col = New1DArray( nSpeciesIn );
	
	//	compute the reciprocal binary diffusion coefficients
	for ( i = 0; i < nSpeciesIn; ++i ) {
		invDij[i][i] = zero;

		for (j = 0; j < i; ++j ) {
			invDij[i][j] = invDij[j][i] =  
#ifdef OPTOMEGAD
					flameNode->OmegaDOverDCoeff[i][j]; // pressureInvT3 is in OmegaDOverDCoeff
#else
					pressureInvT3 * omega_D( temp * omegaCoeff[i][j] ) / dCoeff[i][j];
#endif
		}
	}

	
	//	assemble the matrix
	for ( i = 0; i < nSpeciesIn; ++i ) {
		for ( j = 0; j < nSpeciesIn; ++j ) {
			if ( i == j ) continue;

			for ( k = 0; k < nSpeciesIn; ++k ) {
				if ( k == i ) continue;

				mat[i][j] += Y[k] / mw[k] * invDij[i][k];
			}
			mat[i][j] *= mw[j];
			mat[i][j] += Y[i] * invDij[i][j];		
		}
	}
	
	InverseMatrix( nSpeciesIn, mat, inv, index, col );

	for ( i = 0; i < nSpeciesIn; ++i ) {
		for ( j = 0; j < nSpeciesIn; ++j ) {
			binDij[i][j] = inv[i][j] - mw[i]/mw[j] * inv[i][i];
		}
	}
	
	Free1DArray( col );
	Free1DIntArray( index );
	Free2DArray( inv );
	Free2DArray( mat );
}
#else
void T1DSpecies::Compute_D( TFlameNodePtr flameNode, Double temp, Double *Y, Double pressure, Flag newTemp )
{
	int				j, i;
	int				nSpeciesInSystem = GetNSpeciesInSystem();
//	static int		cont1 = 0, cont2 = 0, contAll = 0, cont0 = 0;
	Double			*molarMass = fMolarMass->vec;
#ifndef OPTOMEGAD
	Double			**dCoeff = fDCoeff->mat;
	Double			**omegaCoeff = fOmegaCoeff->mat;
#endif
	Double			*diffusivity = flameNode->diffusivity;
	Double 			mixtureMolarMass = *flameNode->mixMolarMass;
	Double			**invDij = fInvDij->mat;
	Double			sumYi = 0.0;
	Double			pressureInvT3 = pressure / ( temp * sqrt( temp ) );
	Double			*savedY = flameNode->savedY;
	Double			*sumDiff = flameNode->sumDiff;
	int				changed = -1;	// -1 means nothing changed, -2 means temperature or pressure
									// or more than one Y_i changed, >=0 is changed species, if only
									// one changed
	int				changed2 = -1;	// -1 means nothing changed, -2 means temperature or pressure

	if ( newTemp ) {
		changed = -2;
	}
	else {
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			if ( Y[i] != savedY[i] ) {
				if ( changed2 < 0 ) {
					if ( changed == -1 ) {
						changed = i;
					}
					else {
						changed2 = i;
					}
				}
				else {
					changed = -2;
					break;
				}
			}
		}
	}

		
#ifdef OPTOMEGAD
	if ( flameNode->pressureProp[kCurr] != pressure ) {
		changed = -2;
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			invDij[i][i] = 0.0;
			for ( j = 0; j < i; ++j ) {
				invDij[i][j] = invDij[j][i] = flameNode->OmegaDOverDCoeff[i][j] 
						* pressure / flameNode->pressureProp[kCurr]; // pressureInvT3 is in OmegaDOverDCoeff
	
			}
			flameNode->pressureProp[kCurr] = pressure;
		}
	}
	else {
//			invDij[i][j] = invDij[j][i]  = flameNode->OmegaDOverDCoeff[i][j];
		invDij = flameNode->OmegaDOverDCoeff;
	}
#else
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		for ( j = 0; j < i; ++j ) {
			invDij[i][j] = invDij[j][i]  =  
//					* flameNode->OmegaDOverDCoeff[i][j]; // pressureInvT3 is in OmegaDOverDCoeff
					pressureInvT3 * omega_D( temp * omegaCoeff[i][j] ) / dCoeff[i][j];

		}
	}
#endif
		

	if ( changed >= 0 ) {
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			sumYi += Y[i];
		}
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
	
			sumDiff[i] = sumDiff[i] + invDij[i][changed] / molarMass[changed]
							* ( Y[changed] - savedY[changed] );
			if ( changed2 >= 0 ) {
				sumDiff[i] = sumDiff[i] + invDij[i][changed2] / molarMass[changed2]
							* ( Y[changed2] - savedY[changed2] );
			}


			if ( sumDiff[i] > 1.0e-15 ) {
				diffusivity[i] = ( sumYi - Y[i] ) / ( mixtureMolarMass * sumDiff[i] );
			}
			else {
				for ( j = 1, sumDiff[i] = 0.0; j < nSpeciesInSystem; ++j ) {
					sumDiff[i] += invDij[i][j] / molarMass[j];
				}
				diffusivity[i] = ( ( Double )nSpeciesInSystem - 1.0 ) / ( mixtureMolarMass * sumDiff[i] );
			}
//			if ( isnan( diffusivity[i] ) ) {
//				fprintf( stderr, "diff[%s] is %g\n", fNames[i], diffusivity[i] );
//				diffusivity[i] = 0.0;
//			}
		}
		savedY[changed] = Y[changed];
		if ( changed2 >= 0 ) {
			savedY[changed2] = Y[changed2];
//			++cont2;
		}
//		else {
//			++cont1;
//		}
	}
	else if ( changed == -2 ) {
//		++contAll;
		Double	*YiOverWi = New1DArray( nSpeciesInSystem );
		sumYi = 0.0;
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			YiOverWi[i] = Y[i] / molarMass[i];
			sumYi += Y[i];
		}
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
	
			savedY[i] = Y[i];
	
			sumDiff[i] = YiOverWi[0] * invDij[i][0];
			for ( j = 1; j < nSpeciesInSystem; ++j ) {
				sumDiff[i] += YiOverWi[j] * invDij[i][j];
			}
			
			if ( sumDiff[i] > 1.0e-15 ) {
				diffusivity[i] = ( sumYi - Y[i] ) / ( mixtureMolarMass * sumDiff[i] );
			}
			else {
				for ( j = 1, sumDiff[i] = invDij[i][0] / molarMass[0]; j < nSpeciesInSystem; ++j ) {
					sumDiff[i] += invDij[i][j] / molarMass[j];
				}
				diffusivity[i] = ( ( Double )nSpeciesInSystem - 1.0 ) / ( mixtureMolarMass * sumDiff[i] );
			}
//			if ( isnan( diffusivity[i] ) ) {
//				fprintf( stderr, "diff[%s] is %g\n", fNames[i], diffusivity[i] );
//				fprintf( stderr, "sumDiff[%s] is %g\n", fNames[i], sumDiff[i] );
//				fprintf( stderr, "sumYi is %g\n", sumYi );
//				fprintf( stderr, "Y[%s] is %g\n", fNames[i], Y[i] );
//				fprintf( stderr, "mixtureMolarMass is %g\n", mixtureMolarMass );
//				fprintf( stderr, "temp is %g\n", temp );
//				for ( j = 0; j < nSpeciesInSystem; ++j ) {
//					fprintf( stderr, "YiOverWi[%s] is %g\n", fNames[j], YiOverWi[j] );
//					fprintf( stderr, "invDij[%s][%s] is %g\n", fNames[i], fNames[j], invDij[i][j] );
//				}
//				exit(2);
//			}
	/*		if ( diffusivity[i] < 0.0 ) {*/
	/*			fprintf( stderr, "diff[%s] is %g\n", fNames[i], diffusivity[i] );*/
	/*		}*/
			
		}
		Free1DArray( YiOverWi );
	}
	else {
		// do nothing
//		++cont0;
	}

//	fprintf( stderr, "cont0 = %d\tcont2 = %d\tcont1 = %d\tcontAll = %d\n", cont0, cont1, cont2, contAll );
}
#endif

#ifdef OPTDTHERM
void T1DSpecies::Compute_DTherm( TFlameNodePtr flameNode, Flag calcNewConst )
{
	int				i, j;
	int				nSpecIn = GetNSpeciesInSystem();
	Double			*M = fMolarMass->vec;
	Double			*lambda = flameNode->conductivity;
	Double			*Y = flameNode->Y[kCurr];
	Double			*diffTherm = flameNode->diffTherm[kCurr];
	Double			*deltaI = flameNode->deltaI;
	Double			diff;
	Double			lambdaOverDeltaI;
	const Double	fact = ( 6.0 / 5.0 * 0.93 - 1.0 ) / RGAS;
	Double			**DThermConst = flameNode->DThermConst;
	
	if ( calcNewConst ) {
		for ( i = 0; i < nSpecIn; ++i ) {
			lambdaOverDeltaI = lambda[i] / deltaI[i];
			for ( j = 0, diff = 0.0; j < nSpecIn; ++j ) {
				DThermConst[i][j] = M[j] / ( M[i] + M[j] ) * ( lambda[j] / deltaI[j] - 
						lambdaOverDeltaI );
				diff += Y[j] * DThermConst[i][j];
			}
			diffTherm[i] = diff * fact * Y[i] * M[i];
		}
	}
	else {
		for ( i = 0; i < nSpecIn; ++i ) {
			for ( j = 0, diff = 0.0; j < nSpecIn; ++j ) {
				diff += Y[j] * DThermConst[i][j];
			}
			diffTherm[i] = diff * fact * Y[i] * M[i];
		}
	}
}

#else

void T1DSpecies::Compute_DTherm( TFlameNodePtr flameNode, Flag /*calcNewConst*/ )
{
	int				i, j;
	int				nSpecIn = GetNSpeciesInSystem();
	Double			*M = fMolarMass->vec;
	Double			*lambda = flameNode->conductivity;
	Double			*Y = flameNode->Y[kCurr];
	Double			*diffTherm = flameNode->diffTherm[kCurr];
	Double			*deltaI = flameNode->deltaI;
	Double			diff;
	Double			lambdaOverDeltaI;
	const Double	fact = ( 6.0 / 5.0 * 0.93 - 1.0 ) / RGAS;
	
	for ( i = 0; i < nSpecIn; ++i ) {
		lambdaOverDeltaI = lambda[i] / deltaI[i];
		for ( j = 0, diff = 0.0; j < nSpecIn; ++j ) {
			diff += Y[j] * M[j] / ( M[i] + M[j] ) * ( lambda[j] / deltaI[j] - 
					lambdaOverDeltaI );
		}
		diffTherm[i] = diff * fact * Y[i] * M[i];
	}
}
#endif

Double T1DSpecies::ComputeOneDiffusionCoeff( int number, Double *Y, Double *invDij, Double mixMolarMass )
{
	int		nSpeciesInSystem = GetNSpeciesInSystem();
	Double	*molarMass = fMolarMass->vec;
	Double	sum = Y[0] * invDij[0] / molarMass[0];		/* j = 0 */
	
	for ( int j = 1; j < nSpeciesInSystem; ++j ) 
			sum += Y[j] * invDij[j] / molarMass[j];
	if ( sum ) {
		return ( 1.0 - Y[number] ) / ( sum * mixMolarMass );
	}
	else {
		//return = temp * sqrt( temp ) * dCoeff[i][i] / ( pressure * omega_D( temp * omegaCoeff[i][i] ) );
		cerr << "warning: TSpecies::ComputeOneDiffusionCoeff cannot handle this case ( Y = 1.0 )" << NEWL;
		return 0.0;
	}
}

void T1DSpecies::PrintProdRateTerms( const char *name, T1DFlamePtr flame )
{
	int i;

	if ( ( i = FindSpecies( name ) ) >= 0 ) {
		PrintProdRateTerms( i, flame );
	}
	else {
		fprintf( stderr, "###warning: no species %s\n", name );
		return;
	}
}

void T1DSpecies::PrintProdRateTerms( int i, T1DFlamePtr flame )
{
	TNewtonPtr	bt = flame->GetSolver()->bt;
	int			k;
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		dummy;
	int			counter;

	int				j;
	T1DReactionPtr	reaction = flame->GetReaction();	
	Double			*nu = fNu[i]->vec;
	int				*usedReactions = fUsedReactions[i]->vec;
	int				nUsedReactions = fNOfUsedReactions->vec[i];
	char			**labels = reaction->GetLabels();
	Double			**reactionRate = reaction->GetReactionRate()->mat;
	Double			Mi = fMolarMass->vec[i];
	Double			prodRate;
	char			reac[128];

	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( flame->GetOutFileBuff(), "%s%s_ProdReac_%d.dout", flame->GetOutputPath(), fNames[i], counter );
	if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) {
		cerr << "#warning: unable to open file" << flame->GetOutFileBuff() << NEWL;
		return;
	}

	
	
//	header
	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "y" );
	for ( j = 0; j < nUsedReactions; ++j ) {
//		fprintf( fp, "\t%-12s", labels[usedReactions[j]] );
//		fprintf( fp, "\t%s: ", labels[usedReactions[j]] );
		reaction->PrintReactionEquation( usedReactions[j], flame->GetSpecies(), reac );
//		fprintf( stderr, "%s\n", reac );
		fprintf( fp, "\t%s: %-s", labels[usedReactions[j]], reac );
	}

	for ( k = 0; k < nGridPoints; ++k ) {
		fprintf( fp, "\n%-12E", x[k] );
		for ( j = 0; j < nUsedReactions; ++j ) {
			prodRate = - Mi * nu[j] * reactionRate[k][usedReactions[j]];
			fprintf( fp, "\t%-12E", prodRate );
		}
	}
	
	fclose( fp );
}

void T1DSpecies::PrintProdCons( TNewtonPtr bt, T1DFlamePtr flame )
{
	int			i, k;
	int			nGridPoints = bt->GetCurrentGridPoints();
	int			nSpecies = GetNSpeciesInSystem();
	Double		**reactionRate = flame->GetReaction()->GetReactionRate()->mat;
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	Double		**source = New2DArray( nGridPoints, nSpecies );
	Double		**sink = New2DArray( nGridPoints, nSpecies );
	
	for ( k = 0; k < nGridPoints; ++k ) {
		CompProdCons( source[k], sink[k], reactionRate[k] );
	}
	
	
	sprintf( flame->GetOutFileBuff(), "%sProductionRates.dout", flame->GetOutputPath() );
	if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) { 
		cerr << "#warning: unable to open file " << flame->GetOutFileBuff() << NEWL;
		return;
	}
	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "eta" );
	for ( i = 0; i < nSpecies; ++i ) {
		fprintf( fp, "\tsource_%-12s\tsink_%-12s\tsum_%-12s", fNames[i], fNames[i], fNames[i] );
	}
	for ( k = 0; k < nGridPoints; ++k ) {
		fprintf( fp, "\n%-9E", x[k] );
		for ( i = 0; i < nSpecies; ++i ) {
			fprintf( fp, "\t%-9E\t%-9E\t%-9E", source[k][i]*1.0e-3, sink[k][i]*1.0e-3
						, ( source[k][i] + sink[k][i] ) * 1.0e-3 );
		}
	}
	fclose( fp );
	Free2DArray( sink );
	Free2DArray( source );
}

void T1DSpecies::PrintDiffusionCoeffs( TNewtonPtr bt )
{
	int			i, k;
	int			nSpeciesInSystem = GetNSpeciesInSystem();
	char		fileName[32];
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	Double		**D = fDiffusivity->mat;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	static int	counter = 0;
	
	sprintf( fileName, "diff%d.dout", ++counter );
	if ( counter >= 10 ) counter = 0;
	if ( !( fp = fopen( fileName, "w" ) ) ) { 
		cerr << "#warning: unable to open file " << fileName << NEWL;
		return;
	}
	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "eta" );
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fprintf( fp, "\tD_%-12s", fNames[i] );
	}
	fprintf( fp, "\n%-9E", left );
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fprintf( fp, "\t%-9E", D[-1][i] );
	}
	for ( k = 0; k < nGridPoints; ++k ) {
		fprintf( fp, "\n%-9E", x[k] );
		for ( i = 0; i < nSpeciesInSystem; ++i ) {
			fprintf( fp, "\t%-9E", D[k][i] );
		}
	}
	fprintf( fp, "\n%-9E", right );
	for ( i = 0; i < nSpeciesInSystem; ++i ) {
		fprintf( fp, "\t%-9E", D[nGridPoints][i] );
	}
	fclose( fp );
}

//void T1DSpecies::CompBinDiffCoeff( TFlameNodePtr flameNode, Double temp, Double *Y, Double press )
#endif // ZEROD
