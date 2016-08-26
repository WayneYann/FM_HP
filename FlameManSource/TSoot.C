// history
// cut fCSootStar


#include "FlameMaster.h"

#undef DIFFFLAME

#define ATTENTION 1.0

#define NEWPOLY
#undef NEWPOLYNUCCOND

#define SIMPLESOOT

#undef NOPAH
#undef PAHFROMA4
#undef A4STEADYSTATE

#ifdef SIMPLESOOT
#define PAHFROMA4
#define NOPAH
#define A4STEADYSTATE
#endif

#ifdef DIFFFLAME
#define NODIFF
#define MAGICSOURCE 1.0
#else
#define NODIFF
#define MAGICSOURCE 1.0
#endif

#undef CUTFRACMOM
#undef CUTLOWER

#define NEWSURFGROWTH
#define WITHTHIRDBODY
#define OXMOM0NEW
#undef DEBUGCOAG
#undef NEWINTERPOLATION
#undef FLUXBC
#undef CHECKDIFFUSIONCOEFF
#undef OLDMECH
#define SMALLSOOT 1.0e-20
#define MAGIC 0.33 //Radiation in only one direction: 1/3

#ifndef NOPAH
#undef PAHFROMA4
#endif

void TSoot::InitTSoot( TInputDataPtr input )
{
	CounterPtr	counter = input->GetCounter();
	ReactionPtr	reaction = input->GetPAHReactions();	
	ReactionPtr	sootReac = input->GetSootReactions();	
	int 		nPAHReactions = counter->pahReactions;
	int 		nSootReactions = counter->sootReactions;
	int			i;

//	fprintf( stderr, "###ATTENTION: fAlpha = %g\n", fAlpha );

	fNucleation = input->fNucleation;
	fCondensation = input->fCondensation;
	fCoagulation = input->fCoagulation;
	fSurfaceGrowth = input->fSurfaceGrowth;
	fSurfaceOxidation = input->fSurfaceOxidation;
	fThermoPhoresis = input->fThermoPhoresis;
	fCoagFact = input->fCoagFact;
	fLewis1 = 1.0;
	fSootRadiation = input->fSootRadiation;
	fSootUpdateProdRates = input->fSootUpdateProdRates;
	fSizeDepDiff = input->fSizeDepDiff;
	fSurfDepCoag = input->fSurfDepCoag;

	// init pah reactions
	if ( !nPAHReactions ) {
#ifdef NOPAH
		nPAHReactions = 1;
#else
		cerr << "#error: no reactions for PAH's specified" << NEWL;
		exit( 2 );
#endif
	}
	fA = NewVector( nPAHReactions );
#ifdef NOPAH
	fA->len = 0;
#endif
	fN = NewVector( nPAHReactions );
	fEOverRgas = NewVector( nPAHReactions );

	fSpeciesNumber = new IntVectorPtr[nPAHReactions];
   	if ( !fSpeciesNumber ) FatalError( "memory allocation of TSoot failed" );
	fNu = new VectorPtr[nPAHReactions];
   	if ( !fNu ) FatalError( "memory allocation of TSoot failed" );
#ifdef NOPAH
#	ifdef PAHFROMA4
	fPAHSpeciesIndex = NewIntVector( 2 );
#	else
	fPAHSpeciesIndex = NewIntVector( fNPAHMolecules+3 );
#	endif
#else
	for ( i = 0; i < nPAHReactions; ++i ) {
		fSpeciesNumber[i] = NewIntVector( reaction[i].numberOfSpecies );
	   	fNu[i] = NewVector( reaction[i].numberOfSpecies );
	}
	fPAHSpeciesIndex = NewIntVector( fNPAHMolecules+1 );
#endif


	fRateCoefficients = NewVector( nPAHReactions );
	fRedRateCoeffs = NewVector( kLastRedReacCoeff );

	fZ = NewVector( fNPAHMolecules + 1 );
	fZ->vec[5] = fZ->vec[6] = 1.0;
	fDenom = NewVector( fNPAHMolecules + 1 );
	fF = NewVector( kLastRedReacCoeff );

	fACoeff = NewMatrix( fNPAHMolecules+1, kRestPoly+1, kColumnPointers );
	fBCoeff = NewMatrix( fNPAHMolecules+1, kRestPoly+1, kColumnPointers );
	fBOHCoeff = NewMatrix( fNPAHMolecules+1, kRestPoly+1, kColumnPointers );

	fDelta = NewVector( fNPAHMolecules+1 );

	fFracMom = NewVector( kLastFracMoment );
	fPhi = NewVector( kLastPhi );

	fFracMomPAH = NewVector( kLastFracMoment );
	fPhiPAH = NewVector( kLastPhi );

	// init soot reactions
	if ( !nSootReactions ) {
		cerr << "error: number of soot reactions is zero" << NEWL;
		exit(2);
	}
	fASoot = NewVector( nSootReactions );
	fNSoot = NewVector( nSootReactions );
	fEOverRSoot = NewVector( nSootReactions );
	fSpecNumSoot = new IntVectorPtr[nSootReactions];
   	if ( !fSpecNumSoot ) FatalError( "memory allocation of TSoot failed" );
	fNuSoot = new VectorPtr[nSootReactions];
   	if ( !fNuSoot ) FatalError( "memory allocation of TSoot failed" );
	for ( i = 0; i < nSootReactions; ++i ) {
		fSpecNumSoot[i] = NewIntVector( sootReac[i].numberOfSpecies );
	   	fNuSoot[i] = NewVector( sootReac[i].numberOfSpecies );
	}
	fSootRateCoeffs = NewVector( nSootReactions );
	fSootReactionRate = NewVector( nSootReactions );

/*	// small PAH's*/
/*	fSmallPAHIndTab = new int[kNSmallSoot];*/
/*	fSmallPAHRings = new Double[kNSmallSoot];*/
/*	fSmallPAHNumDens = new Double[kNSmallSoot];*/
/*	fSmallPAHMoments = new Double[fNSootMoments];*/
/*	InitSmallPAHIndTab( input );*/

	if ( !( fLabels = new String[nPAHReactions] ) ) FatalError( "memory allocation of TSoot failed" );
	if ( !( fSootLabels = new String[nSootReactions] ) ) FatalError( "memory allocation of TSoot failed" );
	if ( !( fPAHSymbol = new char[strlen( input->fPAHSymbol )+1] ) ) FatalError( "memory allocation of TSoot failed" );
	strcpy( fPAHSymbol, input->fPAHSymbol );
	if ( !( fSootSymbol = new char[strlen( input->fSootSymbol )+1] ) ) FatalError( "memory allocation of TSoot failed" );
	strcpy( fSootSymbol, input->fSootSymbol );

	FillTSoot( input );
//	fprintf( stderr, "check soot reactions\n" );
	CheckSootReactions();
#ifdef NOPAH
#else
	CheckPAHReactions();
#endif
	
	int	*ind = fPAHSpeciesIndex->vec;

#ifdef OLDMECH
	if ( ( ind[1] = input->FindSpecies( "A4-C18H10" ) ) == -1 ) {
		cerr << "error: can't find the molecule A4-C18H10" << NEWL;
		exit( 2 );
	}
	if ( ( ind[2] = input->FindSpecies( "A4M-C18H9" ) ) == -1 ) {
		cerr << "error: can't find the molecule A4M-C18H9" << NEWL;
		exit( 2 );
	}
	if ( ( ind[3] = input->FindSpecies( "A4C2H-C20H10" ) ) == -1 ) {
		cerr << "error: can't find the molecule A4C2H-C20H10" << NEWL;
		exit( 2 );
	}
	if ( ( ind[4] = input->FindSpecies( "A4C2HS-C20H9" ) ) == -1 ) {
		cerr << "error: can't find the molecule A4C2HS-C20H9" << NEWL;
		exit( 2 );
	}
	if ( ( ind[5] = input->FindSpecies( "A5M-C22H11" ) ) == -1 ) {
		cerr << "error: can't find the molecule A5M-C22H11" << NEWL;
		exit( 2 );
	}
	if ( ( ind[6] = input->FindSpecies( "A5-C22H12" ) ) == -1 ) {
		cerr << "error: can't find the molecule A5-C22H12" << NEWL;
		exit( 2 );
	}
#else
	if ( ( ind[1] = input->FindSpecies( "A4-C18H10" ) ) == -1 ) {
		cerr << "error: can't find the molecule A4-C18H10" << NEWL;
		exit( 2 );
	}
#	ifdef PAHFROMA4
#		ifndef A4STEADYSTATE
	if ( input->GetSpecies()[ind[1]].isSteadyState ) {
		cerr << "error: A4-C18H10 is steady state, while PAHs should be computed from A4" << NEWL;
		exit( 2 );
	}
#		else
	if ( !input->GetSpecies()[ind[1]].isSteadyState ) {
		cerr << "error: A4-C18H10 is not steady state, although it should be" << NEWL;
		exit( 2 );
	}
#		endif
#	else
	if ( ( ind[2] = input->FindSpecies( "A4--C18H9" ) ) == -1 ) {
		cerr << "error: can't find the molecule A4--C18H9" << NEWL;
		exit( 2 );
	}
	if ( ( ind[3] = input->FindSpecies( "A4C2H-C20H10" ) ) == -1 ) {
		cerr << "error: can't find the molecule A4C2H-C20H10" << NEWL;
		exit( 2 );
	}
	if ( ( ind[4] = input->FindSpecies( "A4C2H*-C20H9" ) ) == -1 ) {
		cerr << "error: can't find the molecule A4C2H*-C20H9" << NEWL;
		exit( 2 );
	}
	if ( ( ind[5] = input->FindSpecies( "A5--C22H11" ) ) == -1 ) {
		cerr << "error: can't find the molecule A5--C22H11" << NEWL;
		exit( 2 );
	}
	if ( ( ind[6] = input->FindSpecies( "A5-C22H12" ) ) == -1 ) {
		cerr << "error: can't find the molecule A5-C22H12" << NEWL;
		exit( 2 );
	}
#		ifdef NOPAH
	if ( ( ind[7] = input->FindSpecies( "A6-C24H12" ) ) == -1 ) {
		cerr << "error: can't find the molecule A6-C24H12" << NEWL;
		exit( 2 );
	}
	if ( ( ind[8] = input->FindSpecies( "A6--C24H11" ) ) == -1 ) {
		cerr << "error: can't find the molecule A6--C24H11" << NEWL;
		exit( 2 );
	}
#		endif
#	endif
#endif

	if ( ( f_H = input->FindSpecies( "H" ) ) == -1 ) {
		cerr << "error: can't find the molecule H" << NEWL;
		exit( 2 );
	}
	if ( ( f_H2 = input->FindSpecies( "H2" ) ) == -1 ) {
		cerr << "error: can't find the molecule H2" << NEWL;
		exit( 2 );
	}
	if ( ( f_O2 = input->FindSpecies( "O2" ) ) == -1 ) {
		cerr << "error: can't find the molecule O2" << NEWL;
		exit( 2 );
	}
	if ( ( f_OH = input->FindSpecies( "OH" ) ) == -1 ) {
		cerr << "error: can't find the molecule OH" << NEWL;
		exit( 2 );
	}
	if ( ( f_CO = input->FindSpecies( "CO" ) ) == -1 ) {
		cerr << "error: can't find the molecule f_CO" << NEWL;
		exit( 2 );
	}
	if ( ( f_CH = input->FindSpecies( "CH" ) ) == -1 ) {
		cerr << "error: can't find the molecule f_CH" << NEWL;
		exit( 2 );
	}
	if ( ( f_CHO = input->FindSpecies( "CHO" ) ) == -1 ) {
		if ( ( f_CHO = input->FindSpecies( "HCO" ) ) == -1 ) {
			cerr << "error: can't find the molecule f_CHO or f_HCO" << NEWL;
			exit( 2 );
		}
	}
	if ( ( f_H2O = input->FindSpecies( "H2O" ) ) == -1 ) {
		cerr << "error: can't find the molecule H2O" << NEWL;
		exit( 2 );
	}
	if ( ( f_C2H = input->FindSpecies( "C2H" ) ) == -1 ) {
		cerr << "error: can't find the molecule C2H" << NEWL;
		exit( 2 );
	}
	if ( ( f_C2H2 = input->FindSpecies( "C2H2" ) ) == -1 ) {
		cerr << "error: can't find the molecule C2H2" << NEWL;
		exit( 2 );
	}
#ifdef NEWPOLY
	if ( ( f_A3R5M = input->FindSpecies( "A3R5--C16H9" ) ) == -1 ) {
		cerr << "error: can't find the molecule A3R5--C16H9" << NEWL;
		exit( 2 );
	}
	if ( ( f_A3R5AC = input->FindSpecies( "A3R5AC-C18H11" ) ) == -1 ) {
		cerr << "error: can't find the molecule A3R5AC-C18H11" << NEWL;
	}
#else
	if ( ( f_A3R5AC = input->FindSpecies( "A3R5AC-C18H11" ) ) == -1 ) {
		cerr << "error: can't find the molecule A3R5AC-C18H11" << NEWL;
		exit( 2 );
	}
#endif
	if ( ( f_A1 = input->FindSpecies( "A1-C6H6" ) ) == -1 ) {
		cerr << "error: can't find the molecule A1-C6H6" << NEWL;
		exit( 2 );
	}
	if ( ( fFirstPAH = input->FindSpecies( "A4-C18H10" ) ) == -1 ) {
		cerr << "error: can't find the molecule A4-C18H10" << NEWL;
		exit( 2 );
	}
	for ( i = 1; i < fPAHSpeciesIndex->len ; ++i ) {
		if ( ( ind[i] - fFirstPAH + 1 ) != i ) {
			cerr << "error: PAH molecule no. " << i << TAB << (ind[i] - fFirstPAH + 1) << " is not in right order" << NEWL;
		}
	}
	
	fOffsetMoments = -1; // value of fOffsetMoments is set later

//	InitSootReactions();
}

TSoot::~TSoot( void )
{
	int	i;
	Double	nPAHReactions = GetNPAHReactions();
	Double	nSootReactions = GetNSootReactions();

	delete fSootSymbol;
	for ( i = 0; i < nSootReactions; ++i ) {
		delete fSootLabels[i];
	}
	delete fSootLabels;

	delete fPAHSymbol;
	for ( i = 0; i < nPAHReactions; ++i ) {
		delete fLabels[i];
	}
	delete fLabels;

/*	delete fSmallPAHMoments;*/
/*	delete fSmallPAHNumDens;*/
/*	delete fSmallPAHRings;*/
/*	delete fSmallPAHIndTab;*/
	
	DisposeVector( fSootReactionRate );
	DisposeVector( fSootRateCoeffs );
	for ( i = 0; i < GetNSootReactions(); ++i ) {
		DisposeVector ( fNuSoot[i] );
		DisposeIntVector( fSpecNumSoot[i] );
	}
	delete fNuSoot;
	delete fSpecNumSoot;
	DisposeVector( fEOverRSoot );
	DisposeVector( fNSoot );
	DisposeVector( fASoot );

	DisposeVector( fPhiPAH );
	DisposeVector( fFracMomPAH );

	DisposeVector( fPhi );
	DisposeVector( fFracMom );

	DisposeVector( fDelta );

	DisposeMatrix( fBOHCoeff );
	DisposeMatrix( fBCoeff );
	DisposeMatrix( fACoeff );

	DisposeVector( fF );
	DisposeVector( fDenom );
	DisposeVector( fZ );

	DisposeVector( fRedRateCoeffs );
	DisposeVector( fRateCoefficients );
	
	DisposeIntVector( fPAHSpeciesIndex );

	for ( i = 0; i < nPAHReactions; ++i ) {
		DisposeVector ( fNu[i] );
		DisposeIntVector( fSpeciesNumber[i] );
	}
	delete fNu;
	delete fSpeciesNumber;

	DisposeVector( fEOverRgas );
	DisposeVector( fN );
	DisposeVector( fA );
}

void T0DSoot::InitT0DSoot( TInputDataPtr input )
{
	int 		nPAHReactions = input->GetCounter()->pahReactions;

	fPij = NewMatrix( fNPAHMolecules+1, fNStages, kColumnPointers );

#ifdef NOPAH
	fSumPi = NewVector( fNPAHMolecules+2 );
#else
	fSumPi = NewVector( fNPAHMolecules );
#endif
	fPAHMoments = NewVector( fNPAHMoments );
	fMoments = NewVector( fNSootMoments );

	fReactionRate = NewVector( nPAHReactions );
}

T0DSoot::~T0DSoot( void )
{
	DisposeVector( fReactionRate );

	DisposeVector( fMoments );
	DisposeVector( fPAHMoments );
	DisposeVector( fSumPi );
	DisposeMatrix( fPij );
}

/*void TSoot::InitSmallPAHIndTab( TInputDataPtr input )*/
/*{*/
/*	int count = 0;*/
/**/
/*#ifdef OLDMECH*/
/*	SetSmallPAH( kA1, "A1-C6H6", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA1M, "A1M-C6H5", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA1C2H, "A1C2H-C8H6", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA1C2HM, "A1C2HM-C8H5", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA1C2HS, "A1C2HS-C8H5", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA1C2HAC, "A1C2HAC-C10H7", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA2MX, "A2MX-C10H7", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA2, "A2-C10H8", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA2R5, "A2R5-C12H8", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA2R5M, "A2R5M-C12H7", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA2R5C2H, "A2R5C2H-C14H8", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA2R5C2HS, "A2R5C2HS-C14H7", 2.0, input ); ++count;*/
/*	SetSmallPAH( kANC2HAC, "ANC2HAC-C16H9", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA3R5M, "A3R5M-C16H9", 3.0, input ); ++count;*/
/*	SetSmallPAH( kA3R5, "A3R5-C16H10", 3.0, input ); ++count;*/
/*	SetSmallPAH( kA3R5AC, "A3R5AC-C18H11", 3.0, input ); ++count;*/
/*#else*/
/*	SetSmallPAH( kA1, "A1-C6H6", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA1M, "A1--C6H5", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA1C2H, "A1C2H-C8H6", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA1C2HM, "A1C2H--C8H5", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA1C2HS, "A1C2H*-C8H5", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA1C2HAC, "A1C2HAC-C10H7", 1.0, input ); ++count;*/
/*	SetSmallPAH( kA2MX, "A2-X-C10H7", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA2, "A2-C10H8", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA2R5, "A2R5-C12H8", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA2R5M, "A2R5--C12H7", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA2R5C2H, "A2R5C2H-C14H8", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA2R5C2HS, "A2R5C2H*-C14H7", 2.0, input ); ++count;*/
/*	SetSmallPAH( kANC2HAC, "ANC2HAC-C16H9", 2.0, input ); ++count;*/
/*	SetSmallPAH( kA3R5M, "A3R5--C16H9", 3.0, input ); ++count;*/
/*	SetSmallPAH( kA3R5, "A3R5-C16H10", 3.0, input ); ++count;*/
/*	SetSmallPAH( kA3R5AC, "A3R5AC-C18H11", 3.0, input ); ++count;*/
/*#endif*/
/**/
/*	if ( count != kNSmallSoot ) {*/
/*		cerr << "#error in InitSmallPAHIndTab: count = " << count */
/*					<< TAB << "kNSmallSoot = " << (int)kNSmallSoot << NEWL;*/
/*	}*/
/*}*/
/**/
/*void TSoot::SetSmallPAH( SmallPAH which, char *name, Double nOfRings, TInputDataPtr input )*/
/*{*/
/*	if ( ( fSmallPAHIndTab[which] = input->FindSpecies( name ) ) == -1 ) {*/
/*		cerr << "warning: can't find the molecule " << name << NEWL;*/
/*//		exit( 2 );*/
/*	}*/
/*	else {*/
/*		fSmallPAHRings[which] = nOfRings;*/
/*	}*/
/*}*/

void TSoot::CheckSootReactions( void )
{
	int 	error = -1;
	int		offset = strlen( fSootSymbol );
	
#ifdef NEWSURFGROWTH
	if ( strcmp( fSootLabels[ks7f] + offset, "7f" ) ) error = ks7f; 
	if ( strcmp( fSootLabels[ks7b] + offset, "7b" ) ) error = ks7b; 
	if ( strcmp( fSootLabels[ks9f] + offset, "9f" ) ) error = ks9f; 
	if ( strcmp( fSootLabels[ks9b] + offset, "9b" ) ) error = ks9b; 
	if ( strcmp( fSootLabels[ks111] + offset, "111" ) ) error = ks111; 
	if ( strcmp( fSootLabels[ks112] + offset, "112" ) ) error = ks112; 
	if ( strcmp( fSootLabels[ks13f] + offset, "13f" ) ) error = ks13f; 
	if ( strcmp( fSootLabels[ks13b] + offset, "13b" ) ) error = ks13b; 
#else
	if ( strcmp( fSootLabels[ks9] + offset, "9" ) ) error = ks9; 
	if ( strcmp( fSootLabels[ks11] + offset, "11" ) ) error = ks11; 
#endif

	if ( strcmp( fSootLabels[ks8f] + offset, "8f" ) ) error = ks8f; 
	if ( strcmp( fSootLabels[ks8b] + offset, "8b" ) ) error = ks8b; 

	if ( strcmp( fSootLabels[ks10f] + offset, "10f" ) ) error = ks10f; 
	if ( strcmp( fSootLabels[ks10b] + offset, "10b" ) ) error = ks10b; 

	if ( strcmp( fSootLabels[ks12] + offset, "12" ) ) error = ks12; 

	if ( error > -1 ) {
		cerr << "###error in soot reaction " << fSootLabels[error] << NEWL;
		cerr << "ks12 =  " << (int)ks12 << NEWL;
		cerr << "fSootLabels[ks12] + offset =  " << fSootLabels[ks12] + offset << NEWL;
		exit( 2 );
	}
	
#ifdef NEWSURFGROWTH
//	if ( !fSurfaceGrowth ) {
//		A[ks7f] = A[ks7b] = A[ks8f] = A[ks8b] = A[ks9f] = A[ks9b] = A[ks10f] = A[ks10b] = A[ks13f] = A[ks13b] 
//		= n[ks7f] = n[ks7b] = n[ks8f] = n[ks8b] = n[ks9f] = n[ks9b] = n[ks10f] = n[ks10b] = n[ks13f] = n[ks13b]
//		= E[ks7f] = E[ks7b] = E[ks8f] = E[ks8b] = E[ks9f] = E[ks9b] = E[ks10f] = E[ks10b] = E[ks13f] = E[ks13b]
//		= 0.0;
//	}

//	if ( !fSurfaceOxidation ) {
//		A[ks111] = A[ks112] = A[ks12] 
//		= n[ks111] = n[ks112] = n[ks12]
//		= E[ks111] = E[ks112] = E[ks12] = 0.0;
//	}
#else
//	if ( !fSurfaceGrowth ) {
//		A[ks8f] = A[ks8b] = A[ks9] = A[ks10f] = A[ks10b] 
//				= n[ks8f] = n[ks8b] = n[ks9] = n[ks10f] = n[ks10b]
//				= E[ks8f] = E[ks8b] = E[ks9] = E[ks10f] = E[ks10b] = 0.0;
//	}

//	if ( !fSurfaceOxidation ) {
//		A[ks11] = A[ks12] 
//				= n[ks11] = n[ks12]
//				= E[ks11] = E[ks12] = 0.0;
//	}
#endif
}

void TSoot::PrintSootReactions( TFlamePtr flame, TSpeciesPtr species )
{
	char	**names = species->GetNames();
	int		*ind;
	FILE	*fp = flame->GetOutfile( "sootreacs", flame->kText );

	for ( int i = 0; i < GetNSootReactions(); ++i ) {
		ind = fSpecNumSoot[i]->vec;
		fprintf( fp, "reaction %d:\n", i );
		for ( int j = 0; j < fSpecNumSoot[i]->len; ++j ) {
			fprintf( fp, "\t %g %s", fNuSoot[i]->vec[j], names[ind[j]] );
		}
		fprintf( fp, "\n\nA = %g\nn = %g\nE = %g\n\n\n"
					, fASoot->vec[i], fNSoot->vec[i], fEOverRSoot->vec[i] * RGAS );
	}

	fclose( fp );
}

#ifndef NEWSURFGROWTH
/*void TSoot::InitSootReactions( void )
{
	Double	*a = fASoot->vec;
	Double	*n = fNSoot->vec;
	Double	*eOverR = fEOverRSoot->vec;
	
	a[ks8f] = 7.900e10;
	n[ks8f] = 0.0;
	eOverR[ks8f] = 41.70e6 / RGAS;

	a[ks8b] = 5.683e9;
	n[ks8b] = 0.0;
	eOverR[ks8b] = 20.67e6 / RGAS;

	a[ks9] = 1.000e10;
	n[ks9] = 0.0;
	eOverR[ks9] = 0.0;

	a[ks10f] = 1.000e10;
	n[ks10f] = 0.0;
	eOverR[ks10f] = 0.0;

	a[ks10b] = 2.000e18;
	n[ks10b] = 0.0;
	eOverR[ks10b] = 290.9e6 / RGAS;

	if ( fSurfaceOxidation ) {
		a[ks11] = 1.0e10;
		n[ks11] = 0.0;
		eOverR[ks11] = 0.0;
	
		a[ks12] = 1.3e10;
		n[ks12] = 0.0;
		eOverR[ks12] = 46.0e6 / RGAS;
	}
	else {
		a[ks11] = 0.0;
		n[ks11] = 0.0;
		eOverR[ks11] = 0.0;
	
		a[ks12] = 0.0;
		n[ks12] = 0.0;
		eOverR[ks12] = 0.0;
	}
}
*/
#endif

void TSoot::FillTSoot( TInputDataPtr input )
{
	int				i, j;
	ReactionPtr		reaction = input->GetPAHReactions();
	DimensionPtr	dimension = input->GetDimension();
	int 			numberOfDimensions = input->GetCounter()->dimensions;
	Double 			unitsConverterA = 0.0;
	Double 			unitsConverterE = 0.0;
	Double			reacOrderConv = 0.0;
	Double			tempExpConv = 0.0;
	int				orderOfReaction;
	const Double	moleToKmole = 1000.0;

	// set units converter
	// it is assumed that only A is a function of the order of reaction 
	// and the temperature exponent
	for ( i = 0; i < numberOfDimensions; ++i ) {
		if ( dimension[i].name[0] == 'A' ) {
			if ( unitsConverterA ) {
				cerr << "#error: units converter for 'A' doubly defined or zero" << endl;
				exit( 2 );
			}
			else {
				unitsConverterA = dimension[i].value;
				if ( dimension[i].orderOfReaction ) {
					reacOrderConv = dimension[i].orderOfReacValue;
				}
				if ( dimension[i].tempExponent ) {
					tempExpConv = dimension[i].tempExpValue;
				}
			}
		}
		else if ( dimension[i].name[0] == 'E' ) {
			if ( unitsConverterE ) {
				cerr << "#error: units converter for 'E' doubly defined or zero" << endl;
				exit( 2 );
			}
			else {
				unitsConverterE = dimension[i].value;
			}
		}
		else {
			cerr << "# warning: i don't need units for " << dimension->name << endl;
		}
	}
	
	// set values for A, N, E, including the conversion to SI-units and mole to kmole
	// for PAH reactions
	for ( i = 0; i < GetNPAHReactions(); ++i ) {
		// first get label
		fLabels[i] = new char[strlen( reaction[i].label )+1];
		if ( !fLabels[i] ) FatalError( "memory allocation of TSoot failed" );
		strcpy( fLabels[i], reaction[i].label );
		
	   	// then determine order of reaction
		orderOfReaction = ( reaction[i].withThirdBody ) ? 1: 0;
	   	for ( j = 0; j < reaction[i].numberOfSpecies; ++j ) {
			if ( reaction[i].speciesCoeff[j] > 0.0 ) {
			   	orderOfReaction += ( int )reaction[i].speciesCoeff[j];
			}
		}
		if ( orderOfReaction > 3 && !reaction[i].partialEquilibrium ) {
			cerr << "#warning: order of reaction no. " << i << " is " << orderOfReaction << NEWL;
		}

		// then fill arrays
		// first convert A and E to SI-units, then change mole to kmole
	   	fA->vec[i] = unitsConverterA * myPow( reacOrderConv, orderOfReaction ) 
						* myPow( tempExpConv, reaction[i].n ) * reaction[i].a 
						* myPow( moleToKmole, ( orderOfReaction-1 ) );
		fN->vec[i] = reaction[i].n;
		fEOverRgas->vec[i] = unitsConverterE * reaction[i].e / RGAS * 1000.0;
		
		int	speciesCounter = 0;
		for ( j = 0; j < fNu[i]->len; ++ j ) {
			if ( reaction[i].speciesCoeff[j] > 0.0 ) {
				fSpeciesNumber[i]->vec[speciesCounter] = reaction[i].speciesNumber[j];
		 	  	fNu[i]->vec[speciesCounter] = reaction[i].speciesCoeff[j];
				++speciesCounter;
			}
		}
		for ( j = 0; j < fNu[i]->len; ++ j ) {
			if ( reaction[i].speciesCoeff[j] < 0.0 ) {
				fSpeciesNumber[i]->vec[speciesCounter] = reaction[i].speciesNumber[j];
		 	  	fNu[i]->vec[speciesCounter] = reaction[i].speciesCoeff[j];
				++speciesCounter;
			}
		}
	}
	
	// set values for A, N, E, including the conversion to SI-units and mole to kmole
	// for Soot reactions
	reaction = input->GetSootReactions();
	Double	*a = fASoot->vec;
	Double	*n = fNSoot->vec;
	Double	*eOverR = fEOverRSoot->vec;
	for ( i = 0; i < GetNSootReactions(); ++i ) {
		// first get label
		fSootLabels[i] = new char[strlen( reaction[i].label )+1];
		if ( !fSootLabels[i] ) FatalError( "memory allocation of TSoot failed" );
		strcpy( fSootLabels[i], reaction[i].label );
		
	   	// then determine order of reaction
		orderOfReaction = ( reaction[i].withThirdBody ) ? 1: 0;
	   	for ( j = 0; j < reaction[i].numberOfSpecies; ++j ) {
			if ( reaction[i].speciesCoeff[j] > 0.0 ) {
			   	orderOfReaction += ( int )reaction[i].speciesCoeff[j];
			}
		}
		if ( orderOfReaction > 3 && !reaction[i].partialEquilibrium ) {
			cerr << "#warning: order of reaction no. " << i << " is " << orderOfReaction << NEWL;
		}

		// then fill arrays
		// first convert A and E to SI-units, then change mole to kmole
	   	a[i] = unitsConverterA * myPow( reacOrderConv, orderOfReaction ) 
						* myPow( tempExpConv, reaction[i].n ) * reaction[i].a 
						* myPow( moleToKmole, ( orderOfReaction-1 ) );
		n[i] = reaction[i].n;
		eOverR[i] = unitsConverterE * reaction[i].e / RGAS * 1000.0;
		
		int	speciesCounter = 0;
		for ( j = 0; j < fNuSoot[i]->len; ++ j ) {
			if ( reaction[i].speciesCoeff[j] > 0.0 ) {
				fSpecNumSoot[i]->vec[speciesCounter] = reaction[i].speciesNumber[j];
		 	  	fNuSoot[i]->vec[speciesCounter] = reaction[i].speciesCoeff[j];
				++speciesCounter;
			}
		}
		for ( j = 0; j < fNuSoot[i]->len; ++ j ) {
			if ( reaction[i].speciesCoeff[j] < 0.0 ) {
				fSpecNumSoot[i]->vec[speciesCounter] = reaction[i].speciesNumber[j];
		 	  	fNuSoot[i]->vec[speciesCounter] = reaction[i].speciesCoeff[j];
				++speciesCounter;
			}
		}
	}
}

void TSoot::ComputeSootRateCoeffs( Double *k, Double temp, TReactionPtr reaction )
{
	Double	*a = fASoot->vec;
	Double	*n = fNSoot->vec;
	Double	*eOverR = fEOverRSoot->vec;
	
	for ( int i = 0; i < GetNSootReactions(); ++i ) {
		reaction->ComputeRateCoefficient( temp, k[i], a[i], n[i], eOverR[i] );
	}
#ifdef NEWSURFGROWTH
#	ifdef WITHTHIRDBODY
	k[ks10f] /= 7.0e-3;
	k[ks10b] /= 7.0e-3;
#	endif
#endif
}

Double TSoot::Nucleation( int i, Double temp, Double *pahMoments )
{
	switch( i ) {
		case 0: 
			return 0.5 * GetBeta( temp, pahMoments ) * pahMoments[0] * pahMoments[0];
		case 1:
			return GetBeta( temp, pahMoments ) * pahMoments[0] * pahMoments[1];
		case 2:
			return GetBeta( temp, pahMoments ) * ( pahMoments[0] * pahMoments[2] 
								+ pahMoments[1] * pahMoments[1] );
		case 3:
			return GetBeta( temp, pahMoments ) * ( pahMoments[0] * pahMoments[3] 
								+ 3.0 * pahMoments[1] * pahMoments[2] );
		default:
			cerr << "#error: no need to compute alpha_" << i << NEWL;
			exit( 2 );
	}
	
	return 0.0;
}

Double TSoot::NucleationNow( int i, Double temp, Double *Y, Double rho, Double *molarMass )
{	
	switch( i ) {
		case 0: 
			return 0.5 * 8.16 * GetC( temp ) * rho * Y[fFirstPAH] / molarMass[fFirstPAH] * rho * Y[fFirstPAH] / molarMass[fFirstPAH];
		case 1:
			return 18.0 * 0.5 * 8.16 * GetC( temp ) * rho * Y[fFirstPAH] / molarMass[fFirstPAH] * rho * Y[fFirstPAH] / molarMass[fFirstPAH];
		case 2:
			return 324.0 * 0.5 * 8.16 * GetC( temp ) * rho * Y[fFirstPAH] / molarMass[fFirstPAH] * rho * Y[fFirstPAH] / molarMass[fFirstPAH];
		case 3:
			return 5832.0 * 0.5 * 8.16 * GetC( temp ) * rho * Y[fFirstPAH] / molarMass[fFirstPAH] * rho * Y[fFirstPAH] / molarMass[fFirstPAH];
		default:
			cerr << "#error: no need to compute alpha_" << i << NEWL;
			exit( 2 );
	}
	
	return 0.0;
}

Double TSoot::NucleationNew( int i, Double temp, Double *pahMoments )
{	

#ifdef NEWPOLYNUCCOND
	switch( i ) {
		case 0: 
			return 0.5 * GetBeta( temp, pahMoments ) * pahMoments[0] * pahMoments[0];
		case 1:
			return GetBeta( temp, pahMoments ) * pahMoments[0] * pahMoments[1];
		case 2:
			return GetBeta( temp, pahMoments ) * ( pahMoments[0] * pahMoments[2] 
								+ pahMoments[1] * pahMoments[1] );
		case 3:
			return GetBeta( temp, pahMoments ) * ( pahMoments[0] * pahMoments[3] 
								+ 3.0 * pahMoments[1] * pahMoments[2] );
		default:
			cerr << "#error: no need to compute alpha_" << i << NEWL;
			exit( 2 );
	}
	
	return 0.0;
#else
	switch( i ) {
#	ifdef NOPAH
#		ifdef PAHFROMA4
		case 0: 
			return 8.16 * GetC( temp ) * pahMoments[0] * pahMoments[0];
		case 1:
			return 18.0 * 8.16 * GetC( temp ) * pahMoments[0] * pahMoments[0];
		case 2:
			return 324.0 * 8.16 * GetC( temp ) * pahMoments[0] * pahMoments[0];
		case 3:
			return 5832.0 * 8.16 * GetC( temp ) * pahMoments[0] * pahMoments[0];
#		else
			fprintf( stderr, "#error in TSoot::NucleationNew: not yet implemented\n" );
			exit( 2 );
#		endif
#	else
		case 0: 
		  return 0.5 * GetPhiPAH( 0, 0, temp, pahMoments );//changed by GB
		case 1:
		  return GetPhiPAH( 0, 1, temp, pahMoments );//changed by GB 
		case 2:
			return GetPhi( 0, 2, temp, pahMoments ) + GetPhi( 1, 1, temp, pahMoments );
		case 3:
			return GetPhi( 0, 3, temp, pahMoments ) + 3.0 * GetPhi( 1, 2, temp, pahMoments );
#	endif
		default:
			cerr << "#error: no need to compute alpha_" << i << NEWL;
			exit( 2 );
	}
	
	return 0.0;
#endif
}

Double TSoot::FractionalMoment2( Double fracIndex, int first, int second, Double *moments )
{
/*	if ( fracIndex < 0.0 || fracIndex > 1.0 ) {*/
/*		return 0.0;*/
/*	}*/
	if ( moments[first] < SMALLSOOT || moments[second] < SMALLSOOT ) {
		return 0.0;
	}

#ifdef CUTFRACMOM
#	ifdef CUTLOWER
	return MAX( exp( log( MAX( moments[first], SMALLSOOT ) )
						* GetAlphaI2( first, fracIndex, first, second )
				+ log( MAX( moments[second], SMALLSOOT ) ) 
						* GetAlphaI2( second, fracIndex, first, second ) )
			, MAX( moments[first], SMALLSOOT ) );
#	else
	return MAX( MIN( exp( log( MAX( moments[first], SMALLSOOT ) )
						* GetAlphaI2( first, fracIndex, first, second )
				+ log( MAX( moments[second], SMALLSOOT ) ) 
						* GetAlphaI2( second, fracIndex, first, second ) )
			, moments[second] ), MAX( moments[first], SMALLSOOT ) );
#	endif
#else
	return exp( log( MAX( moments[first], SMALLSOOT ) ) 
						* GetAlphaI2( first, fracIndex, first, second )
				+ log( MAX( moments[second], SMALLSOOT ) ) 
						* GetAlphaI2( second, fracIndex, first, second ) );
#endif
}

Double TSoot::FracMom2( Double fracIndex, Double *moments )
{
#ifdef NEWINTERPOLATION
	int	second;

	if ( fracIndex < 0.0 ) {
		second = 1;
	}
	else if ( fracIndex > fNSootMoments-1 ) {
		second = fNSootMoments-1;
	}
	else {
		second = ( int )fracIndex+1;
	}

	return FractionalMoment2( fracIndex, 0, second, moments );
#else
	int	first;
	
	if ( fracIndex < 0.0 ) {
		first = 0;
	}
	else if ( fracIndex > fNSootMoments-1 ) {
		first = fNSootMoments-2;
	}
	else {
		first = ( int )fracIndex;
	}

	return FractionalMoment2( fracIndex, first, first+1, moments );
#endif
}

/*Double TSoot::GetAlphaI( int i, Double fracIndex )
{
	const Double oneOverSix = -1.0 / 6.0;
	switch ( i ) {
		case 0:
			return oneOverSix * ( 1.0 - fracIndex ) * ( fracIndex - 2.0 ) * ( fracIndex - 3.0 );
		case 1:
			return 0.5 * fracIndex * ( fracIndex - 2.0 ) * ( fracIndex - 3.0 );
		case 2:
			return 0.5 * fracIndex * ( 1.0 - fracIndex ) * ( fracIndex - 3.0 );
		case 3:
			return oneOverSix * fracIndex * ( fracIndex - 1.0 ) * ( fracIndex - 2.0 );
		default:
			cerr << "#error: no need to compute alpha_" << i << NEWL;
			exit( 2 );
	}
}
*/

Double TSoot::GetAlphaI2( int i, Double fracIndex, int first, int second )
{
	if ( i == first ) {
		return ( fracIndex - second ) / ( first - second );
	}
	else if ( i == second ) {
		return ( fracIndex - first ) / ( second - first );
	}
	else {
		return 0.0;
	}
}

Double TSoot::GetAlphaI2( int i, Double fracIndex )
{	
#ifdef NEWINTERPOLATION
	int	second;

	if ( fracIndex < 0.0 ) {
		second = 1;
	}
	else if ( fracIndex > fNSootMoments-1 ) {
		second = fNSootMoments-1;
	}
	else {
		second = ( int )fracIndex+1;
	}

	return GetAlphaI2( i, fracIndex, 0, second );
#else
	int	first;

	if ( fracIndex < 0.0 ) {
		first = 0;
	}
	else if ( fracIndex > fNSootMoments-1 ) {
		first = fNSootMoments-2;
	}
	else {
		first = ( int )fracIndex;
	}

	return GetAlphaI2( i, fracIndex, first, first+1 );
#endif
}

void TSoot::ComputeFractionalMomentsPAH( Double *moments )
{
	Double	*fracMom = fFracMomPAH->vec;

	fracMom[km1_2] = FracMom2( -0.5, moments );
	fracMom[k01_2] = FracMom2( 0.5, moments );
	fracMom[k03_2] = FracMom2( 1.5, moments );
	fracMom[k05_2] = FracMom2( 2.5, moments );

	fracMom[k19_6] = FracMom2( 19.0 / 6.0, moments );
	fracMom[k17_6] = FracMom2( 17.0 / 6.0, moments );
	fracMom[k13_6] = FracMom2( 13.0 / 6.0, moments );
	fracMom[k11_6] = FracMom2( 11.0 / 6.0, moments );

	fracMom[k07_6] = FracMom2( 7.0 / 6.0, moments );
	fracMom[k05_6] = FracMom2( 5.0 / 6.0, moments );
	fracMom[k01_6] = FracMom2( 1.0 / 6.0, moments );
	fracMom[km1_6] = FracMom2( -1.0 / 6.0, moments );
}

void TSoot::ComputeFractionalMoments( Double *moments )
{
	Double	*fracMom = fFracMom->vec;
	
	fracMom[km1_2] = FracMom2( -0.5, moments );
	fracMom[k01_2] = FracMom2( 0.5, moments );
	fracMom[k03_2] = FracMom2( 1.5, moments );
	fracMom[k05_2] = FracMom2( 2.5, moments );

	fracMom[k19_6] = FracMom2( 19.0 / 6.0, moments );
	fracMom[k17_6] = FracMom2( 17.0 / 6.0, moments );
	fracMom[k13_6] = FracMom2( 13.0 / 6.0, moments );
	fracMom[k11_6] = FracMom2( 11.0 / 6.0, moments );

	fracMom[k07_6] = FracMom2( 7.0 / 6.0, moments );
	fracMom[k05_6] = FracMom2( 5.0 / 6.0, moments );
	fracMom[k01_6] = FracMom2( 1.0 / 6.0, moments );
	fracMom[km1_6] = FracMom2( -1.0 / 6.0, moments );
	
	fracMom[km1_24] = FracMom2( -1.0 / 24.0, moments );
	fracMom[k23_24] = FracMom2( 23.0 / 24.0, moments );
	fracMom[k47_24] = FracMom2( 47.0 / 24.0, moments );
	fracMom[k71_24] = FracMom2( 71.0 / 24.0, moments );

	fracMom[km1_3] = FracMom2( -1.0 / 3.0, moments );
	fracMom[k02_3] = FracMom2( 2.0 / 3.0, moments );
	fracMom[k05_3] = FracMom2( 5.0 / 3.0, moments );
	fracMom[k08_3] = FracMom2( 8.0 / 3.0, moments );

}

Double TSoot::GetDerivFracMom( int l, Double which, Double *moments )
{
	// returns dM_which/dM_l
	return GetAlphaI2( l, which ) * FracMom2( which, moments ) / moments[l];
}
Double TSoot::GetDerivFracMom( int l, FracMoments nMom, Double *moments )
{
	// returns dM_r/dM_l

	Double	*fm = fFracMom->vec;

	switch( nMom ) {
		case km1_2:
			return GetAlphaI2( l, -1.0/2.0 ) * fm[nMom] / moments[l];
		case k01_2:
			return GetAlphaI2( l, 1.0/2.0 ) * fm[nMom] / moments[l];
		case k05_6:
			return GetAlphaI2( l, 5.0/6.0 ) * fm[nMom] / moments[l];
		case k01_6:
			return GetAlphaI2( l, 1.0/6.0 ) * fm[nMom] / moments[l];
		case km1_6:
			return GetAlphaI2( l, -1.0/6.0 ) * fm[nMom] / moments[l];

		case k03_2:
			return GetAlphaI2( l, 1.5 ) * fm[nMom] / moments[l];
		case k11_6:
			return GetAlphaI2( l, 11.0/6.0 ) * fm[nMom] / moments[l];
		case k07_6:
			return GetAlphaI2( l, 7.0/6.0 ) * fm[nMom] / moments[l];

		case k05_2:
			return GetAlphaI2( l, 2.5 ) * fm[nMom] / moments[l];
		case k19_6:
			return GetAlphaI2( l, 19.0/6.0 ) * fm[nMom] / moments[l];
		case k17_6:
			return GetAlphaI2( l, 17.0/6.0 ) * fm[nMom] / moments[l];
		case k13_6:
			return GetAlphaI2( l, 13.0/6.0) * fm[nMom] / moments[l];

		case km1_24:
			return GetAlphaI2( l, -1.0/24.0 ) * fm[nMom] / moments[l];
		case k23_24:
			return GetAlphaI2( l, 23.0/24.0 ) * fm[nMom] / moments[l];
		case k47_24:
			return GetAlphaI2( l, 47.0/24.0 ) * fm[nMom] / moments[l];
		case k71_24:
			return GetAlphaI2( l, 71.0/24.0 ) * fm[nMom] / moments[l];

		case km1_3:
			return GetAlphaI2( l, -1.0/3.0 ) * fm[nMom] / moments[l];
		case k02_3:
			return GetAlphaI2( l, 2.0/3.0 ) * fm[nMom] / moments[l];
		case k05_3:
			return GetAlphaI2( l, 5.0/3.0 ) * fm[nMom] / moments[l];
		case k08_3:
			return GetAlphaI2( l, 8.0/3.0 ) * fm[nMom] / moments[l];

		default:
			cerr << "#error: no need to compute derivative of alpha_" << (int)nMom << NEWL;
			exit( 2 );
	}
	return 0.0;
}

Double TSoot::GetDerivPhi( int l, Phi which, Double *moments )
{
	Double	*phi = fPhi->vec;
	Double	*fm = fFracMom->vec;

	switch( which ) {
		case kh00:
			return 0.5 * phi[kh00] 
					* ( GetDerivPhi( l, k000, moments ) / phi[k000]
						+ GetDerivPhi( l, k100, moments ) / phi[k100] );
		case kh11:
			return 0.5 * phi[kh11] 
					* ( GetDerivPhi( l, k011, moments ) / phi[k011]
						+ GetDerivPhi( l, k111, moments ) / phi[k111] );
		case kh12:
			return 0.5 * phi[kh12] 
					* ( GetDerivPhi( l, k012, moments ) / phi[k012]
						+ GetDerivPhi( l, k112, moments ) / phi[k112] );
		case k000:
			return    2.0 * ( fm[k01_6] * GetDerivFracMom( l, km1_2, moments ) 
							+ fm[km1_2] * GetDerivFracMom( l, k01_6, moments )
							+ 2.0 * fm[km1_6] * GetDerivFracMom( l, km1_6, moments ) );
		case k100:
			return    2.0 * ( fm[k07_6] * GetDerivFracMom( l, km1_2, moments ) 
							+ fm[km1_2] * GetDerivFracMom( l, k07_6, moments )
					+ 2.0 * ( fm[k05_6] * GetDerivFracMom( l, km1_6, moments )
							+ fm[km1_6] * GetDerivFracMom( l, k05_6, moments ) )
							+ fm[k01_2] * GetDerivFracMom( l, k01_6, moments )
							+ fm[k01_6] * GetDerivFracMom( l, k01_2, moments ) );
		case k011:
			return    2.0 * ( fm[k07_6] * GetDerivFracMom( l, k01_2, moments ) 
							+ fm[k01_2] * GetDerivFracMom( l, k07_6, moments )
					  + 2.0 * fm[k05_6] * GetDerivFracMom( l, k05_6, moments ) );
		case k111:
			return    2.0 * ( fm[k13_6] * GetDerivFracMom( l, k01_2, moments ) 
							+ fm[k01_2] * GetDerivFracMom( l, k13_6, moments )
					+ 2.0 * ( fm[k11_6] * GetDerivFracMom( l, k05_6, moments )
							+ fm[k05_6] * GetDerivFracMom( l, k11_6, moments ) )
							+ fm[k03_2] * GetDerivFracMom( l, k07_6, moments )
							+ fm[k07_6] * GetDerivFracMom( l, k03_2, moments ) );
		case k012:
			return            fm[k07_6] * GetDerivFracMom( l, k03_2, moments ) 
							+ fm[k03_2] * GetDerivFracMom( l, k07_6, moments )
					+ 2.0 * ( fm[k05_6] * GetDerivFracMom( l, k11_6, moments )
							+ fm[k11_6] * GetDerivFracMom( l, k05_6, moments ) )
							+ fm[k01_2] * GetDerivFracMom( l, k13_6, moments )
							+ fm[k13_6] * GetDerivFracMom( l, k01_2, moments );
		case k112:
			return 		      fm[k13_6] * GetDerivFracMom( l, k03_2, moments ) 
							+ fm[k03_2] * GetDerivFracMom( l, k13_6, moments )
			  + 2.0 * ( 2.0 * fm[k11_6] * GetDerivFracMom( l, k11_6, moments )
							+ fm[k05_6] * GetDerivFracMom( l, k17_6, moments )
							+ fm[k17_6] * GetDerivFracMom( l, k05_6, moments ) )

							+ fm[k03_2] * GetDerivFracMom( l, k13_6, moments )
							+ fm[k13_6] * GetDerivFracMom( l, k03_2, moments )
							+ fm[k07_6] * GetDerivFracMom( l, k05_2, moments )
							+ fm[k05_2] * GetDerivFracMom( l, k07_6, moments )
							+ fm[k01_2] * GetDerivFracMom( l, k19_6, moments )
							+ fm[k19_6] * GetDerivFracMom( l, k01_2, moments );
		default:
			cerr << "#error: no need to compute derivative of phi_" << (int)which << NEWL;
			exit( 2 );
	}
	return 0.0;
}

void TSoot::ComputePhi( Double temp )
{
	Double	c = GetC( temp );
	Double	*fm = fFracMom->vec;
	Double	*phi = fPhi->vec;

	phi[k000] = ( fm[k01_6] * fm[km1_2] + fm[km1_6] * fm[km1_6] ) * 2.0;
	phi[k100] = ( fm[k07_6] * fm[km1_2] + fm[k05_6] * fm[km1_6] * 2.0 
									+ fm[k01_2] * fm[k01_6] ) * 2.0;
	
	phi[k011] = ( fm[k07_6] * fm[k01_2] + fm[k05_6] * fm[k05_6] ) * 2.0;
	phi[k111] = ( fm[k13_6] * fm[k01_2] + fm[k11_6] * fm[k05_6] * 2.0 
									+ fm[k03_2] * fm[k07_6] ) * 2.0;
	
	phi[k012] =   fm[k07_6] * fm[k03_2] + fm[k05_6] * fm[k11_6] * 2.0 
									+ fm[k01_2] * fm[k13_6];
	phi[k112] =   fm[k13_6] * fm[k03_2] 
					+ ( fm[k11_6] * fm[k11_6] + fm[k05_6] * fm[k17_6] ) * 2.0
				 	+ fm[k03_2] * fm[k13_6] + fm[k07_6] * fm[k05_2] + fm[k01_2] * fm[k19_6];
	
	phi[kh00] = sqrt( phi[k000] * phi[k100] ) * c;
	phi[kh11] = sqrt( phi[k011] * phi[k111] ) * c;
	phi[kh12] = sqrt( phi[k012] * phi[k112] ) * c;
	
/*	cerr << "phi[kh00] = " << phi[kh00] << TAB
		 << "phi[kh11] = " << phi[kh11] << TAB
		 << "phi[kh12] = " << phi[kh12] << NEWL;*/
}

Double TSoot::GetPhi( int x, int y, Double temp, Double *moments1, Double *moments2 )
{
	Double	phi0, phi1;
	Double	c = GetC( temp );
	
	phi0 = 		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 );

	phi1 = 		FracMom2( 7.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( 5.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( 1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( 5.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 7.0 / 6.0 + y, moments2 );

	return sqrt( phi0 * phi1 ) * c;
}

Double TSoot::GetPhiPAH( int x, int y, Double temp, Double *moments1, Double *moments2 )
{
	Double	phi0, phi1;
	Double	c = GetCPAH( temp );
	
	phi0 = 		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 );

	phi1 = 		FracMom2( 7.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( 5.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( 1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( 5.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 7.0 / 6.0 + y, moments2 );

	return sqrt( phi0 * phi1 ) * c;
}

//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i// 
Double TSoot::GetPhiPAHFROMA4( Double x, Double temp, Double *moments )
{
	Double	phi0, phi1;
	Double	c = GetC( temp );
	
	phi0 =      pow( 3.0, 1.0 / 3.0 ) * FracMom2(             x, moments )  
		+ 2.0 * pow( 3.0,-1.0 / 3.0 ) * FracMom2( 1.0 / 3.0 + x, moments ) 
		+       1.0 / 3.0             * FracMom2( 2.0 / 3.0 + x, moments );

	phi1 = 	    pow( 3.0, 7.0 / 3.0 ) * FracMom2(             x, moments )  
		+ 2.0 * pow( 3.0, 5.0 / 3.0 ) * FracMom2( 1.0 / 3.0 + x, moments ) 
		+       3.0                   * FracMom2( 2.0 / 3.0 + x, moments ) 
		+       pow( 3.0, 1.0 / 3.0 ) * FracMom2( 1.0       + x, moments )		
		+ 2.0 * pow( 3.0,-1.0 / 3.0 ) * FracMom2( 4.0 / 3.0 + x, moments )
		+       1.0 / 3.0             * FracMom2( 5.0 / 3.0 + x, moments );	

	return sqrt( phi0 * phi1 ) * c;
}
//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//

Double TSoot::GetDerivPhiCond( int l, int x, int y, Double temp, Double *moments1, Double *moments2  )
{
	Double	phi0, phi1, phi, dPhi0dMl, dPhi1dMl;
	Double	c = GetC( temp );

	phi0 = 		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 );

	phi1 = 		FracMom2( 7.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( 5.0 / 6.0 + x, moments1 ) * FracMom2( -1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 2.0 + x, moments1 ) * FracMom2( 1.0 / 6.0 + y, moments2 )
		+		FracMom2( 1.0 / 6.0 + x, moments1 ) * FracMom2( 1.0 / 2.0 + y, moments2 )
		+ 2.0 * FracMom2( -1.0 / 6.0 + x, moments1 ) * FracMom2( 5.0 / 6.0 + y, moments2 )
		+		FracMom2( -1.0 / 2.0 + x, moments1 ) * FracMom2( 7.0 / 6.0 + y, moments2 );

	dPhi0dMl = 		GetDerivFracMom( l, 1.0 / 6.0 + x, moments1 ) 
					* FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * GetDerivFracMom( l, -1.0 / 6.0 + x, moments1 ) 
					* FracMom2( -1.0 / 6.0 + y, moments2 )
		+		GetDerivFracMom( l, -1.0 / 2.0 + x, moments1 ) 
					* FracMom2( 1.0 / 6.0 + y, moments2 );

	dPhi1dMl = 		GetDerivFracMom( l, 7.0 / 6.0 + x, moments1 ) 
					* FracMom2( -1.0 / 2.0 + y, moments2 )
		+ 2.0 * GetDerivFracMom( l, 5.0 / 6.0 + x, moments1 ) 
					* FracMom2( -1.0 / 6.0 + y, moments2 )
		+		GetDerivFracMom( l, 1.0 / 2.0 + x, moments1 ) 
					* FracMom2( 1.0 / 6.0 + y, moments2 )
		+		GetDerivFracMom( l, 1.0 / 6.0 + x, moments1 ) 
					* FracMom2( 1.0 / 2.0 + y, moments2 )
		+ 2.0 * GetDerivFracMom( l, -1.0 / 6.0 + x, moments1 ) 
					* FracMom2( 5.0 / 6.0 + y, moments2 )
		+		GetDerivFracMom( l, -1.0 / 2.0 + x, moments1 ) 
					* FracMom2( 7.0 / 6.0 + y, moments2 );

	phi = sqrt( phi0 * phi1 ) * c;
	
	if ( phi >= c * SMALLSOOT ) {
		return 0.5 * phi * ( dPhi0dMl / phi0 + dPhi1dMl / phi1 );
	}
	else {
		return 0.0;
	}
}

void TSoot::ComputePhiPAH( Double temp )
{
	Double	c = GetC( temp );
	Double	*fm = fFracMomPAH->vec;
	Double	*phi = fPhiPAH->vec;

	phi[k000] = ( fm[k01_6] * fm[km1_2] + fm[km1_6] * fm[km1_6] ) * 2.0;
	phi[k100] = ( fm[k07_6] * fm[km1_2] + fm[k05_6] * fm[km1_6] * 2.0 
									+ fm[k01_2] * fm[k01_6] ) * 2.0;
	
	phi[k011] = ( fm[k07_6] * fm[k01_2] + fm[k05_6] * fm[k05_6] ) * 2.0;
	phi[k111] = ( fm[k13_6] * fm[k01_2] + fm[k11_6] * fm[k05_6] * 2.0 
									+ fm[k03_2] * fm[k07_6] ) * 2.0;
	
	phi[k012] =   fm[k07_6] * fm[k03_2] + fm[k05_6] * fm[k11_6] * 2.0 
									+ fm[k01_2] * fm[k13_6];
	phi[k112] =   fm[k13_6] * fm[k03_2] 
					+ ( fm[k11_6] * fm[k11_6] + fm[k05_6] * fm[k17_6] ) * 2.0
				 	+ fm[k03_2] * fm[k13_6] + fm[k07_6] * fm[k05_2] + fm[k01_2] * fm[k19_6];
	
	phi[kh00] = sqrt( phi[k000] * phi[k100] ) * c;
	phi[kh11] = sqrt( phi[k011] * phi[k111] ) * c;
	phi[kh12] = sqrt( phi[k012] * phi[k112] ) * c;
	
/*	cerr << "phi[kh00] = " << phi[kh00] << TAB
		 << "phi[kh11] = " << phi[kh11] << TAB
		 << "phi[kh12] = " << phi[kh12] << NEWL;*/
}

void TSoot::CheckSolution( Double *theMom )
{
	for ( int i = 0; i < fNSootMoments; ++i ) {
		switch ( i ) {
			case 0:
				theMom[0] = MAX( theMom[0], SMALLSOOT );
				break;
			default:
				theMom[i] = MAX( theMom[i], 9 * theMom[i-1] );
				break;
		}
	}
}

Double TSoot::SourceCoagulation( int i )
{
	switch ( i ) {
		case 0:
			return -0.5 * fPhi->vec[kh00];
		case 1:
			return 0.0;
		case 2:
			return fPhi->vec[kh11];
		case 3:
			return 3.0 * fPhi->vec[kh12];
		default:
			cerr << "#error: no need to compute coagulation source_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

Double TSoot::SourceSurfDepCoag( int i, Double *mom, Double *Y, Double temp
					, Double density, Double *molarMass )
{
	Double	ksg, c, coag, S;
	Double	Cst = 1.0 / ( 6.0 * fCoagFact );

	switch ( i ) {
		case 0:
			ksg = GetSurfGrowthCoeff( Y, density, molarMass );
			c = GetC( temp );
			coag = -0.5 * fPhi->vec[kh00];
			S = MAX( 1.0e-30, GetSCoag( mom, temp ) );
//			fprintf( stderr, "S = %g\tCoag = %g\tksg = %g\tc = %g\tCst = %g\n", S, coag, ksg, c, Cst );

			return coag * ( 1.0 - 1.0 / ( 1.0 + Cst * fPhi->vec[kh00] * ksg / ( c * c * S ) ) );
		case 1:
			return 0.0;
		case 2:
			fprintf( stderr, "###error: SourceSurfDepCoag for higher moments not yet implemented\n" );
			return 0.0;
		case 3:
			fprintf( stderr, "###error: SourceSurfDepCoag for higher moments not yet implemented\n" );
			return 0.0;
		default:
			cerr << "#error: no need to compute coagulation source_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

Double TSoot::GetSCoag( Double *mom, Double /*temp*/ )
{
	Double	S00, S01, S10, S11;
	Double	M_1_6 = FracMom2( 1.0 / 6.0, mom );
	Double	M_m1_6 = FracMom2( -1.0 / 6.0, mom );
	Double	M_1_2 = FracMom2( 1.0 / 2.0, mom );
	Double	M_m1_2 = FracMom2( -1.0 / 2.0, mom );
	Double	M_3_2 = FracMom2( 3.0 / 2.0, mom );
	Double	M_5_2 = FracMom2( 5.0 / 2.0, mom );
	Double	M_5_6 = FracMom2( 5.0 / 6.0, mom );
	Double	M_11_6 = FracMom2( 11.0 / 6.0, mom );
	Double	M_17_6 = FracMom2( 17.0 / 6.0, mom );
	
	S00 = M_m1_6 * M_m1_6 * ( M_m1_2 + 2.0 * M_1_2 + M_3_2 );
	S01 = 2.0 * ( M_11_6 * M_m1_6 * M_m1_2 
					+ M_5_6 * M_5_6 * M_m1_2 
					+ 2.0 * M_5_6 * M_m1_6 * M_1_2 ) 
			+ M_m1_6 * M_m1_6 * M_3_2;
	S10 = 2.0 * M_5_6 * M_m1_6 * M_m1_2
			+ M_m1_6 * M_m1_6 * M_1_2
			+ 4.0 * M_5_6 * M_m1_6 * M_1_2
			+ 2.0 * M_m1_6 * M_m1_6 * M_3_2
			+ 2.0 * M_5_6 * M_m1_6 * M_3_2
			+ M_m1_6 * M_m1_6 * M_5_2;
	S11 = 2.0 * M_17_6 * M_m1_6 * M_m1_2
			+ M_m1_6 * M_m1_6 * M_5_2
			+ 6.0 * M_11_6 * M_5_6 * M_m1_2
			+ 6.0 * M_11_6 * M_m1_6 * M_1_2
			+ 6.0 * M_5_6 * M_m1_6 * M_3_2
			+ 6.0 * M_5_6 * M_5_6 * M_1_2;


	return pow( S00 * S10, 1.0/3.0 ) * pow( S01 * S11, 1.0/6.0 );
}

Double TSoot::SourceCoagulationNew( int i, Double temp, Double *moments )
{
	switch ( i ) {
		case 0:
//			return ( moments[0] > 1.0e-30 ) ? -0.5 * GetPhi( 0, 0, temp, moments ) : 0.0;
			return -/* 1.0 / ( 1.0 + 1.0e5 * SMALLSOOT/ MAX( moments[0], 1.0e-60 ) ) * */ 0.5 * GetPhi( 0, 0, temp, moments );
		case 1:
			return 0.0;
		case 2:
			return GetPhi( 1, 1, temp, moments );
		case 3:
			return 3.0 * GetPhi( 1, 2, temp, moments );
		default:
			cerr << "#error: no need to compute coagulation source_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

Double TSoot::SourceCondensation( int i, Double temp, Double *pahMoments, Double *moments )
{
	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
			return GetBetaCond( temp, moments ) * moments[0] * pahMoments[1];
		case 2:
			return GetBetaCond( temp, moments ) * ( moments[0] * pahMoments[2] 
								+ 2.0 * moments[1] * pahMoments[1] );
		case 3:
			return GetBetaCond( temp, moments ) * ( moments[0] * pahMoments[3] 
								+ 3.0 * ( moments[1] * pahMoments[2] 
										+ moments[2] * pahMoments[1] ) );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

#ifdef NEWPOLYNUCCOND
Double TSoot::SourceCondensationNew( int i, Double temp, Double *pahMoments, Double *moments, Double */*Y*/
                                   , Double /*density*/, Double */*molarMass*/  )
{

	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
			return GetBetaCond( temp, moments ) * moments[0] * pahMoments[1];
		case 2:
			return GetBetaCond( temp, moments ) * ( moments[0] * pahMoments[2] 
								+ 2.0 * moments[1] * pahMoments[1] );
		case 3:
			return GetBetaCond( temp, moments ) * ( moments[0] * pahMoments[3] 
								+ 3.0 * ( moments[1] * pahMoments[2] 
										+ moments[2] * pahMoments[1] ) );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
#else
Double TSoot::SourceCondensationNew( int i, Double temp, Double *pahMom, Double *mom, Double *Y
                                   , Double density, Double *molarMass  )
{	
	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
#	ifdef PAHFROMA4
// changed by hp from
//			return 9.0 * pahMom[0] * GetPhiPAHFROMA4( x, temp, pahMom );
// to
			return 9.0 * pahMom[0] * GetPhiPAHFROMA4( -0.5, temp, mom );
#	else
			return GetPhi( 0, 1, temp, mom, pahMom );
#	endif
		case 2:
			return GetPhi( 0, 2, temp, mom, pahMom ) + 2.0 * GetPhi( 1, 1, temp, mom, pahMom );
		case 3:
			return GetPhi( 0, 3, temp, mom, pahMom ) 
				+ 3.0 * ( GetPhi( 1, 2, temp, mom, pahMom ) + GetPhi( 2, 1, temp, mom, pahMom ) );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
#endif
}

void TSoot::ComputeCSootStar( Double *k, Double *Y, Double density, Double *molarMass
							, Double mixMolarMass )
{
	Double	C_H = density * Y[f_H] / molarMass[f_H];
	Double	C_OH = density * Y[f_OH] / molarMass[f_OH];
	Double	C_H2 = density * Y[f_H2] / molarMass[f_H2];
	Double	C_H2O = density * Y[f_H2O] / molarMass[f_H2O];
	Double	C_O2 = density * Y[f_O2] / molarMass[f_O2];
	Double	C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2];

#ifdef NEWSURFGROWTH
#	ifdef WITHTHIRDBODY
	fThirdBodyConc = density / mixMolarMass;
#	else
	fThirdBodyConc = 1.0;
#	endif

	fFk10 = k[ks13f] / MAX( 1.0e-60, k[ks13f] + k[ks10b] * fThirdBodyConc + k[ks112] * C_O2 );
	
	fCSootStar = ( k[ks7f] * C_OH + k[ks8f] * C_H + k[ks9f] 
				+ k[ks13b] * C_H * ( 1.0 - fFk10 ) + k[ks12] * C_OH ) 
				/ MAX( 1.0e-60, k[ks7b] * ( C_H2O + 1.0e-10 ) + k[ks8b] * ( C_H2 + 1.0e-10 ) + k[ks9b] * C_H + 
						k[ks10f] * fThirdBodyConc * C_C2H2 * fFk10 );

/* 	if ( fCSootStar > 1000.0 ) { */
/* 		return;		 */
/* 	} */
/* 	if ( fCSootStar > 1.0 ) { */
/* 		return;		 */
/* 	} */
	
	fCSootStar = MIN( 0.1, fCSootStar );
#else
	fCSootStar = C_H * k[ks8f] / CatchZero( C_H2 * k[ks8b] + C_H * k[ks9] + 
						C_C2H2 * k[ks10f] + C_O2 * k[ks11] );
#endif
}

void TSoot::ComputeSootReactionRates( Double *w, Double *Y, Double *moments, Double density
						, Double *molarMass, Double /*mixMolarMass*/ )
{
	Double	C_Soot = fAlpha * FracMom2( 2.0/3.0, moments );		//  CSoot = alpha * M_2/3
	Double	C_SootStar = fCSootStar * C_Soot;
	Double	C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2];
	Double	C_H = density * Y[f_H] / molarMass[f_H];
	Double	C_O2 = density * Y[f_O2] / molarMass[f_O2];
	Double	*k = fSootRateCoeffs->vec;
	Double	C_H2 = density * Y[f_H2] / molarMass[f_H2];
	Double	C_OH = density * Y[f_OH] / molarMass[f_OH];

#ifdef NEWSURFGROWTH

	Double	C_H2O = density * Y[f_H2O] / molarMass[f_H2O];
	Double	CSootC2H2 = ( k[ks10f] * fThirdBodyConc * C_C2H2 * C_SootStar 
						+ k[ks13b] * C_Soot * C_H )
				/ CatchZero( k[ks10b] * fThirdBodyConc + k[ks13f] + k[ks112] * C_O2 );

	w[ks7f] = k[ks7f] * C_Soot * C_OH;
	w[ks7b] = k[ks7b] * C_SootStar * C_H2O;
	w[ks8f] = k[ks8f] * C_Soot * C_H;
	w[ks8b] = k[ks8b] * C_SootStar * C_H2;
	w[ks9f] = k[ks9f] * C_Soot;
	w[ks9b] = k[ks9b] * C_SootStar * C_H;
	w[ks10f] = k[ks10f] * fThirdBodyConc * C_SootStar * C_C2H2;
	w[ks10b] = k[ks10b] * fThirdBodyConc * CSootC2H2;
	w[ks111] = k[ks111] * C_SootStar * C_O2;
	w[ks112] = k[ks112] * CSootC2H2 * C_O2;
	w[ks12] = k[ks12] * C_Soot * C_OH;
	w[ks13f] = k[ks13f] * CSootC2H2;
	w[ks13b] = k[ks13b] * C_Soot * C_H;
#else
	w[ks8f] = k[ks8f] * C_Soot * C_H;
	w[ks8b] = k[ks8b] * C_SootStar * C_H2;
	w[ks9] = k[ks9] * C_SootStar * C_H;
	w[ks10f] = k[ks10f] * C_SootStar * C_C2H2;
	w[ks10b] = k[ks10b] * C_Soot * C_H;
	w[ks11] = k[ks11] * C_SootStar * C_O2;
	w[ks12] = k[ks12] * C_Soot * C_OH;
#endif
}

Double TSoot::GetSurfGrowthCoeff( Double *Y, Double density, Double *molarMass )
{
	Double			C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2];
#ifdef NEWSURFGROWTH
	Double			C_H = density * Y[f_H] / molarMass[f_H];
	Double			*k = fSootRateCoeffs->vec;
	Double			locSG = k[ks10f] * fThirdBodyConc * C_C2H2 * fFk10 * fCSootStar 
							- k[ks13b] * C_H * ( 1.0 - fFk10 );

	return fAlpha * locSG;
#else
	return fAlpha * fCSootStar * fSootRateCoeffs->vec[ks10f] * C_C2H2;
#endif
}

Double TSoot::GetSurfGrowthCoeffForw( Double *Y, Double density, Double *molarMass )
{
	Double			C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2];
	Double			C_H = density * Y[f_H] / molarMass[f_H];
#ifdef NEWSURFGROWTH
	//Double			locSG = fSootRateCoeffs->vec[ks10f] * fThirdBodyConc * C_C2H2
	//				* fFk10 * fCSootStar;  changed by GB
	Double			locSG = (fSootRateCoeffs->vec[ks10f] * fThirdBodyConc * C_C2H2 
					* fCSootStar + fSootRateCoeffs->vec[ks13b] * C_H ) * fFk10;

	return fAlpha * locSG;
#else
	return fAlpha * fCSootStar * fSootRateCoeffs->vec[ks10f] * C_C2H2;
#endif
}

Double TSoot::GetSurfGrowthCoeffBackw( Double *Y, Double density, Double *molarMass )
{
#ifdef NEWSURFGROWTH
	Double			C_H = density * Y[f_H] / molarMass[f_H];
	//Double			locSG = fSootRateCoeffs->vec[ks13b] * C_H * ( 1.0 - fFk10 ); changed by GB
	Double			locSG = fSootRateCoeffs->vec[ks13b] * C_H;

	return fAlpha * locSG;
#else
	return 0.0;
#endif
}

Double TSoot::SourceSurfGrowth( int i, Double *Y, Double density, Double *molarMass )
{
	Double	coeff = GetSurfGrowthCoeff( Y, density, molarMass );
	Double	*fracMom = fFracMom->vec;

	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
			return coeff * fracMom[k23_24];
		case 2:
			return coeff * ( fracMom[k23_24] + 2.0 * fracMom[k47_24] );
		case 3:
			return coeff * ( fracMom[k23_24] + 3.0 * ( fracMom[k47_24] + fracMom[k71_24] ) );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );

	}
	return 0.0;
}

Double TSoot::SourceSurfGrowthNew( int i, Double *moments, Double *Y, Double density, Double *molarMass )
{
	Double	fmom23 = FracMom2( 2.0 / 3.0, moments );

	Double	M0Fact = 1.0;// / ( 1.0 + 1.0e5 * SMALLSOOT/ MAX( moments[0], 1.0e-60 ) );
	Double	M1Fact = 1.0;// / ( 1.0 + 1.0e5 * SMALLSOOT/ MAX( moments[1], 1.0e-60 ) );

//	M0Fact = ( moments[0] < SMALLSOOT ) ? 0.0 : M0Fact;
//	M1Fact = ( moments[1] < SMALLSOOT ) ? 0.0 : M1Fact;

	Double	coeffF = GetSurfGrowthCoeffForw( Y, density, molarMass );
	Double	coeffB = M0Fact * M1Fact * GetSurfGrowthCoeffBackw( Y, density, molarMass ); //changed by GB
#ifdef OXMOM0NEW
	Double	beta = MIN( M0Fact * moments[0] / MAX( moments[1], 9.0 * SMALLSOOT ), 1.0 / 9.0 );//changed by GB
//	Double	beta = MIN( M0Fact * moments[0] / MAX( moments[1], 18.0 * SMALLSOOT ), 1.0 / 18.0 );
//	Double	beta = MIN( M0Fact * moments[0] / MAX( moments[1], 18.0 * SMALLSOOT ), 1.0 );
#else
	Double	beta = fBeta;
#endif
//	Double	C_H = density * Y[f_H] / molarMass[f_H];


	switch( i ) {
		case 0: 
//			return 0.0;
		        return - beta * coeffB * FracMom2( -1.0 / 3.0, moments );
		case 1:
			return ATTENTION * ( coeffF - coeffB ) * FracMom2( 2.0 / 3.0, moments);
		case 2:
			return ( coeffF - coeffB ) * ( fmom23 + 2.0 * FracMom2( 5.0 / 3.0, moments ) );
		case 3:
			return ( coeffF - coeffB ) * ( fmom23 + 3.0 * ( FracMom2( 5.0 / 3.0, moments ) 
						+ FracMom2( 8.0 / 3.0, moments ) ) );
		default:
			cerr << "#error: no need to compute source term for surface growth for M_" << i << NEWL;
			exit( 2 );

	}
	return 0.0;
}

Double TSoot::SourceSurfGrowthNew( int i, Double *Y, Double density, Double *molarMass )
{
	Double	coeff = GetSurfGrowthCoeff( Y, density, molarMass );
	Double	*fracMom = fFracMom->vec;

	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
			return coeff * fracMom[k02_3];
		case 2:
			return coeff * ( fracMom[k02_3] + 2.0 * fracMom[k05_3] );
		case 3:
			return coeff * ( fracMom[k02_3] + 3.0 * ( fracMom[k05_3] + fracMom[k08_3] ) );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );

	}
	return 0.0;
}

Double TSoot::GetSootRadiation( Double temp, Double *moments )
{
	const Double	vol1 = 2.46e-10 * 2.46e-10 * 3.51e-10;		// [m^3]
	Double			fvo = vol1 * AVOGADRO * moments[1];

	// alphas taken from:
	//		Hubbard, G. L. , Tien C. L: 
	//		Infrared Mean Absorption Coefficient
	//		of Luminous Flames and Smoke
	//		Journal of Heat Transfer, vol 100, p. 235ff, 1978
	Double			alphas = -3.75e5 + 1735.0 * temp;			// [m^-1] 
	
	return MAX( MAGIC * 4.0 * alphas * fvo * STEFBOLTZ * pow( temp, 4.0 ), 0.0 );
}

Double TSoot::GetSootRadRossCoeff( Double temp, Double *moments, Double cutoff )
{
	const Double	vol1 = 2.46e-10 * 2.46e-10 * 3.51e-10;		// [m^3]
	Double			fvo = vol1 * AVOGADRO * moments[1];

	// alphas taken from:
	//		Hubbard, G. L. , Tien C. L: 
	//		Infrared Mean Absorption Coefficient
	//		of Luminous Flames and Smoke
	//		Journal of Heat Transfer, vol 100, p. 235ff, 1978
	// Rosseland model is del q_R = del(-16/3 sigma T^3/alpha del(T))
	// equal to del(coeff del(T)) = coeff del^2 T + del(coeff) del(T)
	// The coefficient returned here is     coeff = -16/3 sigma T^3/alpha

	Double			alphas = -3.75e5 + 1735.0 * MAX( 273.0, temp );			// [m^-1] 
	
	cutoff = MAX( 1.0e-30, cutoff );

	return -MAX( MAGIC * 16.0/3.0 * STEFBOLTZ * pow( temp, 3.0 ) / ( alphas * MAX( cutoff, fvo ) ), 0.0 );
}

Double TSoot::SourceSootOxidationNew( int i, Double *moments, Double *Y, Double density, Double *molarMass )
{
	Double	coeff = GetSootOxCoeff( Y, density, molarMass );	// F4
//	Double	beta;
	
/*	if ( moments[1] < SMALLSOOT || moments[0] < SMALLSOOT ) {*/
/*		coeff = 0.0;*/
/*	}*/
	
	Double	M0Fact = 1.0;// / ( 1.0 + 1.0e5 * SMALLSOOT/ MAX( moments[0], 1.0e-60 ) );
	Double	M1Fact = 1.0;// / ( 1.0 + 1.0e5 * SMALLSOOT/ MAX( moments[1], 1.0e-60 ) );

//	M0Fact = ( moments[0] < SMALLSOOT ) ? 0.0 : M0Fact;
//	M1Fact = ( moments[1] < SMALLSOOT ) ? 0.0 : M1Fact;

#ifdef OXMOM0NEW
	//beta = MIN( M0Fact * moments[0] / MAX( moments[1], 9.0 * SMALLSOOT ), 1.0 / 9.0 );changed by GB
//	beta = MIN( M0Fact * moments[0] / MAX( moments[1], 18.0 * SMALLSOOT ), 1.0 / 18.0 );
#else
	beta = fBeta;
#endif


	switch( i ) {
		case 0: 
		  //return -M0Fact * M1Fact * beta * coeff * FracMom2( -1.0 / 3.0, moments ); changed by GB
			return -coeff * FracMom2( -1.0 / 3.0, moments );
		case 1:
		  //return -M0Fact * M1Fact * coeff * FracMom2( 2.0 / 3.0, moments ); changed by GB
			return - coeff * FracMom2( 2.0 / 3.0, moments );
		case 2:
			return -coeff * ( -FracMom2( 2.0 / 3.0, moments ) 
								+ 2.0 * FracMom2( 5.0 / 3.0, moments ) );
		case 3:
			return 0.0/*-coeff * ( FracMom2( 2.0 / 3.0, moments ) 
								+ 3.0 * ( -FracMom2( 5.0 / 3.0, moments ) 
								+ FracMom2( 8.0 / 3.0, moments ) ) )*/;
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

Double TSoot::SourceSootOxidationNew( int i, Double *Y, Double density, Double *molarMass )
{
	Double	coeff = GetSootOxCoeff( Y, density, molarMass );	// F4
	Double	*fracMom = fFracMom->vec;

	switch( i ) {
		case 0: 
#ifdef OXMOM0NEW
			fprintf( stderr, "#warning: not yet implemented\n" );
			return -fracMom[km1_3] / fracMom[k02_3] * coeff * fracMom[k02_3];
#else
			return -fBeta * coeff * fracMom[km1_3];
#endif

		case 1:
			return -coeff * fracMom[k02_3];
		case 2:
			return -coeff * ( ( 1.0 - fBeta ) * fracMom[k02_3] + 2.0 * fracMom[k05_3] );
		case 3:
			return -coeff * ( ( 1.0 - fBeta ) * ( fracMom[k02_3] + 3.0 * fracMom[k05_3] ) 
								+ 3.0 * fracMom[k08_3] );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

/*Double TSoot::SourceSootOxidation( int i, Double *Y, Double density, Double *molarMass )
{
	Double	coeff = GetSootOxCoeff( Y, density, molarMass );	// F4
	Double	*fracMom = fFracMom->vec;

	switch( i ) {
		case 0: 
			return -fBeta * coeff * fracMom[km1_24];
		case 1:
			return -coeff * fracMom[k23_24];
		case 2:
			return -coeff * ( ( 1.0 - fBeta ) * fracMom[k23_24] + 2.0 * fracMom[k47_24] );
		case 3:
			return -coeff * ( ( 1.0 - fBeta ) * ( fracMom[k23_24] + 3.0 * fracMom[k47_24] ) 
								+ 3.0 * fracMom[k71_24] );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );

	}
}
*/

Double TSoot::GetSootOxCoeff( Double *Y, Double density, Double *molarMass )
{
	Double			C_OH = density * Y[f_OH] / molarMass[f_OH];
	Double			C_O2 = density * Y[f_O2] / molarMass[f_O2];
	Double                  C_H = density * Y[f_H] / molarMass[f_H]; //added by GB
	Double                  C_C2H2 = density * Y[f_C2H2] / molarMass[f_C2H2]; //added by GB
	Double                  fCSootC2H2;
	Double			*k = fSootRateCoeffs->vec;

#ifdef NEWSURFGROWTH

	fCSootC2H2 = ( k[ks10f] * C_C2H2 * fThirdBodyConc * fCSootStar + k[ks13b] * C_H ) * fFk10 / k[ks13f]; //added by GB
	//return fAlpha * ( k[ks111] * C_O2 * fCSootStar + k[ks12] * C_OH ); changed by GB
	return fAlpha * ( ( k[ks111] * fCSootStar + k[ks112] * fCSootC2H2 ) * C_O2 + k[ks12] * C_OH );

#else
	return fAlpha * ( k[ks10b] * ( density * Y[f_H] / molarMass[f_H] ) 
				+ k[ks12] * C_OH + k[ks11] * C_O2 * fCSootStar );
#endif
}

void TSoot::UpdateSoot( TReactionPtr reaction, TSpeciesPtr species
			, Double *moments, Double temp, Double *Y, Double density, Double mixMolarMass )
{
	Double	*kSoot = fSootRateCoeffs->vec;
	Double	*molarMass = species->GetMolarMass()->vec;

	ComputeFractionalMoments( moments );
	ComputePhi( temp );
	
	ComputeSootRateCoeffs( kSoot, temp, reaction );
	ComputeCSootStar( kSoot, Y, density, molarMass, mixMolarMass );
}

void TSoot::MOverRhoToM( Double *momSource, Double *momDest, int nMom
				, Double *Y, Double temp, Double pressure
				, Double *molarmass, int nSpeciesIn, TPropertiesPtr props )
{
	int		i;
	Double	mixMolarMass;
	Double	rho;
	
	props->ComputeMixtureMolarMass( mixMolarMass, Y, molarmass, nSpeciesIn );
	rho = pressure * mixMolarMass / ( RGAS * temp );
	
	for ( i = 0; i < nMom; ++i ) {
		momDest[i] = momSource[i] * rho;
	}
}

void TSoot::MOverRhoToM( Double *momSource, Double *momDest, int nMom, Double rho )
{
	int		i;
	for ( i = 0; i < nMom; ++i ) {
		momDest[i] = momSource[i] * rho;
	}
}

#ifndef ZEROD
void T1DSoot::InitT1DSoot( TInputDataPtr input )
{
	int			nOfGridPoints = input->fMaxGridPoints;
	int 		nPAHReactions = input->GetCounter()->pahReactions;
	int 		nSpeciesIn = input->GetCounter()->species - input->GetCounter()->steadyStates;
	fReactionRate = NewMatrix( nPAHReactions, nOfGridPoints, kColumnPointers );

	fPij = NewTensor( nOfGridPoints+2, fNPAHMolecules+1, fNStages, kColumnPointers );
	fPij->tensor = &fPij->tensor[kNext];
	fSumPi = NewMatrix( fNPAHMolecules, nOfGridPoints+2, kColumnPointers );
	fSumPi->mat = &fSumPi->mat[kNext];
	fPAHMoments = NewMatrix( fNPAHMoments, nOfGridPoints+2, kColumnPointers );
	fPAHMoments->mat = &fPAHMoments->mat[kNext];

	fMoments = NewMatrix( fNSootMoments, nOfGridPoints+2, kColumnPointers );
	fMoments->mat = &fMoments->mat[kNext];
	fSavedMoments = NewMatrix( fNSootMoments, nOfGridPoints+2, kColumnPointers );
	fSavedMoments->mat = &fSavedMoments->mat[kNext];

	fSootDiff = NewVector( nOfGridPoints+2 );
	fSootDiff->vec = &fSootDiff->vec[kNext];

    fDMdx = NewTensor( nOfGridPoints, fNSootMoments, nSpeciesIn+fNSootMoments+1, kColumnPointers );
	fSource = NewMatrix( fNSootMoments, nOfGridPoints, kColumnPointers );
    fSourceMod = NewVector( fNSootMoments );
}

T1DSoot::~T1DSoot( void )
{
	DisposeVector( fSourceMod );
	DisposeMatrix( fSource );
	DisposeTensor( fDMdx );

	fSootDiff->vec = &fSootDiff->vec[kPrev];
	DisposeVector( fSootDiff );


	fSavedMoments->mat = &fSavedMoments->mat[kPrev];
	DisposeMatrix( fSavedMoments );
	fMoments->mat = &fMoments->mat[kPrev];
	DisposeMatrix( fMoments );

	fPAHMoments->mat = &fPAHMoments->mat[kPrev];
	DisposeMatrix( fPAHMoments );
	fSumPi->mat = &fSumPi->mat[kPrev];
	DisposeMatrix( fSumPi );
	fPij->tensor = &fPij->tensor[kPrev];
	DisposeTensor( fPij );

	DisposeMatrix( fReactionRate );
}

void T1DSoot::UpdateDimensions( int len )
{
	fMoments->cols = len;
}

void T1DSoot::UpdateSolution( T1DFlamePtr flame, Double *y, int gridPoint )
{
	int		i;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		firstSpeciesOff = flame->GetOffsetFirstSpecies();
	Double	*moments = fMoments->mat[gridPoint];
	Double	mixMM, density;
		
	flame->GetProperties()->ComputeMixtureMolarMass( mixMM, &y[firstSpeciesOff]
			, flame->GetSpecies()->GetMolarMass()->vec, nSpeciesInSystem );
	density = flame->GetPressure() * mixMM / ( RGAS * y[flame->GetOffsetTemperature()] );

	for ( i = 0; i < fNSootMoments; ++i ) {
		moments[i] = y[i+fOffsetMoments] * density;
	}
}

void T1DSoot::UpdateSolution( T1DFlamePtr flame, MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec )
{
	T1DPropertiesPtr	prop = flame->GetProperties();
	int		i, k, nGridPoints = yMat->cols;
	int		nSpeciesInSystem = flame->GetSpecies()->GetNSpeciesInSystem();
	int		firstSpeciesOff = flame->GetOffsetFirstSpecies();
	int		tempOff = flame->GetOffsetTemperature();
	Double	**y = yMat->mat;
	Double	*yLeft = yLeftVec->vec;
	Double	*yRight = yRightVec->vec;
	Double	**moments = fMoments->mat;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	mixMM, density;
	Double	pressure = flame->GetPressure();

//	set boundary values
//	left side
	prop->ComputeMixtureMolarMass( mixMM, &yLeft[firstSpeciesOff]
			, molarMass, nSpeciesInSystem );
	density = pressure * mixMM / ( RGAS * yLeft[tempOff] );
	for ( i = 0; i < fNSootMoments; ++i ) {
		moments[kPrev][i] = yLeft[i+fOffsetMoments] * density;
	}

//	right side
	prop->ComputeMixtureMolarMass( mixMM, &yRight[firstSpeciesOff]
			, molarMass, nSpeciesInSystem );
	density = pressure * mixMM / ( RGAS * yRight[tempOff] );
	for ( i = 0; i < fNSootMoments; ++i ) {
		moments[nGridPoints][i] = yRight[i+fOffsetMoments] * density;
	}

//	set inner values
	for ( k = 0; k < nGridPoints; ++k ) {
		prop->ComputeMixtureMolarMass( mixMM, &y[k][firstSpeciesOff]
				, molarMass, nSpeciesInSystem );
		density = pressure * mixMM / ( RGAS * y[k][tempOff] );
		for ( i = 0; i < fNSootMoments; ++i ) {
			moments[k][i] = y[k][i+fOffsetMoments] * density;
		}
	}
}

void T1DSoot::SolutionToSolver( Double **y )
{
	int		i, k, nGridPoints = fMoments->cols;
	Double	**moments = fMoments->mat;

	fprintf( stderr, "check if 'T1DSoot::SolutionToSolver' is correct\n" );

	for ( k = 0; k < nGridPoints; ++k ) {
		for ( i = 0; i < fNSootMoments; ++i ) {
			y[k][i+fOffsetMoments] = moments[k][i];
		}
	}
}

void T1DSoot::SaveSolution( void )
{
	int		i, k;
	int		len = fMoments->cols;
	Double	**moments = fMoments->mat;
	Double	**savedMoments = fSavedMoments->mat;

//	fprintf( stderr, "check if 'T1DSoot::SaveSolution' is correct\n" );

	fSavedMoments->cols = fMoments->cols;
	fSavedMoments->rows = fMoments->rows;

	for ( k = -1; k <= len; ++k ) {
		for ( i = 0; i < fNSootMoments; ++i ) {
			savedMoments[k][i] = moments[k][i];
		}
	}
}

void T1DSoot::RestoreSolution( void )
{
	int		i, k;
	int		len = fSavedMoments->cols;
	Double	**moments = fMoments->mat;
	Double	**savedMoments = fSavedMoments->mat;

	fprintf( stderr, "check if 'T1DSoot::RestoreSolution' is correct\n" );

	for ( k = -1; k <= len; ++k ) {
		for ( i = 0; i < fNSootMoments; ++i ) {
			moments[k][i] = savedMoments[k][i];
		}
	}
}

void T1DSoot::FillJacobi( T1DFlamePtr flame, NodeInfoPtr nodeInfo, CoordType coordinate )
{
	int				i, ioff, l;
	int				tempOff = flame->GetOffsetTemperature();
	int				specOff = flame->GetOffsetFirstSpecies();
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			*pahMoments = flameNode->pahMoments;
	Double			**a = nodeInfo->a;
	Double			hnenn = nodeInfo->hnenn;
	Double			temp = flameNode->temp[kCurr];
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double			density = flameNode->mixDensity[kCurr];
	Double			*moments = flameNode->moments;
	Double			*Y = flameNode->Y[kCurr];
	Double			*kSoot = fSootRateCoeffs->vec;
	
	ComputeFractionalMoments( moments );
	ComputePhi( temp );

	ComputeSootRateCoeffs( kSoot, temp, flame->GetReaction() );
	ComputeCSootStar( kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr] );

	if ( coordinate == kPhysical ) {
		int fVVelocity = flame->GetOffsetVVelocity();
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i+fOffsetMoments;

			//	convection
			FillJacNonlinearConvectUpwind( fVVelocity, ioff, nodeInfo, 1.0 );

			//	diffusion
			if ( fSizeDepDiff ) {
#ifndef NODIFF
				FillJacSootDiffusion( i, -1.0, kPhysical, flame, nodeInfo );
#endif
			}
			else {
				FillJacSootDiffusionNew( i, -1.0, kPhysical, flame, nodeInfo );
			}

			// thermophoresis
			if ( fThermoPhoresis ) {
				FillJacSootThermoPhoresis( i, -1.0, kPhysical, flame, nodeInfo );
			}

			// source term m_i/rho
			for ( l = 0; l < fNSootMoments; ++l ) {
				a[l+fOffsetMoments][ioff] -= flameNode->dMdx[l][i] * hnenn;
			}
			a[tempOff][ioff] -= flameNode->dMdx[fNSootMoments][i] * hnenn;
			int	nSpeciesIn = flame->GetSpecies()->GetNSpeciesInSystem();
			for ( l = 0; l < nSpeciesIn; ++l ) {
				a[l+specOff][ioff] -= flameNode->dMdx[l+1+fNSootMoments][i] * hnenn;
			}

/*			// nucleation
			if ( fNucleation ) {
				a[tempOff][ioff] -= 0.5 / temp * NucleationNew( i, temp, pahMoments )
								 / density * hnenn;
			}

			// coagulation
			if ( fCoagulation ) {
				// 1.0 / density already in 'SootCoagFunc' (since y is used for moments)
				a[tempOff][ioff] -= dfdyUpwind( tempOff, ioff, SootCoagFunc, nodeInfo, flame )
											* hnenn;
//				a[ioff][ioff] -= FillJacSourceCoagulation( i, i, moments )
//											* hnenn;
//				a[ioff][ioff] -= dfdyUpwind( ioff, ioff, SootCoagFunc, nodeInfo, flame )
//											* hnenn;
				for ( l = 0; l < fNSootMoments; ++l ) {
					a[l+fOffsetMoments][ioff] -= dfdyUpwind( l+fOffsetMoments, ioff, SootCoagFunc, nodeInfo, flame )
												* hnenn;
//					a[l+fOffsetMoments][ioff] -= FillJacSourceCoagulation( i, l, moments )
//												* hnenn;
				}
			}

			//	condensation
			if ( fCondensation ) {
				a[tempOff][ioff] -= SourceCondensationNew( i, temp, pahMoments
									, moments, Y, density, molarMass ) * hnenn / ( 2.0 * temp * density );
				for ( l = 0; l < fNSootMoments; ++l ) {
//					a[l+fOffsetMoments][ioff] -= dfdyUpwind( l+fOffsetMoments, ioff, SootCondFunc, nodeInfo, flame )
//												* hnenn;
//					a[l+fOffsetMoments][ioff] -= FillJacCondensationNew( l, i, temp, pahMoments, moments )
//												/ density * hnenn;
				}
			}

			//	surface growth
			if ( fSurfaceGrowth ) {
//				FillJacSurfGrowthNew( i, -1.0 / density, nodeInfo, flame );
			}

			//	surface oxidation
			if ( fSurfaceOxidation ) {
//				FillJacSootOxidationNew( i, -1.0 / density, nodeInfo, flame );
//				for ( l = 0; l < fNSootMoments; ++l ) {
//					a[l+fOffsetMoments][ioff] -= dfdyUpwind( l+fOffsetMoments, ioff, SootOxidationFunc, nodeInfo, flame )
//												* hnenn;
//				}
			}*/
		}
	}
	else if ( coordinate == kSimilarity ){
		int fVVelocity = flame->GetOffsetVVelocity();
		Double	oneOverRhoMuRef = 1.0 / ( flameNode->rhoInf 
						* flameNode->viscosityInf );
		Double	oneOverRhoA = 1.0 / ( flameNode->mixDensity[kCurr] 
						* flame->GetStrainRate() );
						
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i+fOffsetMoments;

			//	convection
			FillJacNonlinearConvectUpwind( fVVelocity, ioff, nodeInfo, 1.0, FALSE );

			//	diffusion
			if ( fSizeDepDiff ) {
				FillJacSootDiffusion( i, oneOverRhoMuRef, kSimilarity, flame, nodeInfo );
			}
			else {
				FillJacSootDiffusionNew( i, oneOverRhoMuRef, kSimilarity, flame, nodeInfo );
			}

			// thermophoresis
			if ( fThermoPhoresis ) {
				FillJacSootThermoPhoresis( i, oneOverRhoMuRef, kSimilarity, flame, nodeInfo );
			}

			// nucleation
			if ( fNucleation ) {
				a[tempOff][ioff] += 0.5 * oneOverRhoA * NucleationNew( i, temp, pahMoments ) 
												/ temp * hnenn;
			}

			// coagulation
			if ( fCoagulation ) {
				a[tempOff][ioff] += oneOverRhoA * dfdyUpwind( tempOff, ioff, SootCoagFunc, nodeInfo, flame )
										* hnenn;
//				a[ioff][ioff] += oneOverRhoA * FillJacSourceCoagulation( i, i, moments )
//											* hnenn;
//				a[ioff][ioff] += oneOverRhoA * dfdyUpwind( ioff, ioff, SootCoagFunc, nodeInfo, flame )
//											* hnenn;
				for ( l = 0; l < fNSootMoments; ++l ) {
					a[l+fOffsetMoments][ioff] += oneOverRhoA * dfdyUpwind( l+fOffsetMoments, ioff, SootCoagFunc, nodeInfo, flame )
												* hnenn;
//					a[l+fOffsetMoments][ioff] += oneOverRhoA * FillJacSourceCoagulation( i, l, moments )
//												* hnenn;
#ifdef DEBUGCOAG
					Double	coagAnalyt = FillJacSourceCoagulation( i, l, moments );
					Double	coagNum = dfdyUpwind( l+fOffsetMoments, ioff, SootCoagFunc, nodeInfo, flame );
					if ( fabs( ( coagAnalyt - coagNum ) / coagNum ) > 0.001 ) {
						fprintf( stderr, "gp = %d\tdCoag_%d/dM_%d Num = %g\tAnalyt= %g\n"
									, nodeInfo->gridPoint, i, l, coagNum, coagAnalyt );
					}
#endif
				}
			}

			//	condensation
			if ( fCondensation ) {
//				FillJacSourceCondensation( i, oneOverRhoA, nodeInfo, flame );
				a[tempOff][ioff] += oneOverRhoA * dfdyUpwind( tempOff, ioff, SootCondFunc, nodeInfo, flame )
											* hnenn;
				for ( l = 0; l < fNSootMoments; ++l ) {
//					a[l+fOffsetMoments][ioff] += oneOverRhoA * dfdyUpwind( l+fOffsetMoments, ioff, SootCondFunc, nodeInfo, flame )
//												* hnenn;
					a[l+fOffsetMoments][ioff] += density * oneOverRhoA * FillJacCondensationNew( l, i, temp, pahMoments, moments )
												* hnenn;
				}
			}

			//	surface growth
			if ( fSurfaceGrowth ) {
				FillJacSurfGrowthNew( i, oneOverRhoA, nodeInfo, flame );
			}

			//	surface oxidation
			if ( fSurfaceOxidation ) {
				FillJacSootOxidationNew( i, oneOverRhoA, nodeInfo, flame );
/*				for ( l = 0; l < fNSootMoments; ++l ) {
					a[l+fOffsetMoments][ioff] += oneOverRhoA * dfdyUpwind( l+fOffsetMoments, ioff, SootOxidationFunc, nodeInfo, flame )
												* hnenn;
				}*/
			}
		}
	}
	else if ( coordinate == kMixtureFraction ) {
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i+fOffsetMoments;

			//	diffusion
			if ( fSizeDepDiff ) {
//				FillJacSecondDerivCentral( ioff, ioff, 1.0, nodeInfo );
//				FillJacSootDiffusion( i, 1.0 / ( 2.0 * GetLewis1() ), kMixtureFraction, flame, nodeInfo );
			}
			else {
//				FillJacSootDiffusionNew( i, 1.0 / ( 2.0 * GetLewis1() ), kMixtureFraction, flame, nodeInfo );
			}

			// nucleation
			if ( fNucleation ) {
				a[tempOff][ioff] += 0.5 / density * NucleationNew( i, temp, pahMoments ) / temp * hnenn;
			}

			// coagulation
			if ( fCoagulation ) {
				a[tempOff][ioff] += dfdyUpwind( tempOff, ioff, SootCoagFunc, nodeInfo, flame )
											 / density * hnenn;
//				a[ioff][ioff] += FillJacSourceCoagulation( i, i, moments )
//											/ density * hnenn;
//				a[ioff][ioff] += dfdyUpwind( ioff, ioff, SootCoagFunc, nodeInfo, flame )
//											/ density * hnenn;
				for ( l = 0; l < fNSootMoments; ++l ) {
					a[l+fOffsetMoments][ioff] += dfdyUpwind( l+fOffsetMoments, ioff, SootCoagFunc, nodeInfo, flame )
												/ density * hnenn;
//					a[l+fOffsetMoments][ioff] += FillJacSourceCoagulation( i, l, moments )
//												/ density * hnenn;
				}
			}

			//	condensation
			if ( fCondensation ) {
				a[tempOff][ioff] += SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass ) 
						/ ( 2.0 * density * temp ) * hnenn;
				for ( l = 0; l < fNSootMoments; ++l ) {
					a[l+fOffsetMoments][ioff] += FillJacCondensationNew( l, i, temp, pahMoments, moments )
												 / density * hnenn;
				}
			}

			//	surface growth
			if ( fSurfaceGrowth ) {
				FillJacSurfGrowthNew( i, 1.0 / density, nodeInfo, flame );
			}

			//	surface oxidation
			if ( fSurfaceOxidation ) {
				FillJacSootOxidationNew( i, 1.0 / density, nodeInfo, flame );
/*				for ( l = 0; l < fNSootMoments; ++l ) {
					a[l+fOffsetMoments][ioff] += dfdyUpwind( l+fOffsetMoments, ioff, SootOxidationFunc, nodeInfo, flame )
												/ density * hnenn;
				}*/
			}
		}
	}
}

void T1DSoot::FillRHS( T1DFlamePtr flame, NodeInfoPtr nodeInfo, CoordType coordinate )
{

        int				i, ioff;
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double			temp = flameNode->temp[kCurr];
	Double			density = flameNode->mixDensity[kCurr];
	Double			*rhs = nodeInfo->rhs;
	Double			*pahMoments = flameNode->pahMoments;
	Double			*moments = flameNode->moments;
	Double			*Y = flameNode->Y[kCurr];
	Double			*kSoot = fSootRateCoeffs->vec;
	Double			*theSource = New1DArray( fNSootMoments );

	ComputeFractionalMoments( moments );
	ComputePhi( temp );
	
	ComputeSootRateCoeffs( kSoot, temp, flame->GetReaction() );
	ComputeCSootStar( kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr] );
	
	// convective term 
	if ( coordinate == kPhysical ) {
		int fVVelocity = flame->GetOffsetVVelocity();
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i + fOffsetMoments;

			//	convection
			rhs[ioff] += NonlinearConvectUpwind( nodeInfo->y[fVVelocity]
					, nodeInfo->yPrev[ioff], nodeInfo->y[ioff]
					, nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h );

			//	diffusion
			if ( fSizeDepDiff ) {
#ifndef NODIFF
				rhs[ioff] -= SootDiffusion( i, kPhysical, flame, nodeInfo );
#endif
			}
			else {
				rhs[ioff] -= SootDiffusionNew( i, kPhysical, flame, nodeInfo );
			}

/*			Double diff = SootDiffusion( i, kPhysical, flame, nodeInfo );*/
/*			Double tp = SootThermoPhoresis( i, kPhysical, flame, nodeInfo );*/
/*			if ( fabs( diff ) * 10.0 < fabs( tp ) ) {*/
/*				cerr << "T = " << temp << TAB*/
/*					<< "Diff = " << diff*/
/*					<< TAB << "ThermoPh = " << tp << NEWL;*/
/*			}*/
/*			if ( fabs( diff ) > fabs( tp ) * 10.0 ) {*/
/*				cerr << TAB << "T = " << temp << TAB*/
/*					<< "Diff = " << diff*/
/*					<< TAB << "ThermoPh = " << tp << NEWL;*/
/*			}*/

			// thermophoresis
			if ( fThermoPhoresis ) {
				rhs[ioff] -= SootThermoPhoresis( i, kPhysical, flame, nodeInfo );
			}
			rhs[ioff] -= flameNode->sootSource[i];

/*			// nucleation*/
/*			if ( fNucleation ) {*/
/*				theSource[i] += NucleationNew( i, temp, pahMoments ) / density;*/
/*				rhs[ioff] -= NucleationNew( i, temp, pahMoments ) / density;*/
/*			rhs[ioff] -= 0.5 * Y[f_A3R5AC] / molarMass[f_A3R5AC];*/
/*			}*/
/**/
/*			//	coagulation*/
/*			if ( fCoagulation ) {*/
/*				theSource[i] += SourceCoagulation( i ) / density;*/
/*				rhs[ioff] -= SourceCoagulation( i ) / density;*/
/*			}*/
/**/
/*			//	condensation*/
/*			if ( fCondensation ) {*/
/*				theSource[i] += SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass ) / density;*/
/*				rhs[ioff] -= SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass ) / density;*/
/*			}*/
/**/
/*			//	surface growth*/
/*			if ( fSurfaceGrowth ) {*/
/*				theSource[i] += SourceSurfGrowthNew( i, moments, Y, density, molarMass ) / density;*/
/*				rhs[ioff] -= SourceSurfGrowthNew( i, moments, Y, density, molarMass ) / density;*/
/*			}*/
/**/
/*			//	surface goxidation*/
/*			if ( fSurfaceOxidation ) {*/
/*				theSource[i] += SourceSootOxidationNew( i, moments, Y, density, molarMass ) / density;*/
/*				rhs[ioff] -= SourceSootOxidationNew( i, moments, Y, density, molarMass ) / density;*/
/*			}*/
/*			if ( theSource[i] > 1.0 && fabs( theSource[i] - flameNode->sootSource[i] ) / MAX( fabs( theSource[i] ), 1.0e-30 ) > 1.0e-8 ) {*/
/*				fprintf( stderr, "gp = %d\ti = %d\ttheSource = %g\tsootSource = %g\n"*/
/*					, nodeInfo->gridPoint, i, theSource[i], flameNode->sootSource[i] );*/
/*			}*/
		}
	}
	else if ( coordinate == kSimilarity ) {
		int fVVelocity = flame->GetOffsetVVelocity();
		Double	oneOverRhoMuRef = 1.0 / ( flameNode->rhoInf 
								* flameNode->viscosityInf );
		Double	oneOverRhoA = 1.0 / ( flameNode->mixDensity[kCurr] 
								* flame->GetStrainRate() );
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i + fOffsetMoments;
			//	convection
			rhs[ioff] += NonlinearConvectUpwind( nodeInfo->y[fVVelocity]
					, nodeInfo->yPrev[ioff], nodeInfo->y[ioff]
					, nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h, FALSE );

			//	diffusion
			if ( fSizeDepDiff ) {
				rhs[ioff] += oneOverRhoMuRef * SootDiffusion( i, kSimilarity, flame, nodeInfo );
			}
			else {
				rhs[ioff] += oneOverRhoMuRef * SootDiffusionNew( i, kSimilarity, flame, nodeInfo );
			}
			// thermophoresis
			if ( fThermoPhoresis ) {
				rhs[ioff] += oneOverRhoMuRef * SootThermoPhoresis( i, kSimilarity, flame, nodeInfo );
			}

			// nucleation
			if ( fNucleation ) {
				rhs[ioff] += oneOverRhoA * NucleationNew( i, temp, pahMoments );
			}
			//	coagulation
			if ( fCoagulation ) {
				fprintf( stderr, "###Attention: old Coagulation used here\n" );
				rhs[ioff] += oneOverRhoA * SourceCoagulation( i );
			}
			//	condensation
			if ( fCondensation ) {
				rhs[ioff] += oneOverRhoA * SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass );
			}
			//	surface growth
			if ( fSurfaceGrowth ) {
				rhs[ioff] += oneOverRhoA * SourceSurfGrowthNew( i, Y, density, molarMass );
			}
			//	surface oxidation
			if ( fSurfaceOxidation ) {
				rhs[ioff] += oneOverRhoA * SourceSootOxidationNew( i, moments, Y, density, molarMass );
			}
		}
	}
	else if ( coordinate == kMixtureFraction ) {
		for ( i = 0; i < fNSootMoments; ++i ) {
			ioff = i + fOffsetMoments;

			//	diffusion
			if ( fSizeDepDiff ) {
//				rhs[ioff] += SecondDeriv( nodeInfo->yPrev[ioff], nodeInfo->y[ioff]
//								, nodeInfo->yNext[ioff], nodeInfo->hm, nodeInfo->h );
//				rhs[ioff] += SootDiffusion( i, kMixtureFraction, flame, nodeInfo )
//							/ ( 2.0 * GetLewis1() );
			}
			else {
//				rhs[ioff] += SootDiffusionNew( i, kMixtureFraction, flame, nodeInfo )
//						/ ( 2.0 * GetLewis1() );
			}

			// nucleation
			if ( fNucleation ) {
//				fprintf( stderr, "Nuc = %g\n", NucleationNew( i, temp, pahMoments ) );
				rhs[ioff] += NucleationNew( i, temp, pahMoments ) / density;
			}

			//	coagulation
			if ( fCoagulation ) {
//				rhs[ioff] += SourceCoagulation( i ) / density;
				rhs[ioff] += SourceCoagulationNew( i, temp, moments ) / density;
			}

			//	condensation
			if ( fCondensation ) {
				rhs[ioff] += SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass ) / density;
			}

			//	surface growth
			if ( fSurfaceGrowth ) {
//				rhs[ioff] += SourceSurfGrowthNew( i, Y, density, molarMass ) / density;
				rhs[ioff] += SourceSurfGrowthNew( i, moments, Y, density, molarMass ) / density;
			}

			//	surface goxidation
			if ( fSurfaceOxidation ) {
				rhs[ioff] += SourceSootOxidationNew( i, moments, Y, density, molarMass ) / density;
			}
		}
	}
	else {
		fprintf( stderr, "#error: invalid coordinate type %d\n", coordinate );
		exit( 2 );
	}
	Free1DArray( theSource );
}

void T1DSoot::FillSource( Double *source, T1DFlamePtr flame )
{
  //  printf("FillSource");
	int				i;
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double			temp = flameNode->temp[kCurr];
	Double			density = flameNode->mixDensity[kCurr];
	Double			*pahMoments = flameNode->pahMoments;
	Double			*moments = flameNode->moments;
	Double			*Y = flameNode->Y[kCurr];
	Double			*kSoot = fSootRateCoeffs->vec;

	ComputeFractionalMoments( moments );
	ComputePhi( temp );
	
	ComputeSootRateCoeffs( kSoot, temp, flame->GetReaction() );
	ComputeCSootStar( kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr] );
	
	for ( i = 0; i < fNSootMoments; ++i ) {
		source[i] = 0.0;

		// nucleation
		if ( fNucleation ) {
//			source[i] += MAGICSOURCE * Nucleation( i, temp, pahMoments );
			source[i] += MAGICSOURCE * NucleationNew( i, temp, pahMoments );
//			Double	C_A3R5AC = density * Y[f_A3R5AC] / molarMass[f_A3R5AC];
//			source[i] += pow( 9.0, (Double) i ) * 1.0e10 * exp( -2405.5 / temp ) * C_A3R5AC;
		}

		//	coagulation
		if ( fCoagulation ) {
//			if ( fSurfDepCoag ) {
//				source[i] += MAGICSOURCE * SourceSurfDepCoag( i, moments, Y, temp, density, molarMass );
//			}
//			else {
				source[i] += MAGICSOURCE * SourceCoagulationNew( i, temp, moments );
			}
//		}

		//	condensation
		if ( fCondensation ) {
			source[i] += MAGICSOURCE * SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass );
		}

		//	surface growth
		if ( fSurfaceGrowth ) {
			source[i] += MAGICSOURCE * SourceSurfGrowthNew( i, moments, Y, density, molarMass );
		}

		//	surface oxidation
		if ( fSurfaceOxidation ) {
			source[i] += MAGICSOURCE * SourceSootOxidationNew( i, moments, Y, density, molarMass );
		}
	}
}

/*void T1DSoot::FillJacSootConvection( int i, T1DFlamePtr flame, NodeInfoPtr nodeInfo, Flag velocityPositive )
{
//	velocityPositive should have the value FALSE, if 'V' has the negative 
//	direction of the physical velocity

// fills the jacobian with     V * d(M_i/rho)/dy

	
	int		iOff = i + fOffsetMoments;
	int		vOff = flame->GetOffsetVVelocity();
	int		tempOff = flame->GetOffsetTemperature();
	Double	*rho = flame->GetFlameNode()->mixDensity;
	Double	coeff;
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*yPrev = nodeInfo->yPrev;
	Double	**a = nodeInfo->a;
	Double	V = y[vOff];

	if ( ( V > 0.0 && velocityPositive ) || ( V < 0.0 && !velocityPositive ) ) {
		coeff = h * ( h + hm );
		a[vOff][iOff] += coeff * ( y[iOff] / rho[kCurr] - yPrev[iOff] / rho[kPrev] );
		a[iOff][iOff] += coeff * V / rho[kCurr];
//		a[tempOff][iOff] += coeff * V * y[iOff] / ( rho[kCurr] * y[tempOff] );
		if ( !nodeInfo->firstPoint ) {
			nodeInfo->c[iOff][iOff] -= coeff * V / rho[kPrev];
//			nodeInfo->c[tempOff][iOff] -= coeff * V * yPrev[iOff] / ( rho[kPrev] * yPrev[tempOff] );
		}
#ifdef FLUXBC
		else if ( !fSizeDepDiff ) {
			Double	c = rho[kPrev] * flame->GetFlameNode()->diffSoot[kPrev] / ( yPrev[vOff] * hm );
			a[iOff][iOff] -= V * coeff / rho[kCurr] * c / ( 1.0 + c ); 
		}
#endif
	}
	else {
		coeff = hm * ( h + hm );
		a[vOff][iOff] += coeff * ( yNext[iOff] / rho[kNext] - y[iOff] / rho[kCurr] );
		a[iOff][iOff] -= coeff * V / rho[kCurr];
//		a[tempOff][iOff] -= coeff * V * y[iOff] / ( rho[kCurr] * y[tempOff] );
		if ( !nodeInfo->lastPoint ) {
			nodeInfo->b[iOff][iOff] += coeff * V / rho[kNext];
//			nodeInfo->b[tempOff][iOff] += coeff * V * yNext[iOff] / ( rho[kNext] * yNext[tempOff] );
		}
#ifdef FLUXBC
		else if ( !fSizeDepDiff ) {
			Double	c = rho[kNext] * flame->GetFlameNode()->diffSoot[kNext] / ( yNext[vOff] * h );
			a[iOff][iOff] -= V * coeff / rho[kCurr] * c / ( 1.0 - c ); 
		}
#endif
	}
}
*/

void T1DSoot::FillJacSootDiffusion( int r, Double constCoeff, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
// fills the jacobian with     constCoeff * d/dy ( rho * diffusivity * d(M_(r-2/3)/rho)/dy )

	int		i;
	int		rOff = r + fOffsetMoments;
	int		tempOff = flame->GetOffsetTemperature();
	TFlameNodePtr	flameNode = flame->GetFlameNode();
/*	Double	diff = flameNode->diffusivity[f_A1];
	Double	diffNext = flameNode->diffusivityNext[f_A1];
	Double	diffPrev = flameNode->diffusivityPrev[f_A1];*/
	Double	diff = flameNode->diffSoot[kCurr];
	Double	diffNext = flameNode->diffSoot[kNext];
	Double	diffPrev = flameNode->diffSoot[kPrev];
	Double	*rho = flameNode->mixDensity;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*yPrev = nodeInfo->yPrev;
	Double	**a = nodeInfo->a;
	Double	**b = nodeInfo->b;
	Double	**c = nodeInfo->c;
	Double	*moments = &y[fOffsetMoments];
	Double	*momentsNext = &yNext[fOffsetMoments];
	Double	*momentsPrev = &yPrev[fOffsetMoments];
	Double	fracIndex = r - 2.0 / 3.0;
	Double	fracMom = FracMom2( fracIndex, moments );
	Double	fracMomPrev = FracMom2( fracIndex, momentsPrev );
	Double	fracMomNext = FracMom2( fracIndex, momentsNext );
	Double	diffPlusHm, diffMinusH, fact, factNext, factPrev;

	if ( coordinate == kPhysical ) {
		diffPlusHm = nodeInfo->hm * constCoeff * ( diff * rho[kCurr]
											+ diffNext * rho[kNext] );
		diffMinusH = nodeInfo->h * constCoeff * ( diffPrev * rho[kPrev]
						+ diff * rho[kCurr] );
	}
	else if ( coordinate == kSimilarity ){
		diffPlusHm = nodeInfo->hm * constCoeff * ( diff * rho[kCurr] * rho[kCurr]
											+ diffNext * rho[kNext] * rho[kNext] );
		diffMinusH = nodeInfo->h * constCoeff * ( diffPrev * rho[kPrev] * rho[kPrev]
						+ diff * rho[kCurr] * rho[kCurr] );
	}
	else if ( coordinate == kMixtureFraction ) {
		diffPlusHm = 2.0 * constCoeff * nodeInfo->hm;
		diffMinusH = 2.0 * constCoeff * nodeInfo->h;
	}

	fact = ( diffPlusHm + diffMinusH ) * fracMom;
	factNext = diffPlusHm * fracMomNext;
	factPrev = diffMinusH * fracMomPrev;
	
#define FULLSOOTDIFFJAC

#ifdef FULLSOOTDIFFJAC
	for ( i = 0; i < fNSootMoments; ++i ) {
		a[i+fOffsetMoments][rOff] -= fact * GetAlphaI2( i, fracIndex ) 
											* rho[kCurr] / moments[i];
		if ( !nodeInfo->lastPoint ) {
			b[i+fOffsetMoments][rOff] += factNext * GetAlphaI2( i, fracIndex ) 
											* rho[kNext] / momentsNext[i];
		}
		if ( !nodeInfo->firstPoint ) {
			c[i+fOffsetMoments][rOff] += factPrev * GetAlphaI2( i, fracIndex ) 
											* rho[kPrev] / momentsPrev[i];
		}
	}
#else
	a[rOff][rOff] -= fact * GetAlphaI2( r, fracIndex ) * rho[kCurr] / moments[r];
	if ( !nodeInfo->lastPoint ) {
		b[rOff][rOff] += factNext * GetAlphaI2( r, fracIndex ) * rho[kNext] / momentsNext[r];
	}
	if ( !nodeInfo->firstPoint ) {
		c[rOff][rOff] += factPrev * GetAlphaI2( r, fracIndex ) * rho[kPrev] / momentsPrev[r];
	}
#endif
}

void T1DSoot::FillJacSootDiffusionNew( int r, Double constCoeff, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
// fills the jacobian with     constCoeff * d/dy ( rho * diffusivity * d(M_r/rho)/dy )

	int		rOff = r + fOffsetMoments;
	int		tempOff = flame->GetOffsetTemperature();
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double	diff = flameNode->diffSoot[kCurr];
	Double	diffNext = flameNode->diffSoot[kNext];
	Double	diffPrev = flameNode->diffSoot[kPrev];
	Double	*rho = flameNode->mixDensity;
	Double	**a = nodeInfo->a;
	Double	**b = nodeInfo->b;
	Double	**c = nodeInfo->c;
	Double	diffPlusHm, diffMinusH, fact, factNext, factPrev;

	if ( coordinate == kPhysical ) {
		diffPlusHm = nodeInfo->hm * constCoeff * ( diff * rho[kCurr]
											+ diffNext * rho[kNext] );
		diffMinusH = nodeInfo->h * constCoeff * ( diffPrev * rho[kPrev]
						+ diff * rho[kCurr] );
	}
	else if ( coordinate == kSimilarity ){
		diffPlusHm = nodeInfo->hm * constCoeff * ( diff * rho[kCurr] * rho[kCurr]
											+ diffNext * rho[kNext] * rho[kNext] );
		diffMinusH = nodeInfo->h * constCoeff * ( diffPrev * rho[kPrev] * rho[kPrev]
						+ diff * rho[kCurr] * rho[kCurr] );
	}
	else if ( coordinate == kMixtureFraction ) {
		diffPlusHm = 2.0 * constCoeff * nodeInfo->hm;
		diffMinusH = 2.0 * constCoeff * nodeInfo->h;
	}

	fact = ( diffPlusHm + diffMinusH );
	factNext = diffPlusHm;
	factPrev = diffMinusH;
	
	a[rOff][rOff] -= fact;
	if ( !nodeInfo->lastPoint ) {
		b[rOff][rOff] += factNext;
	}
#ifdef FLUXBC
	else if ( !fSizeDepDiff ) {
		Double	c = rho[kNext] * diffNext / ( yNext[flame->GetOffsetVVelocity()] * nodeInfo->h );
		a[rOff][rOff] -= factNext * rho[kPrev]/rho[kCurr] * c / ( 1.0 - c ); 
	}
#endif
	if ( !nodeInfo->firstPoint ) {
		c[rOff][rOff] += factPrev;
	}
#ifdef FLUXBC
	else if ( !fSizeDepDiff ) {
		Double	c = rho[kPrev] * diffPrev / ( yPrev[flame->GetOffsetVVelocity()] * nodeInfo->hm );
		a[rOff][rOff] += factPrev * rho[kPrev]/rho[kCurr] * c / ( 1.0 + c ); 
	}
#endif
}

Double T1DSoot::SootConvection( int i, T1DFlamePtr flame, NodeInfoPtr nodeInfo, Flag velocityPositive )
{
//	velocityPositive should have the value FALSE, if 'V' has the negative 
//	direction of the physical velocity

// returns     V * d(M_i/rho)/dy

	int		iOff = i + fOffsetMoments;
	Double	*y = nodeInfo->y;
	Double	*yNext = nodeInfo->yNext;
	Double	*yPrev = nodeInfo->yPrev;
	Double	*rho = flame->GetFlameNode()->mixDensity;
	Double	V = y[flame->GetOffsetVVelocity()];

	if ( ( V > 0.0 && velocityPositive ) || ( V < 0.0 && !velocityPositive ) ) {
		return ( V * ( y[iOff]/rho[kCurr] - yPrev[iOff]/rho[kPrev] ) / nodeInfo->hm );
	}
	else {
		return ( V * ( yNext[iOff]/rho[kNext] - y[iOff]/rho[kCurr] ) / nodeInfo->h );
	}
}

Double T1DSoot::SootDiffusion( int r, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
/*	Double	diff = flameNode->diffusivity[f_A1];
	Double	diffNext = flameNode->diffusivityNext[f_A1];
	Double	diffPrev = flameNode->diffusivityPrev[f_A1];*/
	Double	diff = flameNode->diffSoot[kCurr];
	Double	diffNext = flameNode->diffSoot[kNext];
	Double	diffPrev = flameNode->diffSoot[kPrev];
	Double	*rho = flameNode->mixDensity;
	Double	fracIndex = r - 2.0 / 3.0;
	Double	fracMom = FracMom2( fracIndex, &nodeInfo->y[fOffsetMoments] );
	Double	fracMomPrev = FracMom2( fracIndex, &nodeInfo->yPrev[fOffsetMoments] );
	Double	fracMomNext = FracMom2( fracIndex, &nodeInfo->yNext[fOffsetMoments] );
	Double	diffPlusHm, diffMinusH;
	
	if ( coordinate == kPhysical ) {
		diffPlusHm = nodeInfo->hm * ( diff * rho[kCurr]
											+ diffNext * rho[kNext] );
		diffMinusH = nodeInfo->h * ( diffPrev * rho[kPrev]
						+ diff * rho[kCurr] );
	}
	else if ( coordinate == kSimilarity ) {
		diffPlusHm = nodeInfo->hm * ( diff * rho[kCurr] * rho[kCurr]
											+ diffNext * rho[kNext] * rho[kNext] );
		diffMinusH = nodeInfo->h * ( diffPrev * rho[kPrev] * rho[kPrev]
						+ diff * rho[kCurr] * rho[kCurr] );
	}
	else if ( coordinate == kMixtureFraction ) {
		diffPlusHm = 2.0 * nodeInfo->hm;
		diffMinusH = 2.0 * nodeInfo->h;
	}
	
	return ( diffPlusHm * ( fracMomNext - fracMom ) 
			+ diffMinusH * ( fracMomPrev - fracMom ) ) 
				/ nodeInfo->hnenn;
}

Double T1DSoot::SootDiffusionNew( int r, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
	int		rOff = r + fOffsetMoments;
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double	diff = flameNode->diffSoot[kCurr];
	Double	diffNext = flameNode->diffSoot[kNext];
	Double	diffPrev = flameNode->diffSoot[kPrev];
	Double	*rho = flameNode->mixDensity;
	Double	mom = nodeInfo->y[rOff];
	Double	momPrev = nodeInfo->yPrev[rOff];
	Double	momNext = nodeInfo->yNext[rOff];
	Double	diffPlusHm, diffMinusH;
	
	if ( coordinate == kPhysical ) {
		diffPlusHm = nodeInfo->hm * ( diff * rho[kCurr]
											+ diffNext * rho[kNext] );
		diffMinusH = nodeInfo->h * ( diffPrev * rho[kPrev]
						+ diff * rho[kCurr] );
	}
	else if ( coordinate == kSimilarity ) {
		diffPlusHm = nodeInfo->hm * ( diff * rho[kCurr] * rho[kCurr]
											+ diffNext * rho[kNext] * rho[kNext] );
		diffMinusH = nodeInfo->h * ( diffPrev * rho[kPrev] * rho[kPrev]
						+ diff * rho[kCurr] * rho[kCurr] );
	}
	else if ( coordinate == kMixtureFraction ) {
		diffPlusHm = 2.0 * nodeInfo->hm;
		diffMinusH = 2.0 * nodeInfo->h;
	}

	return ( diffPlusHm * ( momNext - mom ) 
			+ diffMinusH * ( momPrev - mom ) ) 
				/ nodeInfo->hnenn;
}

Double T1DSoot::SootThermoPhoresis( int r, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
	// returns 			d/dy ( 0.55 * nu / T * M_r * dT/dy )	

	const Double	A = 0.9;	// accommodation coefficient following
								// R. J. Santoro: The Transport and Growth
								// of Soot Particles ...
								// Comb. Sci. Tech., 1987, Vol. 53 p. 89
	static Double	pi = 4.0 * atan( 1.0 );
	const Double	fact = 3.0 / ( 4.0 * ( 1.0 + pi * A / 8.0 ) );
	int		rOff = r + fOffsetMoments;
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double	mu = flameNode->mixViscosity[kCurr];
	Double	muNext = flameNode->mixViscosity[kNext];
	Double	muPrev = flameNode->mixViscosity[kPrev];
	Double	M = nodeInfo->y[rOff];
	Double	MNext = nodeInfo->yNext[rOff];
	Double	MPrev = nodeInfo->yPrev[rOff];
	Double	*temp = flameNode->temp;
	Double	*rho = flameNode->mixDensity;
	Double	diffPlusHm, diffMinusH;
	
	if ( coordinate == kPhysical ) {
		diffPlusHm = nodeInfo->hm * fact 
				* ( mu * M / temp[kCurr]
					+ muNext * MNext / temp[kNext] );
		diffMinusH = nodeInfo->h * fact 
				* ( mu * M / temp[kCurr]
					+ muPrev * MPrev / temp[kPrev] );
	}
	else {
		diffPlusHm = nodeInfo->hm * fact 
				* ( rho[kCurr] * mu * M / temp[kCurr]
					+ rho[kNext] * muNext * MNext / temp[kNext] );
		diffMinusH = nodeInfo->h * fact 
				* ( rho[kCurr] * mu * M / temp[kCurr]
					+ rho[kPrev] * muPrev * MPrev /temp[kPrev] );
	}
	return ( diffPlusHm * ( temp[kNext] - temp[kCurr] ) 
			+ diffMinusH * ( temp[kPrev] - temp[kCurr] ) ) 
				/ nodeInfo->hnenn;
}

void T1DSoot::FillJacSootThermoPhoresis( int r, Double constCoeff, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
// fills the jacobian with     constCoeff * d/dy ( 0.55 * nu / T * M_r * dT/dy )

	const Double	A = 0.9;	// accommodation coefficient following
								// R. J. Santoro: The Transport and Growth
								// of Soot Particles ...
								// Comb. Sci. Tech., 1987, Vol. 53 p. 89
	static Double	pi = 4.0 * atan( 1.0 );
	const Double	fact = 3.0 / ( 4.0 * ( 1.0 + pi * A / 8.0 ) );
	constCoeff *= fact;
	int		rOff = r + fOffsetMoments;
	int		tempOff = flame->GetOffsetTemperature();
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double	mu = flameNode->mixViscosity[kCurr];
	Double	muNext = flameNode->mixViscosity[kNext];
	Double	muPrev = flameNode->mixViscosity[kPrev];
	Double	M = nodeInfo->y[rOff];
	Double	MNext = nodeInfo->yNext[rOff];
	Double	MPrev = nodeInfo->yPrev[rOff];
	Double	hm = nodeInfo->hm;
	Double	h = nodeInfo->h;
	Double	*temp = flameNode->temp;
	Double	*rho = flameNode->mixDensity;
	Double	**a = nodeInfo->a;
	Double	**b = nodeInfo->b;
	Double	**c = nodeInfo->c;
	Double	diffPlusHm, diffMinusH;
	Double	coeffPlus, coeffCurr, coeffMinus;

	if ( coordinate == kPhysical ) {
		diffPlusHm = hm * constCoeff 
				* ( mu * M / temp[kCurr]
					+ muNext * MNext / temp[kNext] );
		diffMinusH = h * constCoeff 
				* ( mu * M / temp[kCurr]
					+ muPrev * MPrev / temp[kPrev] );

		coeffPlus = constCoeff * muNext / temp[kNext]
					* hm * ( temp[kNext] - temp[kCurr] );
		coeffMinus = constCoeff * muPrev / temp[kPrev]
					* h * ( temp[kPrev] - temp[kCurr] );
		coeffCurr = constCoeff * mu / temp[kCurr]
					* ( hm * ( temp[kNext] - temp[kCurr] ) 
						+ h * ( temp[kPrev] - temp[kCurr] ) );
	}
	else {
		diffPlusHm = hm * constCoeff 
				* ( rho[kCurr] * mu * M / temp[kCurr]
					+ rho[kNext] * muNext * MNext / temp[kNext] );
		diffMinusH = h * constCoeff 
				* ( rho[kCurr] * mu * M / temp[kCurr]
					+ rho[kPrev] * muPrev * MPrev /temp[kPrev] );

		coeffPlus = constCoeff * rho[kNext] * muNext / temp[kNext]
					* hm * ( temp[kNext] - temp[kCurr] );
		coeffMinus = constCoeff * rho[kPrev] * muPrev / temp[kPrev]
					* h * ( temp[kPrev] - temp[kCurr] );
		coeffCurr = constCoeff * rho[kCurr] * mu / temp[kCurr]
					* ( hm * ( temp[kNext] - temp[kCurr] ) 
						+ h * ( temp[kPrev] - temp[kCurr] ) );
	}

	a[tempOff][rOff] -= diffPlusHm + diffMinusH;
	a[rOff][rOff] += coeffCurr;
	if ( !nodeInfo->lastPoint ) {
		b[tempOff][rOff] += diffPlusHm;
		b[rOff][rOff] += coeffPlus;
	}
	if ( !nodeInfo->firstPoint ) {
		c[tempOff][rOff] += diffMinusH;
		c[rOff][rOff] += coeffMinus;
	}
}

void T1DSoot::PrintFracMom( T1DFlamePtr flame )
{
	TNewtonPtr	bt = flame->GetSolver()->bt;
	int			k;
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		*fm = fFracMom->vec;
	Double		**mo = fMoments->mat;
	Double		dummy;
	int			counter;
	

	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( flame->GetOutFileBuff(), "%sFracMoments_%d.dout", flame->GetOutputPath(), counter );
	if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) {
		cerr << "#warning: unable to open file" << flame->GetOutFileBuff() << NEWL;
		return;
	}

	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "y" );
	fprintf( fp, "\tM_%-6s\tM_%-6s\tM_%-6s\tM_%-6s\tM_%-6s\tM_%-6s\tM_%-6s\tM_%-6s"
				, "km1_3", "0", "k02_3", "1", "k05_3", "2", "k08_3", "3" );
	
	for ( k = 0; k < nGridPoints; ++k ) {
		ComputeFractionalMoments( mo[k] );
		fprintf( fp, "\n%-9E", x[k] );
		fprintf( fp, "\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E"
		, fm[km1_3], mo[k][0], fm[k02_3], mo[k][1], fm[k05_3], mo[k][2], fm[k08_3], mo[k][3] );
	}

	fclose( fp );
}

void T1DSoot::PrintPhi( T1DFlamePtr flame )
{
	TNewtonPtr	bt = flame->GetSolver()->bt;
	int			k;
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		*phi = fPhi->vec;
	Double		**mo = fMoments->mat;
	Double		*temp = flame->GetTemperature()->vec;
	Double		dummy;
	int			counter;
	

	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( flame->GetOutFileBuff(), "%sPhi_%d.dout", flame->GetOutputPath(), counter );
	if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) {
		cerr << "#warning: unable to open file" << flame->GetOutFileBuff() << NEWL;
		return;
	}

	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "y" );
	fprintf( fp, "\tM_%-6s\tPhi_%-6s\tM_%-6s\tPhi_%-6s\tM_%-6s\tPhi_%-6s\tM_%-6s"
				, "0", "kh00", "1", "kh11", "2", "kh12", "3" );
	
	for ( k = 0; k < nGridPoints; ++k ) {
		ComputeFractionalMoments( mo[k] );
		ComputePhi( temp[k] );
		fprintf( fp, "\n%-9E", x[k] );
		fprintf( fp, "\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E\t%-9E"
		, mo[k][0], phi[kh00], mo[k][1], phi[kh11], mo[k][2], phi[kh12], mo[k][3] );
	}

	fclose( fp );
}

void T1DSoot::PrintDiffusivity( T1DFlamePtr flame )
{
	TNewtonPtr	bt = flame->GetSolver()->bt;
	int			k;
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		*diff = fSootDiff->vec;
	Double		**diffA1 = flame->GetSpecies()->GetDiffusivity()->mat;
	Double		dummy;
	int			counter;
	

	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( flame->GetOutFileBuff(), "%sdiff_%d.dout", flame->GetOutputPath(), counter );
	if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) {
		cerr << "#warning: unable to open file" << flame->GetOutFileBuff() << NEWL;
		return;
	}

	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "y" );
	fprintf( fp, "\t%-6s\t%-6s"
				, "D_Soot", "D_A1" );
	
	fprintf( fp, "\n%-9E\t%-9E\t%-9E"
	, left, diff[kPrev], diffA1[kPrev][f_A1] );

	for ( k = 0; k < nGridPoints; ++k ) {
		fprintf( fp, "\n%-9E", x[k] );
		fprintf( fp, "\t%-9E\t%-9E"
		, diff[k], diffA1[k][f_A1] );
	}

	fprintf( fp, "\n%-9E\t%-9E\t%-9E"
	, right, diff[nGridPoints], diffA1[nGridPoints][f_A1] );

	fclose( fp );
}

void T1DSoot::PostIter( T1DFlamePtr flame )
{
	int			i;
	TNewtonPtr	bt = flame->GetSolver()->bt;
	TGridPtr 	currGrid = bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = currGrid->GetNGridPoints();
	Double		*yLeft = currGrid->GetYLeft()->vec;
	Double		*yRight = currGrid->GetYRight()->vec;
	Double		**y = currGrid->GetY()->mat;
#ifdef FLUXBC
	Double		*x = currGrid->GetX()->vec;
	Double		hFirst = x[0] - bt->GetLeft();
	Double		hLast = bt->GetRight() - x[nGridPoints-1];
	Double		*rho0 = &flame->GetProperties()->GetDensity()->vec[0];
	Double		*rhoN = &flame->GetProperties()->GetDensity()->vec[nGridPoints-1];
	Double		cLeft = rho0[kPrev] * GetSootDiff()->vec[kPrev] 
						/ ( yLeft[flame->GetOffsetVVelocity()] * hFirst );
	Double		cRight = rhoN[kNext] * GetSootDiff()->vec[nGridPoints] 
						/ ( yRight[flame->GetOffsetVVelocity()] * hLast );
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			pressure = flame->GetPressure();

	flame->SetFlameNode( 0 );
	flame->ComputeProperties( flameNode, flameNode->temp[kCurr], flameNode->Y[kCurr], pressure );
	flame->SetFlameNode( nGridPoints-1 );
	flame->ComputeProperties( flameNode, flameNode->temp[kCurr], flameNode->Y[kCurr], pressure );
#endif

	//CheckSolution( nGridPoints, y );
#ifdef FLUXBC
	if ( fSizeDepDiff ) {
		int counter = 0;
		Double	normLeft, normRight, val;
		Double	*MLeft = new Double[fNSootMoments];
		Double	*MRight = new Double[fNSootMoments];
		for ( i = 0; i < fNSootMoments; ++i ) {
			yLeft[i+fOffsetMoments] = MAX( yLeft[i+fOffsetMoments], SMALLSOOT );
			yRight[i+fOffsetMoments] = MAX( yRight[i+fOffsetMoments], SMALLSOOT );
			MLeft[i] = yLeft[i+fOffsetMoments];
			MRight[i] = yRight[i+fOffsetMoments];
		}
		
		do {
			++counter;
			normLeft = normRight = 0.0;
			for ( i = 0; i < fNSootMoments; ++i ) {
	/*			yLeft[i+fOffsetMoments] = cLeft * ( rho0[kPrev]/rho0[kCurr] 
						* FracMom2( i-2.0/3.0, &y[0][fOffsetMoments] ) 
						- FracMom2( i-2.0/3.0, &yLeft[fOffsetMoments] ) )
	//										/ ( 1.0 + cLeft )
											;*/
	/*			yRight[i+fOffsetMoments] = cRight * (
							FracMom2( i-2.0/3.0, &yRight[fOffsetMoments] )
							- FracMom2( i-2.0/3.0, &y[nGridPoints-1][fOffsetMoments] ) 
								* rhoN[kNext]/rhoN[kCurr] ) 
	//										/ ( 1.0 - cRight )
											;*/
				Double	c1, c2, f, fPrime;
	
	/*			c1 = cLeft / ( 1.0 + cLeft );
				c2 = rho0[kPrev]/rho0[kCurr] 
						* FracMom2( i-2.0/3.0, &y[0][fOffsetMoments] );
				f = yLeft[i+fOffsetMoments] 
						- c1 * ( c2 - FracMom2( i-2.0/3.0, &yLeft[fOffsetMoments] ) );
				fPrime = 1.0 + c1 * GetDerivFracMom( i, i-2.0/3.0, &yLeft[fOffsetMoments] );
				yLeft[i+fOffsetMoments] -= f / fPrime;*/
	
				c1 = cRight / ( 1.0 - cRight );
				c2 = FracMom2( i-2.0/3.0, &y[nGridPoints-1][fOffsetMoments] );
				f = yRight[i+fOffsetMoments] 
						- c1 * ( FracMom2( i-2.0/3.0, &yRight[fOffsetMoments] ) - c2 );
				fPrime = 1.0 - c1 * GetDerivFracMom( i, i-2.0/3.0, &yRight[fOffsetMoments] );
				yRight[i+fOffsetMoments] -= f / fPrime;
				
				val = ( yLeft[i+fOffsetMoments] - MLeft[i] ) / MLeft[i];
				normLeft += val * val;
				val = ( yRight[i+fOffsetMoments] - MRight[i] ) / MRight[i];
				normRight += val * val;
				if ( i == 1 ) {
	/*				cerr << "i = " << i << NEWL << "left = " << yLeft[i+fOffsetMoments]
						<< TAB << "old = " << MLeft[i]
						<< TAB << "diff = " << yLeft[i+fOffsetMoments] - MLeft[i]
						<< TAB << "normLeft = " << sqrt( normLeft ) << NEWL;*/
					cerr << "right = " << yRight[i+fOffsetMoments] 
						<< TAB << "old = " << MRight[i]
						<< TAB << "diff = " << yRight[i+fOffsetMoments] - MRight[i]
						<< TAB << "normRight = " << sqrt( normRight ) << NEWL;
				}
				MLeft[i] = yLeft[i+fOffsetMoments];
				MRight[i] = yRight[i+fOffsetMoments];
			}
		} while ( ( sqrt( normLeft ) > 0.01 || sqrt( normRight ) > 0.01 ) && counter < 100 );
		cerr << "#counter = " << counter << NEWL;
		if ( counter >=100 ) {
			cerr << "#error: iteration of boundary values not converged" << NEWL;
		}
		for ( i = 0; i < fNSootMoments; ++i ) {
			switch( i ) {
				case 0:
					yLeft[i+fOffsetMoments] = SMALLSOOT;
	//				yRight[i+fOffsetMoments] = SMALLSOOT;
					break;
				default:
					yLeft[i+fOffsetMoments] = 9.0 * yLeft[i+fOffsetMoments-1];
	//				yRight[i+fOffsetMoments] = 9.0 * yRight[i+fOffsetMoments-1];
					break;
			}
		}
		delete MLeft;
		delete MRight;
	else {
		for ( i = 0; i < fNSootMoments; ++i ) {
			yLeft[i+fOffsetMoments] = cLeft * rho0[kPrev]/rho0[kCurr] * y[0][i+fOffsetMoments] 
										/ ( 1.0 + cLeft );
			yRight[i+fOffsetMoments] = - cRight * rhoN[kNext]/rhoN[kCurr] * y[nGridPoints-1][i+fOffsetMoments] 
										/ ( 1.0 - cRight );
		}
	}
#else
	for ( i = 0; i < fNSootMoments; ++i ) {
		switch ( i ) {
			case 0:
				yLeft[i+fOffsetMoments] = SMALLSOOT;
				yRight[i+fOffsetMoments] = SMALLSOOT;
				break;
			default:
				yLeft[i+fOffsetMoments] = 9 * yLeft[i+fOffsetMoments-1];
				yRight[i+fOffsetMoments] = 9 * yRight[i+fOffsetMoments-1];
				break;
		}
	}
#endif
}

void T1DSoot::CheckSolution( int nGridPoints, Double **y )
{
	for ( int k = 0; k < nGridPoints; ++k ) {
		TSoot::CheckSolution( &y[k][fOffsetMoments] );
	}
}

Double T1DSoot::FillJacSourceCoagulation( int i, int l, Double *moments )
{
//	returns dS_i/dM_l
		
	switch ( i ) {
		case 0:
			return -0.5 * GetDerivPhi( l, kh00, moments );
		case 1:
			return 0.0;
		case 2:
			return GetDerivPhi( l, kh11, moments );
		case 3:
			return 3.0 * GetDerivPhi( l, kh12, moments );
		default:
			cerr << "#error: no need to compute derivative of coagulation source_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

void T1DSoot::PrintFlameletFile( int gridPoints, T1DFlamePtr flame, FILE *fp )
{
	int 			i, k;
	Double			**moments = fMoments->mat;
//	Double			fvFact = flame->GetSpecies()->GetMolarMass()->vec[f_A1] / fSootDensity;
	Double			fvFact = fMolarMassSoot / fSootDensity;
	Double			numdensfact = AVOGADRO * 1.0e-6; // kmole/m^3 -> parts/cm^3
	Double			partdiamfact = 3.48e-01; // in [nm]
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double			*density = flame->GetProperties()->GetDensity()->vec;
	Double			*mixMolarMass = flame->GetProperties()->GetMolarMass()->vec;
	Double			**Y = flame->GetMassFracs()->mat;
	Double			*temp = flame->GetTemperature()->vec;
	Double			**pahMoments = fPAHMoments->mat;
	Double			*kSoot = fSootRateCoeffs->vec;

// write moments

	for ( i = 0; i < fNSootMoments; ++i ) {
		fprintf( fp, "conc-SOOT%d\n", i );
		for ( k = 0; k < gridPoints+2; ++k ) {
			fprintf( fp, "\t%-.6e", moments[k-1][i] );
			if ( (k+1) % 5 == 0 ) {
				fprintf( fp, "\n" );
			}
		}
		if ( k % 5 ) {
			fprintf( fp, "\n" );
		}
	}

// write soot mass fraction
	fprintf( fp, "MassFrac_Soot\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", moments[k-1][1] * fMolarMassSoot / density[k-1] );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}

// write soot volume fraction
	fprintf( fp, "fv\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", moments[k-1][1] * fvFact );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}

// write number density
	fprintf( fp, "numdens [cm^-3]\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", moments[k-1][0] * numdensfact );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}

// write particle diameter in [nm]
	fprintf( fp, "partdiam [nm]\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", partdiamfact 
					* pow( MAX( moments[k-1][1] / moments[k-1][0], 1.0e-30 ), 1.0/3.0 ) );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}

// write soot oxidation coefficient
	fprintf( fp, "SootOxCoeff\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", GetSootOxCoeff( Y[k-1], density[k-1], molarMass ) );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}

// write soot surface growth coefficient
	fprintf( fp, "GetSurfGrowthCoeffForw\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", GetSurfGrowthCoeffForw( Y[k-1], density[k-1], molarMass ) );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}

// write soot surface growth coefficient
	fprintf( fp, "GetSurfGrowthCoeffBackw\n" );
	for ( k = 0; k < gridPoints+2; ++k ) {
		fprintf( fp, "\t%-.6e", GetSurfGrowthCoeffBackw( Y[k-1], density[k-1], molarMass ) );
		if ( (k+1) % 5 == 0 ) {
			fprintf( fp, "\n" );
		}
	}
	if ( k % 5 ) {
		fprintf( fp, "\n" );
	}

// write soot source terms
	for ( i = 0; i < fNSootMoments; ++i ) {
	// Nucleation
		if ( fNucleation ) {
			fprintf( fp, "Nucleation%d\n", i );
			for ( k = 0; k < gridPoints+2; ++k ) {
				fprintf( fp, "\t%-.6e", NucleationNew( i, temp[k-1], pahMoments[k-1] ) );
				if ( (k+1) % 5 == 0 ) {
					fprintf( fp, "\n" );
				}
			}
			if ( k % 5 ) {
				fprintf( fp, "\n" );
			}
		}
	// Condensation
		if ( fCondensation ) {
			fprintf( fp, "Condensation%d\n", i );
			for ( k = 0; k < gridPoints+2; ++k ) {
				fprintf( fp, "\t%-.6e", SourceCondensationNew( i, temp[k-1], pahMoments[k-1]
										, moments[k-1], Y[k-1], density[k-1], molarMass ) );
				if ( (k+1) % 5 == 0 ) {
					fprintf( fp, "\n" );
				}
			}
			if ( k % 5 ) {
				fprintf( fp, "\n" );
			}
		}
		// Coagulation
		if ( fCoagulation ) {
			fprintf( fp, "Coagulation%d\n", i );
			for ( k = 0; k < gridPoints+2; ++k ) {
				fprintf( fp, "\t%-.6e", SourceCoagulationNew( i, temp[k-1], moments[k-1] ) );
				if ( (k+1) % 5 == 0 ) {
					fprintf( fp, "\n" );
				}
			}
			if ( k % 5 ) {
				fprintf( fp, "\n" );
			}
		}
		// Surface Growth
		if ( fSurfaceGrowth ) {
			fprintf( fp, "SurfGrowth%d\n", i );
			for ( k = 0; k < gridPoints+2; ++k ) {
				ComputeSootRateCoeffs( kSoot, temp[k-1], flame->GetReaction() );
				ComputeCSootStar( kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1] );
				fprintf( fp, "\t%-.6e", SourceSurfGrowthNew( i, moments[k-1], Y[k-1], density[k-1], molarMass ) );
				if ( (k+1) % 5 == 0 ) {
					fprintf( fp, "\n" );
				}
			}
			if ( k % 5 ) {
				fprintf( fp, "\n" );
			}
		}
	// Oxidation
		if ( fSurfaceOxidation ) {
			fprintf( fp, "SurfOx%d\n", i );
			for ( k = 0; k < gridPoints+2; ++k ) {
				ComputeSootRateCoeffs( kSoot, temp[k-1], flame->GetReaction() );
				ComputeCSootStar( kSoot, Y[k-1], density[k-1], molarMass, mixMolarMass[k-1] );
				fprintf( fp, "\t%-.6e", SourceSootOxidationNew( i, moments[k-1], Y[k-1], density[k-1], molarMass ) );
				if ( (k+1) % 5 == 0 ) {
					fprintf( fp, "\n" );
				}
			}
			if ( k % 5 ) {
				fprintf( fp, "\n" );
			}
		}
	}
}

void T1DSoot::PrintRHSSoot( TNewtonPtr bt, T1DFlamePtr flame )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		i, k;
    int         		N = currentGrid->GetNGridPoints();
	FILE				*fp = NULL;
	Double		dummy;
	int			counter;
	
	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( flame->GetOutFileBuff(), "%srhsSoot%d.dout", flame->GetOutputPath(), counter );
	if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) { 
		cerr << "#warning: unable to open file " << flame->GetOutFileBuff() << NEWL;
		exit(2);
	}
	
	fprintf( fp, "*\n%-12s", "eta" );
	for ( i = 0; i < fNSootMoments; ++i ) {
		fprintf( fp, "\tConv_%-d\tDiff_%-d\tThPh_%-d\tNucl_%-d\tCoag_%-d\tCoaN_%-d\tCond_%-d\tSuGr_%-d\tSuOx_%-d", i, i, i, i, i, i, i, i, i );
	}
	fprintf( fp, "\n" );
		
	for ( k = 0; k < N; ++k ){
		bt->SetNodeInfo( flame, k );
		PrintRHSSoot( flame, nodeInfo, fp );
	}
    fclose( fp );
}

void T1DSoot::PrintRHSSoot( T1DFlamePtr flame, NodeInfoPtr nodeInfo, FILE *fp )
{
	int 			i;
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			temp = flameNode->temp[kCurr];
	Double			*pahMoments = flameNode->pahMoments;
	Double			*moments = flameNode->moments;
	Double			*Y = flameNode->Y[kCurr];
	Double			*kSoot = fSootRateCoeffs->vec;
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double			density = flameNode->mixDensity[kCurr];

	ComputeFractionalMoments( moments );
	ComputePhi( temp );
	ComputeSootRateCoeffs( kSoot, temp, flame->GetReaction() );
	ComputeCSootStar( kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr] );

	fprintf( fp, "%-.6e", *nodeInfo->x );

	for ( i = 0; i < fNSootMoments; ++i ) {

		fprintf( fp, "\t%-.6e", SootConvection( i, flame, nodeInfo, TRUE ) );

		if ( fSizeDepDiff ) {
			fprintf( fp, "\t%-.6e", -SootDiffusion( i, kPhysical, flame, nodeInfo ) );
		}
		else {
			fprintf( fp, "\t%-.6e", -SootDiffusionNew( i, kPhysical, flame, nodeInfo ) );
		}
		fprintf( fp, "\t%-.6e", -SootThermoPhoresis( i, kPhysical, flame, nodeInfo ) );

		fprintf( fp, "\t%-.6e", -NucleationNew( i, temp, pahMoments ) );
		fprintf( fp, "\t%-.6e", -SourceCoagulationNew( i, temp, moments ) );
		fprintf( fp, "\t%-.6e", -SourceSurfDepCoag( i, moments, Y, temp, density, molarMass ) );
		fprintf( fp, "\t%-.6e", -SourceCondensationNew( i, temp, pahMoments, moments, Y, density, molarMass ) );
		fprintf( fp, "\t%-.6e", -SourceSurfGrowthNew( i, moments, Y, density, molarMass ) );
		fprintf( fp, "\t%-.6e", -SourceSootOxidationNew( i, moments, Y, density, molarMass ) );
	}
	fprintf( fp, "\n" );
}

Double T1DSoot::FillJacCondensationNew( int l, int i, Double temp, Double *pahMom, Double *mom )
{
	switch( i ) {
		case 0: 
			return 0.0;
		case 1:
			return GetDerivPhiCond( l, 0, 1, temp, mom, pahMom );
		case 2:
			return GetDerivPhiCond( l, 0, 2, temp, mom, pahMom ) 
				+ 2.0 * GetDerivPhiCond( l, 1, 1, temp, mom, pahMom );
		case 3:
			return GetDerivPhiCond( l, 0, 3, temp, mom, pahMom ) 
				+ 3.0 * ( GetDerivPhiCond( l, 1, 2, temp, mom, pahMom ) 
						+ GetDerivPhiCond( l, 2, 1, temp, mom, pahMom ) );
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
	return 0.0;
}

void T1DSoot::FillJacSurfGrowth( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	int		iOff = i + fOffsetMoments;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	density = flameNode->mixDensity[kCurr];
	Double	**a = nodeInfo->a;
	Double	*moments = flameNode->moments;
	Double	*Y = flameNode->Y[kCurr];

	constCoeff *= GetSurfGrowthCoeff( Y, density, molarMass ) * nodeInfo->hnenn;

	switch( i ) {
		case 0: 
			break;
		case 1:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k23_24, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k23_24, moments );
			break;
		case 2:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k23_24, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k23_24, moments );

			a[fOffsetMoments+1][iOff] += 2.0 * constCoeff * GetDerivFracMom( 1, k47_24, moments );
			a[fOffsetMoments+2][iOff] += 2.0 * constCoeff * GetDerivFracMom( 2, k47_24, moments );
			break;
		case 3:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k23_24, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k23_24, moments );

			a[fOffsetMoments+1][iOff] += 3.0 * constCoeff * GetDerivFracMom( 1, k47_24, moments );
			a[fOffsetMoments+2][iOff] += 3.0 * constCoeff * GetDerivFracMom( 2, k47_24, moments );

			a[fOffsetMoments+2][iOff] += 3.0 * constCoeff * GetDerivFracMom( 2, k71_24, moments );
			a[fOffsetMoments+3][iOff] += 3.0 * constCoeff * GetDerivFracMom( 3, k71_24, moments );
			break;
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
}

void T1DSoot::FillJacSurfGrowthNew( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	int		iOff = i + fOffsetMoments;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	density = flameNode->mixDensity[kCurr];
	Double	**a = nodeInfo->a;
	Double	*moments = flameNode->moments;
	Double	*Y = flameNode->Y[kCurr];

	constCoeff *= GetSurfGrowthCoeff( Y, density, molarMass ) * nodeInfo->hnenn;

	switch( i ) {
		case 0: 
			break;
		case 1:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, 2.0/3.0, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, 2.0/3.0, moments );
/*			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k02_3, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k02_3, moments );*/
			break;
		case 2:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k02_3, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k02_3, moments );

			a[fOffsetMoments+1][iOff] += 2.0 * constCoeff * GetDerivFracMom( 1, k05_3, moments );
			a[fOffsetMoments+2][iOff] += 2.0 * constCoeff * GetDerivFracMom( 2, k05_3, moments );
			break;
		case 3:
			a[fOffsetMoments][iOff] += constCoeff * GetDerivFracMom( 0, k02_3, moments );
			a[fOffsetMoments+1][iOff] += constCoeff * GetDerivFracMom( 1, k02_3, moments );

			a[fOffsetMoments+1][iOff] += 3.0 * constCoeff * GetDerivFracMom( 1, k05_3, moments );
			a[fOffsetMoments+2][iOff] += 3.0 * constCoeff * GetDerivFracMom( 2, k05_3, moments );

			a[fOffsetMoments+2][iOff] += 3.0 * constCoeff * GetDerivFracMom( 2, k08_3, moments );
			a[fOffsetMoments+3][iOff] += 3.0 * constCoeff * GetDerivFracMom( 3, k08_3, moments );
			break;
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
}

void T1DSoot::FillJacSootRadiation( Double fact, T1DFlamePtr flame, NodeInfoPtr nodeInfo )
{
	// alphas taken from:
	//		Hubbard, G. L. , Tien C. L: 
	//		Infrared Mean Absorption Coefficient
	//		of Luminous Flames and Smoke
	//		Journal of Heat Transfer, vol 100, p. 235ff, 1978

	TFlameNodePtr	flameNode = flame->GetFlameNode();
	int				tempOff = flame->GetOffsetTemperature();
	Double			temp = flameNode->temp[kCurr];
	const Double	alpha = -3.75e5, beta = 1735.0;
	Double			*moments = flameNode->moments;
	Double			**a = nodeInfo->a;
	Double			alphas = alpha + beta * temp;			// [m^-1] 
	Double			rad = GetSootRadiation( temp, moments );
	
	fact *= rad * nodeInfo->hnenn;

	a[fOffsetMoments+1][tempOff] += fact * flameNode->mixDensity[kCurr] / moments[1];

	a[tempOff][tempOff] += fact * ( 4.0 / temp 
			+ beta / alphas );
}

void T1DSoot::FillJacSootOxidationNew( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	int		iOff = i + fOffsetMoments;
	Double	*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double	density = flameNode->mixDensity[kCurr];
	Double	**a = nodeInfo->a;
	Double	*moments = flameNode->moments;
	Double	*Y = flameNode->Y[kCurr];
	Double	beta;

	constCoeff *= GetSootOxCoeff( Y, density, molarMass ) * nodeInfo->hnenn;

	switch( i ) {
		case 0: 
#ifdef OXMOM0NEW
			beta = MIN( ( ( moments[0] > SMALLSOOT ) ? moments[0] : 0.0 ) 
						/ MAX( moments[1], SMALLSOOT ), 1.0 / 9.0 );
#else
			beta = fBeta;
#endif
			a[fOffsetMoments][iOff] -= constCoeff * beta * GetDerivFracMom( 0, km1_3, moments );
			a[fOffsetMoments+1][iOff] -= constCoeff * beta * GetDerivFracMom( 1, km1_3, moments );
			break;
		case 1:
			a[fOffsetMoments][iOff] -= constCoeff * GetDerivFracMom( 0, k02_3, moments );
			a[fOffsetMoments+1][iOff] -= constCoeff * GetDerivFracMom( 1, k02_3, moments );
			break;
		case 2:
			a[fOffsetMoments][iOff] -= constCoeff * ( 1.0 - fBeta ) 
											* GetDerivFracMom( 0, k02_3, moments );
			a[fOffsetMoments+1][iOff] -= constCoeff * ( 1.0 - fBeta ) 
											* GetDerivFracMom( 1, k02_3, moments );

			a[fOffsetMoments+1][iOff] -= 2.0 * constCoeff * GetDerivFracMom( 1, k05_3, moments );
			a[fOffsetMoments+2][iOff] -= 2.0 * constCoeff * GetDerivFracMom( 2, k05_3, moments );
			break;
		case 3:
			a[fOffsetMoments][iOff] -= constCoeff * ( 1.0 - fBeta ) 
											* GetDerivFracMom( 0, k02_3, moments );
			a[fOffsetMoments+1][iOff] -= constCoeff * ( 1.0 - fBeta ) 
											* GetDerivFracMom( 1, k02_3, moments );

			a[fOffsetMoments+1][iOff] -= 3.0 * constCoeff * ( 1.0 - fBeta )
											* GetDerivFracMom( 1, k05_3, moments );
			a[fOffsetMoments+2][iOff] -= 3.0 * constCoeff * ( 1.0 - fBeta )
											* GetDerivFracMom( 2, k05_3, moments );

			a[fOffsetMoments+2][iOff] -= 3.0 * constCoeff * GetDerivFracMom( 2, k08_3, moments );
			a[fOffsetMoments+3][iOff] -= 3.0 * constCoeff * GetDerivFracMom( 3, k08_3, moments );
			break;
		default:
			cerr << "#error: no need to compute source term for condensation for M_" << i << NEWL;
			exit( 2 );
	}
}

void T1DSoot::ComputeDiffusivity( T1DFlamePtr flame )
{
	// omega_ij = 1.0

	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			temp = flameNode->temp[kCurr];
	static Double	pi = 4.0 * atan( 1.0 );
	static Double	fact = 1.5 * sqrt( 0.5 / pi );
	static Double	dMin = pow( 6.0 / pi * fMolarMassSoot / ( fSootDensity * AVOGADRO ), 1.0 / 3.0 ); // [m]
	static Double	dMin2 = dMin * dMin; // [m^2]

	flameNode->diffSoot[kCurr] = fact / flameNode->mixDensity[kCurr] 
			* sqrt( flameNode->mixMolarMass[kCurr] * RGAS * temp ) / AVOGADRO / dMin2;
	
}

Double SootCoagFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag /*theFlag*/ )
{
	T1DFlamePtr 	flame = ( T1DFlamePtr ) object;
//	TFlameNodePtr	flameNode = flame->GetFlameNode();
	T1DSootPtr		soot = flame->GetSoot();
	int				momOff = soot->GetOffsetSootMoments();
	Double			*phi = soot->fPhi->vec;
	Double			temp = nodeInfo->y[flame->GetOffsetTemperature()];
//	Double			temp = flameNode->temp[kCurr];
	Double			*moments = &nodeInfo->y[momOff];
	
	soot->ComputeFractionalMoments( moments );
	soot->ComputePhi( temp );
	
	switch ( j - momOff ) {
		case 0:
			return -0.5 * phi[kh00];
		case 1:
			return 0.0;
		case 2:
			return phi[kh11];
		case 3:
			return 3.0 * phi[kh12];
		default:
			cerr << "#error: no need to compute coagulation source_" << j - momOff << NEWL;
			exit( 2 );
	}
	return 0.0;
}

Double SootCondFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag /*theFlag*/ )
{
	T1DFlamePtr 	flame = ( T1DFlamePtr ) object;
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	T1DSootPtr		soot = flame->GetSoot();
	int				momOff = soot->GetOffsetSootMoments();
	Double			density = flameNode->mixDensity[kCurr];
	Double			*Y = flameNode->Y[kCurr];
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double			temp = nodeInfo->y[flame->GetOffsetTemperature()];
//	Double			temp = flameNode->temp[kCurr];
	Double			*mom = &nodeInfo->y[momOff];
	Double			*pahMom = flame->GetFlameNode()->pahMoments;
		
	return soot->SourceCondensationNew( j - momOff, temp, pahMom, mom, Y, density, molarMass );
}

Double SootOxidationFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag /*theFlag*/ )
{
	T1DFlamePtr 	flame = ( T1DFlamePtr ) object;
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			source;
	Double			density = flameNode->mixDensity[kCurr];
	Double			*Y = flameNode->Y[kCurr];
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	T1DSootPtr		soot = flame->GetSoot();
	Double			temp = flameNode->temp[kCurr];
	int				momOff = soot->GetOffsetSootMoments();
	Double			*kSoot = soot->fSootRateCoeffs->vec;
	int				nSootMoments = soot->GetNSootMoments();
	Double			*moments = New1DArray( nSootMoments );

	for ( int i = 0; i < nSootMoments; ++i ) {
		moments[i] *= density;
	}

	soot->ComputeSootRateCoeffs( kSoot, temp, flame->GetReaction() );
	source = soot->SourceSootOxidationNew( j-momOff, moments, Y, density, molarMass );
	Free1DArray( moments );
	return source;
}

void T1DSoot::PrintSurfDepCoag( TNewtonPtr bt, T1DFlamePtr flame )
{
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		k;
    int         		N = currentGrid->GetNGridPoints();
	FILE				*fp = NULL;
	Double		dummy;
	int			counter;
	
	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( flame->GetOutFileBuff(), "%sSurfDepCoag%d.dout", flame->GetOutputPath(), counter );
	if ( !( fp = fopen( flame->GetOutFileBuff(), "w") ) ) { 
		cerr << "#warning: unable to open file " << flame->GetOutFileBuff() << NEWL;
		exit(2);
	}
	
	fprintf( fp, "*\n%-12s", "y" );
	fprintf( fp, "\tratio\tCoagNew/Coag\tksg\tc*phi\tc\tS\ttemp\tcoag\tCst\n" );
		
	for ( k = 0; k < N; ++k ){
		bt->SetNodeInfo( flame, k );
		PrintSurfDepCoag( flame, nodeInfo, fp );
	}
    fclose( fp );
}

void T1DSoot::PrintSurfDepCoag( T1DFlamePtr flame, NodeInfoPtr nodeInfo, FILE *fp )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	Double			temp = flameNode->temp[kCurr];
	Double			*moments = flameNode->moments;
	Double			*Y = flameNode->Y[kCurr];
	Double			*molarMass = flame->GetSpecies()->GetMolarMass()->vec;
	Double			density = flameNode->mixDensity[kCurr];
	Double			*kSoot = fSootRateCoeffs->vec;

	ComputeFractionalMoments( moments );
	ComputePhi( temp );
	ComputeSootRateCoeffs( kSoot, temp, flame->GetReaction() );
	ComputeCSootStar( kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr] );
	Double	ksg = GetSurfGrowthCoeff( Y, density, molarMass );
	Double	c = GetC( temp );
	Double	coag = -0.5 * fPhi->vec[kh00];
	Double	S = GetSCoag( moments, temp );
	Double	Cst = 1.0 / ( 6.0 * fCoagFact );
	Double	newCoag = coag * ( 1.0 - 1.0 / ( 1.0 + Cst * fPhi->vec[kh00] * ksg / ( c * c * S ) ) );
	fprintf( fp, "%-.6e", *nodeInfo->x );


	fprintf( fp, "\t%-.6e", fPhi->vec[kh00]*Cst*ksg/(c*c*S) );
	fprintf( fp, "\t%-.6e", newCoag/coag );
	fprintf( fp, "\t%-.6e", ksg );
	fprintf( fp, "\t%-.6e", fPhi->vec[kh00] );
	fprintf( fp, "\t%-.6e", c );
	fprintf( fp, "\t%-.6e", S );
	fprintf( fp, "\t%-.6e", temp);
	fprintf( fp, "\t%-.6e", coag);
	fprintf( fp, "\t%-.6e", Cst );
	fprintf( fp, "\n" );
}

#endif // ZEROD
