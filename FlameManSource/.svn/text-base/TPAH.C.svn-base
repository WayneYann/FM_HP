#include "FlameMaster.h"

#define EPS 1.0e-30

#define SIMPLESOOT

#undef PAHFROMA4
#undef NOPAH
#undef A4STEADYSTATE

#ifdef SIMPLESOOT
#define PAHFROMA4
#define NOPAH
#define A4STEADYSTATE
#endif

#undef NOCONDENSATION

#define NEWPOLY
#undef ITERATEM0
#undef SIZEDEPBETA
#undef SGINJACOBIAN
#define SGINRHS
#define DEBUGKCO
#undef DEBUG
#undef DEBUGSUM
#define DEBUGSQRT
#define UPDATEPRODRATE
#define MASSMOMENTS

void TSoot::PrintPAHReactions( TFlamePtr flame, TSpeciesPtr species )
{
	FILE	*fp = flame->GetOutfile( "pahreacs", flame->kText );

   	for ( int i = 0; i < GetNPAHReactions(); ++i ) {
		PrintPAHReactions( i, species, fp );
	}
	fclose( fp );
}

int TSoot::PrintPAHReactionEquation( int number, TSpeciesPtr species, FILE *fp )
{
	int 			i;
	Flag			first = TRUE;
	Double			*nu = fNu[number]->vec;
	int				*speciesNumber = fSpeciesNumber[number]->vec;
	int				nSpecPerReac = fSpeciesNumber[number]->len;
	char			**names = species->GetNames();
	int				len = 0;

# if defined (applec) || defined (powerc)
	SpinCursor( 1 );
#endif

/*  print left side of equation  */
	for ( i = 0; i < nSpecPerReac; ++i ) {
		if ( nu[i] > 0.0 ) {
			if ( first ) {
				if ( fabs( nu[i] ) != 1.0 ) {
 					len += fprintf( fp, "%g %s", fabs( nu[i] ), names[speciesNumber[i]] );
				}
				else {
					len += fprintf( fp, "%s", names[speciesNumber[i]] );
				}
				first = FALSE;
			}
			else {
				if ( fabs( nu[i] ) != 1.0 ) {
					len += fprintf( fp, " + %g %s", fabs( nu[i] ), names[speciesNumber[i]] );
				}
				else {
					len += fprintf( fp, " + %s", names[speciesNumber[i]] );
				}
			}
		}
	}

/*  print assignment operator  */
	len += fprintf( fp, " = " );

/*  print right side of equation  */
	first = TRUE;
	for ( i = 0; i < nSpecPerReac; ++i ) {
		if ( nu[i] < 0.0 ) {
			if ( first ) {
				if ( fabs( nu[i] ) != 1.0 ) {
					len += fprintf( fp, "%g %s", fabs( nu[i] ), names[speciesNumber[i]] );
				}
				else {
					len += fprintf( fp, "%s", names[speciesNumber[i]] );
				}
				first = FALSE;
			}
			else {
				if ( fabs( nu[i] ) != 1.0 ) {
					len += fprintf( fp, " + %g %s", fabs( nu[i] ), names[speciesNumber[i]] );
				}
				else {
					len += fprintf( fp, " + %s", names[speciesNumber[i]] );
				}
			}
		}
	}

	return len;
}

void TSoot::PrintPAHReactions( int number, TSpeciesPtr species, FILE *fp )
{
	fprintf( fp, "%s is no. %d\n", fLabels[number], number );

	PrintPAHReactionEquation( number, species, fp );
	fprintf( fp, "\n" );
	fprintf( fp, "A = %g\n", fA->vec[number] );
	fprintf( fp, "n = %g\n", fN->vec[number] );
	fprintf( fp, "E = %g\n", fEOverRgas->vec[number] * RGAS );

	fprintf( fp, "\n\n");
}

void TSoot::CheckPAHReactions( void )
{
	int 	error = -1;
	int		offset = strlen( fPAHSymbol );
	
	if ( strcmp( fLabels[k0f] + offset, "0f" ) ) error = k0f; 
	if ( strcmp( fLabels[k0b] + offset, "0b" ) ) error = k0b; 

	if ( strcmp( fLabels[k1f] + offset, "1f" ) ) error = k1f; 
	if ( strcmp( fLabels[k1b] + offset, "1b" ) ) error = k1b; 

	if ( strcmp( fLabels[k2f] + offset, "2f" ) ) error = k2f; 
	if ( strcmp( fLabels[k2b] + offset, "2b" ) ) error = k2b; 

	if ( strcmp( fLabels[k3f] + offset, "3f" ) ) error = k3f; 
	if ( strcmp( fLabels[k3b] + offset, "3b" ) ) error = k3b; 

	if ( strcmp( fLabels[k4f] + offset, "4f" ) ) error = k4f; 
	if ( strcmp( fLabels[k4b] + offset, "4b" ) ) error = k4b; 

	if ( strcmp( fLabels[k5f] + offset, "5f" ) ) error = k5f; 
	if ( strcmp( fLabels[k5b] + offset, "5b" ) ) error = k5b; 

	if ( strcmp( fLabels[k6f] + offset, "6f" ) ) error = k6f; 
	if ( strcmp( fLabels[k6b] + offset, "6b" ) ) error = k6b; 

	if ( strcmp( fLabels[k7f] + offset, "7f" ) ) error = k7f; 
	if ( strcmp( fLabels[k7b] + offset, "7b" ) ) error = k7b; 

	if ( strcmp( fLabels[k8f] + offset, "8f" ) ) error = k8f; 
	if ( strcmp( fLabels[k8b] + offset, "8b" ) ) error = k8b; 

	if ( strcmp( fLabels[k9f] + offset, "9f" ) ) error = k9f; 
	if ( strcmp( fLabels[k9b] + offset, "9b" ) ) error = k9b; 

	if ( strcmp( fLabels[k10f] + offset, "10f" ) ) error = k10f; 
	if ( strcmp( fLabels[k10b] + offset, "10b" ) ) error = k10b; 

	if ( strcmp( fLabels[k11f] + offset, "11f" ) ) error = k11f; 
	if ( strcmp( fLabels[k11b] + offset, "11b" ) ) error = k11b; 

	if ( strcmp( fLabels[k12f] + offset, "12f" ) ) error = k12f; 
	if ( strcmp( fLabels[k12b] + offset, "12b" ) ) error = k12b; 

	if ( strcmp( fLabels[k13f] + offset, "13f" ) ) error = k13f; 
	if ( strcmp( fLabels[k13b] + offset, "13b" ) ) error = k13b; 

	if ( strcmp( fLabels[k14f] + offset, "14f" ) ) error = k14f; 
	if ( strcmp( fLabels[k14b] + offset, "14b" ) ) error = k14b; 

	if ( strcmp( fLabels[k15f] + offset, "15f" ) ) error = k15f; 
	if ( strcmp( fLabels[k15b] + offset, "15b" ) ) error = k15b; 

	if ( error > -1 ) {
		cerr << "#error in pah reactions" << NEWL;
	}
}

void TSoot::ComputePolymereConcs( Double *Y, Double temp, Double density, Double *molarMass
			, Double **Pij, Double *sumPi, Double *pahMoments, Double *moments, TReactionPtr reaction )
{
	int				j,*pahSpeciesInd = fPAHSpeciesIndex->vec;
	Double			*k = fRateCoefficients->vec;
	Double			*K = fRedRateCoeffs->vec;

#ifdef NOPAH
#	ifdef PAHFROMA4
#		ifdef A4STEADYSTATE

	Double x = -0.5;
	Double RT = RGAS * temp;
	
	Double kPAHOf = 1.000E+10;
	Double kPAHOb = 2.509E+14 * exp( -220790000.0 / RT );
	Double kPAH1 = 3.06e+10 * sqrt( temp );
	Double kPAHOH1 = 1.300E+10 * exp( -46000000.0 / RT );
	
	Double	concH = density * Y[f_H] / molarMass[f_H]; 
	Double	concOH = density * Y[f_OH] / molarMass[f_OH];  
	Double	concA3R5ACXC18H11 = density * Y[f_A3R5AC] / molarMass[f_A3R5AC]; 
	
	Double a = kPAH1 * 2.0;
	Double b = kPAHOb * concH + kPAHOH1 * concOH;
	if ( fCondensation ) {
		b += GetPhiPAHFROMA4( x, temp, moments );
	}
	Double cc = kPAHOf * concA3R5ACXC18H11;
  
	Double concA4 = SolveQuadratic(a,b,cc);

	Y[pahSpeciesInd[1]] =  concA4 * molarMass[pahSpeciesInd[1]] / density;
	if ( fCondensation ) {
//        printf("%g\t%g\t%g\t%g\t%g\t%g\n", kPAHOb * concH + kPAHOH1 * concOH
//        		, GetPhiPAHFROMA4( x, temp, moments ), concA4, moments[0], moments[1], temp );
	}

/*	printf("%e\n",RGAS);	*/
/*	printf("%e\n",kPAHOf);*/
/*	printf("%e\n",kPAHOb);*/
/*	printf("%e\n",kPAH1);*/
/*	printf("%e\n",kPAHOH1);*/
/*	printf("%e\n",molarMass[sH]);*/
/*	printf("%e\n",molarMass[sOH]);*/
/*	printf("%e\n",molarMass[sA3R5ACXC18H11]);*/
/*	exit(1);*/

	Double	r1 = 1.0;
	Double	rand1 = 9.0;

	for ( j = 0; j < fNPAHMoments; ++j, r1 *= rand1 ) {
		pahMoments[j] = r1 * concA4;
	}
#		else
	Double	r1 = 1.0;
	Double	rand1 = 9.0;
/* ingo   */ 
 	Double	conc = density * Y[pahSpeciesInd[1]] / molarMass[pahSpeciesInd[1]]; 
//	Double	conc = density * 10e-7 / molarMass[pahSpeciesInd[1]]; 
//    if( Y[pahSpeciesInd[1]] != 0. ){
//        printf("%e\n",Y[pahSpeciesInd[1]]);
//	}
	for ( j = 0; j < fNPAHMoments; ++j, r1 *= rand1 ) {
		pahMoments[j] = r1 * conc;
	}
#		endif
#	else
	for ( i = 1; i <= fNPAHMolecules; ++i ) {
		sumPi[i-1] = density * Y[pahSpeciesInd[i]] / molarMass[pahSpeciesInd[i]];
	}

	Double	r1, r2, r3, r4;
	Double	rand1 = 9.0;
	Double	rand2 = 10.0;
	Double	rand3 = 11.0;
	Double	rand4 = 12.0;
	r1 = r2 = r3 = r4 = 1.0;
	for ( j = 0; j < fNPAHMoments; ++j, r1 *= rand1, r2 *= rand2, r3 *= rand3, r4 *= rand4 ) {
		pahMoments[j] += r1 * ( sumPi[0] + sumPi[1] )
						+ r2 * ( sumPi[2] + sumPi[3] )
						+ r3 * ( sumPi[4] + sumPi[5] )
						+ r4 * ( sumPi[6] + sumPi[7] );
	}
#	endif
#else
	ComputeRateCoefficients( k, temp, reaction );
#ifdef NEWPOLY
	if ( Y[f_A3R5M] < 1.0e-16 ) { //}
#else
	if ( Y[f_A3R5AC] < 1.0e-20 ) {
#endif
		Clear2DArray( Pij, fNPAHMolecules+1, fNStages );
		Clear1DArray( sumPi, fNPAHMolecules );
		Clear1DArray( pahMoments, fNPAHMoments );

		for ( j = 1; j <= fNPAHMolecules; ++j ) {
			Y[pahSpeciesInd[j]] = 0.0;
		}
	}
	else {
		ComputeRedRateCoeffs( K, Y, k, density, molarMass );
#ifdef ITERATEM0
		for ( int jjj = 0; jjj < 20; ++jjj ) {
#endif
		fKco = ComputeKco( temp, moments, pahMoments );
		ComputeF();
//fprintf( stderr, "\nTemp = %g\tY[f_A3R5M] = %g\n", temp, Y[f_A3R5M] );
		ComputeZ();
#ifdef NEWPOLY
		ComputeNewPoly( Y, density, molarMass, Pij[0] );
		CalcMomentsNewPoly( pahMoments, sumPi, Pij[0] );
#else
		ComputeAB( kFirstPoly );
		ComputeAB( kRestPoly );
		fXInf = ComputeXInf();
		ComputeDelta();
	
		ComputeP1j( Pij );
		ComputePij( Pij );
//	compute concentrations and sumPi
		int 	ind;
		for ( j = 1; j <= fNPAHMolecules; ++j ) {
			ind = pahSpeciesInd[j];
			sumPi[j-1] = 0.0;
			Y[ind] = 0.0;
			for ( int i = 0; i < fNStages; ++i ) {
				sumPi[j-1] += Pij[i][j];
				Y[ind] += Pij[i][j] * ( molarMass[ind] + i * fMolarMassDiff );
			}
			Y[ind] /= density;
		}
#define CHECKPAH
#ifdef CHECKPAH
		Double	sumY = 0.0;
		for ( j = 1; j <= fNPAHMolecules; ++j ) {
			sumY += Y[pahSpeciesInd[j]];
		}
		if ( sumY > 1.0 ) {
			cerr << "sumY = " << sumY
				<< TAB << "Y[A3R5AC] = " << Y[f_A3R5AC]
				<< TAB << "fKco = " << fKco
				<< TAB << "fXInf = " << fXInf << NEWL;
			cerr << "Gamma( kFirstPoly ) = " << ComputeGamma( kFirstPoly )
				<< TAB << "Gamma( kRestPoly ) = " << ComputeGamma( kRestPoly )
				<< TAB << "fDelta[1] = " << fDelta->vec[1]
				<< TAB << "P01 = " << Pij[0][1]
				<< NEWL;
/*			for ( i = 0; i < 2; ++i ) {*/
/*				fprintf( stderr, "\nStage #%d:\n", i );*/
/*				for ( j = 1; j <= fNPAHMolecules; ++j ) {*/
/*					fprintf( stderr, "P_%d = %g\n", j, Pij[i][j] );*/
/*				}*/
/*			}*/
			Clear2DArray( Pij, fNPAHMolecules+1, fNStages );
			Clear1DArray( sumPi, fNPAHMolecules );
			Clear1DArray( pahMoments, fNPAHMoments );
	
			for ( j = 1; j <= fNPAHMolecules; ++j ) {
				Y[pahSpeciesInd[j]] = 0.0;
			}
		}
#endif
#ifdef DEBUGSUM
		Double	sumP = 0.0;
		Double	sumY = 0.0;
		for ( j = 1; j <= fNPAHMolecules; ++j ) {
			ind = pahSpeciesInd[j];
			sumP += sumPi[j-1];
			sumY += Y[ind];
		}
		cerr << "sumP = " << sumP << TAB << "sumY = " << sumY << NEWL;
#endif
		ComputePAHMoments( pahMoments, Pij );
//		fprintf( stderr, "M0Model = %g\n", pahMoments[0] );
#endif
#ifdef ITERATEM0
		}
#endif
//		fprintf( stderr, "%g\n", pahMoments[0] );
	}
#endif
}


#ifdef MASSMOMENTS
void TSoot::ComputePAHMoments( Double *pahMoments, Double **Pij )
{
	int 	i, j;
	Double	rand = 9;
	Double	randPlus, randPlusPlus, randPowR, randPlusPowR, randPlusPlusPowR;

// arrayman function replaced by loop
//	Clear1DArray( pahMoments, fNPAHMoments );
	for ( i = 0; i < fNPAHMoments; ++i ) {
		pahMoments[i] = 0.0;
	}
	
	for ( i = 0; i < fNStages; ++i, rand += 3.0 ) {
		randPlus = rand + 1.0;
		randPlusPlus = randPlus + 1.0;
		randPowR = 1.0;
		randPlusPowR = 1.0;
		randPlusPlusPowR = 1.0;
		for ( j = 0; j < fNPAHMoments; ++j, randPowR *= rand, randPlusPowR *= randPlus, randPlusPlusPowR *= randPlusPlus ) {
			pahMoments[j] += randPowR * ( Pij[i][1] + Pij[i][2] )
							+ randPlusPowR * ( Pij[i][3] + Pij[i][4] )
							+ randPlusPlusPowR * ( Pij[i][5] + Pij[i][6] );
		}
	}
}
#else
void TSoot::ComputePAHMoments( Double *pahMoments, Double **Pij )
{
	int 	i, j;
	Double	rand = fNMinRings;
	Double	randPlus, randPowR, randPlusPowR;

// arrayman function replaced by loop
//	Clear1DArray( pahMoments, fNPAHMoments );
	for ( i = 0; i < fNPAHMoments; ++i ) {
		pahMoments[i] = 0.0;
	}
	
	for ( i = 0; i < fNStages; ++i, rand += 2.0 ) {
		randPlus = rand + 1.0;
		randPowR = 1.0;
		randPlusPowR = 1.0;
		for ( j = 0; j < fNPAHMoments; ++j, randPowR *= rand, randPlusPowR *= randPlus ) {
			pahMoments[j] += randPowR * ( Pij[i][1] + Pij[i][2] 
									+ Pij[i][3] + Pij[i][4] )
							+ randPlusPowR * ( Pij[i][5] + Pij[i][6] );
		}
	}
}
#endif

void TSoot::ComputeRateCoefficients( Double *k, Double temp, TReactionPtr reaction )
{
	Double	nPAHReactions = GetNPAHReactions();
	Double	*a = fA->vec;
	Double	*n = fN->vec;
	Double	*eOverR = fEOverRgas->vec;
	
	for ( int i = 0; i < nPAHReactions; ++i ) {
		reaction->ComputeRateCoefficient( temp, k[i], a[i], n[i], eOverR[i] );
	}
}

void TSoot::ComputeRedRateCoeffs( Double *K, Double *Y, Double *k, Double density, Double *molarMass )
{
	Double	C_H = MAX( 0.0, density * Y[f_H] / molarMass[f_H] );
	Double	C_H2 = MAX( 0.0, density * Y[f_H2] / molarMass[f_H2] );
	Double	C_O2 = MAX( 0.0, density * Y[f_O2] / molarMass[f_O2] );
	Double	C_OH = MAX( 0.0, density * Y[f_OH] / molarMass[f_OH] );
	Double	C_H2O = MAX( 0.0, density * Y[f_H2O] / molarMass[f_H2O] );
	Double	C_C2H = MAX( 0.0, density * Y[f_C2H] / molarMass[f_C2H] );
	Double	C_C2H2 = MAX( 0.0, density * Y[f_C2H2] / molarMass[f_C2H2] );
#ifdef NEWPOLY
	Double	C_A3R5M = MAX( 0.0, density * Y[f_A3R5M] / molarMass[f_A3R5M] );

	K[k01] = k[k0f] * C_A3R5M * C_C2H2;
#else
	Double	C_A3R5AC = MAX( 0.0, density * Y[f_A3R5AC] / molarMass[f_A3R5AC] );

	K[k01] = k[k0f] * C_A3R5AC;
#endif
	K[k10] = k[k0b]*C_H;

	K[k12] = k[k1f] + k[k2f]*C_H + k[k3f]*C_OH + k[k4f]*C_C2H;
	K[k21] = k[k1b]*C_H + k[k2b]*C_H2 + k[k3b]*C_H2O + k[k4b]*C_C2H2;

	K[k23] = k[k5f]*C_C2H2;
	K[k32] = k[k5b]*C_H;

	K[k34] = k[k6f] + k[k7f]*C_H + k[k8f]*C_OH + k[k9f]*C_C2H;
	K[k43] = k[k6b]*C_H + k[k7b]*C_H2 + k[k8b]*C_H2O + k[k9b]*C_C2H2;

	K[k45] = k[k10f]*C_C2H2;
	K[k54] = k[k10b];

	K[k56] = k[k11f]*C_H + k[k12f]*C_H2 + k[k13f]*C_H2O + k[k14f]*C_C2H2;
	K[k65] = k[k11b] + k[k12b]*C_H + k[k13b]*C_OH + k[k14b]*C_C2H;

	K[k51] = k[k15f]*C_C2H2;
	K[k15] = k[k15b]*C_H;

	if ( fOHPAHOxidation ) {
		K[k16] = k[k1OH]*C_OH;
		K[k31] = k[k2OH]*C_OH;
		K[k63] = k[k3OH]*C_OH;
	}
	else {
		K[k16] = K[k31] = K[k63] = 0.0;
	}

	if ( fO2PAHOxidation ) {
		K[k25] = k[k1O2]*C_O2;
		K[k42] = k[k2O2]*C_O2;
		K[k53] = k[k3O2]*C_O2;
	}
	else {
		K[k25] = K[k42] = K[k53] = 0.0;
	}
}

void TSoot::ComputeF( void )
{
	Double	*K = fRedRateCoeffs->vec;
	Double	*F = fF->vec;
	
/*#ifdef NEWPOLY*/
	Double	*nom = fDenom->vec;
	nom[1] = K[k12] + K[k15] + K[k16] + fKco + EPS;
	nom[2] = K[k21] + K[k23] + K[k25] + fKco + EPS;
	nom[3] = K[k32] + K[k34] + K[k31] + fKco + EPS;
	nom[4] = K[k45] + K[k43] + K[k42] + fKco + EPS;
	nom[5] = K[k54] + K[k51] + K[k56] + K[k53] + fKco + EPS;
	nom[6] = K[k65] + K[k63] + fKco + EPS;

	F[k12] = K[k12] / nom[1];
	F[k21] = K[k21] / nom[2];
	F[k23] = K[k23] / nom[2];
	F[k32] = K[k32] / nom[3];
	F[k34] = K[k34] / nom[3];
	F[k43] = K[k43] / nom[4];
	F[k45] = K[k45] / nom[4];
	F[k54] = K[k54] / nom[5];
	F[k15] = K[k15] / nom[1];
	F[k51] = K[k51] / nom[5];
	F[k56] = K[k56] / nom[5];
	F[k65] = K[k65] / nom[6];

/*#else*/

	fF121 = K[k12] / ( K[k12] + K[k10] + K[k16] + fKco + EPS );
	fF151 = K[k15] / ( K[k12] + K[k10] + K[k16] + fKco + EPS );

/*	F[k12] = K[k12] / ( K[k12] + K[k15] + K[k16] + fKco + EPS );*/
/*	F[k21] = K[k21] / ( K[k21] + K[k23] + K[k25] + fKco + EPS );*/
/**/
/*	F[k23] = K[k23] / ( K[k21] + K[k23] + K[k25] + fKco + EPS );*/
/*	F[k32] = K[k32] / ( K[k32] + K[k34] + K[k31] + fKco + EPS );*/
/**/
/*	F[k34] = K[k34] / ( K[k32] + K[k34] + K[k31] + fKco + EPS );*/
/*	F[k43] = K[k43] / ( K[k43] + K[k45] + K[k42] + fKco + EPS );*/
/**/
/*	F[k45] = K[k45] / ( K[k43] + K[k45] + K[k42] + fKco + EPS );*/
/*	F[k54] = K[k54] / ( K[k54] + K[k56] + K[k51] + K[k53] + fKco + EPS );*/
/**/
/**/
/*	F[k15] = K[k15] / ( K[k12] + K[k15] + K[k16] + fKco + EPS );*/
/*	F[k51] = K[k51] / ( K[k51] + K[k54] + K[k56] + K[k53] + fKco + EPS );*/
/**/
/*	F[k56] = K[k56] / ( K[k51] + K[k54] + K[k56] + K[k16] + fKco + EPS );*/
/*	F[k65] = K[k65] / ( K[k65] + K[k63] + fKco + EPS );*/
/**/
/*#endif*/

	if ( fOHPAHOxidation ) {
/*#ifdef NEWPOLY*/
		F[k31] = K[k31]/nom[3];
		F[k63] = K[k63]/nom[6];
		F[k16] = K[k16]/nom[1];
/*#else*/
/*		F[k31] = K[k31] / ( K[k32] + K[k34] + K[k31] + fKco + EPS );*/
/*		F[k16] = K[k16] / ( K[k12] + K[k15] + K[k16] + fKco + EPS );*/
/*		F[k63] = K[k63] / ( K[k65] + K[k63] + fKco + EPS );*/
/*#endif*/
	}
	else {
		F[k16] = F[k31] = F[k63] = 0.0;
	}

	if ( fO2PAHOxidation ) {
/*#ifdef NEWPOLY*/
		F[k25] = K[k25]/nom[2];
		F[k53] = K[k53]/nom[5];
		F[k42] = K[k42]/nom[4];
/*#else*/
/*		F[k25] = K[k25] / ( K[k21] + K[k23] + K[k25] + fKco + EPS );*/
/*		F[k53] = K[k53] / ( K[k54] + K[k56] + K[k51] + K[k53] + fKco + EPS );*/
/*		F[k42] = K[k42] / ( K[k43] + K[k45] + K[k42] + fKco + EPS );*/
/*#endif*/
	}
	else {
		F[k25] = F[k42] = F[k53] = 0.0;
	}
}

void TSoot::ComputeZ( void )
{
	Double	*F = fF->vec;
	Double	*Z = fZ->vec;
	
	Double	F121 = F[k12] * F[k21];
	Double	F232 = F[k23] * F[k32];
	Double	F343 = F[k34] * F[k43];
	Double	F454 = F[k45] * F[k54];
	Double	F565 = F[k56] * F[k65];
	
/*#ifdef NEWPOLY*/
	Z[6] = 1.0;
	Z[5] = 1.0;
	Z[4] = 1.0 -  F565;
	Z[3] = Z[4]   -  F454;
	Z[2] = Z[3]   -  Z[4] *F343 - F[k34]*F[k45]*(F[k53]+F[k56]*F[k63]);
	Z[1] = Z[2]   -  Z[3] *F232 - Z[4]*F[k23]*F[k34]*F[k42] ;
	Z[0] = Z[1]   -  Z[2] *F121 - Z[3]*F[k12]*F[k23]*F[k31];
	fZ53= 1.0 - F343;
	fZ52= 1.0 - (F232 + F[k23]*F[k34]*F[k42] + F343);
	fZ51= 1.0 - (fZ53*F121 + F[k12]*F[k23]*F[k31] + F232 + F[k23]*F[k34]*F[k42] + F343) + EPS;
/*#else*/
/*	Z[4] = Z[5] - F[k56] * F[k65];*/
/*	Z[3] = Z[4] - F[k45] * F[k54];*/
/*	Z[2] = Z[3] - F[k43] * F[k34] * Z[4] */
/*			- F[k34] * ( F[k45] * F[k53] + F[k45] * F[k56] * F[k63] );*/
/*	Z[1] = Z[2] - F[k23] * F[k32] * Z[3] */
/*			- F[k23] * F[k34] * F[k42] * Z[4];*/
/*#endif*/
}

void TSoot::CalcMomentsNewPoly( Double *pahMoments, Double *sumPi, Double *p0 )
{
	int		i, j;
	Double	dn = 3.0;
	Double	fxx, fyy, sum0, sum1, sum2, sum3;
	Double	n0[7];
	n0[1] = 6.0;//9.0;
	n0[2] = 6.0;//9.0;
	n0[3] = 7.0;//10.0;
	n0[4] = 7.0;//10.0;
	n0[5] = 8.0;//11.0;
	n0[6] = 8.0;//11.0;
	
	for ( i = 0; i < fNPAHMoments; ++i ) {
		pahMoments[i] = 0.0;
	}

//    calculate sums_i=0^infty gamma...
//    for sums calculated from p1 values
	fxx  = 1.0-fGamma;
	sum0 = 1.0/fxx;
	fxx  = fxx*(1.0-fGamma);
	sum1 = 1.0/fxx;
	fxx  = fxx*(1.0-fGamma);
	fyy  = 1.0*(1.0+fGamma);
	sum2 = fyy/fxx;
	fxx  = fxx*(1.0-fGamma);
	fyy  = fyy*(1.0+fGamma) + 2.0*fGamma;
	sum3 = fyy/fxx;

//     now calculate SuPij and moment(0)
	for ( j = 1; j <= fNPAHMolecules; ++j ) {
        sumPi[j-1]  = sum0*p0[j];
        pahMoments[0] += sumPi[j-1];
	}

 	for ( j = 1; j <= fNPAHMolecules; ++j ) {
		pahMoments[1] += ( n0[j]*sum0 + dn*sum1  ) *p0[j];
 	}

	if ( fNPAHMoments > 2 ) {
		for ( j = 1; j <= fNPAHMolecules; ++j ) {
			 pahMoments[2] += (  n0[j]*n0[j]          *sum0 
						 + 2.0*n0[j]*dn        *sum1
						 + dn*dn                *sum2  ) *p0[j];
		}
	}

	if ( fNPAHMoments > 3 ) {
		for ( j = 1; j <= fNPAHMolecules; ++j ) {
			 pahMoments[3] += (  n0[j]*n0[j]*n0[j]    *sum0 
						 + 3.0*n0[j]*n0[j]*dn  *sum1
						 + 3.0*n0[j]*dn*dn     *sum2      
						 + dn*dn*dn             *sum3  ) *p0[j];
		}
	}
}

Double TSoot::ComputeNewPoly( Double *Y, Double density, Double *molarMass, Double *p0 )
{
	Double	div2, fxx;
	Double	Ac[7];	
	Double	Bc[7];	
	Double	Cc[7];	
	Double	Dc;	
	Double	*F = fF->vec;
	Double	*Z = fZ->vec;
	Double	*nom = fDenom->vec;
	Double	*K = fRedRateCoeffs->vec;
	Double	F121 = F[k12] * F[k21];
	Double	F232 = F[k23] * F[k32];
	Double	F151 = F[k15] * F[k51];
	Double	F454 = F[k45] * F[k54];
	
	p0[5] = MAX( 0.0, density * Y[f_A3R5M] / molarMass[f_A3R5M] );

	Double	div  = Z[0] - fZ51*(F151+F[k16]*F[k65]*F[k51]) - fZ52*F[k12]*F[k25]*F[k51] 
               - F[k16]*F[k63]*F[k34]*F[k45]*F[k51]; // Nkor

//     calculate Ff
	Double	Ff   = F[k12]*F[k23]*F[k34]*F[k45]*F[k51];

//     calculate Fb
	Double	Fx31 = F[k32]*F[k21] + F[k31];
	Double	Fx32 = F[k32] + F[k31]*F[k12];
	Double	Fx21 = F[k21] + F[k23]*F[k31];
	Double	Fx15 = F[k15] + F[k16]*F[k65];
	Double	Fx53 = F[k53] + F[k56]*F[k63];
	Double	Fy31 = Fx31 + F[k34]*F[k42]*F[k21];
	Double	Fy32 = Fx32 + F[k34]*F[k42];

	Double	Fb   = (F[k54]*(F[k43]*Fx31 + F[k42]*Fx21) + F[k53]*Fy31)*Fx15
		+ F[k56]*F[k63]*Fy31*F[k15]
		+ F[k54]*(F[k43]*Fx32 + F[k42])*F[k25]
		+ Fx53*Fy32*F[k25]
		+ F[k63]*((1.0-F454)*Fx31 + F[k34]*F[k42]*F[k21])*F[k16];

//     calculate Fb2
//	Double	Fb2  = F[k54]*F[k42]*F[k25] * F[k63]*F[k31]*F[k16];

//     calculate gamma

	if (Fb < 1.0e-10) {
		fGamma = Ff/MAX(1e-30,div);
	}
	else {
		div2 = div*div;
		fGamma = ( div - sqrt( div2 - 4.0 * Ff * Fb ) ) / ( 2.0 * Fb );
	}
	if ( fGamma < 1.0e-40 ) {
		fGamma = Ff / MAX(1e-30,div);
	}

//fprintf( stderr, "gamma = %g\tfZ51 = %g\tK[k51] = %g\n", fGamma, fZ51, K[k51] );
//fprintf( stderr, "nom[1] = %g\tfZ52 = %g\tdiv = %g\n", nom[1], fZ52, div );
//     coefficients for p0[i], i=1,2,3,4,6
	fxx    = K[k51]/fZ51;
	Ac[1]  = fZ52*fxx/nom[1];
	fxx    = fxx*F[k12];
	Ac[2]  = fZ53*fxx/nom[2];
	fxx    = fxx*F[k23];
	Ac[3]  = fxx/nom[3];
	fxx    = fxx*F[k34];
	Ac[4]  = fxx/nom[4];
	Ac[5]  = 0.0;
	Ac[6]  = 0.0;

	Bc[1] = (   K[k54]/nom[1] *(F[k43]*Fx31 + F[k42]*Fx21) 
		   + (K[k53]+K[k56]*F[k63])/nom[1] *Fy31           )/fZ51;
	
	Bc[2] = (   K[k51]/nom[2] *F[k16]*F[k63]*(F[k32] + F[k34]*F[k42])
		   + K[k54]/nom[2] *(F[k43]*Fx32 + F[k42])
		   + (K[k53]+K[k56]*F[k63])/nom[2] *Fy32           )/fZ51;
	
	Bc[3] = (   K[k51]/nom[3] *F[k16]*F[k63]
		   + K[k54]/nom[3] *(F[k43]*(1.0-F121)+F[k42]*F[k23])
		   + (K[k53]+K[k56]*F[k63])/nom[3] *(1.0-F121)    )/fZ51;
	
	Bc[4] = (   K[k51]/nom[4] *F[k16]*F[k63]*F[k34]
		   + K[k54]/nom[4] *(1.0-F121-F[k12]*F[k23]*F[k31]-F232)
		   + (K[k53]+K[k56]*F[k63])/nom[4] *F[k34]*(1.0-F121))/fZ51;
	
	Bc[5] = 0.0;
	Bc[6] = K[k56]/nom[6]
		   + K[k51]/nom[6] *F[k16]*fZ52/fZ51;

	Cc[1]  = 0.0;
	Cc[2]  = - K[k54]/nom[2] * F[k42]*F[k63]*F[k31]*F[k16]         /fZ51;
	Cc[3]  =   K[k54]/nom[3] * F[k42]*F[k21]*F[k16]*F[k63]         /fZ51;
	Cc[4]  = - K[k54]/nom[4] * F[k63]*Fx31*F[k16]            /fZ51;
	Cc[5]  = 0.0;
	Cc[6]  = ( K[k54]/nom[6] * (F[k43]*Fx31 + F[k42]*Fx21)
		   +K[k53]/nom[6] * Fy31 )                  /fZ51;

	Dc     = F[k63]*Fy31*F[k16]/fZ51;

//     calculation of p1 instead of p0
	Double	denom = (1.0 - Dc*fGamma);
	p0[1] = (Ac[1] + fGamma*(Bc[1] + fGamma*Cc[1])) / denom * p0[5];
	p0[2] = (Ac[2] + fGamma*(Bc[2] + fGamma*Cc[2])) / denom * p0[5];
	p0[3] = (Ac[3] + fGamma*(Bc[3] + fGamma*Cc[3])) / denom * p0[5];
	p0[4] = (Ac[4] + fGamma*(Bc[4] + fGamma*Cc[4])) / denom * p0[5];
	p0[6] = (Ac[6] + fGamma*(Bc[6] + fGamma*Cc[6])) / denom * p0[5];
	p0[5] = fGamma*p0[5];
	
	fGamma = MIN( fGamma, 0.95 );

	return fGamma;
}

Double TSoot::ComputeN( int poly )
{
	Double	*F = fF->vec;
	Double	*Z = fZ->vec;
	Double	F12 = ( poly == kFirstPoly ) ? fF121: F[k12];
	
	return ( Z[1] - F12 * ( F[k21] * Z[2] + F[k23] * F[k31] * Z[3] ) );
}

void TSoot::ComputeAB( int poly )
{
	Double	*A = fACoeff->mat[poly];
	Double	*B = fBCoeff->mat[poly];
	Double	*BOH = fBOHCoeff->mat[poly];
	Double	*F = fF->vec;
	Double	*Z = fZ->vec;
	Double	*K = fRedRateCoeffs->vec;
	Double	F12 = ( poly == kFirstPoly ) ? fF121: F[k12];
	Double	F15 = ( poly == kFirstPoly ) ? fF151: F[k15];
	Double	N = fNCoeff = ComputeN( poly );
	Double	F12OverN = F12 / ( fNCoeff + EPS );
	
	A[1] = F12OverN;
	A[2] = A[1] * F[k23];
	A[3] = A[2] * F[k34];
	A[4] = A[3] * F[k45];
	A[5] = A[4] * F[k51];
	A[6] = A[5] * F[k65];
	
	A[1] *= Z[1] / ( K[k12] + EPS );
	A[2] *= Z[2] / ( K[k23] + EPS );
	A[3] *= Z[3] / ( K[k34] + EPS );
	A[4] *= Z[4] / ( K[k45] + EPS );
	A[5] *= Z[5] / ( K[k51] + EPS );
	A[6] *= Z[6] / ( K[k65] + EPS );
	
	B[1] = ( ( F[k32] * F[k21] + F[k31] ) * ( F[k43] * F[k54] + F[k53] + F[k63] * F[k56] ) 
			+ F[k21] * F[k42] * ( F[k54] + F[k34] * ( F[k53] + F[k63] * F[k56] ) ) 
			+ F[k23] * F[k31] * F[k42] * F[k54] ) * F15 / ( K[k15] * N + EPS );
	B[2] = ( ( F[k54] * F[k43] + F[k53] + F[k56] * F[k63] ) * ( F[k32] + F[k31] * F12 )
			+ F[k54] * F[k42] + ( F[k53] + F[k56] * F[k63] ) * F[k34] * F[k42] ) 
				* F[k21] / ( K[k21] * N + EPS );
	B[3] = ( ( F[k54] * F[k43] + F[k53] + F[k56] * F[k63] ) * ( 1.0 - F12 * F[k21] )
			+ F[k54] * F[k42] * F[k23] ) * F[k32] / ( K[k32] * N + EPS );
	B[4] = ( F[k54] * ( 1.0 - F12 * F[k21] - F[k23] * F[k32] - F12 * F[k23] * F[k31] )
			+ F[k34] * ( F[k53] + F[k56] * F[k63] ) * ( 1.0 - F12 * F[k21] ) ) 
				* F[k43] / ( K[k43] * N + EPS );
	B[5] = ( 1.0 - F12 * F[k21] - F[k23] * F[k32] - F12 * F[k23] * F[k31]
			- F[k43] * F[k34] * ( 1.0 - F12 * F[k21] ) - F[k23] * F[k34] * F[k42] ) 
				* F[k54] / ( K[k54] * N + EPS );
	B[6] = F[k56] * ( 1.0 - F12 * F[k21] - F[k23] * F[k32] - F12 * F[k23] * F[k31]
			- F[k43] * F[k34] * ( 1.0 - F12 * F[k21] ) - F[k23] * F[k34] * F[k42] ) 
			* F[k65] / ( K[k65] * N + EPS );

	if ( fOHPAHOxidation ) {
		BOH[1] = ( ( F[k32] * F[k21] + F[k31] ) * ( F[k63] * ( 1.0 - F[k45] * F[k54] ) 
				+ F[k65] * ( F[k43] * F[k54] + F[k53] ) )
				+ F[k63] * F[k21] * F[k42] * F[k34] 
				+ F[k65] * F[k42] * F[k21] * ( F[k54] + F[k34] * F[k53] )
				+ F[k23] * F[k31] * F[k54] * F[k65] * F[k42] ) * F15 / ( K[k15] * N + EPS );
		BOH[2] = ( ( F[k63] * ( 1.0 - F[k45] * F[k54] ) + F[k65] * ( F[k54] * F[k43] + F[k53] ) ) 
					* ( F[k32] + F[k31] * F12 )
				+ F[k65] * F[k54] * F[k42] + ( F[k65] * F[k53] + F[k63] ) * F[k34] * F[k42] ) 
					* F[k21] / ( K[k21] * N + EPS );
		BOH[3] = ( ( F[k65] * F[k54] * F[k43] + F[k63] * ( 1.0 - F[k45] * F[k54] ) ) 
					* ( 1.0 - F12 * F[k21] )
				+ F[k65] * F[k54] * F[k42] * F[k23] + F[k65] * F[k53] )
					* F[k32] / ( K[k32] * N + EPS );
		BOH[4] = ( F[k34] * ( 1.0 - F12 * F[k21] ) * ( F[k65] * F[k53] - F[k23] * F[k32] + F[k63] )
				+ F[k65] * F[k54] * ( 1.0 - F12 * F[k21] 
				- F[k23] * F[k32] - F12 * F[k23] * F[k31] ) )
					* F[k43] / ( K[k43] * N + EPS );
		BOH[5] = ( B[5] * F[k65] + F[k63] * F[k34] * F[k45] * ( 1.0 - F12 * F[k21] ) ) 
					* F[k54] / ( K[k54] * N + EPS );
		BOH[6] = ( ( 1.0 - F[k54] * F[k45] ) * ( 1.0 - F12 * F[k21] - F[k23] * F[k32] 
				- F12 * F[k23] * F[k31] ) - F[k23] * F[k34] * F[k42]
				- ( 1.0 - F12 * F[k21] ) * ( F[k45] * F[k34] * F[k53] + F[k43] * F[k34] ) ) 
				* F[k65] / ( K[k65] * N + EPS );
		for ( int i = 1; i <= fNPAHMolecules; ++i ) {
			B[i] += K[k16] / ( K[k15] + EPS ) * BOH[i];
		}
	}
	else {
		BOH[1] = BOH[2] = BOH[3] = BOH[4] = BOH[5] = BOH[6] = 0.0;
	}
}

Double TSoot::ComputeXInf( void )
{
	Double	*K = fRedRateCoeffs->vec;
	Double	*A = fACoeff->mat[kRestPoly];
	Double	*B = fBCoeff->mat[kRestPoly];
	Double	F1;
	Double	F2;
	Double	xInf;
	
	if ( fOHPAHOxidation && fO2PAHOxidation ) {
		cerr << "#warning: pah oxidation including OH && O2 reactions not yet implemented" << NEWL;
		cerr << "          switch off O2 reactions" << NEWL;
		fO2PAHOxidation = FALSE;
	}
	
	if ( fOHPAHOxidation || !( fOHPAHOxidation || fO2PAHOxidation ) ) {
		F1 = A[1] * K[k51] * B[5] * K[k15]; // K[k51] = K[12] K[k15] = K[13]
		F2 = A[5] * K[k51] * B[1] * K[k15];
		if ( F1 * F2 < 1.0e-20 ) {
			xInf = 1.0 / ( 1.0 - F1 - F2 );
		}
		else {
			Double	p = 0.5 * ( 1.0 - F1 - F2 ) / ( F1 * F2 );
			Double	rad = p * p - 1.0 / ( F1 * F2 );
		
#ifdef HP
			if ( rad < 0.0 || isnan( rad ) ) {
#else
			if ( rad < 0.0 ) {
#endif
#							ifdef FOOL_MARKERTOOL
								}
#							endif
				cerr << "#warning: invalid radiand " << rad << " in ComputeXInf" << NEWL;
				fprintf( stderr, "F1 = %g\tF2 = %g\tA[1] = %g\tA[5] = %g\tB[1] = %g\tB[5] = %g\n"
							, F1, F2, A[1], A[5], B[1], B[5] );
				rad = 0.0;
			}
			xInf = p - sqrt( rad );
		}
		xInf = ( 1.0 + F2 * xInf ) / ( 1.0 - F1 * ( 1.0 + F2 * xInf ) );
	}
	else if ( fO2PAHOxidation && !fOHPAHOxidation ) {
		F1 = A[1] * K[k51] * B[5] * K[k15]; // F0
		F2 = A[5] * K[k51] * K[k15] * ( B[1] + B[2] * K[k25] / K[k15] ); // alpha
        Double	cor = 1.0 + A[2]/A[1] * K[k25]/K[k15];
		if ( F1 * F2 < 1.0e-20 ) {
			xInf = cor / ( 1.0 - F1 * cor - F2 );
		}
		else {
			Double	p = 0.5 * ( 1.0 - F1 * cor - F2 ) / ( F1 * F2 );
			Double	rad = p * p - cor / ( F1 * F2 );
		
#			ifdef HP
			if ( rad < 0.0 || isnan( rad ) ) {
#			else
			if ( rad < 0.0 ) {
#			endif
#			ifdef FOOL_MARKERTOOL
			}
#			endif
				rad = 0.0;
				cerr << "#warning: invalid radiand in ComputeXInf" << NEWL;
			}
			xInf = p - sqrt( rad );
		}
		xInf = ( cor + F2 * xInf ) / ( 1.0 - F1 * ( cor + F2 * xInf ) );
	}
	else {
		cerr << "#warning: something wrong in 'TSoot::ComputeXInf'" << NEWL;
		cerr << "#         PAHOHOxidation is " << fOHPAHOxidation << NEWL;
		cerr << "#         PAHO2Oxidation is " << fO2PAHOxidation << NEWL;
	}
	
	return xInf;
}

void TSoot::ComputeDelta( void )
{
	Double	*K = fRedRateCoeffs->vec;
	Double	*A = fACoeff->mat[kRestPoly];
	Double	*B = fBCoeff->mat[kRestPoly];
	Double	*delta = fDelta->vec;

	fGamma = ComputeGamma( kRestPoly );

	if ( fO2PAHOxidation ) {
		for ( int j = 1; j <= fNPAHMolecules; ++j ) {
			delta[j] = K[k51] * ( A[j] + fGamma * A[1] * B[j] * K[k15]
				* ( ( 1.0 + A[5] * K[k51] * B[1] * K[k15] * fXInf ) 
					+ K[k25]/K[k15] * A[2]/A[1] 
					* ( 1.0 + A[1]/A[2] * A[5] * K[k51] * B[2] * K[k15] * fXInf ) ) );
		}
	}
	else {
		for ( int j = 1; j <= fNPAHMolecules; ++j ) {
			delta[j] = K[k51] * ( A[j] + fGamma * A[1] * B[j] * K[k15]
								* ( 1.0 + A[5] * K[k51] * B[1] * K[k15] * fXInf ) );
		}
	}
}

Double TSoot::ComputeGamma( int poly )
{
	Double	*A = fACoeff->mat[kRestPoly];
	Double	*B = fBCoeff->mat[kRestPoly];
	Double	*K = fRedRateCoeffs->vec;
	Double	K51 = ( poly == kFirstPoly ) ? K[k01]: K[k51];
	Double	A5 = ( poly == kFirstPoly ) ? fACoeff->mat[kFirstPoly][5]: A[5];
	Double	B5 = ( poly == kFirstPoly ) ? fBCoeff->mat[kFirstPoly][5]: B[5];

	if ( fO2PAHOxidation ) {
		return A5 * K51 / ( 1.0 - A[1] * K[k51] * B5 * K[k15] * ( ( 1.0 + A[2] / A[1] * K[k25] / K[k15] )
							+ ( B[1] + B[2] * K[k25] / K[k15] ) * K[k15] * A[5] * K[k51] * fXInf ) );
	}
	else {
		return A5 * K51 / ( 1.0 - A[1] * K[k51] * B5 * K[k15] 
								* ( 1.0 + A[5] * K[k51] * B[1] * K[k15] * fXInf ) );
	}
}

Double TSoot::ComputeKco( Double temp, Double *moments, Double *pahMoments/*, Double *M0*/ )
{
	Double	r0 = fRedRateCoeffs->vec[k01];
	
#ifdef DEBUGSQRT
	if ( temp < 0.0 ) {
		fprintf( stderr, "#error in function TSoot::ComputeKco\n" );
		cerr << "temp = " << temp << TAB << "Set temp to 200 K" << NEWL;
		fprintf( stderr, "M0 = %g\tM1 = %g\n", moments[0], moments[1] );
		fprintf( stderr, "MPAH0 = %g\tMPAH1 = %g\n", pahMoments[0], pahMoments[1] );
		temp = 200.0;
	}
#endif

#ifdef NOCONDENSATION
		return sqrt( r0 * GetBeta( temp, pahMoments ) );
#else
	if ( fCondensation ) {
		Double	M0[1];
		Double	betaN = GetBeta( temp, pahMoments );
		Double	betaC = GetBetaCond( temp, moments );
		Double	p = 0.5 * moments[0]/*1e-7*/ * betaC / betaN;
		Double	q = -r0 / betaN;
		
//		M0 = SolveQuadratic( betaN, moments[0] * betaC, r0 );
		*M0 = -p + sqrt( p * p - q );
//		fprintf( stderr, "M0est = %g\t", *M0 );
#	ifdef DEBUGKCO
		if ( *M0 < 0.0 ) {
			cerr << "#error in function ComputeKco" << NEWL;
		}
#	endif
//		fprintf( stderr, "M0 = %g\t", M0 );
#ifdef ITERATEM0
		return betaN * pahMoments[0] + betaC * moments[0];
#else
		return betaN * *M0 + betaC * moments[0];
#endif
	}
	else {
//		fprintf( stderr, "T = %g     M0 = %g\n", temp, sqrt( r0 / GetBeta( temp, pahMoments ) ) );
		return sqrt( r0 * GetBeta( temp, pahMoments ) );
	}
#endif
}

void TSoot::ComputeP1j( Double **Pij )
{
	Double	P21 = 0.0;
	Double	fact;
	Double	*A = fACoeff->mat[kFirstPoly];
	Double	*B = fBCoeff->mat[kFirstPoly];
	Double	*K = fRedRateCoeffs->vec;
	Double	r0 = K[k01];
	
	Pij[0][5] = MAX( ComputeGamma( kFirstPoly ), 0.0 );		// eq. 3.28 in pk
	P21 = fDelta->vec[1] * Pij[0][5];			// eq. 3.30 in pk

	if ( fO2PAHOxidation ) {
		Double	P22 = fDelta->vec[2] * Pij[0][5]; // eq. 3.30 in pk
		fact = K[k15] * P21 + K[k25] * P22;
		for ( int j = 1; j <= fNPAHMolecules; ++j ) {
			Pij[0][j] = A[j] * r0 + B[j] * fact;	// eq. 3.20 in pk
		}
	}
	else {
		fact = K[k15] * P21;
		for ( int j = 1; j <= fNPAHMolecules; ++j ) {
			Pij[0][j] = A[j] * r0 + B[j] * fact;	// eq. 3.20 in pk
		}
	}
}

void TSoot::ComputePij( Double **Pij )
{
	Double	P15 = Pij[0][5];
	Double	powGammaI = 1.0;
	Double	*delta = fDelta->vec;
	
	for ( int i = 1; i < fNStages; ++i ) {
		for ( int j = 1; j <= fNPAHMolecules; ++j ) {
			Pij[i][j] = delta[j] * powGammaI * P15;
		}
		powGammaI *= fGamma;
	}
}

Double CheckPow( Double base, Double ex );
Double CheckPow( Double base, Double ex )
{
#ifdef HP
	if ( isnan( base ) || isnan( ex ) ) {
		cerr << "#pow error: base = " << base << TAB << "exp = " << ex << NEWL;
		exit(2);
	}
#endif
	if ( base <= 0.0 || ex <= 0.0) {
		cerr << "#pow error: base = " << base << TAB << "exp = " << ex << NEWL;
		exit(2);
	}
	return pow( base, ex );
}

void TSoot::UpdateProductionRates( TSpeciesPtr species, TReactionPtr reaction
			, Double *prodRate, Double density, Double *Y, Double temp
			, Double *sumPi, Double *moments, Double *wPAH, Double mixMolarMass )
{

	if ( fSootUpdateProdRates ) {
#ifdef NOPAH
#	ifdef PAHFROMA4
#		ifdef A4STEADYSTATE
			int		*pahSpeciesInd = fPAHSpeciesIndex->vec;
			int		indHCCO = species->FindSpecies( "HCCO" );
			if ( indHCCO == -1 ) {
				cerr << "error: can't find the molecule HCCO" << NEWL;
				exit( 2 );
			}
			int		indA3R5 = species->FindSpecies( "A3R5-C16H10" );
			if ( indHCCO == -1 ) {
				cerr << "error: can't find the molecule A3R5-C16H10" << NEWL;
				exit( 2 );
			}
			Double	*molarMass = species->GetMolarMass()->vec;
			Double	RT = RGAS * temp;
			Double	kPAHOf = 1.000E+10;
			Double	kPAHOb = 2.509E+14 * exp( -220790000.0 / RT );
			Double	kPAHOH1 = 1.300E+10 * exp( -46000000.0 / RT );
			Double	concH = density * Y[f_H] / molarMass[f_H]; 
			Double	concOH = density * Y[f_OH] / molarMass[f_OH];  
			Double	concA3R5AC = density * Y[f_A3R5AC] / molarMass[f_A3R5AC]; 
			Double	concA4 = density * Y[pahSpeciesInd[1]] / molarMass[pahSpeciesInd[1]]; 
			Double	wPAHOf = kPAHOf * concA3R5AC;
			Double	wPAHOb = kPAHOb * concA4 * concH;
			Double	wPAHOH1 = kPAHOH1 * concA4 * concOH;

			prodRate[f_A3R5AC] += molarMass[f_A3R5AC] * ( wPAHOb - wPAHOf );
			prodRate[f_H] -= molarMass[f_H] * ( wPAHOb - wPAHOf );
			prodRate[f_OH] -= molarMass[f_OH] * wPAHOH1;
			prodRate[indHCCO] += molarMass[indHCCO] * wPAHOH1;
			prodRate[indA3R5] += molarMass[indA3R5] * wPAHOH1;
#		else
			int	    *pahSpeciesInd = fPAHSpeciesIndex->vec;
			Double	*molarMass = species->GetMolarMass()->vec;
			Double  x = - 1.0 / 2.0 ;
			Double	conc = density * Y[pahSpeciesInd[1]] / molarMass[pahSpeciesInd[1]]; 
		
			Double  dA4dt = conc *GetPhiPAHFROMA4( x, temp, moments ) ;	
				
//			prodRate[pahSpeciesInd[1]] -= molarMass[pahSpeciesInd[1]] * dA4dt; 
#		endif 
#	endif 
/*  ingo */
#else
		int		i, j;
		int		nSpecies = species->GetNOfSpecies();
		int 	nPAHReactions = GetNPAHReactions();
		int		speciesPerReaction;
		int		thisIndex;
		int		*pahSpeciesInd = fPAHSpeciesIndex->vec;
		Double	*molarMass = species->GetMolarMass()->vec;
		Double	*concs = New1DArray( nSpecies );
		Double	*k = fRateCoefficients->vec;
		int		*speciesIndex;
		Double	*nu;
		
	// set concentrations
		for ( i = 0; i < nSpecies; ++i ) {
			concs[i] = density * Y[i] / molarMass[i];
		}
		for ( j = 1; j <= fNPAHMolecules; ++j ) {
			concs[pahSpeciesInd[j]] = sumPi[j-1];
		}
		
	// compute reactionrates
		for ( i = 0; i < nPAHReactions; ++i ) {
			speciesIndex = fSpeciesNumber[i]->vec;
			nu = fNu[i]->vec;
			speciesPerReaction = fSpeciesNumber[i]->len;
			wPAH[i] = k[i];
			for ( j = 0; j < speciesPerReaction; ++j ) {
				if ( nu[j] > 0 ) {
					if ( nu[j] == 1.0 ) {
						wPAH[i] *= concs[speciesIndex[j]];
					}
					else {
						wPAH[i] *= CheckPow( concs[speciesIndex[j]], nu[j] );
					}
				}
			}
			for ( j = 0; j < speciesPerReaction; ++j ) {
				thisIndex = speciesIndex[j];
				if ( !species->IsSteadyState( thisIndex ) ) {
						prodRate[thisIndex] -= molarMass[thisIndex] * nu[j] * wPAH[i];
				}
			}
		}	
		Free1DArray( concs );
#endif
		if ( fSurfaceGrowth || fSurfaceOxidation ) {
			int 	nSootReactions = GetNSootReactions();
			Double	*kSoot = fSootRateCoeffs->vec;
			Double	*wSoot = fSootReactionRate->vec;
			int		*speciesIndex;
			Double	*nu;
			int		thisIndex;
			int		speciesPerReaction;
			
			ComputeSootRateCoeffs( kSoot, temp, reaction );
			ComputeCSootStar( kSoot, Y, density, molarMass, mixMolarMass );
			ComputeSootReactionRates( wSoot, Y, moments, density, molarMass, mixMolarMass );
			for ( int i = 0; i < nSootReactions; ++i ) {
				speciesIndex = fSpecNumSoot[i]->vec;
				nu = fNuSoot[i]->vec;
				speciesPerReaction = fSpecNumSoot[i]->len;
				
				for ( int j = 0; j < speciesPerReaction; ++j ) {
					thisIndex = speciesIndex[j];
					if ( !species->IsSteadyState( thisIndex ) ) {
						prodRate[thisIndex] -= molarMass[thisIndex] * nu[j] * wSoot[i];
					}
				}
			}
		}
	}
	else {
		return;
	}
}

Double TSoot::GetBeta( Double temp, Double *pahMoments )
{
	//  units of beta are [m^3 / ( kmole s )]

//	const Double	iDouble = 4.0, jDouble = 4.0;
	Double			iDouble, jDouble;
#ifdef SIZEDEPBETA
	if ( pahMoments[1] < 1.0e-30 || pahMoments[0] < 1.0e-30 ) {
#	ifdef MASSMOMENTS
		jDouble = iDouble = 9.0;
#	else
		jDouble = iDouble = 4.0;
#	endif
	}
	else {
		jDouble = iDouble = pahMoments[1] / pahMoments[0];
	}
#else
#	ifdef MASSMOMENTS
		jDouble = iDouble = 9.0;
#	else
		jDouble = iDouble = 4.0;
#	endif
#endif
	const Double	iPowThird = pow( iDouble, 1.0 / 3.0 );
	const Double	jPowThird = pow( jDouble, 1.0 / 3.0 );
	const Double	bRed = sqrt( ( iDouble + jDouble ) / ( iDouble * jDouble ) )
						* ( iPowThird + jPowThird ) * ( iPowThird + jPowThird );

//	cerr << "bRed = " << bRed << NEWL; // O(10^11)
//	cerr << "beta = " << bRed * GetC( temp ) << NEWL; // O(10^11)

	return bRed * GetC( temp );
}

Double TSoot::GetC( Double temp )
{
	//  units of C are [m^3 / ( kmole s )]

#ifdef MASSMOMENTS
	const Double	M1 = 24.0; // molar mass of smallest unit C2
#else
	const Double	M1 = 78.0; // molar mass of Benzene
#endif
	const Double	pi = 4.0 * atan( 1.0 );

// 2.2 is van der Vaals enhancement factor
	const Double	CRed = 2.2 * sqrt( 8.0 * pi * RGAS / M1 )
					* pow( 0.75 * M1 / ( AVOGADRO * pi * fSootDensity ), 2.0 / 3.0 )
					* AVOGADRO;
//	cerr << "CRed = " << CRed << NEWL; // O(10^9)
//	cerr << "C = " << CRed * sqrt( temp ) << NEWL; // O(10^9)

#ifdef DEBUGSQRT
	if ( temp < 0.0 ) {
		fprintf( stderr, "#error in function TSoot::GetC\n" );
		cerr << "temp = " << temp << TAB << "Set temp to 200 K" << NEWL;
		temp = 200.0;
	}
#endif

	return CRed * sqrt( temp );
}

Double TSoot::GetCPAH( Double temp ) //added bu GB
{
	//  units of C are [m^3 / ( kmole s )]

#ifdef MASSMOMENTS
	const Double	M1 = 24.0; // molar mass of smallest unit C2
#else
	const Double	M1 = 78.0; // molar mass of Benzene
#endif
	const Double	pi = 4.0 * atan( 1.0 );

	const Double    fPAHDensity = fSootDensity;

// 2.2 is van der Vaals enhancement factor
	const Double	CRed = 2.2 * sqrt( 8.0 * pi * RGAS / M1 )
					* pow( 0.75 * M1 / ( AVOGADRO * pi * fPAHDensity ), 2.0 / 3.0 )
					* AVOGADRO;
//	cerr << "CRed = " << CRed << NEWL; // O(10^9)
//	cerr << "C = " << CRed * sqrt( temp ) << NEWL; // O(10^9)

#ifdef DEBUGSQRT
	if ( temp < 0.0 ) {
		fprintf( stderr, "#error in function TSoot::GetC\n" );
		cerr << "temp = " << temp << TAB << "Set temp to 200 K" << NEWL;
		temp = 200.0;
	}
#endif

	return CRed * sqrt( temp );
}

Double TSoot::GetBetaCond( Double temp, Double *moments )
{
	//  units of betacond are [m^3 / ( kmole s )]

#ifdef MASSMOMENTS
	const Double	iDouble = 9;
#else
	const Double	iDouble = fNMinRings;
#endif
	const Double	iPowThird = pow( iDouble, 1.0 / 3.0 );
	Double			jDouble;
	if ( moments[0] ) {
		jDouble = moments[1] / moments[0];
	}
	else {
		jDouble = iDouble;
	}
#ifdef MASSMOMENTS
	jDouble = MAX( jDouble, 9 );
#else
	jDouble = MAX( jDouble, fNMinRings );
#endif
	Double			jPowThird = pow( jDouble, 1.0 / 3.0 );
	Double			bRed = sqrt( ( iDouble + jDouble ) / ( iDouble * jDouble ) )
						* ( iPowThird + jPowThird ) * ( iPowThird + jPowThird );

//	cerr << "beta = " << bRed * GetC( temp ) << NEWL; // O(10^9)

	return bRed * GetC( temp );
}

#ifndef ZEROD
void T1DSoot::PrintPathsOfAcetylene( T1DFlamePtr flame )
{
	TFlameNodePtr		flameNode = flame->GetFlameNode();
	TNewtonPtr 			bt = flame->GetSolver()->bt;
	TAdaptiveGridPtr	grid = bt->GetGrid();
    TGridPtr			currentGrid = grid->GetCurrentGrid();
	NodeInfoPtr			nodeInfo = bt->GetNodeInfo();
    int         		k;
    int         		N = currentGrid->GetNGridPoints();
	FILE				*fp = flame->GetOutfile( "C2H2Paths", TFlame::kData );
	Double				*molarMass = flame->GetSpecies()->GetMolarMass()->vec;

	Double	*Y = flameNode->Y[kCurr];
	Double	temp = flameNode->temp[kCurr];	
	ComputePolymereConcs( Y, temp, flameNode->mixDensity[kCurr], molarMass, flameNode->Pij
				, flameNode->sumPi, flameNode->pahMoments, flameNode->moments, flame->GetReaction() );
	fprintf( fp, "*\n%-12s", "x" );
	fprintf( fp, "\t%-30s\t%-30s\t%-30s\t%-30s\n"
			, "m\\dC2H2->Soot\\n [kg/m\\u3\\ns]"
			, "m\\dC2H2->PAH\\n [kg/m\\u3\\ns]"
			, "m\\dA3R5AC->PAH\\n [kg/m\\u3\\ns]"
			, "m\\dPAH->A3R5AC\\n [kg/m\\u3\\ns]" );
		
	for ( k = 0; k < N; ++k ){
		bt->SetNodeInfo( flame, k );
		PrintPathsOfAcetylene( flame, nodeInfo->x[kCurr], fp );
	}
    fclose( fp );
}

void T1DSoot::PrintPathsOfAcetylene( T1DFlamePtr flame, Double x, FILE *fp )
{
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	TSpeciesPtr		species = flame->GetSpecies();
	TReactionPtr	reaction = flame->GetReaction();
	int j;
	int 	nPAHReactions = GetNPAHReactions();
	int 	nSootReactions = GetNSootReactions();
	int		nSteadyStates = species->GetNSteadyStates();
	int		speciesPerReaction;
	int		thisIndex;
	Double	conc;
	Double	mDot;
	Double	prodPAH = 0.0;
	Double	prodSoot = 0.0;
	Double	density = flameNode->mixDensity[kCurr];
	Double	*Y = flameNode->Y[kCurr];
	Double	*k = fRateCoefficients->vec;
	int		*speciesIndex;
	Double	*nu;
	Double	*sumPi = flameNode->sumPi;
	Double	*molarMass = species->GetMolarMass()->vec;
	Double	*w = flameNode->pahReactionRate;
	
	for ( int i = 0; i < nPAHReactions; ++i ) {
		speciesIndex = fSpeciesNumber[i]->vec;
		nu = fNu[i]->vec;
		speciesPerReaction = fSpeciesNumber[i]->len;
		
		// first compute reactionRates
		w[i] = k[i];
		for ( j = 0; j < speciesPerReaction; ++j ) {
			if ( nu[j] > 0 ) {
				thisIndex = speciesIndex[j];
				if ( species->IsSteadyState( thisIndex ) ) {
					conc = sumPi[thisIndex - fFirstPAH];
				} 
				else {
					conc = density * Y[thisIndex] / molarMass[thisIndex];
				}
				if ( nu[j] == 1.0 ) {
					w[i] *= conc;
				}
				else {
					w[i] *= CheckPow( conc, nu[j] );
				}
			}
		}
		for ( j = 0; j < speciesPerReaction; ++j ) {
			thisIndex = speciesIndex[j];
			if ( thisIndex == f_C2H2 ) {
				mDot = nu[j] * w[i];
				prodPAH += molarMass[f_C2H2] * mDot;
			}
		}
	}

	Double	*kSoot = fSootRateCoeffs->vec;
	Double	*wSoot = fSootReactionRate->vec;
	Double	temp = flameNode->temp[kCurr];	
	
	ComputeSootRateCoeffs( kSoot, temp, reaction );
	ComputeSootReactionRates( wSoot, Y, flameNode->moments, density, molarMass, flameNode->mixMolarMass[kCurr] );
	ComputeCSootStar( kSoot, Y, density, molarMass, flameNode->mixMolarMass[kCurr] );
	prodSoot = GetSurfGrowthCoeff( Y, density, molarMass ) 
				* FracMom2( 2.0 / 3.0, flameNode->moments ) * fMolarMassSoot;

	Double	*K = fRedRateCoeffs->vec;
	ComputeRateCoefficients( k, temp, reaction );
	ComputeRedRateCoeffs( K, Y, k, density, molarMass );
	Double	m0f = K[k0f] * 226.0 /* molar mass of A4 in [kg/kmole] */;	
	Double	m0b = K[k0b] * flameNode->Pij[0][1] * 226.0 /* molar mass of A4 in [kg/kmole] */;	

	

	fprintf( fp, "%-.6e\t%-.6e\t%-.6e\t%-.6e\t%-.6e\n", x, prodSoot, prodPAH, m0f, m0b );
}

void T1DSoot::UpdateJacobian( T1DFlamePtr flame )
{
	fprintf( stderr, "T1DSoot::UpdateJacobian called\n" );
	exit(2);
#ifdef UPDATEPRODRATE
	TFlameNodePtr	flameNode = flame->GetFlameNode();
	TSpeciesPtr		species = flame->GetSpecies();
	TReactionPtr	reaction = flame->GetReaction();
	int 	nPAHReactions = GetNPAHReactions();
	int 	nSootReactions = GetNSootReactions();
	int		nSpeciesIn = species->GetNSpeciesInSystem();
	int		speciesPerReaction;
	int		thisIndex, indexL, i, j, l;
	Double	temp = flameNode->temp[kCurr];
	Double	mixMolarMass = flameNode->mixMolarMass[kCurr];
	Double	*Y = flameNode->Y[kCurr];
	Double	**dmdY = flameNode->dMdY;
	Double	*dmdT = flameNode->dMdY[nSpeciesIn];
	Double	coeffi, coeffj;
	Double	dwdT;
	int		*speciesIndex;
	Double	*nu;
	Double	*eOverR = fEOverRgas->vec;
	Double	*molarMass = species->GetMolarMass()->vec;
	Double	*w = flameNode->pahReactionRate;
	Double	sumNu;
	
	for ( i = 0; i < nPAHReactions; ++i ) {
		speciesIndex = fSpeciesNumber[i]->vec;
		nu = fNu[i]->vec;
		speciesPerReaction = fSpeciesNumber[i]->len;
		
		for ( l = 0, sumNu = 0.0; l < speciesPerReaction; ++l ) {
			if ( nu[l] > 0.0 ) {
				sumNu += nu[l];
			}
		}
		
		dwdT = w[i] / temp * ( fN->vec[i] + eOverR[i] / temp - sumNu );

		coeffi = w[i] * sumNu * mixMolarMass;
		for ( j = 0; j < speciesPerReaction; ++j ) {
			thisIndex = speciesIndex[j];
			if ( !species->IsSteadyState( thisIndex ) ) {
				coeffj = molarMass[thisIndex] * coeffi * nu[j];
				for ( l = 0; l < nSpeciesIn; ++l ) {
					dmdY[l][thisIndex] += coeffj / molarMass[l];
				}
				coeffj = molarMass[thisIndex] * nu[j] * w[i];
				for ( l = 0; l < speciesPerReaction; ++l ) {
					indexL = speciesIndex[l];
					if ( nu[l] > 0 && !species->IsSteadyState( indexL ) && Y[indexL] > 0.0 ) {
						dmdY[indexL][thisIndex] -= coeffj * nu[l] / Y[indexL];
					}
				}
				dmdT[thisIndex] -= molarMass[thisIndex] * nu[j] * dwdT;
			}
		}
	}

#ifdef SGINJACOBIAN
	if ( fSurfaceGrowth || fSurfaceOxidation ) {
		Double	density = flameNode->mixDensity[kCurr];
		Double	*kSoot = fSootRateCoeffs->vec;
		Double	*wSoot = fSootReactionRate->vec;
		Double	*a = fASoot->vec;
		Double	*n = fNSoot->vec;
		Double	*eOverR = fEOverRSoot->vec;
		
		ComputeSootRateCoeffs( kSoot, temp, reaction );
		ComputeSootReactionRates( wSoot, Y, flameNode->moments, density, molarMass, mixMolarMass );
		
		for ( i = 0; i < nSootReactions; ++i ) {
			speciesIndex = fSpecNumSoot[i]->vec;
			nu = fNuSoot[i]->vec;
			speciesPerReaction = fSpecNumSoot[i]->len;
			
			for ( l = 0, sumNu = 0.0; l < speciesPerReaction; ++l ) {
				if ( nu[l] > 0.0 ) {
					sumNu += nu[l];
				}
			}
			
			dwdT = wSoot[i] / temp * ( n[i] + eOverR[i] / temp - sumNu );
	
			coeffi = wSoot[i] * sumNu * mixMolarMass;
			for ( j = 0; j < speciesPerReaction; ++j ) {
				thisIndex = speciesIndex[j];
				if ( !species->IsSteadyState( thisIndex ) ) {
					coeffj = molarMass[thisIndex] * coeffi * nu[j];
					for ( l = 0; l < nSpeciesIn; ++l ) {
						dmdY[l][thisIndex] += coeffj / molarMass[l];
					}
					coeffj = molarMass[thisIndex] * nu[j] * wSoot[i];
					for ( l = 0; l < speciesPerReaction; ++l ) {
						indexL = speciesIndex[l];
						if ( nu[l] > 0 && !species->IsSteadyState( indexL ) && Y[indexL] > 0.0 ) {
							dmdY[indexL][thisIndex] -= coeffj * nu[l] / Y[indexL];
						}
					}
					dmdT[thisIndex] -= molarMass[thisIndex] * nu[j] * dwdT;
				}
			}
		}
	}
#endif
/*	if ( fSurfaceGrowth ) {
		Double	*wSoot = fSootReactionRate->vec;
		Double	mSoot;
		const Double	orderOfReaction = 2.0; // all soot reactions are bimolecular
		Double	density = flameNode->mixDensity[kCurr];
		Double	*a = fASoot->vec;
		Double	*n = fNSoot->vec;
		Double	*eOverR = fEOverRSoot->vec;
		
		ComputeSootReactionRates( wSoot, Y, flameNode->moments, density, molarMass, mixMolarMass );
		
		mSoot = molarMass[f_H] * ( -wSoot[ks8f] + wSoot[ks8b] - wSoot[ks9] 
										+ wSoot[ks10f] - wSoot[ks10b] );
		dmdT[f_H] += mSoot / temp * ( n[f_H] + eOverR[f_H] / temp - orderOfReaction );

		mSoot = molarMass[f_H2] * ( wSoot[ks8f] - wSoot[ks8b] );
		dmdT[f_H2] += mSoot / temp * ( n[f_H2] + eOverR[f_H2] / temp - orderOfReaction );

		mSoot = molarMass[f_C2H2] * (  -wSoot[ks10f] + wSoot[ks10b] );
		dmdT[f_C2H2] += mSoot / temp * ( n[f_C2H2] + eOverR[f_C2H2] / temp - orderOfReaction );

		mSoot = molarMass[f_O2] * -wSoot[ks11];
		dmdT[f_O2] += mSoot / temp * ( n[f_O2] + eOverR[f_O2] / temp - orderOfReaction );

		mSoot = molarMass[f_CO] * 2.0 * wSoot[ks11];
		dmdT[f_CO] += mSoot / temp * ( n[f_CO] + eOverR[f_CO] / temp - orderOfReaction );

		mSoot = molarMass[f_OH] * -wSoot[ks12];
		dmdT[f_OH] += mSoot / temp * ( n[f_OH] + eOverR[f_OH] / temp - orderOfReaction );

		mSoot = molarMass[f_CH] * wSoot[ks12];
		dmdT[f_CH] += mSoot / temp * ( n[f_CH] + eOverR[f_CH] / temp - orderOfReaction );

		mSoot = molarMass[f_CHO] * wSoot[ks12];
		dmdT[f_CHO] += mSoot / temp * ( n[f_CHO] + eOverR[f_CHO] / temp - orderOfReaction );
	}*/
#else
	return;
#endif
}

void T1DSoot::PrintPAHj( int j, TNewtonPtr bt )
{
	int			i, k;
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	Double		***Pij = fPij->tensor;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		dummy;
	int			counter;
	
	if ( j < 1 || j > 6 ) {
		cerr << "#error: PAH_" << j << " not existent" << NEWL;
		return;
	}
	
	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( bt->GetOutFileBuff(), "%sPAH%d_%d.dout", bt->GetOutputPath(), j, counter );
	if ( !( fp = fopen( bt->GetOutFileBuff(), "w") ) ) {
		cerr << "#warning: unable to open file" << bt->GetOutFileBuff() << NEWL;
		return;
	}

	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "y" );
	for ( i = 0; i < fNStages; ++i ) {
		fprintf( fp, "\tStage%-6d", i+1 );
	}
	fprintf( fp, "\n%-9E", left );
	for ( i = 0; i < fNStages; ++i ) {
		fprintf( fp, "\t%-9E", Pij[-1][i][j] );
	}
	for ( k = 0; k < nGridPoints; ++k ) {
		fprintf( fp, "\n%-9E", x[k] );
		for ( i = 0; i < fNStages; ++i ) {
			fprintf( fp, "\t%-9E", Pij[k][i][j] );
		}
	}
	fprintf( fp, "\n%-9E", right );
	for ( i = 0; i < fNStages; ++i ) {
		fprintf( fp, "\t%-9E", Pij[nGridPoints][i][j] );
	}

	fclose( fp );
}

void T1DSoot::PrintPAHOneStage( int i, TNewtonPtr bt )
{
	int			j, k;
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	Double		***Pij = fPij->tensor;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		dummy;
	int			counter;
	

	if ( i < 0 || i >= fNStages ) {
		cerr << "#error: PAH stage " << i << " not existent" << NEWL;
		return;
	}
	
	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( bt->GetOutFileBuff(), "%sPAHStage%d_%d.dout", bt->GetOutputPath(), i, counter );
	if ( !( fp = fopen( bt->GetOutFileBuff(), "w") ) ) {
		cerr << "#warning: unable to open file" << bt->GetOutFileBuff() << NEWL;
		return;
	}

	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "y" );
	for ( j = 1; j <= fNPAHMolecules; ++j ) {
		fprintf( fp, "\tPAH%-6d", j );
	}
	fprintf( fp, "\n%-9E", left );
	for ( j = 1; j <= fNPAHMolecules; ++j ) {
		fprintf( fp, "\t%-9E", Pij[-1][i][j] );
	}
	for ( k = 0; k < nGridPoints; ++k ) {
		fprintf( fp, "\n%-9E", x[k] );
		for ( j = 1; j <= fNPAHMolecules; ++j ) {
			fprintf( fp, "\t%-9E", Pij[k][i][j] );
		}
	}
	fprintf( fp, "\n%-9E", right );
	for ( j = 1; j <= fNPAHMolecules; ++j ) {
		fprintf( fp, "\t%-9E", Pij[nGridPoints][i][j] );
	}

	fclose( fp );
}

void T1DSoot::PrintPAHMoments( TNewtonPtr bt )
{
	int			j, k;
	FILE		*fp;
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	Double		**pahMoments = fPAHMoments->mat;
	Double		**sumpi = fSumPi->mat;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		dummy;
	int			counter;
	

	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( bt->GetOutFileBuff(), "%sPAHMoments_%d.dout", bt->GetOutputPath(), counter );
	if ( !( fp = fopen( bt->GetOutFileBuff(), "w") ) ) {
		cerr << "#warning: unable to open file" << bt->GetOutFileBuff() << NEWL;
		return;
	}

	fprintf( fp, "*\n" );
	fprintf( fp, "%-12s", "y" );
	for ( j = 0; j < fNPAHMoments; ++j ) {
		fprintf( fp, "\tM_%-6d", j );
	}
	fprintf( fp, "\tsumPi" );
	fprintf( fp, "\n%-9E", left );
	for ( j = 0; j < fNPAHMoments; ++j ) {
		fprintf( fp, "\t%-9E", pahMoments[-1][j] );
	}
	fprintf( fp, "\t%-9E", sumpi[-1][0] );
	for ( k = 0; k < nGridPoints; ++k ) {
		fprintf( fp, "\n%-9E", x[k] );
		for ( j = 0; j < fNPAHMoments; ++j ) {
			fprintf( fp, "\t%-9E", pahMoments[k][j] );
		}
		fprintf( fp, "\t%-9E", sumpi[k][0] );
	}
	fprintf( fp, "\n%-9E", right );
	for ( j = 0; j < fNPAHMoments; ++j ) {
		fprintf( fp, "\t%-9E", pahMoments[nGridPoints][j] );
	}
	fprintf( fp, "\t%-9E", sumpi[nGridPoints][0] );

	fclose( fp );
}

void T1DSoot::PrintSootInfo( TNewtonPtr bt )
{
	TGridPtr	grid = bt->GetGrid()->GetCurrentGrid();
	Double		*x = grid->GetX()->vec;
	int			nGridPoints = bt->GetCurrentGridPoints();
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	FILE		*fp;
	Double		dummy;
	int			counter;
	
	counter = ( int ) ( modf( bt->GetNIter()/10.0, &dummy ) * 10.0 );
	sprintf( bt->GetOutFileBuff(), "%sSootInfo%d.tout", bt->GetOutputPath(), counter );
	if ( !( fp = fopen( bt->GetOutFileBuff(), "w") ) ) {
		cerr << "#warning: unable to open file" << bt->GetOutFileBuff() << NEWL;
		return;
	}
//	PrintSootReactions( flame );

	fclose( fp );
}
#endif // ZEROD
