//#define FOOL_SOFTBENCH( x ) 

#ifndef __FlameMan__

#define __FlameMan__

#define VERSION "3.3.9"

#ifdef powerc
#define OLDROUTINENAMES 0
#include <fp.h>
#endif

#ifndef __StandardHeaders__
#define __StandardHeaders__
#include "ArrayManager.h"
#include <iostream>
#include <fstream>
#include <cstdio> 
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include <csignal>
#endif

#include "ListTool.h"
#include "Newton.h"
#include "SmallNewton.h"
#include "ScanManStructs.h"

#undef NOXPRODRATES

#define LENOFBUFFER 128

#ifndef __MYIOCONST__
#define __MYIOCONST__
#define TAB "\t"
#define NEWL "\n"
#endif // __MYIOCONST__

enum FlameType				{ kNotSpecified, kStartUpDiffusion, kDiffusionPhys, kDiffPhysEigen, kStartUpDiffusion2, 
								kCountPremPhys, kCountPremSim, kCountDiffMix, kCountDiffCont,
								kHomoIsoChor, kHomoIsoBar, kTransFlamelet, kTrans1DIsoChor, kUnstrPremPhys, kTransCoalParticle };
enum MixtureSpecification	{ kMassFlux = 1, kMassFraction, kMolarFraction };
enum BoundaryCondition		{ kNone = 0, kDirichlet, kNeumann };
enum GridPoint				{ kPrev = -1, kCurr, kNext };
enum Sign					{ kNegative, kPositive };
enum ContinSide				{ kNoSide, kLeft, kRight, kBothSides };
enum CoordType				{ kPhysical, kSimilarity, kMixtureFraction };
enum ContinType				{ kNoCont, kTemperature, kVelocity, kMomentum, kSeshaTemp, kPolytrope, kPressure };
enum EqOfState				{ kDensFromPress, kPressFromDens };
enum CountPremType			{ kBackToBack = 1, kFreshToBurned, kFull};

#define SMALL 1.0e-10

#ifndef RGAS
#define RGAS 8314.34 // [J / kmole K]
#endif
#ifndef AVOGADRO
#define AVOGADRO 6.0221367e26  // [1 / kmole]
#endif
#ifndef STEFBOLTZ
#define STEFBOLTZ 5.67051e-8    // [W / (m^2 K^4)]
#endif

#define MYPOW( a, b ) ( ( (a) <= 0.0 ) ? 0.0 : exp( (b) * log(a) ))

inline
Double CatchZero( Double a )
{
	return ( a == 0.0 ) ? 1.0e-20 : a;
}

typedef const char *const ConstString; 
//typedef const char *const *ConstStringArray; 

#ifndef __CLASSES__
class TReaction;
class TSpecies;
class TProperties;
class TFlame;
class T0DFlame;
class T1DFlame;
struct BoundaryInput;
struct FirstInput;

typedef TReaction *TReactionPtr;
typedef TSpecies *TSpeciesPtr;
typedef TProperties *TPropertiesPtr;
typedef TFlame *TFlamePtr;
typedef T0DFlame *T0DFlamePtr;
typedef T1DFlame *T1DFlamePtr;
typedef BoundaryInput *BoundaryInputPtr;
typedef FirstInput *FirstInputPtr;
#endif // __CLASSES__

//*************************************************************************************************************
class TFlameNode {
public:
	TFlameNode( void ) { InitTFlameNode(); };

// reaction
	Double		*tBodyConc;
	Double		*rateCoeff;
	Double		*tempReaction;
	Double		*pressureReaction;
	Double		*currRateCoeff;
	Flag		*kNewerThanW;
	Double		*reactionRate;
	Double		*YReaction;
	Double		*currReacRate;
	Double		**dMdY;
	
// species
	Double			*viscosity;
	Double			*heatCapacity;
	Double			*conductivity;
	Double			*enthalpy;
	Double			*diffusivityPrev;
	Double			*diffusivity;
	Double			*diffusivityNext;
	Double			**diffTherm;
	Double			*productionRate;
	Double			*savedY;
	Double			*savedDeltaiY;
	Double			*sumDiff;
	Double			*deltaI;
	Double			*tempProp;
	Double			*pressureProp;
	Double			***binDij;
	Double			**GijOverWj;
	Double			**OmegaDOverDCoeff;
	Double			**DThermConst;

// properties
	Double	*mixViscosity;
	Double	*mixDensity;
	Double	*mixConductivity;
	Double	*mixHeatCapacity;
	Double	*mixMolarMass;
	
// flame
	Double				**Y;
	Double				*temp;
	Double				*radiation;
	Double				*diffCorr;
	Double				viscosityInf;	// viscosity at right boundary
	Double				rhoInf;			// density at right boundary

// soot
	Double				**Pij;
	Double				*sumPi;
	Double				*pahMoments;
	Double				*moments;
	Double				*pahReactionRate;
//	Double				*sootReactionRate;
	Double				*diffSoot;
	Double				*sootSource;
	Double				**dMdx;
	
private:
	void InitTFlameNode( void ) {;};
};
typedef TFlameNode *TFlameNodePtr;

//*************************************************************************************************************
struct FirstInput {
	FirstInput( int argc, char *argv[] ) 
		: fVVelocityOffset(0), fUVelocityOffset(1), fTemperatureOffset(2), fPressureOffset(3)
		{ InitFirstInput( argc, argv ); };
	FirstInput( char *fInput ) 
		: fVVelocityOffset(0), fUVelocityOffset(1), fTemperatureOffset(2), fPressureOffset(3)
		{ InitFirstInput( 0, &fInput ); };
	~FirstInput( void );

	void			PrintAdditionalData( void );

	const int		fVVelocityOffset;
	const int		fUVelocityOffset;
	const int		fTemperatureOffset;
	const int		fPressureOffset;

	// Exact Backward Reaction Constants
	Flag                    fExactBackward;

	// 0D-stuff
	Flag			fAdditionalOutput;
	Flag			fEquidistant;
	int				fNOutputs;
	Double			fArtificialSource;
	Double			fArtSourceTime;
	
	// rest
	Flag			fScannerProgress;

	char			*fDateCreated;
	char			*fAuthor;
	int				fFlameType; // uses enum FlameType
	String			fContinType;
	String			fContinSide;
	Double			fContBound;
	Flag			fConstantLewisNumber;
	Flag			fWithRadiation;
	Flag			fArcLengthContin;
	Flag			fStrainRateContin;
	Flag			fThermoDiffusion;
	Flag			fWithSoot;
	int				fNSootMoments;
	Flag			fCompUnPhysChain;
	int				fNUnPhysChain;
	Flag			fIsAxiSymmetric;
	Flag			fClipNegativeConcs;
	Flag			fNoDiffCorr;
	Double			fMinStrainRate;
	Double			fPressure[1000];
	int				fNPressures;
	Double			fPressureComm;	// pressure specified on command line
	Double			fParameterComm;	// parameter specified on command line
	int				fPremConfiguration;

	// Sensitivity Analysis
        Flag		        fSensAnal;
        Flag                    fSensObjAll;     // all species are SensObj
        Flag                    fSensAnalSpec;   // Sens Anal on Species
        Flag                    fSensAnalReac;   // Sens Anal on Reactions
        Double                  fSensAnalFac;    // Factor by which A is multiplied
        Flag                    fSensMax;        // sensitivity on max 
        Flag                    fSensFinal;      // sensitivity on final values
	char			*fSensObj[1000]; // Name of species/variables used as target
	int			fNSensObj;
	Flag			fReactionFluxes;

	Flag			fPrintRHSSpecies;
	Flag			fPrintRHSTemp;
	Double			fStrainRate[1000];
	int			nStrainRates;
	Double			fDissRate[1000];
	int			nDissRates;
	int			fNPhi;
	Double			fPhi[1000];
	Flag			fKeepMassFracs;
	Flag			fLiquidPoolBC;
	Double			fDeltaSCont;
//	Double			fChi;
	char			*fGlobalReaction;
	char			*fFuelName[1000];
	int			fNFuels;
	char			*fFromSpeciesName;
	char			*fToSpeciesName;
	Double			fContInc;
	char			*fOxName;
	Flag			fPrintMolarFractions;
	Flag			fSteadyStatesNumerical;
	int			fVariablesWithoutSpecies;

	Flag			fNucleation;
	Flag			fCondensation;
	Flag			fCoagulation;
	Flag			fSurfaceGrowth;
	Flag			fSurfaceOxidation;
	Flag			fThermoPhoresis;
	Flag			fOHPAHOxidation;
	Flag			fO2PAHOxidation;
	Double			fCoagFact;
	Flag			fSootRadiation;
	Flag			fSootUpdateProdRates;
	Flag			fSizeDepDiff;
	Flag			fSurfDepCoag;

	BoundaryInputPtr	fInitialCond;
	BoundaryInputPtr	leftBoundary;
	BoundaryInputPtr	rightBoundary;
	
	Flag			fWriteBT;
	Flag			fWriteResiduum;
	Flag			fWatchGridding;
	Flag			fWriteEverySolution;
	char			*fOutputPath;
	char			*fOutputPathComm;

	Flag			fUseModifiedNewton;
	Flag			fUseNumericalJac;
        Flag                    fUseSecOrdJac;
	Flag			fUseNumericalDM;
	Flag			fWriteFullRes;
	int			fInitialEquations;
	int			fMaxGridPoints;
	int			fInitialGridPoints;
	Flag			fDampFlag;
	Flag			fTimeDepFlag;
	Flag			fContinFlag;
	int			fDeltaNewGrid;
	Flag			fOneSolOneGrid;
	Flag			fAdjustComputationalDomain;
	Double			fTolRes;
    Double			fTolDy;
    int				fMaxIter;
    Double			fLeft;
    Double			fRight;
    int				fGamma;
	Double			fKappa;
	Double			fTauGrid;
    Double			fR;
    Double			fQ;
	Double			fStart;
	Double			fEnd;
	Double			fAlpha;
    Double			fLambdaMin;
    Double			fDeltaTStart;
    Double			fDeltaTMax;
    int				fContSteps;

	Double			fCl;
	Double			fHv;
	Double			fT_B;

	char			*fLewisNumberFile;
	char			*fReactionFile;
// THIS FILE NAME WAS GIVEN A LENGTH OF 31 BEFORE.  
// THIS LED TO THE FOLLOWING POINTER BEING OVERWRITTEN.
// SHOULD CHECK THAT IT DOESN'T OVERWRITE ANYTHING OR DYNAMICALLY ALLOCATE
	char			fAdditionalInfoFile[LENOFBUFFER];
//	char			*fAdditionalInfoFile;
	char			*fMechanismFileComm;	// name of mechanism file from the command line
	char			*fStartProfileFileComm;	// name of startprofiles file from the command line
	char			*fStartProfileFile;
	char			*fOutFileName;
	char			*fAddFileNo1;
	FILE			*fpR;			// filepointer to fReactionFile
	FILE			*fpA;			// filepointer to fAdditionalInfoFile
	FILE			*fpS;			// filepointer to fStartProfileFile
	FILE			*fpO;			// filepointer to fOutputFile

private:
	void	InitFirstInput( int argc, char *argv[] );
	void	ParseCommandLine( int argc, char *argv[] );
	int		yylex( void );
	
	char	stringBuffer[LENOFBUFFER];
};
//*************************************************************************************************************
struct BoundaryInput {
	BoundaryInput( FirstInputPtr firstInp, int len ) 
		: fVVelocityOffset(firstInp->fVVelocityOffset), fUVelocityOffset(firstInp->fUVelocityOffset)
		, fTemperatureOffset(firstInp->fTemperatureOffset)
		, fPressureOffset(firstInp->fPressureOffset), fNVars(4)
		{ InitBoundaryInput( len ); };
	~BoundaryInput( void );

	BoundaryInput& 	operator=( const BoundaryInput& boundary );
	void			PrintBoundary( FILE *fp );
	
	const int	fVVelocityOffset;
	const int	fUVelocityOffset;
	const int	fTemperatureOffset;
	const int	fPressureOffset;
	const int	fNVars;

	// the following variables contain bc's for the species
	int		fSpecifiedSpeciesBCs;
	int		fMixtureSpecification;	//	uses enum MixtureSpecification
	int		fBcFlagSpecies;			//	uses enum BoundaryCondition
	char	**speciesName;
	Double	*fValueSpecies;			//	contains value of bc for species

	// the following variables contain bc's for f ( or V ), f' ( or U ), T and p
	int		*fBcFlag;				//	specifies kind of bc; uses enum BoundaryCondition
	Double	*fValue;				//	contains value of bc
	
private:
	void	InitBoundaryInput( int len );
	int		fLen;
};
//*************************************************************************************************************

class TInputData {
public:
	TInputData( FirstInputPtr firstInp ) 
		: fVVelocityOffset(firstInp->fVVelocityOffset),fUVelocityOffset(firstInp->fUVelocityOffset),
		fPressureOffset(firstInp->fPressureOffset),
		fTemperatureOffset(firstInp->fTemperatureOffset) { InitInputData( firstInp ); };
	~TInputData( void );
	
	HeaderPtr		GetHeader( void ) { return fHeader; };
	CounterPtr		GetCounter( void ) { return fCounter; };
	ReactionPtr		GetReactions( void ) { return fReaction; };
	ReactionPtr		GetPAHReactions( void ) { return fPAHReaction; };
	ReactionPtr		GetSootReactions( void ) { return fSootReaction; };
	ThirdBodyPtr	GetThirdBody( void ) { return fThirdBody; };
	SpeciesPtr		GetSpecies( void ) { return fSpecies; };
	DimensionPtr	GetDimension( void ) { return fDimension; };
	int				FindSpecies( ConstString species );
	int				FindAtomIndex( ConstString atom );
	int				FindSpeciesCoefficient( int speciesIndex, int reactionIndex );
	void			PrintAdditionalData( void );

	const int		fVVelocityOffset;
	const int		fUVelocityOffset;
	const int		fTemperatureOffset;
	const int		fPressureOffset;

        // Exact Backward Reaction Constants
        Flag                    fExactBackward;

	// 0D-stuff
	Flag			fAdditionalOutput;
	Flag			fEquidistant;
	int				fNOutputs;
	Double			fArtificialSource;
	Double			fArtSourceTime;
	
	char			*fDateCreated;
	char			*fAuthor;
	int				fFlameType; // uses enum FlameType
	ContinType		fContinType;
	ContinSide		fContinSide;
	Double			fContBound;
	Flag			fConstantLewisNumber;
	Flag			fWithRadiation;
	Flag			fArcLengthContin;
	Flag			fStrainRateContin;
	Flag			fThermoDiffusion;
	Flag			fWithSoot;
	int				fNSootMoments;
	Flag			fCompUnPhysChain; // also flag for ConstMassFlux in TUnstrechedPremixed
	int				fNUnPhysChain;
	Flag			fIsAxiSymmetric;
	Flag			fClipNegativeConcs;
	Flag			fNoDiffCorr;
	Double			fMinStrainRate;
	VectorPtr		fPressure;
	Double			fParameterComm;	// parameter specified on command line
	int				fPremConfiguration;

	// Sensitivity Analysis
	Flag		        fSensAnal;
	Flag                    fSensObjAll;         // all species are SensObj
	Flag                    fSensAnalSpec;       // Sens Anal on Species
	Flag                    fSensAnalReac;       // Sens Anal on Reactions
	Double                  fSensAnalFac;        // Factor by which A is multiplied
	Flag                    fSensMax;            // sensitivity on max and location of max
	Flag                    fSensFinal;          // sensitivity on final values
	char			**fSensObj;          // Name of species/variables used as target
	int			fNSensObj;
	Flag			fReactionFluxes;

	Flag			fPrintRHSSpecies;
	Flag			fPrintRHSTemp;
	VectorPtr		fStrainRate;
	VectorPtr		fDissRate;
	VectorPtr		fPhi;
	Flag			fKeepMassFracs;
	Flag			fLiquidPoolBC;
	Double			fDeltaSCont;

	char			*fGlobalReaction;
	IntVectorPtr	fFuelIndex;
	int				fFromSpeciesIndex;
	int				fToSpeciesIndex;
	Double			fContInc;
	int				fOxIndex;
	int				fH2OIndex;
	int				fCO2Index;
	Flag			fPrintMolarFractions;
	Flag			fSteadyStatesNumerical;
	int				fVariablesWithoutSpecies;
	char			*fPAHSymbol;
	char			*fSootSymbol;

	Flag			fNucleation;
	Flag			fCondensation;
	Flag			fCoagulation;
	Flag			fSurfaceGrowth;
	Flag			fSurfaceOxidation;
	Flag			fThermoPhoresis;
	Flag			fOHPAHOxidation;
	Flag			fO2PAHOxidation;
	Double			fCoagFact;
	Flag			fSootRadiation;
	Flag			fSootUpdateProdRates;
	Flag			fSizeDepDiff;
	Flag			fSurfDepCoag;

	BoundaryInputPtr	fInitialCond;
	BoundaryInputPtr	leftBoundary;
	BoundaryInputPtr	rightBoundary;

	Flag			fWriteBT;
	Flag			fWriteResiduum;
	Flag			fWatchGridding;
	Flag			fWriteEverySolution;
	char			*fOutputPath;

	Flag			fWriteFullRes;
	Flag			fUseModifiedNewton;
	Flag			fUseNumericalJac;
        Flag                    fUseSecOrdJac;
	Flag			fUseNumericalDM;
	int			fNVariables;
	int			fInitialEquations;
	int			fMaxGridPoints;
	int			fInitialGridPoints;
	Flag			fDampFlag;
	Flag			fTimeDepFlag;
	Flag			fContinFlag;
	int			fDeltaNewGrid;
	Flag			fOneSolOneGrid;
	Flag			fAdjustComputationalDomain;
	Double			fTolRes;
    Double			fTolDy;
    int				fMaxIter;
    Double			fLeft;
    Double			fRight;
    int				fGamma;
	Double			fKappa;
	Double			fTauGrid;
    Double			fR;
    Double			fQ;
	Double			fStart;
	Double			fEnd;
	Double			fAlpha;
    Double			fLambdaMin;
    Double			fDeltaTStart;
    Double			fDeltaTMax;
    int				fContSteps;

	Double			fCl;
	Double			fHv;
	Double			fT_B;

	char			*fLewisNumberFile;
	char			*fStartProfileFile;
	char			*fAddFileNo1;
	FILE			*fpS;						// filepointer to fStartProfileFile
	char			*fReactionFile;
	char			fAdditionalInfoFile[LENOFBUFFER];
	char			*fOutFileName;
	

private:
	void 			InitInputData( FirstInputPtr firstInp );	// calls ReadInReactions( void ) 
   	void			ReadInAdditionalData( void );
	void			ReadInReactions( void );
	
	HeaderPtr 		NewHeader( void );
	HeaderPtr 		ReadHeader( void );
	CounterPtr		ReadCounter( void );
	AtomsPtr 		NewAtomArray( int len );
	AtomsPtr 		ReadAtoms( void );
	SpeciesPtr 		NewSpeciesArray( int len );
	SpeciesPtr 		ReadSpecies( void );
	ReactionPtr 	NewReactionArray( int len );
	ReactionPtr 	ReadReaction( int len );
	ThirdBodyPtr 	NewThirdBodyArray( int len );
	ThirdBodyPtr 	ReadThirdBody( void );
	DimensionPtr 	NewDimensionArray( int len );
	DimensionPtr 	ReadDimension( void );
	void			FreeAtomsArray( AtomsPtr atoms, int len );
	void			FreeSpeciesArray( SpeciesPtr species, int len );
	void			FreeReactionArray( ReactionPtr reaction, int len );
	void			FreeThirdBodyArray( ThirdBodyPtr thirdBody, int len );
	void			FreeDimensionArray( DimensionPtr dimension, int len );
	int				GetSpeciesNumber( void );

	CounterPtr		fCounter; 	// initialized in InitInputData( void )
	AtomsPtr		fAtoms;
	SpeciesPtr		fSpecies;
	ReactionPtr		fReaction;
	ReactionPtr		fPAHReaction;
	ReactionPtr		fSootReaction;
	ThirdBodyPtr	fThirdBody;
	DimensionPtr	fDimension;
	HeaderPtr		fHeader;
	String			fBuffer; 	// initialized in ReadHeader( void )
		
	FILE	*fpR;						// filepointer to fReactionFile
	FILE	*fpA;						// filepointer to fAdditionalInfoFile
	FILE	*fpO;						// filepointer to fReactionFile
};
typedef TInputData *TInputDataPtr;

class SteadyStateInfo {
public:
	Double	*c;
	Double	*k;
	Double	*M;
	int		speciesIn;
	
	void SetSteadyStateInfo( Double *cIn, Double *kIn, Double *MIn, int speciesIn );
	void GetSteadyStateInfo( Double **cIn, Double **kIn, Double **MIn, int *speciesIn );
};
typedef SteadyStateInfo *SteadyStateInfoPtr;

//*************************************************************************************************************
class TReaction {
/* number of species of each reaction is stored in VectorPtr nu ( nu->len )
 * 'if ( aInf )' is a Flag for the pressure dependence 
 * TReaction needs the following data: T, pressure, Y_{i} 
 */

public:
	TReaction( TInputDataPtr input ) { InitTReaction( input ); };
	~TReaction( void );
	
	int				GetNOfReactions( void ) { return fA->len; };
	int				GetNTBOfSpecies( void ) { return ( fTBSpeciesCoeffs ) ? fTBSpeciesCoeffs->rows : 0; };
	int				GetNOfSpeciesPerReaction( int number ) { return fNu[number]->len; };
	int				GetNOfThirdBodies( void ) { return ( fTBSpeciesCoeffs ) ? fTBSpeciesCoeffs->cols : 0; };
	Flag			IsWithLindemann( int index ) { return ( ( fAInf->vec[index] ) ? TRUE : FALSE ); };
	Flag			IsWithThirdBody( int index ) { return fWithThirdBody[index]; };
	Flag			IsBackwardReaction( int index ) { return ( fForwardReac->vec[index] >= 0 ) ? TRUE : FALSE; };
	Flag			IsForwardReaction( int index ) { return ( fBackwardReac->vec[index] >= 0 ) ? TRUE : FALSE; };
	Flag			IsReducedMech( void ) { return fReducedMech; };
	IntVectorPtr	*GetSpeciesNumber( void ) { return fSpeciesNumber; };
	VectorPtr		*GetNu( void ) { return fNu; };
	VectorPtr		GetA( void ) { return fA; };
	VectorPtr		GetN( void ) { return fN; };
	VectorPtr		GetEOverRgas( void ) { return fEOverRgas; };
	MatrixPtr		GetTBSpeciesCoeffs( void ) { return fTBSpeciesCoeffs; };
	int				GetTBNumber( int index ) { return fThirdBodyNumber->vec[index]; };
	char			**GetLabels( void ) { return fLabels; };
	IntVectorPtr	GetBackwardReacs( void ) { return fBackwardReac; };
	IntVectorPtr	GetForwardReacs( void ) { return fForwardReac; };
	int				PrintReactionEquation( int number, TSpeciesPtr species, FILE *fp );
	void			PrintReactionEquation( int number, TSpeciesPtr species, char *string );
	void                    UpdateReactionRates(Double *w, Double *Ynew, Double *Yold, Double *molarMass, Double mixMolarMass, int index);
	void			ComputeReactionRates( Double *w, Flag &kNewerThanW, Double *currReacRate, Double *k
						, Double *tBConc, Double density, Double *Y
						, Double *currConc, Double *molarMass, TSpeciesPtr species );
	void			ComputeReactionRate( int i, Double &w, Double &k, Double *tBConc, Double density, Double *Y, Double *molarMass );
	void			ComputeReactionRate( int i, Double &w, Double &k, Double *c, Double *tBConc );
	void			ComputeRateCoefficients( Double *k, Double *fCurrRateCoefficients, Flag &kNewerThanW,
				Double temp, Double &currTemp, Double pressure, Double &currPressure, Double *tBConc, TSpeciesPtr fSpecies );
        void                    ComputeBackwardRateCoefficients( Double *k, Double temp, TSpeciesPtr fSpecies);
	void			ComputeRateCoefficient( int number, Double &k, Double temp, Double pressure, Double *tBConc );
	void			ComputeRateCoefficient( Double temp, Double &k, Double a, Double n, Double eOverR );
	void			CompThirdBodyConcs( Double *tBodyConc, Double *Y, Double *molarMass, Double density );
	void			CompThirdBodyConc( int i, Double &tBodyConc, Double *Y, Double *molarMass, Double density );
	void			FilldMdTOnePointAnal( Double *dMdT, Double temp, Double *reactionRate
						, Double *rateCoeff, Double pressure, Double *molarMass, Double *tBConc );
	void			FilldMdYOnePointAnal( Double **dMdY, Double *Y, Double *reactionRate
						, Double mixMolarMass, Double *molarMass, Double *tBodyConc
						, Double *productionRate, Double mixDensity );
	void			CheckSteadyStatesMech( const char *mechName );
	void			CheckBroadeningFactors( const char *mechName );
	void			WriteRoggsMechanismData( TSpeciesPtr species, TInputDataPtr input );
	int				PrintRoggsReactionEquation( int number, TSpeciesPtr species, FILE *fp );
	void			ComputeConcs( Double *c, Double *Y, Double *molarMass, Double rho );
	VectorPtr		GetMolarConcs( void ) { return fMolarConcs; };

	//PP detailed Heat release computation used for Reduction purposes
	Double                  CompDetailedHeatRelease( int number, Double reacRate, Double *enthalpy, Double *W);
	//PP

protected:
	void 			InitTReaction( TInputDataPtr input );
	void			FillReactions( TInputDataPtr input );
	Double			ComputedkdT( int reactionIndex, Double temp, Double pressure, Double rateCoeff, Double *tBConc );
	int				GetNSpeciesForTB( void ) { return ( fTBSpeciesCoeffs ) ? fTBSpeciesCoeffs->rows : 0; };
	int				GetBackwardReaction( int i ) { return GetReverseReaction( i, 'f', 'b' ); };
	int				GetForwardReaction( int i ) { return GetReverseReaction( i, 'b', 'f' ); };
	int				GetReverseReaction( int i, char is, char lookfor );
	int				FindReaction( char *label );
	void			ComputeSteadyStates( Double *c, Double *Y, Double *k, Double *molarMass
							, Double *tBConc, Double density );
//	void			ComputeSteadyStates( Double *k, Double *c, Double *tBConc );
	void			ComputeGlobalReactionRates( Double *w, Double *k, int nOfReactions );

	
	VectorPtr		fA;	// [m, s, kmole, K]					// input; VectorPtr of length nOfReactions
	VectorPtr		fN;										// input; VectorPtr of length nOfReactions
	VectorPtr		fEOverRgas;	// [J / kmole]				// input; VectorPtr of length nOfReactions
	VectorPtr		fAInf;	// [m, s, kmole, K]				// input; VectorPtr of length nOfReactions
	VectorPtr		fNInf;									// input; VectorPtr of length nOfReactions
	VectorPtr		fEInfOverRgas;	// [J / kmole]			// input; VectorPtr of length nOfReactions
	enum			FcCoeff	{ kFca = 0, kFcb, kFcc, kFcTa, kFcTb, kFcTc, kFcLen };
	MatrixPtr		fFcCoeff;								// MatrixPtr of length fFcCoeff[nOfReactions][kFcLen]
	IntVectorPtr	fLindemannNumber;						// input; index for the broadeningfunctionarray
	IntVectorPtr	*fSpeciesNumber;						// input; see fNu; 
	VectorPtr		*fNu;									// nu_product is less than zero; last element is fNu[nOfReactions][speciesPerReaction]
	VectorPtr		fSumEducts;									// last element is fSumEducts[nOfReactions]
															// contains sum_over_educts(nu')+(withThirdBody)?1:0
	IntVectorPtr	fNEducts;	// last element is fNEducts[nOfReactions]
					// contains number of educts for each reaction 
					// input; pointer to nOfReactions VectorPtr; negativ for products 
	// Exact Backward Reaction Constants
        Flag                    fExactBackward;														// 	of length nOfSpeciesPerReaction 
	Flag			*fPartialEquilibrium;					// for later use
	Flag			*fWithThirdBody;
	String			*fThirdBodyName;						//  
	IntVectorPtr	fThirdBodyNumber;						// fThirdBodyNumber can directly be used as index for the following
	IntVectorPtr	*fTBSpeciesNumber;						// input; pointer to nOfThirdBodies IntVectorPtr of length nOfSpecies
	MatrixPtr		fTBSpeciesCoeffs;						// last element is fTBSpeciesCoeffs[nOfThirdBodies][nOfSpecies]
	char			**fLabels;
	IntVectorPtr	fBackwardReac;							// offset of backward reaction IntVectorPtr of length nOfReactions
	IntVectorPtr	fForwardReac;							// offset of forward reaction IntVectorPtr of length nOfReactions
	Flag			fReducedMech;
	Flag			fGlobalMechanism;
	VectorPtr		fMolarConcs;							// length is nOfSpecies
	SteadyStateInfoPtr	fSsInfo;
	VectorPtr		fConcBoxVec;
	NewtonInfoPtr	fNewtonInfo;
	Flag			fSteadyStatesNumerical;

	int				fNSpeciesInSystem;
	int				fNSpecies;
};


//*************************************************************************************************************
class T0DReaction : public TReaction {
public:
	T0DReaction( TInputDataPtr input ) : TReaction( input ) { InitT0DReaction( input ); };
	~T0DReaction( void );
	
	VectorPtr		GetTBConc( void ) { return fTBConc; };
	VectorPtr		GetRateCoefficients( void ) { return fRateCoefficients; };
	VectorPtr		GetReactionRate( void ) { return fReactionRate; };
	Double			*GetCurrRateCoefficients( void ) { return fCurrRateCoefficients->vec; };
	Double			*GetCurrReactionRate( void ) { return fCurrReactionRate->vec; };
	Double			&GetCurrTemp( void ) { return fCurrTemp; };
	Double			&GetCurrPressure( void ) { return fCurrPressure; };
	Double			*GetCurrConc( void ) { return fCurrConc->vec; };
	Flag			&GetKNewerThanW( void ) { return fKNewerThanW; };

private:
	void 		InitT0DReaction( TInputDataPtr input );

	VectorPtr	fTBConc;	// last element is fTBConc[nOfThirdBodies]
								// contains the concentrations C_{M} = rho * Y_{M} / M_{M} = sum( rho * z_{i} * Y_{i} / M_{i} )
								// and can directly be multiplied to w_{k}
	VectorPtr	fRateCoefficients;	// [m, s, kmole]	// k_i; last element is fRateCoefficients[nOfReactions]
	VectorPtr	fCurrRateCoefficients;	// [m, s, kmole]	// k_i; last element is fRateCoefficients[nOfReactions]
	VectorPtr	fReactionRate;	// [kmole / (m^3*s)]	// w_i; last element is fReactionRate[nOfReactions]
	VectorPtr	fCurrReactionRate;	// [kmole / (m^3*s)]	// w_i; last element is fReactionRate[nOfReactions]
	Double		fCurrTemp;
	Double		fCurrPressure;
	VectorPtr	fCurrConc;
	Flag		fKNewerThanW;
};
typedef T0DReaction *T0DReactionPtr;

#ifndef ZEROD
//*************************************************************************************************************
class T1DReaction : public TReaction {
public:
	T1DReaction( TInputDataPtr input ) : TReaction( input )
			{ InitT1DReaction( input ); };
	~T1DReaction( void );
	
	int				GetNSpeciesInSystem( void ) { return fDmdy->phys_rows; };
	MatrixPtr		GetReactionRate( void ) { return fReactionRate; };
	MatrixPtr		GetRateCoefficients( void ) { return fRateCoefficients; };
	MatrixPtr		GetTBConcentrations( void ) { return fTBConcentrations; };
	TensorPtr		GetDmdy( void ) { return fDmdy; };
	void			PrintReactions( int k, TSpeciesPtr species );			// print all
	void			PrintReactions( int number, int gridPoint, TSpeciesPtr species, FILE *fp ); // print one
	void			PrintDmdY( T1DFlamePtr flame );
	void			PrintDmdY( T1DFlamePtr flame, int k, FILE *fp );
	void			PrintReactionRates( T1DFlamePtr flame );
	void			PrintRateCoeffs( T1DFlamePtr flame );
	void			PrintDetailedHeatRelease( T1DFlamePtr flame );
	void			UpdateDMdY( T1DFlamePtr flame, TAdaptiveGridPtr grid, Double pressure );
	void			FilldMdYOnePointAnal( TFlameNodePtr flameNode, Double *molarMass );
	void			FilldMdTOnePointAnal( TFlameNodePtr flameNode, Double pressure, Double *molarMass );
	VectorPtr		GetTempReaction( void ) { return fTempReaction; };
	VectorPtr		GetPressureReaction( void ) { return fPressureReaction; };
	MatrixPtr		GetCurrRateCoeff( void ) { return fCurrRateCoefficients; };
	Flag			*GetKNewerThanW( void ) { return fKNewerThanW; };
	MatrixPtr		GetYReaction( void ) { return fYReaction; };
	MatrixPtr		GetCurrReacRate( void ) { return fCurrReactionRate; };
	void			PrintProdRateGlobalReac( T1DFlamePtr flame );
	Double 			ComputeEIThermalNO( T1DFlamePtr flame );
	Double 			ComputeEIPromptNO( T1DFlamePtr flame );
	Double 			ComputeEINitrousNO( T1DFlamePtr flame );
	Double 			ComputeEIReburnNO( T1DFlamePtr flame );
	
private:
	void 			InitT1DReaction( TInputDataPtr input );
	void			ClearDmdy( void ) { ClearTensor( fDmdy ); };
	
	MatrixPtr		fRateCoefficients;	// [m, s, kmole]	// k_i; last element is fRateCoefficients[nOfGridPoints][nOfReactions]
	MatrixPtr		fCurrRateCoefficients;	// [m, s, kmole]	// k_i; last element is fRateCoefficients[nOfGridPoints][nOfReactions]
	MatrixPtr		fReactionRate;	// [kmole / (m^3*s)]		// w_i; last element is fReactionRate[nOfGridPoints][nOfReactions]
	MatrixPtr		fCurrReactionRate;	// [kmole / (m^3*s)]		// w_i; last element is fReactionRate[nOfGridPoints][nOfReactions]
	TensorPtr		fDmdy;									// last element is fDwdy->tensor[nOfGridPoints][nOfSpecies][nOfSpecies+1]
															// 	dm_i/dy_j = fDmdy->tensor[nOfGridPoints][j][i]
															// contains the derivatives with respect to massfractions and temperature
	MatrixPtr		fTBConcentrations;						// last element is fTBConcentrations[nOfGridPoints][nOfThirdBodies] contains the concentrations 
															// contains the concentrations C_{M} = rho * Y_{M} / M_{M} = sum( rho * z_{i} * Y_{i} / M_{i} )
															// and can directly be multiplied to w_{k}
	VectorPtr		fTempReaction;
	VectorPtr		fPressureReaction;
	MatrixPtr		fYReaction;
	Flag			*fKNewerThanW;
};
typedef T1DReaction *T1DReactionPtr;
#endif // ZEROD

#define NEWSURFGROWTH
//*************************************************************************************************************
enum RateCoeffs				{ k0f, k0b, k1f, k1b, k2f, k2b, k3f, k3b, k4f, k4b, k5f, k5b,
			k6f, k6b, k7f, k7b, k8f, k8b, k9f, k9b, k10f, k10b, 
			k11f, k11b, k12f, k12b, k13f, k13b, k14f, k14b, k15f, k15b, 
			k1OH, k2OH, k3OH, k1O2, k2O2, k3O2 }; // used for TSoot::fRateCoefficients
#ifdef NEWSURFGROWTH
enum SootRateCoeffs			{ ks7f, ks7b, ks8f, ks8b, ks9f, ks9b, ks10f, ks10b, 
			ks111, ks112, ks12, ks13f, ks13b, kLastSootRateCoeff }; // used for TSoot::fSootRateCoeffs
#else
enum SootRateCoeffs			{ ks8f, ks8b, ks9, ks10f, ks10b, 
			ks11, ks12, kLastSootRateCoeff }; // used for TSoot::fSootRateCoeffs
#endif
enum ReducedReactionCoeffs	{ k01, k10, k12, k21, k23, k32, k34, k43, 
			k45, k54, k56, k65, k51, k15,
			k16, k31, k63, k25, k42, k53, kLastRedReacCoeff }; // used for TSoot::fRedRateCoeffs
enum Polymere				{ kFirstPoly, kRestPoly };
enum FracMoments			{ km1_2, k01_2, k03_2, k05_2,
			k19_6, k17_6, k13_6, k11_6, 
			k07_6, k05_6, k01_6, km1_6, 
			km1_24, k23_24, k47_24, k71_24, 
			km1_3, k02_3, k05_3, k08_3, kLastFracMoment }; // used for TSoot::fFracMoms
enum Phi					{ k000, k100, k011, k111, k012, k112, kh00, kh11, 
							kh12, kLastPhi }; // used for TSoot::fPhi
/*enum SmallPAH				{ kA1, kA1M, kA1C2H, kA1C2HM, kA1C2HS, kA1C2HAC, kA2MX*/
/*			, kA2, kA2R5, kA2R5M, kA2R5C2H, kA2R5C2HS, kANC2HAC, kA3R5M*/
/*			, kA3R5, kA3R5AC, kNSmallSoot };*/
enum SootSourcTerms		{ kNuc = 0, kCoag, kCond, kSG, kOx, kTherPhor, kNSootSources };

class TSoot {
public:
	TSoot( TInputDataPtr input ) : fNPAHMolecules( 6 ),fNStages( 50 ),fNMinRings(4)
									,fNSootMoments(input->fNSootMoments)
									,fNPAHMoments(input->fNSootMoments),fAlpha(0.9)
									,fBeta(0.01),fMolarMassSoot(24.0),fSootDensity(1800.0)
									,fMolarMassDiff( 74.02 ) { InitTSoot( input ); };
	~TSoot( void );

	int				GetNPAHReactions( void ) { return fA->len; };
	int				GetNSootReactions( void ) { return fASoot->len; };
	void			ComputePolymereConcs( Double *Y, Double temp, Double density, Double *molarMass
						, Double **Pij, Double *sumPi, Double *pahMoments, Double *moments, TReactionPtr reaction );
	void			UpdateProductionRates( TSpeciesPtr species, TReactionPtr reaction, Double *prodRate
			, Double density, Double *Y, Double temp, Double *sumPi, Double *moments, Double *wPAH, Double mixMolarMass );
	void			UpdateSoot( TReactionPtr reaction, TSpeciesPtr species
			, Double *moments, Double temp, Double *Y, Double density, Double mixMolarMass );
	void			SetMomentsOffset( int offset ) { fOffsetMoments = offset; };
	int				GetOffsetSootMoments( void ) { return fOffsetMoments; };
	int				GetNPAHMoments( void ) { return fNPAHMoments; };
	int				GetNSootMoments( void ) { return fNSootMoments; };
	Double			GetSootRadRossCoeff( Double temp, Double *moments, Double cutoff );
	Double			GetSootRadiation( Double temp, Double *moments );
	Double			GetGamma( void ) { return fGamma; };
	Double			GetCoagFact( void ) { return fCoagFact; };
	Double			GetCCSootStar( void ) { return fCSootStar; };
	Double			GetSCoag( Double *mom, Double temp );
	Double			GetLewis1( void ) { return fLewis1; };
	Double			FracMom2( Double fracIndex, Double *moments );
	Flag			WithNucleation( void ) { return fNucleation; };
	Flag			WithCondensation( void ) { return fCondensation; };
	Flag			WithCoagulation( void ) { return fCoagulation; };
	Flag			WithSurfaceGrowth( void ) { return fSurfaceGrowth; };
	Flag			WithSurfaceOxidation( void ) { return fSurfaceOxidation; };
	Flag			WithThermoPhoresis( void ) { return fThermoPhoresis; };
	Flag			WithSootRadiation( void ) { return fSootRadiation; };
	Flag			WithUpdateProdRates( void ) { return fSootUpdateProdRates; };
	Flag			WithSizeDepDiff( void ) { return fSizeDepDiff; };
	Flag			WithSurfDepCoag( void ) { return fSurfDepCoag; };
	Double			NucleationNow( int i, Double temp, Double *Y, Double rho, Double *molarMass );
	Double			NucleationNew( int i, Double temp, Double *pahMoments );
	Double			SourceCoagulationNew( int i, Double temp, Double *moments );
	Double			SourceSurfDepCoag( int i, Double *mom, Double *Y, Double temp
					, Double density, Double *molarMass );
//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//
/*	Double			SourceCondensationNew( int i, Double temp, Double *pahMoments, Double *moments );*/
	Double			SourceCondensationNew( int i, Double temp, Double *pahMoments, Double *moments
                        , Double *Y, Double density, Double *molarMass );
//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//
	Double			SourceSurfGrowthNew( int i, Double *Y, Double density, Double *molarMass );
	Double			SourceSootOxidationNew( int i, Double *Y, Double density, Double *molarMass );
	Double			SourceSurfGrowthNew( int i, Double *moments, Double *Y, Double density, Double *molarMass );
	Double			SourceSootOxidationNew( int i, Double *moments, Double *Y, Double density, Double *molarMass );
	void			MOverRhoToM( Double *momSource, Double *momDest, int nMom
						, Double *Y, Double temp, Double pressure
						, Double *molarmass, int nSpeciesIn, TPropertiesPtr props );
	void			MOverRhoToM( Double *momSource, Double *momDest, int nMom, Double rho );
	Double			GetAlphaI2( int i, Double fracIndex );
	int				PrintPAHReactionEquation( int number, TSpeciesPtr species, FILE *fp );
	void			PrintPAHReactions( TFlamePtr flame, TSpeciesPtr species );
	void			PrintSootReactions( TFlamePtr flame, TSpeciesPtr species );
	void			ComputeCSootStar( Double *k, Double *Y, Double density, Double *molarMass, Double mixMolarMass );
	VectorPtr		GetSootRateCoeffs( void ) { return fSootRateCoeffs; };
	void			ComputeSootRateCoeffs( Double *k, Double temp, TReactionPtr reaction );
	Double			GetSootOxCoeff( Double *Y, Double density, Double *molarMass );
	
protected:
	void 			InitTSoot( TInputDataPtr input );
	void			FillTSoot( TInputDataPtr input );
	void			CheckPAHReactions( void );
	void			CheckSootReactions( void );
	void			ComputeRateCoefficients( Double *k, Double temp, TReactionPtr reaction );
	void			ComputeRedRateCoeffs( Double *K, Double *Y, Double *k, Double density, Double *molarMass );
	void			ComputeF( void );
	void			ComputeZ( void );
	void			CalcMomentsNewPoly( Double *pahMoments, Double *sumPi, Double *p0 );
	Double			ComputeNewPoly( Double *Y, Double density, Double *molarMass, Double *p0 );
	Double			ComputeN( int poly );
	void			ComputeAB( int poly );
	Double			ComputeXInf( void );
	void			ComputeDelta( void );
	Double			ComputeGamma( int poly );
	Double			ComputeKco( Double temp, Double *moments, Double *pahMoments );
	Double			GetBeta( Double temp, Double *pahMoments );
	Double			GetC( Double temp );
	Double			GetCPAH( Double temp );
	void			ComputeP1j( Double **Pij );
	void			ComputePij( Double **Pij );
	void			ComputePAHMoments( Double *pahMoments, Double **Pij );
	Double			FractionalMoment( Double fracIndex, Double *moments );
	Double			DerivFractionalMoment( int index, Double fracIndex, Double *moments );
	Double			GetAlphaI( int i, Double fracIndex );
	void 			CheckSolution( Double *theMom );
	void			ComputeFractionalMoments( Double *moments );
	void			ComputePhi( Double temp );
	Double			SourceCoagulation( int i );
	Double			GetAlphaI2( int i, Double fracIndex, int first, int second );
	Double			GetDerivFracMom( int l, FracMoments nMom, Double *moments );
	Double			GetDerivFracMom( int l, Double which, Double *moments );
	Double			GetDerivPhi( int l, Phi which, Double *moments );
	Double			GetDerivPhiCond( int l, int x, int y, Double temp, Double *moments1, Double *moments2  );
	Double			FractionalMoment2( Double fracIndex, int first, int second, Double *moments );
	Double			SourceCondensation( int i, Double temp, Double *pahMoments, Double *moments );
	Double			GetBetaCond( Double temp, Double *moments );
	void			InitSootReactions( void );
	void			ComputeSootReactionRates( Double *w, Double *Y, Double *moments, Double density, Double *molarMass, Double mixMolarMass );
	Double			GetSurfGrowthCoeff( Double *Y, Double density, Double *molarMass );
	Double			GetSurfGrowthCoeffForw( Double *Y, Double density, Double *molarMass );
	Double			GetSurfGrowthCoeffBackw( Double *Y, Double density, Double *molarMass );
	Double			SourceSurfGrowth( int i, Double *Y, Double density, Double *molarMass );
	Double			SourceSootOxidation( int i, Double *Y, Double density, Double *molarMass );
	void			ComputeFractionalMomentsPAH( Double *moments );
	void			ComputePhiPAH( Double temp );
	Double			Nucleation( int i, Double temp, Double *pahMoments );
	Double			GetPhi( int x, int y, Double temp, Double *moments ) 
						{ return GetPhi( x, y, temp, moments, moments ); };
	Double			GetPhiPAH( int x, int y, Double temp, Double *moments ) 
						{ return GetPhiPAH( x, y, temp, moments, moments ); };
	Double			GetPhi( int x, int y, Double temp, Double *moments1, Double *moments2 );
	Double			GetPhiPAH( int x, int y, Double temp, Double *moments1, Double *moments2 );
//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//
	Double			GetPhiPAHFROMA4( Double x, Double temp, Double *moments );	
//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//i//
/*	void			InitSmallPAHIndTab( TInputDataPtr input );*/
/*	void			SetSmallPAH( SmallPAH which, char *name, Double nOfRings, TInputDataPtr input );*/
/*	void			ComputeSmallPAHMoments( Double *Y, Double density, TSpeciesPtr species );*/
	void			PrintPAHReactions( int number, TSpeciesPtr species, FILE *fp );

	char			**fLabels;
	char			*fPAHSymbol;
	const int		fNPAHMolecules;
	const int		fNStages;
	const int		fNSootMoments;
	const int		fNPAHMoments;
	const int		fNMinRings;
	const Double	fMolarMassDiff; // M(A_{i+2}) - M(A_{i}) = M( C6H2 ) = 74.02 kg/kmole
	VectorPtr		fA;	// [m, s, kmole, K]					// input; length is nPAHReactions
	VectorPtr		fN;										// input; length is nPAHReactions
	VectorPtr		fEOverRgas;	// [J / kmole]				// input; length is nPAHReactions
	IntVectorPtr	*fSpeciesNumber;						// input; see fNu; 
	VectorPtr		*fNu;									// nu_product is less than zero; last element is fNu[nOfReactions][speciesPerReaction]

	IntVectorPtr	fPAHSpeciesIndex;	// contains indice of pah's in the specieslist; length is fNPAHMolecules+1
	VectorPtr		fRateCoefficients;	// [m, s, kmole]	// k_i; length is nPAHReactions
	VectorPtr		fRedRateCoeffs;	// [m, s, kmole]		// K_i; length is kLastRedReacCoeff
															// K[k01] contains r0 = k0f * [A3R5AC]
	VectorPtr		fZ;										// Z_i; length is fNPAHMolecules + 1; Z[0] not used -> Z[1] means Z1
	VectorPtr		fDenom;										// Z_i; length is fNPAHMolecules + 1; Z[0] not used -> Z[1] means Z1
	Double			fZ53, fZ52, fZ51;
	Double			fKco;									
	Double			fNCoeff;										
	Double			fF121;										
	Double			fF151;										
	VectorPtr		fF;										// F_i; length is kLastRedReacCoeff
	MatrixPtr		fACoeff;		// workspace for the coefficients for the
	MatrixPtr		fBCoeff;		// computation of PAH concentrations; last element is fACoeff[kRestPoly][fNPAHMolecules+1]
	MatrixPtr		fBOHCoeff;		// fACoeff, fBCoeff and fBOHCoeff are using enum Polymere
	VectorPtr		fDelta;			// P_ij = delta_j * gamma^(i-2) * P_15; length is fNPAHMolecules+1
	Double			fGamma;			// P_ij = delta_j * gamma^(i-2) * P_15;
	Double			fXInf;

	VectorPtr		fFracMom;		// contains fractional moments; uses enum FracMoments; length is kLastFracMoment
	VectorPtr		fPhi;			// contains Phi; uses enum Phi; length is kLastPhi

	VectorPtr		fFracMomPAH;		// contains fractional moments; uses enum FracMoments; length is kLastFracMoment
	VectorPtr		fPhiPAH;			// contains Phi; uses enum Phi; length is kLastPhi

	char			*fSootSymbol;
	char			**fSootLabels;
	VectorPtr		fASoot;	// [m, s, kmole, K]					// input; length is nSootReactions
	VectorPtr		fNSoot;										// input; length is nSootReactions
	VectorPtr		fEOverRSoot;	// [J / kmole]				// input; length is nSootReactions
	IntVectorPtr	*fSpecNumSoot;								// input; see fNu; 
	VectorPtr		*fNuSoot;									// nu_product is less than zero; last element is fNu[nOfReactions][speciesPerReaction]
	VectorPtr		fSootRateCoeffs;	// [m, s, kmole]	// k_i; length is kLastSootRateCoeff 
	VectorPtr		fSootReactionRate;	//// k_i; length is kLastSootRateCoeff
	Double			fCSootStar;			// used for the calculation of the surface growth rate
	Double			fFk10;			// used for the calculation of the surface growth rate
	Double			fThirdBodyConc;
	const Double	fAlpha;				// used for the calculation of the surface growth rate
	const Double	fBeta;				// used for the calculation of the soot oxidation rate
	
	int				fOffsetMoments;	//	offset of the moments in the solution vector of the solver

/*	int				*fSmallPAHIndTab;	// length is kNSmallSoot*/
/*	Double			*fSmallPAHRings;	// length is kNSmallSoot*/
/*	Double			*fSmallPAHNumDens;	// last element is fSmallPAHNumDens[kNSmallSoot]*/
/*	Double			*fSmallPAHMoments;	// last element is fSmallPAHMoments[fNSootMoments]*/
	const Double	fMolarMassSoot;	//	24.0 [kg/kmole], difference of two size classes: C2
	const Double	fSootDensity;	//	1800.0 [kg/m^3]

	Flag			fNucleation;
	Flag			fCondensation;
	Flag			fCoagulation;
	Flag			fSurfaceGrowth;
	Flag			fSurfaceOxidation;
	Flag			fThermoPhoresis;
	Flag			fOHPAHOxidation;
	Flag			fO2PAHOxidation;
	Double			fCoagFact;
	Double			fLewis1;
	Flag			fSootRadiation;
	Flag			fSootUpdateProdRates;
	Flag			fSizeDepDiff;
	Flag			fSurfDepCoag;

	int				f_H;
	int				f_H2;
	int				f_O2;
	int				f_OH;
	int				f_CO;
	int				f_CH;
	int				f_CHO;
	int				f_H2O;
	int				f_C2H;
	int				f_C2H2;
	int				f_A3R5AC;
	int				f_A3R5M;
	int				f_A1;
	int				fFirstPAH;
};
typedef TSoot *TSootPtr;

class T0DSoot : public TSoot {
public:
	T0DSoot( TInputDataPtr input ) : TSoot( input ) { InitT0DSoot( input ); };
	~T0DSoot( void );

	MatrixPtr		GetPij( void ) { return fPij; };
	VectorPtr		GetSumPi( void ) { return fSumPi; };
	VectorPtr		GetMoments( void ) { return fMoments; };
	VectorPtr		GetPAHMoments( void ) { return fPAHMoments; };
	VectorPtr		GetReactionRate( void ) { return fReactionRate; };
	
private:
	void 			InitT0DSoot( TInputDataPtr input );

	MatrixPtr		fPij;			// vector of PAH concentrations
									// last element is fPij[nStages][fNPAHMolecules+1]
	VectorPtr		fSumPi;				// last element is fSumPi[fNPAHMolecules]
	VectorPtr		fPAHMoments;		// last element is fMoments[fNPAHMolecules]
	VectorPtr		fMoments;			// last element is fMoments[fNSootMoments]
	VectorPtr		fReactionRate;	// [kmole / (m^3*s)]	// w_i; length is nPAHReactions
};
typedef T0DSoot *T0DSootPtr;

Double SootCoagFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag = FALSE );
Double SootCondFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag = FALSE );
//Double SootCondSmallPAHFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag = FALSE );
Double SootOxidationFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag = FALSE );
class T1DSoot : public TSoot {
friend Double SootCoagFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag );
friend Double SootCondFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag );
//friend Double SootCondSmallPAHFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag );
friend Double SootOxidationFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag );
public:
	T1DSoot( TInputDataPtr input ) : TSoot( input ) { InitT1DSoot( input ); };
	~T1DSoot( void );

	TensorPtr		GetPij( void ) { return fPij; };
	MatrixPtr		GetSumPi( void ) { return fSumPi; };
	MatrixPtr		GetMoments( void ) { return fMoments; };
	MatrixPtr		GetPAHMoments( void ) { return fPAHMoments; };
	MatrixPtr		GetReactionRate( void ) { return fReactionRate; };
	VectorPtr		GetSootDiff( void ) { return fSootDiff; };
	void			UpdateJacobian( T1DFlamePtr flame );
	void			PrintPAHj( int i, TNewtonPtr bt );
	void			PrintPAHOneStage( int i, TNewtonPtr bt );
	void			PrintPAHMoments( TNewtonPtr bt );
	void			PrintSootInfo( TNewtonPtr bt );
	void			PrintSurfDepCoag( TNewtonPtr bt, T1DFlamePtr flame );
	void			PrintSurfDepCoag( T1DFlamePtr flame, NodeInfoPtr nodeInfo, FILE *fp );
	void			UpdateDimensions( int len );
	void			UpdateSolution( T1DFlamePtr flame, MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
	void			UpdateSolution( T1DFlamePtr flame, Double *y, int gridPoint );
	void			SolutionToSolver( Double **y );
	void			SaveSolution( void );
	void			RestoreSolution( void );
	void			FillJacobi( T1DFlamePtr flame, NodeInfoPtr nodeInfo, CoordType coordinate );
	void			FillRHS( T1DFlamePtr flame, NodeInfoPtr nodeInfo, CoordType coordinate );
	void			FillJacSootRadiation( Double fact, T1DFlamePtr flame, NodeInfoPtr nodeInfo );
	void			PostIter( T1DFlamePtr flame );
	void			PrintFlameletFile( int gridPoints, T1DFlamePtr flame, FILE *fp );
	void			PrintRHSSoot( TNewtonPtr bt, T1DFlamePtr flame );
	void			PrintFracMom( T1DFlamePtr flame );
	void			PrintPhi( T1DFlamePtr flame );
	void			PrintDiffusivity( T1DFlamePtr flame );
	void			ComputeDiffusivity( T1DFlamePtr flame );
	void			PrintPathsOfAcetylene( T1DFlamePtr flame );
	void			FillSource( Double *source, T1DFlamePtr flame );
	TensorPtr		GetDMdx( void ) { return fDMdx; };
	MatrixPtr		GetSource( void ) { return fSource; };
	VectorPtr		GetSourceMod( void ) { return fSourceMod; };
	
private:
	void 			InitT1DSoot( TInputDataPtr input );
	void			FillJacSootConvection( int i, T1DFlamePtr flame, NodeInfoPtr nodeInfo, Flag velocityPositive );
	void			FillJacSootDiffusion( int r, Double constCoeff, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo );
	Double			SootConvection( int i, T1DFlamePtr flame, NodeInfoPtr nodeInfo, Flag velocityPositive );
	Double			SootDiffusion( int r, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo );
	void 			CheckSolution( int nGridPoints, Double **y );
	void			FillJacSootDiffusionNew( int r, Double constCoeff, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo );
	Double			SootDiffusionNew( int r, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo );
	Double			FillJacSourceCoagulation( int i, int l, Double *moments );
	void			PrintRHSSoot( T1DFlamePtr flame, NodeInfoPtr nodeInfo, FILE *fp );
	void			FillJacTheDiffusion( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame, Flag sign );
	Double			TheDiffusion( int i, NodeInfoPtr nodeInfo, T1DFlamePtr flame );
	Double			SootThermoPhoresis( int r, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo );
	void			FillJacSourceCondensation( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame );
	Double			FillJacCondensationNew( int l, int i, Double temp, Double *pahMom, Double *mom );
	void			FillJacSurfGrowth( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame );
	void			FillJacSootOxidation( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame );
	void			FillJacSurfGrowthNew( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame );
	void			FillJacSootOxidationNew( int i, Double constCoeff, NodeInfoPtr nodeInfo, T1DFlamePtr flame );
	void			FillJacSootThermoPhoresis( int r, Double constCoeff, CoordType coordinate, T1DFlamePtr flame, NodeInfoPtr nodeInfo );
	void			PrintPathsOfAcetylene( T1DFlamePtr flame, Double x, FILE *fp );

	MatrixPtr		fReactionRate;	// [kmole / (m^3*s)]	// w_i; length is nPAHReactions
	TensorPtr		fPij;			// vector of PAH concentrations
									// last element is fPij[nGridPoints][nStages][fNPAHMolecules+1]
	MatrixPtr		fSumPi;				// last element is fSumPi[nGridPoints][fNPAHMolecules]
	MatrixPtr		fPAHMoments;		// last element is fMoments[nGridPoints][fNPAHMolecules]
	MatrixPtr		fMoments;			// last element is fMoments[nGridPoints][fNSootMoments]
	MatrixPtr		fSavedMoments;		// last element is fMoments[nGridPoints][fNSootMoments]
	VectorPtr		fSootDiff;		//	length is nOfGridPoints+2
	TensorPtr		fDMdx;									// last element is fDMdx->tensor[nOfGridPoints][nMoments][nOfSpecies+1+nMoments]
															// 	d(M_i/rho)/dy_j = fDMdx->tensor[nOfGridPoints][j][i]
															// contains the derivatives with respect to massfractions and temperature and Moments/rho
	MatrixPtr		fSource;
	VectorPtr		fSourceMod;
};
typedef T1DSoot *T1DSootPtr;

//*************************************************************************************************************
class TSpecies {
// TFlame contains a function UpdateThermoProps for the complete calculation of all species and mixture properties
// and the calculation of the production rates

//	compute all properties for steady states except for 
//	steady state species have no storage for fNOfUsedReactions, fUsedReactions, fNu, 
//		fProductionRate

public:
	TSpecies( TInputDataPtr input, TReactionPtr reaction ) { InitSpecies( input, reaction ); };
	~TSpecies( void );

	int				GetNOfSpecies( void ) { return fKOverEps->len; };
	int				GetNSteadyStates( void ) { return fNSteadyStates; };
	int				GetNSpeciesInSystem( void ) { return fKOverEps->len - fNSteadyStates; };
	VectorPtr		GetMolarMass( void ) { return fMolarMass; };
	VectorPtr		GetLewisNumber( void ) { return fLewisNumber; };
	const char		*GetLewisNumberFile( void ) { return fLewisNumberFile; };
	VectorPtr		*GetNu( void ) { return fNu; };
	IntVectorPtr	GetNOfUsedReactions( void ) { return fNOfUsedReactions; };
	IntVectorPtr	*GetUsedReactions( void ) { return fUsedReactions; };
	//	DoublePtr		***GetReactionRate( void ) { return fReactionRate; };
	char			**GetNames( void ) { return fNames; };
	Flag			IsConstantLewisNumber( void ) { return fConstLewisNumber; };
	Flag			IsSteadyState( int i ) { return fIsSteadyState[i]; };
	int				FindSpecies( const char *species );
	void			WriteRoggsSymbolsData( void );
	Double			ComputeDCpiDT( int index, Double temp );
//	Double			CompOnedRhoRhoDkdY( int speciesIndexY, int speciesIndexD, Double *Y, TFlameNodePtr flameNode );
//	Double			CompOnedRhoDkdY( int speciesIndexY, int speciesIndexD, Double *Y, Double rho, Double M, Double *D, Double **invDij );
	void			ComputeProductionRates( Double *productionRate, Double *reactionRate );
	Double			GetPosNegProductionRate( int which, Double *reactionRate, Flag pos );
	void			ComputeNoProductionRates( Double *productionRate );
	void			ComputeProductionRate( int i, Double &productionRate, Double *reactionRate );
	void			ComputeDeltaI( Double *deltaI, Double *Y, Double *viscosity );
	void			CompProdCons( Double *source, Double *sink, Double *reacRate );
	void			CompHeatRelease( Double *heatRel, Double *prodRate, Double *enthalpy );
	int				CompDetailedProdRate( int i, Double *prodRate, Double *reactionRate, TReactionPtr reaction );
	void			ReadLewisNumbers( const char *file, VectorPtr Le );
	MatrixPtr		GetCoeffHigh( void ) { return fCoeffHigh; };
	MatrixPtr		GetCoeffLow( void ) { return fCoeffLow; };
	void			ComputeTheProductionRates( Double *productionRate, Double *reactionRate
						, Double temp, Double pressure, Double *c
						, Double *k, Double *M );
	Double			**GetSqrtMjOverMi( void ) { return sqrtMjOverMi; };
	Double			**GetSqrtMiOverMjPlusOne( void ) { return sqrtMiOverMjPlusOne; };

protected:
	void			InitSpecies( TInputDataPtr input, TReactionPtr reaction );	
	void 			CalcUsedReactions( TReactionPtr reaction );
	void			FillSpecies( TInputDataPtr input, TReactionPtr reaction );
	Double 			omega_mu( Double t );
	Double 			omega_D( Double t );
	Double 			CompDeltaI( int i, int nSpecies, Double *M, Double *mu, Double *Y );
	Double			CompDeltaIOpt( int i, int nSpeciesInSystem, Double **GijOverWj, Double *M, Double *Y );
	Flag			EmptyLine( char *s );

	int				fNSteadyStates;					//	number of steady state species
	IntVectorPtr	fNOfUsedReactions;				//  array of length nOfSpecies
	IntVectorPtr	*fUsedReactions;				//  last element is fUsedReactions[nOfSpecies]->vec[fNOfUsedReactions]
	VectorPtr		*fNu;							//  last element is fNu[nOfSpecies]->vec[fNOfUsedReactions]
													//	nu of educts are > 0
	Flag			fConstLewisNumber;				//
	char			*fLewisNumberFile;
	VectorPtr		fLewisNumber;					//	contains constant Lewis number for every species
													//	initialized in TSpecies::FillSpecies	
	Flag			*fIsSteadyState;				//	signals wether the species is steady state or not
													//	length is nOfSpecies
	
	// constants needed for the calculation of the properties 
	char			**fNames;
	MatrixPtr		fCoeffLow; 			// last element is fCoeffLow[nOfSpecies][7]
	MatrixPtr		fCoeffHigh;			// last element is fCoeffHigh[nOfSpecies][7]
	MatrixPtr		fHCoeffLow; 		// last element is fCoeffLow[nOfSpecies][7]
	MatrixPtr		fHCoeffHigh;		// last element is fCoeffHigh[nOfSpecies][7]
	VectorPtr		fMolarMass;			//  [kg / kmol]
	VectorPtr		fKOverEps;			//	length is nOfSpecies
	VectorPtr		fSigma;				//	length is nOfSpecies
	VectorPtr		fMuCoeff;			//	length is nOfSpecies
	MatrixPtr		fDCoeff;			// fDCoeff[i][j] is the same as fDCoeff[j][i]
	MatrixPtr		fOmegaCoeff;		// fOmegaCoeff[i][j] is the same as fOmegaCoeff[j][i]
	MatrixPtr		fInvDij;			// last element is fInvDij->mat[nOfSpecies][nOfSpecies]
	Double			**sqrtMjOverMi;			//	contains sqrt( Mj/Mi ); last element is sqrtMjOverMi[nOfSpecies][nOfSpecies]
	Double			**sqrtMiOverMjPlusOne;	//	contains 0.3535534 / sqrt( Mi/Mj + 1 ); last element is sqrtMiOverMjPlusOne[nOfSpecies][nOfSpecies]
	Double			fTThermoLimit;
};


//*************************************************************************************************************
class T0DSpecies : public TSpecies {
public:
	T0DSpecies( TInputDataPtr input, TReactionPtr reaction ) 
			: TSpecies( input, reaction )
		{ InitT0DSpecies( input ); };
	~T0DSpecies( void );

	VectorPtr	GetHeatCapacity( void ) { return fHeatCapacity; };
	VectorPtr	GetConductivity( void ) { return fConductivity; };
	VectorPtr	GetViscosity( void ) { return fViscosity; };
	VectorPtr	GetEnthalpy( void ) { return fEnthalpy; };
	VectorPtr	GetProductionRate( void ) { return fProductionRate; };
	VectorPtr	GetDeltaI( void ) { return fDeltaI; };
	Double		&GetCurrTemp( void ) { return fCurrTemp; };

	Flag	ComputeSpeciesProperties( Double temp );
	void	ComputeSpeciesProperties( int number, Double temp );
	void	PrintSpecies( T0DFlamePtr flame );
	void	PrintSpecies( int number, T0DFlamePtr flame, FILE *fp );

private:
	void			InitT0DSpecies( TInputDataPtr input );	

	VectorPtr		fHeatCapacity; // [J / (kg*K)]	//	last element is fHeatCapacity[nOfSpeciesIn]
	VectorPtr		fEnthalpy; // [J / kg]			//	last element is fEnthalpy[nOfSpeciesIn]
	VectorPtr		fViscosity; // [J / (kg*K)]	//	last element is fViscosity[nOfSpeciesIn]
	VectorPtr		fConductivity; // [J / kg]			//	last element is fConductivity[nOfSpeciesIn]
	VectorPtr		fProductionRate;		//  m_i	//	last element is fProductionRate[nOfSpeciesIn]
	VectorPtr		fDeltaI;						//  used by TProperties::CompMixtureProps and TSpecies::Compute_DTherm
	Double			fCurrTemp;
};
typedef T0DSpecies *T0DSpeciesPtr;

#ifndef ZEROD
//*************************************************************************************************************
class T1DSpecies : public TSpecies {
//	steady state species have no storage for fNOfUsedReactions, fUsedReactions, fNu, 
//		fProductionRate
public:
	T1DSpecies( TInputDataPtr input, T1DReactionPtr reaction ) 
			: TSpecies( input, reaction )
	{ 
		InitT1DSpecies( input ); 
	};
	~T1DSpecies( void );

	TensorPtr		GetBinDij( void ) { return fBinDij; };
	TensorPtr		GetGijOverWj( void ) { return fGijOverWj; };
	TensorPtr		GetOmegaDOverDCoeff( void ) { return fOmegaDOverDCoeff; };
	TensorPtr		GetDThermConst( void ) { return fDThermConst; };
	MatrixPtr		GetProductionRate( void ) { return fProductionRate; };
	MatrixPtr		GetViscosity( void ) { return fViscosity; };
	MatrixPtr		GetHeatCapacity( void ) { return fHeatCapacity; };
	MatrixPtr		GetConductivity( void ) { return fConductivity; };
	MatrixPtr		GetEnthalpy( void ) { return fEnthalpy; };
	MatrixPtr		GetDiffusivity( void ) { return fDiffusivity; };
	MatrixPtr		GetDiffTherm( void ) { return fThermoDiffusivity; };
	MatrixPtr		GetDeltaI( void ) { return fDeltaI; };
	void 			PrintSpecies( int k );									// prints all
	void 			PrintSpecies( int number, int gridPoint, FILE *fp );	// prints one
	void			PrintProductionRate( TNewtonPtr bt );
	void			PrintProdCons( TNewtonPtr bt, T1DFlamePtr flame );
	void			PrintDiffusionCoeffs( TNewtonPtr bt );
	void			PrintProdRateTerms( const char *name, T1DFlamePtr flame );
	void			PrintProdRateTerms( int i, T1DFlamePtr flame );
	Flag	 		ComputeSpeciesProperties( TFlameNodePtr flameNode, Double pressure, Double temp ); 				// computes properties of all species
	void 			ComputeSpeciesProperties( TFlameNodePtr flameNode, Double temp, Double pressure, int number );	// computes properties of one species
	void 			Compute_D( TFlameNodePtr flameNode, Double temp, Double *Y, Double pressure, Flag newTemp );	// computes diffusion coefficients
	void			Compute_D( TFlameNodePtr flameNode );
	void			ComputeDeltaIOpt( TFlameNodePtr flameNode, Double *Y, Double **GijOverWj, Flag newTemp );
	void			CompBinDiffCoeff( TFlameNodePtr flameNode, Double temp, Double *Y, Double press );
	void			Compute_DTherm( TFlameNodePtr flameNode, Flag calcNewConst );
	Double			ComputeOneDiffusionCoeff( int number, Double *Y, Double *invDij, Double mixMolarMass );
	MatrixPtr		GetSavedY( void ) { return fSavedY; };
	MatrixPtr		GetSavedDeltaiY( void ) { return fSavedDeltaiY; };
	MatrixPtr		GetSumDiff( void ) { return fSumDiff; };
	VectorPtr		GetTempProp( void ) { return fTempProp; };
	VectorPtr		GetPressureProp( void ) { return fPressureProp; };

private:
	void			InitT1DSpecies( TInputDataPtr input );	

	MatrixPtr		fViscosity; // [kg / (m*s)]		//	last element is fViscosity[nGridPoints][nOfSpecies]
	MatrixPtr		fHeatCapacity; // [J / (kg*K)]	//	last element is fHeatCapacity[nGridPoints][nOfSpecies]
	MatrixPtr		fConductivity; // [W / (m*K)]	//	last element is fConductivity[nGridPoints][nOfSpecies]
	MatrixPtr		fEnthalpy; // [J / kg]			//	last element is fEnthalpy[nGridPoints][nOfSpecies]
	MatrixPtr		fDiffusivity; // [m^2 / s]		//	last element is fDiffusivity[nGridPoints][nOfSpecies]
	MatrixPtr		fThermoDiffusivity; // [m^2 / s]//	last element is fDiffusivity[nGridPoints][nOfSpecies]
	MatrixPtr		fProductionRate;		//  m_i	//	last element is fViscosity[nGridPoints][nOfSpecies - nSteadyStates]
	TensorPtr		fBinDij;			// last element is fBinDij->tensor[nGridPoints][nOfSpeciesIn][nOfSpeciesIn]
	TensorPtr		fGijOverWj;	//	contains 0.3535534 / sqrt( Mi/Mj + 1 ); last element is sqrtMiOverMjPlusOne[nOfSpecies][nOfSpecies]
	TensorPtr		fOmegaDOverDCoeff;
	TensorPtr		fDThermConst;

	MatrixPtr		fSavedY; 
	MatrixPtr		fSavedDeltaiY; 
	MatrixPtr		fSumDiff; 
	MatrixPtr		fDeltaI;						//  used by TProperties::CompMixtureProps and TSpecies::Compute_DTherm
	VectorPtr		fTempProp;
	VectorPtr		fPressureProp;
};
typedef T1DSpecies *T1DSpeciesPtr;
#endif // ZEROD

//*************************************************************************************************************
class TRadiation {
public:
	TRadiation( TInputDataPtr input ) : fExpFuncCoeffCO2(-8.888e-4),fExpFuncCoeffH2O(-1.546e-3)
					,fAlphaCoeffCO2(46.241e-5),fAlphaCoeffH2O(22.6e-5){ InitRadiation( input ); };
	~TRadiation( void );

	void 		ComputeRadiationOnePoint( Double *radiation, Double temp, Double *Y, Double *molarMass, Double	density );
protected:
	void 	InitRadiation( TInputDataPtr input );

	int			fH2OIndex;
	int			fCO2Index;
	int			fCH4Index;
	int			fCOIndex;
	const Double	fExpFuncCoeffCO2; // = -8.888e-4
	const Double	fExpFuncCoeffH2O; // = -1.546e-3
	const Double	fAlphaCoeffCO2; // = 46.241e-5
	const Double	fAlphaCoeffH2O; // = 22.6e-5
};
typedef TRadiation *TRadiationPtr;

//*************************************************************************************************************
class T0DRadiation : public TRadiation {
public:
	T0DRadiation( TInputDataPtr input ) : TRadiation( input ) { InitT0DRadiation( input ); };
	~T0DRadiation( void ) {};

	Double		GetRadiation( void ) { return fRadiation; };
	void		SetRadiation( Double temp, Double *Y, Double *molarMass, Double density ) 
					{ ComputeRadiationOnePoint( &fRadiation, temp, Y, molarMass, density ); };
private:
	void 	InitT0DRadiation( TInputDataPtr /*input*/ ) {};

	Double	fRadiation;
};
typedef T0DRadiation *T0DRadiationPtr;

#ifndef ZEROD
//*************************************************************************************************************
class T1DRadiation : public TRadiation {
public:
	T1DRadiation( TInputDataPtr input ) : TRadiation( input ) { InitT1DRadiation( input ); };
	~T1DRadiation( void );

	VectorPtr	GetRadiation( void ) { return fRadiation; };
	void		FillJacRadiation( Double coeff, T1DFlamePtr flame, NodeInfoPtr nodeInfo );
private:
	void 	InitT1DRadiation( TInputDataPtr input );

	VectorPtr	fRadiation;		//	length is nOfGridPoints
};
typedef T1DRadiation *T1DRadiationPtr;
#endif // ZEROD

//*************************************************************************************************************
class TProperties {	// contains properties of the gas mixture
// TFlame contains a function for the complete calculation of all species and mixture properties
public:
	TProperties( TInputDataPtr input, TSpeciesPtr species ) 
			{ InitProperties( input, species ); };
	~TProperties( void );

	TRadiationPtr 	GetRadiation( void ) { return fRadiation; };
	void 			ComputeMixtureMolarMass( Double &mixMolarMass, Double *Y, Double *molarMass, int nSpecies );
	void			ComputeDCpDT( Double &dCpdT, Double *Y, Double temp, TSpeciesPtr species );
/*	// 0D stuff*/
/*	void			CompMixtureProps( Double &mixHeatCapacity, Double &density, Double *heatCapacity, Double temp, Double mixMolarMass, Double pressure );*/

protected:
	void 			InitProperties( TInputDataPtr input, TSpeciesPtr species );
/*	void			FillProperties( TSpeciesPtr species );*/
	
	TRadiationPtr	fRadiation;			
/*	Double			**sqrtMjOverMi;			//	contains sqrt( Mj/Mi ); last element is sqrtMjOverMi[nOfSpecies][nOfSpecies]*/
/*	Double			**sqrtMiOverMjPlusOne;	//	contains 0.3535534 / sqrt( Mi/Mj + 1 ); last element is sqrtMiOverMjPlusOne[nOfSpecies][nOfSpecies]*/
};
typedef TProperties *TPropertiesPtr;

//*************************************************************************************************************
class T0DProperties : public TProperties {	// contains properties of the gas mixture
public:
	T0DProperties( TInputDataPtr input, TSpeciesPtr species )
		: TProperties( input, species )
			{ InitT0DProperties( input ); };
	~T0DProperties( void );

	Double			GetDensity( void ) { return fDensity; };
	Double			&GetDensityRef( void ) { return fDensity; };
	Double			GetMixHeatCapacity( void ) { return fMixHeatCapacity; };
	Double			GetMixConductivity( void ) { return fMixConductivity; };
	Double	 		GetMixViscosity( void ) { return fMixViscosity; };
	Double			GetMixMolarMass( void ) { return fMixMolarMass; };
	Double			&GetMixMolarMassRef( void ) { return fMixMolarMass; };
	Double			GetPressure( void ) { return fPressure; };
	Double			&GetPressureRef( void ) { return fPressure; };
	T0DRadiationPtr GetRadiation( void ) { return fRadiation; };

	void			CompMixtureProps( Double *heatCapacity, Double *conductivity, Double *mu, Double *Y, Double temp
			, Double &pressure, Double &density, EqOfState what, int nSpeciesIn, T0DSpeciesPtr species );

private:
	void 			InitT0DProperties( TInputDataPtr input );

	Double			fPressure;
	Double			fDensity;
	Double			fMixHeatCapacity;
	Double			fMixViscosity;		// [kg / (m*s)]
	Double			fMixConductivity;
	Double			fMixMolarMass;
	T0DRadiationPtr	fRadiation;
};
typedef T0DProperties *T0DPropertiesPtr;

#ifndef ZEROD
//*************************************************************************************************************
class T1DProperties : public TProperties {	// contains properties of the gas mixture
public:
	T1DProperties( TInputDataPtr input, TSpeciesPtr species ) 
			: TProperties( input, species )
			{ InitT1DProperties( input ); };
	~T1DProperties( void );

	VectorPtr 		GetViscosity( void ) { return fViscosity; };
	VectorPtr 		GetDensity( void ) { return fDensity; };
	VectorPtr 		GetConductivity( void ) { return fConductivity; };
	VectorPtr 		GetHeatCapacity( void ) { return fHeatCapacity; };
	VectorPtr 		GetMolarMass( void ) { return fMolarMass; };
	VectorPtr 		GetDCpdT( void ) { return fdCpdT; };
	T1DRadiationPtr GetRadiation( void ) { return fRadiation; };
	VectorPtr	 	GetDCpDT( void ) { return fdCpdT; };
	void 			CompMixtureProps( TFlameNodePtr flameNode, Double *Y, Double temp
										, Double pressure, TSpeciesPtr species );
	void			PrintProperties( int k );
	void			PrintProperties( TNewtonPtr bt );
private:
	void 			InitT1DProperties( TInputDataPtr input );

	VectorPtr		fViscosity;		// [kg / (m*s)]		//	vector of length maxGridPoints+2, where left boundary is fViscosity[-1] 
	VectorPtr		fDensity;		// [kg / m^3]		//	vector of length maxGridPoints+2, ...
	VectorPtr		fConductivity;		// [W / (m*K)]	//	vector of length maxGridPoints+2, ...
	VectorPtr		fHeatCapacity;		// [J / (kg*K)]	//	vector of length maxGridPoints+2, ...
	VectorPtr		fMolarMass;		// [kg / kmol]		//	vector of length maxGridPoints+2, ...
	T1DRadiationPtr	fRadiation;		
	VectorPtr		fdCpdT;			// 				//	vector of length maxGridPoints, ...
};
typedef T1DProperties *T1DPropertiesPtr;
#endif // ZEROD

//*************************************************************************************************************
class TFlame {
public:
	TFlame( FirstInputPtr firstInp ) { InitTFlame( firstInp ); };
	virtual ~TFlame( void );

/*	virtual TSpeciesPtr		GetSpecies( void ) = 0;
	virtual TPropertiesPtr	GetProperties( void ) = 0;
	virtual TReactionPtr	GetReaction( void ) = 0;*/
//	TSootPtr		GetSoot( void ) { return fSoot; };
	TInputDataPtr	GetInputData( void ) { return fInputData; };
	
	Double			dmdYAnalyt( int mNumber, int yNumber, Double *Y, Double mixDensity
								, Double *reactionRate, Double *theReactionRate );
	void			dmdTAnalyt( int mNumber, Double temp, Double &dmdT, Double *reactionRate );
	virtual Double	GetPressure( void ) { return fPressure->vec[fPressure->len]; };
	void			SetPressure( Double press ) { fPressure->vec[fPressure->len] = press; };
	Flag			ClipNegativeConcs( void ) { return fClipNegativeConcs; };
	int				CheckSolution( Double &temp, Double *Y, int YLen );

	IntVectorPtr	GetFuelIndexVec( void ) { return fFuelIndex; };
	int				GetFuelIndex( int i ) { return fFuelIndex->vec[i]; };
	int				GetFuelIndex( void ) { return fFuelIndex->vec[0]; };
	int				GetNFuels( void ) { return fFuelIndex->len; };
	int				GetFromSpeciesIndex( void ) { return fFromSpeciesIndex; };
	int				GetToSpeciesIndex( void ) { return fToSpeciesIndex; };
	Double			GetNu( char *globalReaction, const char *name );
	Double			GetNuProduct( char *globalReaction, const char *name );
	char			*GetDateCreated( void ) { return fDateCreated; };
	char			*GetAuthor( void ) { return fAuthor; };
	enum			FileType	{ kNone, kData, kText, kBin };
	char			*GetOutputPath( void ) { return fOutputPath; };
	char			*GetOutFileBuff( void ) { return fOutFile; };
	FILE 			*GetOutfile( const char *name, FileType type );
	char			*GetOutfileName( const char *name, FileType type );

protected:
	VectorPtr 		fPressure;
	TInputDataPtr	fInputData;
//	TReactionPtr	fReaction;
//	TSpeciesPtr		fSpecies;
//	TPropertiesPtr	fProperties;
//	TSootPtr		fSoot;
	ContinType		fContinType;
	Double			fContBound;
	Double			fContInc;
	Flag			fUseNumericalJac;
	Flag			fUseNumericalDM;

	// Sensitivity Analysis
	Flag			fSensAnal;
	Flag                    fSensObjAll;         // all species are SensObj
	Flag                    fSensAnalSpec;       // Sens Anal on Species
	Flag                    fSensAnalReac;       // Sens Anal on Reactions
	Double                  fSensAnalFac;        // Factor by which A is multiplied
	Flag                    fSensMax;            // sensitivity on max and location of max
	Flag                    fSensFinal;          // sensitivity on final values
	char			**fSensObj;          // Name of species/variables used as target
	int			fNSensObj;
	Flag			fReactionFluxes;

	Flag			fClipNegativeConcs;

	void			PrintFlameletVector( int len, Double *vec, char *name, FILE *fp, int inc = 1 );
	void			NextPressure( void ) { fPressure->len++; };
	Flag			IsLastPressure( void ) 
					{ 
						return ( fPressure->len+1 == fPressure->phys_len ) 
								? TRUE : FALSE; 
					};

private:
	void	InitTFlame( FirstInputPtr firstInp ); 	// allocate and fill classes

	char			*fOutputPath;
	char			*fOutFile;
	IntVectorPtr	fFuelIndex;
	int				fFromSpeciesIndex;
	int				fToSpeciesIndex;
	char			*fDateCreated;
	char			*fAuthor;
};
//*************************************************************************************************************
class TMassFraction {
public:
	TMassFraction( T1DFlamePtr flame ) { InitTMassFraction( flame ); };
	void	UpdateMixFracDependence( T1DFlamePtr flame, NodeInfoPtr nodeInfo, Flag justTemperature = FALSE );
	void	PrintMassFraction( void );
	int		GetFuelIndex( void ) { return indexFuel; };
	int		GetOxIndex( void ) { return indexOx; };
	void	SetMassFractionsBC( int mixFracIndex, int firstSpeciesIndex, Double *yLeft, Double *yRight );
	Double	GetZStoe( void ) { return fZStoe; };
	
private:
	void 	InitTMassFraction( T1DFlamePtr flame );

	int		indexFuel;
	int		indexOx;
	int		indexH2O;
	int		indexCO2;
	int		aFH;
	int		aFC;
	Double	nuOx;
	Double	nuFuel;
	Double	*fYF1;
	Double	*fYOx2;
	Double	*fYF2;
	Double	*fYOx1;
	Double	fNu;
	Double	fZStoe;
	Double	fYH2OSt;
	Double	fYCO2St;
};
typedef TMassFraction *TMassFractionPtr;

//*************************************************************************************************************
class T0DFlame : public TFlame {
public:
	T0DFlame( FirstInputPtr firstInp ) 
		: TFlame( firstInp ) 
	{ 
		InitT0DFlame();
	};	
	virtual ~T0DFlame( void );

	// access functions
	T0DSpeciesPtr		GetSpecies( void ) { return fSpecies; };
	T0DPropertiesPtr	GetProperties( void ) { return fProperties; };
	T0DReactionPtr		GetReaction( void ) { return fReaction; };
	T0DSootPtr			GetSoot( void ) { return fSoot; };

	virtual int		GetOffsetFirstSpecies( void ) = 0;
	virtual int		GetOffsetTemperature( void ) = 0;
	virtual int		GetVariablesWithoutSpecies( void ) = 0;

	void			CompLewisNumbers( const char *lewisFile );

	void			UpdateThermoProps( Double *Y, Double temp, Double &pressure
								, Double &density, EqOfState what, Double *sootMoments );
	void			ComputeProperties( Double temp, Double *Y, Double &pressure
									, Double &density, EqOfState what, Double *sootMoments );
	Double 			ComputeZBilger( Double *Y, Double *YFuelSide, Double *YOxSide );


protected:

	T0DReactionPtr		fReaction;
	T0DSpeciesPtr		fSpecies;
	T0DPropertiesPtr	fProperties;
	T0DSootPtr			fSoot;

	Double				GetElementMassFraction( Double *Y, const char * const atomName, Double atomMolarMass );

private:
	void 			InitT0DFlame( void );

	// solution
//	VectorPtr		fMassFractions;
//	Double			fTemp;
};

#ifndef ZEROD
//*************************************************************************************************************
class T1DFlame : public TFlame {
typedef void ( *JacFuncPtr )( TFlameNodePtr fFlameNode );
public:
	T1DFlame( FirstInputPtr firstInp ) 
		: TFlame( firstInp ) 
	{ 
		InitT1DFlame();
	};	
	virtual ~T1DFlame( void );

	void				UpdateThermoProps( void );//
	void				UpdateThermoProps( TFlameNodePtr flameNode, NodeInfoPtr nodeInfo );
	void				ComputeProperties( TFlameNodePtr flameNode, Double temp
									, Double *Y, Double pressure );
	void				ComputePropertiesMod( TFlameNodePtr flameNode, Double temp
									, Double *Y, Double pressure );
	TBVPSolverPtr		GetSolver( void ) { return fSolver; };//
	T1DReactionPtr		GetReaction( void ) { return fReaction; };
	T1DSpeciesPtr		GetSpecies( void ) { return fSpecies; };
	T1DPropertiesPtr	GetProperties( void ) { return fProperties; };
	T1DSootPtr			GetSoot( void ) { return fSoot; };
	VectorPtr			GetStrainRateVector( void ) { return fStrainRate; };
	VectorPtr			GetTemperature( void ) { return fSolTemp; };
	MatrixPtr			GetMassFracs( void ) { return fSolMassFracs; };
//	VectorPtr			GetSavedTemperature( void ) { return fSavedTemp; };
//	MatrixPtr			GetSavedMassFracs( void ) { return fSavedMassFracs; };
	Double				GetStrainRate( void ) { return fStrainRate->vec[fStrainRate->len]; };
	void				SetStrainRate( Double strainRate ) { fStrainRate->vec[fStrainRate->len] = strainRate; };
	Flag				UseDiffCorr( void ) { return !fNoDiffCorr; };//
	void				UnSetDiffCorr( void ) { fNoDiffCorr = TRUE; };//
	VectorPtr			GetDiffCorr( void ) { return fDiffusivityCorrection; };//
	Double				ThermoDiffusion( int speciesIndex, CoordType coordinate, NodeInfoPtr nodeInfo );
	void				FillJacThermoDiffusion( int nVariable, Double constCoeff, CoordType coordinate, NodeInfoPtr nodeInfo );
	void				ComputeDiffusivityCorrection( Double **Y, NodeInfoPtr nodeInfo );//
	virtual FILE		*GetOutputFile( char *head, char *tail, FileType type ) = 0;
	void				CompLewisNumbers( const char *lewisFile );
	void 				SetFlameNode( int k );
	void				PostConvergence( void *object );
	void				SetMixtureSpecificationLeft( int num ) { fMixtureSpecificationLeft = num; };
	void				SetMixtureSpecificationRight( int num ) { fMixtureSpecificationRight = num; };
	int					GetMixtureSpecificationLeft( void ) { return fMixtureSpecificationLeft; };
	int					GetMixtureSpecificationRight( void ) { return fMixtureSpecificationRight; };
	TFlameNodePtr		GetFlameNode( void ) { return fFlameNode; };
	TFlameNodePtr		GetFlameNodeSaved( void ) { return fFlameNodeSaved; };
	void				SetBVPInput( TBVPInputPtr input );
	void				FillJacDiffusion( int nVariable, int nEquation, Double constCoeff, Double *diffCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	Double				StandardDiffusion( int nVariable, Double *diffCoeff, NodeInfoPtr nodeInfo );
	void				FillJacMixFracDiffusion( int nVariable, int nEquation, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	Double				SecondDerivMixFracDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double				GetZStoich( void );
	Double				GetZStoich_mf( void );
	Double				GetGeometry( void ) { return fGeometry; };
	void				XToEta( TNewtonPtr bt, VectorPtr etaVec );
	virtual void		EtaToX( TNewtonPtr /*bt*/, VectorPtr /*xPhysVec*/ ) { cerr << "#error: no definition of virtual function EtaToX" << NEWL;};
	void				OriginToZstoich( VectorPtr xVec, VectorPtr mixtureFraction, Double zStoich );
	Double				ComputeEmissionIndex( int speciesIndex, Double *x );
	void				SaveSolution( void );
	void				RestoreSolution( void );
	Flag				AdjustGrid( PostIterFuncPtr PostIter );
	int					SensitivityAnalysis( Double coeffSpecies, Double coeffTemp, CoordType coordinate );
	void				ReactionFluxes( CoordType coordinate );
	Double				dmdYAnalyt( int mNumber, int yNumber, Double *Y, TFlameNodePtr flameNode );
	void				FilldMdYOnePointNum( TFlameNodePtr flameNode );
	void				FilldMdTOnePointNum( TFlameNodePtr flameNode );
	void				FilldMomdMomOnePointNum( TFlameNodePtr flameNode );
	void				FilldMdYOnePointAnal( TFlameNodePtr flameNode );
	void				FilldMdTOnePointAnal( TFlameNodePtr flameNode );
	void				FilldMdYOnePoint( TFlameNodePtr flameNode );
	void				FilldMdTOnePoint( TFlameNodePtr flameNode );
	void				FilldMomdMomOnePoint( TFlameNodePtr flameNode );
	void				(T1DFlame::*FilldMdYOnePointPtr)( TFlameNodePtr flameNode );
	void				(T1DFlame::*FilldMdTOnePointPtr)( TFlameNodePtr flameNode );
	void				(T1DFlame::*FilldMomdMomOnePointPtr)( TFlameNodePtr flameNode );
	void				WriteRoggFiles( TNewtonPtr bt );
	Double 				ComputeZBilger( Double *Y, Double *YFuelSide, Double *YOxSide );

	virtual void	PrintRHSTemp( TNewtonPtr bt ) = 0;
	virtual void	UpdateSolutionOnePoint( Double */*y*/, int /*gridPoint*/ ) { cerr << "#error: wrong instance of function UpdateSolution called" << NEWL; exit(2); };

	virtual int		GetOffsetVVelocity( void ) = 0;
	virtual int		GetOffsetUVelocity( void ) = 0;
	virtual int		GetOffsetTemperature( void ) = 0;
	virtual int		GetOffsetMixFrac( void ) = 0;
	virtual int		GetOffsetFirstSpecies( void ) = 0;
	virtual int		GetVariablesWithoutSpecies( void ) = 0;
	virtual ConstStringArray GetVariableNames( void ) = 0;
	
protected:
	TBVPSolverPtr	fSolver;
	TFlameNodePtr	fFlameNode;
	TFlameNodePtr	fFlameNodeSaved;

	VectorPtr		fSolTemp;
	MatrixPtr		fSolMassFracs;

	VectorPtr		fSavedTemp;
	MatrixPtr		fSavedMassFracs;
	VectorPtr		fSavedGrid;

	VectorPtr		fMmod;

	Flag			fNoDiffCorr;
	VectorPtr		fDiffusivityCorrection; // [m^2/s]	//  this vector contains sum(D_k dY_k/dy)
														//	it has [nGridPoints+2] elements
	VectorPtr		fStrainRate;
	Double			fGeometry;   // 0.0 for planar, 1.0 for axi-symmetric geometry
	Flag			fThermoDiffusion;
	ContinSide		fContinSide;
	Flag			fPrintRHSSpecies;
	Flag				fPrintRHSTemp;
	T1DReactionPtr		fReaction;
	T1DSpeciesPtr		fSpecies;
	T1DPropertiesPtr	fProperties;
	T1DSootPtr			fSoot;

	void			SetInitialBC( TGridPtr grid, TInputDataPtr inp );
	void			CheckBC( void );
	void			ReadStartProfiles( TInputDataPtr inp );
	void			CheckInitialGuess( void );
	void			UpdateDimensions( int len );
	void			UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
	void			UpdateSolution( Double *y, int gridPoint );
	void			SolutionToSolver( void );
	void 			CopyFlameNode( TFlameNodePtr flameNodeSource, TFlameNodePtr flameNodeDest );
	void			CompareFlameNode( TFlameNodePtr flameNodeSource, TFlameNodePtr flameNodeDest );
	Flag			RHSAction( NodeInfoPtr nodeInfo, RHSMode rhsMode );
	Double			GetElementMassFraction( Double *Y, const char * const atomName, Double atomMolarMass );
	
private:
   	void			InitT1DFlame( void );
	TFlameNodePtr 	NewTFlameNodeSaved( void );
	void 			DisposeTFlameNodeSaved( TFlameNodePtr flameNodeSaved );
	int				fMixtureSpecificationLeft; // uses enum MixtureSpecification
	int				fMixtureSpecificationRight;
	virtual void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp ) = 0;
	void			SaveGrid( void );
	void			RestoreGrid( void );
	int				GetNOfX( Double theX, int nGridPoints, Double *x );
	int				CheckFlameLocation( Double xMid, Double deltaXVar );
	Double			CheckGradientLeft( void );
	Double			CheckGradientRight( void );
	void			PrintSensitivityCoefficients( TensorPtr sensitivityCoefficients, IntVectorPtr objects );
	void			PrintReactionFluxesReac( Double **fluxes, char *name, char *header );
	void 			PrintReactionFluxesSpec( Double **fluxes, char *name, char *header );
	void			PrintReduceInfo( Double **fluxes, Double epsilon, char *header );
	void			PrintOneSensMax( Double **sensCoeff, int gridPoints, int nReactions
									, char *fileName );
	Double			GetMaxSens( Double **sensCoeff, int gridPoints, int reaction );
	void			PrintSensMax( TensorPtr sensitivityCoefficients, IntVectorPtr objectsPtr );
};
#endif // ZEROD

class TPremixed {
public:
	Flag			PostConvTPremixed( int isConverged );
	Double			ComputeFlameThickness( Double *temp, Double *x, int nGridPoints );
	Flag			AdjustComputationalDomain( void ) { return fAdjustComputationalDomain; };
protected:
	TPremixed( TInputDataPtr input ) { InitTPremixed( input ); };
	~TPremixed( void );
	Double			GetPhi( void ) { return fPhi->vec[fPhi->len]; };
	Flag			IsLastPhi( void ) 
			{ return ( fPhi->len+1 == fPhi->phys_len ) ? TRUE : FALSE; };
	void			SetEGR( TFlamePtr flame, Double phi, Double *Y, int YLen, Double *molarMass, char **names, Double EGR );
	void			SetPhi( Double phi );
	void			SetMassFracsOfPhi( TFlamePtr flame, Double phi, Double *Y, int YLen, Double *molarMass, char **names, Double EGR = -1 );
	void			SetPhiOfMassFracs( TFlamePtr flame, Double *phi, Double *Y, Double *molarMass, char **names );
	void			NextPhi( void ) { ++fPhi->len; };
	void			InitPhi( void ) { fPhi->len = 0; };
	void			SetPhiVec( VectorPtr phi ) { fPhi = phi; };
	VectorPtr		GetPhiVector( void ) { return fPhi; };
	int				GetNOfPhi( void ) { return fPhi->phys_len; };
private:
	void			InitTPremixed( TInputDataPtr input );
	VectorPtr		fPhi;
	Flag			fAdjustComputationalDomain;
};

/*	to add a new flametype one has to add the class with all initializers for constants
 *	and pure virtual functions, a constructor, a destructor and an initializer
 *	the flametype has to be added to the enumeration list 'FlameType' and an 
 *	identifier for the input specification has to be introduced in file 'ReadAddData.flex' 
 *	in the startcondition 'scReadFlameType' 
 *	flame specific initialization of the class FirstInput has to be done in function 
 *	'FirstInput::InitFirstInput' and for the class TInputData in 
 *	'TInputData::InitInputData'
 *	the print functions of FirstInput and TInputData should be updated
 */

// prototypes

void SkipWhites( char * label );
void SkipIntroWhites( char * label );
void UpperString( char *string );
//void WriteRoggsContinData( TNewtonPtr bt, T1DFlamePtr flame );
void LinearInterpolate( Double *x_old, Double *y_old, int n_old,
						Double *x, Double *y, int n_new );

Double InterpolOne( Double x, Double *x_old, Double *y_old, int n_old );

Double GetXStagnation( int variable, int nGridPoints, Double *x, Double **y );

void FillJacFirstDerivCentral( int nVariable, int nEquation, Double coeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );

void FillJacUpWindConvection( int nVariable, int nEquation, NodeInfoPtr nodeInfo, Flag sign = kPositive );

void FillJacFirstDerivUp( int nVariable, int nEquation, NodeInfoPtr nodeInfo, Flag sign = kPositive );

void FillJacFirstDerivDown( int nVariable, int nEquation, NodeInfoPtr nodeInfo, Flag sign = kPositive );

void FillJacNonlinearConvectCentral( int nVariable1, int nVariable2, NodeInfoPtr nodeInfo, Double coeff = 1.0 );

Double NonlinearConvectCentral( Double y1, Double y2Prev, Double y2, Double y2Next, Double hm, Double h );

void FillJacNonlinearConvectUpwind( int nVariable1, int nVariable2, NodeInfoPtr nodeInfo, Double coeff, Flag velocityPositive = TRUE );

Double NonlinearConvectUpwind( Double y1, Double y2Prev, Double y2, Double y2Next, Double hm, Double h, Flag velocityPositive = TRUE );

void FillJacFTimesSecondDerivF( int nVariable, int nEquation, NodeInfoPtr nodeInfo, Flag sign = kPositive );

void FillJacWithDiffusion( int nVariable, int nEquation, Double constCoeff, Double *coeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );

void FillJacSecondDerivCentral( int nVariable, int nEquation, Double coeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );

void FillJacSpeciesDiffusion( int nVariable, int speciesIndex, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );

Double FirstDerivLeft( Double yPrev, Double y, Double yNext, Double hm, Double h );

Double FirstDerivSq( Double yPrev, Double y, Double yNext, Double hm2, Double h2, Double hnenn );

Double FirstDeriv( Double yPrev, Double y, Double ynext, Double hm, Double h );

Double FirstDerivUpwind( Double y2, Double y1, Double h );

Double FTimesSecondDerivF( int nVariable, NodeInfoPtr nodeInfo );

Double SecondDeriv( Double yPrev, Double y, Double ynext, Double hm, Double h );

Double SecondDerivDiffusion( int nVariable, Double *coeff, NodeInfoPtr nodeInfo );

Double SecondDerivWeightedDiffusion( int nVariable, Double coeff, NodeInfoPtr nodeInfo );

Double SecondDerivMassDiffusion( int nVariable, NodeInfoPtr nodeInfo );

Double ImplicitSpeciesDiffusion( int nVariable, int speciesIndex, NodeInfoPtr nodeInfo );

Double myPow( Double base, Double expo );

Double modPow( Double base, int expo );

Double myLog10( Double arg );

Double Signum( Double val );

void CopyStringArray( char **dest, char **source, int n );

Double Integral( int nPoints, Double *x, Double *y );

Double IntegralFull( int nPoints, Double *x, Double *y );

Double	SolveQuadratic( Double a, Double b, Double c );

void SaveArray( Double **matrix, int rows, int cols, int partition, Double *x, 
		ConstStringArray titles, char *name );
void SaveArrayBin( Double **matrix, int rows, int cols, int partition, Double *x, ConstStringArray titles, char *name, char *path );
		
extern "C" int myCompare( const void *elem1, const void *elem2 );

int LocationOfMax( int len, Double *vec );
int LocationOfMin( int len, Double *vec );
int LocationOfMax( int len, Double *vec, int offset );
int LocationOfMin( int len, Double *vec, int offset );
int LocationOfAbsMax( int len, Double *vec );
int LocationOfMaxSlope( Double *vec, Double *x, int len );

extern "C" int SteadyStatesFunc( const VectorPtr x, VectorPtr fVec, void *object );
void ComputeSteadyStates( Double *k, Double *c, Double *tBConc );
void InverseMatrix( int n, Double **a, Double **inv, int *index, Double *col );

Double LambdaOverCp ( Double temp );

MatrixPtr FortranToCMat( Double *x, int rows, int physRows, int cols, int physCols );
void DisposeFToCMat( MatrixPtr mat );
char *FortranToCChar( char *a, int nChars );
char **FortranToCCharArray( char *a, int nChars, int len );
void CToFortranCharArray( char *aF, char **aC, int nChars, int len );
void CToFortranChar( char *aF, char *aC, int nChars );
void DisposeFToCCharArray( char **newA, int len );
Double TrapIntegrate( int nG, Double *f, Double *x );
Double *TrapIntegrate( int nG, Double *f, Double *x, Double *res );
Double DotProd( int n, const Double *r, int inc_r, const Double *c, int inc_c );
char *MyDataFile( const char* name );


#ifdef HP
int matherr( struct exception *x );
#endif

#endif  /* __FlameMan__ */
