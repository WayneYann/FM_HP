#include "TCountDiffFlameSim.h"

//*************************************************************************************************************
class TSolution {
public:
	VectorPtr	CopyIt( VectorPtr source );
	MatrixPtr	CopyIt( MatrixPtr source );
	char 		**CopyIt( char **source, int len );

	Double		fLeft;
	Double		fRight;
	Double		fPressure;
	VectorPtr	fSolLeft;
	VectorPtr	fSolRight;
	MatrixPtr	fSol;
	VectorPtr	fCoord;
	int			fNOfSolNames;
	char		**fSolNames;
};

typedef TSolution *TSolutionPtr;

void CountDiffContJacRest( void *object, NodeInfoPtr nodeInfo );
void CountDiffContRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
void CountDiffContOutput( void *object, FILE *fp, char* tail );
int CountDiffContPostIter( void *object );
void SetCountDiffContNodeInfo( int k, void *object );
void CountDiffContPostConv( void *object );
void CountDiffContUpdateLeftBoundary( void  *object );
void CountDiffContUpdateRightBoundary( void *object );
Double EntFluxFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag );
Double EntFluxFuncDiffCont( int equation, NodeInfoPtr nodeInfo, void *object, Flag theFlag = FALSE );
ConstStringArray GetCountDiffContVarNames( void *object );

//*************************************************************************************************************
class TCountDiffFlameCont : public T1DFlame {
friend void CountDiffContJacRest( void *object, NodeInfoPtr nodeInfo );
friend void CountDiffContRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
friend void CountDiffContOutput( void *object, FILE *fp, char* tail );
friend int CountDiffContPostIter( void *object );
friend void CountDiffContPostConv( void *object );
friend void CountDiffContUpdateLeftBoundary( void  *object );
friend void CountDiffContUpdateRightBoundary( void *object );
friend Double EntFluxFunc( int j, NodeInfoPtr nodeInfo, void *object, Flag theFlag );
friend Double EntFluxFuncDiffCont( int equation, NodeInfoPtr nodeInfo, void *object, Flag theFlag );

public:
/*	TCountDiffFlameCont( FirstInputPtr firstInp ) : fVVelocity(0), fUVelocity(1), fMixFrac(2), 
								fTemperature(3), fInvStrainRate(4), fFirstSpecies(5),
//								fFirstSpecies(3), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
								fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
								T1DFlame( firstInp ) { InitCountDiffFlameCont();};	*/
	TCountDiffFlameCont( FirstInputPtr firstInp, TSolutionPtr sol, Double strainRate ) : fVVelocity(0), fUVelocity(1), fMixFrac(2), 
								fTemperature(3), fInvStrainRate(4), fFirstSpecies(5),
//								fFirstSpecies(3), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
								fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
								T1DFlame( firstInp ) { InitCountDiffFlameCont( sol, strainRate ); };	
	~TCountDiffFlameCont( void );

	int		GetOffsetVVelocity( void );
	int		GetOffsetUVelocity( void );
	int		GetOffsetTemperature( void );
	int		GetOffsetMixFrac( void );
	int		GetOffsetInvStrainRate( void );
	int		GetOffsetFirstSpecies( void );
	int		GetVariablesWithoutSpecies( void );
	ConstStringArray	GetVariableNames( void );
	
	VectorPtr	GetV( void ) { return fSolV; };
	VectorPtr	GetU( void ) { return fSolU; };
	VectorPtr	GetZ( void ) { return fSolMixFrac; };
	VectorPtr	GetOneOverA( void ) { return fSolInvStrainRate; };

	void		SetdTds( Double dTds ) { fdTds = dTds; };
	void		SetdAInvds( Double dAInvds ) { fdAInvds = dAInvds; };
	void		SetDeltaS( Double ds ) { fDeltaS = ds; };
	Double		GetDeltaS( void ) { return fDeltaS; };
	Double		GetStrainRate( void );
	Double		GetMinStrainRate( void ) { return fMinStrainRate; };
	Flag		ComputeUnphysicalChain( void ) { return fCompUnPhysChain; };

	void	SaveSolution( void );
	void	RestoreSolution( void );
	void	PrintRHSSpecies( TNewtonPtr bt );
	
protected:
	void	PrintRHSTemp( TNewtonPtr bt );

private:
	const int			fVVelocity;
	const int			fUVelocity;
	const int			fMixFrac;
	const int			fFirstSpecies;
	const int			fTemperature;
	const int			fInvStrainRate;
	const int			fVariablesWithoutSpecies;
	char				**fVariableNames;
	TMassFractionPtr	fMassFraction;
	
//	vectors of solution
	VectorPtr			fSolV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSolU;			// len is set to nGridPoints
	VectorPtr			fSolMixFrac;	// first element is fSol->vec[-1]
	VectorPtr			fSolInvStrainRate;	// 

	VectorPtr			fSavedV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSavedU;			// len is set to nGridPoints
	VectorPtr			fSavedMixFrac;		// first element is fSol->vec[-1]
	VectorPtr			fSavedInvStrainRate;		// 

	int		fTmaxLoc;			//	location of maximum temperature
	Double	fArcLengthStart;	//	arclength for the first flamelet
	Double	fAInvStart;			//	inverse of strainrate of the first flamelet
	Double	fTmaxStart;			//	maximum temperature the first flamelet
	Double	fdTds;				//	dT_max/ds, where s is the arclength
	Double	fdAInvds;			//	da_inv/ds, where s is the arclength and a_inv the inverse of the strainrate
	Double	fSaveddTds;			//	saved value of dT_max/ds
	Double	fSaveddAInvds;		//	saved value of da_inv/ds
	Double	fArcLength;			//	current arclength
	Double	fDeltaS;			//	free parameter of the problem: difference of arclength's between two flamelets
	Double	fAInvOld;			//	inverse of strainrate of previous flamelet
	Double	fSConverged;
	Double	fSNotConverged;
	Double	fMinStrainRate;
	Flag	fCompUnPhysChain;
	int		fNUnPhysChain;
	Flag	fUnPhysFound;		//	signals if an unphysical solution has been found so far
	FILE	*fpExtCurve;
	void	ScalarDissipation( TNewtonPtr bt, VectorPtr scalarDissVec, Double *stoechScalarDiss );
//  	void	InitCountDiffFlameCont( void );
  	void	InitCountDiffFlameCont( TSolutionPtr sol, Double strainRate );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	void	SetInitialValues( TInputDataPtr inp, TSolutionPtr sol, Double strainRate );
	void	FillJacSpeciesDiffusion( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	Double	NewDiffCorr( int nVariable, NodeInfoPtr nodeInfo );
	Double	SpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	SecondDerivXDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	NewDiffCorrX( int nVariable, NodeInfoPtr nodeInfo );
	void	PrintRHSSpecies( TNewtonPtr bt, NodeInfoPtr nodeInfo, Double physX, FILE *fp );
	void	PrintRHSTemp( TNewtonPtr bt, NodeInfoPtr nodeInfo, FILE *fp );
	void	FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	void	EtaToX( TNewtonPtr bt, VectorPtr xPhysVec );
	FILE	*GetOutputFile( char *head, char *tail, FileType type );
	void	UpdateDimensions( int len );
	void	UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
	void	UpdateSolutionOnePoint( Double *y, int gridPoint );
	void	SolutionToSolver( void );
};

typedef TCountDiffFlameCont *TCountDiffFlameContPtr;

void CountDiffContPostConvPre( void *object );
class TCountDiffFlameContPre : public TCountDiffFlameSim {
//	storage for fSolPack is not deleted in destructor
friend void CountDiffContPostConvPre( void *object );
public:
	TCountDiffFlameContPre( FirstInputPtr firstInp ) 
			: TCountDiffFlameSim( firstInp ) { InitCountDiffFlameContPre( firstInp ); };	
	~TCountDiffFlameContPre( void );
	Double			GetdTds( void ) { return fdTds; };
	Double			GetdAInvds( void ) { return fdAInvds; };
	TSolutionPtr	GetSolution( void ) { return fSolPack; };

private:
	void	InitCountDiffFlameContPre( FirstInputPtr firstInp );
//	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	void	PostAll( void );

	Double fTmax[2];
	Double fAInv[2];

//	Double	fStrainRateIn;		//	strainrate read from input file
	Double	fdTds;				//	dT_max/ds, where s is the arclength
	Double	fdAInvds;			//	da_inv/ds, where s is the arclength and a_inv the inverse of the strainrate

	TSolutionPtr	fSolPack;
};

typedef TCountDiffFlameContPre *TCountDiffFlameContPrePtr;
