void CountPremSimCheckCompDomain( void *object, NodeInfoPtr nodeInfo );
void CountPremSimJacRest( void *object, NodeInfoPtr nodeInfo );
void CountPremSimRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
void CountPremSimOutput( void *object, FILE *fp, char* tail );
int CountPremSimPostIter( void *object );
void SetCountPremSimNodeInfo( int k, void *object );
void CountPremSimPostConv( void *object );
void CountPremSimUpdateLeftBoundary( void  *object );
void CountPremSimUpdateRightBoundary( void *object );
ConstStringArray GetCountPremSimVarNames( void *object );


//*************************************************************************************************************
class TCountPremFlameSim : public T1DFlame, public TPremixed {
friend void CountPremSimCheckCompDomain( void *object, NodeInfoPtr nodeInfo );
friend void CountPremSimJacRest( void *object, NodeInfoPtr nodeInfo );
friend void CountPremSimRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
friend void CountPremSimOutput( void *object, FILE *fp, char* tail );
friend int CountPremSimPostIter( void *object );
friend void CountPremSimPostConv( void *object );
friend void CountPremSimUpdateLeftBoundary( void  *object );
friend void CountPremSimUpdateRightBoundary( void *object );

public:
	TCountPremFlameSim( FirstInputPtr firstInp ) : fVVelocity(0), fUVelocity(1), 
								fFirstSpecies(2), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
								fLnStrainrate(fTemperature+1),
								fSootMoments(fLnStrainrate+1), fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
								T1DFlame( firstInp ), TPremixed( fInputData ) { InitCountPremFlameSim();};	
	~TCountPremFlameSim( void );

	int		GetOffsetVVelocity( void );
	int		GetOffsetUVelocity( void );
	int		GetOffsetTemperature( void );
	int		GetOffsetFirstSpecies( void );
	int		GetOffsetLnStrainrate( void );
	int		GetVariablesWithoutSpecies( void );
	ConstStringArray	GetVariableNames( void );
// not used by TCountPremFlame but required by TFlame as pure virtual functions
	int		GetOffsetMixFrac( void );
	
	VectorPtr	GetV( void ) { return fSolV; };
	VectorPtr	GetU( void ) { return fSolU; };

	void	SaveSolution( void );
	void	RestoreSolution( void );
	void	PrintRHSSpecies( TNewtonPtr bt );

// continuation stuff
	void		SetTempContStart( Double val ) { fTempContStart = val; };
	void		SetStrainrateContStart( Double val ) { fLnStrainrateContStart = val; };
	void		SetDeltaArcLength( Double val ) { fDeltaArcLength = val; };
	void		SetTMaxLoc( int loc ) { fTmaxLoc = loc; };
	void		SetdTds( Double dTds ) { fdTds = dTds; };
	void		SetdlnStrainrateds( Double dlnStrainrateds ) { fdlnStrainrateds = dlnStrainrateds; };
	void		IncNFlameletCount( void ) { ++fNFlameletsCount; };
	void		ReInitArcLengthCont( void ) { fDeltaArcLength = 0.0; fdTds = 0.0; fdlnStrainrateds = 0.0; fNFlameletsCount = 0; };
	Flag		GetArcLengthCont( void ) { return fArcLengthContin; };
	Flag		GetArcUp( void ) { return fArcUp; };
	int			GetTMaxLoc( void );
	int			GetMaxFlamelets( void ) { return fMaxFlamelets; };
	int			GetNFlameletsCount( void ) { return fNFlameletsCount; };
	Double		GetdTds( void ) { return fdTds; };
	Double		GetdlnStrainrateds( void ) { return fdlnStrainrateds; };
	Double		GetStrainrateContStart( void ) { return fLnStrainrateContStart; };
	Double		GetTempContStart( void ) { return fTempContStart; };
	Double		GetDeltaArcLength( void ) { return fDeltaArcLength; };
	Double		GetDeltaStrainrateref( void ) { return fDeltaStrainrateRef; };
	Double		GetDeltaTref( void ) { return fDeltaTRef; };
	int			GetLnStrainrate( void ) { return fLnStrainrate; };
//	Double		GetStrainRate( void ) { return (fArcLengthContin) ? exp( fSolLnStrainrate->vec[1] ) : T1DFlame::GetStrainRate(); };
	Double		GetStrainRate( void ) { return exp( fSolLnStrainrate->vec[1] ); };

protected:
	void	PrintRHSTemp( TNewtonPtr bt );

	const int			fVVelocity;
	const int			fUVelocity;
	const int			fFirstSpecies;
	const int			fTemperature;
	const int			fVariablesWithoutSpecies;
	const int			fLnStrainrate;
	char				**fVariableNames;
	Flag				fPrintMolarFractions;
	TMassFractionPtr	fMassFraction;
	const int			fSootMoments;
	int					fPremConfiguration; // uses enum CountPremType (kBackToBack, kFreshToBurned, kFull)
	Double				fEnthalpy;
	
//	vectors of solution
	VectorPtr			fSolV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSolU;			// len is set to nGridPoints
	VectorPtr			fSolLnStrainrate;	// 

	VectorPtr			fSavedV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSavedU;			// len is set to nGridPoints
	VectorPtr			fSavedLnStrainrate;		// 

	void	UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
	void	UpdateSolutionOnePoint( Double *y, int gridPoint );

private:
// continuation stuff
	int			fTmaxLoc;			//	location of maximum temperature
	int			fMaxFlamelets;
	int			fNFlameletsCount;
	Flag		fArcLengthContin;
	Flag		fStrainRateContin;
	Flag		fArcUp;
	Double		fTempContStart;
	Double		fLnStrainrateContStart;
	Double		fSavedTempContStart;
	Double		fSavedLnStrainrateContStart;
	Double		fDeltaArcLength;
	Double		fdTds;				//	dT_max/ds, where s is the arclength
	Double		fDeltaStrainrateRef;
	Double		fDeltaTRef;
	Double		fdlnStrainrateds;			//	da_inv/ds, where s is the arclength and a_inv the inverse of the strainrate

  	void	InitCountPremFlameSim( void );
	void	SetInitialBC( TGridPtr grid, TInputDataPtr inp );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	Double	NewDiffCorr( int nVariable, NodeInfoPtr nodeInfo );
	Double	NewDiffCorrX( int nVariable, NodeInfoPtr nodeInfo );
	Double	SecondDerivBinSpecDiff( int nVariable, NodeInfoPtr nodeInfo );
	Double	FullSpeciesDiffusion( int lInd, int kInd, NodeInfoPtr nodeInfo );
	Double	ImplicitSpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	SpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	SecondDerivSpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	SecondDerivXDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	HeatFluxBinSpecDiff( NodeInfoPtr nodeInfo );
	Double	GetBurningVelocity( void );
	int		CheckComputationalDomain( void );
	void	PrintRHSSpecies( int start, int end, NodeInfoPtr nodeInfo, Double physX, FILE *fp );
	void	PrintRHSTemp( TNewtonPtr bt, NodeInfoPtr nodeInfo, Double physX, FILE *fp );
	void	EtaToX( TNewtonPtr bt, VectorPtr xPhysVec );
	FILE	*GetOutputFile( char *head, char *tail, FileType type );
	void	UpdateDimensions( int len );
	void	SolutionToSolver( void );
	void	CompLewis( int which, Double *Le );
};

typedef TCountPremFlameSim *TCountPremFlameSimPtr;
