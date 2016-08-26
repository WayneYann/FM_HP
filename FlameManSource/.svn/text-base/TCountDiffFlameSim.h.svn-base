void CountDiffSimJacRest( void *object, NodeInfoPtr nodeInfo );
void CountDiffSimRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
void CountDiffSimOutput( void *object, FILE *fp, char* tail );
int CountDiffSimPostIter( void *object );
void SetCountDiffSimNodeInfo( int k, void *object );
void CountDiffSimPostConv( void *object );
void CountDiffSimUpdateLeftBoundary( void  *object );
void CountDiffSimUpdateRightBoundary( void *object );
ConstStringArray GetCountDiffSimVarNames( void *object );
// the following function is used as utility function for dfdyUpwind
// the void * is an object of the class T1DFlame
Double Numerical_rhoDk( int equation, NodeInfoPtr nodeInfo, void *object, Flag CalcNewProperties = FALSE );
Double EntFluxFuncDiffSim( int equation, NodeInfoPtr nodeInfo, void *object, Flag theFlag = FALSE );


//*************************************************************************************************************
class TCountDiffFlameSim : public T1DFlame {
friend void CountDiffSimJacRest( void *object, NodeInfoPtr nodeInfo );
friend void CountDiffSimRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
friend void CountDiffSimOutput( void *object, FILE *fp, char* tail );
friend int CountDiffSimPostIter( void *object );
friend void CountDiffSimPostConv( void *object );
friend void CountDiffSimUpdateLeftBoundary( void  *object );
friend void CountDiffSimUpdateRightBoundary( void *object );
friend Double EntFluxFuncDiffSim( int equation, NodeInfoPtr nodeInfo, void *object, Flag theFlag  );

public:
	TCountDiffFlameSim( FirstInputPtr firstInp ) : fVVelocity(0), fUVelocity(1), fMixFrac(2), 
//								fTemperature(3), fFirstSpecies(4),
								fFirstSpecies(3), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
								fSootMoments(fTemperature+1), fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
								T1DFlame( firstInp ) { InitCountDiffFlameSim();};	
	~TCountDiffFlameSim( void );

	int		GetOffsetVVelocity( void );
	int		GetOffsetUVelocity( void );
	int		GetOffsetTemperature( void );
	int		GetOffsetMixFrac( void );
	int		GetOffsetFirstSpecies( void );
	int		GetVariablesWithoutSpecies( void );
	ConstStringArray	GetVariableNames( void );
	
	VectorPtr	GetV( void ) { return fSolV; };
	VectorPtr	GetU( void ) { return fSolU; };
	VectorPtr	GetZ( void ) { return fSolMixFrac; };

	void	SaveSolution( void );
	void	RestoreSolution( void );
	void	PrintRHSSpecies( TNewtonPtr bt );

protected:
	void	PrintRHSTemp( TNewtonPtr bt );

	const int			fVVelocity;
	const int			fUVelocity;
	const int			fMixFrac;
	const int			fFirstSpecies;
	const int			fTemperature;
	const int			fVariablesWithoutSpecies;
	char				**fVariableNames;
	Flag				fPrintMolarFractions;
	TMassFractionPtr	fMassFraction;
	const int			fSootMoments;
	
//	vectors of solution
	VectorPtr			fSolV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSolU;			// len is set to nGridPoints
	VectorPtr			fSolMixFrac;	// first element is fSol->vec[-1]

	VectorPtr			fSavedV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSavedU;			// len is set to nGridPoints
	VectorPtr			fSavedMixFrac;		// first element is fSol->vec[-1]

	void	UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
	void	UpdateSolutionOnePoint( Double *y, int gridPoint );

private:
	Double	GetChiStoichAnal( void );
	Double	GetChiOverChiStoich( Double z );
	void	ScalarDissipation( Double *Z, TNewtonPtr bt, Double *scalarDiss, Double *stoechScalarDiss );
	void	ScalarDissipationBilger( TNewtonPtr bt, VectorPtr scalarDissVec, Double *stoechScalarDiss, Double *Z );
  	void	InitCountDiffFlameSim( void );
//	void	SetInitialBC( TGridPtr grid, TInputDataPtr inp );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	void	FilldTdxdYdxOverTT( Double coeff, NodeInfoPtr nodeInfo );
	void	FilldTdxdYdxOverTTUpwind( Double coeff, NodeInfoPtr nodeInfo );
	void	FillImplicitSpeciesDiffusion( int speciesIndexk, NodeInfoPtr nodeInfo );
	void	FillJacSpeciesDiffusion( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	void	FillJacFullDiffusion( int lInd, int kInd, NodeInfoPtr nodeInfo );
	Double	NewDiffCorr( int nVariable, NodeInfoPtr nodeInfo );
	Double	NewDiffCorrX( int nVariable, NodeInfoPtr nodeInfo );
	Double	SecondDerivBinSpecDiff( int nVariable, NodeInfoPtr nodeInfo );
	Double	FullSpeciesDiffusion( int lInd, int kInd, NodeInfoPtr nodeInfo );
	Double	ImplicitSpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	SpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	SecondDerivSpeciesDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	SecondDerivXDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	SecondDerivXMixFracDiff( int nVariable, NodeInfoPtr nodeInfo );
	Double	SecondDerivConstMixFracDiff( int nVariable, NodeInfoPtr nodeInfo );
	Double	HeatFluxBinSpecDiff( NodeInfoPtr nodeInfo );
	void	PrintRHSSpecies( int start, int end, NodeInfoPtr nodeInfo, Double physX, FILE *fp );
	void	PrintRHSTemp( TNewtonPtr bt, NodeInfoPtr nodeInfo, Double physX, FILE *fp );
	Double	ComputeAbsStrainRate( int k1, Double f1, Double f2, Double fp1, Double fp2, Double h );
	void	FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	Double	Numerical_dRhoDkdeta( int equation, NodeInfoPtr nodeInfo, Flag CalcNewProperties );
	void	EtaToX( TNewtonPtr bt, VectorPtr xPhysVec );
	FILE	*GetOutputFile( char *head, char *tail, FileType type );
	void	UpdateDimensions( int len );
	void	SolutionToSolver( void );
	void	CompLewis( int which, Double *Le );
};

typedef TCountDiffFlameSim *TCountDiffFlameSimPtr;
