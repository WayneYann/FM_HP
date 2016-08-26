void CountDiffPhysJacFirst( void *object, NodeInfoPtr nodeInfo );
void CountDiffPhysJacRest( void *object, NodeInfoPtr nodeInfo );
void CountDiffPhysJacLast( void *object, NodeInfoPtr nodeInfo );
void CountDiffPhysRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
void CountDiffPhysOutput( void *object, FILE *fp, char* tail );
int CountDiffPhysPostIter( void *object );
void SetCountDiffPhysNodeInfo( int k, void *object );
void CountDiffPhysPostConv( void *object );
void CountDiffPhysUpdateRightBoundary( void *object );
void CountDiffPhysUpdateLeftBoundary( void  *object );
ConstStringArray GetCountDiffPhysVarNames( void *object );

//*************************************************************************************************************
class TCountDiffFlamePhys : public T1DFlame {
friend void CountDiffPhysJacFirst( void *object, NodeInfoPtr nodeInfo );
friend void CountDiffPhysJacRest( void *object, NodeInfoPtr nodeInfo );
friend void CountDiffPhysJacLast( void *object, NodeInfoPtr nodeInfo );
friend void CountDiffPhysRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
friend void CountDiffPhysOutput( void *object, FILE *fp, char* tail );
friend int CountDiffPhysPostIter( void *object );
friend void CountDiffPhysUpdateRightBoundary( void *object );
friend void CountDiffPhysUpdateLeftBoundary( void  *object );
friend void CountDiffPhysPostConv( void *object );

public:
	TCountDiffFlamePhys( FirstInputPtr firstInp ) : fVVelocity(0), fUVelocity(1), fMixFrac(2), 
//								fTemperature(3), fFirstSpecies(4),
								fFirstSpecies(3), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
								fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
								T1DFlame( firstInp ) { InitTCountDiffFlamePhys();};	
	~TCountDiffFlamePhys( void );

	int		GetOffsetVVelocity( void );
	int		GetOffsetUVelocity( void );
	int		GetOffsetTemperature( void );
	int		GetOffsetMixFrac( void );
	int		GetOffsetFirstSpecies( void );
	int		GetVariablesWithoutSpecies( void );
	ConstStringArray	GetVariableNames( void );
	Double				GetdYNextdY( int speciesIndex, NodeInfoPtr nodeInfo );
	Double				GetdYPrevdY( int speciesIndex, NodeInfoPtr nodeInfo );
	
	VectorPtr	GetV( void ) { return fSolV; };
	VectorPtr	GetU( void ) { return fSolU; };
	VectorPtr	GetZ( void ) { return fSolMixFrac; };

protected:
	void	PrintRHSTemp( TNewtonPtr bt ) { bt = bt; fprintf( stderr, "nothing happens\n" );};

private:
	const int			fVVelocity;
	const int			fUVelocity;
	const int			fMixFrac;
	const int			fFirstSpecies;
	const int			fTemperature;
	const int			fVariablesWithoutSpecies;
	char				**fVariableNames;
	Flag				fPrintMolarFractions;
	TMassFractionPtr	fMassFraction;

//	vectors of solution
	VectorPtr			fSolV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSolU;			// len is set to nGridPoints
	VectorPtr			fSolMixFrac;	// first element is fSol->vec[-1]

	void	ScalarDissipation( TNewtonPtr bt, VectorPtr scalarDissVec, Double *stoechScalarDiss );
//	void	XToEta( TNewtonPtr bt, VectorPtr etaVec );
  	void	InitTCountDiffFlamePhys( void );
//	void	SetInitialBC( TGridPtr grid, TInputDataPtr inp );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	void	FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	Double	DiffCorr( int nVariable, NodeInfoPtr nodeInfo );
	Double	MixFracDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	SpeciesDiffusion( int nVariable, int speciesIndex, NodeInfoPtr nodeInfo );
	void	FillJacMixFracDiffusion( int nVariable, int nEquation, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	void	FillJacSpeciesDiffusion( int nVariable, int speciesIndex, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	FILE	*GetOutputFile( char *head, char *tail, FileType type );
	void	UpdateDimensions( int len );
	void	UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
};

typedef TCountDiffFlamePhys *TCountDiffFlamePhysPtr;
