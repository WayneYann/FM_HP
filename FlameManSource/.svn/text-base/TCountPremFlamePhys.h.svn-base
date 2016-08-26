void CountPremPhysJacFirst( void *object, NodeInfoPtr nodeInfo );
void CountPremPhysJacRest( void *object, NodeInfoPtr nodeInfo );
void CountPremPhysJacLast( void *object, NodeInfoPtr nodeInfo );
void CountPremPhysRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
void CountPremPhysOutput( void *object, FILE *fp, char* tail );
int CountPremPhysPostIter( void *object );
void CountPremPhysUpdateLeftBoundary( void  *object );
void CountPremPhysUpdateRightBoundary( void  *object );
void SetCountPremPhysNodeInfo( int k, void *object );
void CountPremPhysPostConv( void *object );
ConstStringArray GetCountPremPhysVarNames( void *object );

//*************************************************************************************************************
class TCountPremFlamePhys : public T1DFlame, public TPremixed {
friend void CountPremPhysJacFirst( void *object, NodeInfoPtr nodeInfo );
friend void CountPremPhysJacRest( void *object, NodeInfoPtr nodeInfo );
friend void CountPremPhysJacLast( void *object, NodeInfoPtr nodeInfo );
friend void CountPremPhysRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
friend void CountPremPhysOutput( void *object, FILE *fp, char* tail );
friend int CountPremPhysPostIter( void *object );
friend void CountPremPhysUpdateLeftBoundary( void  *object );
friend void CountPremPhysUpdateRightBoundary( void  *object );

public:
	TCountPremFlamePhys( FirstInputPtr firstInp ) : fVVelocity(0), fUVelocity(1),
								fFirstSpecies(2), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
								fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
								T1DFlame( firstInp ), TPremixed( fInputData ) { InitTCountPremFlamePhys();};	
	~TCountPremFlamePhys( void );

	int		GetOffsetVVelocity( void );
	int		GetOffsetUVelocity( void );
	int		GetOffsetTemperature( void );
	int		GetOffsetFirstSpecies( void );
	int		GetVariablesWithoutSpecies( void );
	ConstStringArray	GetVariableNames( void );
	Double				GetdYNextdY( int speciesIndex, NodeInfoPtr nodeInfo );
	Double				GetdYPrevdY( int speciesIndex, NodeInfoPtr nodeInfo );
// not used by TCountPremFlame but required by TFlame as pure virtual functions
	int		GetOffsetMixFrac( void );
	
	VectorPtr	GetV( void ) { return fSolV; };
	VectorPtr	GetU( void ) { return fSolU; };

protected:
	void	PrintRHSTemp( TNewtonPtr bt ) { bt = bt; fprintf( stderr, "nothing happens\n" );};

private:
	const int			fVVelocity;
	const int			fUVelocity;
	const int			fFirstSpecies;
	const int			fTemperature;
	const int			fVariablesWithoutSpecies;
	char				**fVariableNames;
	Flag				fPrintMolarFractions;
	TMassFractionPtr	fMassFraction;

//	vectors of solution
	VectorPtr			fSolV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSolU;			// len is set to nGridPoints
										// first element is fSol->vec[-1]

//	void	XToEta( TNewtonPtr bt, VectorPtr etaVec );
  	void	InitTCountPremFlamePhys( void );
//	void	SetInitialBC( TGridPtr grid, TInputDataPtr inp );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	void	FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	Double	DiffCorr( int nVariable, NodeInfoPtr nodeInfo );
	Double	SpeciesDiffusion( int nVariable, int speciesIndex, NodeInfoPtr nodeInfo );
	void	FillJacSpeciesDiffusion( int nVariable, int speciesIndex, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	FILE	*GetOutputFile( char *head, char *tail, FileType type );
	void	UpdateDimensions( int len );
	void	UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
	void	UpdateSolutionOnePoint( Double *y, int gridPoint );
};

typedef TCountPremFlamePhys *TCountPremFlamePhysPtr;
