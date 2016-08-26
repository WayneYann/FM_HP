extern "C" int BCLeftNewtonFuncsTUnstr( const VectorPtr x, VectorPtr fVec, void *object );

void UnstrPremPhysJacFirst( void *object, NodeInfoPtr nodeInfo );
void UnstrPremPhysJacRest( void *object, NodeInfoPtr nodeInfo );
void UnstrPremPhysJacLast( void *object, NodeInfoPtr nodeInfo );
void UnstrPremPhysRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
void UnstrPremPhysOutput( void *object, FILE *fp, char* tail );
int UnstrPremPhysPostIter( void *object );
void UnstrPremPhysUpdateLeftBoundary( void  *object );
void UnstrPremPhysUpdateRightBoundary( void *object );
void SetUnstrPremPhysNodeInfo( int k, void *object );
void UnstrPremPhysPostConv( void *object );
ConstStringArray GetUnstrPremPhysVarNames( void *object );

//*************************************************************************************************************
class TUnstrPremFlamePhys : public T1DFlame, public TPremixed {
friend void UnstrPremPhysJacFirst( void *object, NodeInfoPtr nodeInfo );
friend void UnstrPremPhysJacRest( void *object, NodeInfoPtr nodeInfo );
friend void UnstrPremPhysJacLast( void *object, NodeInfoPtr nodeInfo );
friend void UnstrPremPhysRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
friend void UnstrPremPhysOutput( void *object, FILE *fp, char* tail );
friend int UnstrPremPhysPostIter( void *object );
friend void UnstrPremPhysUpdateLeftBoundary( void  *object );
friend void UnstrPremPhysUpdateRightBoundary( void *object );
friend void UnstrPremPhysPostConv( void *object );
friend int BCLeftNewtonFuncsTUnstr( const VectorPtr x, VectorPtr fVec, void *object );

public:
	TUnstrPremFlamePhys( FirstInputPtr firstInp ) : fMassFlowRate(0),
				fFirstSpecies(1),
				fSootMoments(fFirstSpecies + fInputData->GetCounter()->species 
					- fInputData->GetCounter()->steadyStates), 
				fTemperature(fSootMoments+fInputData->fNSootMoments),
				fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
				fDeltaTFix(500.0),
				T1DFlame( firstInp ), TPremixed( fInputData ) 
					{ InitTUnstrPremFlamePhys();};	
	~TUnstrPremFlamePhys( void );

	int		GetOffsetMassFlowRate( void );
	int		GetOffsetTemperature( void );
	int		GetOffsetFirstSpecies( void );
	int		GetVariablesWithoutSpecies( void );
	ConstStringArray	GetVariableNames( void );
// not used by TUnstrPremFlame but required by TFlame as pure virtual functions
	int		GetOffsetVVelocity( void );
	int		GetOffsetUVelocity( void );
	int		GetOffsetMixFrac( void );

	Double				GetdYPrevdY( int speciesIndex, NodeInfoPtr nodeInfo );
	
	VectorPtr	GetM( void ) { return fSolM; };
	void		PrintRHSSpecies( TNewtonPtr bt );
	void		ResetNGridModifications( void ) { fNCutGrid = fNEnlargeGrid = 0; };
	VectorPtr	GetYLeftVec( void ) { return fYLeftVec; };

protected:
	void	PrintRHSTemp( TNewtonPtr bt ) { bt = bt; fprintf( stderr, "nothing happens\n" );};

private:
	const int			fMassFlowRate;
	const int			fFirstSpecies;
	const int			fSootMoments;
	const int			fTemperature;
	const int			fVariablesWithoutSpecies;
	Double				fRadFact;
	char				**fVariableNames;
	Double				fConstMassFlowRate;
	Flag				fConstMassFlux;
	Flag				fPrintMolarFractions;
//	TMassFractionPtr	fMassFraction;
	const Double		fDeltaTFix;		// fixed temperature
	Double				fTfix;			// fixed temperature
	int					fTfixLoc;		// gridpoint of the fixed temperature
	int					fNCutGrid;		// counters for adjusting 
	int					fNEnlargeGrid;	// the computational domain

	NewtonInfoPtr		fNewtonInfoL;
	VectorPtr			fYLeftVec;
	VectorPtr			fRhoY_iV_iPlus;
	Double				*fBC_Left;

//	vectors of solution
	VectorPtr			fSolM;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSavedM;			// length of these vectors is nGridPoints + 2

  	void	InitTUnstrPremFlamePhys( void );
	void	SetInitialBC( TGridPtr grid, TInputDataPtr inp );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	void	SetMassFlux( Double massFlux, MatrixPtr yMat, Double *yLeft, Double *yRight );
	void	FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	Double	DiffCorr( int nVariable, NodeInfoPtr nodeInfo );
	Double	SpeciesDiffusion( int nVariable, int speciesIndex, NodeInfoPtr nodeInfo );
	void	FillJacSpeciesDiffusion( int nVariable, int speciesIndex, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	FILE	*GetOutputFile( char *head, char *tail, FileType type );
	void	UpdateDimensions( int len );
	void	UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
	void	UpdateSolutionOnePoint( Double *y, int gridPoint );
	Double	GetBurningVelocity( void );
	void	PrintRHSSpecies( int start, int end, NodeInfoPtr nodeInfo,  FILE *fp );
	int		CheckComputationalDomain( void );
	Double	SecondDerivXDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	DiffCorrX( int nVariable, NodeInfoPtr nodeInfo );
	void	RestoreSolution( void );
	void	SolutionToSolver( void );
	void	SaveSolution( void );
	void 	CalcAllDiffVeloRhoYiNext( Double *YCurr, Double hNext );
	Double	*GetYLeft( Double *Yguess );
	void	CalcYLeft( void );
	void	CompleteSpeciesDiffusion( NodeInfoPtr nodeInfo );
};

typedef TUnstrPremFlamePhys *TUnstrPremFlamePhysPtr;
