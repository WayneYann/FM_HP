int FlameSheetPostIter( void *object );
void FlameSheetJacRest( void *object, NodeInfoPtr nodeInfo );
void FlameSheetRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
void FlameSheetOutput( void */*object*/, FILE */*fp*/, char* /*tail*/ );
void SetFlameSheetNodeInfo( int k, void *object );
void FlameSheetPostConv( void *object );
ConstStringArray GetFlameSheetVarNames( void *object );

//*************************************************************************************************************
class TFlameSheet : public T1DFlame {
friend void FlameSheetJacRest( void *object, NodeInfoPtr nodeInfo );
friend void FlameSheetRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
friend void FlameSheetOutput( void *object, FILE *fp, char* tail );
friend int FlameSheetPostIter( void *object );

public:
	TFlameSheet( FirstInputPtr firstInp ) : fVVelocity(0), fUVelocity(1), fMixFrac(2), 
//								fTemperature(3), fFirstSpecies(4),
								fFirstSpecies(3), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
								fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
								T1DFlame( firstInp ) { InitTFlameSheet();};	
	~TFlameSheet( void );
	
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

protected:
	void	PrintRHSTemp( TNewtonPtr bt ) { bt = bt;fprintf( stderr, "nothing happens\n" );};


private:
	
	const int			fVVelocity;
	const int			fUVelocity;
	const int			fMixFrac;
	const int			fTemperature;
	const int			fFirstSpecies;
	const int			fVariablesWithoutSpecies;
	char				**fVariableNames;
	TMassFractionPtr	fMassFraction;

//	vectors of solution
	VectorPtr			fSolV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSolU;			// len is set to nGridPoints
	VectorPtr			fSolMixFrac;	// first element is fSol->vec[-1]

 	void	InitTFlameSheet( void );
//	void	SetInitialBC( TGridPtr grid, TInputDataPtr inp );
	void	FillJacMassDiffusion( int nVariable, int nEquation, NodeInfoPtr nodeInfo, Flag sign = kPositive );
	Double 	SecondDerivMassDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	FILE	*GetOutputFile( char *head, char *tail, FileType type );
	void	UpdateDimensions( int len );
	void	UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
};

typedef TFlameSheet *TFlameSheetPtr;
