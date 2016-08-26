 void ResTrans1DIsoChorImpliSolver( Double *T, Double *y, Double *yPrime, Double *delta
				, int *iRes, Double *rPar, int *iPar );

class TTrans1DIsoChorSolver : public T0DFlame, public TPremixed {
friend void ResTrans1DIsoChorImpliSolver( Double *T, Double *y, Double *yPrime, Double *delta
				, int *iRes, Double *rPar, int *iPar );
public:
	TTrans1DIsoChorSolver( FirstInputPtr firstInp ) : fPsi(0), fR(1), fPressure(2), fTemperature(3), 
			fFirstSpecies(4), 
			fVariables(fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates + fInputData->fVariablesWithoutSpecies),
			fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
			fNOfEquations(fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates + fInputData->fVariablesWithoutSpecies ),
			T0DFlame( firstInp ), TPremixed( fInputData ) { InitTTrans1DIsoChorSolver();};	
	virtual ~TTrans1DIsoChorSolver( void );

	int			GetOffsetFirstSpecies( void ) { return fFirstSpecies; };
	int			GetOffsetTemperature( void ) { return fTemperature; };
	int			GetVariablesWithoutSpecies( void ) { return fVariablesWithoutSpecies; };
	ConstStringArray	GetVariableNames( void ) { return fVariableNames; };

	void		Initialize( Double timeStart
					, ConstStringArray names
					, Double **startSolution, int vars
					, Double *grid, int gridPointsA
					, Double pressureStart, Double scalarDissRateStart
					, Double firstTimeStep );
	void		MakeGrid( VectorPtr theGrid, Double left, Double right, Flag equidistant );
	void		MakeGrid( Double *grid, int gridPoints
							, Double left, Double right, Flag equidistant );
	void		PrintSolution( FILE *fp, int nGPoints, int nEq
									, Double **sol, ConstStringArray names );
	Flag		Solve( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, int deltaStepsOut );
	void		SetEndValues( Double timeEnd, Double pressureEnd
					, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd );
	void		WriteFlameletFile( FILE *fp, char *head, char *tail );
	void		WriteFlameletFile( Double time, Double **sol, Double *x, 
											FILE *fp, char *head, char *tail );
	FILE 		*GetOutputFile( Double time, const char *head, const char *tail, FileType type );
	FILE		*GetOutputFile( char *head, char *tail, FileType type );
	Double		GetCurrentTime( void );
	void 		PrintProdRateGlobalReac( Double time );

protected:
	Double	Interpol( Double t, Double valOld, Double tOld, Double valNew, Double tNew );
	Double	ComputeEmissionIndex( int speciesIndex );

	const int	fFirstSpecies;
	const int	fTemperature;
	const int	fPressure;
	const int	fR;
	const int	fPsi;
	const int	fVariablesWithoutSpecies;
	const int	fNOfEquations;
	const int	fVariables;

private:
  	void	InitTTrans1DIsoChorSolver( void );
	void	SetInitialConditions( Double *y, TInputDataPtr inp );
	void 	InitDassl( int nGridPoints );
	void 	InitDasslWorkSpace( int i );
	void 	SetInfo( int i );
	int		GetActualPoint( Double tEnd );
	void	SetInitialPsi( void );
	void	SaveSolution( int k, Double time, Double *y );
	void	SetFDWeights( Double *sol );
	void	ComputeGPDistributionWeights( Double *sol );
	void	SetInitial( Double **startSol, int gridPoints, int vars );
	void	ResTrans1DIsoChorImpliSolver( Double */*T*/, Double *y, Double *yPrime
				, Double *delta
				, int */*iRes*/, Double */*rPar*/, int */*iPar*/ );
	int		GetVariableIndex( const char *name );
	int		GetVariableIndex( const char *name, ConstStringArray array, int len );
	void	SetMaxTimeStep( int kAct, Double t );
	Double	GetTempSource( Double r, Double time );
	Double	GetPressure( void ) ;
	void	PrintProdRate( int speciesIndex, FILE *fp );

	void	DoExit( void );
	void	SetOutSolution( void );
	Flag	FirstImpliStep( void );
	Flag	OneImpliStep( void );

	enum FortranInd	{ kF1, kF2, kF3, kF4, kF5, kF6, kF7, kF8, kF9, kF10, kF11
					, kF12, kF13, kF14, kF15, kF16, kF17, kF18, kF19, kF20 };

	// input
	Double		fTStart;
	Double		fTEnd;
	Double		fRLeft;			// left boundary
	Double		fRRight;		// right boundary
	Double		fPsiLeft;		// left boundary
	Double		fPsiRight;		// right boundary
	Double		fFirstTimeStep;
	Double		*fRTol;
	Double		*fATol;
	int			fNOfSolver;
	int			fNGridPoints;
	Flag		fEquidistant;
	Flag		fPrintMolarFractions;
	Double		fDeltaTMax;

	Flag		fFirstCall;

	FILE		*fOutFilePtr;

	Double		fPressStart;
	Double		fPressEnd;
	Double		fPhiStart;
	Double		fPhiEnd;
	Double		fTempLeftStart;
	Double		fTempRightStart;
	Double		fTempLeftEnd;
	Double		fTempRightEnd;
	Double		fRLeftStart;
	Double		fRRightStart;
	Double		fRLeftEnd;
	Double		fRRightEnd;
	Double		fDPdt;

	int			fActualPoint;
	int			fMaxOrd;
	char		**fVariableNames;
	Double		*fMaxStepTaken;
	Double		fMaxStepSize;
	int			**fNActualOrd;
	int			**fNActualStep;
	Double		fTout;

	Double		*fh;
	Double		*fhm;
	Double		*fhnenn;
	Double		*fFDWCurr;
	Double		*fFDWPlus;
	Double		*fFDWMinus;
	Double		*fWPlus;
	Double		*fWMinus;
	Double		*fGPDistWeight;
	VectorPtr	fGPDens;
	VectorPtr	fSolTime;
	MatrixPtr	fSolMassFracs;
	VectorPtr	fSolTemp;
	VectorPtr	fSolPressure;
	VectorPtr	fSolPsi;
	VectorPtr	fSolR;
	VectorPtr	fSolOldTime;
	MatrixPtr	fSolOldMassFracs;
	VectorPtr	fSolOldTemp;
	VectorPtr	fSolOldPressure;
	VectorPtr	fSolOldPsi;
	VectorPtr	fSolOldR;
	int			fGamma;
	Double		fKappa;
	Double		fTauGrid;

	Double		*fRho;
	Double		*fMixHeatCap;
	Double		fArtificialSource;

	MatrixPtr	fMassFracsWork;
	VectorPtr	fTempWork;
	VectorPtr	fPressureWork;
	VectorPtr	fPsiWork;
	VectorPtr	fRWork;
	VectorPtr	fOutSolWork;

//	variables used for dassl

	IntVectorPtr	*fInfo;		// length is fInfo[nGridPoints]->vec[16]
	VectorPtr		fTime;		// workspace for vector of solution; fTime->vec[fNOfEquations]
	MatrixPtr		fSolution;	// workspace for vector of solution
	MatrixPtr		fSolPrime;	// workspace for partial derivatives
								// lenght is fSolution->mat[nGridPoints][fNOfEquations]
	int				fIdid;
	int				fLRW;
	VectorPtr		*fRWrk;	// workspace for ddassl
								// lenght is fRWrk[nGridPoints]->vec[fLRW]
	int				fLIW;
	IntVectorPtr	*fIWrk;	// workspace for ddassl
								// lenght is fIWrk[nGridPoints]->vec[fLIW]
	int				fDasslNEq;	// = fNOfEquations, but dassl needs non constant
	Double			**fDmDy;	// dMdY[i][j] = dM_j/dY_i; dMdY[nSpeciesIn+1][nSpeciesIn]
};

typedef TTrans1DIsoChorSolver *TTrans1DIsoChorSolverPtr;

