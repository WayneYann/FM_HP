#include <cvode/cvode.h>
#include <cvode/cvode_band.h>
#include <nvector/nvector_serial.h>

#define getflameletsolver getflameletsolver_
#define initflamelet initflamelet_
#define initflameletzref initflameletzref_
#define getspeciesdata getspeciesdata_
#define solveflameletzr solveflameletzr_
#define solveflamelet solveflamelet_
#define getflameletsolution getflameletsolution_
#define getsolutioninfo getsolutioninfo_
#define disposeflamelet disposeflamelet_
#define printsolution printsolution_
#define getdissrate getdissrate_
#define getdissfunc getdissfunc_
#define ADIABFLAMETEMP adiabflametemp_
#define PDFZBETA pdfzbeta_
#define F77BETAPDF f77betapdf_

extern "C" {
	void GETFLAMELETSOLVER( int **object, int *objectSize, char *fInputName, int *nChars );
	void INITFLAMELET( int **object, Double *timeStart
					, char *names, int *nChars
					, Double *startSolution, int *vars, int *offsetVars
					, Double *grid, int *gridPoints, int *offsetGridPoints
					, Double *pressureStart, Double *scalarDissRateStart
					, Double *firstTimeStep );
	void INITFLAMELETZREF( int **object, Double *timeStart
					, char *names, int *nChars
					, Double *startSolution, int *vars, int *offsetVars
					, Double *grid, int *gridPoints, int *offsetGridPoints
					, Double *pressureStart, Double *scalarDissRateStart
					, Double *firstTimeStep, Double *zRef );
	void GETSPECIESDATA( int **object, char *names, int *nChars
					, int *vars, int *offsetVars
					, Double *high, Double *low, Double *molarMass );
	int SOLVEFLAMELETZR( int **object, Double *timeEnd, Double *pressureEnd
						, Double *scalarDissRateEnd
						, Double *temOxEnd, Double *tempFuelEnd, Double *zR );
	int SOLVEFLAMELET( int **object, Double *timeEnd, Double *pressureEnd
					, Double *scalarDissRateEnd
					, Double *temOxEnd, Double *tempFuelEnd );
	void GETFLAMELETSOLUTION( int **object, char *names, int *nChars, Double *outSol, Double *grid
					, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars );
	void GETFLAMELETSOLUTIONEXT( int **object, char *names, int *nChars, Double *outSol, Double *grid
						, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars, Double *density );
	void GETSOLUTIONINFO( int **object, Double *timeStep );
	void DISPOSEFLAMELET( int **object );
	void PRINTSOLUTION( void **object, Double *timeCurr, char *names, int *nChars
			, Double *outSol, Double *grid
			, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars );
	void GETDISSRATE( void **object, Double *dissRateReq, Double *ZReq, Double *dissRate, Double *z );
	int GETSOOTSOURCES( int **object, int *whichMoment, char *names, int *nChars, Double *outSol, Double *grid
					, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars );
	void WRITEFLAMELETFILE( int **object );
	void MAKEGRID( int **object, Double *grid, int *gridPoints
			, Double *left, Double *right, int *equidistant );
	void GETSPECIESNAMES( int **object, char *names, int *nChars, int *len );
	void GETFUEL( int **object, char *name, int *nChars );
	void GETNOFSPECIES( int **object, int *nSpecies, int *nSpeciesIn );
	void GETDISSFUNC( void **object, Double *a, Double *z );
	int GETMAXSPECNAMELEN( int **object );
	void READSTARTPROF( int **object
						, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars
						, Double *grid, Double *startSolution, Double *pressure, Double *chi
						, Double *theTime, Double *currTimeStep
						, int *tempOff, int *speciesOff, int *sootOff );
	void ADIABFLAMETEMP(Double *z, int *komp, Double *c0, Double *temp0, Double *p, char *symbol, 
                                 Double *xmol_i, Double *Enth0,Double *ponalo, Double *ponahi,
                                 int *nndd,
                                 Double *adiaT,Double *adiaY);
	void PDFZBETA( Double *ZMEAN, Double *ZVAR, Double *Z, Double *PDFZ, int *NJ );
	void F77BETAPDF( Double *pdfz, Double *zgrid, Double *zmean, Double *zvari, int *nz, int *nzmax );

	void getflameletsolver( int **object, int *objectSize, char *fInputName, int *nChars );
	void initflamelet( int **object, Double *timeStart
					, char *names, int *nChars
					, Double *startSolution, int *vars, int *offsetVars
					, Double *grid, int *gridPoints, int *offsetGridPoints
					, Double *pressureStart, Double *scalarDissRateStart
					, Double *firstTimeStep );
	void initflameletzref( int **object, Double *timeStart
					, char *names, int *nChars
					, Double *startSolution, int *vars, int *offsetVars
					, Double *grid, int *gridPoints, int *offsetGridPoints
					, Double *pressureStart, Double *scalarDissRateStart
					, Double *firstTimeStep, Double *zRef );
	void getspeciesdata( int **object, char *names, int *nChars
					, int *vars, int *offsetVars
					, Double *high, Double *low, Double *molarMass );
	int solveflameletzr( int **object, Double *timeEnd, Double *pressureEnd
						, Double *scalarDissRateEnd
						, Double *temOxEnd, Double *tempFuelEnd, Double *zR );
	int solveflamelet( int **object, Double *timeEnd, Double *pressureEnd
					, Double *scalarDissRateEnd
					, Double *temOxEnd, Double *tempFuelEnd );
	void getflameletsolution( int **object, char *names, int *nChars, Double *outSol, Double *grid
					, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars );
	void getflameletsolutionext( int **object, char *names, int *nChars, Double *outSol, Double *grid
						, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars, Double *density );
	void getsolutioninfo( int **object, Double *timeStep );
	void disposeflamelet( int **object );
	void printsolution( void **object, Double *timeCurr, char *names, int *nChars
			, Double *outSol, Double *grid
			, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars );
	void getdissrate( void **object, Double *dissRateReq, Double *ZReq, Double *dissRate, Double *z );
	int getsootsources( int **object, int *whichMoment, char *names, int *nChars, Double *sootSources, Double *grid
					, int *gridPoints, int *offsetGridPoints, int *nSources, int *offsetSources );
	void writeflameletfile( int **object );
	void makegrid( int **object, Double *grid, int *gridPoints
			, Double *left, Double *right, int *equidistant );
	void getspeciesnames( int **object, char *names, int *nChars, int *len );
	void getfuel( int **object, char *name, int *nChars );
	void getnofspecies( int **object, int *nSpecies, int *nSpeciesIn );
	void getdissfunc( void **object, Double *a, Double *z );
	int getmaxspecnamelen( int **object );
	void readstartprof( int **object
						, int *gridPoints, int *offsetGridPoints, int *vars, int *offsetVars
						, Double *grid, Double *startSolution, Double *pressure, Double *chi
						, Double *theTime, Double *currTimeStep
						, int *tempOff, int *speciesOff, int *sootOff );
/*	void adiabflametemp(Double *z, int *komp, Double *c0, Double *temp0, Double *p, char *symbol, */
/*                                 Double *xmol_i, Double *Enth0,Double *ponalo, Double *ponahi,*/
/*                                 int *nndd,*/
/*                                 Double *adiaT,Double *adiaY);*/
/*	void pdfzbeta( Double *ZMEAN, Double *ZVAR, Double *Z, Double *PDFZ, int *NJ );*/
}
/*  

	INITFLAMELET	-names				name of variables (e.g.{"T","NO","CO"})
					-nChars				length of each variable name
					-startSolution		start solution (typically air on one side
										fuel on the other side at t=0)
										2-D-Array (number of rows=number of vars, 
										number of columns-1=number of species)
										Fortran: real*8 flamelet(nspec+1,nx)
										(nspec+1) to hold temperature
					-vars				number of vars
					-offsetVars			physical length of var storage (rows) 
					-offsetGridPoints	physical length of grid point storage (columns)
					-firstTimeStep		guess for first time step (=0 don't know)
										

	INITFLAMELETZREF					same as INITFLAMELET except
					-zRef				location in mixture fraction space, where
										scalar dissipation rate is given in calls 
										of INITFLAMELETZREF and SOLVE
										

	GETSPECIESDATA	-names				cf. above
					-nChars				cf. above					
					-vars				cf. above
					-offsetVars			cf. above
					-high				nasa polynomials to calculate enthalpy, cp etc.
					-low				(high and low temperature)
					-molarMass
					
	MAKEGRID( int **object, Double *grid, int *gridPoints
			, Double *left, Double *right, int *equidistant )

				yields grid from 'left' to 'right'
				consisting of 'gridPoints' gridpoints including boundaries
				grid is equidistant, if 'equidistant' is 1,
				nonequidistant, 'equidistant' is 0.
				If 'equidistant'  and 'gridPoints' are consistent with 
				corresponding variables in
				the inputfile, the function yields the grid that is used in 
				the flamelet solver
							


	new flamelet with calling sequence:

	GETFLAMELETSOLVER(...);
	INITFLAMELET(...);	
	
	throw away:
	DISPOSEFLAMELET(...);
			

*/



//*************************************************************************************************************
void ResTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar );
void ResTransFlameImpliSolver( Double *T, Double *y, Double *yPrime, Double *delta
				, int *iRes, Double *rPar, int *iPar );
int ResTransFlameImpliSolverCV( realtype T, N_Vector u, N_Vector udot, void *data );
void JacTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *pd
			, Double *cj, Double *rPar, int *iPar );

class TTransFlameSolver : public T0DFlame {
friend void ResTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar );
friend void ResTransFlameImpliSolver( Double *T, Double *y, Double *yPrime, Double *delta
				, int *iRes, Double *rPar, int *iPar );
friend int ResTransFlameImpliSolverCV( realtype T, N_Vector u, N_Vector udot, void *data );
friend void JacTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *pd
			, Double *cj, Double *rPar, int *iPar );
public:
	TTransFlameSolver( FirstInputPtr firstInp ) : 
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
            fGrid(0),
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
			fFirstSpecies(1), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
			fVariables(fInputData->GetCounter()->species+fInputData->fVariablesWithoutSpecies+1),
			fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies+1),
			fSootMoments(fTemperature+1),
			fNOfEquations(fInputData->GetCounter()->species+fVariablesWithoutSpecies - fInputData->GetCounter()->steadyStates),
			T0DFlame( firstInp ) { InitTTransFlameSolver();};	
	virtual ~TTransFlameSolver( void );

	int			GetOffsetFirstSpecies( void ) { return fFirstSpecies; };
	int			GetOffsetTemperature( void ) { return fTemperature; };
	int			GetVariablesWithoutSpecies( void ) { return fVariablesWithoutSpecies; };
	ConstStringArray	GetVariableNames( void ) { return fVariableNames; };

	void		Initialize( Double timeStart
					, ConstStringArray names
					, Double **startSolution, int vars
					, Double *grid, int gridPointsA
					, Double pressureStart, Double scalarDissRateStart
					, Double firstTimeStep
					, Double ZRef
					, Double ZRStart );
	void		GetSpeciesData( char **newNames, MatrixPtr highMat
					, MatrixPtr lowMat, Double *molarMass );
	void		GetSolution( ConstStringArray names, Double **outSol, Double *grid
					, int gridPointsA, int vars, Double *density = NULL );
	void		GetSolutionInfo( Double *timeStep );
	void		SetSootSources( int whichMoment, MatrixPtr sourcesMat
					, Double theTime, Double *temp, Double **Y
					, Double pressure, Double **moments );
	int			GetSootSources( int whichMoment, ConstStringArray names, Double **sources, Double *grid, int gridPointsA, int vars );
	void		MakeGrid( VectorPtr theGrid, Double left, Double right, Flag equidistant );
	void		MakeGrid( Double *grid, int gridPoints
							, Double left, Double right, Flag equidistant );
	void		PrintSolution( FILE *fp, int nGPoints, int nEq, Double *grid
									, Double **sol, ConstStringArray names );
	Flag		Solve( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, int deltaStepsOut );
	Flag		Solve( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, Double ZREnd, int deltaStepsOut );
	void		SetEndValues( Double timeEnd, Double pressureEnd, Double scalarDissRateEnd
					, Double temOxEnd, Double tempFuelEnd, Double ZREnd );
	void		PrintDissRate( Double t );
	void		WriteFlameletFile( FILE *fp, char *head, char *tail );
	void		WriteFlameletFile( Double time, Double **sol, Double *x, 
											FILE *fp, char *head, char *tail );
	FILE 		*GetOutputFile( Double time, const char *head, const char *tail, FileType type );
	FILE		*GetOutputFile( char *head, char *tail, FileType type );
	Double		GetCurrentTime( void );
	Double 		GetDissRateReq( Double ZReq, Double dissRate, Double z );
	void 		PrintProdRateGlobalReac( Double time );
	Double		GetRandomNumber( void ) { return fRandomNumber; };
	Double		GetDissFunc( Double z );
	void		ReadStartProfiles( TInputDataPtr inp
				, int nGridPoints, int nVars, Double *x, Double **y, Double *pressure
				, Double *chi, Double *theTime, Double *currTimeStep
				, int tempOff, int speciesOff, int sootOff );

protected:
	Double	GetTempOfEnthalpy( Double ent, Double *Y, Double initialGuess );
	Double	Interpol( Double t, Double valOld, Double tOld, Double valNew, Double tNew );
	Double	ComputeEmissionIndex( Double time, int speciesIndex, Double **PDF
								, Double *timePDF );
	Double	ComputeEmissionIndexSoot( Double time, int which, Double **PDF
								, Double *timePDF );
	Double	TurbMeanX( Double t, int speciesIndex, Double ZMean, Double ZVar, Double *specVar = NULL );
	Double	TurbMeanY( Double t, int speciesIndex, Double ZMean, Double ZVar, Double *specVar = NULL );
	Double	TurbMeanTemp( Double t, Double ZMean, Double ZVar, Double *tempVar = NULL );
	Double	TurbMeanTotEnergy( Double t, Double ZMean, Double ZVar );
	Double	TurbMeanZBarlow( Double t, Double ZMean, Double ZVar );
	Double	TurbMeanSoot( Double t, int sootIndex, Double ZMean, Double ZVar, Double *sootVar );

//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	const int	fGrid;
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	const int	fFirstSpecies;
	const int	fTemperature;
	const int	fVariablesWithoutSpecies;
	const int	fNOfEquations;
	const int	fVariables;
	const int	fSootMoments;

private:
  	void	InitTTransFlameSolver( void );
	void	SetInitialConditions( Double *y, TInputDataPtr inp );
	void 	InitDassl( int nGridPoints );
	void 	InitDasslWorkSpace( int i );
	void 	SetInfo( int i );
	Flag	FirstStep( void );
	Flag	OneStep( int k );
	int		GetActualPoint( Double tEnd );
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	void	ComputeGPDistributionWeights( Double *y );
	void	SetFDWeights( Double *y );	
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo 
	void	SaveSolution( int k, Double time, Double *y );
	void	SetWeights( void );
	void	SetMaxVals( void );
	void	SetInitial( ConstStringArray names, Double **startSol, Double *grid
						, int gridPoints, int vars );
	Double	GetTotEnt( Double temp, Double *Y );
	Double	GetTotEnt( int k, Double *Z, Double t );
	Double	GetZStoi( void );
	Double	GetDissRate( Double time, Double z );
	Double	GetDissRate( Double t, Double z, Double rho );
	Double	ExactChi( Double Z );
	Double	DissRateFact( Double t, Double rho );
	Double	GetRefDissRate( Double time );
	Double	GetZR( void ) { return Interpol( GetCurrentTime(), fZRStart, fTStart, fZREnd, fTEnd ); };
	Double	GetDelQ( Double t, Double Z );
	void	SetMolarMassOverRInf( void );
	void 	ResTransFlameSolver( Double *T, Double *y, Double *yPrime, Double *delta
						, int *iRes, Double *rPar, int *iPar );
	void	ResTransFlameImpliSolver( Double */*T*/, Double *y, Double *yPrime
				, Double *delta
				, int */*iRes*/, Double */*rPar*/, int */*iPar*/ );
	void	JacTransFlameSolver( Double *T, Double *y, Double *yPrime
			, Double **pd, Double cj, Double *rPar, int *iPar );
	int		GetVariableIndex( const char *name );
	int		GetVariableIndex( const char *name, ConstStringArray array, int len );
	void	SetMaxTimeStep( int kAct, Double t );
	Double	GetTempSource( Double Z );
	Double	GetPressure( void ) 
			{ cerr << "#error: function GetPressure not allowed in class " << 
					"TTransFlameSolver" << NEWL; return -1.0; };
	Double	GetPressure( Double time );
	void	PrintProdRate( Double time, int speciesIndex, FILE *fp );
	void	WriteSootInfo( Double theTime, Double *temp, Double **Y
				, Double pressure, Double **moments, FILE *fp );

	void	DoExit( void );
	void	SetOutSolution( void );
	Flag	FirstImpliStep( void );
	Flag	OneImpliStep( void );
	void	SetRandomNumber( void );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp
				, int nGridPoints, int nVars, Double *x, Double **y, Double *pressure, Double *chi
				, Double *theTime, Double *currTimeStep
				, int tempOff, int speciesOff, int sootOff );
	void	SetInitialValues( int nGridPoints, int nVars, Double *x, Double **y );
	void	PostIter( Double t );
    Double  GetRosseRadiation( int k, Double *nTemp, Double **moments, Double rfPrev
                                            , Double rfCurr, Double rfNext );
	void	UpdateThermoProps( int k, Double *Y, Double temp, Double &pressure
								, Double &density, EqOfState what, Double *sootMoments );
	void	ComputeDeltaI( int k, Double *Y, Double temp );
	Double	CompOneDeltaI( int i, int k, Double *Y );
	void	CompDeltaIG( int k, Double temp );
	void	SetGrid( void );
	Double	GetU( Double t, Double Z );
	Flag	CheckOutput( Double *t );
	int		InitCVODE( void );
	void	SolutionToCVode( void );
	void	SolutionFromCVode( void );

	enum FortranInd	{ kF1, kF2, kF3, kF4, kF5, kF6, kF7, kF8, kF9, kF10, kF11
					, kF12, kF13, kF14, kF15, kF16, kF17, kF18, kF19, kF20 };

	// input
	Double		fTStart;
	Double		fTEnd;
	Double		fZl;	// left boundary
	Double		fZr;	// right boundary
	Double		fFirstTimeStep;
	Double		*fRTol;
	Double		*fATol;
	int			fNOfSolver;
	int			fNGridPoints;
	Flag		fEquidistant;
	Flag		fPrintMolarFractions;
	Double		fDeltaTMax;
	VectorPtr	fZRin;
	VectorPtr	fTimeIn;
	VectorPtr	fZIn;
	MatrixPtr	fChiIn;
	MatrixPtr	fUstOverU;
	MatrixPtr	fDelQ;
	VectorPtr	fZCount;
	VectorPtr	fChiCount;
	

	Flag		fFirstCall;
	Double		fZRef;
	Double		fWOverRInf;

	FILE		*fOutFilePtr;
	FILE		*fCAinFile;

	Double		fPressStart;
	Double		fChiStart;
	Double		fTempOxStart;
	Double		fTempFuelStart;
	Double		fZRStart;
	Double		fPressEnd;
	Double		fChiEnd;
	Double		fTempOxEnd;
	Double		fTempFuelEnd;
	Double		fZREnd;
	Double		fRandomNumber;
	Double		fDTdtOx;
	Double		fDTdtFu;
	Double		fDPdt;
	Double		fdDeltaZdt;
	VectorPtr	fTotEntStart;
	VectorPtr	fTotEntEnd;

	MatrixPtr	fSpecHeatCp;
	MatrixPtr	fSpecEnthalpy;
	MatrixPtr	fProdRate;
	VectorPtr	fHeatCpMix;
	VectorPtr	fViscosity;
	VectorPtr	fMolarMassMix;
	VectorPtr	fDensity;
	VectorPtr	fLambdaOverCpMix;

	VectorPtr	fDiffTermY;
	VectorPtr	fDiffTermW;
	VectorPtr	fDiffCorrY;
	VectorPtr	fDiffCorrW;
	
	int			fActualPoint;
	int			fMaxOrd;
	char		**fVariableNames;
	int			**fNActualOrd;
	int			**fNActualStep;
	Double		*fActualTimeStepSize;
	Double		fTout;
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	Double		*fh;
	Double		*fhm;
	Double		*fMonFct;
	Double      fKappa;
	Double		fTauGrid;
	VectorPtr	fgpDens;
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo 
	VectorPtr	fSolGrid;
	Double		*fFDWCurr;
	Double		*fFDWPlus;
	Double		*fFDWMinus;
	Double		*fWCurr;
	Double		*fWPlus;
	Double		*fWMinus;
	VectorPtr	fSolTime;
	MatrixPtr	fSolMassFracs;
	MatrixPtr	fSolSootMoments;
	VectorPtr	fSolTemp;
	VectorPtr	fSolOldTime;
	MatrixPtr	fSolOldMassFracs;
	MatrixPtr	fSolOldSootMoments;
	VectorPtr	fSolOldTemp;
	VectorPtr	fMaxVals;

	MatrixPtr	fMassFracsWork;
	MatrixPtr	fSootMomentsWork;
	VectorPtr	fTempWork;
	VectorPtr	fOutSolWork;
	Double		*fTempGSave;
	Double		**fYGSave;
	Double		**fDeltaI;
	Double		***fG_ij;

// cvode stuff
	int			fNActualOrdCV;
	int			fNActualStepCV;
	Double		fActualTimeStepSizeCV;
	void		*fMem;
	N_Vector	fCVY;       // solution vector space
	realtype	*fCYdata;  // pointer to access data in N_Vector

//	variables used for dassl

	IntVectorPtr	*fInfo;		// length is fInfo[nGridPoints]->vec[16]
	VectorPtr		fTime;		// workspace for vector of solution; fTime->vec[fNOfEquations]
	MatrixPtr		fSolution;	// workspace for vector of solution
	MatrixPtr		fSolPrime;	// workspace for partial derivatives
								// lenght is fSolution->mat[nGridPoints][fNOfEquations]
	int				fIdid;
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	int				fML;
//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo//ingo
	int				fLRW;
	VectorPtr		*fRWork;	// workspace for ddassl
								// lenght is fRWork[nGridPoints]->vec[fLRW]
	int				fLIW;
	IntVectorPtr	*fIWork;	// workspace for ddassl
								// lenght is fIWork[nGridPoints]->vec[fLIW]
	int				fDasslNEq;	// = fNOfEquations, but dassl needs non constant
	Double			**fDmDy;	// dMdY[i][j] = dM_j/dY_i; dMdY[nSpeciesIn+1][nSpeciesIn]
};

typedef TTransFlameSolver *TTransFlameSolverPtr;

