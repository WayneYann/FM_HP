extern "C" int BCLeftNewtonFuncs( const VectorPtr x, VectorPtr fVec, void *object );
extern "C" int BCRightNewtonFuncs( const VectorPtr x, VectorPtr fVec, void *object );

void DiffPhysEigenJacFirst( void *object, NodeInfoPtr nodeInfo );
void DiffPhysEigenJacRest( void *object, NodeInfoPtr nodeInfo );
void DiffPhysEigenJacLast( void *object, NodeInfoPtr nodeInfo );
void DiffPhysEigenRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
void DiffPhysEigenOutput( void *object, FILE *fp, char* tail );
int DiffPhysEigenPostIter( void *object );
void DiffPhysEigenUpdateLeftBoundary( void  *object );
void DiffPhysEigenUpdateLeftNewton( void  *object );
void DiffPhysEigenUpdateRightBoundary( void  *object );
void SetDiffPhysEigenNodeInfo( int k, void *object );
void DiffPhysEigenPostConv( void *object );
ConstStringArray GetDiffPhysEigenVarNames( void *object );
//*************************************************************************************************************
class TCountDiffPhysEigen : public T1DFlame {
friend void DiffPhysEigenJacFirst( void *object, NodeInfoPtr nodeInfo );
friend void DiffPhysEigenJacRest( void *object, NodeInfoPtr nodeInfo );
friend void DiffPhysEigenJacLast( void *object, NodeInfoPtr nodeInfo );
friend void DiffPhysEigenRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
friend void DiffPhysEigenOutput( void *object, FILE *fp, char* tail );
friend int DiffPhysEigenPostIter( void *object );
friend void DiffPhysEigenUpdateLeftBoundary( void  *object );
friend void DiffPhysEigenUpdateLeftNewton( void  *object );
friend void DiffPhysEigenUpdateRightBoundary( void  *object );
friend void DiffPhysEigenPostConv( void *object );
friend int BCLeftNewtonFuncs( const VectorPtr x, VectorPtr fVec, void *object );
friend int BCRightNewtonFuncs( const VectorPtr x, VectorPtr fVec, void *object );

public:
	TCountDiffPhysEigen( FirstInputPtr firstInp ) : fVVelocity(0), fUVelocity(1), 
				fMixFrac(2), fPStrain(3),fFirstSpecies(4),
				fSootMoments(fFirstSpecies + fInputData->GetCounter()->species 
					- fInputData->GetCounter()->steadyStates), 
				fTemperature(fSootMoments+fInputData->fNSootMoments),
				fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
				T1DFlame( firstInp ) { InitTCountDiffPhysEigen();};	
	~TCountDiffPhysEigen( void );

	int		GetOffsetVVelocity( void );
	int		GetOffsetUVelocity( void );
	int		GetOffsetTemperature( void );
	int		GetOffsetMixFrac( void );
	int		GetOffsetPStrain( void );
	int		GetOffsetFirstSpecies( void );
	int		GetVariablesWithoutSpecies( void );
	ConstStringArray	GetVariableNames( void );
	Double				GetdYNextdY( int speciesIndex, NodeInfoPtr nodeInfo );
	Double				GetdYPrevdY( int speciesIndex, NodeInfoPtr nodeInfo );
	
/*	void SetNewtonLeftInfo( Double *yLeftIn, Double *bcLeftIn );*/
/*	void GetNewtonLeftInfo( Double **yLeftIn, Double **bcLeftIn );*/
	
	VectorPtr	GetV( void ) { return fSolV; };
	VectorPtr	GetG( void ) { return fSolG; };
	VectorPtr	GetZ( void ) { return fSolMixFrac; };
	VectorPtr	GetP( void ) { return fSolP; };
	
	VectorPtr	GetYLeftVec( void ) { return fYLeftVec; };
	VectorPtr	GetYRightVec( void ) { return fYRightVec; };
//	VectorPtr	GetBCLeftVec( void ) { return fBCLeftVec; };
	
	Double				GetStrainRate( void );
	Double				GetSeshaStrainRate( void );
	Double				GetLiquidPoolStrainRate( void );
	TMassFractionPtr	GetMassFraction( void ) { return fMassFraction; };
	void	SaveSolution( void );
	void	RestoreSolution( void );

protected:
	void	PrintRHSTemp( TNewtonPtr bt );
	void	PrintRHSSpecies( TNewtonPtr bt );

private:
	const int			fVVelocity;
	const int			fUVelocity;
	const int			fMixFrac;
	const int			fPStrain;
	const int			fFirstSpecies;
	const int			fSootMoments;
	const int			fTemperature;
	const int			fVariablesWithoutSpecies;
	char				**fVariableNames;
	Flag				fPrintMolarFractions;
	Flag				fLiquidPoolBC;
	NewtonInfoPtr		fNewtonInfoL;
	NewtonInfoPtr		fNewtonInfoR;
	TMassFractionPtr	fMassFraction;
//	int					fMaxUPoint;
	int					fNLeftofStag;

	VectorPtr			fYLeftVec;
	VectorPtr			fYRightVec;
	VectorPtr			fRhoY_iV_iPlus;
//	VectorPtr			fBCLeftVec;
//	Double				*fY_Left;
	Double				*fBC_Left;
	Double				fKappa;
//	vectors of solution
	VectorPtr			fSolV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSolG;			// len is set to nGridPoints
	VectorPtr			fSolMixFrac;	// first element is fSol->vec[-1]
	VectorPtr			fSolP;			// 

//	saved solution
	VectorPtr			fSavedV;			// length of these vectors is nGridPoints + 2
	VectorPtr			fSavedG;			// len is set to nGridPoints
	VectorPtr			fSavedMixFrac;	// first element is fSaved->vec[-1]
	VectorPtr			fSavedP;			// 

	int		ScalarDissipation( TNewtonPtr bt, VectorPtr scalarDissVec, Double *stoechScalarDiss );
  	void	InitTCountDiffPhysEigen( void );
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
	void	UpdateSolutionOnePoint( Double *y, int gridPoint );
	void	FillJacNonlinearConvectUpwind( int nVariable1, int nVariable2, NodeInfoPtr nodeInfo, Double coeff, Flag velocityPositive = TRUE );
	Double	NonlinearConvectUpwind( int nVariable1, int nVariable2, NodeInfoPtr nodeInfo, Flag velocityPositive = TRUE );
	void	PrintRHSSpecies( int start, int end, NodeInfoPtr nodeInfo, FILE *fp );
	void	PrintRHSTemp( NodeInfoPtr nodeInfo, FILE *fp );
	Flag	UpdateBCLiquidPool( void *object );
	Double	GetPolyDiffVeloRhoYiNext( Double *yLeft, int nVariable, Double hNext );
	void	CalcAllDiffVeloRhoYiNext( Double *YCurr, Double hNext );
	void	CalcAllDiffVeloRhoYiPrev( Double *YCurr, Double hPrev );
	Double	SecondDerivBinSpecDiff( int nVariable, NodeInfoPtr nodeInfo );
	void	CalcYLeft( void );
	void	CalcYRight( void );
	Double	*GetYLeft( Double *Yguess );
	Double	*GetYRight( Double *Yguess );
	Double	SecondDerivXDiffusion( int nVariable, NodeInfoPtr nodeInfo );
	Double	DiffCorrX( int nVariable, NodeInfoPtr nodeInfo );
	Double	GetDiffVeloRhoYiNext( Double *YCurr, int nVariable, Double hNext );
	Double	GetVLiquidPool( void );
	void	ComputePropertiesBound( TFlameNodePtr flameNode, Double temp
									, Double *Y, Double pressure );
	void	SolutionToSolver( void );
	void	CompleteSpeciesDiffusion( NodeInfoPtr nodeInfo );
};

typedef TCountDiffPhysEigen *TCountDiffPhysEigenPtr;
