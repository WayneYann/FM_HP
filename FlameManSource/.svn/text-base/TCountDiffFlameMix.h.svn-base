void CountDiffMixJacRest( void *object, NodeInfoPtr nodeInfo );
void CountDiffMixRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
void CountDiffMixOutput( void *object, FILE *fp, char* tail );
int CountDiffMixPostIter( void *object );
void SetCountDiffMixNodeInfo( int k, void *object );
void CountDiffMixPostConv( void *object );
ConstStringArray GetCountDiffMixVarNames( void *object );
void CountDiffMixUpdateBoundary( void  *object );

//*************************************************************************************************************
class TCountDiffFlameMix : public T1DFlame {
friend void CountDiffMixJacRest( void *object, NodeInfoPtr nodeInfo );
friend void CountDiffMixRHSRest( void *object, NodeInfoPtr nodeInfo, RHSMode rhsMode );
friend void CountDiffMixOutput( void *object, FILE *fp, char* tail );
friend int CountDiffMixPostIter( void *object );
friend void CountDiffMixPostConv( void *object );

public:
	TCountDiffFlameMix( FirstInputPtr firstInp ) : 
								fFirstSpecies(0), fTemperature(fFirstSpecies + fInputData->GetCounter()->species - fInputData->GetCounter()->steadyStates),
								fLnChi(fTemperature+1),
								fSootMoments(fLnChi+1),
								fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
								T1DFlame( firstInp ) { InitTCountDiffFlameMix();};	
	~TCountDiffFlameMix( void );

	int			GetOffsetVVelocity( void );
	int			GetOffsetUVelocity( void );
	int			GetOffsetMixFrac( void );
	int			GetOffsetFirstSpecies( void );
	int			GetOffsetTemperature( void );
	int			GetVariablesWithoutSpecies( void );
	int			GetOffsetLnChi( void );
	ConstStringArray	GetVariableNames( void );
	VectorPtr	GetDissRateVector( void ) { return fDissRate; };
//	Double		GetDissRate( void ) { return fDissRate->vec[fDissRate->len]; };
	Double		GetDissRate( void ) { return (fArcLengthContin) ? exp( fSolLnChi->vec[1] ) : fDissRate->vec[fDissRate->len]; };
	Double		GetDissRate( Double z );
	Double		DissRateFact( Double rho );
//	Double		GetDissRate( Double z, Double rho );
	Double		GetDissRate( Double z, Double rho, Double chiStoich );
	int			GetZRefLoc( VectorPtr xVec );
	Double		GetZRef( void );

// continuation stuff
//	void		IncTempContStart( Double inc ) { fTempContStart += inc; };
	void		SetTempContStart( Double val ) { fTempContStart = val; };
	void		SetChiContStart( Double val ) { fLnChiContStart = val; };
	void		SetDeltaArcLength( Double val ) { fDeltaArcLength = val; };
	void		SetTMaxLoc( int loc ) { fTmaxLoc = loc; };
	void		SetdTds( Double dTds ) { fdTds = dTds; };
	void		SetdlnChids( Double dlnChids ) { fdlnChids = dlnChids; };
	void		IncNFlameletCount( void ) { ++fNFlameletsCount; };
	void		ReInitArcLengthCont( void ) { fDeltaArcLength = 0.0; fdTds = 0.0; fdlnChids = 0.0; fNFlameletsCount = 0; };
	Flag		GetArcLengthCont( void ) { return fArcLengthContin; };
	Flag		GetArcUp( void ) { return fArcUp; };
	int			GetTMaxLoc( void ) { return fTmaxLoc; };
	int			GetMaxFlamelets( void ) { return fMaxFlamelets; };
	int			GetNFlameletsCount( void ) { return fNFlameletsCount; };
	Double		GetdTds( void ) { return fdTds; };
	Double		GetdlnChids( void ) { return fdlnChids; };
	Double		GetChiContStart( void ) { return fLnChiContStart; };
	Double		GetTempContStart( void ) { return fTempContStart; };
	Double		GetDeltaArcLength( void ) { return fDeltaArcLength; };
	Double		GetDeltaChiref( void ) { return fDeltaChiRef; };
	Double		GetDeltaTref( void ) { return fDeltaTRef; };
	int			GetLnChi( void ) { return fLnChi; };

protected:
	void	PrintRHSTemp( TNewtonPtr /*bt*/ ) { cerr << "nothing happens" << NEWL; };

private:
	const int	fFirstSpecies;
	const int	fTemperature;
	const int	fLnChi;
	const int	fSootMoments;
	const int	fVariablesWithoutSpecies;
	char		**fVariableNames;
	Flag		fPrintMolarFractions;
	VectorPtr	fSolLnChi;	// 
	VectorPtr	fSavedLnChi;		// 
	Double		fRhoInf;
	Double		fRhoRef;
	VectorPtr	fZCount;
	VectorPtr	fChiCount;
	VectorPtr	fDissRate;

// continuation stuff
	int			fTmaxLoc;			//	location of maximum temperature
	int			fMaxFlamelets;
	int			fNFlameletsCount;
	Flag		fArcLengthContin;
	Flag		fArcUp;
	Double		fTempContStart;
	Double		fLnChiContStart;
	Double		fSavedTempContStart;
	Double		fSavedLnChiContStart;
	Double		fDeltaArcLength;
	Double		fdTds;				//	dT_max/ds, where s is the arclength
	Double		fDeltaChiRef;
	Double		fDeltaTRef;
	Double		fdlnChids;			//	da_inv/ds, where s is the arclength and a_inv the inverse of the strainrate

// Liquid fuel stuff
	Double		fCl; // cp of liquid
	Double		fHv; // enthalpy of evaporation at 1atm and boiling conditions
	Double		fT_B;// boiling temperature 

  	void	InitTCountDiffFlameMix( void );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	FILE	*GetOutputFile( char *head, char *tail, FileType type );
	void	SetInitialBC( TGridPtr grid, TInputDataPtr inp );
	void	FillJacDiffCorr( int nVariable, Double constCoeff, NodeInfoPtr nodeInfo, Flag sign );
	void	UpdateDimensions( int len );
	void	SaveSolution( void );
	void	RestoreSolution( void );
	void	SolutionToSolver( void );
	void	UpdateSolution( MatrixPtr yMat, VectorPtr yLeftVec, VectorPtr yRightVec );
	void	UpdateSolutionOnePoint( Double *y, int gridPoint );
	Double	ExactChi( Double Z );
	Double	Interpol( Double x, Double xOld, Double yOld, Double xNew, Double yNew );
	void    PrintProdRate( int speciesIndex, FILE *fp );
};

typedef TCountDiffFlameMix *TCountDiffFlameMixPtr;
