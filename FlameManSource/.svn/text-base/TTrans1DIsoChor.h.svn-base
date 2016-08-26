#include "TTrans1DIsoChorSolver.h"

//*************************************************************************************************************
class TTrans1DIsoChor : public TTrans1DIsoChorSolver {
public:
	TTrans1DIsoChor( FirstInputPtr firstInp ) : fPsi(0), fR(1), fPressure(2), fTemperature(3), 
			fFirstSpecies(4),
			fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
			TTrans1DIsoChorSolver( firstInp ) { InitTTrans1DIsoChor();};	
	virtual ~TTrans1DIsoChor( void );


	void				Solve( void );
	ConstStringArray	GetVariableNames( void ) { return fVariableNames; };

private:
  	void	InitTTrans1DIsoChor( void );
	void	SetInitialBC( void );
	void	SetInitialValues( void );
	Double	GetTempEnd( Double time, Double timeStart, Double timeCut
									, Double tempStart, Double tempEnd );
	void	ReadStartProfiles( TInputDataPtr inp );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	void	WriteScalars( Double time );

	const int	fR;
	const int	fPsi;
	const int	fPressure;
	const int	fFirstSpecies;
	const int	fTemperature;
	const int	fVariablesWithoutSpecies;
	char		**fVariableNames;

	FILE		*ffpEI;

	Double		fTStart;
	Double		fTEnd;
	Double		fTCurr;
	Double		fDeltaT;
	Double		fDeltaTStart;
	Double		fTimeCut;
	int			fNOutputs;
	Double		fEquivalenceRatio;
	
	Double		fPressStart;
	Double		fPressEnd;
	Double		fTempLeftStart;
	Double		fTempRightStart;
	Double		fTempLeftEnd;
	Double		fTempRightEnd;
	Double		fRLeftStart;
	Double		fRRightStart;
	Double		fRLeftEnd;
	Double		fRRightEnd;
	Double		fPhiEnd;

	VectorPtr	fGridSol;
	MatrixPtr	fStartSol;
	Flag		fEquidistant;
	int			fNGridPoints;
	int			fMaxGridPoints;
	int			*fBCFlagLeft;
	int			*fBCFlagRight;
};

typedef TTrans1DIsoChor *TTrans1DIsoChorPtr;
