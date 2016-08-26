#include "TTransFlameSolver.h"

//*************************************************************************************************************
class TTransFlamelet : public TTransFlameSolver {
public:
	TTransFlamelet( FirstInputPtr firstInp ) : fTemperature(0), 
			fFirstSpecies(1),
			fVariablesWithoutSpecies(fInputData->fVariablesWithoutSpecies),
			fNOfEquations(fInputData->GetCounter()->species+fVariablesWithoutSpecies - fInputData->GetCounter()->steadyStates),
			TTransFlameSolver( firstInp ) { InitTTransFlamelet();};	
	virtual ~TTransFlamelet( void );


	void				Solve( void );
	ConstStringArray	GetVariableNames( void ) { return fVariableNames; };

private:
  	void	InitTTransFlamelet( void );
	void	SetInitialBC( void );
	void	SetInitialValues( void );
	Double	GetTempEnd( Double time, Double timeStart, Double timeCut
									, Double tempStart, Double tempEnd );
	void	ReadStartProfiles( TInputDataPtr inp );
	void	SetInitialValues( TInputDataPtr inp, StartProfilePtr sp );
	Double	GetScalarDiss( Double theTime );
	void	ReadPDF( void );
	void	WriteScalars( Double time, int i );

	const int	fFirstSpecies;
	const int	fTemperature;
	const int	fVariablesWithoutSpecies;
	const int	fNOfEquations;
	char		**fVariableNames;

	int fDeltaStepsOut;
	int			fPDFGridPoints;
	VectorPtr	fTimePDFIn;
	MatrixPtr	fPDF;
	FILE		*ffpEI;


	Double		fRPM;
	int			fVarsIn;
	VectorPtr	fTimeIn;
	VectorPtr	fTFuelIn;
	VectorPtr	fTOxIn;
	VectorPtr	fPressureIn;
	VectorPtr	fChiIn;
	VectorPtr	fZRIn;
	VectorPtr	fZMeanIn;
	VectorPtr	fZVarIn;
	VectorPtr	fXOverD;
	Flag		fUseInput;
	Double		fTStart;
	Double		fTEnd;
	Double		fTCurr;
	Double		fDeltaT;
	Double		fDeltaTStart;
	Double		fTimeCut;
	int			fNOutputs;
	Double		fScalarDissRate;
	
	Double		fPressStart;
	Double		fTempOxStart;
	Double		fTempFuelStart;
	Double		fPressEnd;
	Double		fChiEnd;
	Double		fTempOxEnd;
	Double		fTempFuelEnd;

	VectorPtr	fGridSol;
	MatrixPtr	fStartSol;
	Flag		fEquidistant;
	int			fNGridPoints;
	int			fMaxGridPoints;
	int			*fBCFlagLeft;
	int			*fBCFlagRight;
};

typedef TTransFlamelet *TTransFlameletPtr;
