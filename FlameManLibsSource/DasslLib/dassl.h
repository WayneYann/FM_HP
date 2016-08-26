#ifdef SUN
#define DDASSL ddassl_
#elif defined SUNGXX
#define DDASSL ddassl_
#elif defined DIGITAL
#define DDASSL ddassl_
#elif defined LINUXGXX
#define DDASSL ddassl_
#endif

typedef void (*ResFunc)( Double *T, Double *y, Double *yPrime, Double *delta
			, int *iRes, Double *rPar, int *iPar );
typedef void (*JacFunc)( Double *T, Double *y, Double *yPrime, Double *pd
			, Double *cj, Double *rPar, int *iPar );
			
extern "C" {
	void DDASSL( ResFunc func, int *fNOfEquations, Double *fTime, Double *y
				, Double *yPrime
				, Double *tout, int *fInfo, Double *fRTol, Double *fATol
				, int *fIdid
				, Double *fRWork, int *fLRW, int *fIWork, int *fLIW
				, Double *rPar, int *object, JacFunc funcer );
	void SDASSL( ResFunc func, int *fNOfEquations, Double *fTime, Double *y
				, Double *yPrime
				, Double *tout, int *fInfo, Double *fRTol, Double *fATol
				, int *fIdid
				, Double *fRWork, int *fLRW, int *fIWork, int *fLIW
				, Double *rPar, int *object, JacFunc funcer );
	void UDASSL( ResFunc func, int *fNOfEquations, Double *fTime, Double *y
				, Double *yPrime, Double *tout, int *fInfo, int *fIsen
				, Double *fRTol, Double *fATol, Double *fRTols, Double *fATols
				, int *fIdid
				, Double *fSWork, int *fLSW, int *fISWork, int *fLISW
				, Double *fRWork, int *fLRW, int *fIWork, int *fLIW
				, Double *rPar, int *object );
}
