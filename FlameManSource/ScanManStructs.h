#if defined (applec) || defined (powerc)
#include <CursorCtl.h>
#endif

typedef char Flag;
typedef char *String;

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define MaxSpeciesOfReaction 100

extern FILE *yyin;


typedef struct CommandLineOption {
	FILE	*inFile;
	char	*inName;
} CommandLineOption, *CommandLineOptionPtr;

typedef struct Header {
	char	version[32];
	int		maxLenOfString;
	Flag	globalMechanism;
} Header, *HeaderPtr;

typedef struct Counter {
	int	atoms;
	int	species;
	int	pahSpecies;
	int	sootSpecies;
	int	steadyStates;
	int	reactions;
	int	pahReactions;
	int	sootReactions;
	int	thirdBodies;
	int dimensions;
} Counter, *CounterPtr;

typedef struct Atoms {
	char			*name;
	int				number;
} Atoms, *AtomsPtr;

#define NOFCOEFFS 7

typedef struct Species {
	char			*name;
	IntVectorPtr	composition;
	int				number;
	int				isSteadyState;	/* signals whether this is a steady state species		*/
	int				aux;			/* auxiliary field needed for reduced mechanisms		*/
	Double			molarMass;
	Double			coeffHot[NOFCOEFFS];
	Double			coeffCold[NOFCOEFFS];
	Double			k_over_eps;			/* temperature factor for evaluating sigma (= k/eps);
									   input is eps/k  (computed) */
	Double			sigma;
	Double			muCoeff;			/* constant factor for evaluation of mu in kg/(m*s)
									   (= 2.6693e-6 * sqrt(M)/sigma)  (computed) */

	VectorPtr		dCoeff;			/* constant factor for evaluation of D in m^2/s 
									   = 0.0188292*sqrt(1/Mi+1/Mj)/(.5*(sigma_i+sigma_j))^2
									   (computed)	*/
	VectorPtr		omegaCoeff;		/* constant factor for evaluation of omega_D 
									   (= sqrt(k_over_eps_i * k_over_eps_j))  (computed) */
} Species, *SpeciesPtr;

typedef struct Reaction {
	char		*label;
	int 		id;	// (jg)
	int			numberOfSpecies; /*  without thirdBody  */
	int			speciesNumber[MaxSpeciesOfReaction];
	Double		speciesCoeff[MaxSpeciesOfReaction];
	Double		a;
	Double		n;
	Double		e;
	Flag		withLindemann;
	Double		aInf;
	Double		nInf;
	Double		eInf;
	int			lindemannNumber;
	Double		fca;
	Double		fcb;
	Double		fcc;
	Double		fcTa;
	Double		fcTb;
	Double		fcTc;
	Flag		withThirdBody;
	int			thirdBodyNumber;
	Flag		partialEquilibrium;
	int			orderOfReaction;	
	Flag		hasForward;
	Flag		hasBackward;
	struct Reaction	*forward;		/*	pointer to forward reaction							*/
	struct Reaction	*backward;		/*	pointer to backward reaction, or NULL				*/
	Flag		isDouble;
    Flag        rateFromEquilibrium;
} Reaction, *ReactionPtr;

typedef struct ThirdBody {
	char			*name;
	int				number;
	IntVectorPtr	speciesNumber;
	VectorPtr		speciesCoeff;
	IntVectorPtr	set;  /* 0 means not set, 1 means set */
} ThirdBody, *ThirdBodyPtr;

typedef struct Link {
	struct Link	*next;
	struct Link	*last;
	void		*item;
} Link, *LinkPtr;

typedef struct List {
	char	*listName;
	LinkPtr	firstItem;
	LinkPtr	lastItem;
	int		items;
} List, *ListPtr;

typedef struct Dimension {
	char			*name;
	Double			value;
	int				kg;
	int				kgParaCoeff;
	char			kgParameter[10];
	int				m;
	int				mParaCoeff;
	char			mParameter[10];
	int				s;
	int				sParaCoeff;
	char			sParameter[10];
	int				K;
	int				KParaCoeff;
	char			KParameter[10];
	int				mole;
	int				moleParaCoeff;
	char			moleParameter[10];
	Flag			orderOfReaction;
	Double			orderOfReacValue;
	Flag			tempExponent;
	Double			tempExpValue;
	ListPtr			parameterValues;
} Dimension, *DimensionPtr;
