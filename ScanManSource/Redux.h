#ifndef __Redux__
#define __Redux__


typedef struct ForwardReaction {
	int			nReacs;			/*	number of all reactions						*/
	int			nForward;		/*	number of forward (all except *b) reactions	*/
	ReactionPtr	*reac;			/*	pointers to forward reactions				*/
} ForwardReaction, *ForwardReactionPtr;


typedef struct Reactants {
	int			nSpecs;
	int			nReactants;
	SpeciesPtr	*spec;
} Reactants, *ReactantsPtr;


typedef struct Reduction {
	int			nReactions;
	int			nSpecies;
	int			nReducs;
	int			*reac;				/* index of reac in ForwardReaction struct	*/
	int			*spec;				/* index of spec in Reactant struct			*/
} Reduction, *ReductionPtr;


typedef struct GlobalReaction {
	int			nGlobalSpecies;		/* number of non-steady-state species		*/
	int			nGlobalReacs;		/* number of global reactions				*/
	int			*spec;				/* index of spec in Reactants struct		*/
	int			*reac;				/* index of reac in ForwardReaction struct	*/
} GlobalReaction, *GlobalReactionPtr;


void FixForwardBackwardsPtrs( void );
void ReduceMechanism( void );
void ReduxCleanup( void );
int SpeciesIsSteadyState( SpeciesPtr sItem );
int NumberOfSteadyStateSpecies( void );
int Column( ReactionPtr reac );
int Row( SpeciesPtr spec );
ReactionPtr GetForwardReaction( int k );

int IndependentEQ( MatrixPtr m );
MatrixPtr MatrixMult( const MatrixPtr mat1, const MatrixPtr mat2 );
MatrixPtr Transpose( const MatrixPtr m );
MatrixPtr Inverse( MatrixPtr m );

/*	Special stuff for the code generator that supports reduced mechanisms.		*/
typedef struct SteadyStateInfo {
	MatrixPtr	nu;
	MatrixPtr	glob;
	MatrixPtr	rate;
} SteadyStateInfo, *SteadyStateInfoPtr;

extern void RedMech( SteadyStateInfoPtr );


extern void FatalError( const char *errorString );
extern const FPUType	kTiny;

#endif	/* __Redux__ */
