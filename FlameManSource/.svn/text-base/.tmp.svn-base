#include "FlameMaster.h"

void ComputeSteadyStates( Double */*k*/, Double */*c*/, Double */*M*/ )
{
	cerr << "#error: no need to use function 'ComputeSteadyStates'" << NEWL;
	exit( 2 );
}

#ifndef MECHANISM 
#define MECHANISM ""
#endif

int SteadyStatesFunc( const VectorPtr /*x*/, VectorPtr /*fVec*/, void */*object*/ )
{
	cerr << "#error: no need to use function 'SteadyStatesFunc'" << NEWL;
	exit( 2 );

	return 0;
}

void TReaction::CheckSteadyStatesMech( const char *mechName )
{
	if ( strncmp( MECHANISM, mechName, strlen( MECHANISM ) ) != 0 ) {
		cerr << "#error: program linked with mechanism " << MECHANISM 
					<< ", input mechanism is " << mechName << NEWL;
	}
}

void ComputeProductionRates( Double */*cdot*/, Double */*w*/, Double */*k*/
								, Double */*c*/, Double */*M*/, Double /*temp*/, Double /*pressure*/ )
{
	cerr << "#error: no need to use function 'ComputeProductionRates'" << NEWL;
	exit( 2 );
}
