/*  The functions in this file compute the broadening factors.
 *  The names of the functions start with "Fc" followed by
 *  the label of the corresponding reaction.
 */

#include"BroadeningFactors.h"

/*	Mechanism file: "nHeptane.allstarnew_oks.mech"	*/

BFFunction gBroadening[4] = { Fca34, Fc36, Fca51, Fca58 };

#ifndef MECHANISM
#define MECHANISM ""
#endif

void TReaction::CheckBroadeningFactors( const char *mechName )
{
	char	*name = new char[strlen( MECHANISM ) + 6];
	sprintf( name, "/%s.pre", MECHANISM );
	if ( strstr( mechName, name ) == NULL ) {
		for ( int i = 0; i < 4; ++i ) {
			gBroadening[i] = FcErr;
		}
	}
}

Double FcErr( Double /*T*/ )
{
	fprintf( stderr, "#error: wrong broadening factors (%s) linked to program\n", MECHANISM );
	exit( 2 );

	return 0;
}

Double Fca34( Double T )
{
#line 80 "nHeptane.allstarnew_oks.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc36( Double T )
{
#line 90 "nHeptane.allstarnew_oks.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fca51( Double T )
{
#line 120 "nHeptane.allstarnew_oks.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fca58( Double T )
{
#line 135 "nHeptane.allstarnew_oks.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

