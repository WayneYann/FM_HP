/*  The functions in this file compute the broadening factors.
 *  The names of the functions start with "Fc" followed by
 *  the label of the corresponding reaction.
 */

#include"BroadeningFactors.h"

/*	Mechanism file: "gri.300.mech"	*/

BFFunction gBroadening[29] = { Fc12, Fc50, Fc52, Fc54, Fc56, Fc57, Fc59, Fc63, Fc70, Fc71, Fc72, Fc74, Fc76, Fc83, Fc85, Fc95, Fc131, Fc140, Fc147, Fc158, Fc174, Fc185, Fc237, Fc241, Fc289, Fc304, Fc312, Fc318, Fc320 };

#ifndef MECHANISM
#define MECHANISM ""
#endif

void TReaction::CheckBroadeningFactors( const char *mechName )
{
	char	*name = new char[strlen( MECHANISM ) + 6];
	sprintf( name, "/%s.pre", MECHANISM );
	if ( strstr( mechName, name ) == NULL ) {
		for ( int i = 0; i < 29; ++i ) {
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

Double Fc12( Double T )
{
#line 20 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc50( Double T )
{
#line 62 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc52( Double T )
{
#line 68 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc54( Double T )
{
#line 74 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc56( Double T )
{
#line 80 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc57( Double T )
{
#line 85 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc59( Double T )
{
#line 91 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc63( Double T )
{
#line 99 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc70( Double T )
{
#line 110 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc71( Double T )
{
#line 115 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc72( Double T )
{
#line 120 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc74( Double T )
{
#line 126 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc76( Double T )
{
#line 132 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc83( Double T )
{
#line 143 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc85( Double T )
{
#line 149 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc95( Double T )
{
#line 163 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc131( Double T )
{
#line 203 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc140( Double T )
{
#line 216 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc147( Double T )
{
#line 227 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc158( Double T )
{
#line 242 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc174( Double T )
{
#line 262 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc185( Double T )
{
#line 276 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc237( Double T )
{
#line 331 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc241( Double T )
{
#line 339 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc289( Double T )
{
#line 391 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc304( Double T )
{
#line 410 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc312( Double T )
{
#line 422 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc318( Double T )
{
#line 432 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

Double Fc320( Double T )
{
#line 438 "gri.300.mech"

	fprintf( stderr, "#error: no broadening function specified\n" );
	exit( 2 );

	return T;
}

