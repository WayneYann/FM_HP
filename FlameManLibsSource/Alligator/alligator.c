/*
	alligator.c: simple allocation package
	
	
	History:
		05.04.99	abort removed from fatalerrormsg (hp+rs)
		24.01.94	general revision (jg)
		27.12.92	1st version (jg)
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alligator.h"


void fatalerrormsg( const char *s, const char *file, int line )
{
	fprintf( stderr, "\n###  Fatal error: %s.\n", s );
	fprintf( stderr, "File \"%s\"; Line %d\n", file, line );
/*	abort();*/
	exit( 2 );
}


void *allocate( unsigned int size, const char *file, int line )
{
	char	*ptr;
	
	if ( size == 0 ) return NULL;
	ptr = (char *)malloc( size );
	if ( ptr != NULL ) {
		memset( ptr, 0, size );
	} else {
		fprintf( stderr, "allocate size %d", size );
		fatalerrormsg( "Memory exhausted", file, line );
	}
	return ptr;
}


void alligator_delete( void *ptr )
{
	if ( !ptr ) return;
	
	free( ptr );
}
