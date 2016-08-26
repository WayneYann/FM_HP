#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "WSS.h"


void call_me( WSSPtr wss, int n )
{
	const int	kAllocSize = 50;
	char		*ptr = WSSPush( wss, kAllocSize );
	
	fprintf( stderr, "# n = %d: ptr = 0x%x\n", n, ptr );
	if ( ptr ) memset( ptr, 0, kAllocSize );
	if ( n ) {
		call_me( wss, n-1 );
	} else {
		WSSTrace( wss, 1 );
		fprintf( stderr, "\t%d free bytes in work space stack 0x%x\n", WSSFreeMem( wss ), wss );
	}
	if ( ptr ) WSSPop( wss );
}


int main( int argc, char *argv[] )
{
	char	*ptr = NULL;
	int		i;
	WSSPtr	wss1 = NewWSS( 1000 );
	WSSPtr	wss2 = NewWSS( 100 );
	
	if ( argc == 1 ) {
		fputs( "# using wss1\n", stderr );
		call_me( wss1, 10 );	/*	should work	*/
		fputs( "# using wss2\n", stderr );
		call_me( wss2, 10 );	/*	should fail	*/
	}
	
	fputs( "# tracing wss1\n", stderr );
	WSSTrace( wss1, 1 );
	fprintf( stderr, "# FreeMem wss1: %d bytes\n", WSSFreeMem( wss1 ) );
	fputs( "# tracing wss2\n", stderr );
	WSSTrace( wss2, 1 );
	fprintf( stderr,"# FreeMem wss2: %d\n", WSSFreeMem( wss2 ) );
	
	if ( argc == 1 ) {
		ptr = WSSPush( wss1, 50 );
		for ( i = 0; i < 10; ++i ) WSSPush( wss1, 50 );
		fputs( "# tracing wss1 (before remove)\n", stderr );
		WSSTrace( wss1, 1 );
		if ( WSSRemove( wss1, ptr ) != 0 ) FATAL( "WSSRemove() failed" );
		fputs( "# tracing wss1 (after remove)\n", stderr );
		WSSTrace( wss1, 1 );
		fprintf( stderr, "# FreeMem wss1: %d bytes\n", WSSFreeMem( wss1 ) );
	}

	DisposeWSS( wss1 );
	DisposeWSS( wss2 );
	
	return 0;
}