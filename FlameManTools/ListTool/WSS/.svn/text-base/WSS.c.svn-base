/*
	WSS.s
	
	work space stack package
	version 1.0§
	11. Nov. 93 (pt)
	
	
	# MPW Library
	Begin
		C WSS.c -r -mc68020 -b3 -s WSS
		C WSS.c -r -mc68020 -b3 -d qMemDebug -o WSS.dbg.c.o -s WSS
		Lib WSS.c.o -o "{myLibs}"WSS.lib
		Lib WSS.dbg.c.o -o "{myLibs}"WSS.dbg.lib
		duplicate -y WSS.h "{myIncludes}"
 		delete WSS.c.o WSS.dbg.c.o
	End ·· "{WorkSheet}"
*/


#include <stdio.h>
#include <stdlib.h>


#include "alligator.h"
#include "WSS.h"


#define kAlignment		8		/* align to 8 byte boundaries	*/
static const int		kAlignMask = kAlignment-1;

#define	ALIGN(n)		n = ( (n & ~kAlignMask) + \
				( (n & kAlignMask) ? kAlignment : 0 ) ) / sizeof(n)

static const char		*corrupted_stack = "corrupted work space stack";

enum {
	kWSSLast = -1,
	kWSSFile = -2,
	kWSSLine = -3
};


/*	
	NewWSS() allocates a work space stack structure sized <size> bytes
	aligned to <kAlignment> byte boundaries.
*/
WSSPtr NewWSS( int size )
{
	WSSPtr	wss = NEW( WSS );
	uint	new_size = size;
	
	ALIGN( new_size );

	/*	the following does the same as the ALIGN() macro
	*/
	/*
	if ( new_size & kAlignMask )
		new_size = (new_size & ~kAlignMask) + kAlignment;
	else 
		new_size = new_size & ~kAlignMask;
	new_size /= sizeof(address);
	*/

	wss->base = wss->cur = NEW2( new_size, address );
	
	wss->top = wss->base + new_size;
	
#ifdef qDebug
	fprintf( stderr, "-> allocated %d bytes for work space stack at 0x%x\n", 
		new_size*sizeof(address), wss );
#endif	
	return wss;
}


/*
	DisposeWSS() deallocates a work space structure previously allocated 
	with a call to NewWSS().
*/
void DisposeWSS( WSSPtr wss )
{
	if ( !wss ) return;
	
	DELETE( wss->base );
	DELETE( wss );
}


/*
	WSSPush() returns a block of <size> bytes aligned to <kAlignment> byte
	boundaries or NULL, if the block cannot be allocated.
*/
void *_WSSPush( WSSPtr wss, int size, const char *file, int line )
{
	uint		new_size = size + 3 * sizeof(address);	/*	add space for bookkeeping.	*/
	address		*save = wss->cur;
	
	ALIGN( new_size );
	
	if ( wss->cur + new_size >= wss->top ) {			/* check for sufficient space.	*/
		fprintf( stderr, "File \"%s\"; Line %d # workspace stack exhausted.\n", 
			file, line );
		return NULL;
	}
	/*	store last block, file and line on the stack.
	*/
	wss->cur += new_size;
	wss->cur[kWSSLast] = (address)save;
	wss->cur[kWSSFile] = (address)file;
	wss->cur[kWSSLine] = (address)line;
	
#ifdef qDebug
	fprintf( stderr, "-> pushed %d bytes on the stack\n", new_size*sizeof(address) );
	fprintf( stderr, "-> cur = 0x%x, last = 0x%x, base = 0x%x\n", 
		wss->cur, wss->cur[kWSSLast], wss->base );
#endif

	return save;
}


/*
	WSSPop() removes the most recent block from the work space stack.
*/
void WSSPop( WSSPtr wss )
{
	if ( wss->cur - wss->base > 0 ) {
		wss->cur = (address *)wss->cur[kWSSLast];
	}
}


/*
	WSSFreeMem() returns the free space on the work space stack. Note that
	a single block of that size cannot be allocated, since the package uses
	12 bytes from the work space stack for each block allocated for bookeeping.
*/
int WSSFreeMem( WSSPtr wss )
{
	int		free_mem = ( wss->top - wss->cur ) * sizeof(address);

	if ( free_mem < 0 ) FATAL( corrupted_stack );

	return free_mem;
}


/*
	WSSTrace() scans through the work space stack and checks for consistency.
	If <print> is not zero, all pointers currently used are printed to stderr.
*/
void WSSTrace( WSSPtr wss, int print )
{
	address		*ptr = wss->cur;
	
	while ( ptr > wss->base ) {
		if ( print ) {
			fprintf( stderr, "File \"%s\"; Line %d # block of %d bytes at 0x%x\n",
				ptr[kWSSFile], ptr[kWSSLine],
				(ptr - (address *)ptr[kWSSLast]) * sizeof(address), 
				ptr[kWSSLast] );
		}
		ptr = (address *)ptr[kWSSLast];
		if ( ptr >= wss->top )  FATAL( corrupted_stack );
	}
	if ( ptr != wss->base ) FATAL( corrupted_stack );
}


/*
	WSSRemove() removes the pointer <to_be_freed> and all pointers which
	have been allocated after it from the work space stack. If successful,
	zero is returned, a non zero value otherwise.
*/
int WSSRemove( WSSPtr wss, void *to_be_freed )
{
	address		*ptr = wss->cur;
	
	while ( ptr > wss->base ) {
		if ( (address *)ptr[kWSSLast] == to_be_freed ) {
			wss->cur = (address *)ptr[kWSSLast];
			return 0;
		}
		ptr = (address *)ptr[kWSSLast];
	}
	if ( ptr != wss->base ) FATAL( corrupted_stack );
	
	return -1;
}

