#include "Interrupt.h"

#ifdef HP
#include <stdio.h>
#include <signal.h>


int gExit;

static void InterruptHandler( int signo )
{
	switch ( signo ) {
		case SIGTERM:
		case SIGINT:
		case SIGABRT:
			fprintf( stderr, "# caught signal no. %d.\n", signo );
			gExit = kExit;
			break;
		case _SIGUSR1:
			fprintf( stderr, "# caught user signal1 no. %d.\n", signo );
			if ( signal( _SIGUSR1, InterruptHandler) == SIG_ERR )
				fprintf( stderr, "# Couldn't install interruptHandler.\n" );
			gExit = kNewGrid;
			break;
		case _SIGUSR2:
			fprintf( stderr, "# caught user signal2 no. %d.\n", signo );
			if ( signal( _SIGUSR2, InterruptHandler) == SIG_ERR )
				fprintf( stderr, "# Couldn't install interruptHandler.\n" );
			break;
						
	}
}


void InstallInterruptHP( void )
{
	gExit = kNoSignal;
    if ( signal( SIGINT, InterruptHandler) == SIG_ERR )
        fprintf( stderr, "# Couldn't install interruptHandler.\n" );
    if ( signal( SIGTERM, InterruptHandler) == SIG_ERR )
        fprintf( stderr, "# Couldn't install interruptHandler.\n" );
    if ( signal( SIGABRT, InterruptHandler) == SIG_ERR )
        fprintf( stderr, "# Couldn't install interruptHandler.\n" );
    if ( signal( _SIGUSR2, InterruptHandler) == SIG_ERR )
        fprintf( stderr, "# Couldn't install interruptHandler.\n" );
    if ( signal( _SIGUSR1, InterruptHandler) == SIG_ERR )
        fprintf( stderr, "# Couldn't install interruptHandler.\n" );
}
#endif
