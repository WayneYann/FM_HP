/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 *
 *	Mathematica.[ch]
 *				Functions related to Mathematica output.
 */

#ifndef __Mathematica__
#define __Mathematica__


void MathematicaSymbol( char dest[80], const char *src );
int MathematicaNumber( char *m, FPUType x );

void WriteMathematicaFile( void );


#endif	/* __Mathematica__ */
