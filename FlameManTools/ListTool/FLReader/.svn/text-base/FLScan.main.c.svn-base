/*
	FLScan.main.c: driver for the flamelet file reader.
	
	© Josef Gšttgens, Peter Terhoeven, Ina Herwono, 1993
	
	History:
		11/3/93		InitFLReader() changed to enable parsing
					more than a single file.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "alligator.h"
#include "ArrayManager.h"

#ifdef applec
#include <CursorCtl.h>
#endif

#include "Stack.h"
#include "regex.h"
#include "List.h"

#include "FLReader.h"

#include "FLScan.h"
#include "FLScan.tab.h"
#include "FLPatch.h"

#undef	qDebugScanner

#define MaxElems	1024

/* regex */

#define BYTEWIDTH 8

/* end of regex */


#define	NEW_PRINT_HEADER

/*
	Globals
*/
FILE		*gBison = NULL;
FILE 		*gFlex = NULL;
FILE		*gFloatBuffer = NULL;
static char	glookup[256];
int 		gSn = 0;

static char			gFormat[256];
static int			gLastDepth = 0;
static const char	*gLastContext = ":";

#ifdef qDebugScanner
YYSTYPE	fllval;
#endif

int				gFBMaxElems,gFBElems;
Double			*gFBuf = NULL, *gFBufhead = NULL;
ListPtr			gFLRSymbolList = NULL, gFLRNumberList = NULL, gFLRIntList = NULL, 
				gFLRFloatList = NULL, gFLRStringList = NULL, gFLRArrayList = NULL;
static ListPtr	gFLRHeaderList = NULL;

struct re_pattern_buffer	gRe_buf = {0};


void ResetFloatBuffer( void )
{
	gFBuf = gFBufhead;
	gFBElems = 0;
}


/* FLRCompare : compare two string with case insensitive 	*/
/* return value = 0 : the strings are the same			*/
/* return value != 0 : different strings				*/

int FLRCompare( const char *ptr1, const char *ptr2 ) 
{
	
	for ( ; glookup[*ptr1] == glookup[*ptr2]; ptr1++, ptr2++ )
		if ( *ptr1 == '\0' ) return 0;
		
	return *ptr1 - *ptr2;
}


/* RFindItem : find the item at position pos in a list with the FLRCompare function */

static void *RFindItem( ListPtr list, int *pos, void *target, 
	int (*cmp)( void *, void * ) )
{
	int		i;
	LinkPtr temp = list->fHead;
	
	if ( *pos > 0 ){
		for ( i = *pos; i > 0 ; --i ) {
			if ( !temp ) return NULL;
			temp = temp->fNext;
		}
	}
	else { 
		*pos = 0;
	}
	
	while ( temp ) {
		++*pos;			/* the next starting position... */
		if ( (*cmp)( temp->fItem, target ) ) {
			return temp->fItem;
		}
		temp = temp->fNext;
	}

	return NULL;
}

/* SortItems : sort the items in the symbollist with sort-function category */

static void SortItems( ListPtr list, void *aux, int (*sort)(void *, void *) )
{
#	ifdef applec
#	pragma unused(aux)
#	endif
	LinkPtr temp = NULL, hPtr = NULL;
	int		changes;
	
	do {
		changes = 0;
		for ( temp = list->fHead; temp && (temp->fNext != NULL); temp = temp->fNext ) {
			if( ((*sort)(temp->fItem,temp->fNext->fItem)) > 0 ){
				hPtr = temp->fItem; 
				temp->fItem = temp->fNext->fItem;
				temp->fNext->fItem = hPtr;
				++changes;
			}
		}
	} while ( changes > 0 );
}

/* SortContext : sort-function for SortItems */

static int SortContext( void *ptr1, void *ptr2 )
{
	FLRSymbolPtr	s = (FLRSymbolPtr)ptr1;
	FLRSymbolPtr	t = (FLRSymbolPtr)ptr2;

	return( FLRCompare( s->context,t->context ) );
}

/* LookupTable: make a lookup table for case insensitive comparing */

static void LookupTable( void ) 
{
	int		i;
	
	for (i = 0; i < 256 ; ++i) {
		glookup[i] = i;
	}
	
	for (i = 65; i < 91 ; ++i) {
		glookup[i] = glookup[i+32];
	}
	
	glookup[128] = glookup[138];
	glookup[129] = glookup[140];
	glookup[130] = glookup[141];
	glookup[131] = glookup[142];
	glookup[132] = glookup[150];
	glookup[133] = glookup[154];
	glookup[134] = glookup[159];
	glookup[203] = glookup[136];
	glookup[204] = glookup[139];
	glookup[205] = glookup[155];
}

/* PrintSymbol : a "ForAll" function to print the list to file pointer fp */

void PrintSymbol( void *ptr, void *aux )
{   
	AMPrintOptions	prnt;
	FILE			*fp = (FILE *)aux;
	FLRSymbolPtr	s = (FLRSymbolPtr)ptr;
 
	if ( !fp ) fp = stdout;
	if (s->type != kFLRArray) fprintf( fp, "\nId = \"%s\"\n  Value = ", s->id );

	switch(s->type) {
		case kFLRInteger:		fprintf(fp,"%d\n",s->val.i); break;
		case kFLRFloat:			fprintf(fp,"%g\n",s->val.f); break;
		case kFLRBoolean:		fputs( (s->val.bool) ? "TRUE\n" : "FALSE\n", fp ); break;
		case kFLRString	:		fprintf( fp, "\"%s\"\n", s->val.s ); break;
		case kFLRArray:			DefaultAMPOpts( &prnt );
								prnt.lineWidth = 80;
								prnt.format = "%12.5g";
								prnt.sep = " ";
								prnt.title =  s->id;
								prnt.rowLabel = prnt.colLabel = NULL;
								PrintVector( s->val.v, &prnt, fp );
								break;
		default			:	/* do nothing */
								break;
	}
	
	fprintf( fp, "  Unit = \"%s\"\n", (s->unit) ? s->unit : "<none>" );
	fprintf( fp, "  Context = \"%s\"\n", s->context );
	fprintf( fp, "  Detached: %s", (s->detached) ? "TRUE" : "FALSE" );
}

/* PrintSymbolList : a special print routine to print out the symbollist */
/*					 to file pointer fp. 								*/

void PrintSymbolList( FILE *fp ) 
{	
	if ( !fp ) fp = stdout;
	fputs( "Symbols:\n", fp );
	ForAll( gFLRSymbolList, PrintSymbol, fp, kListHead );

}

/* DetachItem : function to set the value detached = TRUE of symbol s */
/* aux	: unused */
/* DetachItem could be used as "ForAll" function */

void DetachItem( void *ptr, void *aux )
{
#	ifdef applec
#	pragma unused(aux)
#	endif
	((FLRSymbolPtr)ptr)->detached = TRUE;
}


/* ContextDepth : function to count how deep the context of {all} is. */

static int ContextDepth( const char *all ) 
{
	int		n = 0;
	
	++all;			/* skip first character */
	while( *all ) {
		if ( *all == ':' ) ++n;
		++all;
	}
	
	return n;
}


/* PrintValue : routine used in WriteFlameletFile to print the value  	*/
/*				in flamelet format.										*/

static void PrintValue( FLRSymbolPtr s, FILE *fp )
{
	switch ( s->type ) {
		case kFLRInteger	:	if( s->unit ) {
									fprintf(fp,"%s%s = %d [%s]\n",gFormat,s->id,
											s->val.i,s->unit);
								} else {
									fprintf(fp,"%s%s = %d\n",gFormat,s->id, s->val.i);
								} break;
		case kFLRFloat		:	if ( s->unit ) {
									fprintf(fp,"%s%s = %g [%s]\n",gFormat,s->id,
											s->val.f,s->unit);
								} else {
									fprintf(fp,"%s%s = %g\n",gFormat,s->id,s->val.f);
								} break;
		case kFLRBoolean	:	if ( s->val.bool ) {
									fprintf(fp,"%s%s = TRUE\n",gFormat,s->id);
								} else {
									fprintf(fp,"%s%s = FALSE\n",gFormat,s->id);
								} break;
		case kFLRString		:	fprintf(fp,"%s%s = \"%s\"\n",gFormat,s->id,s->val.s);
								break;
		default				:	break;
	}
}


/* GetContext: function to get the n-th context from a string 			*/
/*				the n-th context will be copied to parameter context	*/

static void GetContext( const char *all, char *context, int n )
{
	char	*buffer = NEW2( strlen(all)+1, char );
	int		i = 1, j = 0, k = 0;
	
	while( all[i] != '\0' ) {
		if( all[i] != ':' ) {
			buffer[k] = all[i];
			++k;
		}
		else {
			++j;
			if (j == n) {		/* found the n'th context */
				buffer[k] = '\0';
				strcpy(context,buffer);
				break;
			}
			else {
				k = 0;
			}
		}
		++i;
	} 
	DELETE( buffer );
}

/* PrintFLHeader : print out the header section of flamelet file 	*/
/*					to file pointer fp							 	*/
/* PrintFLHeader will be called from WriteFlameletFile				*/

#ifndef NEW_PRINT_HEADER
void PrintFLHeader( void *ptr, void *aux )
{
	FLRSymbolPtr	s = (FLRSymbolPtr)ptr;
	int 			n, i = 0, j = 0;	/* gSn : Number of elements in Stack (Global)*/
	int 			counter = 0, done = FALSE, k;
	char 			buffer[80], *contextbuf[10];
	FILE			*fp = (FILE *)aux;
	
#	ifdef applec
	SpinCursor( 1 );
#	endif
	
	if ( !fp ) fp = stdout;
	
	if ( s->type != kFLRArray ) {
		n = ContextDepth( s->context );
		if ( n == 0) {		/* item with global context */
			strcpy( gFormat, "" );
			PrintValue( s, fp );			
		}
		else {				/* item with local context, possibly nested */
			if ( FLRStackEmpty() ) {	/* first item with local context */
				for ( i = 1; i <= n; ++i ) {
					GetContext( s->context, buffer, i );
					fprintf( fp, "\n%s\nBegin\n", buffer );
					strcpy( gFormat, "\t" );
					FLRStackPush( buffer );
					++gSn;
				}
				PrintValue( s, fp );
			}
			else { 			/* stack is not empty, we already printed the context name */
				if( n == gSn ) {		/* this item has the same context depth */
					for ( i = gSn; i > 0; --i ) contextbuf[i] = FLRStackPop();
					for ( i = 1; i <= n; ++i) {
						GetContext( s->context, buffer, i );
						if ( FLRCompare( buffer, contextbuf[i] ) == 0 ) ++counter;
					}	
					if ( counter == n )  { /* and the same context */
						PrintValue( s, fp );
						for ( i = 1; i <= gSn; ++i ) {
							FLRStackPush( contextbuf[i] );
							DELETE( contextbuf[i] );
						}	
					}
					else {				/* item has another context */
						done = FALSE;
						for ( i = 1; i <= n; ++i ) {
							GetContext( s->context, buffer, i ) ;
							if( FLRCompare( buffer, contextbuf[i] ) != 0 ) {
								for ( j = 0; j < (gSn + 1 - i); ++j ) {
									strcpy( gFormat,"" );
									for( k = 0; k < (i-1-j); ++k ) {
										strcat( gFormat,"\t" );
									}
									fprintf ( fp, "%sEnd\n", gFormat );
									if ( i == 1 ) fprintf( fp, "\n" );
									GetContext( s->context, buffer, i+j );
									DELETE( contextbuf[i+j] );
									contextbuf[i+j] = NEW2( strlen(buffer)+1, char );
									strcpy( contextbuf[i+j], buffer );
								}
								for ( j = i; j < ( gSn+1); ++j ) {
									GetContext( s->context, buffer, j );
									strcpy( gFormat, "" );
									for( k = 0; k < (j-1); ++k ) {
										strcat( gFormat, "\t" );
									}
									fprintf( fp, "%s%s\n%sBegin\n", gFormat, buffer,
											gFormat );
								}
								done = TRUE;
							}
							if ( done ) break;
						}
						strcat( gFormat, "\t" );
						PrintValue( s, fp );
						for( i = 1;i < gSn+1; ++i ) {
							FLRStackPush( contextbuf[i] );
							DELETE( contextbuf[i] );
						}	

					}
				} /* end of (n == gSn) */				
				else if ( n > gSn ) {		/* item has a higher context depth */
					for ( i = gSn; i > 0; --i ) contextbuf[i] = FLRStackPop();
					done = FALSE;
					for ( i = 1; i < (gSn+1); ++i ) {
						GetContext( s->context, buffer, i );
						if( FLRCompare( buffer, contextbuf[i] ) != 0 ) {
							for ( j = 0; j< (gSn + 1 - i); ++j ) {
								strcpy( gFormat, "" );
								for ( k = 0; k < (i-1-j); ++k ) {
									strcat( gFormat, "\t" );
								}
								fprintf( fp, "%sEnd\n", gFormat );
								if (i == 1) fprintf( fp, "\n" );
								GetContext( s->context, buffer, i+j );
								DELETE( contextbuf[i+j] );
								contextbuf[i+j] = NEW2( strlen(buffer)+1, char );
								strcpy( contextbuf[i+j], buffer );
							}
							for( j = i; j < ( gSn+1); ++j ) {
								GetContext( s->context, buffer, j );
								strcpy( gFormat, "" );
								for( k = 0; k < (j-1); ++k ) {
									strcat( gFormat, "\t" );
								}
								fprintf( fp, "%s%s\n%sBegin\n", gFormat, buffer, 
											gFormat );
							}
						done = TRUE;
						}
						if ( done ) break;
		
					}
					for ( i = (gSn+1); i < (n+1); ++i ) {
						GetContext( s->context, buffer, i );
						fprintf( fp, "%s%s\n%sBegin\n", gFormat, buffer, gFormat);
						strcat( gFormat, "\t" );
						contextbuf[i] = NEW2( strlen(buffer)+1, char );
						strcpy( contextbuf[i], buffer );
						++gSn;
					}
					PrintValue( s, fp );
					for ( i = 1; i < (n+1); ++i ) {
						FLRStackPush( contextbuf[i] );
						DELETE( contextbuf[i] );
					}
				}	/* end of (n > gSn) */
				
				else { 				/* n < gSn */
					for ( i = 0;i < gSn; ++i) {
						
						strcpy(gFormat,"");
						for ( k = 0; k < (gSn-i-1); ++k ) {
							strcat( gFormat, "\t" );
						}
						fprintf( fp, "%sEnd\n", gFormat );
						contextbuf[i] = FLRStackPop();
						DELETE( contextbuf[i] );
					}
					gSn = 0;
					fprintf( fp, "\n" );
					for( i = 1; i < (n+1); ++i ) {
						GetContext( s->context, buffer, i );
						strcpy( gFormat, "" );
						for ( k = 0; k < (i-1); ++k ) {
							strcat( gFormat, "\t" );
						}
						fprintf( fp, "%s\nBegin\n", buffer );
						FLRStackPush( buffer );
						++gSn;
					}
					strcat( gFormat, "\t" );
					PrintValue( s, fp );
				} /* end of (n < gSn) */
			}
			
			}
		}
}

#else

static void OpenContext( const char *contextName, int curDepth, FILE *fp )
{
	int		i, n = curDepth-1;

	for ( i = 0; i < n; ++i ) fputc( '\t', fp );
	fprintf( fp, "%s Begin\n", contextName );
/*	for ( i = 0; i < n; ++i ) fputc( '\t', fp );
	fputs( "Begin\n", fp );*/
}


static void CloseContext( int curDepth, FILE *fp  )
{
	int		i, n = curDepth;
	
	for ( i = 0; i < n; ++i ) fputc( '\t', fp );
	fputs( "End\n", fp );
}


static int CommonDepth( const char *ptr1, const char *ptr2 )
{
	int		depth;
	
	++ptr1, ++ptr2;		/* skip first char (always ':') */
	for ( depth = 0; glookup[*ptr1] == glookup[*ptr2]; ++ptr1, ++ptr2 ) {
		if ( *ptr1 == ':' ) ++depth;
		if ( *ptr1 == '\0' ) break;
	}	
	
	return depth;
}


void PrintFLHeader( void *ptr, void *aux )
{
	FLRSymbolPtr	sym = (FLRSymbolPtr)ptr;
	FILE			*fp = (FILE *)aux;
	int				i, con;
	char			*buffer = NULL;
	const char		*curContext = sym->context;
	const int		kCurDepth = ContextDepth( curContext );
	const int		kCommonDepth = CommonDepth( gLastContext, curContext );
	
	if ( sym->type == kFLRArray ) return;		/* no arrays in header */
	
	buffer = NEW2( strlen( curContext )+1, char );

	con = gLastDepth;
	
	while ( con > kCommonDepth ) {
		--con;
		CloseContext( con, fp );
	}

	con = kCommonDepth;
	while ( con++ < kCurDepth ) {			/* open new context */
		GetContext( curContext, buffer, con );
		OpenContext( buffer, con, fp );
	}

	for ( i = 0; i < kCurDepth; ++i ) gFormat[i] = '\t';
	gFormat[kCurDepth] = '\0';
	PrintValue( sym, fp );
	
	gLastDepth = kCurDepth;
	gLastContext = curContext;
	
	DELETE( buffer );
}

#endif	/* NEW_PRINT_HEADER */

/* PrintFLBody : print out the body section of flamelet file 	*/
/*					to file pointer fp 							*/
/* PrintFLBody will be called from WriteFlameletFile			*/

void PrintFLBody( void *ptr, void *aux )
{	
	FLRSymbolPtr	s = (FLRSymbolPtr)ptr;
	FILE			*fp = (FILE *)aux;
	Double 			*value = NULL;
	int				len, spalte = 1, i;

#	ifdef applec
	SpinCursor( 1 );
#	endif

	if (!fp) fp = stdout;
	if( s->type = kFLRArray) {
		if ( s->unit )	fprintf(fp,"\n%s [%s]\n",s->id,s->unit);
		else			fprintf(fp,"\n%s\n",s->id);

		value = s->val.v->vec;
		len = s->val.v->len;
		for ( i = 0; i< len; ++i) {
			if (spalte != 5) {
				fprintf(fp,"\t%e",*value);
				++spalte;
			}
			else {
				fprintf(fp,"\t%e\n",*value);
				spalte = 1;
			}
			++value;
		}
	}
}

/* InitFloatBuffers : initialize the floating point buffer 		*/

void InitFloatBuffers( int n ) 
{

	gFBuf = New1DArray(1024);
	gFBufhead = gFBuf;
	gFBMaxElems = n;
	gFBElems = 0;
}

/*	AddNumber: function to add value f to floating point buffer	*/

void AddNumber(Double f)
{
	if ( gFBElems == 0 ) {
		gExpectArray = TRUE;
	}
	if ( gFBElems < MaxElems ) {
		*gFBuf++ = f;
		++gFBElems;
	} else {
		FATAL( "Buffer for floating point numbers too small" );
	}
}

/* 	PrintFloatBuffer: print out the values to a file pointer 	*/
/*						named gFloatBuffer						*/

void PrintFloatBuffer(void)
{
	int		i,j=0, n = gFBElems;
	
	gFBuf = gFBufhead;
	
	for(i=0;i<n;++i,gFBuf++) {
		if(j<4) {
			fprintf(gFloatBuffer,"\t%lf",*gFBuf);
			++j;
		}
		else {
			fprintf(gFloatBuffer,"\t%lf\n",*gFBuf);
			j = 0;
		}
	}
	fputs( "\n", gFloatBuffer );
}



/* PutSymH : put the Symbol in header section in a list (symbollist) 		*/
/* id	: identifier														*/
/* unit : the unit, unit = NULL if the unit doesn't exist					*/
/* type : type of the value : kFLRString = string	(vals)					*/
/*							  kFLRInteger = integer	(vali)					*/
/*							  kFLRFloat   = float	(valf)					*/
/*							  kFLRBoolean = boolean	(vali)					*/
/* Example:																	*/
/* with type = kFLRString , the other value (vali,valf) should be set to NULL	*/

void PutSymH ( char *id,char *unit,int type,int vali,Double valf,char *vals ) 

{
	char			*str1 = NULL, *buffer[10];
	int				i = 0, j, leng;
	FLRSymbolPtr	item = NULL;
	
	item = NEW( FLRSymbol );
	item->type = type;
	item->detached = FALSE;
	switch (type) {
		case kFLRBoolean	: item->val.bool = vali;break;
		case kFLRInteger	: item->val.i    = vali;break;
		case kFLRFloat		: item->val.f    = valf;break;
		case kFLRString		: item->val.s    = vals;break;
		default 			: break;
	}
	
	item->unit = unit;
	item->id = id;

	i = 0;
	if( !FLRStackEmpty() ) {
		for (i = 0 ; i < 10 ; ++i ) {
			buffer[i] = NULL;
		}
		FLRStackReverse();
		i = 0; leng = 0;
		str1 = NEW2(160,char);
		strcat(str1,":");
		while( !FLRStackEmpty() ) {
			buffer[i] = FLRStackPop();
			leng += strlen(buffer[i]) + 1;
			strcat(str1,buffer[i]);
			strcat(str1,":");
			++i;
		}
		item->context = NEW2(leng+i+1,char);
		strcpy(item->context,str1);
		for( j = 0 ; j < i ; ++j ) {
			FLRStackPush( buffer[j] );
			DELETE( buffer[j] );
		}
		DELETE(str1);
	}
	else {
		item->context = NEW2(2,char);
		strcpy(item->context,":");
	}
	AddItem( gFLRSymbolList, item, kListTail );
	AddItem( gFLRHeaderList, item, kListTail );
	switch (type) {
		case kFLRBoolean	: AddItem( gFLRIntList, item, kListTail);break;
		case kFLRInteger	: AddItem( gFLRIntList, item, kListTail);AddItem( gFLRNumberList, item, kListTail);break;
		case kFLRFloat		: AddItem( gFLRFloatList, item, kListTail);AddItem( gFLRNumberList, item, kListTail);break;
		case kFLRString		: AddItem( gFLRStringList, item, kListTail);break;
		default				: break;
	}
}

/* PutSymBody : put the Symbol in body section in symbollist 				*/
/* vector : Double type pointer point to first address from the 1-D Double 	*/
/*			array 															*/
/* len	  : Number of values in the array 									*/

void PutSymBody(char *id,char *unit,Double *vector,int len) 
{
	int				i;
	FLRSymbolPtr	item = NULL;
	
	item = NEW( FLRSymbol );
	item->id = id;
	item->context = NEW2( 2, char );
	item->type = kFLRArray;
	item->detached = FALSE;
	item->unit = unit;
	strcpy( item->context, ":" );
	
	item->val.v = NewVector( len );        /* initialize a new vector */
	for (i = 0 ; i<len ; ++i ) {
		item->val.v->vec[i] = vector[i]; /* copy the contents of the vector to item */
	}
	AddItem( gFLRSymbolList, item , kListTail);
	AddItem( gFLRArrayList, item , kListTail);
}


/* lowstring : set all characters in string to lowercase 	*/
/* lowstring is not used in the program at the moment		*/

static char *lowstring( char *s ) 
{
	char 	*s2 = NULL;
	int 	i, n = strlen( s );
	int		stop = 0;
	
	s2 = NEW2( n+1, char );
	for( i = 0; i < n; ++i ) {
		if ( s[i] == '-' ) stop = 1; 
		if ( stop == 0)	s2[i] = tolower( s[i] );
		else s2[i] = s[i];
	}
	s2[n] = '\0';
	strcpy(s,s2);

	DELETE(s2);
	return s;
}


/* gleichID : FLRCompare-function to FLRCompare the IDs		*/

static int gleichID( void *ptr, void *aux )
{
	FLRSymbolPtr	s = (FLRSymbolPtr)ptr;
	char			*ID = (char *)aux;

	return( ( FLRCompare( s->id, ID ) == 0 ) );
}

/* gleichContext: FLRCompare-function to FLRCompare the contexts */

static int gleichContext( void *ptr, void *aux )
{
	FLRSymbolPtr	s = (FLRSymbolPtr)ptr;
	char			*CON = (char *)aux;


	return( (FLRCompare(s->context,CON) == 0) );
}


/* getSymbol : function to return the symbol with ID */

FLRSymbolPtr getSymbol( char *id ) 
{
	return ( RFindItem( gFLRSymbolList, NULL, lowstring(id), gleichID ) );	
}

/* HasContext : routine to check if {all} has a context   	*/
/* return value = 1 : context exists						*/

static int HasContext( const char *all )
{
	while ( *all ) {
		if ( *all++ == ':' ) return TRUE;
	}
	
	return FALSE;
}

/* GetIDandContext wird nur verwendet nach dem Aufruf von HasContext 	*/
/* GetIDandContext : get the id and the context from {all}				*/

static void GetIDandContext( const char *all, char *id, char *context )

{
	int i, j = 0, k = 0, len;
	
		
	len = strlen(all);
	i = len-1;
	while( all[i] != ':' ) {
		--i;
		++k;
	}
	j = 0;
	++i;
	
	while(i<len) {
		*id = all[i];
		++id;
		++i;
		++j;
	}
	strcat(id,"\0");
	if (all[0] != ':')	strcpy(context,":");
	strncat(context,all,len - k);
	strcat(context,"\0");
}

/* GetString : get a string value from {all} 	*/
/* val : the value which found					*/
/* all : id with(out) context					*/
/* return value = 1: the value found			*/

int GetString( const char *all, char *val ) 
{
	char			*id = NULL, *context = NULL;
	FLRSymbolPtr	item;
	int				pos = 0;
	
	if ( gFLRStringList->fNumItems == 0 ) { 
		fprintf(stderr,"# GetString: item \"%s\" not found (empty stringlist).\n",all);
		return FALSE;
	}
	
	id = NEW2( strlen(all) + 2, char );
	
	if ( HasContext( all ) ) 	{
		context = NEW2( strlen(all)+2, char );
		GetIDandContext( all, id, context );
	}
	else { 
		context = NEW2( 2, char );
		strcpy( context, ":" );
		strcpy( id, all );
	}
	
	item = RFindItem( gFLRStringList, &pos, id, gleichID );
	
	while ( item ) {
		if ( ( FLRCompare(context,item->context ) == 0 ) ) { 
			DELETE(id);
			DELETE(context);
			strcpy(val,item->val.s);		/* return a reference here */
			return TRUE;
		}
		else { 
			item = RFindItem(gFLRStringList, &pos, id, gleichID );
		}	
	}
#	ifdef qDebug
	fprintf(stderr,"# GetString: item \"%s\" not found.\n",id);
#	endif

	/* Clean Up  */
	DELETE(id);
	DELETE(context);

	return FALSE;
}

/*GetNumber : get the number value and unit from {all} 	*/
/* all : id with(out) context							*/
/* return value = 1: the value found					*/

int GetNumber( const char *all, Double *val, char *unit )
{
	char 			*id = NULL,*context = NULL;
	FLRSymbolPtr	item = NULL;
	int  			pos = 0;
	
	if ( gFLRNumberList->fNumItems == 0 ) { 
		fprintf(stderr,"# GetNumber: item \"%s\" not found (empty numberlist).\n",all);
		return FALSE;
	}

	id = NEW2( strlen(all) + 2, char );
	if ( HasContext( all ) ) 	{
		context = NEW2( strlen(all)+2, char );
		GetIDandContext( all, id, context );
	}
	else { 
		context = NEW2( 2, char );
		strcpy( context,":" );
		strcpy( id, all );
	}

	item = RFindItem(gFLRNumberList, &pos, id, gleichID );
	while (item) {
		if ( ( FLRCompare( context,item->context ) == 0 ) ) {
			if ( item->type == kFLRInteger ) {
				DELETE( id );
				DELETE( context );
				if ( unit ) {
					if ( item->unit ) strcpy( unit, item->unit );		/* return a reference here */
					else unit[0] = '\0';
				}
				*val = (Double)item->val.i;
				return TRUE;
			}
			if ( item->type == kFLRFloat )  {
				DELETE( id );
				DELETE( context );
				if ( unit ) {
					if ( item->unit ) strcpy( unit,item->unit );		/* return a reference here */
					else unit[0] = '\0';
				}
				*val = item->val.f;
				return TRUE;
			}
		}
		else { 
			item = RFindItem(gFLRNumberList, &pos, id, gleichID );
		}	
	}

#	ifdef qDebug
	fprintf(stderr,"# GetNumber: item \"%s\" not found.\n",id);
#	endif
	/* Clean Up  */

	DELETE(id);
	DELETE(context);
	
	return FALSE;
}

/*GetInteger : get the integer value and the unitfrom {all}	*/
/* all : id with(out) context								*/
/* return value = 1: the value found						*/

int GetInteger( const char *all, int *val, char *unit )
{
	char			*id = NULL, *context = NULL;
	FLRSymbolPtr	item = NULL;
	int				pos = 0;

	if (gFLRIntList->fNumItems == 0) { 
		fprintf(stderr,"# GetInteger: item \"%s\" not found (empty integerlist).\n",all);
		return FALSE;
	}
	
	id = NEW2(strlen(all) +2, char);
	if ( HasContext(all)) 	{
		context = NEW2(strlen(all)+2, char);
		GetIDandContext(all,id,context);
	}
	else { 
		context = NEW2(2, char);
		strcpy(context,":");
		strcpy(id,all);
	}
	item = RFindItem(gFLRIntList, &pos, id, gleichID );
	while (item) {
		if ( (FLRCompare(context,item->context) == 0) ) {
			DELETE(id);
			DELETE(context);
			if ( unit ) {
				if ( item->unit ) strcpy( unit, item->unit );		/* return a reference here */
				else unit[0] = '\0';
			}
			*val = item->val.i;
			return TRUE;
		}
		else { 
			item = RFindItem(gFLRIntList, &pos, id, gleichID );
		}	
	}

#	ifdef qDebug
	fprintf(stderr,"# GetInteger: item \"%s\" not found.\n",id);
#	endif

	DELETE(id);
	DELETE(context);

	return FALSE;
}

/*GetFloat : get the Double value and the unit	from {all}	*/
/* all : id with(out) context								*/
/* return value = 1: the value found						*/

int GetFloat(const char *all,Double *val,char *unit)
{
	char			*id = NULL, *context = NULL;
	FLRSymbolPtr	item = NULL;
	int				pos = 0;
	
	if (gFLRFloatList->fNumItems == 0) { 
		fprintf(stderr,"# GetFloat: item \"%s\" not found (empty floatlist).\n",all);
		return FALSE;
	}

	id = NEW2(strlen(all) +2, char);
	if ( HasContext(all)) 	{
		context = NEW2(strlen(all)+2, char);
		GetIDandContext(all,id,context);
	}
	else { 
		context = NEW2(2, char);
		strcpy(context,":");
		strcpy(id,all);
	}
	item = RFindItem(gFLRFloatList, &pos, id, gleichID );
	while (item) {
		if ( (FLRCompare(context,item->context) == 0) ) {
			DELETE(id);
			DELETE(context);
			if ( unit ) {
				if ( item->unit ) strcpy( unit, item->unit );		/* return a reference here */
				else unit[0] = '\0';
			}
			*val = item->val.f;
			return TRUE;
		}
		else { 
			item = RFindItem(gFLRFloatList, &pos, id, gleichID );
		}	
	}
#	ifdef qDebug
	fprintf(stderr,"# GetFloat: item \"%s\" not found.\n",id);
#	endif

	DELETE(id);
	DELETE(context);

	return FALSE;
}

/* GetArray : get the array value  and the unit	from {all}	*/
/* all : id with(out) context								*/
/* return value = 1: the value found						*/

int GetArray(const char *all,VectorPtr *val,char *unit)
{
	char 			*id = NULL, *context = NULL;
	FLRSymbolPtr	item = NULL;
	int				pos = 0;
	
	if (gFLRArrayList->fNumItems == 0) { 
		fprintf(stderr,"# GetArray: item \"%s\" not found (empty arraylist).\n",all);
		return FALSE;
	}
	id = NEW2(strlen(all) +2, char);
	if ( HasContext(all)) 	{
		context = NEW2(strlen(all)+2, char);
		GetIDandContext(all,id,context);
	}
	else { 
		context = NEW2(2, char);
		strcpy(context,":");
		strcpy(id,all);
	}

	item = RFindItem(gFLRArrayList, &pos, id, gleichID );
	while (item) {
		if ( (FLRCompare(context,item->context) == 0) ) {
			DELETE(id);
			DELETE(context);
			if ( unit ) {
				if ( item->unit ) strcpy( unit, item->unit );		/* return a reference here */
				else unit[0] = '\0';
			}
			*val = item->val.v;				
			return TRUE;
		}
		else { 
			item = RFindItem(gFLRArrayList, &pos, id, gleichID );
		}	
	}
#	ifdef qDebug
	fprintf(stderr,"# GetArray: item \"%s\" not found.\n",id);
#	endif
	DELETE(id);
	DELETE(context);
	return FALSE;

}


static void FreeSymbol(void *ptr, void *aux)
{	
	FLRSymbolPtr	s = (FLRSymbolPtr)ptr;

#ifdef applec
#pragma unused(aux);
#endif
	DELETE(s->id);
	DELETE(s->unit);
	DELETE(s->context);

	switch ( s->type ) {
		case kFLRString:
			if ( !s->detached ) DELETE(s->val.s);
			break;
		case kFLRArray:
			if ( !s->detached ) DisposeVector( s->val.v );
			break;
		default:
			break;
	}

	DELETE( s );
}


/*GetPatternedItems : return a list of matched items with pattern pat	*/

ListPtr GetPatternedItems( char *pat )
{
  	const char 		*compile_ret = NULL;
	char 			fastmap[(1 << BYTEWIDTH)];
	struct 
	re_registers	regs;
	int 			match_ret_id, match_ret_context, match_ret_2;
	size_t			len1,len2;
	LinkPtr			temp = gFLRSymbolList->fHead;
	FLRSymbolPtr	sym = NULL;
	ListPtr			list = NULL;

	gRe_buf.allocated = 0;
	gRe_buf.buffer = NULL;
	gRe_buf.fastmap = fastmap;
	gRe_buf.translate = glookup;

	re_set_syntax(	RE_BACKSLASH_ESCAPE_IN_LISTS	|
					RE_INTERVALS					|
					RE_NO_BK_BRACES					|
					RE_NO_BK_PARENS					|
					RE_DOT_NEWLINE					|
					RE_NO_BK_VBAR					|
					RE_CONTEXT_INDEP_ANCHORS
				);

	compile_ret = re_compile_pattern( pat, strlen(pat), &gRe_buf );
	if ( compile_ret != NULL ) {
		fprintf( stderr, "### \"%s\": <%s>\n", pat, compile_ret );
	}

	list = NewList("matched items");
	while ( temp ) {
		sym = temp->fItem;
		len1 = strlen(sym->id);
		len2 = strlen(sym->context);
		if ( ( match_ret_id = re_match( &gRe_buf,sym->id,len1,0,&regs ) ) == len1 ) {
			AddItem( list, sym, kListTail );
		}
		else if ( ( match_ret_context = re_match( &gRe_buf,sym->context,len2,0,&regs ) ) 
			== len2 ) {
			AddItem( list, sym, kListTail );
		} 
		else if ( ( match_ret_2 = re_match_2( &gRe_buf, sym->context, len2, 
			sym->id, len1, 0, &regs, len1+len2 ) ) == len1+len2 ) {
			AddItem( list, sym, kListTail );
		}
#		ifdef qDebug
		else {
			fprintf( stderr, "# neither \"%s\" nor \"%s\" matched.\n", 
				sym->id, sym->context );
		}
#		endif

		if ( (match_ret_id == -2 || match_ret_context == -2) && match_ret_2 == -2 ) {
          	fprintf( stderr, "# Get_patterned_Items: re_match failed.\n");
          	ForAll( list, FreeSymbol, NULL, kListHead );
			RemoveList( list );
			return NULL;
		}
		temp = temp->fNext;
	}

	if ( list->fNumItems == 0 ) {	/* don't return an empty list */
		RemoveList( list );
		return NULL;
	}
	
	return list;
}

/*GetItems : is not used */

static int GetItems(ListPtr list)
{
	struct 
	re_registers	regs;
	int 			match_ret_id, match_ret_context, match_ret_2;
	size_t			len1,len2;
	LinkPtr			temp = gFLRSymbolList->fHead;
	FLRSymbolPtr		temp2 = NULL;
	
	while(temp) {
		temp2 = temp->fItem;
		len1 = strlen(temp2->id);
		len2 = strlen(temp2->context);
		match_ret_id  = 
			re_match( &gRe_buf,temp2->id,len1,0,&regs );
		match_ret_context  = 
			re_match( &gRe_buf,temp2->context,len2,0,&regs );
		match_ret_2 = 
			re_match_2( &gRe_buf, temp2->context, len2, temp2->id, len1, 0, &regs, 
			len1+len2 );
		if ( (match_ret_id == -2 || match_ret_context == -2) && match_ret_2 == -2 ) {
          	fprintf (stderr, "# GetItems: re_match failed.\n");
          	return FALSE;
		}
		if ( match_ret_id == len1 || match_ret_context == len2 || match_ret_2 == len1+len2 ) {
			fprintf( stderr, "# matched %d chars in id, %d chars in context, %d chars in full name.\n", 
				match_ret_id, match_ret_context, match_ret_2 );
			AddItem(list,temp2, kListTail);	
		}
		temp = temp->fNext;
	}
	return TRUE;
}



/* InitFLReader : routine to initialize the Flamelet reader */

void InitFLReader(void)
{
#ifdef applec
	SpinCursor(0);
#endif

	FLRStackInit();
	InitFloatBuffers(MaxElems);
	gFLRSymbolList = NewList("Symbol");
	gFLRIntList = NewList("Integer");
	gFLRFloatList = NewList("Float");
	gFLRNumberList = NewList("Numbers");
	gFLRStringList = NewList("Strings");
	gFLRArrayList = NewList("Arrays");
	gFLRHeaderList = NewList( "Header" );
	LookupTable(); /* initialize the lookup table for ASCII character */
	
	gSection = kInitialState;	/*	make sure we can scan the file		*/
	gFLLine = 1;					/*	always start at line 1 with a new file	*/
}

/* ReadFlamelet : read a flamelet file from file pointer fp */

int ReadFlamelet( FILE *fp )
{
	flin = fp;

#	ifdef qDebugScanner
	while ( fllex() );
#	endif

	flrestart( flin );
	return ( flparse() == 0 );
}

/* CleanupFLReader : clean up the memory usage from FLReader  */

void CleanupFLReader(void)
{
	Free1DArray(gFBuf);
	
	/*Delete all items ; Delete all lists	*/

	ForAll(gFLRSymbolList,FreeSymbol,NULL,kListHead);
	
	RemoveList(gFLRSymbolList);
	RemoveList(gFLRNumberList);
	RemoveList(gFLRIntList);
	RemoveList(gFLRFloatList);
	RemoveList(gFLRStringList);
	RemoveList(gFLRArrayList);
	RemoveList( gFLRHeaderList );
	
	FLRStackDelete();

}

/* ChangeString : change string value from an ID 				*/
/* all : id with(out) context (regular expression pattern)		*/
/* newstring	: new string value								*/
/* return value = 1 : value of the id successfull changed 		*/
/* return value = 0 : id is not found in list 					*/

int ChangeString(const char *all,char *newstring) 
{
	char		*id,*context;
	FLRSymbolPtr	item;
	int			pos = 0;
	
	if (gFLRStringList->fNumItems == 0) { 
		fprintf(stderr,"# ChangeString: item \"%s\" not found (empty stringlist).\n",all);
		return 0;
	}
	
	id = NEW2(strlen(all) +2, char);
	
	if ( HasContext(all)) 	{
		context = NEW2(strlen(all)+2, char);
		GetIDandContext(all,id,context);
	}
	else { 
		context = NEW2(2, char);
		strcpy(context,":");
		strcpy(id,all);
	}
	
	item = RFindItem(gFLRStringList, &pos, id, gleichID );
	
	while (item) {
		if ( (FLRCompare(context,item->context) == 0) ){ 
			DELETE(id);
			DELETE(context);
			DELETE(item->val.s);

			item->val.s = NEW2(strlen(newstring)+1,char);
			strcpy(item->val.s,newstring);
			return TRUE;
		}
		else { 
			item = RFindItem(gFLRStringList, &pos, id, gleichID );
		}	
	}
#	ifdef qDebug
	fprintf(stderr,"# ChangeString: item \"%s\" not found.\n",id);
#	endif

	/* Clean Up  */
	DELETE(id);
	DELETE(context);

	return FALSE;
}

/* ChangeNumber : change number value( integer or float) from an ID 	*/
/* all : id with(out) context (regular expression pattern)				*/
/* val	: new value														*/
/* unit : unit = NULL if the unit is to be changed						*/
/* return value = 1 : value of the id successfull changed 				*/
/* return value = 0 : id is not found in list 							*/

int ChangeNumber(const char *all,const Double val, char *unit)
{
	char 			*id = NULL, *context = NULL;
	FLRSymbolPtr	item = NULL;
	int  			pos = 0;
	
	if (gFLRNumberList->fNumItems == 0) { 
		fprintf(stderr,"# ChangeNumber: item \"%s\" not found (empty numberlist).\n",all);
		return 0;
	}

	id = NEW2(strlen(all) +2, char);
	if ( HasContext(all)) 	{
		context = NEW2(strlen(all)+2, char);
		GetIDandContext(all,id,context);
	}
	else { 
		context = NEW2(2, char);
		strcpy(context,":");
		strcpy(id,all);
	}
	item = RFindItem(gFLRNumberList, &pos, id, gleichID );
	while (item) {
		if ( (FLRCompare(context,item->context) == 0) ) {
			if (item->type == kFLRInteger) {
				DELETE(id);
				DELETE(context);
				if ( unit ) {
					DELETE( item->unit );
					item->unit = NEW2(strlen(unit)+1,char);
					strcpy( item->unit, unit );
				}
				item->val.i = (int)val;
				return TRUE;
			}
			if (item->type == kFLRFloat)  {
				DELETE(id);
				DELETE(context);
				if ( unit ) {
					DELETE(item->unit);
					item->unit = NEW2(strlen(unit)+1,char);
					strcpy(item->unit,unit);
				}
				item->val.f = val;
				return TRUE;
			}
		}
		else { 
			item = RFindItem(gFLRNumberList, &pos, id, gleichID );
		}	
	}
#	ifdef qDebug
	fprintf(stderr,"# ChangeNumber: item \"%s\" not found.\n",id);
#	endif

	/* Clean Up  */
	DELETE(id);
	DELETE(context);
	
	return FALSE;
}

/* ChangeArray : change array value( from type VectorPtr) from an ID 	*/
/* all : id with(out) context (regular expression pattern)				*/
/* val	: new value														*/
/* unit : unit = NULL if the unit is to be changed						*/
/* return value = 1 : value of the id successfull changed 				*/
/* return value = 0 : id is not found in list 							*/

int ChangeArray( const char *all, VectorPtr val, char *unit )
{
	char 			*id = NULL, *context = NULL;
	FLRSymbolPtr	item = NULL;
	int				pos = 0;
	
	if (gFLRArrayList->fNumItems == 0) { 
		fprintf(stderr,"# ChangeArray: item \"%s\" not found (empty arraylist).\n",all);
		return 0;
	}
	id = NEW2(strlen(all) +2, char);
	if ( HasContext(all)) 	{
		context = NEW2(strlen(all)+2, char);
		GetIDandContext(all,id,context);
	}
	else { 
		context = NEW2(2, char);
		strcpy(context,":");
		strcpy(id,all);
	}
	item = RFindItem(gFLRArrayList, &pos, id, gleichID );
	while (item) {
		if ( (FLRCompare(context,item->context) == 0) ) {
			DELETE(id);
			DELETE(context);
			if ( unit ) {
				DELETE(item->unit);
				item->unit = NEW2(strlen(unit)+1,char);
				strcpy(item->unit,unit);
			}
			DisposeVector( item->val.v );
			item->val.v = val;
			return TRUE;
		}
		else { 
			item = RFindItem(gFLRArrayList, &pos, id, gleichID );
		}	
	}

#	ifdef qDebug
	fprintf(stderr,"# ChangeArray: item \"%s\" not found.\n",id);
#	endif

	DELETE(id);
	DELETE(context);

	return FALSE;
}

/* DeleteSymbol : delete the symbol "{all}" 			*/
/* all : id with(out) context (regular expression) 		*/
/* return value = 1: the symbol is found and deleted	*/
/* return value = 0: the symbol is not found in list	*/

int DeleteSymbol( const char *all )
{
	char 			*id = NULL, *context = NULL;
	FLRSymbolPtr	item = NULL;
	int  			pos = 0;
	
	if (gFLRSymbolList->fNumItems == 0) { 
		fprintf(stderr,"# DeleteSymbol: item \"%s\" not found (empty Symbol list).\n",all);
		return FALSE;
	}

	id = NEW2(strlen(all) +2, char);
	if ( HasContext(all)) 	{
		context = NEW2(strlen(all)+2, char);
		GetIDandContext(all,id,context);
	}
	else { 
		context = NEW2(2, char);
		strcpy(context,":");
		strcpy(id,all);
	}

	item = RFindItem(gFLRSymbolList, &pos, id, gleichID );
	while (item) {
		if ( (FLRCompare(context,item->context) == 0) ) {
			switch ( item->type ) {
				case kFLRString:
					if ( !item->detached ) DELETE(item->val.s);
					RemoveItem(gFLRStringList,item);
					break;
				case kFLRArray:
					if ( !item->detached ) DisposeVector(item->val.v);
					RemoveItem(gFLRArrayList,item);
					break;
				case kFLRInteger:
					RemoveItem(gFLRIntList,item);
					RemoveItem(gFLRNumberList,item);
					break;
				case kFLRFloat:
					RemoveItem(gFLRFloatList,item);
					RemoveItem(gFLRNumberList,item);
					break;
				default:
					break;
			}

			RemoveItem(gFLRSymbolList,item);
			DELETE(item->id);
			DELETE(item->unit);
			DELETE(item->context);

			DELETE(item);
			DELETE(id);
			DELETE(context);

			return TRUE;
		}
		else { 
			item = RFindItem(gFLRSymbolList, &pos, id, gleichID );
		}	
	}

#	ifdef qDebug
	fprintf(stderr,"# DeleteSymbol: item \"%s\" not found.\n",id);
#	endif

	/* Clean Up  */
	DELETE(id);
	DELETE(context);

	return FALSE;
}

/* WriteFlameletFile:  write the symbollist(global) in flamelet format */
/*						to file pointer fp.								*/

void WriteFlameletFile( FILE *fp ,int section)

{
	char 	*buffer;
	int 	k;
	
/*	SortItems(gFLRSymbolList, NULL, SortContext);*/
	
	if ( (section == kFLRHeader) || (section == kFLRAll) ) {
		fprintf( fp,"Header\n\n");

		ForAll( gFLRHeaderList, PrintFLHeader, fp, kListHead );
#		ifdef NEW_PRINT_HEADER
		while ( gLastDepth  ) {		/* clean up any open context */
			--gLastDepth;
			CloseContext( gLastDepth, fp );
		}
#		endif
		
		while (!FLRStackEmpty()) {
			strcpy(gFormat,"");
			for(k = 0; k< (gSn-1); ++k) {
				strcat(gFormat,"\t");
			}
			fprintf(fp,"%sEnd\n",gFormat);
			buffer = FLRStackPop();
			--gSn;
			DELETE(buffer);
		}
	}

	if ( (section == kFLRBody) || (section == kFLRAll) ) {
		fprintf(fp,"\n\nBody");
		ForAll(gFLRArrayList,PrintFLBody,fp,kListHead);
	}

	if ( (section == kFLRTrailer) || (section == kFLRAll) ) {
		fprintf(fp,"\nTrailer\n");
	}
}
