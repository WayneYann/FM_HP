%{

/*
	FLScan.l
	Scanner file for the Flamelet reader package
	
	© Josef Gšttgens, Peter Terhoeven, Ian Herwono, 1993
	version 1.b1
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef applec
#include <CursorCtl.h>
#endif

#include "ArrayManager.h"
#include "alligator.h"
#include "Stack.h"
#include "regex.h"
#include "List.h"

#include "FLReader.h"
#include "FLScan.h"
#include "FLScan.tab.h"


#undef	YY_DECL
#define	YY_DECL	int fllex( void )

#define yylval fllval


/*
	Global variables
*/
int	gSection = kInitialState;	/* always reinitialized in InitFLReader()	*/
int	gExpectArray = FALSE;
int	gFLLine = 1;



/*#define qDebug*/

#ifdef qDebug
# define PKey(s)	fprintf( gFlex, "Keyword    \"%s\"\n", (s) );
# define PBool(b)	fprintf( gFlex, "Boolean    \"%d\"\n", (b) );
# define PID(id)	fprintf( gFlex, "Identifier \"%s\"\n", (id) );
# define PInt(d)	fprintf( gFlex, "Integer    \"%d\"\n", (d) );
# define PFloat(x)	fprintf( gFlex, "Float      \"%g\"\n", (x) );
# define PUnit(u)	fprintf( gFlex, "Unit       \"%s\"\n", (u) );
# define PString(s)	fprintf( gFlex, "String     \"%s\"\n", (s) );
#else
# define PKey(s)
# define PBool(b)
# define PID(id)
# define PInt(d)
# define PFloat(x)
# define PUnit(u)
# define PString(s)
#endif

#ifdef applec
#undef YY_BREAK
#define YY_BREAK	SpinCursor(1); break;
#endif

static void CanonicalUnits( char *unit );
static void RemoveQuotes( char *string );

%}

ID			[a-zA-Z_][_*/a-zA-Z0-9-]*
INT			[-+]?[0-9]+
FLOAT1		[-+]?[0-9]+"."[0-9]*
FLOAT2		[-+]?"."[0-9]+
EXP			[eE][-+]?[0-9]+
UNIT1		"["[^\n]+"]"
UNIT2		"("[^\n]+")"	
STRING		"\""[^\n]*"\""

%%
[hH][eE][aA][dD][eE][rR]		{	PKey(yytext); return TokHeader;		}
[bB][oO][dD][yY]				{	PKey(yytext); return TokBody;		}
[tT][rR][aA][iI][lL][eE][rR]	{	PKey(yytext); return TokTrailer;	}
([bB][eE][gG][iI][nN])|"{"		{	PKey(yytext); return TokBegin;		}
([eE][nN][dD])|"}"				{	PKey(yytext); return TokEnd;		}
[tT][rR][uU][eE]				{	yylval.bool = 1;
			  						PBool( yylval.bool );
			  						return TokBoolean;
								}
[fF][aA][lL][sS][eE]			{	yylval.bool = 0;
									PBool( yylval.bool );
									return TokBoolean;
								}

[\n]+							{	gFLLine += strlen((char *)yytext);
									if (gExpectArray == FALSE) return '\n';
								}
"="								{	return '=';	}

{ID}							{	yylval.string = NEW2(yyleng+1, char);
									strcpy( yylval.string, (char *)yytext );
									PID( yylval.string );
									return TokID;
								}
{INT}							{	yylval.i = atoi((char *)yytext);
			  						PInt( yylval.i );
									return TokInt;
								}
{FLOAT1}		|
{FLOAT2}		|
{FLOAT1}{EXP}	|
{FLOAT2}{EXP}	|
{INT}{EXP}						{	
									yylval.f = atof((char *)yytext);
									PFloat( yylval.f );
									return TokFloat;
								}
{UNIT1}|{UNIT2}					{	yylval.string = NEW2(yyleng+1,char);
			 						strcpy(yylval.string,(char *)yytext);
									CanonicalUnits(yylval.string);
									PUnit( yylval.string );
									return TokUnit;
								}
{STRING}						{	yylval.string = NEW2(yyleng+1, char);
									strcpy( yylval.string, (char *)yytext );
									RemoveQuotes( yylval.string );
									PString( yylval.string );
									return TokString;
								}
[ \t]+							{	/*  Eat up white space characters */ ;	}

.								{ 	/* Eat up unmatched characters.  */
									fprintf( stderr, "Line %d #  Unmatched character: '%s' (ASCII %d)\n", 
										gFLLine, yytext, (int)yytext[0] );
								}

%%


int flwrap( void )
{
	return 1;
}


static void CanonicalUnits( char *unit )
{
/*	CanonicalUnits() strips off the brackets from the unit.string
    and removes all white-space characters.
*/
    char *ptr = unit;    
    
    if ( (*unit != '[') && (*unit != '(') ) return;
    ++unit;
    while ( (*unit != ']') && (*unit != ')') ) {
        if ( (*unit != ' ') && (*unit != '\t') )	*ptr++ = *unit++;
		else										++unit;
    }
    *ptr = '\0';
}


static void RemoveQuotes( char *string )
{
/*	RemoveQuotes() removes the leading and trailing quotation
    marks from string.
*/
    char *ptr = string;
    
    if ( *string++ != '"' ) return;
    
    while ( *string != '"' ) *ptr++ = *string++;
    *ptr = '\0';
}

