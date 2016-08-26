%{

/*
	FLScan.y
	Parser file for Flamelet reader package
	
	© Josef Gšttgens, Peter Terhoeven, Ian Herwono, 1993
	version 1.b1
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ArrayManager.h"
#include "alligator.h"
#include "Stack.h"
#include "regex.h"
#include "List.h"

#include "FLReader.h"
#include "FLPatch.h"
#include "FLScan.h"


/*#define qOutput*/

#ifdef qOutput
# define PKey(s)				fprintf(gBison,"%s\n",(s))
# define PID_i(id,d)			fprintf(gBison,"\t%s = %d\n",(id),(d))
# define PID_i_unit(id,d,s)		fprintf(gBison,"\t%s = %d [%s]\n",(id),(d),(s))
# define PID_f(id,x)			fprintf(gBison,"\t%s = %lf\n",(id),(x))
# define PID_f_unit(id,x,s)		fprintf(gBison,"\t%s = %lf [%s]\n",(id),(x),(s))
# define PID_string(id,s)		fprintf(gBison,"\t%s = \"%s\"\n",(id),(s))
# define PID_bool(id,d)			fprintf(gBison,"\t%s = %d\n",(id),(d))
# define PID_begin(id)			fprintf(gBison,"\n%s\nBegin\n",(id))
# define PEnd(context)			fprintf(gBison,"End of context %s\n\n", context)
# define PBody_ID(id)			fprintf(gFloatBuffer,"%s\n",(id))
# define PBody_f(x)				fprintf(gFloatBuffer,"%lf\n",(x))
# define PBody_ID_Unit(id,s)	fprintf(gFloatBuffer,"%s [%s]\n",(id),(s))
# define PID_err(id)			fprintf(gBison,"\t%s =  <- not complete\n",(id))
# define PID_err_unit(id,s)		fprintf(gBison,"\t%s %s =  <- not complete\n",(id),(s))
# define fWarning(s)			fprintf(gBison,"%s",(s));
#else
# define PKey(s)				
# define PID_i(id,d)			
# define PID_i_unit(id,d,s)	
# define PID_f(id,x)		
# define PID_f_unit(id,x,s)	
# define PID_string(id,s)	
# define PID_bool(id,d)		
# define PID_begin(id)		
# define PEnd(context)				
# define PBody_ID(id)		
# define PBody_f(x)			
# define PBody_ID_Unit(id,s)	
# define PID_err(id)
# define PID_err_unit(id,s)
# define fWarning(s)		fprintf(stderr,"%s",(s));
#endif
#define Warning(s)		fprintf(stderr,"## Warning: %s near line %d !!\n",(s),gFLLine);

void yyerror(char *s);

static char	*contextBuf;


void CheckSection( int which )
{
	if ( gSection != which ) {
		fprintf( stderr, "###  Syntax error near line %d.\n", gFLLine );
		/*	This should go away one of these days...			*/
		FATAL( "Syntax error" );
	}
}

void checkStack()  /* check whether all contexts have been popped */
{
	if( !FLRStackEmpty() ) {
		Warning("Unterminated context(s)");
	}
}

void checkS(char *s) 
{
	if((strlen(s) == 0)) {
		Warning("empty string");
	}
}

%}

%union {
	char		*string;
	int		 	bool;
	int			i;
	double		f;
}

%token	TokHeader
%token	TokBody
%token	TokTrailer
%token	TokBegin
%token	TokEnd
%token	<bool>  	TokBoolean
%token	<i>			TokInt
%token	<f>			TokFloat
%token  <string>	TokUnit
%token 	<string>	TokID
%token 	<string> 	TokString

%%

input : /* empty */
		|	input statement
	
;	

statement : '\n'
		|	hstatement
		|	bstatement				
		|	tstatement
;

hstatement : 
			TokHeader		{	PKey("Header");
								if ( gSection == kInitialState ) gSection = kInHeader;
								else	FATAL( "Unexpected keyword \"Header\"" );
							}
		|	begincontext	{	CheckSection(kInHeader);	}
		|	endcontext		{	CheckSection(kInHeader);	}
		|	assignment		{	CheckSection(kInHeader);	}
;	

assignment : 
	  TokID '=' TokInt '\n'				{	PID_i($1,$3);
	  										PutSymH($1,NULL,kFLRInteger,$3,0,NULL);
										}
	| TokID '=' TokInt TokUnit '\n' 	{	PID_i_unit($1,$3,$4);
											PutSymH($1,$4,kFLRInteger,$3,0,NULL);
										}
	| TokID TokUnit '=' TokInt '\n'		{	PID_i_unit($1,$4,$2);
											PutSymH($1,$2,kFLRInteger,$4,0,NULL);
										}
	| TokID '=' TokFloat '\n' 			{	PID_f($1,$3);
											PutSymH($1,NULL,kFLRFloat,0,$3,NULL);
										}
	| TokID '=' TokFloat TokUnit '\n'	{	PID_f_unit($1,$3,$4);
											PutSymH($1,$4,kFLRFloat,0,$3,NULL);
										}
	| TokID TokUnit '=' TokFloat '\n'	{	PID_i_unit($1,$4,$2);
											PutSymH($1,$2,kFLRFloat,0,$4,NULL);
										}
	| TokID '=' TokString '\n'			{	checkS($3); PID_string($1,$3);
											PutSymH($1,NULL,kFLRString,0,0,$3);
										}
	| TokID '=' TokBoolean '\n' 		{	PID_bool($1,$<bool>3);
											PutSymH($1,NULL,kFLRBoolean,$3,0,NULL);
										}
	| error '\n'       					{	yyerrok;	}
;

begincontext: 
			TokID TokBegin	'\n'			{	FLRStackPush($1); PID_begin($1); DELETE($1);	}
		|	TokID '\n' TokBegin '\n'		{	FLRStackPush($1); PID_begin($1); DELETE($1);	}
		|	TokID TokID TokBegin '\n'		{	FLRStackPush($1); PID_begin($1); DELETE($1); DELETE($2);	}
		|	TokID TokID '\n' TokBegin '\n'	{	FLRStackPush($1); PID_begin($1); DELETE($1); DELETE($2);	}
;

endcontext:
			TokEnd '\n'  				{	if ( !FLRStackEmpty() ) {
												contextBuf = FLRStackPop();
												PEnd( contextBuf );
												DELETE(contextBuf);
											}
											else FATAL("No context to FLRStackPop()");
										}		
;

bstatement:
			TokBody				{	checkStack();
									PKey("\nBody");
									if ( gSection == kInHeader ) gSection = kInBody;
									else	FATAL( "Unexpected keyword \"Body\"" );
								}
		|	vassignment 		{	CheckSection(kInBody);
									/*PrintFloatBuffer();*/
									ResetFloatBuffer();
									gExpectArray = FALSE;
								}	
;

vassignment :
			TokID '\n'  vector				{	PBody_ID($1);
												PutSymBody($1,NULL,gFBufhead,gFBElems);
											}
		|	TokID TokUnit '\n' vector		{	PBody_ID_Unit($1,$2);
												PutSymBody($1,$2,gFBufhead,gFBElems);
											}
		| 	TokID '=' '\n' vector			{	PBody_ID($1);
												PutSymBody($1,NULL,gFBufhead,gFBElems);
											}
		|	TokID TokUnit '=' '\n' vector	{	PBody_ID_Unit($1,$2);
												PutSymBody($1,$2,gFBufhead,gFBElems);
											}
;

vector :	/* empty */
		|	vector TokFloat 				{	AddNumber($2);	}
;

tstatement :
			TokTrailer 						{	PKey("Trailer");
												if ( gSection == kInBody ) gSection = kInTrailer;
												else	FATAL( "Unexpected keyword \"Trailer\"" );
												/* do not process trailer contents */
												/*flrestart();*/
												YYACCEPT;
											}
;

%%


void yyerror(char *s)
{
	fprintf(gBison,"%s\n",s);
}

