#ifndef lint
static const char yysccsid[] = "@(#)yaccpar	1.9 (Berkeley) 02/21/93";
#endif

#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define YYPATCH 20101229

#define YYEMPTY        (-1)
#define yyclearin      (yychar = YYEMPTY)
#define yyerrok        (yyerrflag = 0)
#define YYRECOVERING() (yyerrflag != 0)

#define YYPREFIX "yy"

#define YYPURE 0

#line 1 "ReactionScan.y"
  
/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 */

#include "ReactionParser.h"
#line 11 "ReactionScan.y"
#ifdef YYSTYPE
#undef  YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
#endif
#ifndef YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
typedef union {
	Double 			*ptrdouble;
	Double 			typdouble;
	int 			typint;
	char 			typstring[10000];
	ReactionPtr		typreaction;
	Dimension		typdimension;
	DimensionPtr	typdimensionptr;
	Exponent		typexponent;
} YYSTYPE;
#endif /* !YYSTYPE_IS_DECLARED */
#line 46 "y.tab.c"
/* compatibility with bison */
#ifdef YYPARSE_PARAM
/* compatibility with FreeBSD */
# ifdef YYPARSE_PARAM_TYPE
#  define YYPARSE_DECL() yyparse(YYPARSE_PARAM_TYPE YYPARSE_PARAM)
# else
#  define YYPARSE_DECL() yyparse(void *YYPARSE_PARAM)
# endif
#else
# define YYPARSE_DECL() yyparse(void)
#endif

/* Parameters sent to lex. */
#ifdef YYLEX_PARAM
# define YYLEX_DECL() yylex(void *YYLEX_PARAM)
# define YYLEX yylex(YYLEX_PARAM)
#else
# define YYLEX_DECL() yylex(void)
# define YYLEX yylex()
#endif

/* Parameters sent to yyerror. */
#define YYERROR_DECL() yyerror(const char *s)
#define YYERROR_CALL(msg) yyerror(msg)

extern int YYPARSE_DECL();

#define FIRSTTOKEN 257
#define LETTERS 258
#define FACTOR 259
#define SPECIES 260
#define THIRDBODY 261
#define RIGHTARROW 262
#define LEFTARROW 263
#define EQUAL 264
#define LEFTBRACE 265
#define SPECIESCOEFF 266
#define MI_NUMBER 267
#define LABEL 268
#define MI_LITERAL 269
#define E_ALPHANUM 270
#define RIGHTHIGH 271
#define UNARYMINUS 272
#define YYERRCODE 256
static const short yylhs[] = {                           -1,
    0,    0,    0,    0,    0,    0,    0,    0,    4,    4,
    4,    4,    4,    4,    3,    3,    3,    3,    3,    3,
    3,    3,    3,    1,    1,    1,    1,    1,    1,    1,
    1,    2,    2,
};
static const short yylen[] = {                            2,
    0,    4,    4,    4,    5,    4,    4,    3,    3,    3,
    3,    2,    3,    1,    3,    2,    4,    3,    3,    3,
    3,    2,    1,    1,    3,    2,    2,    2,    2,    2,
    2,    2,    3,
};
static const short yydefred[] = {                         1,
    0,    0,    0,   24,    0,    0,    0,    0,    0,    8,
    0,   14,    0,    0,    0,   16,    0,   22,    0,    0,
   26,   27,   29,   30,   28,   31,    0,    0,    0,    0,
    0,    0,    0,    0,    4,   12,    0,    0,    0,    0,
    3,    0,   21,    6,   25,   32,    0,    0,    2,    7,
    0,    0,   19,   15,    9,    0,    0,    0,   17,    5,
   33,
};
static const short yydgoto[] = {                          1,
    8,   29,    9,   15,
};
static const short yysindex[] = {                         0,
  -40,  -21,  -37,    0,  -38,  -39,  -39, -168,  -12,    0,
   -6,    0,  -33,  -33,  -30,    0,  -39,    0,   14, -117,
    0,    0,    0,    0,    0,    0, -226, -224, -121,    3,
  -39,  -39,  -39,  -39,    0,    0,   -3,  -33,  -33,  -39,
    0,   61,    0,    0,    0,    0,  -73, -207,    0,    0,
  -15,  -15,    0,    0,    0,  -31,  -31,  -24,    0,    0,
    0,
};
static const short yyrindex[] = {                         0,
    0,    0,    0,    0,  -32,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
  -23,  -17,    0,    0,    0,    4,   68,   -5,    0,    0,
    0,
};
static const short yygindex[] = {                         0,
    0,    0,   41,   40,
};
#define YYTABLESIZE 244
static const short yytable[] = {                          7,
    7,   17,   14,   49,    6,    6,   14,   44,   23,   23,
   23,   39,   23,   23,   23,   41,   38,   18,   31,   18,
   32,   18,   18,   20,   10,   20,   34,   20,   20,   34,
   31,   33,   32,   45,   33,   13,   13,   55,   39,   35,
   13,   13,   46,   38,   10,   10,   18,   19,   50,   10,
   10,   60,   36,   37,   43,   34,   31,   42,   32,   61,
   33,   23,   40,   40,    0,    0,    0,    0,    0,    0,
   18,   51,   52,   53,   54,    0,   20,   56,   57,    0,
   58,    0,    0,    0,    0,    0,    0,   20,   13,    0,
   40,   21,   22,   23,   24,   25,   26,   27,    0,    0,
   28,   59,   34,   31,    0,   32,    0,   33,   11,   11,
    0,    0,    0,   11,   11,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,   47,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,   48,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    2,    0,    3,   11,    0,
   12,   13,    0,   23,   12,   13,    0,    4,    0,    5,
    5,   16,   18,    0,    0,    0,    0,    0,   20,    0,
    0,    0,    0,   30,
};
static const short yycheck[] = {                         40,
   40,   40,   40,  125,   45,   45,   40,  125,   41,   42,
   43,   42,   45,   46,   47,   46,   47,   41,   43,   43,
   45,   45,   46,   41,   46,   43,   42,   45,   46,   42,
   43,   47,   45,  260,   47,   41,   42,   41,   42,   46,
   46,   47,  267,   47,   41,   42,    6,    7,   46,   46,
   47,  125,   13,   14,   41,   42,   43,   17,   45,  267,
   47,   94,   94,   94,   -1,   -1,   -1,   -1,   -1,   -1,
   94,   31,   32,   33,   34,   -1,   94,   38,   39,   -1,
   40,   -1,   -1,   -1,   -1,   -1,   -1,  256,   94,   -1,
   94,  260,  261,  262,  263,  264,  265,  266,   -1,   -1,
  269,   41,   42,   43,   -1,   45,   -1,   47,   41,   42,
   -1,   -1,   -1,   46,   47,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,  256,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,  269,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,  256,   -1,  258,  256,   -1,
  258,  259,   -1,  256,  258,  259,   -1,  268,   -1,  270,
  270,  270,  256,   -1,   -1,   -1,   -1,   -1,  256,   -1,
   -1,   -1,   -1,  256,
};
#define YYFINAL 1
#ifndef YYDEBUG
#define YYDEBUG 1
#endif
#define YYMAXTOKEN 272
#if YYDEBUG
static const char *yyname[] = {

"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,"'('","')'","'*'","'+'",0,"'-'","'.'","'/'",0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"'^'",0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"'}'",0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,"FIRSTTOKEN","LETTERS","FACTOR","SPECIES","THIRDBODY","RIGHTARROW",
"LEFTARROW","EQUAL","LEFTBRACE","SPECIESCOEFF","MI_NUMBER","LABEL","MI_LITERAL",
"E_ALPHANUM","RIGHTHIGH","UNARYMINUS",
};
static const char *yyrule[] = {
"$accept : statement",
"statement :",
"statement : statement reaction moreinfo '}'",
"statement : statement LETTERS dimension '.'",
"statement : statement LETTERS error '.'",
"statement : statement reaction moreinfo error '}'",
"statement : statement reaction error '}'",
"statement : statement exp error '.'",
"statement : statement error '.'",
"dimension : '(' dimension ')'",
"dimension : dimension '/' dimension",
"dimension : dimension '*' dimension",
"dimension : FACTOR dimension",
"dimension : dimension '^' exp",
"dimension : LETTERS",
"exp : exp '*' exp",
"exp : E_ALPHANUM E_ALPHANUM",
"exp : E_ALPHANUM '(' exp ')'",
"exp : exp '+' exp",
"exp : exp '/' exp",
"exp : exp '-' exp",
"exp : '(' exp ')'",
"exp : '-' exp",
"exp : E_ALPHANUM",
"reaction : LABEL",
"reaction : reaction SPECIESCOEFF SPECIES",
"reaction : reaction SPECIES",
"reaction : reaction THIRDBODY",
"reaction : reaction EQUAL",
"reaction : reaction RIGHTARROW",
"reaction : reaction LEFTARROW",
"reaction : reaction LEFTBRACE",
"moreinfo : MI_LITERAL MI_NUMBER",
"moreinfo : moreinfo MI_LITERAL MI_NUMBER",

};
#endif
/* define the initial stack-sizes */
#ifdef YYSTACKSIZE
#undef YYMAXDEPTH
#define YYMAXDEPTH  YYSTACKSIZE
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 500
#define YYMAXDEPTH  500
#endif
#endif

#define YYINITSTACKSIZE 500

int      yydebug;
int      yynerrs;

typedef struct {
    unsigned stacksize;
    short    *s_base;
    short    *s_mark;
    short    *s_last;
    YYSTYPE  *l_base;
    YYSTYPE  *l_mark;
} YYSTACKDATA;
int      yyerrflag;
int      yychar;
YYSTYPE  yyval;
YYSTYPE  yylval;

/* variables for the parser stack */
static YYSTACKDATA yystack;
#line 328 "ReactionScan.y"

#line 283 "y.tab.c"

#if YYDEBUG
#include <stdio.h>		/* needed for printf */
#endif

#include <stdlib.h>	/* needed for malloc, etc */
#include <string.h>	/* needed for memset */

/* allocate initial stack or double stack size, up to YYMAXDEPTH */
static int yygrowstack(YYSTACKDATA *data)
{
    int i;
    unsigned newsize;
    short *newss;
    YYSTYPE *newvs;

    if ((newsize = data->stacksize) == 0)
        newsize = YYINITSTACKSIZE;
    else if (newsize >= YYMAXDEPTH)
        return -1;
    else if ((newsize *= 2) > YYMAXDEPTH)
        newsize = YYMAXDEPTH;

    i = data->s_mark - data->s_base;
    newss = (short *)realloc(data->s_base, newsize * sizeof(*newss));
    if (newss == 0)
        return -1;

    data->s_base = newss;
    data->s_mark = newss + i;

    newvs = (YYSTYPE *)realloc(data->l_base, newsize * sizeof(*newvs));
    if (newvs == 0)
        return -1;

    data->l_base = newvs;
    data->l_mark = newvs + i;

    data->stacksize = newsize;
    data->s_last = data->s_base + newsize - 1;
    return 0;
}

#if YYPURE || defined(YY_NO_LEAKS)
static void yyfreestack(YYSTACKDATA *data)
{
    free(data->s_base);
    free(data->l_base);
    memset(data, 0, sizeof(*data));
}
#else
#define yyfreestack(data) /* nothing */
#endif

#define YYABORT  goto yyabort
#define YYREJECT goto yyabort
#define YYACCEPT goto yyaccept
#define YYERROR  goto yyerrlab

int
YYPARSE_DECL()
{
    int yym, yyn, yystate;
#if YYDEBUG
    const char *yys;

    if ((yys = getenv("YYDEBUG")) != 0)
    {
        yyn = *yys;
        if (yyn >= '0' && yyn <= '9')
            yydebug = yyn - '0';
    }
#endif

    yynerrs = 0;
    yyerrflag = 0;
    yychar = YYEMPTY;
    yystate = 0;

#if YYPURE
    memset(&yystack, 0, sizeof(yystack));
#endif

    if (yystack.s_base == NULL && yygrowstack(&yystack)) goto yyoverflow;
    yystack.s_mark = yystack.s_base;
    yystack.l_mark = yystack.l_base;
    yystate = 0;
    *yystack.s_mark = 0;

yyloop:
    if ((yyn = yydefred[yystate]) != 0) goto yyreduce;
    if (yychar < 0)
    {
        if ((yychar = YYLEX) < 0) yychar = 0;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, reading %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
    }
    if ((yyn = yysindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: state %d, shifting to state %d\n",
                    YYPREFIX, yystate, yytable[yyn]);
#endif
        if (yystack.s_mark >= yystack.s_last && yygrowstack(&yystack))
        {
            goto yyoverflow;
        }
        yystate = yytable[yyn];
        *++yystack.s_mark = yytable[yyn];
        *++yystack.l_mark = yylval;
        yychar = YYEMPTY;
        if (yyerrflag > 0)  --yyerrflag;
        goto yyloop;
    }
    if ((yyn = yyrindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
        yyn = yytable[yyn];
        goto yyreduce;
    }
    if (yyerrflag) goto yyinrecovery;

    yyerror("syntax error");

    goto yyerrlab;

yyerrlab:
    ++yynerrs;

yyinrecovery:
    if (yyerrflag < 3)
    {
        yyerrflag = 3;
        for (;;)
        {
            if ((yyn = yysindex[*yystack.s_mark]) && (yyn += YYERRCODE) >= 0 &&
                    yyn <= YYTABLESIZE && yycheck[yyn] == YYERRCODE)
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *yystack.s_mark, yytable[yyn]);
#endif
                if (yystack.s_mark >= yystack.s_last && yygrowstack(&yystack))
                {
                    goto yyoverflow;
                }
                yystate = yytable[yyn];
                *++yystack.s_mark = yytable[yyn];
                *++yystack.l_mark = yylval;
                goto yyloop;
            }
            else
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: error recovery discarding state %d\n",
                            YYPREFIX, *yystack.s_mark);
#endif
                if (yystack.s_mark <= yystack.s_base) goto yyabort;
                --yystack.s_mark;
                --yystack.l_mark;
            }
        }
    }
    else
    {
        if (yychar == 0) goto yyabort;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
        yychar = YYEMPTY;
        goto yyloop;
    }

yyreduce:
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: state %d, reducing by rule %d (%s)\n",
                YYPREFIX, yystate, yyn, yyrule[yyn]);
#endif
    yym = yylen[yyn];
    if (yym)
        yyval = yystack.l_mark[1-yym];
    else
        memset(&yyval, 0, sizeof yyval);
    switch (yyn)
    {
case 2:
#line 50 "ReactionScan.y"
	{ }
break;
case 3:
#line 51 "ReactionScan.y"
	{ 	
													DimensionPtr dimension = NULL;
													dimension = AddDimension( gUsedDimensionList, yystack.l_mark[-2].typstring, 0, 0, 0, 0, 0, 0 );
													CopyDimension( dimension, yystack.l_mark[-1].typdimension );
													SetExponentFlags( dimension );
												}
break;
case 4:
#line 59 "ReactionScan.y"
	{ 	
													if ( !gErrorOut ) {
														ScanError( "error in the definition of the units for '%s'", yystack.l_mark[-2].typstring, 0 ); 
													}
													gErrorOut = FALSE;
													fprintf( stderr, "# error in the definition of units recovered\n" );
													fprintf( stderr, "\n" );
												}
break;
case 5:
#line 67 "ReactionScan.y"
	{ 
													KillItem( gReactionList, FreeOneReaction, ListIterator( gReactionList, FindReaction, yystack.l_mark[-3].typreaction->label ) );
													if ( !gErrorOut ) {
														ScanError( "error in the additional information to the reaction labeled '%s'", yystack.l_mark[-3].typreaction->label , 0 ); 
													}
													gErrorOut = FALSE;
													fprintf( stderr,"# error in additional information of reaction recovered\n");
													fprintf( stderr, "\n" );
												}
break;
case 6:
#line 76 "ReactionScan.y"
	{ 
													LinkPtr link = NULL;

													ClearComposition( gComposition );
													
													KillItem( gReactionList, FreeOneReaction, ListIterator( gReactionList, FindReaction, yystack.l_mark[-2].typreaction->label ) );
													
													if ( !gErrorOut ) {
														ScanError( "error in the definition of the reaction labeled '%s'", yystack.l_mark[-2].typreaction->label , 0 ); 
													}
													gErrorOut = FALSE;
													fprintf( stderr,"# error in reaction recovered\n");
													fprintf( stderr, "\n" );
												}
break;
case 7:
#line 90 "ReactionScan.y"
	{
													if ( !gErrorOut ) {
														ScanError( ( yystack.l_mark[-2].typexponent.paraCoeff ) ? "error in an exponent of the definition of units\n# the last thing i matched is %s" : 
																				  "error in an exponent of the definition of units\n# the last thing i matched is %s%d", 
																 yystack.l_mark[-2].typexponent.parameter , yystack.l_mark[-2].typexponent.number );
													}
													gErrorOut = FALSE;
													fprintf( stderr,"# error in an exponent recovered\n");
													fprintf( stderr, "\n" );
												}
break;
case 8:
#line 100 "ReactionScan.y"
	{
										if ( !gErrorOut ) {
											ScanError( "error in declaration", "", 0 );
										}
										gErrorOut = FALSE;
										fprintf( stderr,"# error in declaration recovered\n");
										fprintf( stderr, "\n" );
		   							}
break;
case 9:
#line 115 "ReactionScan.y"
	{ 
												yyval.typdimension = yystack.l_mark[-1].typdimension; 
												
													if ( yydebug ) {
														fprintf( stderr,"reduce dimension with parenthesis\n"); 
													}
											}
break;
case 10:
#line 122 "ReactionScan.y"
	{ 
												if ( DivideDimension( &yyval.typdimension, yystack.l_mark[-2].typdimension, yystack.l_mark[0].typdimension ) ) {
												  YYERROR;
											  	}
												
													if ( yydebug ) {
														fprintf( stderr,"reduce dimension with '/'\n");
													}
											}
break;
case 11:
#line 131 "ReactionScan.y"
	{ 
												if ( MultiplyDimension( &yyval.typdimension, yystack.l_mark[-2].typdimension, yystack.l_mark[0].typdimension ) ) {
												  YYERROR;
											  	}
												
													if ( yydebug ) {
														fprintf( stderr,"reduce dimension with '*'\n");
												
													}
											}
break;
case 12:
#line 141 "ReactionScan.y"
	{ 
													if ( yydebug ) {
														fprintf( stderr,"vorher\n");
														fprintf( stderr, "$2.value = %g\n", yystack.l_mark[0].typdimension.value );
														fprintf( stderr,"factor\n");
														fprintf( stderr, "$1->value = %g\n", yystack.l_mark[-1].typdimensionptr->value );
													}
												yystack.l_mark[0].typdimension.value *= yystack.l_mark[-1].typdimensionptr->value;
												CopyDimension( &yyval.typdimension, yystack.l_mark[0].typdimension );
												
													if ( yydebug ) {
														fprintf( stderr,"factor dimension\n");
														fprintf( stderr, "$$.value = %g\n", yyval.typdimension.value );
													}
													
											}
break;
case 13:
#line 157 "ReactionScan.y"
	{ 
												if ( PowerDimension( &yyval.typdimension, yystack.l_mark[-2].typdimension, yystack.l_mark[0].typexponent ) ) {
												  YYERROR;
											  	}
												
													if ( yydebug ) {
														fprintf( stderr,"dimension '^' exp\n");
														fprintf( stderr, "$3.parameter = %s\n", yystack.l_mark[0].typexponent.parameter );
														fprintf( stderr, "$3.paraCoeff = %d\n", yystack.l_mark[0].typexponent.paraCoeff );
														fprintf( stderr, "$3.number = %d\n", yystack.l_mark[0].typexponent.number );
													}
											}
break;
case 14:
#line 169 "ReactionScan.y"
	{ 
												if ( ShiftLetters( &yyval.typdimension, yystack.l_mark[0].typstring ) ) {
												  YYERROR;
											  	}
												
													if ( yydebug ) {
														fprintf( stderr,"LETTERS in dimension\n");
											  			fprintf( stderr, "$1 = %s\n", yystack.l_mark[0].typstring );			
													}
											}
break;
case 15:
#line 181 "ReactionScan.y"
	{	
												if ( MuliplyExponent( &yyval.typexponent, yystack.l_mark[-2].typexponent, yystack.l_mark[0].typexponent ) ) {
												  YYERROR;
											  	}
												
													if ( yydebug ) {
														fprintf( stderr,"reduce exp with '*'\n");
														fprintf( stderr, "$1.paraCoeff = %d\n", yystack.l_mark[-2].typexponent.paraCoeff );
														fprintf( stderr, "$3.paraCoeff = %d\n", yystack.l_mark[0].typexponent.paraCoeff );
														fprintf( stderr, "$1.parameter = %s\n", yystack.l_mark[-2].typexponent.parameter );
														fprintf( stderr, "$3.parameter = %s\n", yystack.l_mark[0].typexponent.parameter );
													}
											}
break;
case 16:
#line 195 "ReactionScan.y"
	{	
													if ( MuliplyExponent( &yyval.typexponent, yystack.l_mark[-1].typexponent, yystack.l_mark[0].typexponent ) ) {
													  YYERROR;
													}
														if ( yydebug ) {
															fprintf( stderr,"reduce exp with '*'\n");
															fprintf( stderr, "$1.paraCoeff = %d\n", yystack.l_mark[-1].typexponent.paraCoeff );
															fprintf( stderr, "$2.paraCoeff = %d\n", yystack.l_mark[0].typexponent.paraCoeff );
															fprintf( stderr, "$1.parameter = %s\n", yystack.l_mark[-1].typexponent.parameter );
															fprintf( stderr, "$2.parameter = %s\n", yystack.l_mark[0].typexponent.parameter );
														}
												}
break;
case 17:
#line 207 "ReactionScan.y"
	{	
													if ( MuliplyExponent( &yyval.typexponent, yystack.l_mark[-3].typexponent, yystack.l_mark[-1].typexponent ) ) {
													  YYERROR;
													}
														if ( yydebug ) {
															fprintf( stderr,"reduce exp with '*'\n");
															fprintf( stderr, "$1.paraCoeff = %d\n", yystack.l_mark[-3].typexponent.paraCoeff );
															fprintf( stderr, "$3.paraCoeff = %d\n", yystack.l_mark[-1].typexponent.paraCoeff );
															fprintf( stderr, "$1.parameter = %s\n", yystack.l_mark[-3].typexponent.parameter );
															fprintf( stderr, "$3.parameter = %s\n", yystack.l_mark[-1].typexponent.parameter );
														}
												}
break;
case 18:
#line 219 "ReactionScan.y"
	{	
												if ( AddExponent( &yyval.typexponent, yystack.l_mark[-2].typexponent, yystack.l_mark[0].typexponent ) ) {
												  YYERROR;
											  	}
												
													if ( yydebug ) {
														fprintf( stderr,"reduce exp with '+'\n");
													}
											}
break;
case 19:
#line 229 "ReactionScan.y"
	{	
												if ( DivideExponent( &yyval.typexponent, yystack.l_mark[-2].typexponent, yystack.l_mark[0].typexponent ) ) {
												  YYERROR;
											  	}
												
													if ( yydebug ) {
														fprintf( stderr,"reduce exp with '/'\n");
													}
											}
break;
case 20:
#line 239 "ReactionScan.y"
	{	
												if ( SubtractExponent( &yyval.typexponent, yystack.l_mark[-2].typexponent, yystack.l_mark[0].typexponent ) ) {
												  YYERROR;
											  	}
												
													if ( yydebug ) {
														fprintf( stderr,"reduce exp with '-'\n");
													}
											}
break;
case 21:
#line 249 "ReactionScan.y"
	{	
												yyval.typexponent = yystack.l_mark[-1].typexponent;
													if ( yydebug ) {
														fprintf( stderr,"reduce exp with parenthesis\n");
														fprintf( stderr, "$$.number = %d\n", yyval.typexponent.number );
														fprintf( stderr, "$$.parameter = %s\n", yyval.typexponent.parameter );
														fprintf( stderr, "$$.paraCoeff = %d\n", yyval.typexponent.paraCoeff );
													}
											}
break;
case 22:
#line 259 "ReactionScan.y"
	{	
												yyval.typexponent.number = -yystack.l_mark[0].typexponent.number;
												yyval.typexponent.paraCoeff = -yystack.l_mark[0].typexponent.paraCoeff;
												strcpy( yyval.typexponent.parameter, yystack.l_mark[0].typexponent.parameter );
													if ( yydebug ) {
														fprintf( stderr,"reduce exp with unary '-'\n");
														fprintf( stderr, "$$.parameter = %s\n", yyval.typexponent.parameter );
														fprintf( stderr, "$$.number = %d\n", yyval.typexponent.number );
														fprintf( stderr, "$$.paraCoeff = %d\n", yyval.typexponent.paraCoeff );
													}
											}
break;
case 23:
#line 270 "ReactionScan.y"
	{	
												ShiftExponent( &yyval.typexponent, yystack.l_mark[0].typexponent );
												
													if ( yydebug ) {
														fprintf( stderr,"E_ALPHANUM in exp\n");
														fprintf( stderr, "$1.number = %d\n", yystack.l_mark[0].typexponent.number );
														fprintf( stderr, "$1.parameter = %s\n", yystack.l_mark[0].typexponent.parameter );
														fprintf( stderr, "$1.paraCoeff = %d\n", yystack.l_mark[0].typexponent.paraCoeff );
													}
											}
break;
case 24:
#line 284 "ReactionScan.y"
	{ yyval.typreaction = yystack.l_mark[0].typreaction; }
break;
case 25:
#line 285 "ReactionScan.y"
	{ PutSpeciesToReaction( yyval.typreaction, yystack.l_mark[0].typint, yystack.l_mark[-1].typdouble ); 
													if ( yydebug ) {
											  			fprintf( stderr,"SPECIESCOEFF\n"); 
													}
											}
break;
case 26:
#line 290 "ReactionScan.y"
	{ PutSpeciesToReaction( yyval.typreaction, yystack.l_mark[0].typint, 1.0 ); 
													if ( yydebug ) {
											  			fprintf( stderr,"SPECIES\n");
													}
											}
break;
case 27:
#line 295 "ReactionScan.y"
	{ yyval.typreaction->withThirdBody = TRUE;
											  yyval.typreaction->thirdBodyNumber = yystack.l_mark[0].typint; }
break;
case 28:
#line 297 "ReactionScan.y"
	{ yyval.typreaction->partialEquilibrium = TRUE;
											  ChangeSignsOfReactionCoeffs( yyval.typreaction, RIGHT ); }
break;
case 29:
#line 299 "ReactionScan.y"
	{ ChangeSignsOfReactionCoeffs( yyval.typreaction, RIGHT ); }
break;
case 30:
#line 300 "ReactionScan.y"
	{ ChangeSignsOfReactionCoeffs( yyval.typreaction, LEFT ); }
break;
case 31:
#line 301 "ReactionScan.y"
	{ ChangeSignsOfReactionCoeffs( yyval.typreaction, RIGHT ); 
												if ( yystack.l_mark[0].typint ) {
													if ( CheckStoichiometry( yyval.typreaction ) ) { 
															gLineContents[strlen(gLineContents)-1] = '\0';
															ScanError( "error in stoichiometry of the reaction labeled '%s'", yyval.typreaction->label, 0 ); 
															WriteCompoOfReaction( yyval.typreaction );
															YYERROR;
													}
												}
											}
break;
case 32:
#line 314 "ReactionScan.y"
	{ *yystack.l_mark[-1].ptrdouble = yystack.l_mark[0].typdouble; 
													if ( yydebug ) {
											  			fprintf( stderr, "MI_LITERAL\n"); 
													}
											}
break;
case 33:
#line 319 "ReactionScan.y"
	{ *yystack.l_mark[-1].ptrdouble = yystack.l_mark[0].typdouble; 
													if ( yydebug ) {
											  			fprintf( stderr, "MI_LITERAL\n"); 
													}
											}
break;
#line 836 "y.tab.c"
    }
    yystack.s_mark -= yym;
    yystate = *yystack.s_mark;
    yystack.l_mark -= yym;
    yym = yylhs[yyn];
    if (yystate == 0 && yym == 0)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
        yystate = YYFINAL;
        *++yystack.s_mark = YYFINAL;
        *++yystack.l_mark = yyval;
        if (yychar < 0)
        {
            if ((yychar = YYLEX) < 0) yychar = 0;
#if YYDEBUG
            if (yydebug)
            {
                yys = 0;
                if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
                if (!yys) yys = "illegal-symbol";
                printf("%sdebug: state %d, reading %d (%s)\n",
                        YYPREFIX, YYFINAL, yychar, yys);
            }
#endif
        }
        if (yychar == 0) goto yyaccept;
        goto yyloop;
    }
    if ((yyn = yygindex[yym]) && (yyn += yystate) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yystate)
        yystate = yytable[yyn];
    else
        yystate = yydgoto[yym];
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *yystack.s_mark, yystate);
#endif
    if (yystack.s_mark >= yystack.s_last && yygrowstack(&yystack))
    {
        goto yyoverflow;
    }
    *++yystack.s_mark = (short) yystate;
    *++yystack.l_mark = yyval;
    goto yyloop;

yyoverflow:
    yyerror("yacc stack overflow");

yyabort:
    yyfreestack(&yystack);
    return (1);

yyaccept:
    yyfreestack(&yystack);
    return (0);
}
