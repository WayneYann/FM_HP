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
extern YYSTYPE yylval;
