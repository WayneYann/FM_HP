%{  
/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 */

#include "ReactionParser.h"
#include"ReactionScan.tab.h"
#undef YY_USER_ACTION
#define YY_USER_ACTION LineCounter();
%}

ISOMERE        		[A-Za-z0-9\-\_\*\,\(\)]*
ID          		[A-Za-z0-9]*
LETTER				[A-Za-z]
LABEL				[a-zA-Z0-9\. ]+
NUMBER				[0-9\.]+
INT					[-]?[0-9]+
FLOAT1				[-]?[0-9]+"."[0-9]*
FLOAT2				[-]?"."[0-9]+
EXP					[eE][-+]?[0-9]+
WHITESPACE			[ \t\n]+
WHITESPACE5			[ \t\n]5

%x scSaveComment
				/* is called from INITIAL */
%x comment				
				/* is called from */
%x scSpecies			
				/* is called from scCoeff */
%x sclet				
				/* is called from INITIAL */
%x allowatoms			
				/* is called from sclet   */
%x scUnits				
				/* is called from preUnits   */
%x scCoeff				
				/* is called from INITIAL */
%x scMoreInfo			
				/* is called from scCoeff */
%x scFc					
				/* is called from scMoreInfo */
%x scMdefinition		
				/* is called from sctbcoeff */
%x sctbcoeff			
				/* is called from sclet */
%x scAdditionalSpecies
				/* is called from sclet */
%x scSteadyStates
				/* is called from sclet */
%x scSkipM				
				/* is called from scMdefinition */
%x scSkipSpecies		
				/* is called from scSpecies, INITIAL */
%x scSkipLet			
				/* is called from allowatoms, scUnits, sclet */
%x preUnits				
				/* is called from sclet */
%x scexponent			
				/* is called from scUnits */
%x scOrderOfReaction
				/* is called from sclet */
%x scTemperatureExponent
				/* is called from sclet */
%x scCComment
				/* is called from INITIAL */
%x scPAHSymbol
				/* is called from sclet */
%x scSootSymbol
				/* is called from sclet */

%%

		static char 		species[256];
		static char			label[256];
		static ReactionPtr 	reaction = NULL;
		static Flag			openParent = FALSE;
		static Flag			errorRecovery = FALSE;		
		static LinkPtr		commentFound = NULL;		
		static int			ERROR = FIRSTTOKEN - 2;
		char				thirdBodyName[32];
		char 				*lastAtom = NULL;
		Double				coefficient = 0.0;

%{
/* eat up whitespace  */
%}
<scCoeff,sclet,allowatoms,scUnits,INITIAL>[ \t\n]+ 
<scTemperatureExponent,scPAHSymbol,scSootSymbol,scSteadyStates,scOrderOfReaction,scexponent>[ \t\n]+ 

%{
/* search initial condition  */
%}
<INITIAL>"/*"	{
										BEGIN(scCComment);
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr,"BEGIN(scCComment)\n");
											}
									}
<INITIAL>#	{
										BEGIN(comment);
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr,"BEGIN(comment)\n");
											}
									}
<INITIAL>#P	{
										BEGIN(scSaveComment);
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr,"BEGIN(scSaveComment)\n");
											}
									}
<sclet,allowatoms,scUnits,INITIAL>#	{
										BEGIN(comment);
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr,"BEGIN(comment)\n");
											}
									}
<scCoeff,scFc,scMoreInfo,scSpecies>#	{
											BEGIN(comment);
												if (gOptions->debugScanner == TRUE) {
													fprintf( stderr,"BEGIN(comment)\n");
												}
										}
<scMdefinition,sctbcoeff>#	{
								BEGIN(comment);
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(comment)\n" );
									}
							}

<scCComment>"*"+"/"	{
						BEGIN(INITIAL);
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "BEGIN(INITIAL)\n" );
							}
					}
<scCComment>"*"+[^*/]*	{
							/* eat up '*'s not followed by '/' */
						}
<scCComment>[^*]*	{
						/* eat up anything that's not a '*' */
					}
					
					
[Ll][Ee][Tt]		{	
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "BEGIN(sclet)\n" );
							}
						BEGIN(sclet);
					}
{LABEL}/:          	{	char	*followLabel = NULL;
						int		lenOfPAHSymbol = strlen( gPAHReactionList->listName );
						int		lenOfSootSymbol = strlen( gSootReactionList->listName );
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "BEGIN(scCoeff)\n" );
							}
						BEGIN(scCoeff);
						/*gLineContents[strlen( gLineContents ) - 1] = '\0';*/
						if ( yyleng > 255 || strstr( ( char * )yytext, "{" ) ) {
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "BEGIN(scSkipSpecies)\n" );
								}
							BEGIN(scSkipSpecies);
							ScanError( "i think you've forgotten to limit the label with a colon", "", 0 );
							strcpy( label, "error" );
						}
						else {
							strcpy( label, ( char * )yytext );
							SkipWhites( label );
							if ( gOptions->useForwardBackward ) {
								if ( label[strlen(label)-1] == 'F' || label[strlen(label)-1] == 'B' ) {
									label[strlen(label)-1] = tolower( label[strlen(label)-1] );
								}
							}
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "label found: '%s'\n", label);
								}
							if ( ( ListIterator( gReactionList, FindReaction, label ) ) ) {
								ScanError( "the label '%s' has already been used", label, 0 );
								strcpy( label, "error" );
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scSkipSpecies)\n" );
									}
								BEGIN(scSkipSpecies);
							}
							
							if ( commentFound ) {
								followLabel = ( char * ) malloc( sizeof( label ) + 1 );
								strcpy( followLabel, label );
								free( ( char * ) commentFound->item );
								commentFound->item = followLabel;
								commentFound = NULL;
							}
						}
						if ( lenOfPAHSymbol && strncmp( label, gPAHReactionList->listName, lenOfPAHSymbol ) == 0 ) {
							reaction = AddReaction( gPAHReactionList, label );
						}
						else if ( lenOfSootSymbol && strncmp( label, gSootReactionList->listName, lenOfSootSymbol ) == 0 ) {
							reaction = AddReaction( gSootReactionList, label );
						}
						else {
							reaction = AddReaction( gReactionList, label );
						}
						
						gtypeOfYylval = TYPREACTION;
						yylval.typreaction = reaction;
						return LABEL;
					}

%{
/* eat up one-line comments */
%}
<comment>\n     	{ 
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "BEGIN(INITIAL)\n" );
							}
						BEGIN(INITIAL);
                	}
<comment>[^\n]  
%{
/* save one-line comments */
%}
<scSaveComment>\n     	{ 
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "BEGIN(INITIAL)\n" );
							}
						BEGIN(INITIAL);
                	}
<scSaveComment>[^\n]+	{
							AddString( gCommentList, ( char * )yytext );
							AddString( gLabelFollowingCommentList, "" );
							commentFound = gLabelFollowingCommentList->lastItem;
						}  
%{
/* species coefficient */
%}
<scCoeff>"{"					{	BEGIN(scMoreInfo);
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(scMoreInfo)\n" );
											fprintf( stderr, "gLineContents = %s\n", gLineContents );
										}
									gtypeOfYylval = TYPINT;
									yylval.typint = TRUE; /* check stoichiometry */
									return LEFTBRACE;
								}
<scCoeff>[ \t\n]*[Oo][Rr][Dd][Ee][Rr][ \t\n]*=[ \t\n]*[0-9]+[ \t\n]*"{"		{
							char		*tmp;
							BEGIN(scMoreInfo);
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "BEGIN(scMoreInfo)\n" );
									fprintf( stderr, "gLineContents = %s\n", gLineContents );
								}
							CopyString( species, ( char * )yytext );
							SkipWhites( species );
							UpperString( species );
							if ( ( tmp = strstr( species, "ORDER=" ) ) != NULL ) {
								species[strlen(species)-1] = '\0';
								reaction->orderOfReaction = atoi( &tmp[6] );
									if (gOptions->debugScanner == TRUE) {
										fprintf(stderr, "reaction order %d found\n", reaction->orderOfReaction );
									}
							}
							gtypeOfYylval = TYPINT;
							yylval.typint = TRUE; /* check stoichiometry */
							return LEFTBRACE;
						}
<scCoeff>[ \t\n]*[Nn][Oo][Cc][Hh][Ee][Cc][Kk][ \t\n]*"{"					{	BEGIN(scMoreInfo);
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(scMoreInfo)\n" );
											fprintf( stderr, "gLineContents = %s\n", gLineContents );
										}
									gtypeOfYylval = TYPINT;
									yylval.typint = FALSE; /* don't check stoichiometry */
									return LEFTBRACE;
								}
	/*<scCoeff>({INT}|{FLOAT1}|{FLOAT2}|{FLOAT1}{EXP}|{FLOAT2}{EXP}|{INT}{EXP})\-		{*/
<scCoeff>({ISOMERE})\-		{
										if (gOptions->debugScanner == TRUE) {
											fprintf(stderr, "i've found an isomere introduced by %s \n", ( char * )yytext );
											fprintf( stderr, "BEGIN(scSpecies)\n" );
										}
									BEGIN(scSpecies);
									gLineContents[strlen( gLineContents ) - yyleng] = '\0';
									yyless( 0 );
								}

<scCoeff>{INT}|{FLOAT1}|{FLOAT2}|{FLOAT1}{EXP}|{FLOAT2}{EXP}|{INT}{EXP}		{
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(scSpecies)\n" );
										}
									BEGIN(scSpecies);
									coefficient = atof( ( char * )yytext );
									gtypeOfYylval = TYPDOUBLE;
									yylval.typdouble = coefficient;
										if (gOptions->debugScanner == TRUE) {
											fprintf(stderr, "i've found the coefficient %s \n", ( char * )yytext );
										}
									return SPECIESCOEFF;
								}
<scCoeff>{LETTER}				{	
										if (gOptions->debugScanner == TRUE) {
											fprintf(stderr, "yytext = '%s'\n", yytext );
											fprintf( stderr, "BEGIN(scSpecies)\n" );
										}
									BEGIN(scSpecies);
									gLineContents[strlen( gLineContents ) - 1] = '\0';
									yyless( 0 );
								}
<scCoeff>":"
<scCoeff>"=>"|"->"				{	
									return RIGHTARROW;
								}
<scCoeff>"=="					{	
									return EQUAL;
								}
<scCoeff>"<="|"<-"				{	
									return LEFTARROW;
								}
<scCoeff>"+"					{	
								}


%{
/*  handle species  */
%}
	/*  isomeres  */
<scSpecies>({ISOMERE})\-/[^>]	{	
										gLineContents[strlen( gLineContents ) - 1] = '\0';
										AdjustNumberOfLines();
										yymore();
											if (gOptions->debugScanner == TRUE) {
												fprintf(stderr, "i've found an isomere introduced by '%s'\n", ( char * )yytext );
											}
									}
	/*  third bodies  */
<scSpecies>[Mm](("'"*)|([0-9]+))		{	
								LinkPtr 	link = NULL;
								char 		*usedThirdBody = NULL;
								
								CopyString( species, ( char * )yytext );
								SkipWhites( species );
								if ( isdigit( yytext[yyleng-1] ) ) {
									UpperString( species );
								}
								else {
									sprintf( species, "M%d", (int)strlen( species ) - 1);
								}
								if ( !( link = ListIterator( gUsedThirdBodyList, FindString, species ) ) ) {
									usedThirdBody = AddString( gUsedThirdBodyList, species );
								}
								else {
									usedThirdBody = ( char * )link->item;
								}
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scCoeff)\n" );
									}
								BEGIN(scCoeff);
								
								gtypeOfYylval = TYPINT;
								yylval.typint = NumberOfListItem( gUsedThirdBodyList, FindString, usedThirdBody )-1;
								return THIRDBODY;
							}
	/*  species  */
<scSpecies>{LETTER}/{LETTER}	{	lastAtom = MatchAtom( ( char * )yytext, yyleng, 2 );
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, " lastAtom  = %s\n",  lastAtom  );
											fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
										}
									if ( lastAtom == 0 ) {
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "BEGIN(scSkipSpecies)\n" );
											}
										BEGIN(scSkipSpecies);
									}	
									gLineContents[strlen( gLineContents ) - 1] = '\0';
									AdjustNumberOfLines();
									yymore();
								}
<scSpecies>{LETTER}			{	lastAtom = MatchAtom( ( char * )yytext, yyleng, 1 );
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, " lastAtom  = %s\n",  lastAtom  );
										fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
									}
								if ( lastAtom == 0 ) {
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(scSkipSpecies)\n" );
										}
									BEGIN(scSkipSpecies);
								}	
								AdjustNumberOfLines();
								yymore();
							}
<scSpecies>[0-9]+			{	int i, num;
								AtomsPtr atom = NULL;
								LinkPtr link = NULL;
								
								for ( i = yyleng; isdigit( yytext[i-1] ); --i ) {/* null statement */;}
								num = atoi( ( char * )&yytext[i] );
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, " lastAtom  = %s\n",  lastAtom  );
									}
								link = ListIterator( gAtomsList, FindAtom, lastAtom );
								atom = ( AtomsPtr )link->item;
								gComposition->vec[atom->number] += num - 1;
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "i've found the number %d following the atom %s\n", num, lastAtom );
									}
								AdjustNumberOfLines();
								yymore();
							}
<scSpecies>{WHITESPACE5}		{	AdjustNumberOfLines();
								yymore();
							}
<scSpecies>{WHITESPACE}		{	AdjustNumberOfLines();
								yymore();
							}
<scSpecies>([Oo][Rr][Dd][Ee][Rr])|([Nn][Oo][Cc][Hh][Ee][Cc][Kk])|([\+-=<{])	{	
						LinkPtr 	link = NULL;
						SpeciesPtr	speciesPtr = NULL;
						char		*tmp;
						
						CopyString( species, ( char * )yytext );
						SkipWhites( species );
						UpperString( species );
							if (gOptions->debugScanner == TRUE) {
								fprintf(stderr, "yytext in ([\\+-=<{]) = '%s'\n", ( char * )yytext );
							}
						if ( strstr( species, "NOCHECK" ) != NULL ) {
								if (gOptions->debugScanner == TRUE) {
									fprintf(stderr, "noCheck found\n" );
								}
							species[ strlen(species)-7 ] = '\0';
							gLineContents[strlen( gLineContents ) - 7] = '\0';
							yyless( yyleng-7 );
						}
						else if ( ( tmp = strstr( species, "ORDER" ) ) != NULL ) {
								if (gOptions->debugScanner == TRUE) {
									fprintf(stderr, "reaction order found\n" );
								}
							species[ strlen(species)-5 ] = '\0';
							gLineContents[strlen( gLineContents ) - 5] = '\0';
							yyless( yyleng-5 );
						}
						else {
							species[ strlen(species)-1 ] = '\0';
							gLineContents[strlen( gLineContents ) - 1] = '\0';
							unput( yytext[yyleng-1] );
						}
							if (gOptions->debugScanner == TRUE) {
								fprintf(stderr, "Species of length %d found: '%s'\n", (int)strlen(species), species );
							}
						if ( !( link = ListIterator( gSpeciesList, gSpeciesCompareFunc, species ) ) ) {
							int lenPAH = strlen( gPAHReactionList->listName );
							int lenSoot = strlen( gSootReactionList->listName );
							speciesPtr = AddSpecies( gSpeciesList, species, gComposition );
							if ( lenPAH && strncmp( label, gPAHReactionList->listName, lenPAH ) == 0 ) {
								/* species of polymerisation reactions are steady state */
								speciesPtr->isSteadyState = TRUE;
							}
							else if ( lenSoot && strncmp( label, gSootReactionList->listName, lenSoot ) == 0 ) {
								/* species of soot reactions are steady state */
								speciesPtr->isSteadyState = TRUE;
							}
							else if ( ListIterator( gSteadyStateList, FindString, species ) ) {
								/* explicitly declared steady state species */
								speciesPtr->isSteadyState = TRUE;
							}
							gComposition = NewIntVector( gAllowAtomList->items );
						}
						else {
							speciesPtr = ( SpeciesPtr )link->item;
							ClearComposition( gComposition );
						}
							if (gOptions->debugScanner == TRUE) {
								fprintf(stderr, "Species of length %d found: '%s'\n", (int)strlen(species), species );
								PrintSpeciesComposition( speciesPtr, stderr );
								fprintf( stderr, "BEGIN(scCoeff)\n" );
							}
						BEGIN(scCoeff);

						gtypeOfYylval = TYPINT;
						yylval.typint = speciesPtr->number;
						return SPECIES;
					}

%{
/*  more Information about the Reaction  */
%}
<scMoreInfo>"}"		{	
						int		len = 0;
						LinkPtr link = NULL;
						/*char	*contents = NULL;*/
						
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "BEGIN(INITIAL)\n" );
							}
						BEGIN(INITIAL);
						return '}';
					}
<scMoreInfo>[Ee]\_?[Ii][Nn]?[Ff]?	{	
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "E_inf follows \n" );
								fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
							}
						reaction->withLindemann = TRUE;
						gtypeOfYylval = PTRDOUBLE;
						yylval.ptrdouble = &reaction->eInf;
						return MI_LITERAL;
					}
<scMoreInfo>([Aa]|[Bb])\_?[Ii][Nn]?[Ff]?	{	
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "A_inf follows \n" );
								fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
							}
						reaction->withLindemann = TRUE;
						gtypeOfYylval = PTRDOUBLE;
						yylval.ptrdouble = &reaction->aInf;
						return MI_LITERAL;
					}
<scMoreInfo>[Nn]\_?[Ii][Nn]?[Ff]?	{	
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "n_inf follows \n" );
								fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
							}
						reaction->withLindemann = TRUE;
						gtypeOfYylval = PTRDOUBLE;
						yylval.ptrdouble = &reaction->nInf;
						return MI_LITERAL;
					}
<scMoreInfo>[Ee]	{	
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "E follows \n" );
								fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
							}
						gtypeOfYylval = PTRDOUBLE;
						yylval.ptrdouble = &reaction->e;
						return MI_LITERAL;
					}
<scMoreInfo>[Aa]|[Bb]	{	
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "A follows \n" );
								fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
							}
						gtypeOfYylval = PTRDOUBLE;
						yylval.ptrdouble = &reaction->a;
						return MI_LITERAL;
					}
<scMoreInfo>[Nn]	{	
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "n follows \n" );
								fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
							}
						gtypeOfYylval = PTRDOUBLE;
						yylval.ptrdouble = &reaction->n;
						return MI_LITERAL;
					}
<scMoreInfo>[Ff][Cc][AaBbCc]	{	
						int		len;
						char  	name[32];
						LinkPtr	namelink = NULL;
						LinkPtr	contentsLink = NULL;
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "fcCoeff follows \n" );
								fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
							}
						gtypeOfYylval = PTRDOUBLE;
						switch ( yytext[2] ) {
							case 'a':
							case 'A':
								yylval.ptrdouble = &reaction->fca;
								break;
							case 'b':
							case 'B':
								yylval.ptrdouble = &reaction->fcb;
								break;
							case 'c':
							case 'C':
								yylval.ptrdouble = &reaction->fcc;
								break;
							default:
								fprintf( stderr, "that's impossible\n" );
						}

						len = strlen( label );
						strcpy( name, label );
						if ( reaction->hasForward || reaction->hasBackward ) {
							name[len-1] = '\0';
						}
						if ( !( namelink = ListIterator( gFcNameList, FindString, name ) ) ) {
							AddString( gFcNameList, name );
							AddString( gFcContentsList, "" );
							AddInt( gFcLineNumberList, gNumberOfLines );
							reaction->lindemannNumber = NumberOfListItem( gFcNameList, FindString, name ) - 1;						
						}
	
						return MI_LITERAL;
					}
<scMoreInfo>[Ff][Cc][Tt][AaBbCc]	{	
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "fcCoeff follows \n" );
								fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
							}
						gtypeOfYylval = PTRDOUBLE;
						switch ( yytext[3] ) {
							case 'a':
							case 'A':
								yylval.ptrdouble = &reaction->fcTa;
								break;
							case 'b':
							case 'B':
								yylval.ptrdouble = &reaction->fcTb;
								break;
							case 'c':
							case 'C':
								yylval.ptrdouble = &reaction->fcTc;
								break;
							default:
								fprintf( stderr, "that's impossible\n" );
						}
						return MI_LITERAL;
					}
<scMoreInfo>{INT}|{FLOAT1}|{FLOAT2}|{FLOAT1}{EXP}|{FLOAT2}{EXP}|{INT}{EXP}			{	
						Double num;
						num = atof( ( char * )yytext );
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, " num  = %g\n",  num  );
							}
						gtypeOfYylval = TYPDOUBLE;
						yylval.typdouble = num;
						return MI_NUMBER;
					}
<scMoreInfo>[,;.:= \t\n]			/* skip */
<scMoreInfo>[Ff][Cc]{WHITESPACE}*={WHITESPACE}*	{	
														if (gOptions->debugScanner == TRUE) {
															fprintf( stderr, "Fc follows \n" );
															fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
															fprintf( stderr, "BEGIN(scFc)\n" );
														}
													BEGIN(scFc);
												}

 /*  handle broadening factor  */
<scFc>;				{	
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "BEGIN(scMoreInfo)\n" );
							}
						BEGIN(scMoreInfo);
					}
<scFc>[^;]+			{	
						int		len;
						char  	name[32];
						char	*broadening = NULL;
						String	buffer = NULL;
						String	otherBuffer = NULL;
						LinkPtr	namelink = NULL;
						LinkPtr	contentsLink = NULL;
						int		number;
						int		lindemannNumber;
						

						broadening = ( char * ) malloc( strlen( ( char * )yytext ) + 1 );
						strcpy( broadening, ( char * )yytext );
						len = strlen( label );
						strcpy( name, label );
						if ( reaction->hasForward || reaction->hasBackward ) {
							name[len-1] = '\0';
						}
						if ( !( namelink = ListIterator( gFcNameList, FindString, name ) ) ) {
							AddString( gFcNameList, name );
							AddString( gFcContentsList, broadening );
							AddInt( gFcLineNumberList, gNumberOfLines );
							lindemannNumber = NumberOfListItem( gFcNameList, FindString, name ) - 1;						
							reaction->lindemannNumber = lindemannNumber;
						}
						else {
							buffer = ( String ) namelink->item;
							number = NumberOfListItem( gFcNameList, FindString, buffer );
							contentsLink =  Find_n_th_Item( number, gFcContentsList );
							buffer = ( String ) contentsLink->item;
							otherBuffer = ( String ) malloc( strlen( buffer ) + 1 );
							strcpy( otherBuffer, buffer );
							SkipWhites( otherBuffer );
							if ( strcmp( otherBuffer, "" ) == 0 ) {
								free( ( String ) contentsLink->item );
								contentsLink->item = ( String ) malloc( strlen( broadening ) + 1 );
								strcpy( ( String )contentsLink->item, broadening );
							}
							reaction->lindemannNumber = number - 1;
							free( otherBuffer );
						}
							if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
							}
						free( broadening );
					}					

%{
/*  let  */
%}
<sclet>"."				{
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "BEGIN(INITIAL)\n" );
								}
							BEGIN(INITIAL);
						}

<sclet>[ \t\n]*[Mm][Ee][Cc][Hh][Aa][Nn][Ii][Ss][Mm][ \t\n]*[Bb][Ee][ \t\n]*[Gg][Ll][Oo][Bb][Aa][Ll]\.  {
																	gHeader->globalMechanism = TRUE;
																	BEGIN(INITIAL);
																}
<sclet>[Aa][Ll][Ll][Oo][Ww][Ee][Dd][ \t\n]*[Aa][Tt][Oo][Mm][Ss]	{
																		if (gOptions->debugScanner == TRUE) {
																			fprintf( stderr, "BEGIN(allowatoms)\n" );
																		}
																	BEGIN(allowatoms);
																}
<sclet>[Uu][Nn][Ii][Tt][Ss]?	{
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(preUnits)\n" );
										}
									BEGIN(preUnits);
								}
<sclet>\[?[ \t\n]*[Mm]	{
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "BEGIN(sctbcoeff)\n" );
								}
							coefficient = 1.0;
							BEGIN(sctbcoeff);
							gLineContents[strlen( gLineContents ) - 1] = '\0';
							unput( yytext[yyleng-1] );
						}
<sclet>[Aa][Dd][Dd][Ii][Tt][Ii][Oo][Nn][Aa][Ll][ \t\n]*[Ss][Pp][Ee][Cc][Ii][Ee][Ss][ \t\n]*	{
												if (gOptions->debugScanner == TRUE) {
													fprintf( stderr, "BEGIN(scAdditionalSpecies)\n" );
												}
											BEGIN(scAdditionalSpecies);
										}

<sclet>[Oo][Rr][Dd][Ee][Rr][ \t\n]*[Oo][Ff][ \t\n]*[Rr][Ee][Aa][Cc][Tt][Ii][Oo][Nn]	{
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scOrderOfReaction)\n" );
									}
								BEGIN(scOrderOfReaction);
							}
<sclet>[Tt][Ee][Mm][Pp][Ee][Rr][Aa][Tt][Uu][Rr][Ee][ \t\n]*[Ee][Xx][Pp][Oo][Nn][Ee][Nn][Tt]	{
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scTemperatureExponent)\n" );
									}
								BEGIN(scTemperatureExponent);
							}
<sclet>[Ss][Tt][Ee][Aa][Dd][Yy][ \t\n]*[Ss][Tt][Aa][Tt][Ee][Ss]	{
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scSteadyStates)\n" );
									}
								BEGIN(scSteadyStates);
							}
<sclet>[Ss][Yy][Mm][Bb][Oo][Ll][ \t\n]*[Ff][Oo][Rr][ \t\n]*[Pp][Aa][Hh]	{
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scPAHSymbol)\n" );
									}
								BEGIN(scPAHSymbol);
							}
<sclet>[Ss][Yy][Mm][Bb][Oo][Ll][ \t\n]*[Ff][Oo][Rr][ \t\n]*[Ss][Oo][Oo][Tt]	{
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scSootSymbol)\n" );
									}
								BEGIN(scSootSymbol);
							}

%{
/*   scPAHSymbol   */
%}
<scPAHSymbol>([Bb][Ee])|=		/* skip */
<scPAHSymbol>"."			{
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "BEGIN(INITIAL)\n" );
											}
										BEGIN(INITIAL);
									}
<scPAHSymbol>[a-zA-Z\_]+	{
								LinkPtr link = NULL;
								
								if ( strcmp( gPAHReactionList->listName, "undefined" ) != 0 ) {
									ScanError( "it is not allowed to declare more than one symbol for pah reactions", "", 0 );
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scSkipLet)\n" );
									}
									BEGIN(scSkipLet);
								}
								else {
									if ( yyleng > gHeader->maxLenOfString ) {
										ScanError( "pah symbol too long", "", 0 );
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(scSkipLet)\n" );
										}
										BEGIN(scSkipLet);
									}
									else {
										char	*symbol = NULL;
										
										if ( !( symbol = ( char * )malloc( strlen( ( char * )yytext ) + 1 ) ) ) {
											FATAL( "memory allocation failed" );
										}
										strcpy( symbol, ( char * )yytext );
										free( gPAHReactionList->listName );
										gPAHReactionList->listName = symbol;
									}
								}
							}
%{
/*   scSootSymbol   */
%}
<scSootSymbol>([Bb][Ee])|=		/* skip */
<scSootSymbol>"."			{
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "BEGIN(INITIAL)\n" );
											}
										BEGIN(INITIAL);
									}
<scSootSymbol>[a-zA-Z\_]+	{
								LinkPtr link = NULL;
								
								if ( strcmp( gSootReactionList->listName, "undefined" ) != 0 ) {
									ScanError( "it is not allowed to declare more than one symbol for soot reactions", "", 0 );
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scSkipLet)\n" );
									}
									BEGIN(scSkipLet);
								}
								else {
									if ( yyleng > gHeader->maxLenOfString ) {
										ScanError( "soot symbol too long", "", 0 );
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(scSkipLet)\n" );
										}
										BEGIN(scSkipLet);
									}
									else {
										char	*symbol = NULL;
										
										if ( !( symbol = ( char * )malloc( strlen( ( char * )yytext ) + 1 ) ) ) {
											FATAL( "memory allocation failed" );
										}
										strcpy( symbol, ( char * )yytext );
										free( gSootReactionList->listName );
										gSootReactionList->listName = symbol;
									}
								}
							}

%{
/*   scOrderOfReaction   */
%}
<scOrderOfReaction>([Bb][Ee])|=		/* skip */
<scOrderOfReaction>"."				{
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "BEGIN(INITIAL)\n" );
											}
										BEGIN(INITIAL);
									}
<scOrderOfReaction>[a-zA-Z\_]+		{
										static Flag alreadyDefined = FALSE;
										LinkPtr link = NULL;
										
										if ( alreadyDefined ) {
											ScanError( "it is not allowed to declare more than one parameter for the order of reaction", "", 0 );
												if (gOptions->debugScanner == TRUE) {
													fprintf( stderr, "BEGIN(scSkipLet)\n" );
												}
											BEGIN(scSkipLet);
										}
										else {
											if ( link = ListIterator( gAllowReacExpList, FindString, "n" ) ) {
												KillItem( gAllowReacExpList, FreeOneString, link );
											}
											AddString( gAllowReacExpList, ( char * )yytext );
											alreadyDefined = TRUE;
										}
									}
									
%{
/*   scTemperatureExponent   */
%}
<scTemperatureExponent>([Bb][Ee])|=		/* skip */
<scTemperatureExponent>"."			{
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "BEGIN(INITIAL)\n" );
											}
										BEGIN(INITIAL);
									}
<scTemperatureExponent>[a-zA-Z\_]+	{
										static Flag alreadyDefined = FALSE;
										LinkPtr link = NULL;
										
										if ( alreadyDefined ) {
											ScanError( "it is not allowed to declare more than one parameter for the order of reaction", "", 0 );
												if (gOptions->debugScanner == TRUE) {
													fprintf( stderr, "BEGIN(scSkipLet)\n" );
												}
											BEGIN(scSkipLet);
										}
										else {
											if ( link = ListIterator( gAllowTempExpList, FindString, "n_k" ) ) {
												KillItem( gAllowTempExpList, FreeOneString, link );
											}
											if ( link = ListIterator( gAllowTempExpList, FindString, "nk" ) ) {
												KillItem( gAllowTempExpList, FreeOneString, link );
											}
											AddString( gAllowTempExpList, ( char * )yytext );
										}
									}

%{
/*   scSteadyStates   */
%}
<scSteadyStates>([Bb][Ee])|=		/* skip */
<scSteadyStates>#[^\n]*\n	/* skip */
<scSteadyStates>"."			{
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "BEGIN(INITIAL)\n" );
											}
										BEGIN(INITIAL);
									}
	/*<scSteadyStates>([0-9a-z,A-Z]+\-[a-zA-Z0-9]+)|[a-zA-Z0-9]+	{*/
<scSteadyStates>(({ISOMERE})\-[a-zA-Z0-9]+)|[a-zA-Z0-9]+	{
										LinkPtr link = NULL;
										
										CopyString( species, ( char * ) yytext );
										UpperString( species );
										if ( link = ListIterator( gSteadyStateList, FindString, species ) ) {
											fprintf( stderr, "#warning: multiple appearance of species %s in steady state list\n", ( char * ) yytext );
										}
										else {
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "steady state species %s found\n", species );
											}
											AddString( gSteadyStateList, species );
										}
									}

%{
/*   preUnits   */
%}
<preUnits>[fF][oO][rR]		/*  skip  */
<preUnits>([Bb][Ee])|=		{
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scUnits)\n" );
									}
								BEGIN(scUnits);
							}
<preUnits>[Aa]|[Ee]			{
								gtypeOfYylval = TYPSTRING;
						
								strcpy( yylval.typstring, ( char * )yytext );
								UpperString( yylval.typstring );
								return LETTERS;
							}
<preUnits>{WHITESPACE5}
<preUnits>{WHITESPACE}

%{
/*   scUnits   */
%}
<scUnits>\[					/*  skip  */
<scUnits>\]					/*  skip  */
<scUnits>\.				{	
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "BEGIN(INITIAL)\n" );
								}
							BEGIN(INITIAL);
							return '.';
						}

<scUnits>kg|mole|cal|mol|atm		{
									gtypeOfYylval = TYPSTRING;
									strcpy( yylval.typstring, ( char * )yytext );
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "LETTERS found\n" );
								}
									return LETTERS;
								}
<scUnits>m/[ \t\n\^\/\*\]\)\.]	{
									gLineContents[strlen( gLineContents ) - 1] = '\0';
									gtypeOfYylval = TYPSTRING;
									strcpy( yylval.typstring, ( char * )yytext );
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "LETTERS found\n" );
								}
									return LETTERS;
								}
<scUnits>m				{
							DimensionPtr dimension = NewDimension( "", 0.0, 0, 0, 0, 0, 0 );
							
							if ( ShiftLetters( dimension, "milli" ) ) {
								ScanError( "%s is not an allowed dimension", ( char * )yytext, 0 );
							}
							gtypeOfYylval = TYPDIMENSIONPTR;
							yylval.typdimensionptr = dimension;
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "FACTOR found\n" );
								}
							return FACTOR;
						}
<scUnits>[kMGc]			{
							DimensionPtr dimension = NewDimension( "", 0.0, 0, 0, 0, 0, 0 );
							
							if ( ShiftLetters( dimension, ( char * )yytext ) ) {
								ScanError( "%s is not an allowed dimension", ( char * )yytext, 0 );
							}
							gtypeOfYylval = TYPDIMENSIONPTR;
							yylval.typdimensionptr = dimension;
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "FACTOR found\n" );
								}
							return FACTOR;
						}
<scUnits>[a-bd-jln-zA-FH-LN-Z\_]+		{
									gtypeOfYylval = TYPSTRING;
									CopyString( yylval.typstring, ( char * )yytext );
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "LETTERS found\n" );
								}
									return LETTERS;
								}
<scUnits>[\/\*\(\)]		{
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "/ * ( or ) found\n" );
								}
							return yytext[0];
						}
<scUnits>\^				{
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "BEGIN(scexponent)\n" );
								}
							BEGIN(scexponent);							
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "^ found\n" );
								}
							return '^';
						}

%{
/*   scexponent   */
%}
<scexponent>\(			{
							openParent = TRUE;
							return '(';
						}
<scexponent>[\*\/\+\-]	{
							return yytext[0];
						}
<scexponent>[0-9]+	{
							gtypeOfYylval = TYPEXPONENT;
							yylval.typexponent.number = atoi( ( char * )yytext );
							yylval.typexponent.parameter[0] = '\0';
							yylval.typexponent.paraCoeff = 0;
							if ( !openParent ) {
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scUnits)\n" );
									}
								BEGIN(scUnits);
							}
							return E_ALPHANUM;
						}
<scexponent>[a-zA-Z\_]+	{
								if ( !( ListIterator( gAllowTempExpList, FindString, ( char * )yytext ) ) && 
										!( ListIterator( gAllowReacExpList, FindString, ( char * )yytext ) ) ) {
									ScanError( "error in exponent of the definition of units: parameter %s is not declared", ( char * )yytext, 0 );
								}
								gtypeOfYylval = TYPEXPONENT;
								yylval.typexponent.number = 0;
								CopyString( yylval.typexponent.parameter, ( char * )yytext );
								yylval.typexponent.paraCoeff = 1;
								if ( !openParent ) {
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(scUnits)\n" );
										}
									BEGIN(scUnits);
								}
								return E_ALPHANUM;
						}
<scexponent>\)			{
							openParent = FALSE;
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "BEGIN(scUnits)\n" );
								}
							BEGIN(scUnits);
							return ')';
						}

%{
/*   definition of M   */
%}

<scMdefinition>\*					/* skip */
<scMdefinition>\[					/* skip */

	/*  isomeres  */
<scMdefinition>({ISOMERE})\-	{	
										AdjustNumberOfLines();
										yymore();
											if (gOptions->debugScanner == TRUE) {
												fprintf(stderr, "i've found an isomere introduced by '%s'\n", ( char * )yytext );
											}
									}
<scMdefinition>[\+\.\] \t\n]	{	
						LinkPtr 	link;
						int			speciesNumber;
						SpeciesPtr	speciesPtr = NULL;

						strcpy( species, ( char * )yytext );
						species[ strlen(species)-1 ] = '\0';
						SkipWhites( species );
						UpperString( species );
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
								fprintf(stderr, "Species of length %d found: '%s'\n", (int)strlen(species), species );
							}
						if ( strcmp( species, "" ) == 0 ) {
							AdjustNumberOfLines();
							yymore();
						}
						else {
							gLineContents[strlen( gLineContents ) - 1] = '\0';
							unput( yytext[yyleng-1] );
							if ( link = ListIterator( gSpeciesList, gSpeciesCompareFunc, species ) ) {
								speciesPtr = ( SpeciesPtr )link->item;
								speciesNumber = speciesPtr->number;
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, " thirdBodyName  = %s\n",  thirdBodyName  );
									}
								FillThirdBodyArray( speciesNumber, coefficient );
								ClearComposition( gComposition );
							}
							else {
								ClearComposition( gComposition );
								ScanError( "'%s' is not a valid species", species, 0 );
							}
	
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "BEGIN(sctbcoeff)\n" );
								}
							coefficient = 1.0;
							BEGIN(sctbcoeff);
						}
					}
<scMdefinition>[Oo][Tt][Hh][Ee][Rr][|Ss]?  {
							FillOtherThirdBodyArray( coefficient );
							
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "BEGIN(sctbcoeff)\n" );
								}
							coefficient = 1.0;
							BEGIN(sctbcoeff);
						}
<scMdefinition>{LETTER}/{LETTER}  {	
							lastAtom = MatchAtom( ( char * )yytext, yyleng, 2 );
							if ( lastAtom == 0 ) {
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scSkipM)\n" );
									}
								BEGIN(scSkipM);
							}	
							gLineContents[strlen( gLineContents ) - 1] = '\0';
							AdjustNumberOfLines();
							yymore();
						}
<scMdefinition>{LETTER}  {	
							lastAtom = MatchAtom( ( char * )yytext, yyleng, 1 );
							if ( lastAtom == 0 ) {
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "BEGIN(scSkipM)\n" );
									}
								BEGIN(scSkipM);
							}	
							AdjustNumberOfLines();
							yymore();
						}
<scMdefinition>[0-9]+  {	
							int i, num;
							AtomsPtr atom = NULL;
							LinkPtr link = NULL;
							
							for ( i = yyleng; isdigit( yytext[i-1] ); --i ) {/* null statement */;}
							num = atoi( ( char * )&yytext[i] );
							link = ListIterator( gAtomsList, FindAtom, lastAtom );
							atom = ( AtomsPtr )link->item;
							gComposition->vec[atom->number] += num - 1;
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "&yytext[i] = %s\n", ( char * )&yytext[i] );
									fprintf( stderr, " lastAtom  = %s\n",  lastAtom  );
									fprintf( stderr, "i've found the number %d following the atom %s\n", num, lastAtom );
								}
							AdjustNumberOfLines();
							yymore();
						}

%{
/*  coefficient of third body  */
%}
<sctbcoeff>{WHITESPACE}
<sctbcoeff>=|[Bb][Ee]
<sctbcoeff>\]
<sctbcoeff>\+
<sctbcoeff>[Mm][ \t\n]*(("'"*)|([0-9]+))		  	{
								CopyString( species, ( char * )yytext );
								SkipWhites( species );
								if ( isdigit( yytext[yyleng-1] ) ) {
									UpperString( species );
								}
								else {
									sprintf( species, "M%d", (int)strlen( species ) - 1);
								}
								strcpy( thirdBodyName, species);
								if ( ListIterator( gUsedThirdBodyList, FindString, species ) ) {
									if ( !( ListIterator( gThirdBodyList, FindThirdBody, species ) ) ) {
										AddThirdBody( gThirdBodyList, species );
									}
									else {
										ScanError( "the third body efficiency '%s' is already defined", species, 0 );
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "BEGIN(scSkipM)\n" );
											}
										BEGIN(scSkipM);
									}
								}
								else {
									ScanError( "you don't need to define the third body efficiency '%s'\n, because you don't use it", species, 0 );
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(scSkipM)\n" );
										}
									BEGIN(scSkipM);
								}
								
						}
	/*<sctbcoeff>[Mm][ \t\n]*[0-9]*  {	
								CopyString( species, ( char * )yytext );
								SkipWhites( species );
								UpperString( species );
								strcpy( thirdBodyName, species);
								if ( ListIterator( gUsedThirdBodyList, FindString, species ) ) {
									if ( !( ListIterator( gThirdBodyList, FindThirdBody, species ) ) ) {
										AddThirdBody( gThirdBodyList, species );
									}
									else {
										ScanError( "the third body efficiency '%s' is already defined", species, 0 );
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "BEGIN(scSkipM)\n" );
											}
										BEGIN(scSkipM);
									}
								}
								else {
									ScanError( "you don't need to define the third body efficiency '%s', because you don't use it", species, 0 );
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(scSkipM)\n" );
										}
									BEGIN(scSkipM);
								}
								
						}*/
<sctbcoeff>"."					{	
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(INITIAL)\n" );
										}
									BEGIN(INITIAL);
								}
<sctbcoeff>{INT}|{FLOAT1}|{FLOAT2}|{FLOAT1}{EXP}|{FLOAT2}{EXP}|{INT}{EXP}		{	
										if (gOptions->debugScanner == TRUE) {
											fprintf(stderr, "i've found the third body coefficient '%s': %g\n", ( char * )yytext, coefficient );
											fprintf( stderr, "BEGIN(scMdefinition)\n" );
										}
									BEGIN(scMdefinition);
									coefficient = atof( ( char * )yytext );
								}
<sctbcoeff>{LETTER}				{	
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(scMdefinition)\n" );
										}
									BEGIN(scMdefinition);
									yyless( 0 );
								}
<sctbcoeff>\[			

%{
/*  scSkipM  */
%}
<scSkipM>[bB][eE][ \t\n]*\.			{
										if ( !errorRecovery ) {
											errorRecovery = TRUE;
											return ERROR; /* error */
										}
									}
<scSkipM>[A-Za-z]+[0-9]*[ \t\n]*\.	{
										errorRecovery = FALSE;
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "BEGIN(INITIAL)\n" );
											}
										BEGIN(INITIAL);
										return '.';
									}
<scSkipM>.|\n

%{
/*  scSkipLet  */
%}
<scSkipLet>\.		{
						errorRecovery = FALSE;
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "BEGIN(INITIAL)\n" );
							}
						BEGIN(INITIAL);
						return '.';
					}
<scSkipLet>.|\n		{		
						if ( !errorRecovery ) {
							errorRecovery = TRUE;
							return ERROR; /* error */
						}
					}


%{
/*  scSkipSpecies  */
%}
<scSkipSpecies>\}					{
										errorRecovery = FALSE;
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "BEGIN(INITIAL)\n" );
											}
										BEGIN(INITIAL);
										return '}';
									}
<scSkipSpecies>.|\n					{
										if ( !errorRecovery ) {
											errorRecovery = TRUE;
											return ERROR; /* error */
										}
									}


%{
/*  allowed atoms  */
%}
<allowatoms>=|[Bb][Ee]|","	/* skip */
<allowatoms>{ID}	{	char	at[2];
						if ( yyleng < 3 ) {
							strcpy( at, ( char * )yytext );
						}
						else {
							ScanError( "name of atom '%s' too long\n# maximum size of atomnames is 2 chars", ( char * )yytext, 0 );
							exit( 2 );
						}
						UpperString( at );
						if ( !( ListIterator( gAllowAtomList, FindString, at ) ) ) {
							AddString( gAllowAtomList, at );
						}

					}
<allowatoms>"."		{	gComposition = NewIntVector( gAllowAtomList->items );
						
							if (gOptions->debugScanner == TRUE) {
								fprintf( stderr, "BEGIN(INITIAL)\n" );
							}
						BEGIN(INITIAL);
					}

%{
/*  Additional Species  */
%}
<scAdditionalSpecies>=|[Bb][Ee]    /* skip */
	/*  isomeres  */
<scAdditionalSpecies>({ISOMERE})\-	{	
										AdjustNumberOfLines();
										yymore();
											if (gOptions->debugScanner == TRUE) {
												fprintf(stderr, "i've found an isomere introduced by '%s'\n", ( char * )yytext );
											}
									}
	/*  species  */
<scAdditionalSpecies>{LETTER}/{LETTER}	{	
									lastAtom = MatchAtom( ( char * )yytext, yyleng, 2 );
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, " lastAtom  = %s\n",  lastAtom  );
											fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
										}
									if ( lastAtom == 0 ) {  /* error handling, a message has already printed out in MatchAtom  */
											if (gOptions->debugScanner == TRUE) {
												fprintf( stderr, "BEGIN(scSkipSpecies)\n" );
											}
										BEGIN(scSkipSpecies);
									}	
									gLineContents[strlen( gLineContents ) - 1] = '\0';
									AdjustNumberOfLines();
									yymore();
								}
<scAdditionalSpecies>{LETTER}	{	
								lastAtom = MatchAtom( ( char * )yytext, yyleng, 1 );
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, " lastAtom  = %s\n",  lastAtom  );
										fprintf( stderr, "yytext  = %s\n",  ( char * )yytext  );
									}
								if ( lastAtom == 0 ) {  /* error handling, a message has already printed out in MatchAtom  */
										if (gOptions->debugScanner == TRUE) {
											fprintf( stderr, "BEGIN(scSkipSpecies)\n" );
										}
									BEGIN(scSkipSpecies);
								}	
								AdjustNumberOfLines();
								yymore();
							}
<scAdditionalSpecies>[0-9]+	{	
								int i, num;
								AtomsPtr atom = NULL;
								LinkPtr link = NULL;
								
								for ( i = yyleng; isdigit( yytext[i-1] ); --i ) { /* null statement */; }
								num = atoi( ( char * )&yytext[i] );
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, " lastAtom  = %s\n",  lastAtom  );
									}
								link = ListIterator( gAtomsList, FindAtom, lastAtom );
								atom = ( AtomsPtr )link->item;
								gComposition->vec[atom->number] += num - 1;
									if (gOptions->debugScanner == TRUE) {
										fprintf( stderr, "i've found the number %d following the atom %s\n", num, lastAtom );
									}
								AdjustNumberOfLines();
								yymore();
							}
<scAdditionalSpecies>[,.]	{	
						LinkPtr 	link = NULL;
						SpeciesPtr	speciesPtr = NULL;	

						CopyString( species, ( char * )yytext );
						species[ strlen(species)-1 ] = '\0'; /*  skip the last char of species, that is ',' or '.'  */
						SkipWhites( species );
						UpperString( species );
							if (gOptions->debugScanner == TRUE) {
								fprintf(stderr, "Species of length %d found: '%s'\n", (int)strlen(species), species );
							}
						if ( !( link = ListIterator( gSpeciesList, gSpeciesCompareFunc, species ) ) ) {
							speciesPtr = AddSpecies( gSpeciesList, species, gComposition );
							gComposition = NewIntVector( gAllowAtomList->items );
						}
						else {
							speciesPtr = ( SpeciesPtr )link->item;
							ClearComposition( gComposition );
						}
						if ( yytext[yyleng-1] == '.' ) {
								if (gOptions->debugScanner == TRUE) {
									fprintf(stderr, "end of declation of additional species\n" );
									fprintf( stderr, "BEGIN(INITIAL)\n" );
								}
							BEGIN(INITIAL);
						}
					}
<scAdditionalSpecies>{WHITESPACE}	{	
								AdjustNumberOfLines();
								yymore();
							}

%{
/*  anything else  */
%}
<sclet>.   		{ 
					ScanError( "i can't match the token '%s'", ( char * )yytext, 0);
						if (gOptions->debugScanner == TRUE) {
							fprintf( stderr, "BEGIN(scSkipLet)\n" );
						}
					BEGIN(scSkipLet);
				}
<scCoeff,scSpecies>.	{ 
							ScanError( "i can't match the token '%s'", ( char * )yytext, 0);
								if (gOptions->debugScanner == TRUE) {
									fprintf( stderr, "BEGIN(scSkipSpecies)\n" );
								}
							BEGIN(scSkipSpecies);
						}
<scMoreInfo>.   { 
					ScanError( "i can't match the token '%s'", ( char * )yytext, 0);
						if (gOptions->debugScanner == TRUE) {
							fprintf( stderr, "BEGIN(scSkipSpecies)\n" );
						}
					BEGIN(scSkipSpecies);
				}
<allowatoms>.   {
					ScanError( "i can't match the token '%s'", ( char * )yytext, 0);
						if (gOptions->debugScanner == TRUE) {
							fprintf( stderr, "BEGIN(scSkipLet)\n" );
						}
					BEGIN(scSkipLet);
				}
<scUnits>.   	{
					if ( !gErrorOut ) {
						ScanError( "i can't match the token '%s', maybe you have forgotten to delimit an exponent", ( char * )yytext, 0);
					}
						if (gOptions->debugScanner == TRUE) {
							fprintf( stderr, "BEGIN(scSkipLet)\n" );
						}
					BEGIN(scSkipLet);
				}
<INITIAL>.   	{
					ScanError( "a '%s' is a complete surprise for me in this context\n# i think you've forgotten to limit the label with a colon", ( char * )yytext, 0);
						if (gOptions->debugScanner == TRUE) {
							fprintf( stderr, "BEGIN(comment)\n" );
						}
					BEGIN(comment);
				}
<scFc,scMdefinition,sctbcoeff,scAdditionalSpecies,scSteadyStates,scSkipM,scSkipSpecies>.   	{
					fprintf( stderr, "a '%s' is a complete surprise for me in this context\n", ( char * )yytext );
				}
<scSkipLet,preUnits,scexponent,scOrderOfReaction,scTemperatureExponent,scCComment>.   	{
					fprintf( stderr, "a '%s' is a complete surprise for me in this context\n", ( char * )yytext );
				}
<scPAHSymbol,scSootSymbol>.   	{
					fprintf( stderr, "a '%s' is a complete surprise for me in this context\n", ( char * )yytext );
				}

%{
/*  end of file  */
%}
<<EOF>>             { return 0; }

%%

int yywrap()
	{
	return 1;
	}
