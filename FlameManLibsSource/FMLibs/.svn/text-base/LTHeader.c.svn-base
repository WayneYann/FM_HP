#define FOOL_SOFTBENCH( x ) 

/*
	LTHeader.c: header functions for ListTool
*/

#include "ListTool.h"

#undef DEBUG

/* Globals */

extern int			gNumOfParas;
extern parameter	*gParas;

char	gTokenList[kMaxTokens];
char	*gDefaultContext = "*";
char	*gHeaderTooLarge = "Header too large, \"Body\" forgotten?";
int		gNumOfTokens = 0;
int		gInGroup = FALSE;		/* is TRUE if inside of a group */


#ifdef DEBUG
void printTableEntry(parameter *p, int i)
{
	fprintf(stderr, "(%d) identifier=%s, context=%s, tag=%d, ",
		i, p->identifier, p->context, p->tag);
	if ( p->tag == kTextType )
		fprintf(stderr, "string=%s\n", p->what.string);
	else
		fprintf(stderr, "value=%g, unit=%s\n",
			p->what.quantity.value, p->what.quantity.unit);
}
#endif


static int isp1(int i)			/* e.g. pressure (atm) = 5.0 */
{
	if ( i + 3 > gNumOfTokens ) return ( FALSE );
	if (	gTokenList[i] == kSymbol && gTokenList[i+1] == kUnit
		 && gTokenList[i+2] == kAssignment && gTokenList[i+3] == kNumber )
		 return ( TRUE );
	else
		return ( FALSE );
}

static int isp2(int i)			/* e.g. pressure = 5.0 (atm) */
{
	if ( i + 3 > gNumOfTokens ) return ( FALSE );
	if (	gTokenList[i] == kSymbol && gTokenList[i+1] == kAssignment
		 && gTokenList[i+2] == kNumber && gTokenList[i+3] == kUnit )
		 return ( TRUE );
	else
		return ( FALSE );
}

int isp3(int i)			/* e.g. pressure = 5.0 */
{
	if ( i + 2 > gNumOfTokens ) return ( FALSE );
	if (	gTokenList[i] == kSymbol && gTokenList[i+1] == kAssignment
		 && gTokenList[i+2] == kNumber )
		 return ( TRUE );
	else
		return ( FALSE );
}

int isp4(int i)			/* e.g. author = "claire" */
{
	if ( i + 2 > gNumOfTokens ) return ( FALSE );
	if (	gTokenList[i] == kSymbol && gTokenList[i+1] == kAssignment
		 && gTokenList[i+2] == kString )
		 return ( TRUE );
	else
		return ( FALSE );
}

int isp5(int i)			/* begin of group, e.g. unburnt begin */
{
	if ( i + 2 > gNumOfTokens ) return ( FALSE );
	if ( gTokenList[i] == kSymbol && gTokenList[i+1] == kBegin )
		return ( TRUE );
	else
		return ( FALSE );
}

int isp6(int i)			/* begin of group, e.g. unburnt state begin */
{
	if ( i + 3 > gNumOfTokens ) return ( FALSE );
	if ( gTokenList[i] == kSymbol && gTokenList[i+1] == kSymbol
		&& gTokenList[i+2] == kBegin )
		return ( TRUE );
	else
		return ( FALSE );
}


int isp7(int i)			/* gInGroup must be true to signal the end of a group */
{
	if ( gTokenList[i] == kEnd )
		return ( TRUE );
	else
		return ( FALSE );
}


void makeSymbolTable(char *stream)
{
	static	char currentContext[kMaxLenOfContext];
	static	char noUnit[] = "n/a";
	char	*token;
	int		i = 0;			/* index to current token */
	int		items = 0;		/* items read by sscanf */
	double	x;
	
	token = strtok(stream, kTokenSepStr);
#	ifdef DEBUG
	fprintf(stderr, "%d>  token -> %s\n", i, token);
#	endif
	strcpy(currentContext, gDefaultContext);
	for (i = 1; i < gNumOfTokens; ) {
		if ( isp1(i) ) {				/* e.g. pressure (atm) = 5 */
			token = strtok(NULL, kTokenSepStr);			/* read identifier */
			items = sscanf(token, "%s", gParas[gNumOfParas].identifier);
			token = strtok(NULL, kTokenSepStr);			/* read unit */
			items = sscanf(token, "%s", gParas[gNumOfParas].what.quantity.unit);
			/*strcpy(gParas[gNumOfParas].what.quantity.unit, token);*/
			token = strtok(NULL, kTokenSepStr);			/* read value */
			items = sscanf(token, "%lf", &x);
			gParas[gNumOfParas].what.quantity.value = x;
			strcpy(gParas[gNumOfParas].context, currentContext);
			gParas[gNumOfParas].tag = kPhysicalQuantity;
			++gNumOfParas;
			i += 4;
#			ifdef DEBUG
			printTableEntry(&gParas[gNumOfParas-1], gNumOfParas);
#			endif
		} else if ( isp2(i) ) {				/* e.g. pressure = 5 (atm) */
			token = strtok(NULL, kTokenSepStr);			/* read identifier */
			items = sscanf(token, "%s", gParas[gNumOfParas].identifier);
			/*strcpy(gParas[gNumOfParas].what.quantity.unit, token);*/
			token = strtok(NULL, kTokenSepStr);			/* read value */
			items = sscanf(token, "%lf", &x);
			token = strtok(NULL, kTokenSepStr);			/* read unit */
			items = sscanf(token, "%s", gParas[gNumOfParas].what.quantity.unit);
			gParas[gNumOfParas].what.quantity.value = x;
			strcpy(gParas[gNumOfParas].context, currentContext);
			gParas[gNumOfParas].tag = kPhysicalQuantity;
			++gNumOfParas;
			i += 4;
#			ifdef DEBUG
			printTableEntry(&gParas[gNumOfParas-1], gNumOfParas);
#			endif
		} else if ( isp3(i) ) {			/* e.g. pressure = 5 */
			token = strtok(NULL, kTokenSepStr);
			items = sscanf(token, "%s", gParas[gNumOfParas].identifier);
			token = strtok(NULL, kTokenSepStr);
			items = sscanf(token, "%lf", &x);
			gParas[gNumOfParas].what.quantity.value = x;
			strcpy(gParas[gNumOfParas].context, currentContext);
			gParas[gNumOfParas].tag = kPhysicalQuantity;
			strcpy(gParas[gNumOfParas].what.quantity.unit, noUnit);
			++gNumOfParas;
			i += 3;
#			ifdef DEBUG
			printTableEntry(&gParas[gNumOfParas-1], gNumOfParas);
#			endif
		} else if ( isp4(i) ) {			/* e.g. title = "xyz" */
			token = strtok(NULL, kTokenSepStr);
			items = sscanf(token, "%s", gParas[gNumOfParas].identifier);
			token = strtok(NULL, kTokenSepStr);
			strcpy(gParas[gNumOfParas].what.string, token);
			strcpy(gParas[gNumOfParas].context, currentContext);
			gParas[gNumOfParas].tag = kTextType;
			++gNumOfParas;
			i += 3;
#			ifdef DEBUG
			printTableEntry(&gParas[gNumOfParas-1], gNumOfParas);
#			endif
		} else if ( isp5(i) ) {			/* e.g. unburnt begin */
			token = strtok(NULL, kTokenSepStr);
			items = sscanf(token, "%s", currentContext);
			gInGroup = TRUE;
			i += 2;
#			ifdef DEBUG
			fprintf(stderr, "new context: %s\n", currentContext);
#			endif
		} else if ( isp6(i) ) {	/* e.g. unburnt state begin */
			token = strtok(NULL, kTokenSepStr);
			items = sscanf(token, "%s", currentContext);
			token = strtok(NULL, kTokenSepStr);
			gInGroup = TRUE;
			i += 3;
#			ifdef DEBUG
			fprintf(stderr, "new context: %s\n", currentContext);
#			endif
		} else if ( isp7(i) && gInGroup) {	/* end of group
			strcpy(currentContext, gDefaultContext); */
			gInGroup = FALSE;
			++i;
#			ifdef DEBUG
			fprintf(stderr, "¥ back to global context\n");
#			endif
		} else FatalError("Syntax or semantic error in header");
	}
#	ifdef DEBUG
	fprintf(stderr, "¥ %d assignment statements parsed.\n", gNumOfParas);
#	endif
}


char *trim(char *str)	/* remove leading and trailing spaces from strings */
{
	register char *ptr, *s = str;
	
	ptr = s + strlen(s) - 2;			/* skip to 2nd last character */
	while ( isspace(*ptr) ) --ptr;		/* find last nonspace character */
	*(ptr+1) = '"';						/* and terminate string properly */
	*(ptr+2) = '\0';
	ptr = ++s;							/* skip to 2nd character */
	if ( !isspace(*ptr) ) return (str);	/* we don't need to move characters */
	while ( isspace(*ptr) ) ++ptr;		/* find 1st nonspace character */
	do {								/* move characters, trailing '\0' included */
		*s++ = *ptr;
	} while ( *ptr++ );
	return (str);
}


int getTokens(char *stream, FILE *fp)
{
	register char *ptr = stream;
			 char *maxPtr = stream + kMaxHeader;
			 char str[256];
	register int  j, k, c;
		   	 int  inHeader = TRUE;
	long	 int  tokCount = 0;				/* tokCount counts the number of tokens */
	
	do {
		do {
			c = getc(fp);
		} while ( isspace(c) || iscntrl(c) );	/* skip over white space characters */
#		ifdef DEBUG2
		fprintf(stderr, "> <%c> (%d)\n", c, c);
#		endif
		if ( c == EOF ) return ( c );
		switch ( c ) {
			case '=':
				gTokenList[tokCount++] = kAssignment;
				break;
			case '(':
			case '[':
			case '{':
				while ( (c = getc(fp)) != EOF ) {
					if ( iscntrl(c) ) continue;		/* always skip control chars */
					if (c == ')' || c == '}' || c == ']') break;
					/* check if room for 2 chars (c and token separator) */
					if ( !isspace(c) && ptr + 2 <= maxPtr ) *ptr++ = c;
					else	FatalError(gHeaderTooLarge);
				}
				if ( c == EOF ) return ( c );
				*ptr++ = kTokenSeparator;
				gTokenList[tokCount++] = kUnit;
				break;
			case '"':
				j = 0;
				str[j++] = c;
				while ( (c = getc(fp)) != EOF ) {
					if ( iscntrl(c) ) continue;
					if ( j < 255 )	str[j++] = c;
					else			break;
					if (c == '"' )	break;
				}
				if ( c == EOF ) return ( c );
				str[j] = '\0';		/* make it a regular C string */
				trim(str);			/* remove leading and trailing spaces */
				k = strlen(str);
				str[k] = kTokenSeparator;
				str[++k] = '\0';
				if ( ptr + k <= maxPtr )	strcpy(ptr, str);
				else						FatalError(gHeaderTooLarge);
				while ( *ptr ) ++ptr;
				gTokenList[tokCount++] = kString;
				break;
			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':
			case '+':
			case '-':
			case '.':
				if ( ungetc(c, fp) != c ) FatalError("ungetc() didn't work");
				if ( ptr + kSafety <= maxPtr ) {
					if ( fscanf(fp, "%s", ptr) != 1 ) return ( EOF );
				} else
					FatalError(gHeaderTooLarge);
				while ( *ptr ) ++ptr;
				*ptr++ = kTokenSeparator;
				gTokenList[tokCount++] = kNumber;
				break;
			default:	/* token may be either a reserved word or an identifier */

				(void)ungetc(c, fp);
				if ( fscanf(fp, "%s", str) != 1 ) return ( EOF );
				lowercase(str);

				if ( streq(str, "begin") ) {
					gTokenList[tokCount++] = kBegin;
				} else if ( streq(str, "end") ) {
					gTokenList[tokCount++] = kEnd;
				} else if ( streq(str, "header") ) {
					gTokenList[tokCount++] = kHeader;
					strcpy(ptr, str);
					while ( *ptr ) ++ptr;
					*ptr++ = kTokenSeparator;
				} else if ( streq(str, "body") ) {
				/*gTokenList[tokCount++] = kBody; */
					inHeader = FALSE;
				} else {
					gTokenList[tokCount++] = kSymbol;
					strcpy(ptr, str);
					while ( *ptr ) ++ptr;
					*ptr++ = kTokenSeparator;
				}
		}
	} while ( inHeader );
	*ptr = '\0';
	gTokenList[tokCount] = kEndOfList;
	gNumOfTokens = tokCount;
#	ifdef DEBUG
	fprintf(stderr, "Content of stream:\n%s\n", stream);
	fprintf(stderr, "gNumOfTokens = %d\ngTokenList:\n", gNumOfTokens);
	for ( tokCount = 0; tokCount < gNumOfTokens; ++tokCount )
		fprintf(stderr, "%d> %d\n", tokCount, gTokenList[tokCount]);
#	endif
	return ( noErr );
}

void parseHeader(char *stream, FILE *fp)
{
	if ( getTokens(stream, fp) == noErr && gTokenList[0] == kHeader)
		makeSymbolTable(stream);
	else
		FatalError("Couldn't assemble symbol table");
}
