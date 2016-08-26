#define FOOL_SOFTBENCH( x ) 

/*
	LTBody.c: body functions
*/

#include "ListTool.h"

#undef DEBUG

/* Globals */
extern char		*gLabels;
extern double	*gData;
extern int		gRows;
extern int		gVectors;


static void addLabel(char *str)
{
	static char	*ptr;
	static int	len, init = FALSE;
	
	if (!init) {
		ptr = gLabels;
		init = TRUE;
	}
	len += strlen(str) + 1;
	if (len > kMaxLabelChars) FatalError("Too many label characters");
	strcpy(ptr, str);
#	ifdef DEBUG
	fprintf(stderr, "Label: \"%s\", current length: %d\n", ptr, len);
#	endif
	ptr += strlen(str) + 1;
}


static int expandArray(void)
{
	static size_t size = 0;
	
	size += kPool * sizeof(double);
	if ( !(gData = (doublePtr)realloc(gData, size)) )
		FatalError("Not enough room in heap zone");
#	ifdef DEBUG
	fprintf(stderr, "¥ New size of gData: %u\n", size);
#	endif
	return ( kPool );
}

/*static int expandArray(void)
{
	static size_t size = 0;
	
	size += kPool * sizeof(double);
	if ( !gData ) {
		if ( !( gData = (doublePtr)malloc( gData, size ) ) )
			FatalError("malloc of gData failed");
	}
	else {
		if ( !(gData = (doublePtr)realloc(gData, size)) )
			FatalError("Not enough room in heap zone");
#ifdef DEBUG
		fprintf(stderr, "¥ New size of gData: %u\n", size);
#endif
	}
	return ( kPool );
}
*/

static void addNumber(char *str)
{
	static int	left = 0,		/* number of empty slots */
				counter = 0;	/* total number of fp values stored so far */
	
	if (!left) left = expandArray();
	sscanf(str, "%lf", &gData[counter++]);
	--left;
#	if defined (applec) || defined (powerc)
	RotateCursor(counter);
#	endif
}


static int readArray(FILE *fp, int *error)
{
	static int	first = TRUE;
	char		buffer[128], unit[128];
	char		*bptr, *uptr;
	int 		done = FALSE;
	int			len = 0;

	*error = noErr;
	do {
		*error = fscanf(fp, "%s", buffer);
		if ( *error == EOF ) return ( len );
		if ( streq(lowercase(buffer), "trailer") ) {
			*error = kAtTrailer;
			return ( len );
		}
#		ifdef DEBUG
		/*fprintf(stderr, "> buffer = %s, *error = %d\n", buffer, *error); */
#		endif
		switch ( buffer[0] ) {
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
				addNumber(buffer);
				++len;
				break;
			case '=':
			/* ignore equal signs */
				break;
			case '"':
				FatalError("strings are not allowed within body");
				break;
			default:
				if ( isReserved(buffer) ) FatalError("Syntax error");
				if (!first) done = TRUE;
				else		first = FALSE;
				/*if  ( fscanf(fp, "%[^+-.0-9]", unit) == 1 ) {*/		/* handle units */
				if  ( fscanf(fp, " [%s]", unit) == 1 || fscanf(fp, "%[^+-.0-9]", unit) == 1 ) {		/* handle units */
					/*fprintf(stderr, "matching '%s' in '%s'\n", unit, buffer);*/
					strip(unit);
					if ( unit[0] ) {
						bptr = buffer + strlen(buffer);
						uptr = unit;
						*bptr++ = ' ';
						while ( *uptr ) {
						  if ( !isspace(*uptr) && *uptr != '=' ) *bptr++ = tolower(*uptr);
						  ++uptr;
						}
						*bptr = '\0';
					}
				}
				addLabel(buffer);
				break;
		}
	} while (!done);
	return ( len );
}


void parseBody(FILE *fp)
{
	int		vLength,		/* length of a vector */
			error;			/* error flag */
	
	/* read first vector */
	gRows = readArray(fp, &error);
	if ( (error == EOF && gRows == 0) || (error == kAtTrailer && gRows == 0) ) {
		Warning("Nothing in body");
		return;
	}
	gVectors = 1;
#	ifdef DEBUG
	fprintf(stderr, "¥ gRows = %d\n", gRows);
#	endif
	
	/* read the rest and check vector lengths	 */
	while ( vLength = readArray(fp, &error) ) {
		if ( vLength != gRows ) {
			fprintf(stderr, "# last correct vector was no. %d\n", gVectors);
			FatalError("input vectors have different lengths!");
		}
		++gVectors;
		if (error == EOF || error == kAtTrailer) break;
	}
#	ifdef DEBUG
	fprintf(stderr, "¥ gRows = %d, gVectors = %d, error = %d\n", gRows, gVectors, error);
#	endif
}


void saveData(FILE *fp)
{
	char	*ptr = gLabels;
	int		i, j, increment = 0;
	
	/* save column names first, ... */
	fprintf(fp, "*\n");
	for ( j = 0; j < gVectors - 1; ++j) {
		fprintf(fp, "%s\t", ptr);
		ptr += strlen(ptr) + 1;
	}
	fprintf(fp, "%s\n", ptr);			/* last string terminates with a newline char */
	/* ... and then the data itself */
	for ( i = 0; i < gRows; ++i ) {
		for ( j = 0; j < gVectors - 1; ++j ) {
			fprintf(fp, "%g\t", gData[i + j * gRows]);
#			if defined (applec) || defined (powerc)
			RotateCursor(--increment);
#			endif
		}
		fprintf(fp, "%g\n", gData[i + j * gRows]);
	}
}
