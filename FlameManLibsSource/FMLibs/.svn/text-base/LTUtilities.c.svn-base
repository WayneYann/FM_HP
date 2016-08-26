#define FOOL_SOFTBENCH( x ) 

/*
	LTUtilities.c: utilities for the ListTool program
*/

#include "ListTool.h"


char *lowercase(register char *str)
{
	register char *s = str;
	
	do {
		*s = tolower(*s);
	} while ( *s++ );
	return ( str );
}


void FatalError(const char *str)
{
	fprintf(stderr, "### Error: %s.\n", str);
	exit(2);
}


void Warning(char *str)
{
	fprintf(stderr, "### Warning: %s.\n", str);
}


void Message(char *str)
{
	fprintf(stderr, "# %s", str);
}


void printHelp( char *tname )
{
	fprintf(stdout, "Usage:\t%s -H, or\n", tname);
	fprintf(stdout, "\t\t%s -i <input file> [options]\n", tname);
	fprintf(stdout, "\t\t\tOptions:\n");
	fprintf(stdout, "\t\t\t-i <input file>: the name of the mandatory input file\n");
	fprintf(stdout, "\t\t\t-o <output file>: save tabulator separated data in <output file>\n");
	fprintf(stdout, "\t\t\t-s <symbol file>: compare symbols with the set given in <symbol file>\n");
	fprintf(stdout, "\t\t\t-h: print the header in a tabular form to stdout\n");
	fprintf(stdout, "\t\t\t-t: print the content of the trailer to stdout\n");
	fprintf(stdout, "\t\t\t-c: check the format of the input file (no output is produced)\n");
	fprintf(stdout, "\t\t\t-p: write some progress information to stderr\n");
	fprintf(stdout, "\t\t\t-T: write timing information to stderr\n");
	fprintf(stdout, "\t\t\t-H: write this help info to stdout\n\n");
}


int moveto(char *str, FILE *fp)
{
	char	buffer[50];
	int		error;
	
	lowercase(str);
	do {
		error = fscanf(fp, "%s", buffer);
		if (error == 0)
			FatalError("Unable to read string");
		else if (error == EOF) {
			Warning("Couldn't find search string");
			return ( 1 );
		}
	} while ( strcmp(str, lowercase(buffer)) );
	return ( noErr );
}


int isReserved(char *str)
{
	if ( streq(str, "header") ||  streq(str, "body") || streq(str, "trailer")
		|| streq(str, "begin") || streq(str, "end") )
		return ( TRUE );
	else
		return ( FALSE );
}


void strip(register char *str)	/* removes all space characters from a string */
{
	register char *ptr = str;
	
	do {
		if ( !isspace(*str) ) *ptr++ = *str;
	} while ( *str++ );
	*ptr = '\0';
}





