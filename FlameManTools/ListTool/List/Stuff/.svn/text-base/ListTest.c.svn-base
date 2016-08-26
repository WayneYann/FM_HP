/*
	ListTest.c: program to test the list package.

	HP: cc -Aa -DHP -z -I/jgLib ListTest.c list.a /jgLib/alligator.a -o test
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef applec
#include <CursorCtl.h>
#endif


#include "alligator.h"
#include "List.h"


void print_int( void *item, void *aux )
{
	fprintf( (FILE *)aux, "   %d", *(int *)item );
}


void print_string( void *item, void *aux )
{
	fprintf( (FILE *)aux, "   %s", (char *)item );
}


int Test1( void )
{
	int	a[] = { 3, 65, 3, 38, 1, 15, 3, 95, 8, 42, 9, 4, -1 };
	int	*ptr = NULL;
	ListPtr	intList = NewList( "integer list" );
	
	ptr = a;
	while ( *ptr > 0 ) AddItem( intList, ptr++, kListHead );
	fprintf( stdout, "List of %d integers:\n", ItemsInList(intList) );
	ForAll( intList, print_int, stdout, kListHead );
	fprintf( stdout, "\n\n" );
	fflush( stdout );
	fprintf( stdout, "List of %d integers (reversed):\n", ItemsInList(intList) );
	ForAll( intList, print_int, stdout, kListTail );
	fprintf( stdout, "\n\n" );
	fflush( stdout );
	RemoveList( intList );
	
	fprintf( stdout, "Allocate in reverse order:\n" );
	intList = NewList( "integer list #2" );
	ptr = a;
	while ( *ptr > 0 ) AddItem( intList, ptr++, kListTail );
	fprintf( stdout, "Reversed list of %d integers:\n", ItemsInList(intList) );
	ForAll( intList, print_int, stdout, kListHead );
	fprintf( stdout, "\n\n" );
	fflush( stdout );
	RemoveList( intList );
	
	return 0;
}


int Test2( void )
{
	char	*s[] = { "eins", "zwei", "drei", "vier", "fünf", "sechs", 
					 "sieben", "acht", "neun", "zehn", NULL };
	char	*ptr = NULL;
	int		i;
	ListPtr	list = NewList( "string list" );
	ListPtr clone = NULL;

	for ( i = 0; s[i] != NULL; ++i ) {
		ptr = NEW2( strlen(s[i])+1, char );
		strcpy( ptr, s[i] );
		AddItem( list, ptr, kListHead );
	}
	fprintf( stdout, "List of strings:\n" );
	ForAll( list, print_string, stdout, kListTail );
	fprintf( stdout, "\n\n" );
	fflush( stdout );
	
	clone = CloneList( list, "cloned list" );
	fprintf( stdout, "Cloned list of strings:\n" );
	ForAll( list, print_string, stdout, kListTail );
	fprintf( stdout, "\n\n" );
	fflush( stdout );
	
	DeleteList( list );
	RemoveList( clone );
	
	return 0;
}


int main( void )
{
#ifdef applec
	InitCursorCtl( NULL );
#endif

	fprintf( stderr, "\nTest #1: Simple allocation and deallocation of int list.\n" );
	if ( Test1() != 0 )	FATAL( "Test 1 failed" );
	else				fprintf( stderr, "Test 1 successful.\n\n" );
	
	fprintf( stderr, "\nTest #2: CloneList and DeleteList.\n" );
	if ( Test2() != 0 )	FATAL( "Test 2 failed" );
	else				fprintf( stderr, "Test 2 successful.\n\n" );
	
	return 0;
}
