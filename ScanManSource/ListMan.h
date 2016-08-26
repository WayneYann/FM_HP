/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 */

#ifndef __ListMan__
#define __ListMan__

#ifndef stdin
#include<stdio.h>
#endif
#include<stdlib.h>
#include<string.h>

/*  debugging tools  */
/*#define prints( a ) fprintf( stderr, #a " = '%s'\n", ( char * )a )*/
/*#define print( a ) fprintf( stderr, #a "\n" )*/
/*#define printg( a ) fprintf( stderr, #a " = %g\n", a )*/
/*#define printd( a ) fprintf( stderr, #a " = %d\n", a )*/

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

typedef char Flag;

typedef struct Link {
	void		*item;
	struct Link	*next;
	struct Link	*last;
} Link, *LinkPtr;

typedef struct List {
	LinkPtr	firstItem;
	LinkPtr	lastItem;
	int		items;
	char	*listName;
} List, *ListPtr;


typedef int (*ListFuncPtr)( LinkPtr link, void *var );


void RemoveLinkFromList( ListPtr list, LinkPtr link );

void AddItem	( void *item, ListPtr list );

void AppendLink( LinkPtr Link, ListPtr list );

LinkPtr NewLink	( void *item );

ListPtr NewList	( char* name );

LinkPtr ListIterator( ListPtr list, int (*func)( LinkPtr link, void *var ), void *var );

void ListInfo	( ListPtr list, FILE *fp );

void FreeList( ListPtr list, int (*func)( LinkPtr link, void *var ) );

LinkPtr Find_n_th_Item( int n, ListPtr list );

int PrintLinks( LinkPtr link, void *var );

int FreeString( LinkPtr link, void *var );

int FreeOneString( LinkPtr link );

int FindString	( LinkPtr link, void *var );

int PrintString	( LinkPtr link, void *var );

char *AddString( ListPtr list, char *name );

int FreeInt( LinkPtr link, void *var );

int FindInt( LinkPtr link, void *var );

int PrintInt( LinkPtr link, void *var );

int *AddInt( ListPtr list, int ivalue );

int NumberOfListItem( ListPtr list, int (*func)( LinkPtr link, void *var ), char *name );

int NumberOfItems( ListPtr list );

int AllOnesInTwo( ListPtr one, ListPtr two, int (*func)( LinkPtr link, void *var ), char **error );

void KillItem( ListPtr list, int (*func)( LinkPtr link ), LinkPtr link ) ;

void RemoveLink( ListPtr list, LinkPtr link );

LinkPtr InsertItem( ListPtr list, void *(*func)( ListPtr list, void *var ), void *var, int appendTo );

void SortList( ListPtr list, Flag (*CompFunc)( LinkPtr link1, LinkPtr link2 ) );

void SwapItems( ListPtr list, LinkPtr link );

/*	jg stuff	*/
typedef void* Pointer;
typedef void (*actionFunc)(Pointer, Pointer);

void Iterate( ListPtr list, actionFunc, Pointer );
int CountItems( const ListPtr list );
Pointer FindListItem( ListPtr list, Pointer, int (*)(Pointer, Pointer) );
Pointer ListItem( ListPtr list, int item );

#endif /* __ListMan__ */
