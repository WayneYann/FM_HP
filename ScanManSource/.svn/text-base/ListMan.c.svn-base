/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 */


#include "ListMan.h"
#undef qDebug


void AddItem( void *item, ListPtr list )
{
	LinkPtr	link = NULL;
	
	link = NewLink( item );
	
	link->last = list->lastItem;
	if ( list->items > 0 ){
		list->lastItem->next = link;
	}
	else {
		list->firstItem = link;
	}
	list->lastItem = link;
	++list->items;
	
}

void AppendLink( LinkPtr link, ListPtr list )
{
	if ( link == NULL ) {
		fprintf( stderr, "# error in AppendLink: address of specified Item is NULL\n" );
		exit( 2 );
	}
	if ( list->items == 0 ) {
		list->firstItem = link;
	}
	else {
		list->lastItem->next = link;
	}
	link->last = list->lastItem;
	list->lastItem = link;
	link->next = NULL;
	++list->items;
}

LinkPtr NewLink( void *item )
{
	LinkPtr	link;
	
	if ( (link = (LinkPtr) malloc( sizeof(Link))) == NULL){
		fprintf( stderr, "# memory allocation of Link_struct failed\n");
		exit(2);
	}

	link->next = NULL;
	link->last = NULL;
	link->item = item;

	return link;
}

ListPtr NewList( char* name )
{
	ListPtr	list;
	char	*listName = NULL;
	
	if ( (list = (ListPtr) malloc( sizeof(List))) == NULL){
		fprintf( stderr, "# memory allocation of List_struct failed\n");
		exit(2);
	}

	listName = (char *) malloc( strlen( name ) + 1 );
	strcpy(listName, name);
	list->listName = listName;
	list->firstItem = NULL;
	list->lastItem = NULL;
	list->items = 0;
	
	return list;
}

LinkPtr ListIterator( ListPtr list, int (*func)( LinkPtr link, void *var ), void *var )
{	
	LinkPtr	link = NULL;
	int 	err = 0;
	
	for ( link = list->firstItem; link != NULL && err == 0; ) {
		err = func( link, var);
		if ( !err && link ) {
			link = link->next;
		}
	}

	return link;
}


void ListInfo( ListPtr list, FILE *fp )
{
	fprintf( fp, "\nThe List '%s' has %d Items.\n", list->listName, list->items );
}

void FreeList( ListPtr list, int (*func)( LinkPtr link, void *var ) )
{
	if ( func ) {
		ListIterator( list, func, list );
	}
	free( list->listName );
	free( list );
}

LinkPtr Find_n_th_Item( int n, ListPtr list )
{
	int i;
	LinkPtr link = list->firstItem;
	
	if ( n == 0 ) {
		return 0;
	}
	
	for ( i = 0; i < n-1; ++i, link = link->next ) { /* null statement */; }

	return link;
}

int PrintLinks( LinkPtr link, void *var )
{
	ListPtr 		list = var;
	if ( link == list->firstItem ) {
		fprintf(stdout, "%s\n", list->listName );
		fprintf(stdout, "firstItem = %d lastItem = %d\n", (int)list->firstItem, (int)list->lastItem );
	}
	fprintf( stdout, "link->last = %d link = %d link->next = %d\n", (int)link->last, (int)link, (int)link->next );
	
	return 0;
}

int FreeString( LinkPtr link, void *var )
{
	ListPtr list = var;
	int		err = 0;
	
	if ( !(list->firstItem == link ) ) {
		link = link->last;
		FreeOneString( link );
		free( link );
	}
	if (link == list->lastItem ) {
		FreeOneString( link );
		free( link );
		err = 1;
	}
	return err;
}

int FreeOneString( LinkPtr link )
{
	char	*string = NULL;
	
	string = ( char * )link->item;
	free( string );

	return 0;
}

int FindString( LinkPtr link, void *var )
{
	char		*name = var;
	char		*string = link->item;
	
	if ( strcmp( name, string ) == 0 ) {
		return 1;
	} 
	else {
		return 0;
	}
}

int PrintString( LinkPtr link, void *var )
{
	ListPtr 		list = var;	

	if ( link == list->firstItem ) {
		ListInfo( list, stdout );	
		fprintf( stdout, "%s are ", list->listName );
	}

	if ( link->next ) {
		fprintf( stdout, "%s, ", link->item );
	}
	else {
		fprintf( stdout, "%s.\n\n", link->item );
	}
	return 0;
}

char *AddString( ListPtr list, char *name )
{
	char	*string = NULL;
	
	string = ( char* ) malloc( strlen( name ) + 1 );
	strcpy( string, name );
	AddItem( string, list );
	
	return string;
}

int FindInt( LinkPtr link, void *var )
{
	int	*compNumber = var;
	int	*inumber = link->item;
	
	if ( *compNumber == *inumber ) {
		return 1;
	} 
	else {
		return 0;
	}
}

int PrintInt( LinkPtr link, void *var )
{
	ListPtr 		list = var;	

	if ( link == list->firstItem ) {
		ListInfo( list, stdout );	
		fprintf( stdout, "%s are ", list->listName );
	}

	if ( link->next ) {
		fprintf( stdout, "%d, ", (int)link->item );
	}
	else {
		fprintf( stdout, "%d.\n\n", (int)link->item );
	}
	return 0;
}

int *AddInt( ListPtr list, int ivalue )
{
	int	*inumber = NULL;
	
	inumber = ( int * ) malloc( sizeof( int ) );
	*inumber = ivalue;
	AddItem( inumber, list );
	
	return inumber;
}

int FreeInt( LinkPtr link, void *var )
{
	ListPtr list = var;
	int		err = 0;
	
	if ( !(list->firstItem == link ) ) {
		link = link->last;
		FreeOneString( link );
		free( link );
	}
	if (link == list->lastItem ) {
		FreeOneString( link );
		free( link );
		err = 1;
	}
	return err;
}

int FreeOneInt( LinkPtr link )
{
	int		*inumber = NULL;
	
	inumber = ( int * )link->item;
	free( inumber );

	return 0;
}

int NumberOfListItem( ListPtr list, int (*func)( LinkPtr link, void *var ), char *name )
{
	int	number = 0;
	LinkPtr link = NULL;
	
	for ( link = list->firstItem; link != NULL; link = link->next ) {
		++number;
		if ( func( link, name ) ) {
			return number;
		}
	}
	
	return 0;
}

int NumberOfItems( ListPtr list )
{
	LinkPtr link = NULL;
	int 	number = 0;
	
	for ( link = list->firstItem; link != NULL; link = link->next ) {
		++number;
	}

	if ( number != list->items ) {
		fprintf( stderr, "# i count %d items for the list '%s', but the listinternal\n# value is %d\n", number, list->listName, list->items );
	}

	return number;
}

/*  AllOnesInTwo returns TRUE if all items of list one appears in list two, 
 *  otherwise FALSE
 *  list one has to be a list of strings
 *  criterium is the equality of the char *name
 */
int AllOnesInTwo( ListPtr one, ListPtr two, int (*func)( LinkPtr link, void *var ), char **error )
{
	LinkPtr link1 = NULL;
	LinkPtr link2 = NULL;
	char	*stringOne = NULL;
	int 	i;
	
	for ( i = 0; i < NumberOfItems( one ); ++i ) {
		link1 = Find_n_th_Item( i+1, one );
		stringOne = ( char * )link1->item;
		if ( !ListIterator( two, func, stringOne ) ) {
			*error = stringOne;
			return FALSE;
		}
	}
	
	return TRUE;
}

void KillItem( ListPtr list, int (*func)( LinkPtr link ), LinkPtr link ) 
{
	RemoveLinkFromList( list, link );

/*  free link and item  */
	func( link );
	free( link );
}

void RemoveLinkFromList( ListPtr list, LinkPtr link )
{
	LinkPtr workLink = NULL;
	
/*  test specified link  */
	if ( link == NULL ) {
		fprintf( stderr, "# error in RemoveLinkFromList: address of specified Item is NULL\n" );
		exit( 2 );
	}
	
	for ( workLink = list->firstItem; workLink != link && workLink != NULL; workLink = workLink->next ) {
		if ( workLink == NULL ) {
			fprintf( stderr, "# error in RemoveLinkFromList: can't find specified item\n" );
			exit( 2 );
		}
	}
	
/*  update list  */
	if ( link == list->firstItem ) {
		if ( list->items == 1 ) {
			list->firstItem = NULL;
			list->lastItem = NULL;
		}
		else {
			list->firstItem = link->next;
			link->next->last = NULL;
		}
	}
	else if ( link == list->lastItem ) {
		list->lastItem = link->last;
		link->last->next = NULL;
	}
	else {
		link->last->next = link->next;
		link->next->last = link->last;
	}
	--list->items;
}

void RemoveLink( ListPtr list, LinkPtr link )
{
	RemoveLinkFromList( list, link );
	free( link );
}

/* the function InsertItem has not the best performance, but is easy to
 * implement, because it needs only the addfunction as an itemspecific 
 * function; it first adds the item, then moves its link to the specified
 * location
 */
LinkPtr InsertItem( ListPtr list, void *(*func)( ListPtr list, void *var ), void *var, int appendTo )
{
	LinkPtr link = NULL;
	LinkPtr newLink = NULL;

/* check */
	if ( appendTo + 1 > list->items ) {
		fprintf( stderr, "#error: something wrong in function InsertItem" );
		exit(2);
	}

/* addItem */
	func( list, var );
	newLink = list->lastItem;

/* find link */
	link = Find_n_th_Item( appendTo+1, list );

/*  update list  */
	/* first remove link from end of list */
	newLink->last->next = NULL;
	list->lastItem = newLink->last;
	
	/* then place it to the specified location */
	if ( link == list->lastItem ) {
		link->next = newLink;
		newLink->last = link;
		list->lastItem = newLink;
		newLink->next = NULL;
	}
	else {
		link->next->last = newLink;
		newLink->next = link->next;
		
		link->next = newLink;
		newLink->last = link;
	}
	
	return newLink;
}


void SortList( ListPtr list, Flag (*CompFunc)( LinkPtr link1, LinkPtr link2 ) )
{
	int		i;
	int 	items;
	Flag	sorted = TRUE;
	LinkPtr link = NULL;

/* check */
	if ( list->items == 0 ) {
		fprintf( stderr, "#error: something wrong in function SortList\n" );
		exit(2);
	}

	items = NumberOfItems( list );

	do {
		sorted = TRUE;
		link = list->firstItem;
		for ( i = 0; i < items-1; ++i ) {
			if ( CompFunc( link, link->next ) ) {
				SwapItems( list, link );
				sorted = FALSE;
			}
			else {
				link = link->next;
			}
		}
	} while ( !sorted );
}

void SwapItems( ListPtr list, LinkPtr link )
{
	LinkPtr link2 = link->next;
	
	if ( link->next == NULL ) {
		return;
	}

	if ( link == list->firstItem && link2 == list->lastItem ) {
		list->firstItem = link2;
		list->lastItem = link;
		link->next = NULL;
		link2->last = NULL;
	}
	else if ( link == list->firstItem ) {
		list->firstItem = link2;
		link2->next->last = link;
		link->next = link2->next;
		link2->last = NULL;
	}
	else if ( link2 == list->lastItem ) {
		list->lastItem = link;
		link->last->next = link2;
		link2->last = link->last;
		link->next = NULL;
	}
	else {
		link2->next->last = link;
		link->last->next = link2;
		link2->last = link->last;
		link->next = link2->next;
	}
	link2->next = link;
	link->last = link2;
}

int CountItems( const ListPtr list )
{
	extern void Warning( const char * );
	LinkPtr	link = list->firstItem;
	int		items = 0;

	while ( link ) {
		++items;
		link = link->next;
	}
	
	if ( list->items != items )
		Warning( "Ooops, number of list items is wrong" );
	
	return items;
}


void Iterate( ListPtr list, actionFunc func, void *aux )
{
	LinkPtr	link = list->firstItem;
	
	while ( link ) {
		(*func)( link->item, aux );
		link = link->next;
	}
}


Pointer FindListItem( ListPtr list, Pointer ptr, int (*compare)(Pointer, Pointer) )
{
	LinkPtr	link = list->firstItem;

	while ( link ) {
		if ( (*compare)( link->item, ptr ) ) return link->item; 
		link = link->next;
	}
	return NULL;
}


Pointer ListItem( ListPtr list, int item )
{
	int		i;
	LinkPtr	link = list->firstItem;
	
	for ( i = 0; i < item; ++i ) {
		link = link->next;
		if ( link == NULL ) return NULL;
	}
	
	return link->item;
}
