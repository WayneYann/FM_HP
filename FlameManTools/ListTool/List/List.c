/*
	List.[ch]: singly linked list package.

	Requires the "alligator" package.
	
	Library:
	
		Begin
			C -b3 -r -mc68020 -mc68881 -elems881 List.c
			Lib List.c.o -o "{myLibs}"List.lib
			C -b3 -r -mc68020 -mc68881 -d qMemDebug -elems881 List.c
			Lib List.c.o -o "{myLibs}"List.dbg.lib
			Duplicate -y List.h "{myIncludes}"
		End …… "{Worksheet}"

*/

#include <string.h>
#include <stdlib.h>
#include <string.h>

#include "List.h"

#ifdef applec
#pragma segment List
#endif


ListPtr NewList( const char *label )
{
	ListPtr	list = NEW( List );
	if (label != NULL)
		strncpy( list->fLabel, label, 39 );
	return list;
}


LinkPtr NewLink( LinkPtr prev, LinkPtr next, void *item )
{
	LinkPtr	link = NEW( Link );
	link->fNext = next;
	link->fPrev = prev;
	link->fItem = item;
	return link;
}


void AddItem( ListPtr list, void *item, int where )
{
	LinkPtr	temp = NULL;
	if ( list->fHead == NULL ) {		/*	list is empty	*/
		temp = NewLink( NULL, NULL, item );
		list->fTail = list->fHead = temp;
	} else if ( where == kListTail ) {	/*	insert at tail	*/
		temp = NewLink( list->fTail, NULL, item );
		list->fTail->fNext = temp;
		list->fTail = temp;
	} else {							/*	insert at head	*/
		temp = NewLink( NULL, list->fHead, item );
		list->fHead->fPrev = temp;
		list->fHead = temp;
	}
	++list->fNumItems;	
}


void RemoveItem( ListPtr list, void *item )
{
	/*	Remove the link, but not the object itself.						*/
	LinkPtr	temp = list->fHead;
	while ( temp != NULL ) {
		if ( temp->fItem == item ) {
			if ( temp == list->fHead )	{			/*	remove head		*/
				list->fHead = temp->fNext;
				temp->fNext->fPrev = temp->fPrev;	/*	i.e. NULL		*/
			} else if ( temp == list->fTail ) {	/*	remove tail		*/
				list->fTail = temp->fPrev;
				temp->fPrev->fNext = temp->fNext;	/*	i.e. NULL		*/
			} else {
				temp->fPrev->fNext = temp->fNext;
				temp->fNext->fPrev = temp->fPrev;
			}
			DELETE( temp );
			--list->fNumItems;
			return;
		}
		temp = temp->fNext;
	}
}


void ForAll( ListPtr list, void (*DoIt)(void *, void *), void *aux, int start )
{
	LinkPtr temp = NULL;
	if ( start == kListTail ) {
		temp = list->fTail;
		while ( temp ) {
			(*DoIt)( temp->fItem, aux );
			temp = temp->fPrev;
		}
	} else {
		temp = list->fHead;
		while ( temp ) {
			(*DoIt)( temp->fItem, aux );
			temp = temp->fNext;
		}
	}
}


void RemoveList( ListPtr list )
{
	/*	Remove all links, but not the objects.							*/
	LinkPtr temp = list->fHead, next = NULL;
	while ( temp ) {
		next = temp->fNext;
		DELETE( temp );
		temp = next;
	}
	DELETE( list );
}



void DeleteList( ListPtr list )
{
	/*	Remove all links and all objects.  This assumes that all
		objects have been allocated with a single NEW or NEW2 call.		*/
	LinkPtr temp = list->fHead, next = NULL;
	while ( temp ) {
		next = temp->fNext;
		DELETE( temp->fItem );
		DELETE( temp );
		temp = next;
	}
	DELETE( list );
}



void *FindItem( ListPtr list, void *target, int (*cmp )(void *, void *) )
{
	LinkPtr temp = list->fHead;
	while ( temp ) {
		if ( (*cmp)( temp->fItem, target ) ) return temp->fItem;
		temp = temp->fNext;
	}
	return NULL;
}


int ItemsInList(ListPtr list)
{
	return list->fNumItems;
}


ListPtr CloneList(ListPtr list, const char *new_label)
{
	LinkPtr	temp = list->fHead;
	ListPtr	clone = NewList( new_label );
	while ( temp ) {
		AddItem( clone, temp->fItem, kListHead );
		temp = temp->fNext;
	}
	return clone;
}


void DisposeList(ListPtr list)
{
	RemoveList( list );
}

