/*
	Stack.c: package for a stack of strings
	
	To do:
		-	check for correct access
		-	visit function (starting at the head)
		-	print function
		-	search function
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alligator.h"
#include "Stack.h"


typedef struct Node {
	char		*key;
	struct Node	*next;
} Node, *NodePtr;


static NodePtr head;				/*	first element of stack	*/
static NodePtr z;					/*	sentinel element		*/
static NodePtr t;					/*	temporary node			*/

static int	   NumElems;			/* current number of elements in stack */

void FLRStackInit(void) 
{
	head = NEW( Node );
	z = NEW( Node );
	head->next = z; 
	head->key = NULL;
	z->next = z;
	NumElems = 0;
}

void FLRStackDelete(void)
{
	DELETE(head);
	DELETE(z);
}

void FLRStackPush(char *v)
{
	char *v1;
	
	t = NEW( Node );
	v1 = NEW2(strlen(v)+1,char);
	strcpy(v1,v);
	t->key = v1;
	t->next = head->next;
	head->next = t;
	++NumElems;
}

char *FLRStackPop(void)
{
	char *x;
	
	t = head->next;
	head->next = t->next;
	x = t->key;
	DELETE(t);
	--NumElems;
	return x;
	
}

int FLRStackEmpty(void)
{
	return (head->next == z);
}

char **FLRStackSearch(char *v)
{
	int		i = 0;
	char	*x;
	
	t = head->next;
	while( i++ < NumElems ) {
		x = t->key;
		if( strcmp(x,v) == 0 ) 	return &x;
		else 					t = t->next;
	}
	return NULL;
}

void FLRStackPrint(void)
{
	int	i;
	
	t = head->next;
	for(i=0;i< NumElems; ++i) {
		fprintf(stderr,"%s\n",t->key);
		t = t->next;
	}
	
}

void FLRStackReverse(void)
{
	NodePtr temp1,temp2,temp3;
	
	temp1 = head->next;temp2 = head->next;temp3 = head->next;

	if (NumElems >= 3) {
		temp1 = temp3->next;
		temp2 = temp1->next;
		while((temp2 != z->next)) {
			temp1->next = temp3;
			temp3 = temp1;
			temp1 = temp2;
			temp2 = temp2->next;
		}
		head->next->next = temp2;
		head->next = temp1;
		temp1->next = temp3;
		
	}
	else if(NumElems == 2) {
		temp2 = temp1->next;
		temp2->next = temp1;
		temp1->next = z;
		head->next = temp2;
	}
}


