/*
	MapMan.cp:	package for manipulating sets of profiles.
	
	History:
		09/06/93		first version (jg)
*/


#ifndef __dump__
#include <iostream> // GB
#ifdef __ZTC__
#include <strstrea.h>
#else
#include <sstream>  // GB
#endif
#include <fstream>  // GB
#include <cstdio>   // GB
#include <cstdlib>  // GB
#include <cstring>  // GB
#include <cmath>    // GB
#include <cassert>  // GB
#if defined (applec) || defined (powerc)
# include <CursorCtl.h>
# include "getopt.h"
#elif defined(__ZTC__)
# include "getopt.h"
#else
# include <unistd.h>
#endif
#include "alligator.h"
#include "ArrayManager.h"
#endif
#include "MapMan.h"

#define SKIPOMPARSE

#if defined (applec) || defined (powerc)
# pragma segment MapMan
#endif

//		Class MMDataSet
//--------------------------------------------------------------------------------


MMDataSet::MMDataSet()
{
	fVec = NULL;
	fLen = 0;
	fOffset = 0;
	//fName[0] = '\0';
	//ostrstream(fName, sizeof(fName)) << "<anonymous>" << ends; C++ Standard : Deprecated
	strcpy( fName, "<anonymous>" );
	fExtrapolate = lastValue;
	fBag = NULL;

#ifdef qDebug
	cerr << "-> MMDataSet: anonymous allocator called.\n";
#endif
}


MMDataSet::MMDataSet( Double* v, int n, int off, const char* name,
	MMDataBag* bag )
{
	assert( (v != NULL) && (n > 0) && (off > 0) );
	fVec = v;
	fLen = n;
	fOffset = off;
	if ( name != NULL )
	  //ostrstream(fName, sizeof(fName)) << name << ends;  C++ Standard : Deprecated
	        strcpy( fName, name );
	else
	        fName[0] = '\0';
	fExtrapolate = lastValue;
	fBag = bag;

#ifdef qDebug
	cerr << "-> MMDataSet: normal allocator called.\n";
#endif
}


MMDataSet::MMDataSet( VectorPtr v, const char* name, MMDataBag* bag )
{
	assert( (v != NULL) && (v->len > 0) );
	fVec = v->vec;
	fLen = v->len;
	fOffset = 1;
	if ( name != NULL )
	        //ostrstream(fName, sizeof(fName)) << name << ends;  C++ Standard : Deprecated
	        strcpy( fName, name );
	else
	        fName[0] = '\0';
	fExtrapolate = lastValue;
	fBag = bag;

#ifdef qDebug
	cerr << "-> MMDataSet: VectorPtr allocator called.\n";
#endif
}


void
MMDataSet::FreeMemory( void )
{
	if ( IsValid() ) {
		DELETE( fVec );
		fVec = NULL;
		fLen = 0;
	}
}


int
MMDataSet::IsSameDataSet( const Double* /*v*/, int len, int off, const char* name ) const
{
	if ( (len == fLen) && (off == fOffset) && (strcmp(name,fName) == 0) )
		return TRUE;
	else
		return FALSE;
}


void
MMDataSet::SetXType( MMDataSet::XType type )
{
	fExtrapolate = type;
}

#ifdef SKIPOMPARSE
MMDataSet::XType
#endif
MMDataSet::GetXType( void )
{
	return fExtrapolate;
}



#if 0
MMDataSet::~MMDataSet()
{
#ifdef qDebug
	cerr << "About to destruct " << fName << endl;
#endif
}
#endif


ostream& operator<<( ostream& os, MMDataSet& ds )
{
	if ( ds.IsValid() ) {
		os << ds.fName << " (" << ds.fLen << ", ";
		char buf[40];
		switch ( ds.fExtrapolate ) {
			case MMDataSet::zeroValue:
			        //ostrstream( buf, sizeof(buf) ) << "zero" << ends; C++ Standard : Deprecated
			        strcpy( buf, "zero" );
				break;
			case MMDataSet::lastValue:
			        //ostrstream( buf, sizeof(buf) ) << "last value" << ends; C++ Standard : Deprecated
			        strcpy( buf, "last value" );
				break;
			case MMDataSet::zeroGradient:
			        //ostrstream( buf, sizeof(buf) ) << "zero gradient" << ends; C++ Standard : Deprecated
			        strcpy( buf, "zero gradient" );
				break;
			case MMDataSet::lastGradient:
			        //ostrstream( buf, sizeof(buf) ) << "last gradient" << ends; C++ Standard : Deprecated
			        strcpy( buf, "last gradient" );
				break;
			default:
			        //ostrstream( buf, sizeof(buf) ) << "<invalid type>" << ends; C++ Standard : Deprecated
			        strcpy( buf, "<invalid type>" );
		}
		os << buf << "):";
		Double* ptr = ds.fVec;
		for ( int i = 0; i < ds.fLen; ++i, ptr += ds.fOffset ) {
			if ( (i % 8) == 0 ) os << '\n';
			os << "  " << *ptr;
		}
	}
	os << endl;

	return os;
}


const Double*
MMDataSet::GetData( void ) const
{
	return (const Double*)fVec;
}


void
MMDataSet::Map( Double* y_new, int len, int offset )
{
	//	Do some error checking.
	assert( fBag->GetNewInpedVar()->Length() == len );

	const Double*	x_new = fBag->GetNewInpedVar()->GetData();
	Double*		y_old = fVec;		// offset is fOffset
	const Double*	x_old = fBag->GetOldInpedVar()->GetData();
	int				i;					// simple loop counter
	
	//	Find the number of left extrapolation, interpolation,
	//	and right extrapolation points.
	int i1 = 0;
	int i2 = len;
	{	const Double* ptr = x_new;
		while ( *ptr++ < x_old[0] ) ++i1;
		ptr = x_new + len - 1;
		while ( *ptr-- > x_old[fLen-1] ) --i2;
	}
#ifdef qDebug
	cerr << "i1 = " << i1 << ", i2 = " << i2 << ", len = " << len << endl;
#endif

	
	//	Left side extrapolation.
	for ( i = 0; i < i1; ++i, ++x_new, y_new += offset ) {
		switch ( fExtrapolate ) {
			case zeroValue:
				*y_new = 0.0;
				break;
			case lastValue:
			case zeroGradient:
				*y_new = y_old[0];
				break;
			case lastGradient: {
				const FPUType m = (y_old[2*fOffset] - y_old[0]) / (x_old[2] - x_old[0]);
					// we are using 2 nodes to approximate the derivative
				const FPUType a = y_old[0] - m * x_old[0];
				*y_new = a + m * *x_new;
				} break;
			default:
				FATAL( "unknown extrapolation type" );				
		}
	}
	
	//	Linear interpolation.
	const FPUType kTiny = 1.0e-16;
	for ( i = i1; i < i2; ++i, ++x_new, y_new += offset ) {
		while ( *x_old < *x_new ) { ++x_old; y_old += fOffset; }
		if ( fabs(*x_old - *x_new) < kTiny ) {
			*y_new = *y_old;
		} else {
			register FPUType factor = (*x_new - x_old[-1]) / (x_old[0] - x_old[-1]);
			assert( factor >= 0.0 && factor <= 1.0 );
			*y_new = y_old[-fOffset] + factor * (y_old[0] - y_old[-fOffset]);
		}
#ifdef qDebug
		cerr << y_old[-fOffset] << "  " << *y_new << "  " << y_old[0] 
			 << "  (" << factor << ")\n";
#endif
	}

	//	Right side extrapolation

	//	y_old and x_old must be reset, because the linear interpolation
	//	algorithm modifies these pointers.
	y_old = fVec;
	x_old = fBag->GetOldInpedVar()->GetData();

	for ( i = i2; i < len; ++i, ++x_new, y_new += offset ) {
		switch ( fExtrapolate ) {
			case zeroValue:
				*y_new = 0.0;
				break;
			case lastValue:
			case zeroGradient:
				*y_new = y_old[(fLen-1)*fOffset];
				break;
			case lastGradient: {
				const FPUType m = (y_old[(fLen-1)*fOffset] - y_old[(fLen-3)*fOffset]) / (x_old[fLen-1] - x_old[fLen-3]);
					// we are using 2 nodes to approximate the derivative
				const FPUType a = y_old[(fLen-3)*fOffset] - m * x_old[fLen-3];
				*y_new = a + m * *x_new;
				} break;
			default:
				FATAL( "unknown extrapolation type" );				
		}
	}

}


VectorPtr
MMDataSet::Map( void )
{
	const int len = fBag->GetNewInpedVar()->Length();
	VectorPtr v = NewVector( len );
	Map( v->vec, v->len, 1 );
	return v;
}


void
MMDataSet::Map( VectorPtr v )
{
	Map( v->vec, v->len, 1 );
}



//		Class MMDataBag
//--------------------------------------------------------------------------------



MMDataBag::MMDataBag( int maxItems )
{
	fElems = maxItems;
	fCurrent = 0;
	fPool = NULL;
	fArrayLength = 0;
#ifdef qDebug
	cerr << "About to allocate an MMDataBag.\n";
#endif
}


MMDataBag::~MMDataBag()
{
#ifdef qDebug
	cerr << "About to deallocate an MMDataBag.\n";
#endif
	delete [] fPool;
}


void
MMDataBag::Initialize( void )
{
	fPool = new MMDataSet[fElems];	// calls MMDataSet() constructor
	assert( fPool != NULL );
}


ostream& operator<<( ostream& os, MMDataBag& bag )
{
	os << "Bag with " << bag.fCurrent << " items of length " 
	   << bag.fArrayLength << ":\n";
	
	if ( bag.fOldIndepVar.IsValid() ) {
		os << "Old independent variable " << bag.fOldIndepVar;
	}
	if ( bag.fNewIndepVar.IsValid() ) {
		os << "New independent variable " << bag.fNewIndepVar;
	}
	
	for ( int i = 0; i < bag.fCurrent; ++i )
		os << " [" << i << "]  " << bag.fPool[i];
	
	return os;
}


MMDataSet&
MMDataBag::operator[]( int index )
{
	assert( (index >= 0) && (index < fElems) );
	return fPool[index];
}


int
MMDataBag::Insert(  Double* v, int len, int offset, const char *name )
{
	assert( fCurrent < fElems );
	
	//	The first array determines the length of all arrays.
	if ( fArrayLength == 0 ) {
		fArrayLength = len;
	} else {
		assert( fArrayLength == len );
	}
	
	if ( fOldIndepVar.IsSameDataSet( v, len, offset, name ) ) {
		cerr << "# Warning: independent variable \"" << name
			 << "\" will not be added to list of dependent variables.\n";
		return -1;
	} else {
		fPool[fCurrent++] = MMDataSet( v, len, offset, name, this );
		return 0;
	}
}


#define qNewVersion

int
MMDataBag::Insert( VectorPtr v, const char *name )
{
#ifdef qNewVersion

	return Insert( v->vec, v->len, 1, name );

#else
	assert( fCurrent < fElems );
	
	//	The first array determines the length of all arrays.
	if ( fArrayLength == 0 ) {
		fArrayLength = v->len;
	} else {
		assert( fArrayLength == v->len );
	}
	
	if ( fOldIndepVar.IsSameDataSet( v->vec, v->len, 1, name ) ) {
		cerr << "# Warning: independent variable \"" << name
			 << "\" will not be added to list of dependent variables.\n";
		return -1;
	} else {
		fPool[fCurrent++] = MMDataSet( v, name, this );
		return 0;
	}
#endif
}



void
MMDataBag::SetNewInpedVar( Double *v, int len, int off, const char *name )
{
	fNewIndepVar = MMDataSet( v, len, off, name, this );
}

void
MMDataBag::SetNewInpedVar( VectorPtr v, const char *name )
{
	fNewIndepVar = MMDataSet( v, name, this );
}

void
MMDataBag::SetNewInpedVar( MMDataSet& ds )
{
	fNewIndepVar = ds;
}

void
MMDataBag::SetOldInpedVar( Double *v, int len, int off, const char *name )
{
	fOldIndepVar = MMDataSet( v, len, off, name, this );
}

void
MMDataBag::SetOldInpedVar( VectorPtr v, const char *name )
{
	fOldIndepVar = MMDataSet( v, name, this );
}

void
MMDataBag::SetOldInpedVar( MMDataSet& ds )
{
	fOldIndepVar = ds;
}


void
MMDataBag::FreeMemory( void )
{
	fOldIndepVar.FreeMemory();
	fNewIndepVar.FreeMemory();
	for ( int i = 0; i < fElems; ++i ) fPool[i].FreeMemory();
}





//		Class NodeMover
//--------------------------------------------------------------------------------

NodeMover::NodeMover()
{
	fLen = 0;
	fNumToMove = 0;
	fFrom = fromLeftSide;
	fOld = NULL;
	fNew = NULL;
}


NodeMover::NodeMover( Double* oldVar, int len )
{
	SetOldVar( oldVar, len );
}


NodeMover::NodeMover( VectorPtr oldVar )
{
	SetOldVar( oldVar->vec, oldVar->len );
}


NodeMover::NodeMover( Double* oldVar, Double* newVar, int len, 
	int n, fromType from )
{
	SetOldVar( oldVar, len );
	SetNewVar( newVar, len );
	assert( n < len );
	fNumToMove = n;
	fFrom = from;
}


NodeMover::NodeMover( VectorPtr oldVar, VectorPtr newVar, 
	int n, fromType from )
{
	SetOldVar( oldVar->vec, oldVar->len );
	SetNewVar( newVar->vec, newVar->len );
	fNumToMove = n;
	fFrom = from;

}


void
NodeMover::SetOldVar( Double* oldVar, int len )
{
	assert( (oldVar != NULL) && (len > 0) );
	fLen = len;
	fOld = oldVar;	
}


void
NodeMover::SetNewVar( Double* newVar, int len )
{
	assert( (newVar != NULL) && (len == fLen) );
	fNew = newVar;	
}


int
NodeMover::BadMembers( void )
{
	if ( (fLen < 0) || (fLen <= fNumToMove) || (fOld == NULL)
			|| (fNew == NULL) )
		return -1;
	else
		return 0;
}


int
NodeMover::MoveIt( void )
{
	if ( BadMembers() ) FATAL( "NodeMover with invalid data" );
	
	if ( fFrom == fromLeftSide )
		return MoveNodesFromLeft();
	else if ( fFrom == fromRightSide )
		return MoveNodesFromRight();
	else {
		FATAL( "Unknown fromType value" );
		return -1;
	}
}


int
NodeMover::MoveIt( Double* newVar, int len, int n, fromType from )
{
	assert( (newVar != NULL) && (len == fLen) && (n < fLen) );
	fNew = newVar;
	fNumToMove = n;
	fFrom = from;
	return MoveIt();
}


int
NodeMover::MoveIt( VectorPtr newVar, int n, fromType from )
{
	return MoveIt( newVar->vec, newVar->len, n, from );
}


int
NodeMover::MoveNodesFromLeft( void )
{
	int	i;
	const FPUType	dx = fOld[fLen-1] - fOld[fLen-2];
		// use last cell size for the extension
	const int		n = fLen -  fNumToMove;	// number of values to shift
	
	for ( i = 0; i < n; ++i ) fNew[i] = fOld[i+fNumToMove];
	for ( i = n; i < fLen; ++i ) fNew[i] = fNew[i-1] + dx;
	
	return 0;
}


int
NodeMover::MoveNodesFromRight( void )
{
	int	i;
	const FPUType	dx = fOld[1] - fOld[0];
		// use last cell size for the extension
	const int		n = fLen -  fNumToMove;	// number of values to shift
	
	for ( i = 0; i < n; ++i ) fNew[i+fNumToMove] = fOld[i];
	for ( i = fNumToMove-1; i > -1; --i ) fNew[i] = fNew[i+1] - dx;
	
	return 0;
}



ostream& operator<<( ostream& os, NodeMover& nm )
{
	int		i;
	
	os << "NodeMover at $" << &nm << ":\n";
	os << "  fLen       = " << nm.fLen << '\n'
	   << "  fNumToMove = " << nm.fNumToMove << '\n'
	   << "  move nodes from "
	   << ((nm.fFrom == NodeMover::fromLeftSide) ? "left" : "right") << " side\n";

	if ( nm.fOld != NULL ) {
		os << "Old variable:";
		for ( i = 0; i < nm.fLen; ++i ) {
			if ( (i % 8) == 0 ) os << '\n';
			os << "  " << nm.fOld[i];
		}
		os << endl;
	}

	if ( nm.fNew != NULL ) {
		os << "New variable:";
		for ( i = 0; i < nm.fLen; ++i ) {
			if ( (i % 8) == 0 ) os << '\n';
			os << "  " << nm.fNew[i];
		}
		os << endl;
	}

	return os;
}
