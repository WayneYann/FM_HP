/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 */


#ifdef qUseDump
#pragma load "ScanMan.dump"
#else
#include "ScanMan.h"
#endif


AtomsPtr NewAtoms( char *name )
{
	AtomsPtr			atoms;
	
	if ( (atoms = (AtomsPtr) malloc( sizeof(Atoms))) == NULL){
		fprintf( stderr, "# memory allocation of Atoms_struct failed\n");
		exit(2);
	}

	atoms->name = ( char * ) malloc( strlen( name ) + 1 );
	strcpy( atoms->name, name );

	return atoms;
}

SpeciesPtr NewSpecies( char *name )
{
	SpeciesPtr			species;
	
	if ( (species = (SpeciesPtr) malloc( sizeof(Species))) == NULL){
		fprintf( stderr, "# memory allocation of Species_struct failed\n");
		exit(2);
	}

	species->name = ( char * ) malloc( strlen( name ) + 1 );
	strcpy( species->name, name );
	species->isSteadyState = FALSE;

	return species;
}

ReactionPtr NewReaction( char *label )
{
	ReactionPtr			reaction;
	int					i;
	
	if ( (reaction = (ReactionPtr) malloc( sizeof(Reaction))) == NULL){
		fprintf( stderr, "# memory allocation of Reaction_struct failed\n");
		exit(2);
	}
	memset( reaction, 0, sizeof(*reaction) );

	reaction->label = ( char * ) malloc( strlen( label ) + 1 );
	strcpy( reaction->label, label );

	reaction->numberOfSpecies = 0;

	for ( i = 0; i < MaxSpeciesOfReaction; ++i) {
		reaction->speciesNumber[i] = 0;
		reaction->speciesCoeff[i] = 0.0;
	}

	reaction->a = 0.0;
	reaction->n = 0.0;
	reaction->e = 0.0;

	reaction->withLindemann = FALSE;

	reaction->aInf = 0.0;
	reaction->nInf = 0.0;
	reaction->eInf = 0.0;

	reaction->lindemannNumber = -1;

	reaction->fca = 0.0;
	reaction->fcb = 0.0;
	reaction->fcc = 0.0;
	reaction->fcTa = 0.0;
	reaction->fcTb = 0.0;
	reaction->fcTc = 0.0;
	
	reaction->withThirdBody = FALSE;
	reaction->thirdBodyNumber = -1;

	reaction->partialEquilibrium = FALSE;

	reaction->orderOfReaction = 0;
	
	reaction->forward = reaction->backward = NULL;
	reaction->isDouble = FALSE;
	reaction->rateFromEquilibrium = FALSE;


	return reaction;
}

ThirdBodyPtr NewThirdBody( char *name )
{
	ThirdBodyPtr			thirdBody;
	int						i;
	
	if ( (thirdBody = (ThirdBodyPtr) malloc( sizeof(ThirdBody))) == NULL){
		fprintf( stderr, "# memory allocation of ThirdBody_struct failed\n");
		exit(2);
	}

	thirdBody->name = ( char * ) malloc( strlen( name ) + 1 );
	strcpy( thirdBody->name, name );
	thirdBody->speciesNumber = NewIntVector( gSpeciesList->items );
	thirdBody->speciesCoeff = NewVector( gSpeciesList->items );
	thirdBody->set = NewIntVector( gSpeciesList->items );
	for ( i = 0; i < thirdBody->speciesNumber->len; ++i ) {
		thirdBody->speciesNumber->vec[i] = i;
	}

	return thirdBody;
}

DimensionPtr NewDimension( char *name, Double value, int kg, int m, int s, int K, int mole )
{
	DimensionPtr			dimension;
	
	if ( (dimension = (DimensionPtr) malloc( sizeof( Dimension ) ) ) == NULL ){
		fprintf( stderr, "# memory allocation of Dimension_struct failed\n");
		exit(2);
	}

	dimension->name = ( char * ) malloc( strlen( name ) + 1 );
	strcpy( dimension->name, name );
	dimension->value = value;
	dimension->kg = kg;
	dimension->m = m;
	dimension->s = s;
	dimension->K = K;
	dimension->mole = mole;
	
	dimension->kgParaCoeff = 0;
	dimension->mParaCoeff = 0;
	dimension->sParaCoeff = 0;
	dimension->KParaCoeff = 0;
	dimension->moleParaCoeff = 0;

	dimension->orderOfReaction = FALSE;
	dimension->tempExponent = FALSE;

	return dimension;
}

ParameterPtr NewParameter( char *name )
{
	ParameterPtr			parameter;
	
	if ( (parameter = (ParameterPtr) malloc( sizeof(Parameter))) == NULL){
		fprintf( stderr, "# memory allocation of Parameter_struct failed\n");
		exit(2);
	}

	parameter->name = ( char * ) malloc( strlen( name ) + 1 );
	strcpy( parameter->name, name );
	
	parameter->value = 1.0;

	return parameter;
}

int FreeAtoms( LinkPtr link, void *var )
{
	ListPtr list = var;
	int		err = 0;
	
	if ( !(list->firstItem == link ) ) {
		link = link->last;
		FreeOneAtoms( link );
		free( link );
	}
	if ( list->lastItem == link ) {
		FreeOneAtoms( link );
		free( link );
		err = 1;
	}
	
	return err;
}

int FreeSpecies( LinkPtr link, void *var )
{
	ListPtr list = var;
	int		err = 0;
	
	if ( !(list->firstItem == link ) ) {
		link = link->last;
		FreeOneSpecies( link );
		free( link );
	}
	if ( list->lastItem == link ) {
		FreeOneSpecies( link );
		free( link );
		err = 1;
	}
	
	return err;
}

int FreeReaction( LinkPtr link, void *var )
{
	ListPtr list = var;
	int		err = 0;
	
	if ( !(list->firstItem == link ) ) {
		link = link->last;
		FreeOneReaction( link );
		free( link );
	}
	if (list->lastItem == link ) {
		FreeOneReaction( link );
		free( link);
		err = 1;
	}
	
	return err;
}

int FreeThirdBody( LinkPtr link, void *var )
{
	ListPtr list = var;
	int		err = 0;
	
	if ( link != list->firstItem ) {
		link = link->last;
		FreeOneThirdBody( link );
		free( link );
	}
	if ( link == list->lastItem ) {
		FreeOneThirdBody( link );
		free( link );
		err = 1;
	}
	
	return err;
}

int FreeDimension( LinkPtr link, void *var )
{
	ListPtr list = var;
	int		err = 0;
	
	if ( !(list->firstItem == link ) ) {
		link = link->last;
		FreeOneDimension( link );
		free( link );
	}
	if ( list->lastItem == link ) {
		FreeOneDimension( link );
		free( link );
		err = 1;
	}
	
	return err;
}

int FreeParameter( LinkPtr link, void *var )
{
	ListPtr list = var;
	int		err = 0;
	
	if ( !(list->firstItem == link ) ) {
		link = link->last;
		FreeOneParameter( link );
		free( link );
	}
	if ( list->lastItem == link ) {
		FreeOneParameter( link );
		free( link );
		err = 1;
	}
	
	return err;
}

int FreeOneAtoms( LinkPtr link )
{
	AtomsPtr atom;
	
	atom = ( AtomsPtr )link->item;
	free( atom->name);
	free( atom );

	return 0;
}

int FreeOneSpecies( LinkPtr link )
{
	SpeciesPtr species;
	
	species = ( SpeciesPtr )link->item;
	free( species->name );
	if ( species->omegaCoeff ) {
		DisposeVector( species->omegaCoeff );
	}
	if ( species->dCoeff ) {
		DisposeVector( species->dCoeff );
	}
	DisposeIntVector( species->composition );
	free( species );

	return 0;
}

int FreeOneReaction( LinkPtr link )
{
	ReactionPtr reaction = NULL;
	
	reaction = ( ReactionPtr )link->item;
	free( reaction->label );
	free( reaction );

	return 0;
}

int FreeOneThirdBody( LinkPtr link )
{
	ThirdBodyPtr thirdBody = NULL;
	
	thirdBody = ( ThirdBodyPtr )link->item;
	free( thirdBody->name );
	DisposeVector( thirdBody->speciesCoeff );
	DisposeIntVector( thirdBody->speciesNumber );
	DisposeIntVector( thirdBody->set );
	free( thirdBody );

	return 0;
}

int FreeOneDimension( LinkPtr link )
{
	DimensionPtr dimension;
	
	dimension = ( DimensionPtr )link->item;
	free( dimension->name);
	free( dimension );

	return 0;
}

int FreeOneParameter( LinkPtr link )
{
	ParameterPtr parameter;
	
	parameter = ( ParameterPtr )link->item;
	free( parameter->name);
	free( parameter );

	return 0;
}

AtomsPtr AddAtom( ListPtr list, char *name )
{
	AtomsPtr atom = NULL;
	
	atom = NewAtoms( name );
	atom->number = list->items;
	AddItem( atom, list );

	return atom;
}

SpeciesPtr AddSpecies( ListPtr list, char *name, IntVectorPtr composition )
{
	SpeciesPtr species = NULL;
	
	species = NewSpecies( name );
	species->composition = composition;
	species->number = list->items;
	species->molarMass = 0.0;
	species->dCoeff = NULL;
	species->omegaCoeff = NULL;
	AddItem( species, list );

	return species;
}

ReactionPtr AddReaction( ListPtr list, char *name )
{
	static int	id = 0;					/* (jg) this is new!		*/
	ReactionPtr	reaction = NULL;
	int			len;
	
	reaction = NewReaction( name );
	reaction->id = id++;
	AddItem( reaction, list );

	reaction->hasForward = FALSE;
	reaction->hasBackward = FALSE;
	len = strlen( reaction->label );
	if ( gOptions->useForwardBackward ) {
		if ( reaction->label[len-1] == 'f' ) {
			reaction->hasBackward = TRUE;
		}
		else if ( reaction->label[len-1] == 'b' ) {
			reaction->hasForward = TRUE;
		}
	}

	return reaction;
}

void *AddReactionVoid( ListPtr list, void *var )
{
	char *name = ( char * ) var;
	
	return AddReaction( list, name );
}

ReactionPtr InsertReaction( ListPtr list, char *name, int appendTo )
{
	LinkPtr		link = NULL;
	ReactionPtr	reaction = NULL;
	
	link = InsertItem( list, AddReactionVoid, name, appendTo );
	reaction = ( ReactionPtr )link->item;

	return reaction;
}

ThirdBodyPtr AddThirdBody( ListPtr list, char *name )
{
	static int id = 0;
	ThirdBodyPtr thirdBody = NULL;
	
	thirdBody = NewThirdBody( name );
	/*thirdBody->id = list->items; (jg) This could mean trouble! */
	/*thirdBody->id = id++;*/
	thirdBody->id = NumberOfListItem( gUsedThirdBodyList, FindString, name ) - 1;
	AddItem( thirdBody, list );

	return thirdBody;
}

DimensionPtr AddDimension( ListPtr list, char *name, Double value, int kg, int m, int s, int K, int mole )
{
	DimensionPtr dimension = NULL;
	
	dimension = NewDimension( name, value, kg, m, s, K, mole );
	AddItem( dimension, list );
	return dimension;
}

ParameterPtr AddParameter( ListPtr list, char *name )
{
	ParameterPtr parameter = NULL;
	
	parameter = NewParameter( name );
	AddItem( parameter, list );

	return parameter;
}

