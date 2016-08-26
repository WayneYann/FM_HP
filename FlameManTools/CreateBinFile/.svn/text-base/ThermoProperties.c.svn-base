/*
	ThermoProperties.c: functions that compute thermodynamic and transport
						properties.
	
	Revision history:	7/19/91	first version by pt
						7/30/91 cleaned up (jg)
						10/14/91 final library version together with 
								Enthalpy.c MassFractions.c

   changed by go 12/20/94: BR and F added
*/

#undef DEBUGCOMPOSITION
#undef DEBUG

#ifdef DEBUG
#include <stdio.h>
#endif
#ifndef __STDLIB__
#include <stdlib.h>
#endif
#ifndef __CTYPE__
#include <ctype.h>
#endif
#ifndef __MATH__
#include <math.h>
#endif
#ifndef __STRING__
#include <string.h>
#endif
#ifndef __ArrayManager__
#include "ArrayManager.h"
#endif

#include "ThermoProperties.h"


/*
#ifdef applec
#pragma segment ThTrProps
#endif
*/


/*
	NewSpeciesRecordArray allocates a one dimensional Array of SpeciesRecords.
*/

SpeciesRecordPtr NewSpeciesRecordArray( int n )
{
	SpeciesRecordPtr	p = NULL;
	
	if ( !(p = (SpeciesRecordPtr)calloc( n, sizeof(SpeciesRecord) )) )
		FatalError( "Couldn't allocate memory for SpeciesRecordArray" );
	
	return p;
}

/*
	FreeSpeciesRecordArray releases the memory which was previously
	allocated with NewSpeciesRecordArray.
*/

void FreeSpeciesRecordArray( SpeciesRecordPtr ptr )
{
	free( ptr );
}


/*
	NewThTrPropertiesArray allocates a one dimensional array of ThTrProperties.
*/

ThTrPropertiesPtr NewThTrPropertiesArray( int n )
{
	ThTrPropertiesPtr	p = NULL;
	
	if ( ( p = (ThTrPropertiesPtr)calloc( n, sizeof( ThTrProperties ) ) ) == NULL ) {
		FatalError( "Couldn't allocate memory for tht-propertiesArray!" );
	}
	
	return p;
}


/*
	FreeThTrPropertiesArray releases the momory which was previously allocated with 
	NewThTrPropertiesArray
*/

void FreeThTrPropertiesArray( ThTrPropertiesPtr ptr )
{
	free( ptr );
}


/*
	omega_mu() returns the ½_µ for a given dimensionless temperature t/(eps/k).
*/

static Double omega_mu( Double t )
{
	static Double m1 = 3.3530622607;
	static Double m2 = 2.53272006;
	static Double m3 = 2.9024238575;
	static Double m4 = 0.11186138893;
	static Double m5 = 0.8662326188;			/* = -0.1337673812 + 1.0 */
	static Double m6 = 1.3913958626;
	static Double m7 = 3.158490576;
	static Double m8 = 0.18973411754;
	static Double m9 = 0.00018682962894;
	Double num, den;

	num = m1 + t*(m2 + t*(m3 + t*m4));
	den = m5 + t*(m6 + t*(m7 + t*(m8 + t*m9)));
	return num / den;
}


/*
	ComputeThTrEnthalpy computes the enthalpy of the species pointed to by s
	at the temperature temp
*/
void ComputeThTrEnthalpy( Double temp, SpeciesRecordPtr s, ThTrPropertiesPtr t )
{
	Double		*coeff;						/* pointer to nasa coefficients */
	Double		r_over_m = kR / s->M;		/* in J / (kg K) */
	
	if ( temp > 1000.0 )	coeff = s->hot;
	else					coeff = s->cold;
	
	t->h = temp*(coeff[0]+temp*(coeff[1]/2.0+temp*(coeff[2]/3.0+
			temp*(coeff[3]/4.0+temp*coeff[4]/5.0))))+coeff[5];
	t->h *= r_over_m;
#ifdef DEBUG
	fprintf( stderr, "h = %e\n", t->h );
#endif
}

/*
	ComputeThTrProperties computes the properties cp, h, mu and
	lambda of the species pointed to by s at the temperature temp
*/
void ComputeThTrProperties( Double temp, SpeciesRecordPtr s, ThTrPropertiesPtr t )
{
	Double		*coeff;						/* pointer to nasa coefficients */
	Double		r_over_m = kR / s->M;		/* in J / (kg K) */
	
	if ( temp > 1000.0 )	coeff = s->hot;
	else					coeff = s->cold;
	
	t->cp = coeff[0]+temp*(coeff[1]+temp*(coeff[2]+temp*(coeff[3]+temp*coeff[4])));
	t->cp *= r_over_m;
	
	t->h = temp*(coeff[0]+temp*(coeff[1]/2.0+temp*(coeff[2]/3.0+
			temp*(coeff[3]/4.0+temp*coeff[4]/5.0))))+coeff[5];
	t->h *= r_over_m;
					
	t->mu = s->mu_coeff * sqrt( temp ) / omega_mu( temp * s->k_over_eps );

	t->lambda = t->mu * ( t->cp + 1.2 * r_over_m );

#ifdef DEBUG
	fprintf( stderr, "cp = %e\th = %e\nmu = %e\tlambda = %e\n",
		t->cp, t->h, t->mu, t->lambda );
#endif
}


/*
	omega_D() returns the Sto§integral ½_¶ for a given dimensionless 
	Temperature t/(eps/k)
*/

static Double omega_D( Double t )
{					
	static Double m1 = 6.8728271691;
	static Double m2 = 9.4122316321;
	static Double m3 = 7.7442359037;
	static Double m4 = 0.23424661229;
	static Double m5 = 1.45337701568;			/* = 1.0 + 0.45337701568 */
	static Double m6 = 5.2269794238;
	static Double m7 = 9.7108519575;
	static Double m8 = 0.46539437353;
	static Double m9 = 0.00041908394781;
	Double	num, den;
	
	num = m1 + t * (m2 + t * (m3 + t * m4));
	den = m5 + t * (m6 + t * (m7 + t * (m8 + t * m9)));
	return num / den;
}	 


/*
	Compute_D computes a mixture averaged diffusion coefficent for each species
*/

void Compute_D( Double temp, Double pressure, Double *X, Double *Y,
				SpeciesRecordPtr s, int nspecies, ThTrPropertiesPtr t )
{
	static Double	**inverse_Dij, zero = 0.0;
	static int		init = FALSE;
	int				j, i;
	Double			sum;
	
	if ( !init ) {
		init = TRUE;
		inverse_Dij = New2DArray( nspecies, nspecies );
		for ( i = 0; i < nspecies; ++i ) {
			for ( j = 0; j < nspecies; ++j ) inverse_Dij[i][j] = zero;
		}
	}

	for ( i = 0; i < nspecies; ++i ) {

		inverse_Dij[i][i] = zero;

		for (j = 0; j < i; ++j ) {

			inverse_Dij[j][i]  = pressure * omega_D( temp * s[i].omega_coeff[j] );
			inverse_Dij[j][i] /= s[i].D_coeff[j] * temp * sqrt( temp );
			inverse_Dij[i][j] = inverse_Dij[j][i];

#			ifdef DEBUG
			fprintf( stderr, "inverse_Dij = %g\t", inverse_Dij[i][j] );
#			endif
		}

#	ifdef DEBUG
	fprintf( stderr, "\n" );
#	endif
	}
	
	for ( i = 0; i < nspecies; ++i ) {

		sum = X[0] * inverse_Dij[i][0];		/* j = 0 */
		for ( j = 1; j < nspecies; ++j ) sum += X[j] * inverse_Dij[i][j];
		t[i].D = ( 1.0 - Y[i] ) / sum;

#		ifdef DEBUG
		fprintf( stderr, "D = %e\n", t[i].D);
#		endif
	}

}



void ComputeMuCoeff( SpeciesRecordPtr sr )
{
	/*
		Computes the constant part of the viscosity approximation
		
			µ = µ_coeff * ÃT / ½_µ				in  [kg / (m s)]
			µ_coeff = 2.6693e-6 ÃM / sigma^2
	*/
	sr->mu_coeff = 2.6693e-6 * sqrt( sr->M ) / (sr->sigma * sr->sigma);
}


/**go** old version of 'ComputeComposition' renamed to '_ComputeComposition' *
	void ComputeComposition( SpeciesRecordPtr sr  ) --- OLD VERSION ---
	
	Given a pointer to a SpeciesRecord sr, this functions parses the 
	name string of the species and counts the number of C, H, N, and
	O atoms. Other atoms are not considered. The algorithm currently
	works only for names that can be described by the following pro-
	ductions:
	
		name  -> list
		list  -> list desc | desc
		desc  -> id num | id
		id	  -> C | H | N | O
		num	  -> num digit | digit
		digit -> 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
*/

void _ComputeComposition( SpeciesRecordPtr sr  )
{
	CompositionPtr	cPtr = &sr->composition_elems;
	char			theChar, *ptr = sr->name;
	int				num_of_atoms;
	char			*isoPtr = NULL;
	Flag			err = FALSE;
	
	/* initialize everything to zero */
	cPtr->C = (unsigned char)0, cPtr->H = (unsigned char)0,
	cPtr->N = (unsigned char)0, cPtr->O = (unsigned char)0;
	
	if ( isoPtr = strrchr( ptr, '-' ) ) {
#ifdef DEBUGCOMPOSITION
		fprintf( stderr, "ptr = %s\t", ptr );
#endif
		ptr = &isoPtr[1];
#ifdef DEBUGCOMPOSITION
		fprintf( stderr, "ptr = %s\n", ptr );
#endif
	}
	while ( *ptr ) {
		
		theChar = *ptr++;	/* save this character and move to the next one */
		
		/*
			Get the number of atoms
		*/
		num_of_atoms = 0;
		while ( *ptr && isdigit(*ptr) ) {
			num_of_atoms *= 10;
			num_of_atoms += *ptr++ - '0';
		}
		if ( !num_of_atoms ) num_of_atoms = 1;	/* there is always at least one */
		
		switch ( toupper(theChar) ) {
			case 'C': cPtr->C += num_of_atoms; break;
			case 'H': cPtr->H += num_of_atoms; break;
			case 'N': cPtr->N += num_of_atoms; break;
			case 'O': cPtr->O += num_of_atoms; break;
			default: cPtr->C = cPtr->H = cPtr->N = cPtr->O = 0;
				err = TRUE;
				break;					/* don't do anything right now */
		}	/* switch */
		if ( err ) {
			break;
		}
	}	/* while */

#	ifdef DEBUG
	fprintf( stderr, "# \"%s\": ", sr->name );
	fprintf( stderr, "%d C atoms, %d H atoms, %d N atoms, %d O atoms\n",
		cPtr->C, cPtr->H, cPtr->N, cPtr->O );
#	endif
}

/**go** new version of 'ComputComposition', F and BR added, extended to 2 letter elements */
/*
	void ComputeComposition( SpeciesRecordPtr sr  )  --- NEW VERSION ---
	
	Given a pointer to a SpeciesRecord sr, this functions parses the 
	name string of the species and counts the number of C, H, N, BR, F, and
	O atoms. Other atoms are not considered. The algorithm currently
	works only for names that can be described by the following pro-
	ductions:
	
		name  -> list
		list  -> list desc | desc
		desc  -> id num | id
		id	  -> C | H | N | O | F | BR | AR
		num	  -> num digit | digit
		digit -> 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
*/

void ComputeComposition( SpeciesRecordPtr sr  )
{
	CompositionPtr	cPtr = &sr->composition_elems;
	char			theChar, *ptr = sr->name;
	int				num_of_atoms;
	char			*isoPtr = NULL;
	Flag			err = FALSE;
	
	/* initialize everything to zero */
	cPtr->C = (unsigned char)0, cPtr->H = (unsigned char)0,
	cPtr->N = (unsigned char)0, cPtr->O = (unsigned char)0,
	cPtr->F = (unsigned char)0, cPtr->BR = (unsigned char)0,
	cPtr->AR = (unsigned char)0, cPtr->HE = (unsigned char)0;/*PP*/
	
	if ( isoPtr = strrchr( ptr, '-' ) ) {
#ifdef DEBUGCOMPOSITION
		fprintf( stderr, "ptr = %s\t", ptr );
#endif
		ptr = &isoPtr[1];
#ifdef DEBUGCOMPOSITION
		fprintf( stderr, "ptr = %s\n", ptr );
#endif
	}
	while ( *ptr ) {
		
		theChar = *ptr++;	/* save this character and move to the next one */

/**go** if theChar is a 'B', it must be a BR atom, if theChar is a 'A', 
	it must be a AR atom (at least in this implementation), */
/**go** but let's take a look at the next letter */
		if ((theChar == 'A') ||(theChar == 'B')) {
/**go** if next letter is not 'R', I dunno know this element */ 
			if (*ptr != 'R') {
				cPtr->C = cPtr->H = cPtr->N = cPtr->O  = cPtr->F = cPtr->BR = cPtr->AR = 0;
				break;
			}
			else 
/**go** if it is 'R', goto next character */
			  ptr++;
		}
/**go** if it is 'H', species can be H or HE: check next character   */ /*PP*/
		else if (theChar == 'H') {
/**go** if next letter is 'E', goto next character */
		  if (*ptr == 'E') {
		    theChar = 'E'; /* Because H is already used in switch below*/
		    ptr++;
		  }
/**go** if next letter is not 'E' or one of known atoms, I dunno know this element */ 
		  else if ( (*ptr != 'C') && (*ptr != 'N') && (*ptr != 'F') && (*ptr != 'A') && (*ptr != 'H') && (*ptr != 'O') && (*ptr != 'B')){
		    cPtr->C = cPtr->H = cPtr->N = cPtr->O  = cPtr->F = cPtr->BR = cPtr->AR = 0;
		    break;
		  }
		}

		/*
			Get the number of atoms
		*/
		num_of_atoms = 0;
		while ( *ptr && isdigit(*ptr) ) {
			num_of_atoms *= 10;
			num_of_atoms += *ptr++ - '0';
		}
		if ( !num_of_atoms ) num_of_atoms = 1;	/* there is always at least one */
		
		switch ( toupper(theChar) ) {
			case 'C': cPtr->C  += num_of_atoms; break;
			case 'H': cPtr->H  += num_of_atoms; break;
			case 'N': cPtr->N  += num_of_atoms; break;
			case 'O': cPtr->O  += num_of_atoms; break;
			case 'F': cPtr->F  += num_of_atoms; break;
			case 'B': cPtr->BR += num_of_atoms; break;
			case 'A': cPtr->AR += num_of_atoms; break;
		        case 'E': cPtr->HE += num_of_atoms; break;
			default: cPtr->C = cPtr->H = cPtr->N = cPtr->O = cPtr->F = cPtr->BR = cPtr->AR = cPtr->HE = 0;
				err = TRUE;
				break;					/* don't do anything right now */
		}	/* switch */
		if ( err ) {
			break;
		}
	}	/* while */

#	ifdef DEBUG
	fprintf( stderr, "# \"%s\": ", sr->name );
	fprintf( stderr, "%d C atoms, %d H atoms, %d N atoms, %d O atoms, %d F atoms, %d BR, %d AR, %d HE atoms\n",
		cPtr->C, cPtr->H, cPtr->N, cPtr->O, cPtr->F, cPtr->BR , cPtr->AR, cPtr->HE );
#	endif
}



