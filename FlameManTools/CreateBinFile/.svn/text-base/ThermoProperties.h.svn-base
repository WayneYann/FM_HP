/* 
	Thermoproperties.h: header file for the thermodynamical properties
						functions.
	
	Revision history:	7/19/91	first version by pt
						7/30/91 cleaned up (jg)
   changed by go 12/20/94: BR and F added
	 */

#ifndef __ThermoProperties__
#define __ThermoProperties__

#ifndef __ArrayManager__
#include "ArrayManager.h"
#endif

#ifndef FALSE
#define	FALSE	0
#define TRUE	1
#endif

typedef char Flag;

#define kR			(8314.4)		/* universal gas constant in J/(kmol K)*/
#define kNameLen	32				/* length identifiers for chemical species */

/**go** F and BR added to struct 'Composition' */
typedef struct Composition {
  unsigned char	C, H, N, O, BR, F, AR, HE;	
} Composition, *CompositionPtr;

typedef struct SpeciesRecord {
	char		name[kNameLen];		/* name of component  (input) */
	Double		M;					/* molecular weight in kg/kmol  (computed) */
	Composition	composition_elems;	/* element composition record  (computed) */
	Double		hot[7];				/* nasa coefficients for T>1000 K  (input) */
	Double		cold[7];			/* nasa coefficients for T<1000 K (input) */
	Double		k_over_eps;			/* temperature factor for evaluating sigma (= k/eps);
									   input is eps/k  (computed) */
	Double		sigma;
	Double		mu_coeff;			/* constant factor for evaluation of mu in kg/(m*s)
									   (= 2.6693e-6 * sqrt(M)/sigma)  (computed) */

	Double		*D_coeff;			/* constant factor for evaluation of D in m^2/s 
									   = 0.0188292*sqrt(1/Mi+1/Mj)/(.5*(sigma_i+sigma_j))^2
									   (computed)	*/
	Double		*omega_coeff;		/* constant factor for evaluation of omega_D 
									   (= sqrt(k_over_eps_i * k_over_eps_j))  (computed) */
	int			id;					/* not used by preprocessor (set to -1) */
} SpeciesRecord, *SpeciesRecordPtr;


typedef struct ThTrProperties {			/* THermodynamic and TRansport properties */
	Double		mu;						/* dynamic viscosity in kg/(m*s) */
	Double		lambda;					/* thermal conductivity in W/(m*K) */
	Double		D;						/* diffusion coefficient in m^2/s */
	Double		cp;						/* specific heat in J/(kg K) */
	Double		h;						/* specific enthalpy in J/kg */
} ThTrProperties, *ThTrPropertiesPtr;


/*
	Prototypes
*/

SpeciesRecordPtr NewSpeciesRecordArray( int n );
void FreeSpeciesRecordArray( SpeciesRecordPtr ptr );
ThTrPropertiesPtr NewThTrPropertiesArray( int n );
void FreeThTrPropertiesArray( ThTrPropertiesPtr ptr );
void ComputeThTrEnthalpy( Double temp, SpeciesRecordPtr s, ThTrPropertiesPtr t );
void ComputeThTrProperties( Double temp, SpeciesRecordPtr s, ThTrPropertiesPtr t );
void Compute_D( Double temp, Double pressure, Double *X, Double *Y,
				SpeciesRecordPtr s, int nspecies, ThTrPropertiesPtr t );
void ComputeMuCoeff( SpeciesRecordPtr sr );
void ComputeComposition( SpeciesRecordPtr sr  );
void UpperString( char *string );

extern void FatalError( char *errorString );


#endif	/* __ThermoProperties__ */
