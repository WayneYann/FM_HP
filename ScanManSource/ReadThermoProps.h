/*
 *	ScanMan:	Written by Josef Goettgens and Heinz Pitsch
 *				Copyright 1991-92, Josef Goettgens and Heinz Pitsch
 *				All rights reserved.
 */

typedef struct Composition {
  unsigned char	C, H, N, O, BR, F, AR, HE;	
} Composition, *CompositionPtr;

typedef struct SpeciesRecord {
	char		name[32];		/* name of component  (input) */
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

typedef struct SpeciesRecordHigh {
	SpeciesRecordPtr 	speciesRecordPtr;
	int					len;
} SpeciesRecordHigh;


SpeciesRecordPtr NewSpeciesRecordArray( int n );

void FreeSpeciesRecordArray( SpeciesRecordPtr ptr );

Flag MatchIsomereName( LinkPtr link, SpeciesRecordHigh	*recordStruct );

void CopySpecies( SpeciesPtr species, SpeciesRecordPtr record );
