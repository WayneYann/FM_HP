#ifndef __MAGICFILE__
#define __MAGICFILE__
#define MECHANISM "H2.allstar3"

typedef enum SpeciesLabel {

	/* Computed species s.. */
	/* Steady-state species ss.. */
	sN2 = 0,
	sO2 = 1,
	sH = 2,
	sOH = 3,
	sO = 4,
	sH2 = 5,
	sH2O = 6,
	sHO2 = 7,
	sH2O2 = 8,
	sEnd
} SpeciesLabel;


typedef enum ReactionLabel {
	/* Reactions */
	r1f = 0,
	r1b = 1,
	r2f = 2,
	r2b = 3,
	r3f = 4,
	r3b = 5,
	r4f = 6,
	r4b = 7,
	r5f = 8,
	r5b = 9,
	r6f = 10,
	r6b = 11,
	r7f = 12,
	r7b = 13,
	r8f = 14,
	r8b = 15,
	r9f = 16,
	r9b = 17,
	r10f = 18,
	r10b = 19,
	r11f = 20,
	r11b = 21,
	r12f = 22,
	r12b = 23,
	r13f = 24,
	r13b = 25,
	r14f = 26,
	r14b = 27,
	r15f = 28,
	r15b = 29,
	r16f = 30,
	r16b = 31,
	r17f = 32,
	r17b = 33,
	r18f = 34,
	r18b = 35,
	r19f = 36,
	r19b = 37,
	/* PAHReactions */
	/* SootReactions */
	rEnd
} ReactionLabel;


typedef enum TirdBodyLabel {

	mM1 = 0,
	mEnd
} TirdBodyLabel;

#endif
