#define MECHANISM "nHeptane.allstarnew_oks"
#include "FlameMaster.h"

/*	Mechanism file: "nHeptane.allstarnew_oks.mech"	*/

typedef Double (*BFFunction)(Double T);

/* prototypes */
Double Fca34( Double T );
Double Fc36( Double T );
Double Fca51( Double T );
Double Fca58( Double T );
Double FcErr( Double T );


extern BFFunction gBroadening[4];
