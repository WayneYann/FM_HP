#define MECHANISM "CH4.Igni73"
#include "FlameMaster.h"

/*	Mechanism file: "CH4.Igni73.mech"	*/

typedef Double (*BFFunction)(Double T);

/* prototypes */
Double Fc34( Double T );
Double Fc36( Double T );
Double Fc51( Double T );
Double Fc58( Double T );
Double FcErr( Double T );


extern BFFunction gBroadening[4];
