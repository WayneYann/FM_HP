#define MECHANISM "gri.300"
#include "FlameMaster.h"

/*	Mechanism file: "gri.300.mech"	*/

typedef Double (*BFFunction)(Double T);

/* prototypes */
Double Fc12( Double T );
Double Fc50( Double T );
Double Fc52( Double T );
Double Fc54( Double T );
Double Fc56( Double T );
Double Fc57( Double T );
Double Fc59( Double T );
Double Fc63( Double T );
Double Fc70( Double T );
Double Fc71( Double T );
Double Fc72( Double T );
Double Fc74( Double T );
Double Fc76( Double T );
Double Fc83( Double T );
Double Fc85( Double T );
Double Fc95( Double T );
Double Fc131( Double T );
Double Fc140( Double T );
Double Fc147( Double T );
Double Fc158( Double T );
Double Fc174( Double T );
Double Fc185( Double T );
Double Fc237( Double T );
Double Fc241( Double T );
Double Fc289( Double T );
Double Fc304( Double T );
Double Fc312( Double T );
Double Fc318( Double T );
Double Fc320( Double T );
Double FcErr( Double T );


extern BFFunction gBroadening[29];
