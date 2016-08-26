COMPFLAGS     = #-O3  

FFLAGS        = $(COMPFLAGS) $(FlameManF77Opts)

FF            = $(FlameManFF)

HDRS	      = dassl.h

LIBRARY	      = ddasslLib.a

MAKEFILE	    = Makefile

OBJS	      = ddassl.o daux.o dlinpk.o

SRCS	      = ddassl.f daux.f dlinpk.f

$(LIBRARY):	$(OBJS)
		ar cr $(LIBRARY) $(OBJS)
		$(FlameManRANLIB)
		mv $(LIBRARY) $(FlameManLib)
		cp $(HDRS) $(FlameManLib)

$(OBJS) : $(HDRS) $(MAKEFILE) 

.SUFFIXES: $(.SUFFIXES) .f

.f.o:
	$(FF) $(FFLAGS) -c $<
