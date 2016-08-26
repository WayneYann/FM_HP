#
# makefile for the list package
#

# -g for debugging, -O for optimization
CC            = $(FlameManCC)
CFLAGS	      = $(FlameManCCOpts) -c -D$(FlameManMach) -I$(FlameManLib)
dbg  = #.dbg
defs = #-DqMemDebug
Sources = List.c
Objects = $(Sources:.c=.o)
Headers = List.h
LIBRARY	      = list$(dbg).a

install : $(Objects)     
	ar cr $(LIBRARY) $(Objects)
	$(FlameManRANLIB)
	$(FlameManRANLIB)
	cp $(Headers) $(FlameManLib)
	mv $(LIBRARY) $(FlameManLib)
