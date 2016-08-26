COMPFLAGS     =   

OBJECTS = alligator.o
COpts = $(COMPFLAGS) $(FlameManCCOpts) -c -D$(FlameManMach)
CC = $(FlameManCC)
LIBRARY = alligator.a

$(LIBRARY): $(OBJECTS) Makefile
	ar cr $(LIBRARY) $(OBJECTS)
	$(FlameManRANLIB)
	mv $(LIBRARY) $(FlameManLib)
	cp alligator.h $(FlameManLib)
alligator.o:
	$(CC) $(COpts) alligator.c
