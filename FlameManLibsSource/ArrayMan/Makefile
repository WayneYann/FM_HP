COMPFLAGS     =  

OBJECTS = AM_alloc.o \
			ArrayManager.o \
			AM_blas.o \
			AM_Printing.o \
			AM_Gauss.o \
			AM_BlockTriDiagSolver.o

COpts = $(COMPFLAGS) $(FlameManCCOpts) -c -I$(FlameManLib) -D$(FlameManMach)
CC = $(FlameManCC)

LIBRARY	      = libAM.a

$(LIBRARY): $(OBJECTS) Makefile
	ar cr $(LIBRARY) $(OBJECTS)
	$(FlameManRANLIB)
	mv $(LIBRARY) $(FlameManLib)
	cp ArrayManager.h jgTypes.h $(FlameManLib)
ArrayManager.o:
	$(CC) $(COpts) ArrayManager.c
AM_Printing.o:
	$(CC) $(COpts) AM_Printing.c
AM_Gauss.o:
	$(CC) $(COpts) AM_Gauss.c
AM_BlockTriDiagSolver.o:
	$(CC) $(COpts) AM_BlockTriDiagSolver.c
AM_blas.o:
	$(CC) $(COpts) AM_blas.c
AM_alloc.o:
	$(CC) $(COpts) AM_alloc.c
