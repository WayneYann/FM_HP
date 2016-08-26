DBG = -UqMemDebug 

CFLAGS	      = $(FlameManCCOpts) -c -D$(FlameManMach) -I$(FlameManLib) $(DBG)

CC            = $(FlameManCC)
YACC          = $(FlameManYACC)

OBJECTS = \
	FLScan.tab.o \
	FLScan.main.o \
	lex.fl.o \
	Stack.o \
	regex.o \
	xmalloc.o \
	alloca.o 

MAKEFILE = Makefile

LIBRARY = FLReader.a


$(LIBRARY): $(MAKEFILE) $(OBJECTS)
	ar cru $(LIBRARY) $(OBJECTS) 
	$(FlameManRANLIB)
	mv $(LIBRARY) $(FlameManLib)
	cp FLReader.h $(FlameManLib)
	rm *.o lex.fl.c *.tab.[ch]

FLScan.main.o : $(MAKEFILE) FLScan.main.c FLScan.h
	$(CC) -c FLScan.main.c $(CFLAGS)

FLScan.tab.c FLScan.tab.h: FLScan.y $(MAKEFILE) FLScan.h
	$(YACC) -p fl -dl -b FLScan FLScan.y
	sed -n -e '/#define	yy/p' FLScan.tab.c > tmp
	sed -n -e '/extern/!p' FLScan.tab.h >> tmp
	sed -n -e '/extern/p' FLScan.tab.h >> tmp
	mv tmp FLScan.tab.h

FLScan.tab.o : $(MAKEFILE) FLScan.tab.c FLScan.h 
	$(CC) -c FLScan.tab.c $(CFLAGS)

lex.fl.c: FLScan.l $(MAKEFILE) FLScan.tab.h FLScan.h
	$(FLEXPATH)flex -s -Pfl -t FLScan.l > lex.fl.c

lex.fl.o : $(MAKEFILE) FLScan.l FLScan.h lex.fl.c
	$(CC) -c lex.fl.c $(CFLAGS)

Stack.o : $(MAKEFILE) Stack.c Stack.h
	$(CC) -c Stack.c $(CFLAGS)

regex.o : $(MAKEFILE) regex.c regex.h
	$(CC) -c regex.c $(CFLAGS)

xmalloc.o : $(MAKEFILE) xmalloc.c
	$(CC) -c xmalloc.c $(CFLAGS)

alloca.o : $(MAKEFILE) alloca.c
	$(CC) -c alloca.c $(CFLAGS)

clean:;	@rm $(OBJS) lex.fl.* *.tab.*
