DBG = -DqMemDebug

SYMOPTS = -O

CFLAGS = $(DBG) $(SYMOPTS) $(COPTS)

OBJECTS = \
	FLScan.tab.o \
	FLScan.main.o \
	lex.fl.o \
	Stack.o \
	regex.o \
	xmalloc.o \
	alloca.o 

MAKEFILE = debug.Lib.make

TARGET = FLReader.dbg.a


$(TARGET): $(MAKEFILE) $(OBJECTS)
	ar cru $(TARGET) $(OBJECTS) 
	mv $(TARGET) $(myLibs)
	cp FLReader.h $(myLibs)
	rm *.o lex.fl.c *.tab.[ch]
	FixPermissions

FLScan.main.o : $(MAKEFILE) FLScan.main.c FLScan.h
	cc -c FLScan.main.c $(CFLAGS)

FLScan.tab.c\
FLScan.tab.h: FLScan.y $(MAKEFILE) FLScan.h
	bison -dp fl FLScan.y

FLScan.tab.o : $(MAKEFILE) FLScan.tab.c FLScan.h 
	cc -c FLScan.tab.c $(CFLAGS)

lex.fl.c: FLScan.l $(MAKEFILE) FLScan.tab.h FLScan.h
	flex -s -Pfl FLScan.l

lex.fl.o : $(MAKEFILE) FLScan.l FLScan.h
	cc -c lex.fl.c $(CFLAGS)

Stack.o : $(MAKEFILE) Stack.c Stack.h
	cc -c Stack.c $(CFLAGS)

regex.o : $(MAKEFILE) regex.c regex.h
	cc -c regex.c $(CFLAGS)

xmalloc.o : $(MAKEFILE) xmalloc.c
	cc -c xmalloc.c $(CFLAGS)

alloca.o : $(MAKEFILE) alloca.c
	cc -c alloca.c $(CFLAGS)

clean:;	@rm $(OBJS) lex.fl.* *.tab.*
