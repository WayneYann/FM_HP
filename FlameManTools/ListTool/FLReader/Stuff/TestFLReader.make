DBG = -DqMemDebug

SYMOPTS = -O

COPTS = $(DBG) $(SYMOPTS) $(COPTS)

OBJECTS = TestFLReader.o

LIBS =	$(myLibs)/FLReader.dbg.a \
		/usr/local/lib/jgLib/list.dbg.a \
		/usr/local/lib/jgLib/libAM.a \
		/usr/local/lib/jgLib/alligator.a \
		/lib/libm.a

MAKEFILE = TestFLReader.make

TARGET = TestFLReader

$(TARGET): $(MAKEFILE) $(OBJECTS)
	cc -o $(TARGET) $(SYMOPTS) $(OBJECTS) $(LIBS)

TestFLReader.o : $(MAKEFILE) TestFLReader.c
	cc -c TestFLReader.c $(COPTS)

