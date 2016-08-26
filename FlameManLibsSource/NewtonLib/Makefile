COMPFLAGS     =  
CXXFLAGS      = $(COMPFLAGS) $(FlameManCXXOpts) -I$(FlameManLib) -D$(FlameManMach)
 
ARFLAGS       = cr

DEST	      = .

CXX           = $(FlameManCXX)

# Uncomment this if you want system header files to be expanded
#
# SYSHDRS       =

EXTHDRS	      =

libLoc        = $(FlameManLib)
HDRS	      = Newton.h

INSTALL	      = install

LIBNAME       = NewtonLib

LIBRARY	      = $(LIBNAME).a

LINTFLAGS     = -u $(CFLAGS) -o $(LIBNAME)

LINTLIB       = llib-l$(LIBNAME).ln

LINTLIBS      = 

MAKEFILE      = Makefile

OBJS	      = Newton.o \
		NewtonUt.o \
		AdaptiveGrid.o

PRINT	      = pr

PRINTFLAGS    =

LP	      = lp

LPFLAGS       = 

SHELL	      = /bin/sh

SRCS	      = Newton.C \
		NewtonUt.C \
		AdaptiveGrid.C

all:		$(LIBRARY) 
			@echo 

$(LIBRARY):	$(OBJS) $(MAKEFILE)
		ar $(ARFLAGS) $(LIBRARY) $(OBJS)
		$(FlameManRANLIB)
		@cp $(HDRS) $(libLoc)
		@mv $(LIBRARY) $(libLoc)

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(LIBRARY) core tags

depend:;	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)

echo:;		@echo $(HDRS) $(SRCS)

extract:;	@ar x $(DEST)/$(LIBRARY)

-1:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(LIBRARY)
		@echo Installing $(LIBRARY) in $(DEST)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(LIBRARY); $(INSTALL) -f $(DEST) $(LIBRARY)); fi

lint:		$(LINTLIB)

$(LINTLIB):	$(LINTLIBS) $(SRCS) $(HDRS) $(EXTHDRS)
		lint $(LINTFLAGS) $(LINTLIBS) $(SRCS)

print:;		@$(PRINT) $(PRINTFLAGS) $(HDRS) $(SRCS) | $(LP) $(LPFLAGS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:         $(DEST)/$(LIBRARY)

touch:;		@touch $(SRCS) $(HDRS)
###
Newton.o: Newton.h $(MAKEFILE)
NewtonUt.o: Newton.h $(MAKEFILE)
AdaptiveGrid.o: Newton.h $(MAKEFILE)

.SUFFIXES: $(.SUFFIXES) .cp .C

.cp.C:
	mv $< $*.C
.C.o:
	$(CXX) $(CXXFLAGS) -c $<

