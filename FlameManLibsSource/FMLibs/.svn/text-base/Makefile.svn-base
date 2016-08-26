COMPFLAGS     =  

CFLAGS        = $(COMPFLAGS) $(FlameManCCOpts) -c -D$(FlameManMach) -I$(FlameManLib)
CXXFLAGS      = $(COMPFLAGS) $(FlameManCXXOpts) -I$(FlameManLib) -D$(FlameManMach)
FFLAGS        = $(COMPFLAGS) $(FlameManF77Opts)

CC            = $(FlameManCC)
CXX           = $(FlameManCXX)
FF            = $(FlameManFF)


ARFLAGS       = cr

DEST	      = .

# Uncomment this if you want system header files to be expanded
#
# SYSHDRS       =

EXTHDRS	      =

HDRS	      = Interrupt.h \
		ListTool.h \
		Spline.h \
		MapMan.h \
		TofZ.h \
		BetaPDF.h \
		SmallNewton.h \
		paramtr.h

INSTALL	      = install

LIBNAME       = FMLib

LIBRARY	      = $(LIBNAME).a

LINTFLAGS     = -u $(CFLAGS) -o $(LIBNAME)

LINTLIB       = llib-l$(LIBNAME).ln

LINTLIBS      = 

MAKEFILE      = Makefile

OBJS	      = Interrupt.o \
		LTBody.o \
		LTHeader.o \
		LTUtilities.o \
		ReadStartProfile.o \
		Spline.o \
		MapMan.o \
		TofZ.o \
		BetaPDF.o \
		SmallNewton.o
		
#        Adiab_flam.o \
#        PDFF.o

PRINT	      = pr

PRINTFLAGS    =

LP	      = lp

LPFLAGS       = 

SHELL	      = /bin/sh

SRCS	      = Interrupt.c \
		LTBody.c \
		LTHeader.c \
		LTUtilities.c \
		ReadStartProfile.c \
		Spline.c \
		MapMan.C \
		TofZ.c \
		BetaPDF.c \
		SmallNewton.c

#        Adiab_flam.f \
#        PDFF.f

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS) $(MAKEFILE)
		ar $(ARFLAGS) $(LIBRARY) $(OBJS)
		$(FlameManRANLIB)
		@cp $(HDRS) $(FlameManLib)
		@mv $(LIBRARY) $(FlameManLib)

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
Interrupt.o: Interrupt.h
LTBody.o: ListTool.h
LTHeader.o: ListTool.h
LTUtilities.o: ListTool.h
ReadStartProfile.o: ListTool.h
Spline.o: Spline.h
TofZ.o: TofZ.h
BetaPDF.o: BetaPDF.h
SmallNewton.o: SmallNewton.h
MapMan.o: MapMan.h
$(OBJS) : $(MAKEFILE)

.SUFFIXES: $(.SUFFIXES) .cp .C .f

.cp.C:
	mv $< $*.C
.C.o:
	$(CXX) $(CXXFLAGS) -c $<
.f.o:
	$(FF) $(FFLAGS) -c $<
