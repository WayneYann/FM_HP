COMPFLAGS     =   
CFLAGS	      = $(FlameManCCOpts) -c -D$(FlameManMach) -I$(FlameManLib) $(COMPFLAGS)
PFLAGS	      = $(COMPFLAGS)
FFLAGS	      = $(COMPFLAGS)
CCFLAGS       = $(COMPFLAGS)
#COpts = -ansi -c -I$(FlameManLib) -DHP

CXXFLAGS      = $(COMPFLAGS)

DEST	      = .

CC            = $(FlameManCC)
YACC          = $(FlameManYACC)

# Uncomment this if you want system header files to be expanded
#
# SYSHDRS       =

EXTHDRS	      =

HDRS	      = ListMan.h \
		ReactionParser.h \
		ReactionScan.tab.h \
		ReadThermoProps.h \
		ScanMan.h \
		Redux.h \
		Mathematica.h \
		svd.h \
		FitPak.h

INSTALL	      = install

LD	      = $(CC)

LDFLAGS	      = $(COMPFLAGS) $(FlameManLDCCOpts) -lm

LIBS	      = $(FlameManLib)/libAM.a $(FlameManLib)/alligator.a

LINTLIBS      =

LINTFLAGS     = -u $(CFLAGS)

MAKEFILE      = Makefile

OBJS	      = CheckUtilities.o \
		ListMan.o \
		ListUtilities1.o \
		ListUtilities2.o \
		ParseUtilities.o \
		ReactionScan.tab.o \
		ScanMan.o \
		FitPakLinearSystems.o \
		FitPakFits.o \
		FitPakUtil.o \
		ScanUtilities.o \
		Labels.o \
		Redux.o \
		Mathematica.o \
		svd.o \
		Matrix.o \
		IndependentEQ.o \
		lex.yy.o

PRINT	      = pr

PRINTFLAGS    =

LP	      = lp

LPFLAGS       = 

PROGRAM       = ScanMan

SHELL	      = /bin/sh

SRCS	      = CheckUtilities.c \
		ListMan.c \
		ListUtilities1.c \
		ListUtilities2.c \
		ParseUtilities.c \
		ReactionScan.y \
		ReactionScan.tab.c \
		ScanMan.c \
		FitPakLinearSystems.c \
		FitPakFits.c \
		FitPakUtil.c \
		ScanUtilities.c \
		Labels.c \
		Redux.c \
		Mathematica.c \
		svd.c \
		Matrix.c \
		IndependentEQ.c \
		lex.yy.c

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS) $(MAKEFILE)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(OBJS) $(LIBS) $(LDFLAGS) -o $(PROGRAM)
		cp $(PROGRAM) $(FlameManLib)
		@echo "done"

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

depend:;	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)

echo:;		@echo $(HDRS) $(SRCS)

-1:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi

lint:		$(LINTLIBS) $(SRCS) $(HDRS) $(EXTHDRS)
		lint $(LINTFLAGS) $(LINTLIBS) $(SRCS)

print:;		@$(PRINT) $(PRINTFLAGS) $(HDRS) $(SRCS) | $(LP) $(LPFLAGS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

touch:;		@touch $(SRCS) $(HDRS); rm lex.yy.c ReactionScan.tab.[ch]
###
lex.yy.c: ReactionScan.flex $(MAKEFILE)
			$(FLEXPATH)flex -8 -s ReactionScan.flex
ReactionScan.tab.c: ReactionScan.y $(MAKEFILE) ReactionParser.h ScanMan.h ListMan.h
			$(YACC) -dt ReactionScan.y
			mv y.tab.c ReactionScan.tab.c
			mv y.tab.h ReactionScan.tab.h
CheckUtilities.o: ScanMan.h ListMan.h $(MAKEFILE)
#CheckUtilities.o: ScanMan.h ListMan.h $(MAKEFILE)
#	$(CC)  $(COpts) CheckUtilities.c
ListMan.o: ListMan.h $(MAKEFILE)
ListUtilities1.o: ScanMan.h ListMan.h $(MAKEFILE)
ListUtilities2.o: ScanMan.h ListMan.h $(MAKEFILE)
ParseUtilities.o: ReactionParser.h ReactionScan.tab.c ScanMan.h ListMan.h $(MAKEFILE)
ReactionScan.o: ReactionParser.h ScanMan.h ListMan.h $(MAKEFILE)
ScanMan.o: ScanMan.h ListMan.h ReadThermoProps.h $(MAKEFILE)
FitPakLinearSystems.o: FitPak.h $(MAKEFILE)
FitPakFits.o: FitPak.h $(MAKEFILE)
FitPakUtil.o: FitPak.h $(MAKEFILE)
ScanUtilities.o: ScanMan.h ListMan.h $(MAKEFILE)
Labels.o: ScanMan.h ListMan.h $(MAKEFILE)
lex.yy.o: ReactionParser.h ScanMan.h ListMan.h $(MAKEFILE)
