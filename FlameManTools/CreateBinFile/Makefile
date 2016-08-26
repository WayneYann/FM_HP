COMPFLAGS     = 

CFLAGS	      = $(FlameManCCOpts) -c -D$(FlameManMach) -I$(FlameManLib) $(COMPFLAGS)

CC = $(FlameManCC)

DEST	      = .

# Uncomment this if you want system header files to be expanded
#
# SYSHDRS       =

EXTHDRS	      =

HDRS	      = ThermoProperties.h

INSTALL	      = install

LD	      = $(CC)

LDFLAGS	      = $(COMPFLAGS) $(FlameManLDCCOpts) -lm

LIBS	      = $(FlameManLib)/libAM.a $(FlameManLib)/alligator.a

LINTLIBS      =

LINTFLAGS     = -u $(CFLAGS)

MAKEFILE      = Makefile

OBJS	      = CreateBinFile.o \
		ThermoProperties.o

PRINT	      = pr

PRINTFLAGS    =

LP	      = lp

LPFLAGS       = 

PROGRAM       = CreateBinFile

SHELL	      = /bin/sh

SRCS	      = CreateBinFile.c \
		ThermoProperties.c

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

touch:;		@touch $(SRCS) $(HDRS)

###
CreateBinFile.o: ThermoProperties.h
ThermoProperties.o: ThermoProperties.h

