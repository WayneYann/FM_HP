CFLAGS	      = $(FlameManCCOpts) -c -D$(FlameManMach) -I$(FlameManLib)

DEST	      = .
CC            = $(FlameManCC)
EXTHDRS	      =

HDRS	      =

INSTALL	      = /etc/install

LD	      = $(CC)

LDFLAGS	      = $(SYMOPTS) -lm

LIBS	      = $(FlameManLib)/FLReader.a \
		$(FlameManLib)/list.a \
		$(FlameManLib)/libAM.a \
		$(FlameManLib)/alligator.a

#		/usr/local/lib/libdf.a \
#		/usr/local/lib/ptLib/WSSLib.a \

MAKEFILE      = Makefile

OBJS	      = ListTool.o 

PRINT	      = pr

PROGRAM       = ListTool

SHELL	      = /bin/sh

SRCS	      = ListTool.c

SYSHDRS	      =

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
		mv $(PROGRAM) $(FlameManLib)
		chmod 755 $(FlameManLib)/$(PROGRAM)

clean:;		@rm -f $(OBJS) 

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

depend:;	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)

echo:;		@echo $(HDRS) $(SRCS)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		@echo Installing $(PROGRAM) in $(DEST)
		@-strip $(PROGRAM)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(PROGRAM); $(INSTALL) -f $(DEST) $(PROGRAM)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

$(OBJS) : $(HDRS) $(MAKEFILE) 

