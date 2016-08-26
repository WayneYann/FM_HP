COMPFLAGS     =

CFLAGS	      = $(FlameManCCOpts) -c -D$(FlameManMach) -I$(FlameManLib)

DEST	      = .

EXTHDRS	      =

HDRS	      = WSS.h

INSTALL	      = /etc/install

LIBRARY	      = WSSLib.a

MAKEFILE      = Makefile

OBJS	      = \
		WSS.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = \
		WSS.c

CC            = $(FlameManCC)

SYSHDRS	      =

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS)
		ar cru $(LIBRARY) $(OBJS)
		$(FlameManRANLIB)
		mv $(LIBRARY) $(FlameManLib)
		cp $(HDRS) $(FlameManLib)
		rm $(OBJS)

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(LIBRARY) core tags

depend:;	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)

echo:;		@echo $(HDRS) $(SRCS)

extract:;	@ar x $(DEST)/$(LIBRARY)

index:;		@ctags -wx $(HDRS) $(SRCS)

install:	$(LIBRARY)
		@echo Installing $(LIBRARY) in $(DEST)
		@if [ $(DEST) != . ]; then \
		(rm -f $(DEST)/$(LIBRARY); $(INSTALL) -f $(DEST) $(LIBRARY)); fi

print:;		@$(PRINT) $(HDRS) $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags $(HDRS) $(SRCS)

update:         $(DEST)/$(LIBRARY)

$(OBJS) : $(HDRS) $(MAKEFILE) 

