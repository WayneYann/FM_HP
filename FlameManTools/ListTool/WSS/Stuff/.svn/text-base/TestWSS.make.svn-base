CFLAGS	      = $(COPTS) -O -DqMemDebug -z

DEST	      = .

EXTHDRS	      =

HDRS	      =

INSTALL	      = /etc/install

LD	      = cc

LDFLAGS	      =

LIBS	      = \
		/jgLib/alligator.a \
		/lib/libm.a

MAKEFILE      = TestWSS.make

OBJS	      = TestWSS.o \
		WSS.o

PRINT	      = pr

PROGRAM       = TestWSS

SHELL	      = /bin/sh

SRCS	      = TestWSS.c \
		WSS.c

SYSHDRS	      =

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS) tags
		$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

clean:;		@rm -f $(OBJS) core

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

$(DEST)/$(PROGRAM): $(SRCS) $(LIBS) $(HDRS) $(EXTHDRS)
		@$(MAKE) -f $(MAKEFILE) ROOT=$(ROOT) DEST=$(DEST) install


$(OBJS) : $(HDRS) $(MAKEFILE) 

