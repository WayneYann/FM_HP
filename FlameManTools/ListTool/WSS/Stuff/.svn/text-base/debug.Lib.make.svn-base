CFLAGS	      = $(COPTS) -O -DqMemDebug -z

DEST	      = .

EXTHDRS	      =

HDRS	      = WSS.h

INSTALL	      = /etc/install

LIBRARY	      = WSSLib.dbg.a

MAKEFILE      = debug.Lib.make

OBJS	      = \
		WSS.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = \
		WSS.c

SYSHDRS	      =

all:		$(LIBRARY)

$(LIBRARY):	$(OBJS)
		ar cru $(LIBRARY) $(OBJS)
		mv $(LIBRARY) $(myLibs) 
		cp $(HDRS) $(myLibs)
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

$(DEST)/$(LIBRARY): $(SRCS) $(HDRS) $(EXTHDRS)
		@$(MAKE) -f $(MAKEFILE) ROOT=$(ROOT) DEST=$(DEST) install


$(OBJS) : $(HDRS) $(MAKEFILE) 

