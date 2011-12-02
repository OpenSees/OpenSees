include ../../Makefile.def

OBJS       = commands.o TclFeViewer.o TclVideoPlayer.o

all: $(OBJS)

tcl: $(OBJS)  tclMain.o tclAppInit.o

tk: $(OBJS)  tkMain.o tkAppInit.o

FEAP_ELE_LIBS = $(FE)/element/feap/elmt01.o \
	$(FE)/element/feap/elmt02.o \
	$(FE)/element/feap/elmt03.o \
	$(FE)/element/feap/elmt04.o \
	$(FE)/element/feap/elmt05.o \
	$(FE)/element/feap/fillCommon.o \
	$(FE)/element/feap/getCommon.o \
	$(FE)/element/feap/dummyFeap.o \
	$(FE)/element/feap/fElement.o \
	$(FE)/element/feap/fElmt02.o 

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.





