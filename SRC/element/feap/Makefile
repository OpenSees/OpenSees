include ../../../Makefile.def

OBJS       = elmt01.o elmt02.o elmt03.o elmt04.o elmt05.o \
	fillCommon.o fElement.o fElmt02.o dummyFeap.o getCommon.o \
	TclFeapElementCommand.o

all:         $(OBJS)


# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
