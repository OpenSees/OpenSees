include ../../../../Makefile.def

OBJS       = SingleDomEleIter.o SingleDomNodIter.o \
	SingleDomSP_Iter.o SingleDomMP_Iter.o SingleDomAllSP_Iter.o \
	SingleDomParamIter.o SingleDomPC_Iter.o

# Compilation control

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
