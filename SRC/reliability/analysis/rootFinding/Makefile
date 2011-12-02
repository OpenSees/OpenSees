include ../../../../Makefile.def

OBJS       = RootFinding.o \
	SecantRootFinding.o \
	SafeGuardedZeroFindingAlgorithm.o \
	ZeroFindingAlgorithm.o



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
