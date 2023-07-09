include ../../../../Makefile.def

OBJS       = 	PerformanceFunction.o PerformanceFunctionIter.o \
	LimitStateFunction.o LimitStateFunctionIter.o \
	ConstraintFunction.o

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
