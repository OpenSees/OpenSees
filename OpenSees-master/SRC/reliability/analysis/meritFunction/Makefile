include ../../../../Makefile.def

OBJS       = 	AdkZhangMeritFunctionCheck.o \
	CriteriaReductionMeritFunctionCheck.o \
	MeritFunctionCheck.o


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
