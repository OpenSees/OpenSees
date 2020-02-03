include ../../../../Makefile.def

OBJS       = 	UnloadingRule.o \
	ConstantUnloadingRule.o EnergyUnloadingRule.o \
	KarsanUnloadingRule.o TakedaUnloadingRule.o \
	TclModelBuilderUnloadingRuleCommand.o

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
