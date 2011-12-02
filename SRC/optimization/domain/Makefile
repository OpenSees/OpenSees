include ../../../Makefile.def

OBJS       = Optimization.o \
	OptimizationDomainComponent.o \
	OptimizationDomain.o

################### TARGETS ########################
all: $(OBJS)
	@$(CD) $(FE)/optimization/domain/component; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core test

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core test
	@$(RM) $(RMFLAGS) test.\$\$\$ 
	@$(RM) $(RMFLAGS) *.ti 
	@$(RM) $(RMFLAGS) tca.map 
	@$(RM) $(RMFLAGS) .inslog*

spotless: clean

wipe: spotless
	@$(CD) $(FE)/optimization/domain/component; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
