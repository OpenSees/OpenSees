include ../../Makefile.def

ifndef OPTIMIZATION_LIBRARY
  OPTIMIZATION_LIBRARY = $(HOME)/lib/libOptimization.a
endif

OPTIMIZATION_OBJS = $(FE)/optimization/domain/Optimization.o \
	$(FE)/optimization/domain/OptimizationDomain.o \
	$(FE)/optimization/domain/OptimizationDomainComponent.o \
	$(FE)/optimization/domain/component/ConstraintFunction.o \
	$(FE)/optimization/domain/component/DesignVariable.o \
	$(FE)/optimization/domain/component/DesignVariableIter.o \
	$(FE)/optimization/domain/component/DesignVariablePositioner.o \
	$(FE)/optimization/domain/component/DesignVariablePositionerIter.o \
	$(FE)/optimization/domain/component/ObjectiveFunction.o \
	$(FE)/optimization/tcl/TclOptimizationBuilder.o


all:
	@$(CD) $(FE)/optimization/tcl; $(MAKE);
	@$(CD) $(FE)/optimization/domain; $(MAKE);
	@$(CD) $(FE)/optimization/SNOPT; $(MAKE);
	$(AR) $(ARFLAGS) $(OPTIMIZATION_LIBRARY) $(OPTIMIZATION_OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core test

clean: tidy

spotless: clean

wipe: spotless
	@$(CD) $(FE)/optimization/tcl; $(MAKE) wipe;
	@$(CD) $(FE)/optimization/domain; $(MAKE) wipe;
	@$(CD) $(FE)/optimization/SNOPT; $(MAKE) wipe;
	@$(RM) $(RMFLAGS) $(OPTIMIZATION_LIBRARY)

# DO NOT DELETE THIS LINE -- make depend depends on it.
