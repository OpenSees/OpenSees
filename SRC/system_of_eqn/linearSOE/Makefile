include ../../../Makefile.def

OBJS       = LinearSOE.o DomainSolver.o LinearSOESolver.o

ifeq ($(PROGRAMMING_MODE), PARALLEL)

all:         $(OBJS)
	@$(CD) $(FE)/system_of_eqn/linearSOE/bandSPD; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/bandGEN; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/fullGEN; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/profileSPD; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/sparseGEN; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/sparseSYM; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/sparseSYM; $(MAKE) law;
	@$(CD) $(FE)/system_of_eqn/linearSOE/umfGEN; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/cg; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/BJsolvers; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/diagonal; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/petsc; $(MAKE);

else

all:         $(OBJS)
	@$(CD) $(FE)/system_of_eqn/linearSOE/bandSPD; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/bandGEN; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/fullGEN; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/profileSPD; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/sparseGEN; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/sparseSYM; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/sparseSYM; $(MAKE) law;
	@$(CD) $(FE)/system_of_eqn/linearSOE/umfGEN; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/cg; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/BJsolvers; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/linearSOE/diagonal; $(MAKE);

endif



# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/system_of_eqn/linearSOE/bandSPD; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn/linearSOE/bandGEN; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn/linearSOE/fullGEN; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn/linearSOE/profileSPD; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn/linearSOE/sparseGEN; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn/linearSOE/sparseSYM; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn/linearSOE/umfGEN; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn/linearSOE/cg; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn/linearSOE/BJsolvers; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn/linearSOE/diagonal; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn/linearSOE/petsc; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
