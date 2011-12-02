include ../../Makefile.def

all: 
	@$(CD) $(FE)/../PACKAGES/NewMaterial/cpp; $(MAKE);
	@$(CD) $(FE)/../PACKAGES/NewMaterial/c; $(MAKE);
	@$(CD) $(FE)/../PACKAGES/NewMaterial/fortran; $(MAKE);

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(CD) $(FE)/../PACKAGES/NewMaterial/cpp; $(MAKE) wipe;
	@$(CD) $(FE)/../PACKAGES/NewMaterial/c; $(MAKE) wipe;
	@$(CD) $(FE)/../PACKAGES/NewMaterial/fortran; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
