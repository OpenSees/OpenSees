include ../../Makefile.def

all: 
	@$(CD) $(FE)/../PACKAGES/NewElement/cpp; $(MAKE);
	@$(CD) $(FE)/../PACKAGES/NewElement/c; $(MAKE);
	@$(CD) $(FE)/../PACKAGES/NewElement/fortran; $(MAKE);

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(CD) $(FE)/../PACKAGES/NewElement/cpp; $(MAKE) wipe;
	@$(CD) $(FE)/../PACKAGES/NewElement/c; $(MAKE) wipe;
	@$(CD) $(FE)/../PACKAGES/NewElement/fortran; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
