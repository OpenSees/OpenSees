include ../../../Makefile.def

OBJS       = FE_Element.o

# Compilation control

all:         $(OBJS)
	@$(CD) $(FE)/analysis/fe_ele/penalty; $(MAKE);
	@$(CD) $(FE)/analysis/fe_ele/lagrange; $(MAKE);
	@$(CD) $(FE)/analysis/fe_ele/transformation; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	@$(CD) $(FE)/analysis/fe_ele/penalty; $(MAKE) wipe;
	@$(CD) $(FE)/analysis/fe_ele/lagrange; $(MAKE) wipe;
	@$(CD) $(FE)/analysis/fe_ele/transformation; $(MAKE) wipe;


# DO NOT DELETE THIS LINE -- make depend depends on it.
