include ../../Makefile.def

all:        
	     @$(CD) $(FE)/analysis/analysis; $(MAKE);
	     @$(CD) $(FE)/analysis/integrator; $(MAKE);
	     @$(CD) $(FE)/analysis/model; $(MAKE);
	     @$(CD) $(FE)/analysis/numberer; $(MAKE);
	     @$(CD) $(FE)/analysis/handler; $(MAKE);
	     @$(CD) $(FE)/analysis/algorithm; $(MAKE);
	     @$(CD) $(FE)/analysis/dof_grp; $(MAKE);
	     @$(CD) $(FE)/analysis/fe_ele; $(MAKE);

# Miscellaneous
tidy:	
	@rm -f Makefile.bak *~ #*# core

clean: tidy
	@rm -f $(OBJS) *.o

spotless: clean

wipe: spotless
	     @$(CD) $(FE)/analysis/analysis; $(MAKE) wipe;
	     @$(CD) $(FE)/analysis/integrator; $(MAKE) wipe;
	     @$(CD) $(FE)/analysis/model; $(MAKE) wipe;
	     @$(CD) $(FE)/analysis/numberer; $(MAKE) wipe;
	     @$(CD) $(FE)/analysis/handler; $(MAKE) wipe;
	     @$(CD) $(FE)/analysis/algorithm; $(MAKE) wipe;
	     @$(CD) $(FE)/analysis/dof_grp; $(MAKE) wipe;
	     @$(CD) $(FE)/analysis/fe_ele; $(MAKE) wipe;


# DO NOT DELETE THIS LINE -- make depend depends on it.
