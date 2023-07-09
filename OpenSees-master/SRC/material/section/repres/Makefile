include ../../../../Makefile.def

OBJS       = 

# Compilation control

all:	
	@$(CD) $(FE)/material/section/repres/cell; $(MAKE);
	@$(CD) $(FE)/material/section/repres/patch; $(MAKE);
	@$(CD) $(FE)/material/section/repres/reinfBar; $(MAKE);
	@$(CD) $(FE)/material/section/repres/reinfLayer; $(MAKE);
	@$(CD) $(FE)/material/section/repres/section; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy

spotless: clean

wipe: spotless
	@$(CD) $(FE)/material/section/repres/cell; $(MAKE) wipe;
	@$(CD) $(FE)/material/section/repres/patch; $(MAKE) wipe;
	@$(CD) $(FE)/material/section/repres/reinfBar; $(MAKE) wipe;
	@$(CD) $(FE)/material/section/repres/reinfLayer; $(MAKE) wipe;
	@$(CD) $(FE)/material/section/repres/section; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
