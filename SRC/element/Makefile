include ../../Makefile.def

OBJS       = Element.o ElementalLoad.o  Information.o TclElementCommands.o 
# Compilation control

all:         $(OBJS)
	@$(CD) $(FE)/element/beam2d; $(MAKE);
	@$(CD) $(FE)/element/beam3d; $(MAKE);
	@$(CD) $(FE)/element/truss; $(MAKE);
	@$(CD) $(FE)/element/zeroLength; $(MAKE);
	@$(CD) $(FE)/element/feap; $(MAKE);
	@$(CD) $(FE)/element/nonlinearBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/fourNodeQuad; $(MAKE);
	@$(CD) $(FE)/element/beamWithHinges; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) Makefile.bak $(OBJS) *.o *~ #*# core

spotless: clean
	@$(CD) $(FE)/element/beam2d; $(MAKE) wipe;
	@$(CD) $(FE)/element/beam3d; $(MAKE) wipe;
	@$(CD) $(FE)/element/truss; $(MAKE) wipe;
	@$(CD) $(FE)/element/zeroLength; $(MAKE) wipe;
	@$(CD) $(FE)/element/feap; $(MAKE) wipe;
	@$(CD) $(FE)/element/nonlinearBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/fourNodeQuad; $(MAKE) wipe;
	@$(CD) $(FE)/element/beamWithHinges; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
