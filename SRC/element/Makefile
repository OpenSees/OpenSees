include ../../Makefile.def

OBJS       = Element.o ElementalLoad.o  Information.o TclElementCommands.o NewElement.o
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
	@$(CD) $(FE)/element/8nbrick; $(MAKE);
	@$(CD) $(FE)/element/20nbrick; $(MAKE);
	@$(CD) $(FE)/element/27nbrick; $(MAKE);
	@$(CD) $(FE)/element/TotalLagrangianFD20NodeBrick; $(MAKE);
	@$(CD) $(FE)/element/upU; $(MAKE);
	@$(CD) $(FE)/element/shell; $(MAKE);
	@$(CD) $(FE)/element/brick; $(MAKE);
	@$(CD) $(FE)/element/joint; $(MAKE);
	@$(CD) $(FE)/element/dispBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/elasticBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/updatedLagrangianBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/forceBeamColumn; $(MAKE);

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
	@$(CD) $(FE)/element/8nbrick; $(MAKE) wipe;
	@$(CD) $(FE)/element/20nbrick; $(MAKE) wipe;
	@$(CD) $(FE)/element/27nbrick; $(MAKE) wipe;
	@$(CD) $(FE)/element/TotalLagrangianFD20NodeBrick; $(MAKE) wipe;
	@$(CD) $(FE)/element/upU; $(MAKE) wipe;
	@$(CD) $(FE)/element/shell; $(MAKE) wipe;
	@$(CD) $(FE)/element/brick; $(MAKE) wipe;
	@$(CD) $(FE)/element/joint; $(MAKE) wipe;
	@$(CD) $(FE)/element/dispBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/elasticBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/updatedLagrangianBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/forceBeamColumn; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
