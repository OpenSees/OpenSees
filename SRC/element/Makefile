include ../../Makefile.def

OBJS       = Element.o ElementalLoad.o  Information.o TclElementCommands.o NewElement.o WrapperElement.o

# Compilation control
#	@$(CD) $(FE)/element/8nbrick; $(MAKE);
#	@$(CD) $(FE)/element/20nbrick; $(MAKE);
#	@$(CD) $(FE)/element/27nbrick; $(MAKE);
#	@$(CD) $(FE)/element/UP_ucdavis; $(MAKE);
#	@$(CD) $(FE)/element/upU; $(MAKE);
#	@$(CD) $(FE)/element/TotalLagrangianFD20NodeBrick; $(MAKE);

all:         $(OBJS)
	@$(CD) $(FE)/element/beam2d; $(MAKE);
	@$(CD) $(FE)/element/beam3d; $(MAKE);
	@$(CD) $(FE)/element/truss; $(MAKE);
	@$(CD) $(FE)/element/zeroLength; $(MAKE);
	@$(CD) $(FE)/element/feap; $(MAKE);
	@$(CD) $(FE)/element/nonlinearBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/fourNodeQuad; $(MAKE);
	@$(CD) $(FE)/element/beamWithHinges; $(MAKE);
	@$(CD) $(FE)/element/UP-ucsd; $(MAKE);
	@$(CD) $(FE)/element/shell; $(MAKE);
	@$(CD) $(FE)/element/brick; $(MAKE);
	@$(CD) $(FE)/element/joint; $(MAKE);
	@$(CD) $(FE)/element/triangle; $(MAKE);
	@$(CD) $(FE)/element/PFEMElement; $(MAKE);
	@$(CD) $(FE)/element/PML; $(MAKE);
	@$(CD) $(FE)/element/dispBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/mixedBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/elasticBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/updatedLagrangianBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/forceBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/dispBeamColumnInt; $(MAKE);
	@$(CD) $(FE)/element/generic; $(MAKE);
	@$(CD) $(FE)/element/elastomericBearing; $(MAKE);
	@$(CD) $(FE)/element/frictionBearing; $(MAKE);
	@$(CD) $(FE)/element/adapter; $(MAKE);
	@$(CD) $(FE)/element/twoNodeLink; $(MAKE);
	@$(CD) $(FE)/element/pyMacro; $(MAKE);
	@$(CD) $(FE)/element/surfaceLoad; $(MAKE);
	@$(CD) $(FE)/element/UWelements; $(MAKE);
	@$(CD) $(FE)/element/XMUelements; $(MAKE);
	@$(CD) $(FE)/element/HUelements; $(MAKE);
	@$(CD) $(FE)/element/IGA; $(MAKE);
	@$(CD) $(FE)/element/XMUelements; $(MAKE);
	@$(CD) $(FE)/element/mvlem; $(MAKE);
	@$(CD) $(FE)/element/PFEMElement; $(MAKE);
	@$(CD) $(FE)/element/componentElement; $(MAKE);
	@$(CD) $(FE)/element/catenaryCable; $(MAKE);
	@$(CD) $(FE)/element/tetrahedron; $(MAKE);
	@$(CD) $(FE)/element/absorbentBoundaries; $(MAKE);
	@$(CD) $(FE)/element/gradientInelasticBeamColumn; $(MAKE);
	@$(CD) $(FE)/element/RockingBC; $(MAKE);
	@$(CD) $(FE)/element/masonry; $(MAKE);
	@$(CD) $(FE)/element/CEqElement; $(MAKE);
	@$(CD) $(FE)/element/mefi; $(MAKE);

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
	@$(CD) $(FE)/element/UP-ucsd; $(MAKE) wipe;
	@$(CD) $(FE)/element/shell; $(MAKE) wipe;
	@$(CD) $(FE)/element/brick; $(MAKE) wipe;
	@$(CD) $(FE)/element/joint; $(MAKE) wipe;
	@$(CD) $(FE)/element/triangle; $(MAKE) wipe;
	@$(CD) $(FE)/element/PFEMElement; $(MAKE) wipe;
	@$(CD) $(FE)/element/PML; $(MAKE) wipe;
	@$(CD) $(FE)/element/dispBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/mixedBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/elasticBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/updatedLagrangianBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/forceBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/dispBeamColumnInt; $(MAKE) wipe;
	@$(CD) $(FE)/element/generic; $(MAKE) wipe;
	@$(CD) $(FE)/element/elastomericBearing; $(MAKE) wipe;
	@$(CD) $(FE)/element/adapter; $(MAKE) wipe;
	@$(CD) $(FE)/element/twoNodeLink; $(MAKE) wipe;
	@$(CD) $(FE)/element/frictionBearing; $(MAKE) wipe;
	@$(CD) $(FE)/element/pyMacro; $(MAKE) wipe;
	@$(CD) $(FE)/element/surfaceLoad; $(MAKE) wipe;
	@$(CD) $(FE)/element/UWelements; $(MAKE) wipe;
	@$(CD) $(FE)/element/XMUelements; $(MAKE) wipe;
	@$(CD) $(FE)/element/HUelements; $(MAKE) wipe;
	@$(CD) $(FE)/element/IGA; $(MAKE) wipe;
	@$(CD) $(FE)/element/XMUelements; $(MAKE) wipe;
	@$(CD) $(FE)/element/mvlem; $(MAKE) wipe;
	@$(CD) $(FE)/element/PFEMElement; $(MAKE) wipe;
	@$(CD) $(FE)/element/componentElement; $(MAKE) wipe;
	@$(CD) $(FE)/element/catenaryCable; $(MAKE) wipe;
	@$(CD) $(FE)/element/tetrahedron; $(MAKE) wipe;
	@$(CD) $(FE)/element/absorbentBoundaries; $(MAKE) wipe;
	@$(CD) $(FE)/element/gradientInelasticBeamColumn; $(MAKE) wipe;
	@$(CD) $(FE)/element/RockingBC; $(MAKE) wipe;
	@$(CD) $(FE)/element/masonry; $(MAKE) wipe;
	@$(CD) $(FE)/element/CEqElement; $(MAKE) wipe;
	@$(CD) $(FE)/element/mefi; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
