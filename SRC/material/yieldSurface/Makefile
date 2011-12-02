include ../../../Makefile.def


all:
	@$(CD) $(FE)/material/yieldSurface/yieldSurfaceBC; $(MAKE);
	@$(CD) $(FE)/material/yieldSurface/plasticHardeningMaterial; $(MAKE);
	@$(CD) $(FE)/material/yieldSurface/evolution; $(MAKE);

# Miscellaneous

tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/material/yieldSurface/yieldSurfaceBC; $(MAKE) wipe;
	@$(CD) $(FE)/material/yieldSurface/plasticHardeningMaterial; $(MAKE) wipe;
	@$(CD) $(FE)/material/yieldSurface/evolution; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
