include ../../Makefile.def

OBJS       = Material.o

all:         $(OBJS)
	@$(CD) $(FE)/material/uniaxial; $(MAKE);
	@$(CD) $(FE)/material/nD; $(MAKE);
	@$(CD) $(FE)/material/section; $(MAKE);

# Miscellaneous

tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/material/uniaxial; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD; $(MAKE) wipe;
	@$(CD) $(FE)/material/section; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
