include ../../Makefile.def

OBJS       = TaggedObject.o

# Compilation control

all:         $(OBJS)
	@$(CD) $(FE)/tagged/storage; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/tagged/storage; $(MAKE) clean;

wipe: spotless
	@$(CD) $(FE)/tagged/storage; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.


