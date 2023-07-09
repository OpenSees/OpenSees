include ../../../Makefile.def

OBJS       = 

# Compilation control
all:         
	     @$(CD) $(FE)/reliability/domain/components; $(MAKE);
	     @$(CD) $(FE)/reliability/domain/distributions; $(MAKE);
	     @$(CD) $(FE)/reliability/domain/filter; $(MAKE);
	     @$(CD) $(FE)/reliability/domain/functionEvaluator; $(MAKE);
	     @$(CD) $(FE)/reliability/domain/performanceFunction; $(MAKE);
	     @$(CD) $(FE)/reliability/domain/modulatingFunction; $(MAKE);
	     @$(CD) $(FE)/reliability/domain/spectrum; $(MAKE);
	     @$(CD) $(FE)/reliability/domain/storage; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	     @$(CD) $(FE)/reliability/domain/components; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/domain/distributions; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/domain/filter; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/domain/functionEvaluator; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/domain/performanceFunction; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/domain/modulatingFunction; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/domain/spectrum; $(MAKE) wipe;
	     @$(CD) $(FE)/reliability/domain/storage; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
