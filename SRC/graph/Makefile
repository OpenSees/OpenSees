include ../../Makefile.def

OBJS       = 

# Compilation control

all:         $(OBJS)
	@$(CD) graph; $(MAKE);
	@$(CD) numberer; $(MAKE);
	@$(CD) partitioner; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core 

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	@$(CD) graph; $(MAKE) wipe;
	@$(CD) numberer; $(MAKE) wipe;
	@$(CD) partitioner; $(MAKE) wipe;


# DO NOT DELETE THIS LINE -- make depend depends on it.


