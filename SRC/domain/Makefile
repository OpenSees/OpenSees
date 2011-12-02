include ../../Makefile.def

OBJS       = 

#	@$(CD) $(FE)/domain/partitioner; $(MAKE);
#	@$(CD) $(FE)/domain/loadBalancer; $(MAKE);

all:     
	@$(CD) $(FE)/domain/domain; $(MAKE);
	@$(CD) $(FE)/domain/component; $(MAKE);
	@$(CD) $(FE)/domain/load; $(MAKE);
	@$(CD) $(FE)/domain/node; $(MAKE);
	@$(CD) $(FE)/domain/subdomain; $(MAKE);
	@$(CD) $(FE)/domain/constraints; $(MAKE);
	@$(CD) $(FE)/domain/pattern; $(MAKE);
	@$(CD) $(FE)/domain/groundMotion; $(MAKE);
	@$(CD) $(FE)/domain/partitioner; $(MAKE);
	@$(CD) $(FE)/domain/region; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy

spotless: clean

wipe: spotless
	@$(CD) $(FE)/domain/domain; $(MAKE) wipe;
	@$(CD) $(FE)/domain/component; $(MAKE) wipe;
	@$(CD) $(FE)/domain/load; $(MAKE) wipe;
	@$(CD) $(FE)/domain/subdomain; $(MAKE) wipe;
	@$(CD) $(FE)/domain/node; $(MAKE) wipe;
	@$(CD) $(FE)/domain/constraints; $(MAKE) wipe;
	@$(CD) $(FE)/domain/partitioner; $(MAKE) wipe;
	@$(CD) $(FE)/domain/loadBalancer; $(MAKE) wipe;
	@$(CD) $(FE)/domain/pattern; $(MAKE) wipe;
	@$(CD) $(FE)/domain/groundMotion; $(MAKE) wipe;
	@$(CD) $(FE)/domain/partitioner; $(MAKE) wipe;
	@$(CD) $(FE)/domain/region; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
