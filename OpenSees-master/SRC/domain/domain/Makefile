include ../../../Makefile.def

OBJS       = Domain.o DomainModalProperties.o

# Compilation control

all:         $(OBJS)
	     @$(CD) $(FE)/domain/domain/single; $(MAKE);
	     @$(CD) $(FE)/domain/domain/partitioned; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(CD) $(FE)/domain/domain/single; $(MAKE) wipe clean;
	@$(CD) $(FE)/domain/domain/partitioned; $(MAKE) wipe clean;
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	@$(CD) $(FE)/domain/domain/single; $(MAKE) wipe;
	@$(CD) $(FE)/domain/domain/partitioned; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
