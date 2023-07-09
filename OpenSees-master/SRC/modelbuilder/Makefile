include ../../Makefile.def

#	PartitionedModelBuilder.o PartitionedQuick2dFrame.o

OBJS       = ModelBuilder.o PlaneFrame.o

# Compilation control

all:    $(OBJS)
	@$(CD) $(FE)/modelbuilder/tcl; $(MAKE) all;

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	@$(CD) $(FE)/modelbuilder/tcl; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
