include ../../../Makefile.def

OBJS       = SectionForceDeformation.o \
	ElasticSection2d.o ElasticSection3d.o \
	FiberSection.o \
	GenericSection1d.o GenericSectionNd.o \
	SectionAggregator.o \
	TclModelBuilderSectionCommand.o

all:         $(OBJS)
	@$(CD) $(FE)/material/section/repres; $(MAKE);
	@$(CD) $(FE)/material/section/fiber; $(MAKE);

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	@$(CD) $(FE)/material/section/repres; $(MAKE) wipe;
	@$(CD) $(FE)/material/section/fiber; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depends on it.
