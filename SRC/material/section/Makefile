include ../../../Makefile.def

OBJS       = SectionForceDeformation.o \
	ElasticSection2d.o \
	ElasticSection3d.o \
	GenericSection1d.o \
	SectionAggregator.o \
	ElasticPlateSection.o \
	ElasticMembranePlateSection.o \
	MembranePlateFiberSection.o \
	TclModelBuilderSectionCommand.o \
	FiberSection2d.o \
	FiberSection3d.o \
	Bidirectional.o \
	FiberSectionGJ.o \
	Isolator2spring.o

all:         $(OBJS)
	@$(CD) $(FE)/material/section/repres; $(MAKE);
	@$(CD) $(FE)/material/section/fiber; $(MAKE);
	@$(CD) $(FE)/material/section/yieldSurface; $(MAKE);

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	@$(CD) $(FE)/material/section/repres; $(MAKE) wipe;
	@$(CD) $(FE)/material/section/fiber; $(MAKE) wipe;
	@$(CD) $(FE)/material/section/yieldSurface; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depends on it.
