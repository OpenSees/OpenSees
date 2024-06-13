include ../../../Makefile.def

OBJS       = SectionForceDeformation.o \
	FiberSection.o \
	ElasticSection2d.o ElasticSection3d.o \
	ElasticShearSection2d.o ElasticShearSection3d.o \
	ElasticBDShearSection2d.o \
	ElasticWarpingShearSection2d.o \
	ElasticTubeSection3d.o \
	GenericSection1d.o \
	SectionAggregator.o \
	ParallelSection.o \
	ElasticPlateSection.o \
	ElasticMembranePlateSection.o \
	MembranePlateFiberSection.o \
	DoubleMembranePlateFiberSection.o \
	TclModelBuilderSectionCommand.o \
	FiberSection2d.o \
	NDFiberSection2d.o \
	NDFiberSection3d.o \
	NDFiberSectionWarping2d.o \
	FiberSection2dThermal.o \
	FiberSection3d.o \
	FiberSectionWarping3d.o \
	FiberSectionAsym3d.o \
	Bidirectional.o \
	Elliptical.o \
	Elliptical2.o \
	Isolator2spring.o \
	BiaxialHysteretic.o \
	McftSection2dfiber.o \
	LayeredShellFiberSection.o \
	LayeredShellFiberSectionThermal.o \
	FiberSection3dThermal.o \
	MembranePlateFiberSectionThermal.o \
	FiberSectionGJThermal.o

all:         $(OBJS)
	@$(CD) $(FE)/material/section/repres; $(MAKE);
	@$(CD) $(FE)/material/section/fiber; $(MAKE);
	@$(CD) $(FE)/material/section/yieldSurface; $(MAKE);
	@$(CD) $(FE)/material/section/integration; $(MAKE);
	@$(CD) $(FE)/material/section/LayeredMembraneSection; $(MAKE);

test: main.o
	@$(LINKER) $(LINKFLAGS) main.o $(INTERPRETER_LIBS_MPI) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o section_tst

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
	@$(CD) $(FE)/material/section/integration; $(MAKE) wipe;
	@$(CD) $(FE)/material/section/LayeredMembraneSection; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depends on it.
