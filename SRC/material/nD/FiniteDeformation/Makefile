include ../../../../Makefile.def

OBJS       = FiniteDeformationElastic3D.o \
             NeoHookeanCompressible3D.o \
             FDdecoupledElastic3D.o \
             TclFiniteDeformationElastic3DCommand.o \
             W.o \
             LogWEnergy.o \
             MooneyRivlinWEnergy.o \
             NeoHookeanWEnergy.o \
             OgdenWEnergy.o \
             SimoPisterWEnergy.o \
             OgdenSimoWEnergy.o \
             MooneyRivlinSimoWEnergy.o \
             FDEPState.o \
             FiniteDeformationEP3D.o \
             TclFiniteDeformationEP3DCommand.o

OBJS1      = FiniteDeformationElastic3D.o \
             NeoHookeanCompressible3D.o \
             FDdecoupledElastic3D.o \
             TclFiniteDeformationElastic3DCommand.o \
             W.o \
             LogWEnergy.o \
             MooneyRivlinWEnergy.o \
             NeoHookeanWEnergy.o \
             OgdenWEnergy.o \
             SimoPisterWEnergy.o \
             OgdenSimoWEnergy.o \
             MooneyRivlinSimoWEnergy.o \
             FDEPState.o \
             FiniteDeformationEP3D.o \
             TclFiniteDeformationEP3DCommand.o
 

all: $(OBJS)
	@$(CD) $(FE)/material/nD/FiniteDeformation/fdYield; $(MAKE);
	@$(CD) $(FE)/material/nD/FiniteDeformation/fdFlow; $(MAKE);
	@$(CD) $(FE)/material/nD/FiniteDeformation/fdEvolution; $(MAKE);

fd_test: $(OBJS)   fd_test.o
	$(LINKER) $(LINKFLAGS) fd_test.o $(OBJS1)  \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(TCL_LIBRARY)\
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o fd_test

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS)  *.o *.bak
	@$(CD) $(FE)/material/nD/FiniteDeformation/fdYield; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/FiniteDeformation/fdFlow; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/FiniteDeformation/fdEvolution; $(MAKE) wipe;
	
spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.




