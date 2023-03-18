# Makefile for model objects

include ../../../Makefile.def

OBJS       = TclModelBuilder.o myCommands.o TclUniaxialMaterialTester.o \
	Block2D.o Block3D.o TclSectionTester.o TclPlaneStressMaterialTester.o

# Compilation control

all: $(OBJS) 

tcl: $(OBJS) 
	@$(CD) $(FE)/tcl; $(MAKE) tcl;
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_TCL) \
	$(TCL_COMMAND_LIBS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM);
	@$(ECHO) @$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_TCL) \
	$(TCL_COMMAND_LIBS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM);



tk: $(OBJS)  
	@$(CD) $(FE)/tcl; $(MAKE) tk;
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_TK) \
	$(TCL_COMMAND_LIBS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)Tk

ifeq ($(PROGRAMMING_MODE), PARALLEL)

tcl: $(OBJS)  
	@$(CD) $(FE)/tcl; $(MAKE) tcl;
	$(ECHO)	$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_MPI) \
	$(TCL_COMMAND_LIBS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_MPI) \
	$(TCL_COMMAND_LIBS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)

tk: $(OBJS)  
	@$(CD) $(FE)/tcl; $(MAKE) tk;
	$(ECHO)	$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_MPI) \
	$(INTERPRETER_LIBS_TK) \
	$(TCL_COMMAND_LIBS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)Tk
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_MPI) \
	$(TCL_COMMAND_LIBS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)Tk

endif

ifeq ($(PROGRAMMING_MODE), PARALLEL_INTERPRETERS)

tcl: $(OBJS)  
	@$(CD) $(FE)/tcl; $(MAKE) tcl;
	@$(CD) $(FE)/system_of_eqn/diagonal; $(MAKE) mpi;
	@$(CD) $(FE)/actor/channel; $(MAKE) mpi;
	$(ECHO)	$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_MPI_PARAMETERS) \
	$(TCL_COMMAND_LIBS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_MPI_PARAMETERS) \
	$(TCL_COMMAND_LIBS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)

tk: $(OBJS)  
	@$(CD) $(FE)/tcl; $(MAKE) tk;
	@$(CD) $(FE)/system_of_eqn/diagonal; $(MAKE) mpi;
	@$(CD) $(FE)/actor/channel; $(MAKE) mpi;
	$(ECHO)	$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_MPI_PARAMETERS) \
	$(INTERPRETER_LIBS_TK) \
	$(TCL_COMMAND_LIBS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)Tk
	@$(LINKER) $(LINKFLAGS) $(INTERPRETER_LIBS_MPI_PARAMETERS) \
	$(TCL_COMMAND_LIBS) $(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY)  \
	$(MACHINE_SPECIFIC_LIBS) -o $(OpenSees_PROGRAM)Tk

endif

GRAPHIC_LIBS_TCL = $(FE)/tcl/TclFeViewer.o \
	$(FE)/tcl/TclVideoPlayer.o 

ifeq ($(GRAPHICS), NONE)

GRAPHIC_LIBS_TCL = $(FE)/tcl/TclFeViewer.o \
	$(FE)/renderer/PlainMap.o

endif


INTERPRETER_LIBS_TCL = $(FE)/tcl/tclMain.o \
	$(FE)/tcl/tclAppInit.o \
	$(FE)/tcl/commands.o \
	$(FE)/modelbuilder/tcl/myCommands.o $(GRAPHIC_LIBS_TCL)

TCL_COMMAND_LIBS = $(FE)/domain/pattern/TclPatternCommand.o \
	$(FE)/domain/pattern/TclSeriesCommand.o \
	$(FE)/domain/pattern/TclSeriesIntegratorCommand.o \
	$(FE)/domain/groundMotion/TclGroundMotionCommand.o \
	$(FE)/material/uniaxial/limitState/TclLimitState.o \
	$(FE)/element/zeroLength/TclZeroLength.o \
	$(FE)/element/forceBeamColumn/TclForceBeamColumnCommand.o \
	$(FE)/coordTransformation/TclGeomTransfCommand.o \
	$(FE)/damping/TclDampingCommand.o \
	$(FE)/element/joint/TclJoint2dCommand.o \
	$(FE)/element/joint/TclJoint3dCommand.o \
	$(FE)/element/joint/TclBeamColumnJointCommand.o \
	$(FE)/element/updatedLagrangianBeamColumn/TclElement2dYS.o \
	$(FE)/element/updatedLagrangianBeamColumn/TclElement2dGNL.o \
	$(FE)/element/updatedLagrangianBeamColumn/TclCyclicModelCommands.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/TclModelBuilderYieldSurfaceBCCommand.o \
	$(FE)/material/yieldSurface/evolution/TclModelBuilderYS_EvolutionCommand.o \
	$(FE)/material/yieldSurface/plasticHardeningMaterial/TclModelBuilderYSPlasticMaterialCommand.o \
	$(FE)/damage/TclModelBuilderDamageModelCommand.o \
	$(FE)/material/uniaxial/PY/TclPyTzQzMaterialCommand.o \
	$(FE)/material/nD/soil/TclUpdateMaterialStageCommand.o \
	$(FE)/element/UP-ucsd/TclFourNodeQuadUPCommand.o \
	$(FE)/recorder/TclRecorderCommands.o \
	$(FE)/database/TclDatabaseCommands.o \
	$(FE)/element/TclElementCommands.o \
	$(FE)/element/dispBeamColumnInt/TclDispBeamColumnIntCommand.o \
	$(FE)/element/elasticBeamColumn/TclElasticBeamCommand.o \
	$(FE)/element/fourNodeQuad/TclFourNodeQuadCommand.o \
	$(FE)/element/brick/TclTwenty_Node_BrickCommand.o \
	$(FE)/element/brick/TclBrickCommand.o \
	$(FE)/element/feap/TclFeapElementCommand.o \
	$(FE)/material/nD/feap/TclFeapMaterialCommand.o \
	$(FE)/element/frictionBearing/frictionModel/TclModelBuilderFrictionModelCommand.o \
	$(FE)/element/PFEMElement/TclModelBuilder_addPFEMElement.o \
	$(FE)/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.o \
	$(FE)/material/uniaxial/backbone/TclModelBuilderBackboneCommand.o \
	$(FE)/material/uniaxial/stiffness/TclModelBuilderStiffnessDegradationCommand.o \
	$(FE)/material/uniaxial/strength/TclModelBuilderStrengthDegradationCommand.o \
	$(FE)/material/uniaxial/unloading/TclModelBuilderUnloadingRuleCommand.o \
	$(FE)/material/nD/TclModelBuilderNDMaterialCommand.o \
	$(FE)/material/section/yieldSurface/TclModelBuilderYS_SectionCommand.o \
	$(FE)/material/section/TclModelBuilderSectionCommand.o \
	$(FE)/modelbuilder/tcl/TclUniaxialMaterialTester.o \
	$(FE)/modelbuilder/tcl/TclPlaneStressMaterialTester.o \
	$(FE)/modelbuilder/tcl/TclSectionTester.o \
	$(FE)/modelbuilder/tcl/TclModelBuilder.o \
	$(FE)/domain/component/TclUpdateMaterialCommand.o \
	$(FE)/domain/component/TclParameterCommands.o \
	$(FE)/domain/region/TclRegionCommands.o \
	$(FE)/material/uniaxial/snap/TclSnapMaterialCommand.o \
	$(FE)/material/uniaxial/fedeas/TclFedeasMaterialCommand.o \
	$(FE)/material/uniaxial/drain/TclDrainMaterialCommand.o \
	$(FE)/api/elementAPI_TCL.o


INTERPRETER_LIBS_MPI = $(FE)/tcl/tclMain.o \
	$(FE)/tcl/mpiMain.o \
	$(FE)/tcl/commands.o \
	$(FE)/modelbuilder/tcl/myCommands.o \
	$(FE)/actor/channel/MPI_Channel.o \
	$(FE)/actor/machineBroker/MPI_MachineBroker.o \
	$(FE)/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.o \
	$(FE)/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.o \
	$(FE)/domain/subdomain/ActorSubdomain.o \
	$(FE)/domain/subdomain/ShadowSubdomain.o $(GRAPHIC_LIBS_TCL)

INTERPRETER_LIBS_MPI_PARAMETERS = $(FE)/tcl/tclMain.o \
	$(FE)/tcl/mpiParameterMain.o \
	$(FE)/tcl/commands.o \
	$(FE)/modelbuilder/tcl/myCommands.o \
	$(FE)/actor/channel/MPI_Channel.o \
	$(FE)/actor/machineBroker/MPI_MachineBroker.o \
	$(FE)/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.o \
	$(FE)/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.o \
	$(FE)/domain/subdomain/ActorSubdomain.o \
	$(FE)/domain/subdomain/ShadowSubdomain.o $(GRAPHIC_LIBS_TCL)

INTERPRETER_LIBS_TK = $(FE)/tcl/tkMain.o \
	$(FE)/tcl/tkAppInit.o \
	$(FE)/tcl/commands.o \
	$(FE)/modelbuilder/tcl/myCommands.o $(GRAPHIC_LIBS_TCL)

GRAPHIC_LIBS_TCL = $(FE)/tcl/TclFeViewer.o \
	$(FE)/tcl/TclVideoPlayer.o 


# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o g3

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
