# Makefile for fe objects
ifndef SHARED

include ../Makefile.def

endif


ifeq ($(RELIABILITY), YES_RELIABILITY)

all:
	@$(RM) $(RMFLAGS) $(FE_LIBRARY);
	@$(CD) $(FE)/matrix; $(MAKE);
	@$(CD) $(FE)/database; $(MAKE);
	@$(CD) $(FE)/system_of_eqn; $(MAKE);
	@$(CD) $(FE)/modelbuilder; $(MAKE);
	@$(CD) $(FE)/tagged; $(MAKE);
	@$(CD) $(FE)/material; $(MAKE);
	@$(CD) $(FE)/analysis; $(MAKE);
	@$(CD) $(FE)/graph; $(MAKE);
	@$(CD) $(FE)/utility; $(MAKE);
	@$(CD) $(FE)/convergenceTest; $(MAKE)
	@$(CD) $(FE)/recorder; $(MAKE);
	@$(CD) $(FE)/tcl; $(MAKE);
	@$(CD) $(FE)/renderer; $(MAKE);
	@$(CD) $(FE)/domain; $(MAKE);
	@$(CD) $(FE)/handler; $(MAKE);
	@$(CD) $(FE)/actor; $(MAKE);
	@$(CD) $(FE)/machine; $(MAKE);
	@$(CD) $(FE)/element; $(MAKE);
	@$(CD) $(FE)/coordTransformation; $(MAKE);
	@$(CD) $(FE)/api; $(MAKE);
	@$(CD) $(FE)/damage; $(MAKE);
	@$(CD) $(FE)/api; $(MAKE);
	@$(CD) $(FE)/reliability; $(MAKE);
	@$(CD) $(FE)/optimization; $(MAKE);
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(MATRIX_LIBS) 
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(ELE_LIBS) 
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(MATERIAL_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(SECTION_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(DATABASE_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(FEDEAS_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(MODEL_BUILDER_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(DOMAIN_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(PY_SJB_RWB_BJ_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(ANALYSIS_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(SysOfEqn_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(GRAPH_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(TAGGED_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(UTILITY_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(ACTOR_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(RECORDER_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(LAW_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(GGD_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(RENDERER_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(REMO_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(GLF_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(ZHY_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(ERRORHANDLER_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(DRAIN_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(SNAP_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(FMK_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(KJE_LIBS)
	@$(RANLIB) $(RANLIBFLAGS) $(FE_LIBRARY)
	@$(AR) $(ARFLAGS) $(HOME)/lib/libOpenSeesTclCommands.a $(TCL_COMMAND_LIBS)
	@$(RANLIB) $(RANLIBFLAGS) $(HOME)/lib/libOpenSeesTclCommands.a 

else

all:
	@$(RM) $(RMFLAGS) $(FE_LIBRARY);
	@$(CD) $(FE)/matrix; $(MAKE);
	@$(CD) $(FE)/database; $(MAKE);
	@$(CD) $(FE)/system_of_eqn; $(MAKE);
	@$(CD) $(FE)/modelbuilder; $(MAKE);
	@$(CD) $(FE)/tagged; $(MAKE);
	@$(CD) $(FE)/material; $(MAKE);
	@$(CD) $(FE)/analysis; $(MAKE);
	@$(CD) $(FE)/graph; $(MAKE);
	@$(CD) $(FE)/utility; $(MAKE);
	@$(CD) $(FE)/convergenceTest; $(MAKE);
	@$(CD) $(FE)/recorder; $(MAKE);
	@$(CD) $(FE)/tcl; $(MAKE);
	@$(CD) $(FE)/renderer; $(MAKE);
	@$(CD) $(FE)/domain; $(MAKE);
	@$(CD) $(FE)/handler; $(MAKE);
	@$(CD) $(FE)/actor; $(MAKE);
	@$(CD) $(FE)/machine; $(MAKE);
	@$(CD) $(FE)/element; $(MAKE);
	@$(CD) $(FE)/coordTransformation; $(MAKE);
	@$(CD) $(FE)/api; $(MAKE);
	@$(CD) $(FE)/damage; $(MAKE);
	@$(CD) $(FE)/api; $(MAKE);
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(MATRIX_LIBS) 
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(ELE_LIBS) 
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(MATERIAL_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(SECTION_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(DATABASE_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(FEDEAS_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(MODEL_BUILDER_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(DOMAIN_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(PY_SJB_RWB_BJ_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(ANALYSIS_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(SysOfEqn_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(GRAPH_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(TAGGED_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(UTILITY_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(ACTOR_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(RECORDER_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(LAW_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(GGD_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(RENDERER_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(REMO_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(GLF_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(ZHY_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(ERRORHANDLER_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(DRAIN_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(SNAP_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(FMK_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(KJE_LIBS)
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(ACTOR_LIBS)
	@$(RANLIB) $(RANLIBFLAGS) $(FE_LIBRARY)
	@$(AR) $(ARFLAGS) $(HOME)/lib/libOpenSeesTclCommands.a $(TCL_COMMAND_LIBS)
	@$(RANLIB) $(RANLIBFLAGS) $(HOME)/lib/libOpenSeesTclCommands.a 

endif

TCL_COMMAND_LIBS = $(FE)/domain/pattern/TclPatternCommand.o \
	$(FE)/domain/pattern/TclSeriesCommand.o \
	$(FE)/domain/pattern/TclSeriesIntegratorCommand.o \
	$(FE)/domain/groundMotion/TclGroundMotionCommand.o \
	$(FE)/material/uniaxial/limitState/TclLimitState.o \
	$(FE)/element/zeroLength/TclZeroLength.o \
	$(FE)/element/forceBeamColumn/TclForceBeamColumnCommand.o \
	$(FE)/element/elasticBeamColumn/TclElasticBeamCommand.o \
	$(FE)/coordTransformation/TclGeomTransfCommand.o \
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
	$(FE)/element/fourNodeQuad/TclFourNodeQuadCommand.o \
	$(FE)/element/brick/TclTwenty_Node_BrickCommand.o \
	$(FE)/element/brick/TclBrickCommand.o \
	$(FE)/element/adapter/TclActuatorCommand.o \
	$(FE)/element/adapter/TclActuatorCorotCommand.o \
	$(FE)/element/adapter/TclAdapterCommand.o \
	$(FE)/element/generic/TclGenericClientCommand.o \
	$(FE)/element/generic/TclGenericCopyCommand.o \
	$(FE)/element/elastomericBearing/TclElastomericBearingBoucWenCommand.o \
	$(FE)/element/feap/TclFeapElementCommand.o \
	$(FE)/material/nD/feap/TclFeapMaterialCommand.o \
	$(FE)/element/beamWithHinges/TclBeamWithHingesBuilder.o \
	$(FE)/element/frictionBearing/frictionModel/TclModelBuilderFrictionModelCommand.o \
	$(FE)/element/PFEMElement/TclModelBuilder_addPFEMElement.o \
	$(FE)/material/uniaxial/backbone/TclModelBuilderBackboneCommand.o \
	$(FE)/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.o \
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
	$(FE)/element/twoNodeLink/TclTwoNodeLinkCommand.o \
	$(FE)/element/frictionBearing/TclRJWatsonEQSCommand.o \
	$(FE)/material/uniaxial/snap/TclSnapMaterialCommand.o \
	$(FE)/element/frictionBearing/TclFlatSliderCommand.o \
	$(FE)/element/elastomericBearing/TclElastomericBearingPlasticityCommand.o \
	$(FE)/element/elastomericBearing/TclElastomericBearingUFRPCommand.o \
	$(FE)/element/frictionBearing/TclSingleFPCommand.o \
	$(FE)/material/uniaxial/fedeas/TclFedeasMaterialCommand.o \
	$(FE)/material/uniaxial/drain/TclDrainMaterialCommand.o \
	$(FE)/api/elementAPI_TCL.o


skypack:
	@$(CD) $(HOME)/lib; rm -f $(HOME)/lib/libskypack.a;
	@$(CD) $(FE)/system_of_eqn/linearSOE/profileSPD/skypack; creator;

metis:
	@$(CD) $(HOME)/lib; rm -f $(METIS_LIBRARY);
	@$(CD) $(FE)/graph/partitioner/metis-2.0; $(MAKE) clean;
	@$(CD) $(FE)/graph/partitioner/metis-2.0; $(MAKE) lib;

FMK_LIBS = $(FE)/database/FE_Datastore.o \
	$(FE)/api/packages.o 


KJE_LIBS = $(FE)/material/uniaxial/limitState/LimitStateMaterial.o \
	$(FE)/material/uniaxial/limitState/limitCurve/LimitCurve.o \
	$(FE)/material/uniaxial/limitState/limitCurve/AxialCurve.o \
	$(FE)/material/uniaxial/limitState/limitCurve/WrapperLimitCurve.o \
	$(FE)/material/uniaxial/limitState/limitCurve/RotationShearCurve.o \
	$(FE)/material/uniaxial/limitState/limitCurve/ShearCurve.o \
	$(FE)/material/uniaxial/limitState/limitCurve/ThreePointCurve.o 

# the files that make up the FE_LIBRARY
GLF_LIBS = $(FE)/element/zeroLength/ZeroLength.o \
	$(FE)/element/zeroLength/ZeroLengthSection.o \
	$(FE)/element/zeroLength/CoupledZeroLength.o \
	$(FE)/element/zeroLength/ZeroLengthND.o \
	$(FE)/element/zeroLength/ZeroLengthVG_HG.o \
	$(FE)/element/zeroLength/ZeroLengthImpact3D.o

REMO_LIBS = $(FE)/element/nonlinearBeamColumn/matrixutil/MatrixUtil.o \
	$(FE)/element/forceBeamColumn/ForceBeamColumn2d.o \
	$(FE)/element/forceBeamColumn/ForceBeamColumn2dThermal.o \
	$(FE)/element/forceBeamColumn/ForceBeamColumnCBDI2d.o \
	$(FE)/element/forceBeamColumn/ForceBeamColumnWarping2d.o \
	$(FE)/element/forceBeamColumn/ForceBeamColumn3d.o \
	$(FE)/element/forceBeamColumn/ElasticForceBeamColumn2d.o \
	$(FE)/element/forceBeamColumn/ElasticForceBeamColumnWarping2d.o \
	$(FE)/element/forceBeamColumn/ElasticForceBeamColumn3d.o \
	$(FE)/element/forceBeamColumn/BeamIntegration.o \
	$(FE)/element/forceBeamColumn/LobattoBeamIntegration.o \
	$(FE)/element/forceBeamColumn/CompositeSimpsonBeamIntegration.o \
	$(FE)/element/forceBeamColumn/LegendreBeamIntegration.o \
	$(FE)/element/forceBeamColumn/RadauBeamIntegration.o \
	$(FE)/element/forceBeamColumn/NewtonCotesBeamIntegration.o \
	$(FE)/element/forceBeamColumn/HingeMidpointBeamIntegration.o \
	$(FE)/element/forceBeamColumn/HingeEndpointBeamIntegration.o \
	$(FE)/element/forceBeamColumn/HingeRadauBeamIntegration.o \
	$(FE)/element/forceBeamColumn/HingeRadauTwoBeamIntegration.o \
	$(FE)/element/forceBeamColumn/UserDefinedBeamIntegration.o \
	$(FE)/element/forceBeamColumn/UserDefinedHingeIntegration.o \
	$(FE)/element/forceBeamColumn/DistHingeIntegration.o \
	$(FE)/element/forceBeamColumn/RegularizedHingeIntegration.o \
	$(FE)/element/forceBeamColumn/MidDistanceBeamIntegration.o \
	$(FE)/element/forceBeamColumn/TrapezoidalBeamIntegration.o \
	$(FE)/element/forceBeamColumn/FixedLocationBeamIntegration.o \
	$(FE)/element/forceBeamColumn/LowOrderBeamIntegration.o \
	$(FE)/coordTransformation/CrdTransf.o \
	$(FE)/coordTransformation/LinearCrdTransf2d.o \
	$(FE)/coordTransformation/PDeltaCrdTransf2d.o \
	$(FE)/coordTransformation/CorotCrdTransf2d.o \
	$(FE)/coordTransformation/CorotCrdTransfWarping2d.o \
	$(FE)/coordTransformation/LinearCrdTransf3d.o \
	$(FE)/coordTransformation/PDeltaCrdTransf3d.o \
	$(FE)/coordTransformation/CorotCrdTransf3d.o

LAW_LIBS = $(FE)/system_of_eqn/linearSOE/sparseSYM/SymSparseLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/SymSparseLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/grcm.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/nest.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/nmat.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/symbolic.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/utility.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/newordr.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/nnsim.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/genmmd.o \
	$(FE)/system_of_eqn/eigenSOE/EigenSOE.o \
	$(FE)/system_of_eqn/eigenSOE/EigenSolver.o \
	$(FE)/system_of_eqn/eigenSOE/ArpackSOE.o \
	$(FE)/system_of_eqn/eigenSOE/ArpackSolver.o \
	$(FE)/system_of_eqn/eigenSOE/SymBandEigenSOE.o \
	$(FE)/system_of_eqn/eigenSOE/SymBandEigenSolver.o \
	$(FE)/analysis/analysis/EigenAnalysis.o \
	$(FE)/analysis/integrator/EigenIntegrator.o 


GGD_LIBS = $(FE)/element/joint/MP_Joint2D.o \
	$(FE)/element/joint/MP_Joint3D.o \
	$(FE)/element/joint/Joint2D.o \
	$(FE)/element/joint/Joint3D.o \
	$(FE)/element/joint/BeamColumnJoint2d.o \
	$(FE)/element/joint/BeamColumnJoint3d.o \
	$(FE)/element/joint/ElasticTubularJoint.o \
	$(FE)/element/updatedLagrangianBeamColumn/Elastic2DGNL.o \
	$(FE)/element/updatedLagrangianBeamColumn/Inelastic2DYS02.o \
	$(FE)/element/updatedLagrangianBeamColumn/InelasticYS2DGNL.o \
	$(FE)/element/updatedLagrangianBeamColumn/Inelastic2DYS01.o \
	$(FE)/element/updatedLagrangianBeamColumn/Inelastic2DYS03.o \
	$(FE)/element/updatedLagrangianBeamColumn/UpdatedLagrangianBeam2D.o \
	$(FE)/element/updatedLagrangianBeamColumn/BilinearCyclic.o \
	$(FE)/element/updatedLagrangianBeamColumn/CyclicModel.o \
	$(FE)/element/updatedLagrangianBeamColumn/LinearCyclic.o \
	$(FE)/element/updatedLagrangianBeamColumn/QuadraticCyclic.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/Attalla2D.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/ElTawil2D.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/ElTawil2DUnSym.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/Hajjar2D.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/Orbison2D.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/YieldSurface_BC2D.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/YieldSurface_BC.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/NullYS2D.o \
	$(FE)/material/yieldSurface/evolution/BkStressLimSurface2D.o \
	$(FE)/material/yieldSurface/evolution/CombinedIsoKin2D01.o \
	$(FE)/material/yieldSurface/evolution/Isotropic2D01.o \
	$(FE)/material/yieldSurface/evolution/Kinematic2D02.o \
	$(FE)/material/yieldSurface/evolution/PeakOriented2D01.o \
	$(FE)/material/yieldSurface/evolution/PlasticHardening2D.o \
	$(FE)/material/yieldSurface/evolution/YS_Evolution.o \
	$(FE)/material/yieldSurface/evolution/BoundingSurface2D.o \
	$(FE)/material/yieldSurface/evolution/CombinedIsoKin2D02.o \
	$(FE)/material/yieldSurface/evolution/Kinematic2D01.o \
	$(FE)/material/yieldSurface/evolution/NullEvolution.o \
	$(FE)/material/yieldSurface/evolution/PeakOriented2D02.o \
	$(FE)/material/yieldSurface/evolution/YS_Evolution2D.o \
	$(FE)/material/yieldSurface/plasticHardeningMaterial/ExponReducing.o \
	$(FE)/material/yieldSurface/plasticHardeningMaterial/MultiLinearKp.o \
	$(FE)/material/yieldSurface/plasticHardeningMaterial/NullPlasticMaterial.o \
	$(FE)/material/yieldSurface/plasticHardeningMaterial/PlasticHardeningMaterial.o \
	$(FE)/damage/DamageModel.o \
	$(FE)/damage/Kratzig.o \
	$(FE)/damage/ParkAng.o \
	$(FE)/damage/DamageResponse.o \
	$(FE)/damage/HystereticEnergy.o \
	$(FE)/damage/Mehanny.o \
	$(FE)/damage/NormalizedPeak.o

ERRORHANDLER_LIBS = $(FE)/handler/StandardStream.o \
	$(FE)/handler/FileStream.o \
	$(FE)/handler/OPS_Stream.o \
	$(FE)/handler/DataFileStream.o \
	$(FE)/handler/DataFileStreamAdd.o \
	$(FE)/handler/XmlFileStream.o \
	$(FE)/handler/BinaryFileStream.o \
	$(FE)/handler/DummyStream.o \
	$(FE)/handler/TCP_Stream.o \
	$(FE)/handler/DatabaseStream.o 


PY_SJB_RWB_BJ_LIBS = $(FE)/material/uniaxial/PY/PySimple1.o \
	$(FE)/material/uniaxial/PY/QzSimple1.o \
	$(FE)/material/uniaxial/PY/TzSimple1.o \
	$(FE)/material/uniaxial/PY/PySimple2.o \
	$(FE)/material/uniaxial/PY/PySimple3.o \
	$(FE)/material/uniaxial/PY/QzSimple2.o \
	$(FE)/material/uniaxial/PY/TzSimple2.o \
	$(FE)/material/uniaxial/PY/ShallowFoundationGen.o \
	$(FE)/material/uniaxial/PY/PyLiq1.o \
	$(FE)/material/uniaxial/PY/TzLiq1.o \
	$(FE)/material/uniaxial/PY/PySimple1Gen.o \
	$(FE)/material/uniaxial/PY/TzSimple1Gen.o

ZHY_LIBS = $(FE)/material/nD/soil/FluidSolidPorousMaterial.o \
	$(FE)/material/nD/soil/MultiYieldSurface.o \
	$(FE)/material/nD/soil/MultiYieldSurfaceClay.o \
	$(FE)/material/nD/soil/PressureDependMultiYield.o \
	$(FE)/material/nD/soil/PressureDependMultiYield02.o \
	$(FE)/material/nD/soil/PressureIndependMultiYield.o \
	$(FE)/material/nD/soil/T2Vector.o \
	$(FE)/element/UP-ucsd/BrickUP.o \
	$(FE)/element/UP-ucsd/BBarBrickUP.o \
	$(FE)/element/UP-ucsd/FourNodeQuadUP.o \
	$(FE)/element/UP-ucsd/BBarFourNodeQuadUP.o \
	$(FE)/element/UP-ucsd/Nine_Four_Node_QuadUP.o \
	$(FE)/element/UP-ucsd/Twenty_Eight_Node_BrickUP.o \
	$(FE)/element/UP-ucsd/shp3dv.o

ifeq ($(GRAPHICS), UsingOpenGL)

RENDERER_LIBS = $(FE)/renderer/OpenGlDevice.o \
	$(FE)/renderer/OpenGlRenderer.o \
	$(FE)/renderer/Renderer.o \
	$(FE)/renderer/PlainMap.o \
	$(FE)/recorder/FilePlotter.o \
	$(FE)/recorder/AlgorithmIncrements.o \

else

RENDERER_LIBS = $(FE)/renderer/db.o \
	$(FE)/renderer/View.o \
	$(FE)/renderer/Projection.o \
	$(FE)/renderer/Viewport.o \
	$(FE)/renderer/Scan.o \
	$(FE)/renderer/Clipping.o \
	$(FE)/renderer/X11Device.o \
	$(FE)/renderer/WindowDevice.o \
	$(FE)/renderer/Device.o \
	$(FE)/renderer/gMatrix.o \
	$(FE)/renderer/DofColorMap.o \
	$(FE)/renderer/Renderer.o \
	$(FE)/renderer/WindowRenderer.o \
	$(FE)/renderer/X11Renderer.o \
	$(FE)/renderer/PlainMap.o \
	$(FE)/recorder/FilePlotter.o \
	$(FE)/recorder/AlgorithmIncrements.o


endif

ifeq ($(GRAPHICS), NONE)

RENDERER_LIBS = $(FE)/renderer/PlainMap.o

endif

RECORDER_LIBS = $(FE)/recorder/Recorder.o \
	$(FE)/recorder/DatastoreRecorder.o \
	$(FE)/recorder/NodeRecorder.o \
	$(FE)/recorder/EnvelopeNodeRecorder.o \
	$(FE)/recorder/DriftRecorder.o \
	$(FE)/recorder/EnvelopeDriftRecorder.o \
	$(FE)/recorder/PatternRecorder.o \
	$(FE)/recorder/ElementRecorder.o \
	$(FE)/recorder/NormEnvelopeElementRecorder.o \
	$(FE)/recorder/NormElementRecorder.o \
	$(FE)/recorder/EnvelopeElementRecorder.o \
	$(FE)/recorder/response/Response.o \
	$(FE)/recorder/response/CompositeResponse.o \
	$(FE)/recorder/response/ElementResponse.o \
	$(FE)/recorder/response/MaterialResponse.o \
	$(FE)/recorder/response/FiberResponse.o \
	$(FE)/recorder/DamageRecorder.o \
	$(FE)/recorder/RemoveRecorder.o \
	$(FE)/recorder/PVDRecorder.o \
	$(FE)/recorder/GmshRecorder.o \
	$(FE)/recorder/MPCORecorder.o 


DATABASE_LIBS = $(FE)/database/FileDatastore.o \
	$(FE)/database/NEESData.o

MATRIX_LIBS   = $(FE)/matrix/Matrix.o \
	$(FE)/matrix/Vector.o \
	$(FE)/matrix/ID.o

TAGGED_LIBS =   $(FE)/tagged/TaggedObject.o \
	$(FE)/tagged/storage/ArrayOfTaggedObjects.o \
	$(FE)/tagged/storage/ArrayOfTaggedObjectsIter.o \
	$(FE)/tagged/storage/MapOfTaggedObjects.o \
	$(FE)/tagged/storage/MapOfTaggedObjectsIter.o

UTILITY_LIBS = $(FE)/utility/Timer.o \
	$(FE)/utility/SimulationInformation.o \
	$(FE)/utility/File.o \
	$(FE)/utility/FileIter.o \
	$(FE)/utility/PeerNGA.o \
	$(FE)/utility/StringContainer.o 


GRAPH_LIBS = $(FE)/graph/graph/DOF_Graph.o \
	$(FE)/graph/graph/VertexIter.o \
	$(FE)/graph/graph/Vertex.o \
	$(FE)/graph/graph/Graph.o \
	$(FE)/graph/graph/DOF_GroupGraph.o \
	$(FE)/graph/numberer/RCM.o \
	$(FE)/graph/numberer/AMDNumberer.o \
	$(FE)/graph/numberer/MyRCM.o \
	$(FE)/graph/numberer/GraphNumberer.o \
	$(FE)/graph/numberer/SimpleNumberer.o \
	$(FE)/graph/partitioner/Metis.o

ELE_LIBS   =  $(FE)/element/Element.o \
	$(FE)/element/WrapperElement.o \
	$(FE)/element/Information.o \
	$(FE)/element/ElementalLoad.o \
	$(FE)/element/truss/Truss.o \
	$(FE)/element/truss/TrussSection.o \
	$(FE)/element/truss/CorotTruss.o \
	$(FE)/element/truss/CorotTrussSection.o \
	$(FE)/element/truss/Truss2.o \
	$(FE)/element/truss/CorotTruss2.o \
	$(FE)/element/truss/N4BiaxialTruss.o \
	$(FE)/element/zeroLength/ZeroLengthContact2D.o \
	$(FE)/element/zeroLength/ZeroLengthContact3D.o \
	$(FE)/element/zeroLength/ZeroLengthContactNTS2D.o \
	$(FE)/element/zeroLength/ZeroLengthInterface2D.o \
	$(FE)/element/zeroLength/ZeroLengthRocking.o \
	$(FE)/element/surfaceLoad/SurfaceLoad.o \
	$(FE)/element/surfaceLoad/TriSurfaceLoad.o \
	$(FE)/element/UWelements/SimpleContact2D.o \
	$(FE)/element/UWelements/SimpleContact3D.o \
	$(FE)/element/UWelements/BeamContact2D.o \
	$(FE)/element/UWelements/BeamContact2Dp.o \
	$(FE)/element/UWelements/BeamContact3D.o \
	$(FE)/element/UWelements/BeamContact3Dp.o \
	$(FE)/element/UWelements/BeamEndContact3D.o \
	$(FE)/element/UWelements/BeamEndContact3Dp.o \
	$(FE)/element/UWelements/PileToe3D.o \
	$(FE)/element/HUelements/MultipleShearSpring.o \
	$(FE)/element/HUelements/MultipleNormalSpring.o \
	$(FE)/element/HUelements/KikuchiBearing.o \
	$(FE)/element/HUelements/DBESI0.o \
	$(FE)/element/HUelements/DBESI1.o \
	$(FE)/element/HUelements/YamamotoBiaxialHDR.o \
	$(FE)/element/XMUelements/AC3D8HexWithSensitivity.o \
	$(FE)/element/XMUelements/ASI3D8QuadWithSensitivity.o \
	$(FE)/element/XMUelements/AV3D4QuadWithSensitivity.o \
	$(FE)/element/XMUelements/VS3D4QuadWithSensitivity.o \
	$(FE)/element/feap/elmt01.o \
	$(FE)/element/feap/elmt02.o \
	$(FE)/element/feap/elmt03.o \
	$(FE)/element/feap/elmt04.o \
	$(FE)/element/feap/elmt05.o \
	$(FE)/element/feap/fillCommon.o \
	$(FE)/element/feap/getCommon.o \
	$(FE)/element/feap/dummyFeap.o \
	$(FE)/element/feap/fElement.o \
	$(FE)/element/feap/fElmt02.o \
	$(FE)/element/dispBeamColumn/AxEqDispBeamColumn2d.o \
	$(FE)/element/dispBeamColumn/TimoshenkoBeamColumn2d.o \
	$(FE)/element/dispBeamColumn/DispBeamColumn2d.o \
	$(FE)/element/dispBeamColumn/DispBeamColumnNL2d.o \
	$(FE)/element/dispBeamColumn/DispBeamColumn2dThermal.o \
	$(FE)/element/dispBeamColumn/DispBeamColumn3d.o \
	$(FE)/element/dispBeamColumn/DispBeamColumn3dThermal.o \
	$(FE)/element/dispBeamColumn/DispBeamColumn2dWithSensitivity.o \
	$(FE)/element/dispBeamColumn/DispBeamColumn3dWithSensitivity.o \
	$(FE)/element/dispBeamColumnInt/DispBeamColumn2dInt.o \
	$(FE)/element/dispBeamColumnInt/FiberSection2dInt.o \
	$(FE)/element/dispBeamColumnInt/LinearCrdTransf2dInt.o \
	$(FE)/element/elasticBeamColumn/ElasticBeam2d.o \
	$(FE)/element/elasticBeamColumn/WheelRail.o \
	$(FE)/element/elasticBeamColumn/ModElasticBeam2d.o \
	$(FE)/element/elasticBeamColumn/ElasticBeam3d.o \
	$(FE)/element/elasticBeamColumn/ElasticTimoshenkoBeam2d.o \
	$(FE)/element/elasticBeamColumn/ElasticTimoshenkoBeam3d.o \
	$(FE)/element/elasticBeamColumn/WheelRail.o \
	$(FE)/element/fourNodeQuad/FourNodeQuad.o \
	$(FE)/element/fourNodeQuad/FourNodeQuad3d.o \
	$(FE)/element/fourNodeQuad/FourNodeQuadWithSensitivity.o \
	$(FE)/element/fourNodeQuad/ConstantPressureVolumeQuad.o \
	$(FE)/element/fourNodeQuad/EnhancedQuad.o \
	$(FE)/element/componentElement/ComponentElement2d.o \
	$(FE)/element/fourNodeQuad/NineNodeMixedQuad.o \
	$(FE)/element/shell/ShellMITC4.o \
	$(FE)/element/shell/ShellMITC4Thermal.o \
	$(FE)/element/shell/ShellMITC9.o \
	$(FE)/element/shell/ShellDKGQ.o \
	$(FE)/element/shell/ShellNLDKGQ.o \
	$(FE)/element/shell/ShellDKGT.o \
	$(FE)/element/shell/ShellNLDKGT.o \
	$(FE)/element/shell/ShellNLDKGQThermal.o \
	$(FE)/element/shell/R3vectors.o \
	$(FE)/element/shell/ShellANDeS.o \
	$(FE)/element/brick/Brick.o \
	$(FE)/element/brick/BbarBrick.o \
	$(FE)/element/brick/BbarBrickWithSensitivity.o \
	$(FE)/element/brick/shp3d.o \
	$(FE)/element/brick/Twenty_Node_Brick.o \
	$(FE)/element/generic/GenericClient.o \
	$(FE)/element/generic/GenericCopy.o \
	$(FE)/element/adapter/ActuatorCorot.o \
	$(FE)/element/adapter/Actuator.o \
	$(FE)/element/adapter/Adapter.o \
	$(FE)/element/elastomericBearing/ElastomericBearingBoucWen2d.o \
	$(FE)/element/elastomericBearing/ElastomericBearingBoucWen3d.o \
	$(FE)/element/elastomericBearing/ElastomericBearingBoucWenMod3d.o \
	$(FE)/element/elastomericBearing/ElastomericBearingUFRP2d.o \
	$(FE)/element/elastomericBearing/ElastomericBearingPlasticity2d.o \
	$(FE)/element/elastomericBearing/ElastomericBearingPlasticity3d.o \
	$(FE)/element/elastomericBearing/HDR.o \
	$(FE)/element/elastomericBearing/LeadRubberX.o \
	$(FE)/element/elastomericBearing/ElastomericX.o \
	$(FE)/element/frictionBearing/MultiFP2d.o \
	$(FE)/element/frictionBearing/TripleFrictionPendulum.o \
	$(FE)/element/frictionBearing/TFP_Bearing.o \
	$(FE)/element/frictionBearing/TFP_Bearing2d.o \
	$(FE)/element/frictionBearing/TPB1D.o \
	$(FE)/element/frictionBearing/FPBearingPTV.o \
	$(FE)/element/frictionBearing/FlatSliderSimple2d.o \
	$(FE)/element/frictionBearing/FlatSliderSimple3d.o \
	$(FE)/element/frictionBearing/SingleFPSimple2d.o \
	$(FE)/element/frictionBearing/SingleFPSimple3d.o \
	$(FE)/element/frictionBearing/frictionModel/Coulomb.o \
	$(FE)/element/frictionBearing/frictionModel/FrictionModel.o \
	$(FE)/element/frictionBearing/frictionModel/FrictionResponse.o \
	$(FE)/element/frictionBearing/frictionModel/VelNormalFrcDep.o \
	$(FE)/element/frictionBearing/RJWatsonEQS2d.o \
	$(FE)/element/frictionBearing/RJWatsonEQS3d.o \
	$(FE)/element/frictionBearing/frictionModel/VelPressureDep.o \
	$(FE)/element/frictionBearing/frictionModel/VelDependent.o \
	$(FE)/element/frictionBearing/frictionModel/VelDepMultiLinear.o \
	$(FE)/element/twoNodeLink/TwoNodeLink.o \
	$(FE)/element/pyMacro/PY_Macro2D.o \
	$(FE)/element/triangle/Tri31.o \
	$(FE)/element/PFEMElement/PFEMElement2D.o \
	$(FE)/element/PFEMElement/PFEMElement3D.o \
	$(FE)/element/PFEMElement/PFEMElement2DCompressible.o \
	$(FE)/element/PFEMElement/PFEMElement2DBubble.o \
	$(FE)/element/PFEMElement/PFEMElement3DBubble.o \
	$(FE)/element/PFEMElement/PFEMElement2Dmini.o \
	$(FE)/element/PFEMElement/PFEMElement2DFIC.o \
	$(FE)/element/PFEMElement/PFEMElement2DQuasi.o \
	$(FE)/element/PFEMElement/MINI.o \
	$(FE)/element/PFEMElement/LineMesh.o \
	$(FE)/element/PFEMElement/TriMesh.o \
	$(FE)/element/PFEMElement/TetMesh.o \
	$(FE)/element/PFEMElement/QuadMesh.o \
	$(FE)/element/PFEMElement/TriangleMeshGenerator.o \
	$(FE)/element/PFEMElement/TetMeshGenerator.o \
	$(FE)/element/PFEMElement/LineMeshGenerator.o \
	$(FE)/element/PFEMElement/QuadMeshGenerator.o \
	$(FE)/element/PFEMElement/tetgen.o \
	$(FE)/element/PFEMElement/predicates.o \
	$(FE)/element/PFEMElement/triangle.o \
	$(FE)/element/PFEMElement/BackgroundMesh.o \
	$(FE)/element/PFEMElement/ParticleGroup.o \
	$(FE)/element/PFEMElement/Particle.o \
	$(FE)/element/PFEMElement/BackgroundDef.o \
	$(FE)/element/PFEMElement/TaylorHood2D.o \
	$(FE)/element/PFEMElement/TriGaussPoints.o \
	$(FE)/element/PFEMElement/HigherOrder.o \
	$(FE)/element/PFEMElement/Mesh.o \
	$(FE)/element/UWelements/SSPquad.o \
	$(FE)/element/UWelements/SSPquadUP.o \
	$(FE)/element/UWelements/SSPbrick.o \
	$(FE)/element/UWelements/SSPbrickUP.o \
	$(FE)/element/mvlem/MVLEM.o \
	$(FE)/element/mvlem/SFI_MVLEM.o \
	$(FE)/element/UWelements/Quad4FiberOverlay.o \
	$(FE)/element/UWelements/Brick8FiberOverlay.o \
	$(FE)/element/UWelements/QuadBeamEmbedContact.o \
	$(FE)/element/catenaryCable/CatenaryCable.o \
	$(FE)/element/tetrahedron/FourNodeTetrahedron.o \
	$(FE)/element/absorbentBoundaries/LysmerTriangle.o \
	$(FE)/element/UWelements/EmbeddedBeamInterfaceL.o \
	$(FE)/element/UWelements/EmbeddedBeamInterfaceP.o \
	$(FE)/element/UWelements/EmbeddedEPBeamInterface.o \
	$(FE)/element/UWelements/Tcl_generateInterfacePoints.o

#	$(FE)/material/nD/Damage2p.o \
#	$(FE)/material/nD/Damage2p3D.o \
#	$(FE)/material/nD/Damage2ppstress.o \


ifdef DAMAGE2P_FLAG
	DAMAGE_2P = $(FE)/material/nD/Damage2p.o $(FE)/material/nD/Damage2p3D.o $(FE)/material/nD/Damage2ppstress.o \
	$(FE)/material/nD/Concrete.o $(FE)/material/nD/CPlaneStress.o $(FE)/material/nD/CPlaneStress2d.o $(FE)/material/nD/CPlaneStrain.o $(FE)/material/nD/CThreeDimensional.o 

else
	DAMAGE_2P =
endif

ifdef DAMAGE2P_FLAG
	DAMAGE_2P = $(FE)/material/nD/Damage2p.o $(FE)/material/nD/Damage2p3D.o $(FE)/material/nD/Damage2ppstress.o \
	$(FE)/material/nD/Concrete.o $(FE)/material/nD/CPlaneStress.o $(FE)/material/nD/CPlaneStress2d.o $(FE)/material/nD/CPlaneStrain.o $(FE)/material/nD/CThreeDimensional.o 
else
	DAMAGE_2P =
endif

ifdef PLAINCONCRETE_FLAG
	PLAINCONCRETE = $(FE)/material/nD/PlainConcrete.o $(FE)/material/nD/forumat.o
else
	PLAINCONCRETE = 
endif

#	$(FE)/material/nD/Faria1998PlaneStress.o \
#	$(FE)/material/nD/Faria1998PlaneStress2d.o \
#	$(FE)/material/nD/Faria1998PlaneStrain.o \
#	$(FE)/material/nD/Faria1998.o \
#	$(FE)/material/nD/Faria1998ThreeDimensional.o \


MATERIAL_LIBS   =  $(FE)/material/Material.o \
	$(FE)/material/uniaxial/UniaxialMaterial.o \
	$(FE)/material/uniaxial/UniaxialJ2Plasticity.o \
	$(FE)/material/uniaxial/WrapperUniaxialMaterial.o \
	$(FE)/material/uniaxial/ElasticMaterial.o \
	$(FE)/material/uniaxial/Elastic2Material.o \
	$(FE)/material/uniaxial/ElasticMultiLinear.o \
	$(FE)/material/uniaxial/ReinforcingSteel.o \
	$(FE)/material/uniaxial/ElasticPPMaterial.o \
	$(FE)/material/uniaxial/ElasticBilin.o \
	$(FE)/material/uniaxial/MultiLinear.o \
	$(FE)/material/uniaxial/HookGap.o \
	$(FE)/material/uniaxial/HardeningMaterial.o \
	$(FE)/material/uniaxial/ParallelMaterial.o \
	$(FE)/material/uniaxial/SeriesMaterial.o \
	$(FE)/material/uniaxial/Neoprene.o \
	$(FE)/material/uniaxial/Concrete01.o \
	$(FE)/material/uniaxial/Concrete02.o \
	$(FE)/material/uniaxial/Concrete04.o \
	$(FE)/material/uniaxial/Concrete06.o \
	$(FE)/material/uniaxial/Concrete07.o \
	$(FE)/material/uniaxial/ConcreteSakaiKawashima.o \
	$(FE)/material/uniaxial/FRPConfinedConcrete.o \
	$(FE)/material/uniaxial/FRPConfinedConcrete02.o \
	$(FE)/material/uniaxial/ResilienceMaterialHR.o \
	$(FE)/material/uniaxial/ResilienceLow.o \
	$(FE)/material/uniaxial/CFSSSWP.o \
	$(FE)/material/uniaxial/CFSWSWP.o \
	$(FE)/material/uniaxial/TriMatrix.o \
	$(FE)/material/uniaxial/CubicSpline.o \
	$(FE)/material/uniaxial/c14-SK-M.o \
	$(FE)/material/uniaxial/Steel01.o \
	$(FE)/material/uniaxial/SteelMPF.o \
	$(FE)/material/uniaxial/ConcreteCM.o \
	$(FE)/material/uniaxial/RambergOsgoodSteel.o \
	$(FE)/material/uniaxial/Steel02.o \
	$(FE)/material/uniaxial/Steel4.o \
	$(FE)/material/uniaxial/Steel2.o \
	$(FE)/material/uniaxial/OriginCentered.o \
	$(FE)/material/uniaxial/Steel03.o \
	$(FE)/material/uniaxial/SteelMP.o \
	$(FE)/material/uniaxial/BWBN.o \
	$(FE)/material/uniaxial/Steel01Thermal.o \
	$(FE)/material/uniaxial/Steel02Thermal.o \
	$(FE)/material/uniaxial/Concrete02Thermal.o \
	$(FE)/material/uniaxial/Concrete02Thermal.o \
	$(FE)/material/uniaxial/ConcreteECThermal.o \
	$(FE)/material/uniaxial/ElasticMaterialThermal.o \
	$(FE)/material/uniaxial/StainlessECThermal.o \
	$(FE)/material/uniaxial/SteelECThermal.o \
	$(FE)/material/uniaxial/limitState/PinchingLimitStateMaterial.o \
	$(FE)/material/uniaxial/ModIMKPinching.o \
	$(FE)/material/uniaxial/ModIMKPinching02.o \
	$(FE)/material/uniaxial/ModIMKPeakOriented.o \
	$(FE)/material/uniaxial/ModIMKPeakOriented02.o \
	$(FE)/material/uniaxial/SmoothPSConcrete.o \
	$(FE)/material/uniaxial/Concrete01WithSITC.o \
	$(FE)/material/uniaxial/ConcretewBeta.o \
	$(FE)/material/uniaxial/ConcreteD.o \
	$(FE)/material/uniaxial/KikuchiAikenHDR.o \
	$(FE)/material/uniaxial/KikuchiAikenLRB.o \
	$(FE)/material/uniaxial/AxialSpHD.o \
	$(FE)/material/uniaxial/AxialSp.o \
	$(FE)/material/uniaxial/ECC01.o \
	$(FE)/material/uniaxial/HystereticMaterial.o \
	$(FE)/material/uniaxial/EPPGapMaterial.o \
	$(FE)/material/uniaxial/ViscousMaterial.o \
	$(FE)/material/uniaxial/Maxwell.o \
	$(FE)/material/uniaxial/ViscousDamper.o \
	$(FE)/material/uniaxial/DamperMaterial.o \
	$(FE)/material/uniaxial/BilinearOilDamper.o \
	$(FE)/material/uniaxial/Cast.o \
	$(FE)/material/uniaxial/DoddRestrepo.o \
	$(FE)/material/uniaxial/Dodd_Restrepo.o \
	$(FE)/material/uniaxial/DoddRestr.o \
	$(FE)/material/uniaxial/STEELDR.o \
	$(FE)/material/uniaxial/PathIndependentMaterial.o \
	$(FE)/material/uniaxial/BackboneMaterial.o \
	$(FE)/material/uniaxial/MinMaxMaterial.o \
	$(FE)/material/uniaxial/TensionOnlyMaterial.o \
	$(FE)/material/uniaxial/SimpleFractureMaterial.o \
	$(FE)/material/uniaxial/InitStrainMaterial.o \
	$(FE)/material/uniaxial/InitStressMaterial.o \
	$(FE)/material/uniaxial/pyUCLA.o \
	$(FE)/material/uniaxial/FatigueMaterial.o \
	$(FE)/material/uniaxial/SAWSMaterial.o \
	$(FE)/material/uniaxial/snap/Bilin.o \
	$(FE)/material/uniaxial/snap/Bilin02.o \
	$(FE)/material/uniaxial/ConfinedConcrete01.o \
	$(FE)/material/uniaxial/CableMaterial.o \
	$(FE)/material/uniaxial/ENTMaterial.o \
	$(FE)/material/uniaxial/BoucWenMaterial.o \
	$(FE)/material/uniaxial/BoucWenOriginal.o \
	$(FE)/material/uniaxial/Pinching4Material.o \
	$(FE)/material/uniaxial/ShearPanelMaterial.o \
	$(FE)/material/uniaxial/BarSlipMaterial.o \
	$(FE)/material/uniaxial/Bond_SP01.o \
	$(FE)/material/uniaxial/HyperbolicGapMaterial.o \
	$(FE)/material/uniaxial/ImpactMaterial.o \
	$(FE)/material/uniaxial/SteelBRB.o \
	$(FE)/material/uniaxial/GNGMaterial.o \
	$(FE)/material/uniaxial/SelfCenteringMaterial.o \
	$(FE)/material/uniaxial/backbone/HystereticBackbone.o \
	$(FE)/material/uniaxial/backbone/ArctangentBackbone.o \
	$(FE)/material/uniaxial/backbone/TrilinearBackbone.o \
	$(FE)/material/uniaxial/backbone/ReeseSoftClayBackbone.o \
	$(FE)/material/uniaxial/backbone/ReeseStiffClayBelowWS.o \
	$(FE)/material/uniaxial/backbone/ReeseSandBackbone.o \
	$(FE)/material/uniaxial/backbone/ManderBackbone.o \
	$(FE)/material/uniaxial/backbone/RaynorBackbone.o \
	$(FE)/material/nD/NDMaterial.o \
	$(FE)/material/nD/PlaneStressLayeredMaterial.o \
	$(FE)/material/nD/PlaneStressRebarMaterial.o \
	$(DAMAGE_2P) $(FE)/material/nD/WrapperNDMaterial.o \
	$(FE)/material/nD/J2Plasticity.o \
	$(FE)/material/nD/J2PlaneStress.o \
	$(FE)/material/nD/J2PlaneStrain.o \
	$(FE)/material/nD/J2AxiSymm.o \
	$(FE)/material/nD/J2ThreeDimensional.o \
	$(FE)/material/nD/J2PlateFiber.o \
	$(FE)/material/nD/J2PlateFibre.o \
	$(FE)/material/nD/PlaneStressRebarMaterial.o \
	$(FE)/material/nD/PlaneStressLayeredMaterial.o \
	$(FE)/material/nD/J2BeamFiber2d.o \
	$(FE)/material/nD/J2BeamFiber3d.o \
	$(FE)/material/nD/FSAM.o \
	$(FE)/material/nD/ConcreteMcftNonLinear5.o \
	$(FE)/material/nD/ConcreteMcftNonLinear7.o \
	$(FE)/material/nD/InitStressNDMaterial.o \
	$(FE)/material/nD/matCMM/MaterialCMM.o \
	$(FE)/material/nD/PlasticDamageConcrete3d.o \
	$(FE)/material/nD/PlasticDamageConcretePlaneStress.o \
	$(FE)/material/nD/eig3.o \
	$(FE)/material/nD/CycLiqCP.o \
	$(FE)/material/nD/CycLiqCP3D.o \
	$(FE)/material/nD/CycLiqCPPlaneStrain.o \
	$(FE)/material/nD/CycLiqCPSP.o \
	$(FE)/material/nD/CycLiqCPSP3D.o \
	$(FE)/material/nD/CycLiqCPSPPlaneStrain.o \
	$(FE)/material/nD/UWmaterials/DruckerPrager.o \
	$(FE)/material/nD/UWmaterials/DruckerPrager3D.o \
	$(FE)/material/nD/UWmaterials/DruckerPragerPlaneStrain.o \
	$(FE)/material/nD/UWmaterials/ContactMaterial3D.o \
	$(FE)/material/nD/UWmaterials/ContactMaterial2D.o \
	$(FE)/material/nD/UWmaterials/BoundingCamClay.o \
	$(FE)/material/nD/UWmaterials/BoundingCamClay3D.o \
	$(FE)/material/nD/UWmaterials/BoundingCamClayPlaneStrain.o \
	$(FE)/material/nD/UWmaterials/ManzariDafalias.o \
	$(FE)/material/nD/UWmaterials/ManzariDafalias3D.o \
	$(FE)/material/nD/UWmaterials/ManzariDafaliasPlaneStrain.o \
	$(FE)/material/nD/UWmaterials/ManzariDafaliasRO.o \
	$(FE)/material/nD/UWmaterials/ManzariDafalias3DRO.o \
	$(FE)/material/nD/UWmaterials/ManzariDafaliasPlaneStrainRO.o \
	$(FE)/material/nD/UWmaterials/InitialStateAnalysisWrapper.o \
	$(FE)/material/nD/UWmaterials/J2CyclicBoundingSurface.o \
	$(FE)/material/nD/UWmaterials/PM4Sand.o \
	$(FE)/material/nD/UWmaterials/PM4Silt.o \
	$(FE)/material/nD/LinearCap.o \
	$(FE)/material/nD/AcousticMedium.o \
	$(FE)/material/nD/FeapMaterial.o \
	$(FE)/material/nD/feap/FeapMaterial01.o \
	$(FE)/material/nD/feap/matl01.o \
	$(FE)/material/nD/feap/FeapMaterial02.o \
	$(FE)/material/nD/feap/matl02.o \
	$(FE)/material/nD/feap/FeapMaterial03.o \
	$(FE)/material/nD/feap/matl03.o \
	$(FE)/material/nD/feap/feapCommon.o \
	$(FE)/material/nD/PressureDependentElastic3D.o \
	$(FE)/material/nD/ElasticIsotropicMaterial.o \
	$(FE)/material/nD/ElasticIsotropicThreeDimensional.o \
	$(FE)/material/nD/ElasticIsotropicPlaneStress2D.o \
	$(FE)/material/nD/ElasticIsotropicPlaneStrain2D.o \
	$(FE)/material/nD/ElasticIsotropicAxiSymm.o \
	$(FE)/material/nD/ElasticIsotropicPlateFiber.o \
	$(FE)/material/nD/ElasticIsotropicBeamFiber.o \
	$(FE)/material/nD/ElasticIsotropicBeamFiber2d.o \
	$(FE)/material/nD/ElasticOrthotropicMaterial.o \
	$(FE)/material/nD/ElasticOrthotropicThreeDimensional.o \
	$(FE)/material/nD/PlaneStressMaterial.o \
	$(FE)/material/nD/PlaneStrainMaterial.o \
	$(FE)/material/nD/PlateFiberMaterial.o \
	$(FE)/material/nD/ElasticIsotropic3DThermal.o \
	$(FE)/material/nD/ElasticIsotropicMaterialThermal.o \
	$(FE)/material/nD/J2PlasticityThermal.o \
	$(FE)/material/nD/J2ThreeDimensionalThermal.o \
	$(FE)/material/nD/DruckerPrager3DThermal.o \
	$(FE)/material/nD/DruckerPragerThermal.o \
	$(FE)/material/nD/J2ThreeDimensionalThermal.o \
	$(FE)/material/nD/PlateFiberMaterialThermal.o \
	$(FE)/material/nD/PlateFromPlaneStressMaterialThermal.o \
	$(FE)/material/nD/PlateRebarMaterialThermal.o \
	$(FE)/material/nD/BeamFiberMaterial.o \
	$(FE)/material/nD/BeamFiberMaterial2d.o \
	$(FE)/material/nD/cyclicSoil/MultiaxialCyclicPlasticity.o \
	$(FE)/material/nD/cyclicSoil/MultiaxialCyclicPlasticityPlaneStrain.o \
	$(FE)/material/nD/cyclicSoil/MultiaxialCyclicPlasticity3D.o \
	$(FE)/material/nD/cyclicSoil/MultiaxialCyclicPlasticityAxiSymm.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/FAReinforcedConcretePlaneStress.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/ReinforcedConcretePlaneStress.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/FAFourSteelRCPlaneStress.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/RAFourSteelRCPlaneStress.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/PrestressedConcretePlaneStress.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/FAPrestressedConcretePlaneStress.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/FAFourSteelPCPlaneStress.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/RAFourSteelPCPlaneStress.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/ConcreteZ01.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/ConcreteL01.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/SteelZ01.o \
	$(FE)/material/nD/reinforcedConcretePlaneStress/TendonL01.o \
	$(FE)/material/nD/CapPlasticity.o \
	$(FE)/material/nD/ConcreteS.o \
	$(FE)/material/nD/SimplifiedJ2.o \
	$(FE)/material/nD/PlaneStressSimplifiedJ2.o \
	$(FE)/material/nD/PlaneStressUserMaterial.o $(PLAINCONCRETE) \
	$(FE)/material/nD/PlateFromPlaneStressMaterial.o \
	$(FE)/material/nD/PlateRebarMaterial.o \
	$(FE)/material/nD/stressDensityModel/stressDensity.o \
	$(FE)/material/nD/stressDensityModel/SDM-UC.o

FEDEAS_LIBS = $(FE)/material/uniaxial/FedeasMaterial.o \
	$(FE)/material/uniaxial/fedeas/FedeasHardeningMaterial.o \
	$(FE)/material/uniaxial/fedeas/Hard1.o \
	$(FE)/material/uniaxial/fedeas/PlasticDamageMaterial.o \
	$(FE)/material/uniaxial/fedeas/PD1.o \
	$(FE)/material/uniaxial/fedeas/FedeasBond1Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasBond2Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasConcr1Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasConcr2Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasConcr3Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasHyster1Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasHyster2Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasSteel1Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasSteel2Material.o


DRAIN_LIBS = $(FE)/material/uniaxial/DrainMaterial.o \
	$(FE)/material/uniaxial/drain/DrainHardeningMaterial.o \
	$(FE)/material/uniaxial/drain/anal00.o \
	$(FE)/material/uniaxial/drain/common00.o \
	$(FE)/material/uniaxial/drain/DrainBilinearMaterial.o \
	$(FE)/material/uniaxial/drain/DrainClough1Material.o \
	$(FE)/material/uniaxial/drain/DrainClough2Material.o \
	$(FE)/material/uniaxial/drain/DrainPinch1Material.o

SNAP_LIBS = \
	$(FE)/material/uniaxial/snap/Pinching.o \
	$(FE)/material/uniaxial/snap/Bilinear.o \
	$(FE)/material/uniaxial/snap/Clough.o \
	$(FE)/material/uniaxial/snap/CloughHenry.o \
	$(FE)/material/uniaxial/snap/CloughDamage.o \
	$(FE)/material/uniaxial/snap/PinchingDamage.o \


SECTION_LIBS = $(FE)/material/section/SectionForceDeformation.o \
	$(FE)/material/section/fiber/Fiber.o \
	$(FE)/material/section/fiber/UniaxialFiber2d.o \
	$(FE)/material/section/fiber/UniaxialFiber3d.o \
	$(FE)/material/section/fiber/NDFiber2d.o \
	$(FE)/material/section/fiber/NDFiber3d.o \
	$(FE)/material/section/ElasticSection2d.o \
	$(FE)/material/section/ElasticSection3d.o \
	$(FE)/material/section/ElasticShearSection2d.o \
	$(FE)/material/section/ElasticShearSection3d.o \
	$(FE)/material/section/ElasticWarpingShearSection2d.o \
	$(FE)/material/section/ElasticTubeSection3d.o \
	$(FE)/material/section/FiberSection2d.o \
	$(FE)/material/section/NDFiberSection2d.o \
	$(FE)/material/section/NDFiberSectionWarping2d.o \
	$(FE)/material/section/FiberSection2dThermal.o \
	$(FE)/material/section/LayeredShellFiberSectionThermal.o \
	$(FE)/material/section/FiberSection3dThermal.o \
	$(FE)/material/section/MembranePlateFiberSectionThermal.o \
	$(FE)/material/section/FiberSectionGJThermal.o \
	$(FE)/material/section/FiberSection3d.o \
	$(FE)/material/section/NDFiberSection3d.o \
	$(FE)/material/section/SectionAggregator.o \
	$(FE)/material/section/ParallelSection.o \
	$(FE)/material/section/GenericSection1d.o \
	$(FE)/material/section/ElasticPlateSection.o \
	$(FE)/material/section/ElasticMembranePlateSection.o \
	$(FE)/material/section/MembranePlateFiberSection.o \
	$(FE)/material/section/Bidirectional.o \
	$(FE)/material/section/Elliptical2.o \
	$(FE)/material/section/Isolator2spring.o \
	$(FE)/material/section/BiaxialHysteretic.o \
	$(FE)/material/section/LayeredShellFiberSection.o \
	$(FE)/material/section/integration/SectionIntegration.o \
	$(FE)/material/section/integration/RCSectionIntegration.o \
	$(FE)/material/section/integration/RCCircularSectionIntegration.o \
	$(FE)/material/section/integration/RCTunnelSectionIntegration.o \
	$(FE)/material/section/integration/WideFlangeSectionIntegration.o \
	$(FE)/material/section/integration/RCTBeamSectionIntegration.o \
	$(FE)/material/section/integration/TubeSectionIntegration.o \
	$(FE)/material/section/repres/patch/Patch.o \
	$(FE)/material/section/repres/patch/QuadPatch.o \
	$(FE)/material/section/repres/patch/CircPatch.o \
	$(FE)/material/section/repres/cell/Cell.o \
	$(FE)/material/section/repres/cell/QuadCell.o \
	$(FE)/material/section/repres/cell/CircSectionCell.o \
	$(FE)/material/section/repres/reinfBar/ReinfBar.o \
	$(FE)/material/section/repres/reinfLayer/ReinfLayer.o \
	$(FE)/material/section/repres/reinfLayer/StraightReinfLayer.o \
	$(FE)/material/section/repres/reinfLayer/CircReinfLayer.o \
	$(FE)/material/section/repres/section/SectionRepres.o \
	$(FE)/material/section/repres/section/FiberSectionRepr.o \
	$(FE)/material/section/yieldSurface/YieldSurfaceSection2d.o \
	$(FE)/material/section/yieldSurface/YS_Section2D01.o \
	$(FE)/material/section/yieldSurface/YS_Section2D02.o \
	$(FE)/material/section/yieldSurface/SoilFootingSection2d.o


SUPER_LU_OBJ = $(FE)/system_of_eqn/linearSOE/sparseGEN/SuperLU.o 

ifeq ($(PROGRAMMING_MODE), THREADS)
SUPER_LU_OBJ = $(FE)/system_of_eqn/linearSOE/sparseGEN/ThreadedSuperLU.o
endif

ifeq ($(PROGRAMMING_MODE), PARALLEL)
SUPER_LU_OBJ = $(FE)/system_of_eqn/linearSOE/sparseGEN/SuperLU.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/DistributedSuperLU.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenColLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.o \
	$(FE)/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.o \
	$(FE)/system_of_eqn/linearSOE/mumps/MumpsSOE.o \
	$(FE)/system_of_eqn/linearSOE/mumps/MumpsSolver.o
endif

# ifeq ($(PROGRAMMING_MODE), PARALLEL_INTERPRETERS)
# SUPER_LU_OBJ = $(FE)/system_of_eqn/linearSOE/sparseGEN/SuperLU.o \
# 	$(FE)/system_of_eqn/linearSOE/sparseGEN/DistributedSuperLU.o \
# 	$(FE)/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenColLinSOE.o 
# endif

CUDA_CLASSES = 

SequentialSysOfEqn_LIBS =	$(FE)/system_of_eqn/linearSOE/LinearSOE.o \
	$(FE)/system_of_eqn/linearSOE/LinearSOESolver.o \
	$(FE)/system_of_eqn/linearSOE/DomainSolver.o \
	$(FE)/system_of_eqn/linearSOE/bandGEN/BandGenLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/bandGEN/DistributedBandGenLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/bandGEN/BandGenLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/bandGEN/BandGenLinLapackSolver.o \
	$(FE)/system_of_eqn/linearSOE/fullGEN/FullGenLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/fullGEN/FullGenLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/fullGEN/FullGenLinLapackSolver.o \
	$(FE)/system_of_eqn/linearSOE/bandSPD/BandSPDLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/bandSPD/BandSPDLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/bandSPD/DistributedBandSPDLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/bandSPD/BandSPDLinLapackSolver.o \
	$(FE)/system_of_eqn/linearSOE/itpack/ItpackLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/itpack/ItpackLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/diagonal/DiagonalSOE.o \
	$(FE)/system_of_eqn/linearSOE/diagonal/DiagonalSolver.o \
	$(FE)/system_of_eqn/linearSOE/diagonal/DiagonalDirectSolver.o \
	$(FE)/system_of_eqn/linearSOE/diagonal/DistributedDiagonalSOE.o \
	$(FE)/system_of_eqn/linearSOE/diagonal/DistributedDiagonalSolver.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/SProfileSPDLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/SProfileSPDLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/DistributedProfileSPDLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectSolver.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSubstrSolver.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenColLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMSolver.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Umfpack.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleSolver.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMQuasiLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMQuasiSolver.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenColLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenRowLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenRowLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMDiaLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMDiaSolver.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Laplace.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_LumpM.o \
	$(SUPER_LU_OBJ) \
	$(FE)/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.o \
	$(FE)/system_of_eqn/eigenSOE/FullGenEigenSOE.o \
	$(FE)/system_of_eqn/eigenSOE/FullGenEigenSolver.o

ifeq ($(PROGRAMMING_MODE), PARALLEL_INTERPRETERS)
SysOfEqn_LIBS = $(SequentialSysOfEqn_LIBS) \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Mumps.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleSolver_Mumps.o
else
SysOfEqn_LIBS = $(SequentialSysOfEqn_LIBS)
endif

MODEL_BUILDER_LIBS = $(FE)/modelbuilder/ModelBuilder.o \
	$(FE)/modelbuilder/PlaneFrame.o \
	$(FE)/modelbuilder/tcl/Block2D.o \
	$(FE)/modelbuilder/tcl/Block3D.o

DOMAIN_LIBS =  $(FE)/domain/component/DomainComponent.o \
	$(FE)/domain/component/Parameter.o \
	$(FE)/domain/component/ElementParameter.o \
	$(FE)/domain/component/InitialStateParameter.o \
	$(FE)/domain/component/ElementStateParameter.o \
	$(FE)/domain/component/RVParameter.o \
	$(FE)/domain/component/DVParameter.o \
	$(FE)/domain/component/NodeResponseParameter.o \
	$(FE)/domain/component/LoadFactorParameter.o \
	$(FE)/domain/component/MaterialStageParameter.o \
	$(FE)/domain/component/MatParameter.o \
	$(FE)/domain/domain/Domain.o \
	$(FE)/domain/domain/single/SingleDomEleIter.o \
	$(FE)/domain/domain/single/SingleDomNodIter.o \
	$(FE)/domain/domain/single/SingleDomSP_Iter.o \
	$(FE)/domain/domain/single/SingleDomPC_Iter.o \
	$(FE)/domain/domain/single/SingleDomAllSP_Iter.o \
	$(FE)/domain/domain/single/SingleDomMP_Iter.o \
	$(FE)/domain/domain/single/SingleDomParamIter.o \
	$(FE)/domain/domain/partitioned/PartitionedDomain.o \
	$(FE)/domain/domain/partitioned/PartitionedDomainEleIter.o \
	$(FE)/domain/domain/partitioned/PartitionedDomainSubIter.o \
	$(FE)/domain/partitioner/DomainPartitioner.o \
	$(FE)/domain/region/MeshRegion.o \
	$(FE)/domain/node/Node.o \
	$(FE)/domain/node/NodalLoad.o \
	$(FE)/domain/constraints/SP_Constraint.o \
	$(FE)/domain/constraints/MP_Constraint.o \
	$(FE)/domain/constraints/Pressure_Constraint.o \
	$(FE)/domain/constraints/ImposedMotionSP.o \
	$(FE)/domain/constraints/ImposedMotionSP1.o \
	$(FE)/domain/constraints/RigidDiaphragm.o \
	$(FE)/domain/constraints/RigidRod.o \
	$(FE)/domain/constraints/RigidBeam.o \
	$(FE)/domain/load/NodalLoadIter.o \
	$(FE)/domain/load/ElementalLoadIter.o \
	$(FE)/domain/load/Load.o \
	$(FE)/domain/load/Beam2dUniformLoad.o \
	$(FE)/domain/load/Beam2dPartialUniformLoad.o \
	$(FE)/domain/load/Beam2dPointLoad.o \
	$(FE)/domain/load/Beam2dTempLoad.o \
	$(FE)/domain/load/Beam2dThermalAction.o \
	$(FE)/domain/load/NodalThermalAction.o \
	$(FE)/domain/load/ThermalActionWrapper.o \
	$(FE)/domain/load/Beam3dThermalAction.o \
	$(FE)/domain/load/ShellThermalAction.o \
	$(FE)/domain/load/Beam3dUniformLoad.o \
	$(FE)/domain/load/Beam3dPartialUniformLoad.o \
	$(FE)/domain/load/Beam3dPointLoad.o \
	$(FE)/domain/load/BrickSelfWeight.o \
	$(FE)/domain/load/SurfaceLoader.o \
	$(FE)/domain/load/LysmerVelocityLoader.o \
	$(FE)/domain/load/SelfWeight.o \
	$(FE)/domain/pattern/LoadPattern.o \
	$(FE)/domain/pattern/FireLoadPattern.o \
	$(FE)/domain/pattern/EarthquakePattern.o \
	$(FE)/domain/pattern/MultiSupportPattern.o \
	$(FE)/domain/pattern/UniformExcitation.o \
	$(FE)/domain/pattern/LoadPatternIter.o \
	$(FE)/domain/pattern/TimeSeries.o \
	$(FE)/domain/pattern/LinearSeries.o \
	$(FE)/domain/pattern/RectangularSeries.o \
	$(FE)/domain/pattern/ConstantSeries.o \
	$(FE)/domain/pattern/TrigSeries.o \
	$(FE)/domain/pattern/PathSeries.o \
	$(FE)/domain/pattern/PeerMotion.o \
	$(FE)/domain/pattern/PeerNGAMotion.o \
	$(FE)/domain/pattern/PathTimeSeries.o \
	$(FE)/domain/pattern/PathTimeSeriesThermal.o \
	$(FE)/domain/pattern/PulseSeries.o \
	$(FE)/domain/pattern/TriangleSeries.o \
	$(FE)/domain/pattern/TimeSeriesIntegrator.o \
	$(FE)/domain/pattern/TrapezoidalTimeSeriesIntegrator.o \
	$(FE)/domain/pattern/SimpsonTimeSeriesIntegrator.o \
	$(FE)/domain/pattern/drm/Mesh3DSubdomain.o \
	$(FE)/domain/pattern/drm/GeometricBrickDecorator.o \
	$(FE)/domain/pattern/drm/DRMBoundaryLayerDecorator.o \
	$(FE)/domain/pattern/drm/DRMLoadPattern.o \
	$(FE)/domain/pattern/drm/DRMLoadPatternWrapper.o \
	$(FE)/domain/pattern/drm/DRMInputHandler.o \
	$(FE)/domain/pattern/drm/PlaneDRMInputHandler.o \
	$(FE)/domain/groundMotion/GroundMotion.o \
	$(FE)/domain/groundMotion/GroundMotionRecord.o \
	$(FE)/domain/groundMotion/InterpolatedGroundMotion.o \
	$(FE)/domain/subdomain/Subdomain.o \
	$(FE)/domain/subdomain/ShadowSubdomain.o \
	$(FE)/domain/subdomain/ActorSubdomain.o \
	$(FE)/domain/subdomain/SubdomainNodIter.o \
	$(FE)/analysis/analysis/DomainUser.o

ANALYSIS_LIBS = $(FE)/analysis/analysis/Analysis.o \
	$(FE)/analysis/analysis/StaticAnalysis.o \
	$(FE)/analysis/analysis/TransientAnalysis.o \
	$(FE)/analysis/analysis/DirectIntegrationAnalysis.o \
	$(FE)/analysis/analysis/VariableTimeStepDirectIntegrationAnalysis.o \
	$(FE)/analysis/analysis/PFEMAnalysis.o \
	$(FE)/analysis/analysis/DomainDecompositionAnalysis.o \
	$(FE)/analysis/analysis/StaticDomainDecompositionAnalysis.o \
	$(FE)/analysis/analysis/TransientDomainDecompositionAnalysis.o \
	$(FE)/analysis/analysis/SubstructuringAnalysis.o \
	$(FE)/analysis/algorithm/SolutionAlgorithm.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/EquiSolnAlgo.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/Linear.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/NewtonRaphson.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/ModifiedNewton.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/NewtonLineSearch.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/NewtonHallM.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/Broyden.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/BFGS.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/KrylovNewton.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/PeriodicNewton.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/LineSearch.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/BisectionLineSearch.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/SecantLineSearch.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/RegulaFalsiLineSearch.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/InitialInterpolatedLineSearch.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/AcceleratedNewton.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/accelerator/Accelerator.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/accelerator/RaphsonAccelerator.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/accelerator/PeriodicAccelerator.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/accelerator/KrylovAccelerator.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/accelerator/SecantAccelerator1.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/accelerator/SecantAccelerator2.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/accelerator/SecantAccelerator3.o \
	$(FE)/convergenceTest/ConvergenceTest.o \
	$(FE)/convergenceTest/CTestNormUnbalance.o \
	$(FE)/convergenceTest/CTestNormDispIncr.o \
	$(FE)/convergenceTest/NormDispAndUnbalance.o \
	$(FE)/convergenceTest/NormDispOrUnbalance.o \
	$(FE)/convergenceTest/CTestEnergyIncr.o \
	$(FE)/convergenceTest/CTestRelativeNormUnbalance.o \
	$(FE)/convergenceTest/CTestRelativeNormDispIncr.o \
	$(FE)/convergenceTest/CTestRelativeEnergyIncr.o \
	$(FE)/convergenceTest/CTestRelativeTotalNormDispIncr.o \
	$(FE)/convergenceTest/CTestFixedNumIter.o \
	$(FE)/convergenceTest/CTestPFEM.o \
	$(FE)/analysis/algorithm/domainDecompAlgo/DomainDecompAlgo.o \
	$(FE)/analysis/integrator/Integrator.o \
	$(FE)/analysis/integrator/IncrementalIntegrator.o \
	$(FE)/analysis/integrator/StaticIntegrator.o \
	$(FE)/analysis/integrator/LoadControl.o \
	$(FE)/analysis/integrator/LoadPath.o \
	$(FE)/analysis/integrator/ArcLength.o \
	$(FE)/analysis/integrator/HSConstraint.o \
	$(FE)/analysis/integrator/ArcLength1.o \
	$(FE)/analysis/integrator/MinUnbalDispNorm.o \
	$(FE)/analysis/integrator/DisplacementControl.o \
	$(FE)/analysis/integrator/DistributedDisplacementControl.o \
	$(FE)/analysis/integrator/TransientIntegrator.o \
	$(FE)/analysis/integrator/Newmark.o \
	$(FE)/analysis/integrator/PFEMIntegrator.o \
	$(FE)/analysis/integrator/TRBDF2.o \
	$(FE)/analysis/integrator/TRBDF3.o \
	$(FE)/analysis/integrator/HHT.o \
	$(FE)/analysis/integrator/AlphaOS_TP.o \
	$(FE)/analysis/integrator/AlphaOSGeneralized_TP.o \
	$(FE)/analysis/integrator/HHT_TP.o \
	$(FE)/analysis/integrator/HHTExplicit_TP.o \
	$(FE)/analysis/integrator/HHTGeneralized_TP.o \
	$(FE)/analysis/integrator/HHTGeneralizedExplicit_TP.o \
	$(FE)/analysis/integrator/HHTHSFixedNumIter_TP.o \
	$(FE)/analysis/integrator/HHTHSIncrLimit_TP.o \
	$(FE)/analysis/integrator/HHTHSIncrReduct_TP.o \
	$(FE)/analysis/integrator/GeneralizedAlpha.o \
	$(FE)/analysis/integrator/Newmark1.o \
	$(FE)/analysis/integrator/CentralDifference.o \
	$(FE)/analysis/integrator/CentralDifferenceAlternative.o \
	$(FE)/analysis/integrator/CentralDifferenceNoDamping.o \
	$(FE)/analysis/integrator/WilsonTheta.o \
	$(FE)/analysis/integrator/ExplicitDifference.o \
	$(FE)/analysis/integrator/NewmarkExplicit.o \
	$(FE)/analysis/integrator/NewmarkHSIncrReduct.o \
	$(FE)/analysis/integrator/NewmarkHSIncrLimit.o \
	$(FE)/analysis/integrator/NewmarkHSFixedNumIter.o \
	$(FE)/analysis/integrator/HHTExplicit.o \
	$(FE)/analysis/integrator/HHTGeneralized.o \
	$(FE)/analysis/integrator/HHTGeneralized_TP.o \
	$(FE)/analysis/integrator/HHTGeneralizedExplicit.o \
	$(FE)/analysis/integrator/HHTHSIncrReduct.o \
	$(FE)/analysis/integrator/HHTHSIncrLimit.o \
	$(FE)/analysis/integrator/HHTHSFixedNumIter.o \
	$(FE)/analysis/integrator/AlphaOS.o \
	$(FE)/analysis/integrator/AlphaOS_TP.o \
	$(FE)/analysis/integrator/KRAlphaExplicit.o \
	$(FE)/analysis/integrator/KRAlphaExplicit_TP.o \
	$(FE)/analysis/integrator/AlphaOSGeneralized.o \
	$(FE)/analysis/integrator/AlphaOSGeneralized_TP.o \
	$(FE)/analysis/integrator/Collocation.o \
	$(FE)/analysis/integrator/CollocationHSIncrReduct.o \
	$(FE)/analysis/integrator/CollocationHSIncrLimit.o \
	$(FE)/analysis/integrator/CollocationHSFixedNumIter.o \
	$(FE)/analysis/integrator/KRAlphaExplicit_TP.o \
	$(FE)/analysis/integrator/HHT_TP.o \
	$(FE)/analysis/integrator/Houbolt.o \
	$(FE)/analysis/integrator/ParkLMS3.o \
	$(FE)/analysis/integrator/BackwardEuler.o \
	$(FE)/analysis/model/AnalysisModel.o \
	$(FE)/analysis/model/FE_EleIter.o \
	$(FE)/analysis/model/DOF_GrpIter.o \
	$(FE)/analysis/handler/ConstraintHandler.o \
	$(FE)/analysis/handler/PlainHandler.o \
	$(FE)/analysis/handler/PenaltyConstraintHandler.o \
	$(FE)/analysis/handler/LagrangeConstraintHandler.o \
	$(FE)/analysis/handler/TransformationConstraintHandler.o \
	$(FE)/analysis/numberer/DOF_Numberer.o \
	$(FE)/analysis/numberer/PlainNumberer.o \
	$(FE)/analysis/numberer/ParallelNumberer.o \
	$(FE)/analysis/dof_grp/DOF_Group.o \
	$(FE)/analysis/dof_grp/LagrangeDOF_Group.o \
	$(FE)/analysis/dof_grp/TransformationDOF_Group.o \
	$(FE)/analysis/fe_ele/FE_Element.o \
	$(FE)/analysis/fe_ele/penalty/PenaltySP_FE.o \
	$(FE)/analysis/fe_ele/penalty/PenaltyMP_FE.o \
	$(FE)/analysis/fe_ele/lagrange/LagrangeSP_FE.o \
	$(FE)/analysis/fe_ele/lagrange/LagrangeMP_FE.o \
	$(FE)/analysis/fe_ele/transformation/TransformationFE.o


ACTOR_LIBS = $(FE)/actor/channel/Channel.o \
	$(FE)/actor/channel/TCP_Socket.o \
	$(FE)/actor/channel/UDP_Socket.o \
	$(FE)/actor/channel/Socket.o \
	$(FE)/actor/channel/HTTP.o \
	$(FE)/actor/message/Message.o \
	$(FE)/actor/machineBroker/MachineBroker.o \
	$(FE)/actor/objectBroker/FEM_ObjectBroker.o \
	$(FE)/actor/objectBroker/FEM_ObjectBrokerAllClasses.o \
	$(FE)/actor/actor/Actor.o \
	$(FE)/actor/actor/MovableObject.o \
	$(FE)/actor/shadow/Shadow.o \
	$(FE)/actor/address/ChannelAddress.o

# Miscellaneous

tidy:
	@rm -f Makefile.bak *~ #*# core

clean:  tidy
	@rm -f $(OBJS) *.o

spotless: clean
	@rm -f $(PROGRAM) fake

wipe: spotless
	@$(CD) $(FE)/material; $(MAKE) wipe;
	@$(CD) $(FE)/analysis; $(MAKE) wipe;
	@$(CD) $(FE)/graph; $(MAKE) wipe;
	@$(CD) $(FE)/actor; $(MAKE) wipe;
	@$(CD) $(FE)/utility; $(MAKE) wipe;
	@$(CD) $(FE)/tcl; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn; $(MAKE) wipe;
	@$(CD) $(FE)/tagged; $(MAKE) wipe;
	@$(CD) $(FE)/domain; $(MAKE) wipe;
	@$(CD) $(FE)/modelbuilder; $(MAKE) wipe;
	@$(CD) $(FE)/matrix; $(MAKE) wipe;
	@$(CD) $(FE)/convergenceTest; $(MAKE) wipe;
	@$(CD) $(FE)/database; $(MAKE) wipe;
	@$(CD) $(FE)/recorder; $(MAKE) wipe;
	@$(CD) $(FE)/renderer; $(MAKE) wipe;
	@$(CD) $(FE)/handler; $(MAKE) wipe;
	@$(CD) $(FE)/machine; $(MAKE) wipe;
	@$(CD) $(FE)/coordTransformation; $(MAKE) wipe;
	@$(CD) $(FE)/element; $(MAKE) wipe;
	@$(CD) $(FE)/reliability; $(MAKE) wipe;
	@$(CD) $(FE)/optimization; $(MAKE) wipe;
	@$(CD) $(FE)/api; $(MAKE) wipe;
	@$(CD) $(FE)/interpreter; $(MAKE) wipe;
	@$(CD) $(FE)/damage; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.






