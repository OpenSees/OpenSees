# Makefile for fe objects

include ../Makefile.def

ifeq ($(RELIABILITY), YES_RELIABILITY)

all:
	@$(RM) $(RMFLAGS) $(FE_LIBRARY);
	@$(CD) $(FE)/matrix; $(MAKE);
	@$(CD) $(FE)/nDarray; $(MAKE);
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
	@$(CD) $(FE)/reliability; $(MAKE);
	@$(CD) $(FE)/coordTransformation; $(MAKE);
	@$(CD) $(FE)/damage; $(MAKE);
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(MATRIX_LIBS) $(ELE_LIBS) \
	$(MATERIAL_LIBS) $(HYSTERETIC_LIBS) $(SECTION_LIBS) $(DATABASE_LIBS) \
	$(FEDEAS_LIBS) $(MODEL_BUILDER_LIBS) $(DOMAIN_LIBS)\
	$(MATMOD_BJ_LIBS) $(NDARRAY_BJ_LIBS) $(ELEMENT_BJ_LIBS) \
	$(FINITE_DEFORMATION_BJ_LIBS) $(PY_SJB_RWB_BJ_LIBS) \
	$(ANALYSIS_LIBS) $(SysOfEqn_LIBS) $(GRAPH_LIBS) $(TAGGED_LIBS) \
	$(UTILITY_LIBS) $(ACTOR_LIBS) $(RECORDER_LIBS) $(LAW_LIBS) $(GGD_LIBS)\
	$(RENDERER_LIBS) $(REMO_LIBS) $(GLF_LIBS) $(ZHY_LIBS) $(ERRORHANDLER_LIBS) \
	$(DRAIN_LIBS) $(SNAP_LIBS) $(FMK_LIBS)
	@$(RANLIB) $(RANLIBFLAGS) $(FE_LIBRARY)

else

all:
	@$(RM) $(RMFLAGS) $(FE_LIBRARY);
	@$(CD) $(FE)/matrix; $(MAKE);
	@$(CD) $(FE)/nDarray; $(MAKE);
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
	@$(CD) $(FE)/damage; $(MAKE);
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(MATRIX_LIBS) $(ELE_LIBS) \
	$(MATERIAL_LIBS) $(HYSTERETIC_LIBS) $(SECTION_LIBS) $(DATABASE_LIBS) \
	$(FEDEAS_LIBS) $(MODEL_BUILDER_LIBS) $(DOMAIN_LIBS)\
	$(MATMOD_BJ_LIBS) $(NDARRAY_BJ_LIBS) $(ELEMENT_BJ_LIBS) \
	$(FINITE_DEFORMATION_BJ_LIBS) $(PY_SJB_RWB_BJ_LIBS) \
	$(ANALYSIS_LIBS) $(SysOfEqn_LIBS) $(GRAPH_LIBS) $(TAGGED_LIBS) \
	$(UTILITY_LIBS) $(ACTOR_LIBS) $(RECORDER_LIBS) $(LAW_LIBS) $(GGD_LIBS) \
	$(RENDERER_LIBS) $(REMO_LIBS) $(GLF_LIBS) $(ZHY_LIBS) $(ERRORHANDLER_LIBS) \
	$(DRAIN_LIBS) $(SNAP_LIBS) $(FMK_LIBS)
	@$(RANLIB) $(RANLIBFLAGS) $(FE_LIBRARY)

endif

skypack:
	@$(CD) $(HOME)/lib; rm -f $(HOME)/lib/libskypack.a;
	@$(CD) $(FE)/system_of_eqn/linearSOE/profileSPD/skypack; creator;

metis:
	@$(CD) $(HOME)/lib; rm -f $(METIS_LIBRARY);
	@$(CD) $(FE)/graph/partitioner/metis-2.0; $(MAKE) clean;
	@$(CD) $(FE)/graph/partitioner/metis-2.0; $(MAKE) lib;


FMK_LIBS = $(FE)/database/FE_Datastore.o


# the files that make up the FE_LIBRARY
GLF_LIBS = $(FE)/element/zeroLength/ZeroLength.o \
	$(FE)/element/zeroLength/ZeroLengthSection.o \
	$(FE)/element/zeroLength/TclZeroLength.o



REMO_LIBS = $(FE)/element/nonlinearBeamColumn/quadrule/QuadRule.o \
	$(FE)/element/nonlinearBeamColumn/quadrule/QuadRule1d.o \
	$(FE)/element/nonlinearBeamColumn/quadrule/QuadRule1d01.o \
	$(FE)/element/nonlinearBeamColumn/quadrule/GaussQuadRule1d.o \
	$(FE)/element/nonlinearBeamColumn/quadrule/GaussQuadRule1d01.o \
	$(FE)/element/nonlinearBeamColumn/quadrule/GaussLobattoQuadRule1d01.o \
	$(FE)/element/nonlinearBeamColumn/matrixutil/MatrixUtil.o \
	$(FE)/element/nonlinearBeamColumn/element/NLBeamColumn2d.o \
	$(FE)/element/nonlinearBeamColumn/element/NLBeamColumn3d.o \
	$(FE)/element/nonlinearBeamColumn/element/TclNLBeamColumnCommand.o \
	$(FE)/element/forceBeamColumn/ForceBeamColumn2d.o \
	$(FE)/element/forceBeamColumn/ForceBeamColumn3d.o \
	$(FE)/element/forceBeamColumn/TclForceBeamColumnCommand.o \
	$(FE)/element/forceBeamColumn/BeamIntegration.o \
	$(FE)/element/forceBeamColumn/LobattoBeamIntegration.o \
	$(FE)/element/forceBeamColumn/UserDefinedBeamIntegration.o \
	$(FE)/element/forceBeamColumn/HingeMidpointBeamIntegration2d.o \
	$(FE)/element/forceBeamColumn/HingeMidpointBeamIntegration3d.o \
	$(FE)/element/forceBeamColumn/HingeRadauBeamIntegration2d.o \
	$(FE)/element/forceBeamColumn/HingeRadauBeamIntegration3d.o \
	$(FE)/element/forceBeamColumn/HingeRadauTwoBeamIntegration2d.o \
	$(FE)/element/forceBeamColumn/HingeRadauTwoBeamIntegration3d.o \
	$(FE)/element/forceBeamColumn/UserDefinedHingeIntegration2d.o \
	$(FE)/element/forceBeamColumn/UserDefinedHingeIntegration3d.o \
	$(FE)/coordTransformation/CrdTransf.o \
	$(FE)/coordTransformation/CrdTransf2d.o \
	$(FE)/coordTransformation/CrdTransf3d.o \
	$(FE)/coordTransformation/LinearCrdTransf2d.o \
	$(FE)/coordTransformation/PDeltaCrdTransf2d.o \
	$(FE)/coordTransformation/LinearCrdTransf3d.o \
	$(FE)/coordTransformation/PDeltaCrdTransf3d.o \
	$(FE)/coordTransformation/TclGeomTransfCommand.o

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
	$(FE)/system_of_eqn/eigenSOE/BandArpackSOE.o \
	$(FE)/system_of_eqn/eigenSOE/BandArpackSolver.o \
	$(FE)/system_of_eqn/eigenSOE/SymArpackSOE.o \
	$(FE)/system_of_eqn/eigenSOE/SymArpackSolver.o \
	$(FE)/system_of_eqn/eigenSOE/SymBandEigenSOE.o \
	$(FE)/system_of_eqn/eigenSOE/SymBandEigenSolver.o \
	$(FE)/analysis/analysis/EigenAnalysis.o \
	$(FE)/analysis/integrator/EigenIntegrator.o \
	$(FE)/analysis/algorithm/eigenAlgo/EigenAlgorithm.o \
	$(FE)/analysis/algorithm/eigenAlgo/FrequencyAlgo.o \
	$(FE)/analysis/algorithm/eigenAlgo/StandardEigenAlgo.o


GGD_LIBS = $(FE)/element/joint/MP_Joint2D.o \
	$(FE)/element/joint/MP_Joint3D.o \
	$(FE)/element/joint/Joint2D.o \
	$(FE)/element/joint/Joint3D.o \
	$(FE)/element/joint/TclJoint2dCommand.o \
	$(FE)/element/joint/TclJoint3dCommand.o \
	$(FE)/element/joint/BeamColumnJoint2d.o \
	$(FE)/element/joint/BeamColumnJoint3d.o \
	$(FE)/element/joint/TclBeamColumnJointCommand.o \
	$(FE)/element/updatedLagrangianBeamColumn/Elastic2DGNL.o \
	$(FE)/element/updatedLagrangianBeamColumn/Inelastic2DYS02.o \
	$(FE)/element/updatedLagrangianBeamColumn/InelasticYS2DGNL.o \
	$(FE)/element/updatedLagrangianBeamColumn/TclElement2dYS.o \
	$(FE)/element/updatedLagrangianBeamColumn/Inelastic2DYS01.o \
	$(FE)/element/updatedLagrangianBeamColumn/Inelastic2DYS03.o \
	$(FE)/element/updatedLagrangianBeamColumn/TclElement2dGNL.o \
	$(FE)/element/updatedLagrangianBeamColumn/UpdatedLagrangianBeam2D.o \
	$(FE)/element/updatedLagrangianBeamColumn/BilinearCyclic.o \
	$(FE)/element/updatedLagrangianBeamColumn/CyclicModel.o \
	$(FE)/element/updatedLagrangianBeamColumn/LinearCyclic.o \
	$(FE)/element/updatedLagrangianBeamColumn/QuadraticCyclic.o \
	$(FE)/element/updatedLagrangianBeamColumn/TclCyclicModelCommands.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/Attalla2D.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/ElTawil2D.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/ElTawil2DUnSym.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/Hajjar2D.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/Orbison2D.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/YieldSurface_BC2D.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/YieldSurface_BC.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/NullYS2D.o \
	$(FE)/material/yieldSurface/yieldSurfaceBC/TclModelBuilderYieldSurfaceBCCommand.o \
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
	$(FE)/material/yieldSurface/evolution/TclModelBuilderYS_EvolutionCommand.o \
	$(FE)/material/yieldSurface/plasticHardeningMaterial/ExponReducing.o \
	$(FE)/material/yieldSurface/plasticHardeningMaterial/MultiLinearKp.o \
	$(FE)/material/yieldSurface/plasticHardeningMaterial/NullPlasticMaterial.o \
	$(FE)/material/yieldSurface/plasticHardeningMaterial/PlasticHardeningMaterial.o \
	$(FE)/material/yieldSurface/plasticHardeningMaterial/TclModelBuilderYSPlasticMaterialCommand.o \
	$(FE)/recorder/YsVisual.o \
	$(FE)/damage/DamageModel.o \
	$(FE)/damage/Kratzig.o \
	$(FE)/damage/ParkAng.o \
	$(FE)/damage/DamageResponse.o \
	$(FE)/damage/HystereticEnergy.o \
	$(FE)/damage/Mehanny.o \
	$(FE)/damage/NormalizedPeak.o \
	$(FE)/damage/TclModelBuilderDamageModelCommand.o \

ERRORHANDLER_LIBS = $(FE)/handler/StandardStream.o \
	$(FE)/handler/FileStream.o \
	$(FE)/handler/OPS_Stream.o


NDARRAY_BJ_LIBS = $(FE)/nDarray/basics.o \
                  $(FE)/nDarray/nDarray.o \
                  $(FE)/nDarray/BJmatrix.o \
                  $(FE)/nDarray/BJvector.o \
                  $(FE)/nDarray/BJtensor.o \
                  $(FE)/nDarray/straint.o  \
                  $(FE)/nDarray/stresst.o


MATMOD_BJ_LIBS = $(FE)/material/nD/ElasticIsotropic3D.o  \
                 $(FE)/material/nD/PressureDependentElastic3D.o \
                 $(FE)/material/nD/ElasticCrossAnisotropic.o \
                 $(FE)/material/nD/Template3Dep/MatPoint3D.o \
                 $(FE)/material/nD/Template3Dep/TclTemplate3DepCommand.o \
                 $(FE)/material/nD/Template3Dep/Template3Dep.o \
                 $(FE)/material/nD/Template3Dep/EPState.o \
                 $(FE)/material/nD/Template3Dep/PS.o \
                 $(FE)/material/nD/Template3Dep/YS.o \
                 $(FE)/material/nD/Template3Dep/MD_YS.o  \
                 $(FE)/material/nD/Template3Dep/MD_PS.o \
                 $(FE)/material/nD/Template3Dep/CAM_YS.o  \
                 $(FE)/material/nD/Template3Dep/CAM_PS.o \
                 $(FE)/material/nD/Template3Dep/DP_YS.o  \
                 $(FE)/material/nD/Template3Dep/DP_PS.o \
                 $(FE)/material/nD/Template3Dep/VM_YS.o  \
                 $(FE)/material/nD/Template3Dep/VM_PS.o  \
                 $(FE)/material/nD/Template3Dep/EL_T.o   \
                 $(FE)/material/nD/Template3Dep/EL_S.o   \
                 $(FE)/material/nD/Template3Dep/EL_NLEijMD.o \
                 $(FE)/material/nD/Template3Dep/EL_NLEp.o  \
                 $(FE)/material/nD/Template3Dep/EL_NLEeq.o  \
                 $(FE)/material/nD/Template3Dep/EL_LEeq.o  \
                 $(FE)/material/nD/Template3Dep/EL_NLEij.o \
                 $(FE)/material/nD/Template3Dep/EL_LEij.o \
                 $(FE)/material/nD/Template3Dep/RMC01.o  \
                 $(FE)/material/nD/Template3Dep/RMC01_YS.o  \
                 $(FE)/material/nD/Template3Dep/RMC01_PS.o



ELEMENT_BJ_LIBS = $(FE)/element/8nbrick/EightNodeBrick.o \
                  $(FE)/element/8nbrick/TclEightNodeBrickCommand.o \
                  $(FE)/element/20nbrick/TwentyNodeBrick.o \
                  $(FE)/element/20nbrick/TclTwentyNodeBrickCommand.o \
                  $(FE)/element/27nbrick/TwentySevenNodeBrick.o \
                  $(FE)/element/27nbrick/TclTwentySevenNodeBrickCommand.o \
                  $(FE)/element/TotalLagrangianFD20NodeBrick/TotalLagrangianFD20NodeBrick.o \
                  $(FE)/element/TotalLagrangianFD20NodeBrick/TclTLFD20NodeBrickCommand.o \
                  $(FE)/element/upU/EightNodeBrick_u_p_U.o \
                  $(FE)/element/upU/TwentyNodeBrick_u_p_U.o \
                  $(FE)/element/upU/TclEightNodeBrick_u_p_U.o \
                  $(FE)/element/upU/TclTwentyNodeBrick_u_p_U.o


FINITE_DEFORMATION_BJ_LIBS = $(FE)/material/nD/FiniteDeformation/NeoHookeanWEnergy.o \
                             $(FE)/material/nD/FiniteDeformation/OgdenWEnergy.o \
                             $(FE)/material/nD/FiniteDeformation/SimoPisterWEnergy.o \
                             $(FE)/material/nD/FiniteDeformation/FiniteDeformationElastic3D.o \
                             $(FE)/material/nD/FiniteDeformation/LogWEnergy.o \
                             $(FE)/material/nD/FiniteDeformation/MooneyRivlinWEnergy.o \
                             $(FE)/material/nD/FiniteDeformation/OgdenSimoWEnergy.o \
                             $(FE)/material/nD/FiniteDeformation/MooneyRivlinSimoWEnergy.o \
			     $(FE)/material/nD/FiniteDeformation/W.o \
			     $(FE)/material/nD/FiniteDeformation/NeoHookeanCompressible3D.o \
			     $(FE)/material/nD/FiniteDeformation/FDdecoupledElastic3D.o \
			     $(FE)/material/nD/FiniteDeformation/TclFiniteDeformationElastic3DCommand.o \
			     $(FE)/material/nD/FiniteDeformation/FDEPState.o \
			     $(FE)/material/nD/FiniteDeformation/FiniteDeformationEP3D.o \
			     $(FE)/material/nD/FiniteDeformation/fdEvolution/fdEvolution_S.o \
                             $(FE)/material/nD/FiniteDeformation/fdEvolution/fdEvolution_SLS.o \
			     $(FE)/material/nD/FiniteDeformation/fdEvolution/fdEvolution_T.o \
			     $(FE)/material/nD/FiniteDeformation/fdEvolution/fdEvolution_TL.o \
                             $(FE)/material/nD/FiniteDeformation/fdFlow/fdFlow.o \
                             $(FE)/material/nD/FiniteDeformation/fdFlow/fdFlowVM.o \
			     $(FE)/material/nD/FiniteDeformation/fdFlow/fdFlowDP.o \
                             $(FE)/material/nD/FiniteDeformation/fdYield/fdYield.o \
			     $(FE)/material/nD/FiniteDeformation/fdYield/fdYieldVM.o \
			     $(FE)/material/nD/FiniteDeformation/fdYield/fdYieldDP.o \
			     $(FE)/material/nD/FiniteDeformation/TclFiniteDeformationEP3DCommand.o





PY_SJB_RWB_BJ_LIBS = $(FE)/material/uniaxial/PY/PySimple1.o \
                     $(FE)/material/uniaxial/PY/QzSimple1.o \
                     $(FE)/material/uniaxial/PY/TzSimple1.o \
                     $(FE)/material/uniaxial/PY/PyLiq1.o \
                     $(FE)/material/uniaxial/PY/TzLiq1.o \
                     $(FE)/material/uniaxial/PY/TclPyTzQzMaterialCommand.o \
                     $(FE)/material/uniaxial/PY/PySimple1Gen.o \
                     $(FE)/material/uniaxial/PY/TzSimple1Gen.o



ZHY_LIBS =	$(FE)/material/nD/soil/FluidSolidPorousMaterial.o \
	$(FE)/material/nD/soil/MultiYieldSurface.o \
	$(FE)/material/nD/soil/PressureDependMultiYield.o \
	$(FE)/material/nD/soil/PressureDependMultiYield02.o \
	$(FE)/material/nD/soil/TclUpdateMaterialStageCommand.o \
	$(FE)/material/nD/soil/PressureIndependMultiYield.o \
	$(FE)/material/nD/soil/T2Vector.o \
	$(FE)/element/upU/BrickUP.o \
	$(FE)/element/upU/FourNodeQuadUP.o \
	$(FE)/element/upU/Nine_Four_Node_QuadUP.o \
	$(FE)/element/upU/Twenty_Eight_Node_BrickUP.o \
	$(FE)/element/upU/shp3dv.o \
	$(FE)/element/upU/TclFourNodeQuadUPCommand.o \
	$(FE)/element/20nbrick/TclTwenty_Node_BrickCommand.o \
	$(FE)/element/20nbrick/Twenty_Node_Brick.o 


ifeq ($(GRAPHICS), UsingOpenGL)

RENDERER_LIBS = $(FE)/renderer/OpenGlDevice.o \
	$(FE)/renderer/Renderer.o \
	$(FE)/renderer/PlainMap.o \
	$(FE)/renderer/OpenGlRenderer.o

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
	$(FE)/renderer/PlainMap.o

endif

RECORDER_LIBS = $(FE)/recorder/MaxNodeDispRecorder.o \
	$(FE)/recorder/DatastoreRecorder.o \
	$(FE)/recorder/NodeRecorder.o \
	$(FE)/recorder/EnvelopeNodeRecorder.o \
	$(FE)/recorder/DriftRecorder.o \
	$(FE)/recorder/PatternRecorder.o \
	$(FE)/recorder/ElementRecorder.o \
	$(FE)/recorder/EnvelopeElementRecorder.o \
	$(FE)/recorder/FilePlotter.o \
	$(FE)/recorder/GSA_Recorder.o \
	$(FE)/recorder/AlgorithmIncrements.o \
	$(FE)/recorder/TclRecorderCommands.o \
	$(FE)/recorder/response/Response.o \
	$(FE)/recorder/response/ElementResponse.o \
	$(FE)/recorder/response/MaterialResponse.o \
	$(FE)/recorder/response/FiberResponse.o \
	$(FE)/recorder/DamageRecorder.o 



DATABASE_LIBS = $(FE)/database/FileDatastore.o \
	$(FE)/database/MySqlDatastore.o \
	$(FE)/database/NEESData.o \
	$(FE)/database/TclDatabaseCommands.o

MATRIX_LIBS   = $(FE)/matrix/Matrix.o \
	$(FE)/matrix/Vector.o \
	$(FE)/matrix/ID.o  \

TAGGED_LIBS =   $(FE)/tagged/TaggedObject.o \
	$(FE)/tagged/storage/ArrayOfTaggedObjects.o \
	$(FE)/tagged/storage/ArrayOfTaggedObjectsIter.o  \
	$(FE)/tagged/storage/MapOfTaggedObjects.o \
	$(FE)/tagged/storage/MapOfTaggedObjectsIter.o

UTILITY_LIBS = $(FE)/utility/Timer.o


GRAPH_LIBS = $(FE)/graph/graph/DOF_Graph.o \
	$(FE)/graph/graph/VertexIter.o \
	$(FE)/graph/graph/Vertex.o  \
	$(FE)/graph/graph/Graph.o \
	$(FE)/graph/graph/DOF_GroupGraph.o  \
	$(FE)/graph/numberer/RCM.o \
	$(FE)/graph/numberer/MyRCM.o \
	$(FE)/graph/numberer/GraphNumberer.o \
	$(FE)/graph/numberer/SimpleNumberer.o \
	$(FE)/graph/partitioner/Metis.o

ELE_LIBS   =  $(FE)/element/Element.o \
	$(FE)/element/Information.o \
	$(FE)/element/ElementalLoad.o \
	$(FE)/element/beam2d/beam2d03.o \
	$(FE)/element/beam2d/beam2d04.o \
	$(FE)/element/beam3d/beam3d01.o \
	$(FE)/element/beam3d/beam3d02.o \
	$(FE)/element/truss/Truss.o \
	$(FE)/element/truss/TrussSection.o \
	$(FE)/element/truss/CorotTruss.o \
	$(FE)/element/truss/CorotTrussSection.o \
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
	$(FE)/element/TclElementCommands.o \
	$(FE)/element/feap/TclFeapElementCommand.o \
	$(FE)/element/truss/TclTrussCommand.o \
	$(FE)/element/beamWithHinges/TclBeamWithHingesBuilder.o \
	$(FE)/element/beamWithHinges/BeamWithHinges2d.o \
	$(FE)/element/beamWithHinges/BeamWithHinges3d.o \
	$(FE)/element/dispBeamColumn/DispBeamColumn2d.o \
	$(FE)/element/dispBeamColumn/DispBeamColumn3d.o \
	$(FE)/element/dispBeamColumn/TclDispBeamColumnCommand.o \
	$(FE)/element/elasticBeamColumn/ElasticBeam2d.o \
	$(FE)/element/elasticBeamColumn/ElasticBeam3d.o \
	$(FE)/element/elasticBeamColumn/TclElasticBeamCommand.o \
	$(FE)/element/fourNodeQuad/TclFourNodeQuadCommand.o \
	$(FE)/element/fourNodeQuad/FourNodeQuad.o \
	$(FE)/element/fourNodeQuad/ConstantPressureVolumeQuad.o \
	$(FE)/element/fourNodeQuad/EnhancedQuad.o \
	$(FE)/element/fourNodeQuad/NineNodeMixedQuad.o \
	$(FE)/element/shell/ShellMITC4.o \
	$(FE)/element/shell/R3vectors.o \
	$(FE)/element/shell/TclShellCommand.o \
	$(FE)/element/brick/Brick.o \
	$(FE)/element/brick/BbarBrick.o \
	$(FE)/element/brick/shp3d.o \
	$(FE)/element/brick/TclBrickCommand.o


MATERIAL_LIBS   =  $(FE)/material/Material.o \
	$(FE)/material/uniaxial/UniaxialMaterial.o \
	$(FE)/material/uniaxial/ElasticMaterial.o \
	$(FE)/material/uniaxial/ElasticPPMaterial.o \
	$(FE)/material/uniaxial/HardeningMaterial.o \
	$(FE)/material/uniaxial/ParallelMaterial.o \
	$(FE)/material/uniaxial/SeriesMaterial.o \
	$(FE)/material/uniaxial/Concrete01.o \
	$(FE)/material/uniaxial/Steel01.o \
	$(FE)/material/uniaxial/HystereticMaterial.o \
	$(FE)/material/uniaxial/EPPGapMaterial.o \
	$(FE)/material/uniaxial/ViscousMaterial.o \
	$(FE)/material/uniaxial/PathIndependentMaterial.o \
	$(FE)/material/uniaxial/MinMaxMaterial.o \
	$(FE)/material/uniaxial/FatigueMaterial.o \
	$(FE)/material/uniaxial/CableMaterial.o \
	$(FE)/material/uniaxial/ENTMaterial.o \
	$(FE)/material/uniaxial/BoucWenMaterial.o \
	$(FE)/material/uniaxial/Pinching4Material.o \
	$(FE)/material/uniaxial/BarSlipMaterial.o \
	$(FE)/material/nD/NDMaterial.o \
	$(FE)/material/nD/J2Plasticity.o \
	$(FE)/material/nD/J2PlaneStress.o \
	$(FE)/material/nD/J2PlaneStrain.o \
	$(FE)/material/nD/J2AxiSymm.o \
	$(FE)/material/nD/J2ThreeDimensional.o \
	$(FE)/material/nD/J2PlateFiber.o \
	$(FE)/material/nD/FeapMaterial.o \
	$(FE)/material/nD/feap/TclFeapMaterialCommand.o \
	$(FE)/material/nD/feap/FeapMaterial01.o \
	$(FE)/material/nD/feap/matl01.o \
	$(FE)/material/nD/feap/FeapMaterial02.o \
	$(FE)/material/nD/feap/matl02.o \
	$(FE)/material/nD/feap/FeapMaterial03.o \
	$(FE)/material/nD/feap/matl03.o \
	$(FE)/material/nD/feap/feapCommon.o \
	$(FE)/material/nD/ElasticIsotropicMaterial.o \
	$(FE)/material/nD/ElasticIsotropicPlaneStress2D.o \
	$(FE)/material/nD/ElasticIsotropicPlaneStrain2D.o \
	$(FE)/material/nD/ElasticIsotropicAxiSymm.o \
	$(FE)/material/nD/ElasticIsotropicPlateFiber.o \
	$(FE)/material/nD/ElasticIsotropicBeamFiber.o \
	$(FE)/material/nD/PlaneStressMaterial.o \
	$(FE)/material/nD/PlateFiberMaterial.o \
	$(FE)/material/nD/BeamFiberMaterial.o \
	$(FE)/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.o \
	$(FE)/material/nD/TclModelBuilderNDMaterialCommand.o

FEDEAS_LIBS = $(FE)/material/uniaxial/FedeasMaterial.o \
	$(FE)/material/uniaxial/fedeas/FedeasHardeningMaterial.o \
	$(FE)/material/uniaxial/fedeas/Hard1.o \
	$(FE)/material/uniaxial/fedeas/FedeasBond1Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasBond2Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasConcr1Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasConcr2Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasConcr3Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasHyster1Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasHyster2Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasSteel1Material.o \
	$(FE)/material/uniaxial/fedeas/FedeasSteel2Material.o \
	$(FE)/material/uniaxial/fedeas/TclFedeasMaterialCommand.o


DRAIN_LIBS = $(FE)/material/uniaxial/DrainMaterial.o \
	$(FE)/material/uniaxial/drain/DrainHardeningMaterial.o \
	$(FE)/material/uniaxial/drain/anal00.o \
	$(FE)/material/uniaxial/drain/common00.o \
	$(FE)/material/uniaxial/drain/DrainBilinearMaterial.o \
	$(FE)/material/uniaxial/drain/DrainClough1Material.o \
	$(FE)/material/uniaxial/drain/DrainClough2Material.o \
	$(FE)/material/uniaxial/drain/DrainPinch1Material.o \
	$(FE)/material/uniaxial/drain/TclDrainMaterialCommand.o

SNAP_LIBS = $(FE)/material/uniaxial/snap/TclSnapMaterialCommand.o \
	$(FE)/material/uniaxial/snap/Pinching.o \
	$(FE)/material/uniaxial/snap/Bilinear.o \
	$(FE)/material/uniaxial/snap/Clough.o \
	$(FE)/material/uniaxial/snap/CloughDamage.o \
	$(FE)/material/uniaxial/snap/PinchingDamage.o \

SECTION_LIBS = $(FE)/material/section/SectionForceDeformation.o \
	$(FE)/material/section/ElasticSection2d.o \
	$(FE)/material/section/ElasticSection3d.o \
	$(FE)/material/section/FiberSection2d.o \
	$(FE)/material/section/FiberSection3d.o \
	$(FE)/material/section/FiberSectionGJ.o \
	$(FE)/material/section/SectionAggregator.o \
	$(FE)/material/section/GenericSection1d.o \
	$(FE)/material/section/ElasticPlateSection.o \
	$(FE)/material/section/ElasticMembranePlateSection.o \
	$(FE)/material/section/MembranePlateFiberSection.o \
	$(FE)/material/section/Bidirectional.o \
	$(FE)/material/section/fiber/Fiber.o \
	$(FE)/material/section/fiber/UniaxialFiber2d.o \
	$(FE)/material/section/fiber/UniaxialFiber3d.o \
	$(FE)/material/section/repres/patch/Patch.o \
	$(FE)/material/section/repres/patch/QuadPatch.o \
	$(FE)/material/section/repres/patch/CircPatch.o \
	$(FE)/material/section/repres/cell/Cell.o \
	$(FE)/material/section/repres/cell/QuadCell.o \
	$(FE)/material/section/repres/reinfBar/ReinfBar.o \
	$(FE)/material/section/repres/reinfLayer/ReinfLayer.o \
	$(FE)/material/section/repres/reinfLayer/StraightReinfLayer.o \
	$(FE)/material/section/repres/reinfLayer/CircReinfLayer.o \
	$(FE)/material/section/repres/section/SectionRepres.o \
	$(FE)/material/section/repres/section/FiberSectionRepr.o \
	$(FE)/material/section/yieldSurface/YieldSurfaceSection2d.o \
	$(FE)/material/section/yieldSurface/YS_Section2D01.o \
	$(FE)/material/section/yieldSurface/YS_Section2D02.o \
	$(FE)/material/section/yieldSurface/TclModelBuilderYS_SectionCommand.o \
	$(FE)/material/section/TclModelBuilderSectionCommand.o


ifeq ($(PROGRAMMING_MODE), THREADS)
SUPER_LU_OBJ = $(FE)/system_of_eqn/linearSOE/sparseGEN/ThreadedSuperLU.o
else
SUPER_LU_OBJ = $(FE)/system_of_eqn/linearSOE/sparseGEN/SuperLU.o
endif

SysOfEqn_LIBS = $(FE)/system_of_eqn/SystemOfEqn.o \
	$(FE)/system_of_eqn/Solver.o \
	$(FE)/system_of_eqn/linearSOE/LinearSOE.o \
	$(FE)/system_of_eqn/linearSOE/LinearSOESolver.o \
	$(FE)/system_of_eqn/linearSOE/DomainSolver.o \
	$(FE)/system_of_eqn/linearSOE/bandGEN/BandGenLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/bandGEN/BandGenLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/bandGEN/BandGenLinLapackSolver.o \
	$(FE)/system_of_eqn/linearSOE/fullGEN/FullGenLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/fullGEN/FullGenLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/fullGEN/FullGenLinLapackSolver.o \
	$(FE)/system_of_eqn/linearSOE/bandSPD/BandSPDLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/bandSPD/BandSPDLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/bandSPD/BandSPDLinLapackSolver.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectSolver.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectSkypackSolver.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectBlockSolver.o \
	$(FE)/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinSubstrSolver.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenColLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/sparseGEN/SparseGenColLinSolver.o \
	$(SUPER_LU_OBJ) \
	$(FE)/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.o


MODEL_BUILDER_LIBS = $(FE)/modelbuilder/ModelBuilder.o \
	$(FE)/modelbuilder/PlaneFrame.o \
	$(FE)/modelbuilder/tcl/TclUniaxialMaterialTester.o \
	$(FE)/modelbuilder/tcl/TclModelBuilder.o \
	$(FE)/modelbuilder/tcl/Block2D.o \
	$(FE)/modelbuilder/tcl/Block3D.o

DOMAIN_LIBS =  $(FE)/domain/component/DomainComponent.o \
	$(FE)/domain/domain/Domain.o \
	$(FE)/domain/domain/single/SingleDomEleIter.o \
	$(FE)/domain/domain/single/SingleDomNodIter.o \
	$(FE)/domain/domain/single/SingleDomSP_Iter.o \
	$(FE)/domain/domain/single/SingleDomAllSP_Iter.o \
	$(FE)/domain/domain/single/SingleDomMP_Iter.o \
	$(FE)/domain/domain/partitioned/PartitionedDomain.o \
	$(FE)/domain/domain/partitioned/PartitionedDomainEleIter.o \
	$(FE)/domain/domain/partitioned/PartitionedDomainSubIter.o \
	$(FE)/domain/partitioner/DomainPartitioner.o \
	$(FE)/domain/region/MeshRegion.o \
	$(FE)/domain/region/TclRegionCommands.o \
	$(FE)/domain/node/Node.o \
	$(FE)/domain/node/NodalLoad.o \
	$(FE)/domain/constraints/SP_Constraint.o \
	$(FE)/domain/constraints/MP_Constraint.o \
	$(FE)/domain/constraints/ImposedMotionSP.o \
	$(FE)/domain/constraints/ImposedMotionSP1.o \
	$(FE)/domain/constraints/RigidDiaphragm.o \
	$(FE)/domain/constraints/RigidRod.o \
	$(FE)/domain/constraints/RigidBeam.o \
	$(FE)/domain/load/NodalLoadIter.o \
	$(FE)/domain/load/ElementalLoadIter.o \
	$(FE)/domain/load/Load.o \
	$(FE)/domain/load/Beam2dUniformLoad.o \
	$(FE)/domain/load/Beam2dPointLoad.o \
	$(FE)/domain/load/Beam2dTempLoad.o \
	$(FE)/domain/load/Beam3dUniformLoad.o \
	$(FE)/domain/load/Beam3dPointLoad.o \
	$(FE)/domain/load/BrickSelfWeight.o \
	$(FE)/domain/pattern/LoadPattern.o \
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
	$(FE)/domain/pattern/PathTimeSeries.o \
	$(FE)/domain/pattern/TclPatternCommand.o \
	$(FE)/domain/pattern/TclSeriesCommand.o \
	$(FE)/domain/pattern/TclSeriesIntegratorCommand.o \
	$(FE)/domain/pattern/TimeSeriesIntegrator.o \
	$(FE)/domain/pattern/TrapezoidalTimeSeriesIntegrator.o \
	$(FE)/domain/groundMotion/GroundMotion.o \
	$(FE)/domain/groundMotion/GroundMotionRecord.o \
	$(FE)/domain/groundMotion/InterpolatedGroundMotion.o \
	$(FE)/domain/groundMotion/TclGroundMotionCommand.o \
	$(FE)/domain/subdomain/Subdomain.o \
	$(FE)/domain/subdomain/SubdomainNodIter.o \
	$(FE)/analysis/analysis/DomainUser.o \
	$(FE)/domain/pattern/PBowlLoading.o


ANALYSIS_LIBS = $(FE)/analysis/analysis/Analysis.o \
	$(FE)/analysis/analysis/StaticAnalysis.o \
	$(FE)/analysis/analysis/TransientAnalysis.o \
	$(FE)/analysis/analysis/DirectIntegrationAnalysis.o \
	$(FE)/analysis/analysis/VariableTimeStepDirectIntegrationAnalysis.o \
	$(FE)/analysis/analysis/DomainDecompositionAnalysis.o \
	$(FE)/analysis/analysis/SubstructuringAnalysis.o \
	$(FE)/analysis/algorithm/SolutionAlgorithm.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/EquiSolnAlgo.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/Linear.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/NewtonRaphson.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/ModifiedNewton.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/NewtonLineSearch.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/Broyden.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/BFGS.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/KrylovNewton.o \
        $(FE)/analysis/algorithm/equiSolnAlgo/PeriodicNewton.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/LineSearch.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/BisectionLineSearch.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/SecantLineSearch.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/RegulaFalsiLineSearch.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/InitialInterpolatedLineSearch.o \
	$(FE)/convergenceTest/ConvergenceTest.o \
	$(FE)/convergenceTest/CTestNormUnbalance.o \
	$(FE)/convergenceTest/CTestNormDispIncr.o \
	$(FE)/convergenceTest/CTestEnergyIncr.o \
	$(FE)/convergenceTest/CTestRelativeNormUnbalance.o \
	$(FE)/convergenceTest/CTestRelativeNormDispIncr.o \
	$(FE)/convergenceTest/CTestRelativeEnergyIncr.o \
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
	$(FE)/analysis/integrator/TransientIntegrator.o \
	$(FE)/analysis/integrator/Newmark.o \
	$(FE)/analysis/integrator/HHT.o \
	$(FE)/analysis/integrator/Newmark1.o \
	$(FE)/analysis/integrator/HHT1.o \
	$(FE)/analysis/integrator/CentralDifference.o \
	$(FE)/analysis/integrator/WilsonTheta.o \
	$(FE)/analysis/model/AnalysisModel.o \
	$(FE)/analysis/model/simple/SimpleFE_Iter.o \
	$(FE)/analysis/model/simple/SimpleDOF_Iter.o \
	$(FE)/analysis/handler/ConstraintHandler.o \
	$(FE)/analysis/handler/PlainHandler.o \
	$(FE)/analysis/handler/PenaltyConstraintHandler.o \
	$(FE)/analysis/handler/LagrangeConstraintHandler.o \
	$(FE)/analysis/handler/TransformationConstraintHandler.o \
	$(FE)/analysis/numberer/DOF_Numberer.o \
	$(FE)/analysis/numberer/PlainNumberer.o \
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
	$(FE)/actor/message/Message.o \
	$(FE)/actor/machineBroker/AlphaMachineBroker.o \
	$(FE)/actor/machineBroker/MillMachineBroker.o \
	$(FE)/actor/machineBroker/DecMachineBroker.o \
	$(FE)/actor/objectBroker/FEM_ObjectBroker.o \
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
	@$(CD) $(FE)/nDarray; $(MAKE) wipe;
	@$(CD) $(FE)/convergenceTest; $(MAKE) wipe;
	@$(CD) $(FE)/database; $(MAKE) wipe;
	@$(CD) $(FE)/recorder; $(MAKE) wipe;
	@$(CD) $(FE)/renderer; $(MAKE) wipe;
	@$(CD) $(FE)/handler; $(MAKE) wipe;
	@$(CD) $(FE)/machine; $(MAKE) wipe;
	@$(CD) $(FE)/coordTransformation; $(MAKE) wipe;
	@$(CD) $(FE)/element; $(MAKE) wipe;
	@$(CD) $(FE)/reliability; $(MAKE) wipe;
	@$(CD) $(FE)/damage; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.






