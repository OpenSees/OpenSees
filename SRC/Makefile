# Makefile for fe objects

include ../Makefile.def

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
	@$(AR) $(ARFLAGS) $(FE_LIBRARY) $(MATRIX_LIBS) $(ELE_LIBS)\
	$(MATERIAL_LIBS) $(HYSTERETIC_LIBS) $(SECTION_LIBS) $(DATABASE_LIBS) \
	$(MODEL_BUILDER_LIBS) $(DOMAIN_LIBS)  \
	$(ANALYSIS_LIBS) $(SysOfEqn_LIBS) $(GRAPH_LIBS) $(TAGGED_LIBS) \
	$(UTILITY_LIBS) $(ACTOR_LIBS) $(RECORDER_LIBS) $(LAW_LIBS) \
	$(RENDERER_LIBS) $(REMO_LIBS) $(GLF_LIBS) $(ERRORHANDLER_LIBS)
	@$(RANLIB) $(FE_LIBRARY)

skypack:
	@$(CD) $(HOME)/lib; rm -f $(HOME)/lib/libskypack.a;     
	@$(CD) $(FE)/system_of_eqn/linearSOE/profileSPD/skypack; creator;

metis:
	@$(CD) $(HOME)/lib; rm -f $(METIS_LIBRARY);     
	@$(CD) $(FE)/graph/partitioner/metis-2.0; $(MAKE) clean;
	@$(CD) $(FE)/graph/partitioner/metis-2.0; $(MAKE) lib;

# the files that make up the FE_LIBRARY
GLF_LIBS = $(FE)/element/zeroLength/ZeroLength.o \
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
	$(FE)/coordTransformation/CrdTransf.o \
	$(FE)/coordTransformation/CrdTransf2d.o \
	$(FE)/coordTransformation/CrdTransf3d.o \
	$(FE)/coordTransformation/LinearCrdTransf2d.o \
	$(FE)/coordTransformation/LinearCrdTransf3d.o \
	$(FE)/element/nonlinearBeamColumn/tcl/TclElmtBuilder.o

LAW_LIBS = $(FE)/system_of_eqn/linearSOE/sparseSYM/SymSparseLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/SymSparseLinSolver.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/grcm.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/nest.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/nmat.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/symbolic.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/utility.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/newordr.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/nnsim.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/tim.o \
	$(FE)/system_of_eqn/linearSOE/sparseSYM/genmmd.o \
	$(FE)/system_of_eqn/eigenSOE/EigenSOE.o \
	$(FE)/system_of_eqn/eigenSOE/EigenSolver.o \
	$(FE)/system_of_eqn/eigenSOE/BandArpackSOE.o \
	$(FE)/system_of_eqn/eigenSOE/BandArpackSolver.o \
	$(FE)/system_of_eqn/eigenSOE/SymArpackSOE.o \
	$(FE)/system_of_eqn/eigenSOE/SymArpackSolver.o \
	$(FE)/analysis/analysis/EigenAnalysis.o \
	$(FE)/analysis/integrator/EigenIntegrator.o \
	$(FE)/analysis/algorithm/eigenAlgo/EigenAlgorithm.o \
	$(FE)/analysis/algorithm/eigenAlgo/FrequencyAlgo.o 


ERRORHANDLER_LIBS = $(FE)/handler/ErrorHandler.o \
	$(FE)/handler/ConsoleErrorHandler.o 

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


RECORDER_LIBS = $(FE)/recorder/MaxNodeDispRecorder.o \
	$(FE)/recorder/DatastoreRecorder.o \
	$(FE)/recorder/NodeRecorder.o \
	$(FE)/recorder/ElementRecorder.o \
	$(FE)/recorder/FilePlotter.o \
	$(FE)/recorder/TclRecorderCommands.o

DATABASE_LIBS = $(FE)/database/FileDatastore.o \
	$(FE)/database/FE_Datastore.o \
	$(FE)/database/TclDatabaseCommands.o

MATRIX_LIBS   = $(FE)/matrix/Matrix.o \
	$(FE)/matrix/Vector.o \
	$(FE)/matrix/ID.o  \
	$(FE)/matrix/nDarray.o \
	$(FE)/matrix/BJtensor.o \
	$(FE)/matrix/BJmatrix.o \
	$(FE)/matrix/BJvector.o \
	$(FE)/matrix/basics.o

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
	$(FE)/element/beam2d/ElasticBeam2d.o \
	$(FE)/element/beam2d/beam2d03.o \
	$(FE)/element/beam2d/beam2d04.o \
	$(FE)/element/beam3d/ElasticBeam3d.o \
	$(FE)/element/beam3d/beam3d01.o \
	$(FE)/element/beam3d/beam3d02.o \
	$(FE)/element/truss/Truss.o \
	$(FE)/element/truss/TrussSection.o \
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
	$(FE)/element/fourNodeQuad/TclFourNodeQuadCommand.o \
	$(FE)/element/fourNodeQuad/FourNodeQuad.o 


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
	$(FE)/material/nD/NDMaterial.o \
	$(FE)/material/nD/ElasticIsotropicMaterial.o \
	$(FE)/material/nD/ElasticIsotropicPlaneStress2D.o \
	$(FE)/material/nD/ElasticIsotropicPlaneStrain2D.o \
	$(FE)/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.o \
	$(FE)/material/nD/TclModelBuilderNDMaterialCommand.o 


SECTION_LIBS = $(FE)/material/section/SectionForceDeformation.o \
	$(FE)/material/section/ElasticSection2d.o \
	$(FE)/material/section/ElasticSection3d.o \
	$(FE)/material/section/FiberSection.o \
	$(FE)/material/section/SectionAggregator.o \
	$(FE)/material/section/GenericSection1d.o \
	$(FE)/material/section/GenericSectionNd.o \
	$(FE)/material/section/fiber/Fiber.o \
	$(FE)/material/section/fiber/UniaxialFiber.o \
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
	$(FE)/material/section/TclModelBuilderSectionCommand.o


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
	$(FE)/system_of_eqn/linearSOE/sparseGEN/SuperLU.o \
	$(FE)/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSOE.o \
	$(FE)/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.o 


MODEL_BUILDER_LIBS = $(FE)/modelbuilder/ModelBuilder.o \
	$(FE)/modelbuilder/PlaneFrame.o \
	$(FE)/modelbuilder/tcl/TclModelBuilder.o 


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
	$(FE)/domain/node/Node.o \
	$(FE)/domain/node/NodalLoad.o \
	$(FE)/domain/constraints/SP_Constraint.o \
	$(FE)/domain/constraints/MP_Constraint.o \
	$(FE)/domain/constraints/RigidDiaphragm.o \
	$(FE)/domain/constraints/RigidRod.o \
	$(FE)/domain/constraints/RigidBeam.o \
	$(FE)/domain/load/NodalLoadIter.o \
	$(FE)/domain/load/ElementalLoadIter.o \
	$(FE)/domain/load/Load.o \
	$(FE)/domain/pattern/LoadPattern.o \
	$(FE)/domain/pattern/EarthquakePattern.o \
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
	$(FE)/domain/groundMotion/GroundMotion.o \
	$(FE)/domain/groundMotion/GroundMotionIntegrator.o \
	$(FE)/domain/groundMotion/GroundMotionRecord.o \
	$(FE)/domain/groundMotion/TrapezoidalGroundMotionIntegrator.o \
	$(FE)/domain/subdomain/Subdomain.o \
	$(FE)/domain/subdomain/SubdomainNodIter.o \
	$(FE)/analysis/analysis/DomainUser.o 


ANALYSIS_LIBS = $(FE)/analysis/analysis/Analysis.o \
	$(FE)/analysis/analysis/StaticAnalysis.o \
	$(FE)/analysis/analysis/TransientAnalysis.o \
	$(FE)/analysis/analysis/DirectIntegrationAnalysis.o \
	$(FE)/analysis/analysis/DomainDecompositionAnalysis.o \
	$(FE)/analysis/analysis/SubstructuringAnalysis.o \
	$(FE)/analysis/algorithm/SolutionAlgorithm.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/EquiSolnAlgo.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/Linear.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/NewtonRaphson.o \
	$(FE)/analysis/algorithm/equiSolnAlgo/ModifiedNewton.o \
	$(FE)/convergenceTest/ConvergenceTest.o \
	$(FE)/convergenceTest/CTestNormUnbalance.o \
	$(FE)/convergenceTest/CTestNormDispIncr.o \
	$(FE)/convergenceTest/CTestEnergyIncr.o \
	$(FE)/analysis/algorithm/domainDecompAlgo/DomainDecompAlgo.o \
	$(FE)/analysis/integrator/Integrator.o \
	$(FE)/analysis/integrator/IncrementalIntegrator.o \
	$(FE)/analysis/integrator/StaticIntegrator.o \
	$(FE)/analysis/integrator/LoadControl.o \
	$(FE)/analysis/integrator/LoadPath.o \
	$(FE)/analysis/integrator/ArcLength.o \
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
	@$(CD) $(FE)/convergenceTest; $(MAKE) wipe;
	@$(CD) $(FE)/database; $(MAKE) wipe;
	@$(CD) $(FE)/recorder; $(MAKE) wipe;
	@$(CD) $(FE)/renderer; $(MAKE) wipe;
	@$(CD) $(FE)/handler; $(MAKE) wipe;
	@$(CD) $(FE)/machine; $(MAKE) wipe;
	@$(CD) $(FE)/coordTransformation; $(MAKE) wipe;
	@$(CD) $(FE)/element; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.





