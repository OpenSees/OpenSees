#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

# add_library(OPS_Tcl OBJECT)
# Uncomment once all tcl-related files are moved under the SRC/tcl
# directory

#target_sources(OPS_Tcl PRIVATE tclMain.cpp)

target_sources(OPS_InterpTcl PRIVATE 
    "${OPS_SRC_DIR}/coordTransformation/TclGeomTransfCommand.cpp"
    "${OPS_SRC_DIR}/database/TclDatabaseCommands.cpp"
    "${OPS_SRC_DIR}/database/TclMySQL.cpp"
    #"${OPS_SRC_DIR}/database/TclBerkeleyDB.cpp"
    "${OPS_SRC_DIR}/domain/component/TclParameterCommands.cpp"
    "${OPS_SRC_DIR}/domain/component/TclUpdateMaterialCommand.cpp"
    "${OPS_SRC_DIR}/domain/groundMotion/TclGroundMotionCommand.cpp"
    "${OPS_SRC_DIR}/domain/pattern/TclPatternCommand.cpp"
    "${OPS_SRC_DIR}/domain/pattern/TclSeriesCommand.cpp"
    "${OPS_SRC_DIR}/domain/pattern/TclSeriesIntegratorCommand.cpp"
    #"${OPS_SRC_DIR}/domain/pattern/drm/TclPatternCommand.cpp"
    "${OPS_SRC_DIR}/domain/region/TclRegionCommands.cpp"
    "${OPS_SRC_DIR}/element/TclElementCommands.cpp"
    "${OPS_SRC_DIR}/element/UP-ucsd/TclFourNodeQuadUPCommand.cpp"
    "${OPS_SRC_DIR}/element/UWelements/Tcl_generateInterfacePoints.cpp"
    "${OPS_SRC_DIR}/element/adapter/TclActuatorCommand.cpp"
    "${OPS_SRC_DIR}/element/adapter/TclActuatorCorotCommand.cpp"
    "${OPS_SRC_DIR}/element/adapter/TclAdapterCommand.cpp"
    "${OPS_SRC_DIR}/element/beamWithHinges/TclBeamWithHingesBuilder.cpp"
    "${OPS_SRC_DIR}/element/brick/TclBrickCommand.cpp"
    "${OPS_SRC_DIR}/element/brick/TclTwenty_Node_BrickCommand.cpp"
    "${OPS_SRC_DIR}/element/dispBeamColumnInt/TclDispBeamColumnIntCommand.cpp"
    "${OPS_SRC_DIR}/element/elasticBeamColumn/TclElasticBeamCommand.cpp"
    "${OPS_SRC_DIR}/element/elastomericBearing/TclElastomericBearingBoucWenCommand.cpp"
    "${OPS_SRC_DIR}/element/elastomericBearing/TclElastomericBearingPlasticityCommand.cpp"
    "${OPS_SRC_DIR}/element/elastomericBearing/TclElastomericBearingUFRPCommand.cpp"
    "${OPS_SRC_DIR}/element/forceBeamColumn/TclForceBeamColumnCommand.cpp"
    "${OPS_SRC_DIR}/element/fourNodeQuad/TclFourNodeQuadCommand.cpp"
    "${OPS_SRC_DIR}/element/frictionBearing/TclFlatSliderCommand.cpp"
    "${OPS_SRC_DIR}/element/frictionBearing/TclRJWatsonEQSCommand.cpp"
    "${OPS_SRC_DIR}/element/frictionBearing/TclSingleFPCommand.cpp"
    "${OPS_SRC_DIR}/element/frictionBearing/frictionModel/TclModelBuilderFrictionModelCommand.cpp"
    "${OPS_SRC_DIR}/element/generic/TclGenericClientCommand.cpp"
    "${OPS_SRC_DIR}/element/generic/TclGenericCopyCommand.cpp"
    "${OPS_SRC_DIR}/element/gradientInelasticBeamColumn/TclGradientInelasticBeamColumnCommand.cpp"
    "${OPS_SRC_DIR}/element/joint/TclBeamColumnJointCommand.cpp"
    "${OPS_SRC_DIR}/element/joint/TclJoint2dCommand.cpp"
    "${OPS_SRC_DIR}/element/joint/TclJoint3dCommand.cpp"
    "${OPS_SRC_DIR}/element/twoNodeLink/TclTwoNodeLinkCommand.cpp"
    "${OPS_SRC_DIR}/element/updatedLagrangianBeamColumn/TclCyclicModelCommands.cpp"
    "${OPS_SRC_DIR}/element/updatedLagrangianBeamColumn/TclElement2dGNL.cpp"
    "${OPS_SRC_DIR}/element/updatedLagrangianBeamColumn/TclElement2dYS.cpp"
    "${OPS_SRC_DIR}/element/zeroLength/TclZeroLength.cpp"
    #"${OPS_SRC_DIR}/element/PFEMElement/TclModelBuilder_addPFEMElement.cpp"
    #"${OPS_SRC_DIR}/element/PFEMElement/TclPFEMCommands.cpp"
    #"${OPS_SRC_DIR}/element/feap/TclFeapElementCommand.cpp"
    #"${OPS_SRC_DIR}/element/nonlinearBeamColumn/tcl/TclElmtBuilder.cpp"
    # Materials
    "${OPS_SRC_DIR}/material/nD/TclModelBuilderNDMaterialCommand.cpp"
    "${OPS_SRC_DIR}/material/nD/feap/TclFeapMaterialCommand.cpp"
    "${OPS_SRC_DIR}/material/nD/soil/TclUpdateMaterialStageCommand.cpp"
    "${OPS_SRC_DIR}/material/section/TclModelBuilderSectionCommand.cpp"
    "${OPS_SRC_DIR}/material/section/yieldSurface/TclModelBuilderYS_SectionCommand.cpp"
    "${OPS_SRC_DIR}/material/uniaxial/PY/TclPyTzQzMaterialCommand.cpp"
    "${OPS_SRC_DIR}/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp"
    "${OPS_SRC_DIR}/material/uniaxial/backbone/TclModelBuilderBackboneCommand.cpp"
    "${OPS_SRC_DIR}/material/uniaxial/drain/TclDrainMaterialCommand.cpp"
    "${OPS_SRC_DIR}/material/uniaxial/fedeas/TclFedeasMaterialCommand.cpp"
    "${OPS_SRC_DIR}/material/uniaxial/limitState/TclLimitState.cpp"
    "${OPS_SRC_DIR}/material/uniaxial/snap/TclSnapMaterialCommand.cpp"
    "${OPS_SRC_DIR}/material/uniaxial/stiffness/TclModelBuilderStiffnessDegradationCommand.cpp"
    "${OPS_SRC_DIR}/material/uniaxial/strength/TclModelBuilderStrengthDegradationCommand.cpp"
    "${OPS_SRC_DIR}/material/uniaxial/unloading/TclModelBuilderUnloadingRuleCommand.cpp"
    #"${OPS_SRC_DIR}/material/uniaxial/TclNewUniaxialMaterial.cpp"
    "${OPS_SRC_DIR}/material/yieldSurface/evolution/TclModelBuilderYS_EvolutionCommand.cpp"
    "${OPS_SRC_DIR}/material/yieldSurface/plasticHardeningMaterial/TclModelBuilderYSPlasticMaterialCommand.cpp"
    "${OPS_SRC_DIR}/material/yieldSurface/yieldSurfaceBC/TclModelBuilderYieldSurfaceBCCommand.cpp"
    # Model Builders
    "${OPS_SRC_DIR}/modelbuilder/tcl/TclModelBuilder.cpp"
    "${OPS_SRC_DIR}/modelbuilder/tcl/TclPlaneStressMaterialTester.cpp"
    "${OPS_SRC_DIR}/modelbuilder/tcl/TclSectionTester.cpp"
    "${OPS_SRC_DIR}/modelbuilder/tcl/TclUniaxialMaterialTester.cpp"
    "${OPS_SRC_DIR}/optimization/tcl/TclOptimizationBuilder.cpp"
    "${OPS_SRC_DIR}/recorder/TclRecorderCommands.cpp"
    "${OPS_SRC_DIR}/damage/TclModelBuilderDamageModelCommand.cpp"
    #"${OPS_SRC_DIR}/api/elementAPI_TCL.cpp"
    "${OPS_SRC_DIR}/tcl/commands.cpp"
  "${OPS_SRC_DIR}/tcl/TclFeViewer.cpp"
    
    #"${OPS_SRC_DIR}/interpreter/TclInterpreter.cpp"
    "${OPS_SRC_DIR}/interpreter/TclWrapper.cpp"

    #"${OPS_SRC_DIR}/reliability/analysis/gFunction/TclGFunEvaluator.cpp"
    #"${OPS_SRC_DIR}/reliability/analysis/gFunction/TclMatlabGFunEvaluator.cpp"
    #"${OPS_SRC_DIR}/reliability/domain/functionEvaluator/TclEvaluator.cpp"
    #"${OPS_SRC_DIR}/reliability/tcl/TclReliabilityBuilder.cpp"
    #"${OPS_SRC_DIR}/system_of_eqn/linearSOE/petsc/TclPetsc.cpp"
    #"${OPS_SRC_DIR}/tcl/TclVideoPlayer.cpp"
)


