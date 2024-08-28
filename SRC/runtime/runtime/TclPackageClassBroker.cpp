//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Purpose: This file contains the class definition for TclPackageClassBroker.
// TclPackageClassBroker is is an object broker class that is meant to become
// a threadsafe replacement for the BrokerAllClasses class.
// All methods are virtual to allow for subclasses; which can be
// used by programmers when introducing new subclasses of the main objects.
//
#ifdef _PARALLEL_PROCESSING
#  include <mpi.h>
#endif

#ifdef _PARALLEL_INTERPRETERS
#  include <mpi.h>
#endif

#include <Hash.h>
using namespace OpenSees::Hash;
using namespace OpenSees::Hash::literals;
#define DISPATCH(symbol) case  hasher<std::string>()(#symbol): return new symbol();
//
// case hasher<std::string>()(Truss::class_name):  return new Truss();
//
#include "packages.h"
#include <TclPackageClassBroker.h>

// ActorTypes
#include "domain/subdomain/ActorSubdomain.h"

// Convergence tests
#include "analysis/criteria/CTestNormUnbalance.h"
#include "analysis/criteria/CTestRelativeNormUnbalance.h"
#include "analysis/criteria/CTestNormDispIncr.h"
#include "analysis/criteria/CTestRelativeNormDispIncr.h"
#include "analysis/criteria/CTestRelativeTotalNormDispIncr.h"
#include "analysis/criteria/CTestEnergyIncr.h"
#include "analysis/criteria/CTestRelativeEnergyIncr.h"
#include "analysis/criteria/CTestFixedNumIter.h"

// graph numbering schemes
#include "graph/numberer/RCM.h"
#include "graph/numberer/SimpleNumberer.h"

// uniaxial material model header files
#include "BoucWen/BoucWenMaterial.h"
#include "SPSW02.h"
#include "ElasticMaterial.h"
#include "ElasticMultiLinear.h"
#include "ElasticPowerFunc.h"
#include "Elastic2Material.h"
#include "ElasticPPMaterial.h"
#include "ParallelMaterial.h"
#include "ASD_SMA_3K.h"
#include "Concrete01.h"
#include "Concrete02.h"
#include "Concrete04.h"
#include "Concrete06.h"
#include "Concrete07.h"
#include "ConcretewBeta.h"
#include "OriginCentered.h"
#include "Steel01.h"
#include "Steel02.h"
#include "Steel2.h"
#include "Steel4.h"
#include "FatigueMaterial.h"
#include "ReinforcingSteel/ReinforcingSteel.h"
#include "HardeningMaterial.h"
#include "HystereticMaterial.h"
#include "EPPGapMaterial.h"
#include "ViscousMaterial.h"
#include "ViscousDamper.h"
#include "PathIndependentMaterial.h"
#include "SeriesMaterial.h"
#include "CableMaterial.h"
#include "ENTMaterial.h"
#include "MinMaxMaterial.h"
#include "ModIMKPeakOriented.h"
#include "snap/Clough.h"
#include "limitState/LimitStateMaterial.h"
#include "InitStressMaterial.h"
#include "InitStrainMaterial.h"
#include "Bond_SP01.h"
#include "SimpleFractureMaterial.h"
#include "ConfinedConcrete01.h"
#include <HystereticPoly.h> // Salvatore Sessa 14-Jan-2021

// PY springs: RWBoulanger and BJeremic
#include "PY/PySimple1.h"
#include "PY/TzSimple1.h"
#include "PY/QzSimple1.h"
#include "PY/PySimple2.h"
#include "PY/TzSimple2.h"
#include "PY/QzSimple2.h"
#include "PY/PyLiq1.h"
#include "PY/TzLiq1.h"
#include "PY/QzLiq1.h"

#include "fedeas/FedeasBond1Material.h"
#include "fedeas/FedeasBond2Material.h"
#include "fedeas/FedeasConcr1Material.h"
#include "fedeas/FedeasConcr2Material.h"
#include "fedeas/FedeasConcr3Material.h"
#include "fedeas/FedeasHardeningMaterial.h"
#include "fedeas/FedeasHyster1Material.h"
#include "fedeas/FedeasHyster2Material.h"
#include "fedeas/FedeasSteel1Material.h"
#include "fedeas/FedeasSteel2Material.h"

#include "Bilin.h"
#include "drain/DrainBilinearMaterial.h"
#include "drain/DrainClough1Material.h"
#include "drain/DrainClough2Material.h"
#include "drain/DrainPinch1Material.h"
#include "HyperbolicGapMaterial.h"
#include "ImpactMaterial.h"

// Sections
#include "ElasticSection2d.h"
#include "ElasticSection3d.h"
#include "ElasticShearSection2d.h"
#include "ElasticShearSection3d.h"
#include "GenericSection1d.h"
//#include "GenericSectionNd.h"
#include "SectionAggregator.h"
//#include "FiberSection.h"
#include "FiberSection2d.h"
#include "FiberSection3d.h"
#include "FiberSectionAsym3d.h" //Xinlong Du
#include "ElasticPlateSection.h"
#include "ElasticMembranePlateSection.h"
#include "MembranePlateFiberSection.h"
#include "Bidirectional.h"
#include "LayeredShellFiberSection.h" // Yuli Huang & Xinzheng Lu

// NDMaterials
#include "ElasticIsotropicPlaneStrain2D.h"
#include "ElasticIsotropicPlaneStress2D.h"
#include "ElasticIsotropicPlateFiber.h"
#include "ElasticIsotropicAxiSymm.h"
#include "ElasticIsotropicThreeDimensional.h"
#include "J2PlaneStrain.h"
#include "J2PlaneStress.h"
#include "J2PlateFiber.h"
#include "J2AxiSymm.h"
#include "J2ThreeDimensional.h"
#include "PlaneStressMaterial.h"
#include "PlateFiberMaterial.h"
// start Yuli Huang & Xinzheng L
#include "PlateRebarMaterial.h"
#include "PlateFromPlaneStressMaterial.h"
//#include "ConcreteS.h"
#include "PlaneStressUserMaterial.h"
// end Yuli Huang & Xinzheng Lu
#include "feap/FeapMaterial03.h"
#include "CycLiqCP3D.h"
#include "CycLiqCPPlaneStrain.h"
#include "CycLiqCPSP3D.h"
#include "CycLiqCPSPPlaneStrain.h"

#include "soil/FluidSolidPorousMaterial.h"
#include "soil/PressureDependMultiYield.h"
#include "soil/PressureDependMultiYield02.h"
#include "soil/PressureIndependMultiYield.h"

#include "UWmaterials/ContactMaterial2D.h"
#include "UWmaterials/ContactMaterial3D.h"
#include "UWmaterials/DruckerPrager3D.h"
#include "UWmaterials/DruckerPragerPlaneStrain.h"
#include "UWmaterials/BoundingCamClay.h"
#include "UWmaterials/BoundingCamClay3D.h"
#include "UWmaterials/BoundingCamClayPlaneStrain.h"
#include "UWmaterials/ManzariDafalias.h"
#include "UWmaterials/ManzariDafalias3D.h"
#include "UWmaterials/ManzariDafaliasPlaneStrain.h"
#include "UWmaterials/ManzariDafaliasRO.h"
#include "UWmaterials/ManzariDafalias3DRO.h"
#include "UWmaterials/ManzariDafaliasPlaneStrainRO.h"
#include "UWmaterials/PM4Sand.h"
#include "UWmaterials/PM4Silt.h"
#include "J2CyclicBoundingSurface.h"
#include "J2CyclicBoundingSurface3D.h"
#include "J2CyclicBoundingSurfacePlaneStrain.h"
#include "UWmaterials/InitialStateAnalysisWrapper.h"
#include "stressDensityModel/stressDensity.h"
#include "InitStressNDMaterial.h"

//
// element header files
//
#include "Element.h"
#include "Truss/Truss.h"
#include "Truss/Truss2.h"
#include "Truss/TrussSection.h"
#include "Truss/CorotTruss.h"
#include "Truss/CorotTrussSection.h"
#include "Truss/InertiaTruss.h"
#include "Point/ZeroLength.h"
#include "Point/ZeroLengthSection.h"
#include "Point/ZeroLengthContact2D.h"
#include "Point/ZeroLengthContact3D.h"
#include "Point/ZeroLengthContactNTS2D.h"
#include "Point/ZeroLengthInterface2D.h"
#include "Point/ZeroLengthContactASDimplex.h"
//#include "ZeroLengthND.h"

#include "Plane/FourNodeQuad.h"
#include "Plane/EnhancedQuad.h"
#include "Plane/NineNodeMixedQuad.h"
#include "Plane/NineNodeQuad.h"
#include "Plane/EightNodeQuad.h"
#include "Plane/ConstantPressureVolumeQuad.h"
#include "Plane/Nine_Four_Node_QuadUP.h"
#include "Plane/BBarFourNodeQuadUP.h"
#include "Plane/FourNodeQuadUP.h"
#include "Plane/Tri31.h"

#include "Frame/Elastic/ElasticBeam2d.h"
#include "Frame/Elastic/ElasticBeam3d.h"
#include "Frame/Elastic/ModElasticBeam2d.h" //SAJalali
#include "Frame/Elastic/ModElasticBeam3d.h"
#include "Frame/Elastic/ElasticTimoshenkoBeam2d.h"
#include "Frame/Elastic/ElasticTimoshenkoBeam3d.h"
#include "Frame/Other/Force/ForceBeamColumn2d.h"
#include "Frame/Other/Force/ForceBeamColumn3d.h"
#include "Frame/Other/Displ/DispBeamColumn2d.h"
#include "Frame/Other/Displ/DispBeamColumn3d.h"
#include "Frame/Other/Displ/DispBeamColumnAsym3d.h"   // Xinlong Du
#include "Frame/Other/Mixed/MixedBeamColumnAsym3d.h"  // Xinlong Du


#include "UWelements/SSPquad.h"
#include "UWelements/SSPquadUP.h"
#include "UWelements/SSPbrick.h"
#include "UWelements/SSPbrickUP.h"
#include "UWelements/BeamContact2D.h"
#include "UWelements/BeamContact2Dp.h"
#include "UWelements/BeamContact3D.h"
#include "UWelements/BeamContact3Dp.h"
#include "UWelements/BeamEndContact3D.h"
#include "UWelements/BeamEndContact3Dp.h"
#include "UWelements/QuadBeamEmbedContact.h"

#include "Other/PML/PML2D.h"
#include "Other/PML/PML3D.h"

#include "Shell/ShellMITC4.h"
#include "Shell/ShellMITC9.h"
#include "Shell/ShellDKGQ.h"   // Added by Lisha Wang, Xinzheng Lu, Linlin Xie, Song Cen & Quan Gu
#include "Shell/ShellNLDKGQ.h" // Added by Lisha Wang, Xinzheng Lu, Linlin Xie, Song Cen & Quan Gu
#include "Shell/ASDShellQ4.h"  // Massimo Petracca
#include "Brick/Brick.h"
#include "Brick/BbarBrick.h"
#include "Brick/BrickUP.h"
#include "Brick/BBarBrickUP.h"
#include "Brick/Twenty_Eight_Node_BrickUP.h"
#include "Joint/Joint2D.h" // Arash
#include "Link/TwoNodeLink.h"
#include "Link/LinearElasticSpring.h"
#include "Link/Inerter.h"

// Bearings
#include "Bearing/elastomeric/ElastomericBearingBoucWen2d.h"
#include "Bearing/elastomeric/ElastomericBearingBoucWen3d.h"
#include "Bearing/elastomeric/ElastomericBearingPlasticity2d.h"
#include "Bearing/elastomeric/ElastomericBearingPlasticity3d.h"
#include "Bearing/elastomeric/ElastomericBearingUFRP2d.h"
#include "Bearing/elastomeric/ElastomericX.h"
#include "Bearing/elastomeric/HDR.h"
#include "Bearing/elastomeric/LeadRubberX.h"
#include "Bearing/friction/FlatSliderSimple2d.h"
#include "Bearing/friction/FlatSliderSimple3d.h"
#include "Bearing/friction/FPBearingPTV.h"
#include "Bearing/friction/RJWatsonEQS2d.h"
#include "Bearing/friction/RJWatsonEQS3d.h"
#include "Bearing/friction/SingleFPSimple2d.h"
#include "Bearing/friction/SingleFPSimple3d.h"
#include "Bearing/friction/TripleFrictionPendulum.h"
// friction models
#include "Bearing/friction/frictionModel/Coulomb.h"
#include "Bearing/friction/frictionModel/VelDependent.h"
#include "Bearing/friction/frictionModel/VelPressureDep.h"
#include "Bearing/friction/frictionModel/VelDepMultiLinear.h"
#include "Bearing/friction/frictionModel/VelNormalFrcDep.h"


#include "mvlem/MVLEM.h"        // Kristijan Kolozvari
#include "mvlem/SFI_MVLEM.h"    // Kristijan Kolozvari
#include "mvlem/MVLEM_3D.h"     // Kristijan Kolozvari
#include "mvlem/SFI_MVLEM_3D.h" // Kristijan Kolozvari

#include "Boundary/RockingBC.h"

#include "ASDEA/CEqElement/ASDEmbeddedNodeElement.h"
#include "ASDEA/absorbentBoundaries/ASDAbsorbingBoundary2D.h"
#include "ASDEA/absorbentBoundaries/ASDAbsorbingBoundary3D.h"

#include "LinearCrdTransf2d.h"
#include "LinearCrdTransf3d.h"
#include "PDeltaCrdTransf2d.h"
#include "PDeltaCrdTransf3d.h"
#include "CorotCrdTransf2d.h"
#include "CorotCrdTransf3d.h"

#include "quadrature/Frame/HingeMidpointBeamIntegration.h"
#include "quadrature/Frame/HingeEndpointBeamIntegration.h"
#include "quadrature/Frame/HingeRadauBeamIntegration.h"
#include "quadrature/Frame/HingeRadauTwoBeamIntegration.h"
#include "quadrature/Frame/UserDefinedHingeIntegration.h"
#include "quadrature/Frame/DistHingeIntegration.h"
#include "quadrature/Frame/RegularizedHingeIntegration.h"

#include "quadrature/Frame/LobattoBeamIntegration.h"
#include "quadrature/Frame/LegendreBeamIntegration.h"
#include "quadrature/Frame/RadauBeamIntegration.h"
#include "quadrature/Frame/NewtonCotesBeamIntegration.h"
#include "quadrature/Frame/TrapezoidalBeamIntegration.h"
#include "quadrature/Frame/UserDefinedBeamIntegration.h"
#include "quadrature/Frame/FixedLocationBeamIntegration.h"
#include "quadrature/Frame/LowOrderBeamIntegration.h"
#include "quadrature/Frame/MidDistanceBeamIntegration.h"
#include "quadrature/Frame/CompositeSimpsonBeamIntegration.h"

// node header files
#include "Node.h"
#ifdef HEAP_NODE
#include "HeapNode.h"
#endif

#include "FileStream.h"
#include "StandardStream.h"
#include "XmlFileStream.h"
#include "DataFileStream.h"
#include "DataFileStreamAdd.h"
#include "BinaryFileStream.h"
#include "DatabaseStream.h"
#include "DummyStream.h"

#include "NodeRecorder.h"
#include "ElementRecorder.h"
#include "EnvelopeNodeRecorder.h"
#include "EnvelopeElementRecorder.h"
#include "DriftRecorder.h"
//#include "MPCORecorder.h"
#include "VTK_Recorder.h"
#include "GmshRecorder.h"

// mp_constraint header files
#include "MP_Constraint.h"
#include "Joint/MP_Joint2D.h"

// sp_constraint header files
#include "SP_Constraint.h"
#include "SP_Constraint.h"
#include "ImposedMotionSP.h"
#include "ImposedMotionSP1.h"

// Pressure_Constraint header file
#include "Pressure_Constraint.h"

// nodal load header files
#include "NodalLoad.h"

// elemental load header files
#include "ElementalLoad.h"
#include "Beam2dUniformLoad.h"
#include "Beam2dPointLoad.h"
#include "Beam3dUniformLoad.h"
#include "Beam3dPointLoad.h"
#include "BrickSelfWeight.h"
#include "SelfWeight.h"
#include "SurfaceLoader.h"

// matrix, vector & id header files
#include "Matrix.h"
#include "Vector.h"
#include "ID.h"

// subdomain header files
#include "Subdomain.h"

// constraint handler header files
#include "ConstraintHandler.h"
#include "PlainHandler.h"
#include "PenaltyConstraintHandler.h"
#include "LagrangeConstraintHandler.h"
#include "TransformationConstraintHandler.h"

// dof numberer header files
#include "DOF_Numberer.h"
#include "PlainNumberer.h"

// analysis model header files
#include "AnalysisModel.h"

// equi soln algo header files
#include "EquiSolnAlgo.h"
#include "Linear.h"
#include "NewtonRaphson.h"
#include "Broyden.h"
#include "NewtonLineSearch.h"
#include "KrylovNewton.h"
#include "AcceleratedNewton.h"
#include "ModifiedNewton.h"

#include "accelerator/KrylovAccelerator.h"
#include "accelerator/RaphsonAccelerator.h"

#include "BisectionLineSearch.h"
#include "InitialInterpolatedLineSearch.h"
#include "RegulaFalsiLineSearch.h"
#include "SecantLineSearch.h"

// domain decomp soln algo header files
#include "DomainDecompAlgo.h"

// integrator header files
#include "ArcLength.h"
#include "DisplacementControl.h"
#ifdef _PARALLEL_PROCESSING
#include "DistributedDisplacementControl.h"
#endif
#include "LoadControl.h"
// #include "StagedLoadControl.h"

#include "TransientIntegrator.h"
#include "AlphaOS.h"
#include "AlphaOS_TP.h"
#include "AlphaOSGeneralized.h"
#include "AlphaOSGeneralized_TP.h"
#include "CentralDifference.h"
#include "CentralDifferenceAlternative.h"
#include "CentralDifferenceNoDamping.h"
#include "Collocation.h"
#include "CollocationHSFixedNumIter.h"
#include "CollocationHSIncrLimit.h"
#include "CollocationHSIncrReduct.h"
#include "HHT.h"
#include "HHT_TP.h"
#include "HHTExplicit.h"
#include "HHTExplicit_TP.h"
#include "HHTGeneralized.h"
#include "HHTGeneralized_TP.h"
#include "HHTGeneralizedExplicit.h"
#include "HHTGeneralizedExplicit_TP.h"
#include "HHTHSFixedNumIter.h"
#include "HHTHSFixedNumIter_TP.h"
#include "HHTHSIncrLimit.h"
#include "HHTHSIncrLimit_TP.h"
#include "HHTHSIncrReduct.h"
#include "HHTHSIncrReduct_TP.h"
#include "KRAlphaExplicit.h"
#include "KRAlphaExplicit_TP.h"
#include "Newmark.h"
// #include "StagedNewmark.h"
#include "NewmarkExplicit.h"
#include "NewmarkHSFixedNumIter.h"
#include "NewmarkHSIncrLimit.h"
#include "NewmarkHSIncrReduct.h"
#include "TRBDF2.h"
#include "TRBDF3.h"
#include "WilsonTheta.h"

// system of eqn header files
#include "LinearSOE.h"
#include "DomainSolver.h"
#include "fullGEN/FullGenLinSOE.h"
#include "bandGEN/BandGenLinSOE.h"
#include "bandSPD/BandSPDLinSOE.h"
#include "profileSPD/ProfileSPDLinSOE.h"
#include "profileSPD/ProfileSPDLinSubstrSolver.h"
#include "sparseGEN/SparseGenColLinSOE.h"
#include "DomainDecompositionAnalysis.h"

// load patterns
#include "LoadPattern.h"
#include "UniformExcitation.h"
#include "MultiSupportPattern.h"
#include "GroundMotion.h"
#include "InterpolatedGroundMotion.h"
#ifdef OPSDEF_DRM
#  include "drm/DRMLoadPatternWrapper.h"
#endif // OPSDEF_DRM

#ifdef _H5DRM
#  include "drm/H5DRM.h"
#endif

#include "Parameter.h"
#include "ElementParameter.h"
#include "MaterialStageParameter.h"
#include "MatParameter.h"
#include "InitialStateParameter.h"
#include "ElementStateParameter.h"

// time series
#include "LinearSeries.h"
#include "PathSeries.h"
#include "PathTimeSeries.h"
#include "RectangularSeries.h"
#include "ConstantSeries.h"
#include "TrigSeries.h"
#include "TriangleSeries.h"

// time series integrators
#include "TrapezoidalTimeSeriesIntegrator.h"

#include "eigenSOE/ArpackSOE.h"

#ifdef _PETSC
#  include "PetscSOE.h"
#  include "PetscSolver.h"
#  include "SparseGenColLinSOE.h"
#endif

#ifdef _MUMPS
#  include "MumpsSOE.h"
#  ifdef _PARALLEL_PROCESSING
#  include "MumpsParallelSOE.h"
#  endif
#endif

#ifdef _PARALLEL_PROCESSING
#  include "DistributedBandSPDLinSOE.h"
#  include "DistributedProfileSPDLinSOE.h"
#  include "DistributedSparseGenColLinSOE.h"
#  include "DistributedSparseGenRowLinSOE.h"
#  include "DistributedBandGenLinSOE.h"
#  include "DistributedSuperLU.h"
#  include "ParallelNumberer.h"
#  include "StaticDomainDecompositionAnalysis.h"
#  include "TransientDomainDecompositionAnalysis.h"
#  include "DistributedDiagonalSOE.h"
#endif


typedef struct uniaxialPackage {
  int classTag;
  char *libName;
  char *funcName;
  UniaxialMaterial *(*funcPtr)(void);
  struct uniaxialPackage *next;
} UniaxialPackage;

static UniaxialPackage *theUniaxialPackage = NULL;

TclPackageClassBroker::TclPackageClassBroker() : lastDomainSolver(0) {}

TclPackageClassBroker::~TclPackageClassBroker() {}

Actor *
TclPackageClassBroker::getNewActor(int classTag, Channel *theChannel)
{
  switch (classTag) {

#ifdef _PARALLEL_PROCESSING
  case ACTOR_TAGS_SUBDOMAIN:
    return new ActorSubdomain(*theChannel, *this);
#endif

  default:
    opserr << "TclPackageClassBroker::getNewActor - ";
    opserr << " - no ActorType type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

PartitionedModelBuilder *
TclPackageClassBroker::getPtrNewPartitionedModelBuilder(Subdomain &theSubdomain,
                                                        int classTag)
{
  switch (classTag) {
    /*
    case PartitionedModelBuilder_TAGS_PartitionedQuick2dFrameModel:
         return new PartitionedQuick2dFrame(theSubdomain);
         */

  default:
    opserr << "TclPackageClassBroker::getPtrNewPartitionedModelBuilder - ";
    opserr << " - no PartitionedModelBuilder type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

GraphNumberer *
TclPackageClassBroker::getPtrNewGraphNumberer(int classTag)
{
  switch (classTag) {
  case GraphNUMBERER_TAG_RCM:
    return new RCM();

  case GraphNUMBERER_TAG_SimpleNumberer:
    return new SimpleNumberer();

  default:
    opserr << "TclPackageClassBroker::getPtrNewGraphNumberer - ";
    opserr << " - no GraphNumberer type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

/*****************************************
 *
 * METHODS TO GET NEW MODELLING CLASSES
 *
 *****************************************/

Element *
TclPackageClassBroker::getNewElement(int classTag)
{
  switch ((std::size_t)classTag) {

    DISPATCH(Truss);
    DISPATCH(Truss2);
    DISPATCH(TrussSection);
    DISPATCH(CorotTruss);
    DISPATCH(CorotTrussSection);
    DISPATCH(InertiaTruss);

    // case ELE_TAG_ZeroLengthND:
    // return new ZeroLengthND();

    DISPATCH(FourNodeQuadUP);
    DISPATCH(FourNodeQuad);
    DISPATCH(Tri31);
    DISPATCH(ElasticBeam2d);
    DISPATCH(ModElasticBeam2d);
    DISPATCH(ModElasticBeam3d);
    DISPATCH(ElasticBeam3d);
    DISPATCH(ElasticTimoshenkoBeam2d);
    DISPATCH(ElasticTimoshenkoBeam3d);
    DISPATCH(ForceBeamColumn2d);
    DISPATCH(ForceBeamColumn3d);
    DISPATCH(DispBeamColumn2d);
    DISPATCH(DispBeamColumn3d);
    DISPATCH(DispBeamColumnAsym3d);
    DISPATCH(MixedBeamColumnAsym3d);

// Quads
    DISPATCH(EnhancedQuad);
    DISPATCH(NineNodeMixedQuad);
    DISPATCH(NineNodeQuad);
    DISPATCH(EightNodeQuad);
    DISPATCH(ConstantPressureVolumeQuad);
    DISPATCH(BBarFourNodeQuadUP);
    DISPATCH(NineFourNodeQuadUP);

#if defined(OPSDEF_Elements_UW)
    DISPATCH(SSPquad);

  case ELE_TAG_SSPquadUP:
    return new SSPquadUP;

    DISPATCH(SSPbrick);
    DISPATCH(SSPbrickUP);
#endif
    DISPATCH(PML2D);

    DISPATCH(PML3D);

// Bricks
    DISPATCH(Brick);
    DISPATCH(BbarBrick);
    DISPATCH(BBarBrickUP);
    DISPATCH(BrickUP);
    DISPATCH(TwentyEightNodeBrickUP);

// Shells
    DISPATCH(ShellMITC4);
    DISPATCH(ShellMITC9);
    DISPATCH(ShellDKGQ);
    DISPATCH(ShellNLDKGQ);
    DISPATCH(ASDShellQ4); // Massimo Petracca

#if defined(OPSDEF_Elements_UW)
    DISPATCH(BeamContact2D);
    DISPATCH(BeamContact2Dp);
    DISPATCH(BeamContact3D);
    DISPATCH(BeamContact3Dp);
    DISPATCH(BeamEndContact3D);
    DISPATCH(BeamEndContact3Dp);
    DISPATCH(QuadBeamEmbedContact);
#endif // UW elements

    DISPATCH(ZeroLength);
    DISPATCH(ZeroLengthSection);
    DISPATCH(ZeroLengthContact2D);
    DISPATCH(ZeroLengthContact3D);
    DISPATCH(ZeroLengthInterface2D);
    DISPATCH(ZeroLengthContactNTS2D);
    DISPATCH(ZeroLengthContactASDimplex);



    DISPATCH(Joint2D); // Arash
    DISPATCH(TwoNodeLink);
    DISPATCH(LinearElasticSpring);
    DISPATCH(Inerter);

    DISPATCH(MVLEM); // Kristijan Kolozvari
    DISPATCH(SFI_MVLEM); // Kristijan Kolozvari
    DISPATCH(MVLEM_3D); // Kristijan Kolozvari
    DISPATCH(SFI_MVLEM_3D); // Kristijan Kolozvari


    DISPATCH(ElastomericBearingBoucWen2d);
    DISPATCH(ElastomericBearingBoucWen3d);
    DISPATCH(ElastomericBearingPlasticity2d);
    DISPATCH(ElastomericBearingPlasticity3d);
    DISPATCH(ElastomericBearingUFRP2d);
    DISPATCH(ElastomericX);
    DISPATCH(HDR);
    DISPATCH(LeadRubberX);
    DISPATCH(FlatSliderSimple2d);
    DISPATCH(FlatSliderSimple3d);
    DISPATCH(FPBearingPTV);
    DISPATCH(RJWatsonEQS2d);
    DISPATCH(RJWatsonEQS3d);
    DISPATCH(SingleFPSimple2d);
    DISPATCH(SingleFPSimple3d);
    DISPATCH(TripleFrictionPendulum);


    DISPATCH(RockingBC);
    DISPATCH(ASDEmbeddedNodeElement);
    DISPATCH(ASDAbsorbingBoundary2D);
    DISPATCH(ASDAbsorbingBoundary3D);

  default:
    opserr << "TclPackageClassBroker::getNewElement - ";
    opserr << " - no Element type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

Node *
TclPackageClassBroker::getNewNode(int classTag)
{
  switch (classTag) {
  case NOD_TAG_Node:
#ifdef HEAP_NODE
    return new HeapNode(classTag);
#else
    return new Node(classTag);
#endif

  default:
    opserr << "TclPackageClassBroker::getNewNode - ";
    opserr << " - no Node type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

MP_Constraint *
TclPackageClassBroker::getNewMP(int classTag)
{
  switch (classTag) {
  case CNSTRNT_TAG_MP_Constraint:
    return new MP_Constraint(classTag);

  case CNSTRNT_TAG_MP_Joint2D: // Arash
    return new MP_Joint2D();   // Arash

  default:
    opserr << "TclPackageClassBroker::getNewMP - ";
    opserr << " - no MP_Constraint type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

SP_Constraint *
TclPackageClassBroker::getNewSP(int classTag)
{
  switch (classTag) {
  case CNSTRNT_TAG_SP_Constraint:
    return new SP_Constraint(classTag);

  case CNSTRNT_TAG_ImposedMotionSP:
    return new ImposedMotionSP();

  case CNSTRNT_TAG_ImposedMotionSP1:
    return new ImposedMotionSP1();

  default:
    opserr << "TclPackageClassBroker::getNewSP - ";
    opserr << " - no SP_Constraint type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

Pressure_Constraint *
TclPackageClassBroker::getNewPC(int classTag)
{
  switch (classTag) {
  case CNSTRNT_TAG_Pressure_Constraint:
    return new Pressure_Constraint(classTag);

  default:
    opserr << "TclPackageClassBroker::getNewPC - ";
    opserr << " - no Pressure_Constraint type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}


ElementalLoad *
TclPackageClassBroker::getNewElementalLoad(int classTag)
{
  switch (classTag) {

  case LOAD_TAG_Beam2dUniformLoad:
    return new Beam2dUniformLoad();

  case LOAD_TAG_Beam2dPointLoad:
    return new Beam2dPointLoad();

  case LOAD_TAG_Beam3dUniformLoad:
    return new Beam3dUniformLoad();

  case LOAD_TAG_Beam3dPointLoad:
    return new Beam3dPointLoad();

  case LOAD_TAG_BrickSelfWeight:
    return new BrickSelfWeight();

  case LOAD_TAG_SelfWeight:
    return new SelfWeight();

  case LOAD_TAG_SurfaceLoader:
    return new SurfaceLoader();

  default:
    opserr << "TclPackageClassBroker::getNewNodalLoad - ";
    opserr << " - no NodalLoad type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }

  return 0;
}

CrdTransf *
TclPackageClassBroker::getNewCrdTransf(int classTag)
{
  switch (classTag) {
  case CRDTR_TAG_LinearCrdTransf2d:
    return new LinearCrdTransf2d();
  case CRDTR_TAG_PDeltaCrdTransf2d:
    return new PDeltaCrdTransf2d();
  case CRDTR_TAG_CorotCrdTransf2d:
    return new CorotCrdTransf2d();
  case CRDTR_TAG_LinearCrdTransf3d:
    return new LinearCrdTransf3d();
  case CRDTR_TAG_PDeltaCrdTransf3d:
    return new PDeltaCrdTransf3d();
  case CRDTR_TAG_CorotCrdTransf3d:
    return new CorotCrdTransf3d();
  default:
    opserr << "TclPackageClassBroker::getCrdTransf - ";
    opserr << " - no CrdTransf type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

BeamIntegration *
TclPackageClassBroker::getNewBeamIntegration(int classTag)
{
  switch (classTag) {
  case BEAM_INTEGRATION_TAG_Lobatto:
    return new LobattoBeamIntegration();

  case BEAM_INTEGRATION_TAG_Legendre:
    return new LegendreBeamIntegration();

  case BEAM_INTEGRATION_TAG_Radau:
    return new RadauBeamIntegration();

  case BEAM_INTEGRATION_TAG_NewtonCotes:
    return new NewtonCotesBeamIntegration();

  case BEAM_INTEGRATION_TAG_Trapezoidal:
    return new TrapezoidalBeamIntegration();

  case BEAM_INTEGRATION_TAG_UserDefined:
    return new UserDefinedBeamIntegration();

  case BEAM_INTEGRATION_TAG_FixedLocation:
    return new FixedLocationBeamIntegration();

  case BEAM_INTEGRATION_TAG_LowOrder:
    return new LowOrderBeamIntegration();

  case BEAM_INTEGRATION_TAG_MidDistance:
    return new MidDistanceBeamIntegration();

  case BEAM_INTEGRATION_TAG_CompositeSimpson:
    return new CompositeSimpsonBeamIntegration();

  case BEAM_INTEGRATION_TAG_HingeMidpoint:
    return new HingeMidpointBeamIntegration();

  case BEAM_INTEGRATION_TAG_HingeRadau:
    return new HingeRadauBeamIntegration();

  case BEAM_INTEGRATION_TAG_HingeRadauTwo:
    return new HingeRadauTwoBeamIntegration();

  case BEAM_INTEGRATION_TAG_HingeEndpoint:
    return new HingeEndpointBeamIntegration();

  case BEAM_INTEGRATION_TAG_UserHinge:
    return new UserDefinedHingeIntegration();

  case BEAM_INTEGRATION_TAG_DistHinge:
    return new DistHingeIntegration();

  case BEAM_INTEGRATION_TAG_RegularizedHinge:
    return new RegularizedHingeIntegration();

  default:
    opserr << "TclPackageClassBroker::getBeamIntegration - ";
    opserr << " - no BeamIntegration type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

UniaxialMaterial *
TclPackageClassBroker::getNewUniaxialMaterial(int classTag)
{
  switch (classTag) {
  case MAT_TAG_SPSW02:
    return new SPSW02(); // SAJalali
  case MAT_TAG_BoucWen:
    return new BoucWenMaterial(); // SAJalali
  case MAT_TAG_ElasticMaterial:
    return new ElasticMaterial();

  case MAT_TAG_Elastic2Material:
    return new Elastic2Material();

  case MAT_TAG_ElasticPPMaterial:
    return new ElasticPPMaterial();

  case MAT_TAG_ElasticMultiLinear:
    return new ElasticMultiLinear();

  case MAT_TAG_ElasticPowerFunc:
    return new ElasticPowerFunc();

  case MAT_TAG_ParallelMaterial:
    return new ParallelMaterial();

  case MAT_TAG_ASD_SMA_3K:
    return new ASD_SMA_3K();

  case MAT_TAG_Concrete01:
    return new Concrete01();

  case MAT_TAG_Concrete02:
    return new Concrete02();

  case MAT_TAG_Concrete04:
    return new Concrete04();

  case MAT_TAG_Concrete06:
    return new Concrete06();

  case MAT_TAG_Concrete07:
    return new Concrete07();

  case MAT_TAG_ConcretewBeta:
    return new ConcretewBeta();

  case MAT_TAG_Steel01:
    return new Steel01();

  case MAT_TAG_Steel02:
    return new Steel02();

  case MAT_TAG_Steel2:
    return new Steel2();

  case MAT_TAG_Steel4:
    return new Steel4();

  case MAT_TAG_OriginCentered:
    return new OriginCentered();

  case MAT_TAG_ReinforcingSteel:
    return new ReinforcingSteel(0);

  case MAT_TAG_Hardening:
    return new HardeningMaterial();

  case MAT_TAG_PySimple1:
    return new PySimple1();

  case MAT_TAG_PyLiq1:
    return new PyLiq1();

  case MAT_TAG_TzSimple1:
    return new TzSimple1();

  case MAT_TAG_PySimple2:
    return new PySimple2();

  case MAT_TAG_TzSimple2:
    return new TzSimple2();

  case MAT_TAG_Fatigue:
    return new FatigueMaterial();

  case MAT_TAG_TzLiq1:
    return new TzLiq1();

  case MAT_TAG_QzSimple1:
    return new QzSimple1();

  case MAT_TAG_QzSimple2:
    return new QzSimple2();

  case MAT_TAG_QzLiq1:
    return new QzLiq1();

  case MAT_TAG_Hysteretic:
    return new HystereticMaterial();

  case MAT_TAG_ModIMKPeakOriented:
    return new ModIMKPeakOriented();

  case MAT_TAG_SnapClough:
    return new Clough();

  case MAT_TAG_LimitState:
    return new LimitStateMaterial();

  case MAT_TAG_EPPGap:
    return new EPPGapMaterial();

  case MAT_TAG_Viscous:
    return new ViscousMaterial();

  case MAT_TAG_ViscousDamper:
    return new ViscousDamper();

  case MAT_TAG_PathIndependent:
    return new PathIndependentMaterial();

  case MAT_TAG_SeriesMaterial:
    return new SeriesMaterial();

  case MAT_TAG_CableMaterial:
    return new CableMaterial();

  case MAT_TAG_ENTMaterial:
    return new ENTMaterial();
#if defined(OPSDEF_UNIAXIAL_FEDEAS)
  case MAT_TAG_FedeasBond1:
    return new FedeasBond1Material();

  case MAT_TAG_FedeasBond2:
    return new FedeasBond2Material();

  case MAT_TAG_FedeasConcrete1:
    return new FedeasConcr1Material();

  case MAT_TAG_FedeasConcrete2:
    return new FedeasConcr2Material();

  case MAT_TAG_FedeasConcrete3:
    return new FedeasConcr3Material();

  case MAT_TAG_FedeasHardening:
    return new FedeasHardeningMaterial();

  case MAT_TAG_FedeasHysteretic1:
    return new FedeasHyster1Material();

  case MAT_TAG_FedeasHysteretic2:
    return new FedeasHyster2Material();

  case MAT_TAG_FedeasSteel1:
    return new FedeasSteel1Material();

  case MAT_TAG_FedeasSteel2:
    return new FedeasSteel2Material();
#endif // OPSDEF_UNIAXIAL_FEDEAS
  case MAT_TAG_DrainBilinear:
    return new DrainBilinearMaterial();

  case MAT_TAG_HyperbolicGapMaterial:
    return new HyperbolicGapMaterial();

  case MAT_TAG_ImpactMaterial:
    return new ImpactMaterial();

  case MAT_TAG_Bilin:
    return new Bilin();

  case MAT_TAG_DrainClough1:
    return new DrainClough1Material();

  case MAT_TAG_DrainClough2:
    return new DrainClough2Material();

  case MAT_TAG_DrainPinch1:
    return new DrainPinch1Material();

  case MAT_TAG_MinMax:
    return new MinMaxMaterial();

  case MAT_TAG_InitStrain:
    return new InitStrainMaterial();

  case MAT_TAG_InitStress:
    return new InitStressMaterial();

  case MAT_TAG_Bond_SP01:
    return new Bond_SP01();

  case MAT_TAG_SimpleFractureMaterial:
    return new SimpleFractureMaterial();

  case MAT_TAG_ConfinedConcrete01:
    return new ConfinedConcrete01();

  case MAT_TAG_HystereticPoly: // Salvatore Sessa
    return new HystereticPoly();

  default:

    UniaxialPackage *matCommands = theUniaxialPackage;
    bool found = false;
    while (matCommands != NULL && found == false) {
      if ((matCommands->classTag == classTag) && (matCommands->funcPtr != 0)) {
        UniaxialMaterial *result = (*(matCommands->funcPtr))();
        return result;
      }
      matCommands = matCommands->next;
    }

    opserr << "TclPackageClassBroker::getNewUniaxialMaterial - ";
    opserr << " - no UniaxialMaterial type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

SectionForceDeformation *
TclPackageClassBroker::getNewSection(int classTag)
{
  switch (classTag) {
  case SEC_TAG_Elastic2d:
    return new ElasticSection2d();

  case SEC_TAG_Elastic3d:
    return new ElasticSection3d();

  case SEC_TAG_ElasticShear2d:
    return new ElasticShearSection2d();

  case SEC_TAG_ElasticShear3d:
    return new ElasticShearSection3d();

  case SEC_TAG_Generic1d:
    return new GenericSection1d();

    // case SEC_TAG_GenericNd:
    // return new GenericSectionNd();

  case SEC_TAG_Aggregator:
    return new SectionAggregator();

    // case SEC_TAG_Fiber:
    // return new FiberSection();

  case SEC_TAG_FiberSection2d:
    return new FiberSection2d();

  case SEC_TAG_FiberSection3d:
    return new FiberSection3d();

  case SEC_TAG_FiberSectionAsym3d:
    return new FiberSectionAsym3d(); // Xinlong Du

  case SEC_TAG_ElasticPlateSection:
    return new ElasticPlateSection();

  case SEC_TAG_ElasticMembranePlateSection:
    return new ElasticMembranePlateSection();

  case SEC_TAG_MembranePlateFiberSection:
    return new MembranePlateFiberSection();

  // start Yuli Huang & Xinzheng Lu LayeredShellFiberSection
  case SEC_TAG_LayeredShellFiberSection:
    return new LayeredShellFiberSection();
    // end Yuli Huang & Xinzheng Lu LayeredShellFiberSection

  case SEC_TAG_Bidirectional:
    return new Bidirectional();

  default:
    opserr << "TclPackageClassBroker::getNewSection - ";
    opserr << " - no section type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

NDMaterial *
TclPackageClassBroker::getNewNDMaterial(int classTag)
{
  switch (classTag) {
  case ND_TAG_ElasticIsotropicPlaneStrain2d:
    return new ElasticIsotropicPlaneStrain2D();

  case ND_TAG_ElasticIsotropicPlaneStress2d:
    return new ElasticIsotropicPlaneStress2D();

  case ND_TAG_ElasticIsotropicAxiSymm:
    return new ElasticIsotropicAxiSymm();

  case ND_TAG_ElasticIsotropicPlateFiber:
    return new ElasticIsotropicPlateFiber();

  case ND_TAG_ElasticIsotropicThreeDimensional:
    return new ElasticIsotropicThreeDimensional();

  case ND_TAG_J2PlaneStrain:
    return new J2PlaneStrain();

  case ND_TAG_J2PlaneStress:
    return new J2PlaneStress();

  case ND_TAG_J2AxiSymm:
    return new J2AxiSymm();

  case ND_TAG_J2PlateFiber:
    return new J2PlateFiber();

  case ND_TAG_J2ThreeDimensional:
    return new J2ThreeDimensional();

  case ND_TAG_PlaneStressMaterial:
    return new PlaneStressMaterial();

  // start Yuli Huang & Xinzheng
  case ND_TAG_PlateRebarMaterial:
    return new PlateRebarMaterial();

  case ND_TAG_PlateFromPlaneStressMaterial:
    return new PlateFromPlaneStressMaterial();

    // case ND_TAG_ConcreteS:
    //    return new ConcreteS();

  case ND_TAG_PlaneStressUserMaterial:
    return new PlaneStressUserMaterial();
    // end Yuli Huang & Xinzheng Lu

  case ND_TAG_PlateFiberMaterial:
    return new PlateFiberMaterial();

  case ND_TAG_FluidSolidPorousMaterial:
    return new FluidSolidPorousMaterial();

  case ND_TAG_PressureDependMultiYield:
    return new PressureDependMultiYield();

  case ND_TAG_PressureDependMultiYield02:
    return new PressureDependMultiYield02();

  case ND_TAG_PressureIndependMultiYield:
    return new PressureIndependMultiYield();

#if defined(OPSDEF_ELEMENT_FEAP)
  case ND_TAG_FeapMaterial03:
    return new FeapMaterial03();
#endif // OPSDEF_ELEMENT_FEAP

  case ND_TAG_ContactMaterial2D:
    return new ContactMaterial2D();

  case ND_TAG_ContactMaterial3D:
    return new ContactMaterial3D();

  case ND_TAG_DruckerPrager3D:
    return new DruckerPrager3D();

  case ND_TAG_DruckerPragerPlaneStrain:
    return new DruckerPragerPlaneStrain();

  case ND_TAG_BoundingCamClay:
    return new BoundingCamClay();

  case ND_TAG_BoundingCamClay3D:
    return new BoundingCamClay3D();

  case ND_TAG_BoundingCamClayPlaneStrain:
    return new BoundingCamClayPlaneStrain();

  case ND_TAG_ManzariDafalias:
    return new ManzariDafalias();

  case ND_TAG_ManzariDafalias3D:
    return new ManzariDafalias3D();

  case ND_TAG_ManzariDafaliasPlaneStrain:
    return new ManzariDafaliasPlaneStrain();

  case ND_TAG_ManzariDafaliasRO:
    return new ManzariDafaliasRO();

  case ND_TAG_ManzariDafalias3DRO:
    return new ManzariDafalias3DRO();

  case ND_TAG_ManzariDafaliasPlaneStrainRO:
    return new ManzariDafaliasPlaneStrainRO();

  case ND_TAG_PM4Sand:
    return new PM4Sand();

  case ND_TAG_PM4Silt:
    return new PM4Silt();

  case ND_TAG_J2CyclicBoundingSurface:
    return new J2CyclicBoundingSurface();

  case ND_TAG_J2CyclicBoundingSurface3D:
    return new J2CyclicBoundingSurface3D();

  case ND_TAG_J2CyclicBoundingSurfacePlaneStrain:
    return new J2CyclicBoundingSurfacePlaneStrain();

  case ND_TAG_InitialStateAnalysisWrapper:
    return new InitialStateAnalysisWrapper();
  case ND_TAG_stressDensity:
    return new stressDensity();
  case ND_TAG_CycLiqCP3D:
    return new CycLiqCP3D();

  case ND_TAG_CycLiqCPPlaneStrain:
    return new CycLiqCPPlaneStrain();

  case ND_TAG_CycLiqCPSP3D:
    return new CycLiqCPSP3D();

  case ND_TAG_CycLiqCPSPPlaneStrain:
    return new CycLiqCPSPPlaneStrain();

  case ND_TAG_InitStressNDMaterial:
    return new InitStressNDMaterial();

  default:
    opserr << "TclPackageClassBroker::getNewNDMaterial - ";
    opserr << " - no NDMaterial type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

Fiber *
TclPackageClassBroker::getNewFiber(int classTag)
{
  switch (classTag) {
  default:
    opserr << "TclPackageClassBroker::getNewFiber - ";
    opserr << " - no Fiber type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

FrictionModel *
TclPackageClassBroker::getNewFrictionModel(int classTag)
{
  switch (classTag) {
  case FRN_TAG_Coulomb:
    return new Coulomb();

  case FRN_TAG_VelDependent:
    return new VelDependent();

  case FRN_TAG_VelPressureDep:
    return new VelPressureDep();

  case FRN_TAG_VelDepMultiLinear:
    return new VelDepMultiLinear();

  case FRN_TAG_VelNormalFrcDep:
    return new VelNormalFrcDep();

  default:
    opserr << "TclPackageClassBroker::getNewFrictionModel - ";
    opserr << " - no FrictionModel type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

ConvergenceTest *
TclPackageClassBroker::getNewConvergenceTest(int classTag)
{
  switch (classTag) {
  case CONVERGENCE_TEST_CTestNormUnbalance:
    return new CTestNormUnbalance();

  case CONVERGENCE_TEST_CTestRelativeNormUnbalance:
    return new CTestRelativeNormUnbalance();

  case CONVERGENCE_TEST_CTestNormDispIncr:
    return new CTestNormDispIncr();

  case CONVERGENCE_TEST_CTestRelativeNormDispIncr:
    return new CTestRelativeNormDispIncr();

  case CONVERGENCE_TEST_CTestRelativeTotalNormDispIncr:
    return new CTestRelativeTotalNormDispIncr();

  case CONVERGENCE_TEST_CTestEnergyIncr:
    return new CTestEnergyIncr();

  case CONVERGENCE_TEST_CTestRelativeEnergyIncr:
    return new CTestRelativeEnergyIncr();

  case CONVERGENCE_TEST_CTestFixedNumIter:
    return new CTestFixedNumIter();

  default:
    opserr << "TclPackageClassBroker::getNewConvergenceTest - ";
    opserr << " - no ConvergenceTest type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

LoadPattern *
TclPackageClassBroker::getNewLoadPattern(int classTag)
{
  switch (classTag) {
  case PATTERN_TAG_LoadPattern:
    return new LoadPattern();

  case PATTERN_TAG_UniformExcitation:
    return new UniformExcitation();

  case PATTERN_TAG_MultiSupportPattern:
    return new MultiSupportPattern();

#ifdef OPSDEF_DRM
  case PATTERN_TAG_DRMLoadPattern:
    return new DRMLoadPatternWrapper();
#endif // OPSDEF_DRM

#ifdef _H5DRM
  case PATTERN_TAG_H5DRM:
    return new H5DRM();
#endif
  default:
    opserr << "TclPackageClassBroker::getPtrLoadPattern - ";
    opserr << " - no Load type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

GroundMotion *
TclPackageClassBroker::getNewGroundMotion(int classTag)
{
  switch (classTag) {

  case GROUND_MOTION_TAG_GroundMotion:
    return new GroundMotion(GROUND_MOTION_TAG_GroundMotion);

  case GROUND_MOTION_TAG_InterpolatedGroundMotion:
    return new GroundMotion(GROUND_MOTION_TAG_InterpolatedGroundMotion);

  default:
    opserr << "TclPackageClassBroker::getPtrGroundMotion - ";
    opserr << " - no Load type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

TimeSeries *
TclPackageClassBroker::getNewTimeSeries(int classTag)
{
  switch (classTag) {
  case TSERIES_TAG_LinearSeries:
    return new LinearSeries;

  case TSERIES_TAG_RectangularSeries:
    return new RectangularSeries;

  case TSERIES_TAG_PathTimeSeries:
    return new PathTimeSeries;

  case TSERIES_TAG_PathSeries:
    return new PathSeries;

  case TSERIES_TAG_ConstantSeries:
    return new ConstantSeries;

  case TSERIES_TAG_TriangleSeries:
    return new TriangleSeries;

  case TSERIES_TAG_TrigSeries:
    return new TrigSeries;

  default:
    opserr << "TclPackageClassBroker::getPtrTimeSeries - ";
    opserr << " - no Load type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

TimeSeriesIntegrator *
TclPackageClassBroker::getNewTimeSeriesIntegrator(int classTag)
{
  switch (classTag) {
  case TIMESERIES_INTEGRATOR_TAG_Trapezoidal:
    return new TrapezoidalTimeSeriesIntegrator();

  default:
    opserr << "TclPackageClassBroker::getPtrTimeSeriesIntegrator - ";
    opserr << " - no Load type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

Matrix *
TclPackageClassBroker::getPtrNewMatrix(int classTag, int noRows, int noCols)
{
  switch (classTag) {
  case MATRIX_TAG_Matrix:
    return new Matrix(noRows, noCols);

  default:
    opserr << "TclPackageClassBroker::getPtrNewMatrix - ";
    opserr << " - no NodalLoad type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

Vector *
TclPackageClassBroker::getPtrNewVector(int classTag, int size)
{
  switch (classTag) {
  case VECTOR_TAG_Vector:
    return new Vector(size);

  default:
    opserr << "TclPackageClassBroker::getPtrNewVector - ";
    opserr << " - no Vector type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

ID *
TclPackageClassBroker::getPtrNewID(int classTag, int size)
{
  switch (classTag) {
  case ID_TAG_ID:
    return new ID(size);

  default:
    opserr << "TclPackageClassBroker::getPtrNewID - ";
    opserr << " - no ID type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

/*****************************************
 *
 * METHODS TO GET NEW OUTPUT CLASS OBJECTS
 *
 *****************************************/

OPS_Stream *
TclPackageClassBroker::getPtrNewStream(int classTag)
{
  switch (classTag) {
  case OPS_STREAM_TAGS_StandardStream:
    return new StandardStream();

  case OPS_STREAM_TAGS_FileStream:
    return new FileStream();

  case OPS_STREAM_TAGS_XmlFileStream:
    return new XmlFileStream();

  case OPS_STREAM_TAGS_DataFileStream:
    return new DataFileStream();

  case OPS_STREAM_TAGS_DataFileStreamAdd:
    return new DataFileStreamAdd();

  case OPS_STREAM_TAGS_BinaryFileStream:
    return new BinaryFileStream();

  case OPS_STREAM_TAGS_DatabaseStream:
    return new DatabaseStream();

  case OPS_STREAM_TAGS_DummyStream:
    return new DummyStream();

  default:
    opserr << "TclPackageClassBroker::getPtrNewStream - ";
    opserr << " - no DataOutputHandler type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

Recorder *
TclPackageClassBroker::getPtrNewRecorder(int classTag)
{
  switch (classTag) {
  case RECORDER_TAGS_ElementRecorder:
    return new ElementRecorder();

  case RECORDER_TAGS_NodeRecorder:
    return new NodeRecorder();

  case RECORDER_TAGS_EnvelopeNodeRecorder:
    return new EnvelopeNodeRecorder();

  case RECORDER_TAGS_EnvelopeElementRecorder:
    return new EnvelopeElementRecorder();

  case RECORDER_TAGS_VTK_Recorder:
    return new VTK_Recorder();

  case RECORDER_TAGS_DriftRecorder:
    return new DriftRecorder();

  case RECORDER_TAGS_TclFeViewer:
    return 0;

  case RECORDER_TAGS_GmshRecorder:
    return new GmshRecorder();

    //        case RECORDER_TAGS_MPCORecorder:
    //          return new MPCORecorder();

  default:
    opserr << "TclPackageClassBroker::getNewRecordr - ";
    opserr << " - no Recorder type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

/*****************************************
 *
 * METHODS TO GET NEW ANALYSIS CLASSES
 *
 *****************************************/

ConstraintHandler *
TclPackageClassBroker::getNewConstraintHandler(int classTag)
{
  switch (classTag) {
  case HANDLER_TAG_PlainHandler:
    return new PlainHandler();

  case HANDLER_TAG_PenaltyConstraintHandler:
    return new PenaltyConstraintHandler(1.0e12, 1.0e12);

  case HANDLER_TAG_LagrangeConstraintHandler:
    return new LagrangeConstraintHandler(1.0, 1.0);

  case HANDLER_TAG_TransformationConstraintHandler:
    return new TransformationConstraintHandler();

  default:
    opserr << "TclPackageClassBroker::getNewConstraintHandler - ";
    opserr << " - no ConstraintHandler type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

DOF_Numberer *
TclPackageClassBroker::getNewNumberer(int classTag)
{
  switch (classTag) {
  case NUMBERER_TAG_DOF_Numberer:
    return new DOF_Numberer();

  case NUMBERER_TAG_PlainNumberer:
    return new PlainNumberer();

#ifdef _PARALLEL_PROCESSING
  case NUMBERER_TAG_ParallelNumberer:
    return new ParallelNumberer();
#endif

  default:
    opserr << "TclPackageClassBroker::getNewConstraintHandler - ";
    opserr << " - no ConstraintHandler type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

AnalysisModel *
TclPackageClassBroker::getNewAnalysisModel(int classTag)
{
  switch (classTag) {
  case AnaMODEL_TAGS_AnalysisModel:
    return new AnalysisModel();

  default:
    opserr << "TclPackageClassBroker::getNewAnalysisModel - ";
    opserr << " - no AnalysisModel type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

EquiSolnAlgo *
TclPackageClassBroker::getNewEquiSolnAlgo(int classTag)
{
  switch (classTag) {
  case EquiALGORITHM_TAGS_Linear:
    return new Linear();

  case EquiALGORITHM_TAGS_NewtonRaphson:
    return new NewtonRaphson();

  case EquiALGORITHM_TAGS_NewtonLineSearch:
    return new NewtonLineSearch();

  case EquiALGORITHM_TAGS_KrylovNewton:
    return new KrylovNewton();

  case EquiALGORITHM_TAGS_AcceleratedNewton:
    return new AcceleratedNewton();

  case EquiALGORITHM_TAGS_ModifiedNewton:
    return new ModifiedNewton(CURRENT_TANGENT);

  case EquiALGORITHM_TAGS_Broyden:
    return new Broyden();

  default:
    opserr << "TclPackageClassBroker::getNewEquiSolnAlgo - ";
    opserr << " - no EquiSolnAlgo type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

Accelerator *
TclPackageClassBroker::getAccelerator(int classTag)
{
  switch (classTag) {

  case ACCELERATOR_TAGS_Krylov:
    return new KrylovAccelerator;
  case ACCELERATOR_TAGS_Raphson:
    return new RaphsonAccelerator;

  default:
    opserr << "TclPackageClassBroker::getAccelerator - ";
    opserr << " - no EquiSolnAlgo type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

LineSearch *
TclPackageClassBroker::getLineSearch(int classTag)
{
  switch (classTag) {

  case LINESEARCH_TAGS_BisectionLineSearch:
    return new BisectionLineSearch();

  case LINESEARCH_TAGS_InitialInterpolatedLineSearch:
    return new InitialInterpolatedLineSearch();

  case LINESEARCH_TAGS_RegulaFalsiLineSearch:
    return new RegulaFalsiLineSearch();

  case LINESEARCH_TAGS_SecantLineSearch:
    return new SecantLineSearch();
  default:
    opserr << "TclPackageClassBroker::getNewEquiSolnAlgo - ";
    opserr << " - no EquiSolnAlgo type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

DomainDecompAlgo *
TclPackageClassBroker::getNewDomainDecompAlgo(int classTag)
{
  switch (classTag) {
  case DomDecompALGORITHM_TAGS_DomainDecompAlgo:
    return new DomainDecompAlgo();

  default:
    opserr << "TclPackageClassBroker::getNewDomainDecompAlgo - ";
    opserr << " - no DomainDecompAlgo type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

StaticIntegrator *
TclPackageClassBroker::getNewStaticIntegrator(int classTag)
{
  switch (classTag) {
  case INTEGRATOR_TAGS_LoadControl:
    return new LoadControl(1.0, 1, 1.0, .10); // must recvSelf
#if 0
  case INTEGRATOR_TAGS_StagedLoadControl:
    return new StagedLoadControl(1.0, 1, 1.0, .10); // must recvSelf
#endif
#ifdef _PARALLEL_PROCESSING
  case INTEGRATOR_TAGS_DistributedDisplacementControl:
    return new DistributedDisplacementControl(); // must recvSelf
#endif

  case INTEGRATOR_TAGS_ArcLength:
    return new ArcLength(1.0); // must recvSelf

  default:
    opserr << "TclPackageClassBroker::getNewStaticIntegrator - ";
    opserr << " - no StaticIntegrator type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

TransientIntegrator *
TclPackageClassBroker::getNewTransientIntegrator(int classTag)
{
  switch (classTag) {
  case INTEGRATOR_TAGS_AlphaOS:
    return new AlphaOS();

  case INTEGRATOR_TAGS_AlphaOS_TP:
    return new AlphaOS_TP();

  case INTEGRATOR_TAGS_AlphaOSGeneralized:
    return new AlphaOSGeneralized();

  case INTEGRATOR_TAGS_AlphaOSGeneralized_TP:
    return new AlphaOSGeneralized_TP();

  case INTEGRATOR_TAGS_CentralDifference:
    return new CentralDifference(); // must recvSelf

  case INTEGRATOR_TAGS_CentralDifferenceAlternative:
    return new CentralDifferenceAlternative(); // must recvSelf

  case INTEGRATOR_TAGS_CentralDifferenceNoDamping:
    return new CentralDifferenceNoDamping(); // must recvSelf

  case INTEGRATOR_TAGS_Collocation:
    return new Collocation();

  case INTEGRATOR_TAGS_CollocationHSFixedNumIter:
    return new CollocationHSFixedNumIter();

  case INTEGRATOR_TAGS_CollocationHSIncrLimit:
    return new CollocationHSIncrLimit();

  case INTEGRATOR_TAGS_CollocationHSIncrReduct:
    return new CollocationHSIncrReduct();

  case INTEGRATOR_TAGS_HHT:
    return new HHT();

  case INTEGRATOR_TAGS_HHT_TP:
    return new HHT_TP();

  case INTEGRATOR_TAGS_HHTExplicit:
    return new HHTExplicit();

  case INTEGRATOR_TAGS_HHTExplicit_TP:
    return new HHTExplicit_TP();

  case INTEGRATOR_TAGS_HHTGeneralized:
    return new HHTGeneralized();

  case INTEGRATOR_TAGS_HHTGeneralized_TP:
    return new HHTGeneralized_TP();

  case INTEGRATOR_TAGS_HHTGeneralizedExplicit:
    return new HHTGeneralizedExplicit();

  case INTEGRATOR_TAGS_HHTGeneralizedExplicit_TP:
    return new HHTGeneralizedExplicit_TP();

  case INTEGRATOR_TAGS_HHTHSFixedNumIter:
    return new HHTHSFixedNumIter();

  case INTEGRATOR_TAGS_HHTHSFixedNumIter_TP:
    return new HHTHSFixedNumIter_TP();

  case INTEGRATOR_TAGS_HHTHSIncrLimit:
    return new HHTHSIncrLimit();

  case INTEGRATOR_TAGS_HHTHSIncrLimit_TP:
    return new HHTHSIncrLimit_TP();

  case INTEGRATOR_TAGS_HHTHSIncrReduct:
    return new HHTHSIncrReduct();

  case INTEGRATOR_TAGS_HHTHSIncrReduct_TP:
    return new HHTHSIncrReduct_TP();

  case INTEGRATOR_TAGS_KRAlphaExplicit:
    return new KRAlphaExplicit();

  case INTEGRATOR_TAGS_KRAlphaExplicit_TP:
    return new KRAlphaExplicit_TP();

  case INTEGRATOR_TAGS_Newmark:
    return new Newmark();
#if 0
  case INTEGRATOR_TAGS_StagedNewmark:
    return new StagedNewmark();
#endif
  case INTEGRATOR_TAGS_NewmarkExplicit:
    return new NewmarkExplicit();

  case INTEGRATOR_TAGS_NewmarkHSFixedNumIter:
    return new NewmarkHSFixedNumIter();

  case INTEGRATOR_TAGS_NewmarkHSIncrLimit:
    return new NewmarkHSIncrLimit();

  case INTEGRATOR_TAGS_NewmarkHSIncrReduct:
    return new NewmarkHSIncrReduct();
#if 0
  case INTEGRATOR_TAGS_PFEMIntegrator:
    return new PFEMIntegrator();
#endif
  case INTEGRATOR_TAGS_TRBDF2:
    return new TRBDF2();

  case INTEGRATOR_TAGS_TRBDF3:
    return new TRBDF3();

  case INTEGRATOR_TAGS_WilsonTheta:
    return new WilsonTheta();

  default:
    opserr << "TclPackageClassBroker::getNewTransientIntegrator - ";
    opserr << " - no TransientIntegrator type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

IncrementalIntegrator *
TclPackageClassBroker::getNewIncrementalIntegrator(int classTag)
{
  switch (classTag) {
  case INTEGRATOR_TAGS_LoadControl:
    return new LoadControl(1.0, 1, 1.0, 1.0); // must recvSelf

  case INTEGRATOR_TAGS_ArcLength:
    return new ArcLength(1.0); // must recvSelf

  case INTEGRATOR_TAGS_Newmark:
    return new Newmark();

#ifdef _PARALLEL_PROCESSING
  case INTEGRATOR_TAGS_DistributedDisplacementControl:
    return new DistributedDisplacementControl(); // must recvSelf
#endif

  default:
    opserr << "TclPackageClassBroker::getNewIncrementalIntegrator - ";
    opserr << " - no IncrementalIntegrator type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

LinearSOE *
TclPackageClassBroker::getNewLinearSOE(int classTagSOE)
{
  LinearSOE *theSOE = nullptr;

  switch (classTagSOE) {

  case LinSOE_TAGS_SparseGenColLinSOE:
    theSOE = new SparseGenColLinSOE();
    return theSOE;

#ifdef _PARALLEL_PROCESSING

#ifdef _MUMPS
  case LinSOE_TAGS_MumpsParallelSOE:
    theSOE = new MumpsParallelSOE();
    return theSOE;
#endif

#ifdef _PETSC
  case LinSOE_TAGS_PetscSOE:
    theSOE = new PetscSOE(*(new PetscSolver()));
    return theSOE;
#endif

  case LinSOE_TAGS_DistributedBandGenLinSOE:

    theSOE = new DistributedBandGenLinSOE();
    return theSOE;

  case LinSOE_TAGS_DistributedBandSPDLinSOE:

    theSOE = new DistributedBandSPDLinSOE();
    return theSOE;

  case LinSOE_TAGS_DistributedProfileSPDLinSOE:

    theSOE = new DistributedProfileSPDLinSOE();
    return theSOE;

  case LinSOE_TAGS_DistributedDiagonalSOE:

    theSOE = new DistributedDiagonalSOE();
    return theSOE;

  case LinSOE_TAGS_DistributedSparseGenColLinSOE:

    theSOE = new DistributedSparseGenColLinSOE();
    return theSOE;

#endif

  default:
    opserr << "TclPackageClassBroker::getNewLinearSOE - ";
    opserr << " - no LinearSOE type exists for class tag ";
    opserr << classTagSOE << endln;
    return 0;
  }
}

EigenSOE *
TclPackageClassBroker::getNewEigenSOE(int classTagSOE)
{
  EigenSOE *theSOE = 0;

  switch (classTagSOE) {

  case EigenSOE_TAGS_ArpackSOE:
    theSOE = new ArpackSOE();
    return theSOE;

  default:
    opserr << "TclPackageClassBroker::getNewEigenSOE - ";
    opserr << " - no EigenSOE type exists for class tag ";
    opserr << classTagSOE << endln;
    return 0;
  }
}

DomainSolver *
TclPackageClassBroker::getNewDomainSolver(void)
{
  return lastDomainSolver;
}

LinearSOE *
TclPackageClassBroker::getPtrNewDDLinearSOE(int classTagSOE,
                                            int classTagDDSolver)
{
  ProfileSPDLinSubstrSolver *theProfileSPDSolver = 0;

  switch (classTagSOE) {
  case LinSOE_TAGS_ProfileSPDLinSOE:

    if (classTagDDSolver == SOLVER_TAGS_ProfileSPDLinSubstrSolver) {
      theProfileSPDSolver = new ProfileSPDLinSubstrSolver();
      LinearSOE *theSOE = new ProfileSPDLinSOE(*theProfileSPDSolver);
      lastDomainSolver = theProfileSPDSolver;
      return theSOE;
    } else {
      opserr << "TclPackageClassBroker::getNewLinearSOE - ";
      opserr << " - no ProfileSPD Domain Solver type exists for class tag ";
      opserr << classTagDDSolver << endln;
      return 0;
    }

  default:
    opserr << "TclPackageClassBroker::getNewLinearSOE - ";
    opserr << " - no LinearSOE type exists for class tag ";
    opserr << classTagSOE << endln;
    return 0;
  }
}

DomainDecompositionAnalysis *
TclPackageClassBroker::getNewDomainDecompAnalysis(int classTag,
                                                  Subdomain &theSubdomain)
{
  switch (classTag) {
//case DomDecompANALYSIS_TAGS_DomainDecompositionAnalysis:
//  return new DomainDecompositionAnalysis(theSubdomain);

#ifdef _PARALLEL_PROCESSING
  case ANALYSIS_TAGS_StaticDomainDecompositionAnalysis:
    return new StaticDomainDecompositionAnalysis(theSubdomain);

  case ANALYSIS_TAGS_TransientDomainDecompositionAnalysis:
    return new TransientDomainDecompositionAnalysis(theSubdomain);
#endif

  default:
    opserr << "TclPackageClassBroker::getNewDomainDecompAnalysis ";
    opserr << " - no DomainDecompAnalysis type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

Subdomain *
TclPackageClassBroker::getSubdomainPtr(int classTag)
{
  opserr << "TclPackageClassBroker: NOT IMPLEMENTED YET";
  return 0;
}

int
TclPackageClassBroker::addUniaxialMaterial(int classTag, const char *lib,
                                           const char *funcName,
                                           UniaxialMaterial *(*funcPtr)(void))
{
  // check to see if it's already added

  UniaxialPackage *matCommands = theUniaxialPackage;
  bool found = false;
  while (matCommands != NULL && found == false) {
    if ((strcmp(lib, matCommands->libName) == 0) &&
        (strcmp(funcName, matCommands->funcName) == 0)) {
      return 0;
    }
  }

  //
  // if funPtr == 0; go get the handle
  //

  void *libHandle;
  if (funcPtr == 0) {
    if (getLibraryFunction(lib, funcName, &libHandle, (void **)&funcPtr) != 0) {
      opserr << "TclPackageClassBroker::addUniaxialMaterial - could not find "
                "function\n";
      return -1;
    }
  }

  //
  // add the new funcPtr
  //

  char *libNameCopy = new char[strlen(lib) + 1];
  char *funcNameCopy = new char[strlen(funcName) + 1];
  UniaxialPackage *theMat = new UniaxialPackage;
  if (libNameCopy == 0 || funcNameCopy == 0 || theMat == 0) {
    opserr << "TclPackageClassBroker::addUniaxialMaterial - could not add lib, "
              "out of memory\n";
    return -1;
  }
  strcpy(libNameCopy, lib);
  strcpy(funcNameCopy, funcName);

  theMat->classTag = classTag;
  theMat->funcName = funcNameCopy;
  theMat->libName = libNameCopy;
  theMat->funcPtr = funcPtr;
  theMat->next = theUniaxialPackage;
  theUniaxialPackage = theMat;

  return 0;
}

Parameter *
TclPackageClassBroker::getParameter(int classTag)
{
  Parameter *theRes = 0;

  switch (classTag) {
  case PARAMETER_TAG_Parameter:
    theRes = new Parameter;
    break;

  case PARAMETER_TAG_ElementParameter:
    theRes = new ElementParameter;
    break;

  case PARAMETER_TAG_MaterialStageParameter:
    theRes = new MaterialStageParameter();
    break;

  case PARAMETER_TAG_MatParameter:
    theRes = new MatParameter();
    break;

  case PARAMETER_TAG_InitialStateParameter:
    theRes = new InitialStateParameter();
    break;

  case PARAMETER_TAG_ElementStateParameter:
    theRes = new ElementStateParameter();
    break;

  default:;
  }

  return theRes;
}
