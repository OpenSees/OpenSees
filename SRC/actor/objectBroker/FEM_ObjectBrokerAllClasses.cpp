/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// Written: fmk
// Revision: A
//
// Purpose: This file contains the class definition for FEM_ObjectBrokerAllClasses.
// FEM_ObjectBrokerAllClasses is is an object broker class for the finite element
// method. All methods are virtual to allow for subclasses; which can be
// used by programmers when introducing new subclasses of the main objects.

#ifdef _PARALLEL_PROCESSING
#include <mpi.h>
#endif

#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#endif

#include <FEM_ObjectBrokerAllClasses.h>

// ActorTypes
#include "domain/subdomain/ActorSubdomain.h"

// Convergence tests
#include "convergenceTest/CTestNormUnbalance.h"
#include "convergenceTest/CTestRelativeNormUnbalance.h"
#include "convergenceTest/CTestNormDispIncr.h"
#include "convergenceTest/CTestRelativeNormDispIncr.h"
#include "convergenceTest/CTestRelativeTotalNormDispIncr.h"
#include "convergenceTest/CTestEnergyIncr.h" 
#include "convergenceTest/CTestRelativeEnergyIncr.h"
#include "convergenceTest/CTestFixedNumIter.h"

// graph numbering schemes
#include "graph/numberer/RCM.h"
#include "graph/numberer/MyRCM.h"
#include "graph/numberer/SimpleNumberer.h"


// uniaxial material model header files
#include "ElasticBilin.h"
#include "BoucWenMaterial.h"
#include "BWBN.h"
#include "BoucWenOriginal.h"
#include "BoucWenInfill.h"
#include "SPSW02.h"			//SAJalali
#include "ElasticMaterial.h"
#include "ElasticMultiLinear.h"
#include "ElasticPowerFunc.h"
#include "Elastic2Material.h"
#include "ElasticPPMaterial.h"
#include "ParallelMaterial.h"
#include "DamperMaterial.h"
#include "PenaltyMaterial.h"
#include "MultiplierMaterial.h"
#include "TensionOnlyMaterial.h"
#include "ASD_SMA_3K.h"
#include "ASDConcrete1DMaterial.h"
#include "Concrete01.h"
#include "Concrete01WithSITC.h"
#include "Concrete02.h"
#include "Concrete02IS.h"
#include "Concrete04.h"
#include "Concrete06.h" 
#include "Concrete07.h"
#include "ConcretewBeta.h"
#include "TDConcrete.h"
#include "TDConcreteNL.h"
#include "TDConcreteEXP.h"
#include "TDConcreteMC10.h"
#include "TDConcreteMC10NL.h"
#include "CreepMaterial.h"
#include "OriginCentered.h"
#include "Steel01.h"
#include "Steel02.h"
#include "SteelMPF.h"
#include "Steel2.h"
#include "Steel4.h"
#include "RambergOsgoodSteel.h"
#include "AxialSp.h"
#include "AxialSpHD.h"
#include "SelfCenteringMaterial.h"
#include "SteelFractureDI.h"
#include "FatigueMaterial.h"
#include "ReinforcingSteel.h"
#include "HardeningMaterial.h"
#include "HystereticMaterial.h"
#include "HystereticSMMaterial.h"
#include "OOHystereticMaterial.h"
#include "EPPGapMaterial.h"
#include "HookGap.h"
#include "ViscousMaterial.h"
#include "ViscousDamper.h"
#include "APDVFD.h"
#include "APDMD.h"
#include "APDFMD.h"
#include "BilinearOilDamper.h"
#include "Maxwell.h"
#include "ContinuumUniaxial.h"
#include "PathIndependentMaterial.h"
#include "BackboneMaterial.h"
#include "SeriesMaterial.h"
#include "CableMaterial.h"
#include "ENTMaterial.h"
#include "GNGMaterial.h"
#include "Ratchet.h"
#include "MinMaxMaterial.h"
#include "ModIMKPeakOriented.h"
#include "IMKBilin.h"
#include "IMKPeakOriented.h"
#include "IMKPinching.h"
#include "UVCuniaxial.h"
#include "SLModel.h"
#include "SteelDRC.h"
#include "Cast.h"
#include "snap/Clough.h"
#include "limitState/LimitStateMaterial.h"
#include "InitStressMaterial.h"
#include "InitStrainMaterial.h"
#include "Bond_SP01.h"
#include "FRCC.h"
#include "SimpleFractureMaterial.h"
#include "ConfinedConcrete01.h"
#include "FRPConfinedConcrete.h"
#include "ConcreteCM.h"
#include "BarSlipMaterial.h"

#include <HystereticPoly.h>					// Salvatore Sessa 14-Jan-2021
#include <HystereticSmooth.h>					// Salvatore Sessa Apr-19-2022
#include <HystereticAsym.h>					// Salvatore Sessa Apr-21-2022


#include "DowelType.h"
#include "DuctileFracture.h" // Kuanshi Zhong
#include <GMG_CyclicReinforcedConcrete.h>      // Rasool Ghorbani

//PY springs: RWBoulanger and BJeremic
#include "PY/PySimple1.h"
#include "PY/TzSimple1.h"
#include "PY/QzSimple1.h"
#include "PY/PySimple2.h"
#include "PY/TzSimple2.h"
#include "PY/QzSimple2.h"
#include "PY/PyLiq1.h"
#include "PY/TzLiq1.h"
#include "PY/QzLiq1.h"
// Unified CPT based method
#include "TzSandCPT.h"
#include "QbSandCPT.h"

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
#include "Hertzdamp.h"
#include "JankowskiImpact.h"
#include "ViscoelasticGap.h"
#include "Pinching4Material.h"
#include "SAWSMaterial.h"
#include "KikuchiAikenHDR.h"
#include "KikuchiAikenLRB.h"
#include "CFSSSWP.h"
#include "CFSWSWP.h"

// Sections
#include "ElasticSection2d.h"
#include "ElasticSection3d.h"
#include "ElasticShearSection2d.h"
#include "ElasticBDShearSection2d.h"
#include "ElasticShearSection3d.h"
#include "ElasticTubeSection3d.h"
#include "ElasticWarpingShearSection2d.h"
#include "GenericSection1d.h"
//#include "GenericSectionNd.h"
#include "SectionAggregator.h"
#include "ParallelSection.h"
//#include "FiberSection.h"
#include "FiberSection2d.h"
#include "FiberSection3d.h"
#include "FiberSection2dThermal.h"
#include "FiberSection3dThermal.h"
#include "NDFiberSection2d.h"
#include "NDFiberSectionWarping2d.h"
#include "NDFiberSection3d.h"
#include "FiberSectionAsym3d.h" //Xinlong Du
#include "FiberSectionWarping3d.h"
#include "ElasticPlateSection.h"
#include "ElasticMembranePlateSection.h"
#include "MembranePlateFiberSection.h"
#include "DoubleMembranePlateFiberSection.h"
#include "Bidirectional.h"
#include "Elliptical2.h"
#include "Isolator2spring.h"
#include "LayeredShellFiberSection.h" // Yuli Huang & Xinzheng Lu 
#include "LayeredMembraneSection/ReinforcedConcreteLayeredMembraneSection.h" // M. J. Nunez
#include "LayeredMembraneSection/LayeredMembraneSection.h" // M. J. Nunez
#include "LayeredMembraneSection/ElasticMembraneSection.h" // M. J. Nunez

// NDMaterials
#include "ElasticIsotropicPlaneStrain2D.h"
#include "ElasticIsotropicPlaneStress2D.h"
#include "ElasticIsotropicPlateFiber.h"
#include "ElasticIsotropicBeamFiber.h"
#include "ElasticIsotropicBeamFiber2d.h"
#include "ElasticIsotropicAxiSymm.h"
#include "ElasticIsotropicThreeDimensional.h"
#include "ElasticOrthotropicThreeDimensional.h"
#include "ElasticOrthotropicPlaneStress.h"
#include "J2PlaneStrain.h"
#include "J2PlaneStress.h"
#include "J2PlateFiber.h"
#include "J2PlateFibre.h"
#include "J2BeamFiber2d.h"
#include "J2BeamFiber3d.h"
#include "J2AxiSymm.h"
#include "J2ThreeDimensional.h"
#include "SimplifiedJ2.h"
#include "PlaneStressSimplifiedJ2.h"
#include "PlaneStrainMaterial.h"
#include "PlaneStressMaterial.h"
#include "PlateFiberMaterial.h"
#include "OrthotropicMaterial.h"
#include "Series3DMaterial.h"
#include "Parallel3DMaterial.h"
#include "PlaneStressRebarMaterial.h"
#include "PlaneStressLayeredMaterial.h"
//start Yuli Huang & Xinzheng L
#include "PlateRebarMaterial.h"
#include "PlateFromPlaneStressMaterial.h"
//#include "ConcreteS.h"
#include "PlaneStressUserMaterial.h"
//end Yuli Huang & Xinzheng Lu
#include "feap/FeapMaterial03.h"
#include "CycLiqCP3D.h"
#include "CycLiqCPPlaneStrain.h"
#include "CycLiqCPSP3D.h"
#include "CycLiqCPSPPlaneStrain.h"


#include "soil/FluidSolidPorousMaterial.h"
#include "soil/PressureDependMultiYield.h"
#include "soil/PressureDependMultiYield02.h"
#include "soil/PressureDependMultiYield03.h"
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
#include "UANDESmaterials/SAniSandMS.h"
#include "UANDESmaterials/SAniSandMS3D.h"
#include "UANDESmaterials/SAniSandMSPlaneStrain.h"
#include "UWmaterials/PM4Sand.h"
#include "UWmaterials/PM4Silt.h"
#include "J2CyclicBoundingSurface.h"
#include "J2CyclicBoundingSurface3D.h"
#include "J2CyclicBoundingSurfacePlaneStrain.h"
#include "UWmaterials/InitialStateAnalysisWrapper.h"
#include "stressDensityModel/stressDensity.h"
#include "InitStressNDMaterial.h"
#include "InitStrainNDMaterial.h"
#include "MinMaxNDMaterial.h"
#include "ASDConcrete3DMaterial.h"
#include "PlasticDamageConcrete3d.h"
#include "PlasticDamageConcretePlaneStress.h"
#include "ConcreteS.h"
#include "OrthotropicRotatingAngleConcreteT2DMaterial01/OrthotropicRotatingAngleConcreteT2DMaterial01.h" // M. J. Nunez
#include "SmearedSteelDoubleLayerT2DMaterial01/SmearedSteelDoubleLayerT2DMaterial01.h" // M. J. Nunez

// Fibers
#include "fiber/UniaxialFiber2d.h"
#include "fiber/UniaxialFiber3d.h"

// friction models
#include "frictionBearing/frictionModel/Coulomb.h"
#include "frictionBearing/frictionModel/VelDependent.h"
#include "frictionBearing/frictionModel/VelPressureDep.h"
#include "frictionBearing/frictionModel/VelDepMultiLinear.h"
#include "frictionBearing/frictionModel/VelNormalFrcDep.h"

#include "CappedBackbone.h"
#include "LinearCappedBackbone.h"
#include "MaterialBackbone.h"
#include "ArctangentBackbone.h"
#include "TrilinearBackbone.h"
#include "MultilinearBackbone.h"
#include "ManderBackbone.h"
//#include "KentParkBackbone.h"
#include "RaynorBackbone.h"
#include "ReeseStiffClayBelowWS.h"
#include "ReeseStiffClayAboveWS.h"
#include "ReeseSoftClayBackbone.h"
#include "ReeseSandBackbone.h"
#include "VuggyLimestone.h"
#include "CementedSoil.h"
#include "WeakRock.h"
#include "LiquefiedSand.h"

#include "ConstantUnloadingRule.h"
#include "TakedaUnloadingRule.h"
#include "EnergyUnloadingRule.h"
#include "KarsanUnloadingRule.h"

#include "ConstantStiffnessDegradation.h"
#include "DuctilityStiffnessDegradation.h"
#include "EnergyStiffnessDegradation.h"
#include "PincheiraStiffnessDegradation.h"

#include "ConstantStrengthDegradation.h"
#include "ACIStrengthDegradation.h"
#include "DuctilityStrengthDegradation.h"
#include "PetrangeliStrengthDegradation.h"
#include "EnergyStrengthDegradation.h"
#include "SectionStrengthDegradation.h"


// element header files
#include "Element.h"
#include "truss/Truss.h"
#include "truss/Truss2.h"
#include "truss/TrussSection.h"
#include "truss/CorotTruss.h"
#include "truss/CorotTrussSection.h"
#include "truss/InertiaTruss.h"
#include "zeroLength/ZeroLength.h"
#include "zeroLength/ZeroLengthSection.h"
#include "zeroLength/CoupledZeroLength.h"
#include "zeroLength/ZeroLengthContact2D.h"
#include "zeroLength/ZeroLengthContact3D.h"
#include "zeroLength/ZeroLengthContactNTS2D.h"
#include "zeroLength/ZeroLengthInterface2D.h"
#include "zeroLength/ZeroLengthContactASDimplex.h"
//#include "ZeroLengthND.h"

#include "fourNodeQuad/FourNodeQuad.h"
#include "fourNodeQuad/FourNodeQuad3d.h"
#include "fourNodeQuad/EnhancedQuad.h"
#include "fourNodeQuad/NineNodeMixedQuad.h"
#include "fourNodeQuad/NineNodeQuad.h"
#include "fourNodeQuad/EightNodeQuad.h"
#include "fourNodeQuad/ConstantPressureVolumeQuad.h"
#include "elasticBeamColumn/ElasticBeam2d.h"
#include "elasticBeamColumn/ElasticBeam3d.h"
#include "componentElement/ComponentElement2d.h"
#include "componentElement/ComponentElement3d.h"
#include "elasticBeamColumn/ModElasticBeam2d.h"			//SAJalali
#include "elasticBeamColumn/ModElasticBeam3d.h"
#include "elasticBeamColumn/ElasticTimoshenkoBeam2d.h"
#include "elasticBeamColumn/ElasticTimoshenkoBeam3d.h"
#include "forceBeamColumn/ForceBeamColumn2d.h"
#include "forceBeamColumn/ForceBeamColumn3d.h"
#include "dispBeamColumn/TimoshenkoBeamColumn2d.h"
#include "dispBeamColumn/TimoshenkoBeamColumn3d.h"
#include "gradientInelasticBeamColumn/GradientInelasticBeamColumn2d.h"
#include "gradientInelasticBeamColumn/GradientInelasticBeamColumn3d.h"
#include "triangle/Tri31.h"
#include "fourNodeQuad/SixNodeTri.h"

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
#include "UWelements/Quad4FiberOverlay.h"
#include "UWelements/Brick8FiberOverlay.h"
#include "EmbeddedBeamInterfaceL.h"
#include "SurfaceLoad.h"
#include "TriSurfaceLoad.h"

#include "PML/PML2D.h"
#include "PML/PML3D.h"
#include "PML/PML2D_3.h"
#include "PML/PML2D_5.h"
#include "PML/PML2D_12.h"
#include "PML/PML2DVISCOUS.h"


#include "UP-ucsd/Nine_Four_Node_QuadUP.h"
#include "UP-ucsd/BrickUP.h"
#include "UP-ucsd/BBarBrickUP.h"
#include "UP-ucsd/BBarFourNodeQuadUP.h"
#include "UP-ucsd/Twenty_Eight_Node_BrickUP.h"
#include "UP-ucsd/FourNodeQuadUP.h"

#include "dispBeamColumn/DispBeamColumn2d.h"
#include "dispBeamColumn/DispBeamColumn3d.h"
#include "mixedBeamColumn/MixedBeamColumn2d.h"
#include "mixedBeamColumn/MixedBeamColumn3d.h"
#include "dispBeamColumn/DispBeamColumnNL2d.h"
#include "dispBeamColumn/DispBeamColumnNL3d.h"
#include "dispBeamColumn/DispBeamColumnAsym3d.h"    //Xinlong Du
#include "mixedBeamColumn/MixedBeamColumnAsym3d.h"  //Xinlong Du
#include "dispBeamColumn/DispBeamColumn2dThermal.h"
#include "dispBeamColumn/DispBeamColumn3dThermal.h"
#include "shell/ShellMITC4.h"
#include "shell/ShellMITC9.h"
#include "shell/ShellDKGQ.h"   //Added by Lisha Wang, Xinzheng Lu, Linlin Xie, Song Cen & Quan Gu
#include "shell/ShellNLDKGQ.h" //Added by Lisha Wang, Xinzheng Lu, Linlin Xie, Song Cen & Quan Gu
#include "shell/ShellDKGT.h"
#include "shell/ShellNLDKGT.h"
#include "shell/ASDShellQ4.h" // Massimo Petracca
#include "shell/ASDShellT3.h" // Massimo Petracca
#include "brick/Brick.h"
#include "brick/BbarBrick.h"
#include "joint/Joint2D.h"		// Arash
#include "joint/Inno3DPnPJoint.h" // Cristian Miculas
#include "twoNodeLink/TwoNodeLink.h"
#include "twoNodeLink/TwoNodeLinkSection.h"
#include "twoNodeLink/LinearElasticSpring.h"
#include "twoNodeLink/Inerter.h"
#include "tetrahedron/FourNodeTetrahedron.h"

#include "mvlem/MVLEM.h"		// Kristijan Kolozvari
#include "mvlem/SFI_MVLEM.h"	// Kristijan Kolozvari
#include "mvlem/MVLEM_3D.h"		// Kristijan Kolozvari
#include "mvlem/SFI_MVLEM_3D.h"	// Kristijan Kolozvari
#include "mvlem/E_SFI_MVLEM_3D.h"	// Kristijan Kolozvari
#include "mvlem/E_SFI.h"		// C. N. Lopez

#include "mefi/MEFI.h"		// C. N. Lopez

#include "elastomericBearing/ElastomericBearingBoucWen2d.h"
#include "elastomericBearing/ElastomericBearingBoucWen3d.h"
#include "elastomericBearing/ElastomericBearingPlasticity2d.h"
#include "elastomericBearing/ElastomericBearingPlasticity3d.h"
#include "elastomericBearing/ElastomericBearingUFRP2d.h"
#include "elastomericBearing/ElastomericX.h"
#include "elastomericBearing/HDR.h"
#include "elastomericBearing/LeadRubberX.h"

#include "frictionBearing/FlatSliderSimple2d.h"
#include "frictionBearing/FlatSliderSimple3d.h"
#include "frictionBearing/FPBearingPTV.h"
#include "frictionBearing/RJWatsonEQS2d.h"
#include "frictionBearing/RJWatsonEQS3d.h"
#include "frictionBearing/SingleFPSimple2d.h"
#include "frictionBearing/SingleFPSimple3d.h"
#include "frictionBearing/TripleFrictionPendulum.h"
#include "frictionBearing/TripleFrictionPendulumX.h"

#include "PFEMElement/PFEMElement2D.h"
#include "RockingBC/RockingBC.h"

#include "CEqElement/ASDEmbeddedNodeElement.h"

#include "absorbentBoundaries/ASDAbsorbingBoundary2D.h"
#include "absorbentBoundaries/ASDAbsorbingBoundary3D.h"
#include "absorbentBoundaries/FSIFluidElement2D.h"
#include "absorbentBoundaries/FSIInterfaceElement2D.h"
#include "absorbentBoundaries/FSIFluidBoundaryElement2D.h"

#include "LinearCrdTransf2d.h"
#include "LinearCrdTransf3d.h"
#include "PDeltaCrdTransf2d.h"
#include "PDeltaCrdTransf3d.h"
#include "CorotCrdTransf2d.h"
#include "CorotCrdTransf3d.h"

#include "UniformDamping.h"
#include "SecStifDamping.h"
#include "URDDamping.h"
#include "URDDampingbeta.h"

#include "HingeMidpointBeamIntegration.h"
#include "HingeEndpointBeamIntegration.h"
#include "HingeRadauBeamIntegration.h"
#include "HingeRadauTwoBeamIntegration.h"
#include "UserDefinedHingeIntegration.h"
#include "DistHingeIntegration.h"
#include "RegularizedHingeIntegration.h"

#include "LobattoBeamIntegration.h"
#include "LegendreBeamIntegration.h"
#include "RadauBeamIntegration.h"
#include "NewtonCotesBeamIntegration.h"
#include "TrapezoidalBeamIntegration.h"
#include "UserDefinedBeamIntegration.h"
#include "FixedLocationBeamIntegration.h"
#include "LowOrderBeamIntegration.h"
#include "MidDistanceBeamIntegration.h"
#include "CompositeSimpsonBeamIntegration.h"
#include "SimpsonBeamIntegration.h"
#include "ChebyshevBeamIntegration.h"

#include "ConcentratedPlasticityBeamIntegration.h"
#include "ConcentratedCurvatureBeamIntegration.h"

#include "RCCircularSectionIntegration.h"
#include "RCSectionIntegration.h"
#include "RCTBeamSectionIntegration.h"
#include "RCTunnelSectionIntegration.h"
#include "TubeSectionIntegration.h"
#include "WideFlangeSectionIntegration.h"
#include "HSSSectionIntegration.h"

// node header files
#include "Node.h"


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
#ifdef _HDF5
#include "MPCORecorder.h"
#endif // _HDF5
#include "VTK_Recorder.h"
#include "GmshRecorder.h"

// mp_constraint header files
#include "MP_Constraint.h"
#include "joint/MP_Joint2D.h"

// sp_constraint header files
#include "SP_Constraint.h"
#include "SP_Constraint.h"
#include "ImposedMotionSP.h"
#include "ImposedMotionSP1.h"

// Pressure_Constraint header file
#include "Pressure_Constraint.h"

#include "EQ_Constraint.h"

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
#include "AutoConstraintHandler.h"

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
#include "StagedLoadControl.h"

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
#include "StagedNewmark.h"
#include "NewmarkExplicit.h"
#include "NewmarkHSFixedNumIter.h"
#include "NewmarkHSIncrLimit.h"
#include "NewmarkHSIncrReduct.h"
#include "PFEMIntegrator.h"
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
#include "drm/DRMLoadPatternWrapper.h"

#ifdef _H5DRM
#include "drm/H5DRMLoadPattern.h"
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
#include "MPAccSeries.h"   //Tang.S

// time series integrators
#include "TrapezoidalTimeSeriesIntegrator.h"

#include "eigenSOE/ArpackSOE.h"

#ifdef _PETSC
#include "PetscSOE.h"
#include "PetscSolver.h"
#include "SparseGenColLinSOE.h"
#endif


#ifdef _MUMPS
#include "MumpsSOE.h"
#ifdef _PARALLEL_PROCESSING
#include "MumpsParallelSOE.h"
#endif
#endif

#ifdef _PARALLEL_PROCESSING
#include "DistributedBandSPDLinSOE.h"
#include "DistributedProfileSPDLinSOE.h"
#include "DistributedSparseGenColLinSOE.h"
#include "DistributedSparseGenRowLinSOE.h"
#include "DistributedBandGenLinSOE.h"
#include "DistributedSuperLU.h"
#include "ParallelNumberer.h"
#include "StaticDomainDecompositionAnalysis.h"
#include "TransientDomainDecompositionAnalysis.h"
#include "DistributedDiagonalSOE.h"
#endif

//#include "TclFeViewer.h"

#include "packages.h"

typedef struct uniaxialPackage {
  int classTag;
  char *libName;
  char *funcName;
  UniaxialMaterial *(*funcPtr)(void);
  struct uniaxialPackage *next;
} UniaxialPackage;

static UniaxialPackage *theUniaxialPackage = NULL;



FEM_ObjectBrokerAllClasses::FEM_ObjectBrokerAllClasses()
:lastDomainSolver(0)
{

}


FEM_ObjectBrokerAllClasses::~FEM_ObjectBrokerAllClasses()
{

}


Actor *
FEM_ObjectBrokerAllClasses::getNewActor(int classTag, Channel *theChannel)
{
  switch(classTag) {

#ifdef _PARALLEL_PROCESSING
  case ACTOR_TAGS_SUBDOMAIN:  
    return new ActorSubdomain(*theChannel, *this);
#endif

  default:
    opserr << "FEM_ObjectBrokerAllClasses::getNewActor - ";
    opserr << " - no ActorType type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}


PartitionedModelBuilder          *
FEM_ObjectBrokerAllClasses::getPtrNewPartitionedModelBuilder(Subdomain &theSubdomain,
						   int classTag)
{
    switch(classTag) {
	/*
	case PartitionedModelBuilder_TAGS_PartitionedQuick2dFrameModel:  
	     return new PartitionedQuick2dFrame(theSubdomain);
	     */

	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getPtrNewPartitionedModelBuilder - ";
	     opserr << " - no PartitionedModelBuilder type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }    
}


GraphNumberer *
FEM_ObjectBrokerAllClasses::getPtrNewGraphNumberer(int classTag)
{
    switch(classTag) {
	case GraphNUMBERER_TAG_RCM:  
 	     return new RCM();
	     
	     
	case GraphNUMBERER_TAG_MyRCM:  
	     return new MyRCM();
	     	     
	     
	case GraphNUMBERER_TAG_SimpleNumberer:  
	     return new SimpleNumberer();				
	     
	     
	default:
	     opserr << "ObjectBrokerAllClasses::getPtrNewGraphNumberer - ";
	     opserr << " - no GraphNumberer type exists for class tag " ;
	     opserr << classTag << endln;
	     return 0;
	     
	 }
}

/*****************************************
 *
 * METHODS TO GET NEW MODELLING CLASSES
 *
 *****************************************/



Element       *
FEM_ObjectBrokerAllClasses::getNewElement(int classTag)
{
    switch(classTag) {
	     
    case ELE_TAG_Truss:  
      return new Truss(); 
      
    case ELE_TAG_Truss2:  
      return new Truss2(); 
      
    case ELE_TAG_TrussSection:  
      return new TrussSection(); 	     
      
    case ELE_TAG_CorotTruss:  
      return new CorotTruss(); 
      
    case ELE_TAG_CorotTrussSection:  
      return new CorotTrussSection(); 

	case ELE_TAG_InertiaTruss:
		return new InertiaTruss();
      
    case ELE_TAG_ZeroLength:  
      return new ZeroLength();

    case ELE_TAG_CoupledZeroLength:
      return new CoupledZeroLength();
      
    case ELE_TAG_ZeroLengthSection:  
      return new ZeroLengthSection(); 	     
      
    case ELE_TAG_ZeroLengthContact2D:  
      return new ZeroLengthContact2D(); 	     
      
    case ELE_TAG_ZeroLengthContact3D:  
      return new ZeroLengthContact3D(); 	     
      
    case ELE_TAG_ZeroLengthInterface2D:  
      return new ZeroLengthInterface2D(); 	     
      
    case ELE_TAG_ZeroLengthContactNTS2D:  
      return new ZeroLengthContactNTS2D(); 	     
    
    case ELE_TAG_ZeroLengthContactASDimplex:
      return new ZeroLengthContactASDimplex();
      
      //case ELE_TAG_ZeroLengthND:  
      //return new ZeroLengthND(); 	     
      
    case ELE_TAG_FourNodeQuadUP:  
      return new FourNodeQuadUP(); 	     
      
    case ELE_TAG_FourNodeQuad:  
      return new FourNodeQuad(); 

    case ELE_TAG_FourNodeQuad3d:  
      return new FourNodeQuad3d(); 
      
    case ELE_TAG_Tri31:  
      return new Tri31();

    case ELE_TAG_SixNodeTri:  
      return new SixNodeTri();      
      
    case ELE_TAG_ElasticBeam2d:
      return new ElasticBeam2d();
      
	  //SAJalali
	case ELE_TAG_ModElasticBeam2d:
		return new ModElasticBeam2d();

	case ELE_TAG_ModElasticBeam3d:
		return new ModElasticBeam3d();

	case ELE_TAG_ElasticBeam3d:
      return new ElasticBeam3d();

    case ELE_TAG_ComponentElement2d:
      return new ComponentElement2d();

    case ELE_TAG_ComponentElement3d:
      return new ComponentElement3d();      
      
    case ELE_TAG_ElasticTimoshenkoBeam2d:
      return new ElasticTimoshenkoBeam2d();
      
    case ELE_TAG_ElasticTimoshenkoBeam3d:
      return new ElasticTimoshenkoBeam3d();

    case ELE_TAG_TimoshenkoBeamColumn2d:
      return new TimoshenkoBeamColumn2d();
      
    case ELE_TAG_TimoshenkoBeamColumn3d:
      return new TimoshenkoBeamColumn3d();

    case ELE_TAG_GradientInelasticBeamColumn2d:
      return new GradientInelasticBeamColumn2d();
      
    case ELE_TAG_GradientInelasticBeamColumn3d:
      return new GradientInelasticBeamColumn3d();
      
    case ELE_TAG_ForceBeamColumn2d:  
      return new ForceBeamColumn2d();					     
      
    case ELE_TAG_ForceBeamColumn3d:  
      return new ForceBeamColumn3d();  
      
    case ELE_TAG_DispBeamColumn2d:  
      return new DispBeamColumn2d();					     
      
    case ELE_TAG_DispBeamColumn3d:  
      return new DispBeamColumn3d();

    case ELE_TAG_MixedBeamColumn2d:  
      return new MixedBeamColumn2d();					     
      
    case ELE_TAG_MixedBeamColumn3d:  
      return new MixedBeamColumn3d();      

    case ELE_TAG_DispBeamColumnNL2d:  
      return new DispBeamColumnNL2d();					     
      
    case ELE_TAG_DispBeamColumnNL3d:  
      return new DispBeamColumnNL3d();       

	case ELE_TAG_DispBeamColumnAsym3d:
		return new DispBeamColumnAsym3d();    //Xinlong Du

	case ELE_TAG_MixedBeamColumnAsym3d:
		return new MixedBeamColumnAsym3d();   //Xinlong Du

    case ELE_TAG_DispBeamColumn2dThermal:  
      return new DispBeamColumn2dThermal();					     
      
    case ELE_TAG_DispBeamColumn3dThermal:  
      return new DispBeamColumn3dThermal();
      
    case ELE_TAG_EnhancedQuad:
      return new EnhancedQuad();
      
    case ELE_TAG_NineNodeMixedQuad:
      return new NineNodeMixedQuad();
      
    case ELE_TAG_NineNodeQuad:
      return new NineNodeQuad();
      
    case ELE_TAG_EightNodeQuad:
      return new EightNodeQuad();
      
    case ELE_TAG_ConstantPressureVolumeQuad:
      return new ConstantPressureVolumeQuad();
      
    case ELE_TAG_Brick:
      return new Brick();
      
    case ELE_TAG_SSPquad:          
      return new SSPquad();
      
    case ELE_TAG_SSPquadUP:     
      return new SSPquadUP;
      
    case ELE_TAG_SSPbrick:  
      return new SSPbrick();
      
    case ELE_TAG_SSPbrickUP:
      return new SSPbrickUP();

    case ELE_TAG_SurfaceLoad:
      return new SurfaceLoad();

    case ELE_TAG_TriSurfaceLoad:
      return new TriSurfaceLoad();      
      
    case ELE_TAG_Quad4FiberOverlay:
      return new Quad4FiberOverlay(); //Amin Pakzad
	
	case ELE_TAG_Brick8FiberOverlay:
      return new Brick8FiberOverlay(); //Amin Pakzad

	case ELE_TAG_EmbeddedBeamInterfaceL:
	  return new EmbeddedBeamInterfaceL(); //Amin Pakzad

    case ELE_TAG_FourNodeTetrahedron:
      return new FourNodeTetrahedron();
      
	case ELE_TAG_PML2D:
	  return new PML2D();

	case ELE_TAG_PML3D:
	  return new PML3D();

	case ELE_TAG_PML2D_3:
	  return new PML2D_3(); // Amin Pakzad
	
	case ELE_TAG_PML2D_5:
	  return new PML2D_5(); // Amin Pakzad

	case ELE_TAG_PML2D_12:
	  return new PML2D_12(); // Amin Pakzad

	case ELE_TAG_PML2DVISCOUS:
	  return new PML2DVISCOUS(); // Amin Pakzad
	
    case ELE_TAG_BeamContact2D:
      return new BeamContact2D();
      
    case ELE_TAG_BeamContact2Dp:
      return new BeamContact2Dp();
      
    case ELE_TAG_BeamContact3D:
      return new BeamContact3D();
      
    case ELE_TAG_BeamContact3Dp:
      return new BeamContact3Dp();
      
    case ELE_TAG_BeamEndContact3D:
      return new BeamEndContact3D();
      
    case ELE_TAG_BeamEndContact3Dp:
      return new BeamEndContact3Dp();
	  
    case ELE_TAG_QuadBeamEmbedContact:
      return new QuadBeamEmbedContact();
      
    case ELE_TAG_ShellMITC4:
      return new ShellMITC4();

    case ELE_TAG_ShellMITC9:
      return new ShellMITC9();

    case ELE_TAG_ShellDKGQ:      //Added by Lisha Wang, Xinzheng Lu, Linlin Xie, Song Cen & Quan Gu
      return new ShellDKGQ();  //Added by Lisha Wang, Xinzheng Lu, Linlin Xie, Song Cen & Quan Gu
      
    case ELE_TAG_ShellNLDKGQ:      //Added by Lisha Wang, Xinzheng Lu, Linlin Xie, Song Cen & Quan Gu
      return new ShellNLDKGQ();  //Added by Lisha Wang, Xinzheng Lu, Linlin Xie, Song Cen & Quan Gu

    case ELE_TAG_ShellDKGT:
      return new ShellDKGT();
      
    case ELE_TAG_ShellNLDKGT:
      return new ShellNLDKGT();

    case ELE_TAG_ASDShellQ4:   // Massimo Petracca
      return new ASDShellQ4(); // Massimo Petracca
    
    case ELE_TAG_ASDShellT3:   // Massimo Petracca
      return new ASDShellT3(); // Massimo Petracca
	    
    case ELE_TAG_BbarBrick:
      return new BbarBrick();
            
    case ELE_TAG_Joint2D:				// Arash
      return new Joint2D();			// Arash
	  
    case ELE_TAG_Inno3DPnPJoint:	// Cristian Miculas
      return new Inno3DPnPJoint();	// Cristian Miculas
      
    case ELE_TAG_TwoNodeLink:				
      return new TwoNodeLink();

    case ELE_TAG_TwoNodeLinkSection:				
      return new TwoNodeLinkSection();			      
      
    case ELE_TAG_LinearElasticSpring:
        return new LinearElasticSpring();

    case ELE_TAG_Inerter:
        return new Inerter();

	case ELE_TAG_MVLEM:				// Kristijan Kolozvari
		return new MVLEM();	// Kristijan Kolozvari

	case ELE_TAG_SFI_MVLEM:			// Kristijan Kolozvari
		return new SFI_MVLEM();	// Kristijan Kolozvari

	case ELE_TAG_MVLEM_3D:		// Kristijan Kolozvari
		return new MVLEM_3D();	// Kristijan Kolozvari

	case ELE_TAG_SFI_MVLEM_3D:		// Kristijan Kolozvari
		return new SFI_MVLEM_3D();	// Kristijan Kolozvari

	case ELE_TAG_E_SFI_MVLEM_3D:		// Kristijan Kolozvari
		return new E_SFI_MVLEM_3D();	// Kristijan Kolozvari
		
	case ELE_TAG_E_SFI:			// C. N. Lopez
		return new E_SFI();		// C. N. Lopez	
		
	case ELE_TAG_MEFI:			// C. N. Lopez
		return new MEFI();		// C. N. Lopez		

    case ELE_TAG_BBarFourNodeQuadUP:
      return new BBarFourNodeQuadUP();			
      
    case ELE_TAG_BBarBrickUP:
      return new BBarBrickUP();			
      
    case ELE_TAG_Nine_Four_Node_QuadUP:
      return new NineFourNodeQuadUP();
      
    case ELE_TAG_BrickUP:
      return new BrickUP();
      
    case ELE_TAG_Twenty_Eight_Node_BrickUP:
      return new TwentyEightNodeBrickUP();
      
    case ELE_TAG_ElastomericBearingBoucWen2d:
      return new ElastomericBearingBoucWen2d();
      
    case ELE_TAG_ElastomericBearingBoucWen3d:
      return new ElastomericBearingBoucWen3d();
      
    case ELE_TAG_ElastomericBearingPlasticity2d:
      return new ElastomericBearingPlasticity2d();
      
    case ELE_TAG_ElastomericBearingPlasticity3d:
      return new ElastomericBearingPlasticity3d();
      
    case ELE_TAG_ElastomericBearingUFRP2d:
      return new ElastomericBearingUFRP2d();
      
    case ELE_TAG_ElastomericX:
      return new ElastomericX();
      
    case ELE_TAG_HDR:
      return new HDR();
      
    case ELE_TAG_LeadRubberX:
      return new LeadRubberX();
      
    case ELE_TAG_FlatSliderSimple2d:
      return new FlatSliderSimple2d();
      
    case ELE_TAG_FlatSliderSimple3d:
      return new FlatSliderSimple3d();
      
    case ELE_TAG_FPBearingPTV:
      return new FPBearingPTV();
      
    case ELE_TAG_RJWatsonEQS2d:
      return new RJWatsonEQS2d();
      
    case ELE_TAG_RJWatsonEQS3d:
      return new RJWatsonEQS3d();
      
    case ELE_TAG_SingleFPSimple2d:
      return new SingleFPSimple2d();
      
    case ELE_TAG_SingleFPSimple3d:
      return new SingleFPSimple3d();
      
    case ELE_TAG_TripleFrictionPendulum:
      return new TripleFrictionPendulum();

    case ELE_TAG_TripleFrictionPendulumX:
      return new TripleFrictionPendulumX();
		    
    case ELE_TAG_PFEMElement2D:
      return new PFEMElement2D();

    case ELE_TAG_RockingBC:
      return new RockingBC();

    case ELE_TAG_ASDEmbeddedNodeElement:
      return new ASDEmbeddedNodeElement();

    case ELE_TAG_ASDAbsorbingBoundary2D:
      return new ASDAbsorbingBoundary2D();

    case ELE_TAG_ASDAbsorbingBoundary3D:
      return new ASDAbsorbingBoundary3D();

    case ELE_TAG_FSIFluidElement2D:
      return new FSIFluidElement2D();

    case ELE_TAG_FSIInterfaceElement2D:
      return new FSIInterfaceElement2D();

    case ELE_TAG_FSIFluidBoundaryElement2D:
      return new FSIFluidBoundaryElement2D();

    default:
      opserr << "FEM_ObjectBrokerAllClasses::getNewElement - ";
      opserr << " - no Element type exists for class tag " ;
      opserr << classTag << endln;
      return 0;
      
    }
}

Node          *
FEM_ObjectBrokerAllClasses::getNewNode(int classTag)
{
    switch(classTag) {
	case NOD_TAG_Node:  
	     return new Node(classTag);
	     
	     
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewNode - ";
	     opserr << " - no Node type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }    
}


MP_Constraint *
FEM_ObjectBrokerAllClasses::getNewMP(int classTag)
{
    switch(classTag) {
	case CNSTRNT_TAG_MP_Constraint:  
	     return new MP_Constraint(classTag);

 	case CNSTRNT_TAG_MP_Joint2D:			// Arash
	     return new MP_Joint2D();			// Arash
	
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewMP - ";
	     opserr << " - no MP_Constraint type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }    
}


EQ_Constraint *
FEM_ObjectBrokerAllClasses::getNewEQ(int classTag)
{
    switch(classTag) {
	case CNSTRNT_TAG_EQ_Constraint:  
	     return new EQ_Constraint(classTag);

	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewEQ - ";
	     opserr << " - no EQ_Constraint type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }    
}


SP_Constraint *
FEM_ObjectBrokerAllClasses::getNewSP(int classTag)
{
    switch(classTag) {
	case CNSTRNT_TAG_SP_Constraint:  
	     return new SP_Constraint(classTag);

	case CNSTRNT_TAG_ImposedMotionSP:  
	     return new ImposedMotionSP();

	case CNSTRNT_TAG_ImposedMotionSP1:  
	     return new ImposedMotionSP1();
	     
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewSP - ";
	     opserr << " - no SP_Constraint type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }    
}

Pressure_Constraint *
FEM_ObjectBrokerAllClasses::getNewPC(int classTag)
{
    switch(classTag) {
    case CNSTRNT_TAG_Pressure_Constraint:  
        return new Pressure_Constraint(classTag);
	
    default:
        opserr << "FEM_ObjectBrokerAllClasses::getNewPC - ";
        opserr << " - no Pressure_Constraint type exists for class tag ";
        opserr << classTag << endln;
        return 0;
	
    }    
}

NodalLoad     *
FEM_ObjectBrokerAllClasses::getNewNodalLoad(int classTag)
{
    switch(classTag) {
	case LOAD_TAG_NodalLoad:  
	     return new NodalLoad(classTag);
	     
	     
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewNodalLoad - ";
	     opserr << " - no NodalLoad type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }    
}


ElementalLoad *
FEM_ObjectBrokerAllClasses::getNewElementalLoad(int classTag)
{
  switch(classTag) {
    
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
    opserr << "FEM_ObjectBrokerAllClasses::getNewNodalLoad - ";
    opserr << " - no NodalLoad type exists for class tag ";
    opserr << classTag << endln;
    return 0;
    
  }    
  
  return 0;
}

CrdTransf*
FEM_ObjectBrokerAllClasses::getNewCrdTransf(int classTag)
{
	switch(classTag) {
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
	  opserr << "FEM_ObjectBrokerAllClasses::getCrdTransf - ";
	  opserr << " - no CrdTransf type exists for class tag ";
	  opserr << classTag << endln;
	  return 0;
	}

}

Damping*
FEM_ObjectBrokerAllClasses::getNewDamping(int classTag)
{
	switch(classTag) {
	case DMP_TAG_UniformDamping:
		return new UniformDamping();
	case DMP_TAG_SecStifDamping:
		return new SecStifDamping();
	case DMP_TAG_URDDamping:
		return new URDDamping();
	case DMP_TAG_URDDampingbeta:
		return new URDDampingbeta();
	default:
	  opserr << "FEM_ObjectBrokerAllClasses::getNewDamping - ";
	  opserr << " - no Damping type exists for class tag ";
	  opserr << classTag << endln;
	  return 0;
	}

}

BeamIntegration *
FEM_ObjectBrokerAllClasses::getNewBeamIntegration(int classTag)
{
  switch(classTag) {
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

  case BEAM_INTEGRATION_TAG_Simpson:        
    return new SimpsonBeamIntegration();

  case BEAM_INTEGRATION_TAG_Chebyshev:
    return new ChebyshevBeamIntegration();

  case BEAM_INTEGRATION_TAG_ConcentratedPlasticity:
	  return new ConcentratedPlasticityBeamIntegration();

  case BEAM_INTEGRATION_TAG_ConcentratedCurvature:
	  return new ConcentratedCurvatureBeamIntegration();

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
    opserr << "FEM_ObjectBrokerAllClasses::getBeamIntegration - ";
    opserr << " - no BeamIntegration type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

SectionIntegration *
FEM_ObjectBrokerAllClasses::getNewSectionIntegration(int classTag)
{
  switch(classTag) {
  case SECTION_INTEGRATION_TAG_WideFlange:        
    return new WideFlangeSectionIntegration();

  case SECTION_INTEGRATION_TAG_RC:        
    return new RCSectionIntegration();

  case SECTION_INTEGRATION_TAG_RCT:        
    return new RCTBeamSectionIntegration();        

  case SECTION_INTEGRATION_TAG_RCCIRCULAR:        
    return new RCCircularSectionIntegration();

  case SECTION_INTEGRATION_TAG_RCTUNNEL:        
    return new RCTunnelSectionIntegration();

  case SECTION_INTEGRATION_TAG_Tube:        
    return new TubeSectionIntegration();    

  case SECTION_INTEGRATION_TAG_HSS:        
    return new HSSSectionIntegration();
	  
  default:
    opserr << "FEM_ObjectBrokerAllClasses::getSectionIntegration - ";
    opserr << " - no SectionIntegration type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}


HystereticBackbone *
FEM_ObjectBrokerAllClasses::getNewHystereticBackbone(int classTag)
{
  switch (classTag) {
  case BACKBONE_TAG_Capped:
    return new CappedBackbone();

  case BACKBONE_TAG_LinearCapped:
    return new LinearCappedBackbone();    

  case BACKBONE_TAG_Material:
    return new MaterialBackbone();
    
  case BACKBONE_TAG_Arctangent:
    return new ArctangentBackbone();

  case BACKBONE_TAG_Trilinear:
    return new TrilinearBackbone();

  case BACKBONE_TAG_Multilinear:
    return new MultilinearBackbone();

  case BACKBONE_TAG_Mander:
    return new ManderBackbone();

    //  case BACKBONE_TAG_KentPark:
    //return new KentParkBackbone();

  case BACKBONE_TAG_Raynor:
    return new RaynorBackbone();

  case BACKBONE_TAG_ReeseStiffClayBelowWS:
    return new ReeseStiffClayBelowWS();

  case BACKBONE_TAG_ReeseStiffClayAboveWS:
    return new ReeseStiffClayAboveWS();

  case BACKBONE_TAG_ReeseSoftClay:
    return new ReeseSoftClayBackbone();

  case BACKBONE_TAG_ReeseSand:
    return new ReeseSandBackbone();

  case BACKBONE_TAG_VuggyLimestone:
    return new VuggyLimestone();

  case BACKBONE_TAG_CementedSoil:
    return new CementedSoil();

  case BACKBONE_TAG_WeakRock:
    return new WeakRock();

  case BACKBONE_TAG_LiquefiedSand:
    return new LiquefiedSand();
    
  default:
    opserr << "FEM_ObjectBrokerAllClasses::getNewHystereticBackbone - ";
    opserr << " - no HystereticBackbone type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

UnloadingRule *
FEM_ObjectBrokerAllClasses::getNewUnloadingRule(int classTag)
{
  switch (classTag) {
  case DEG_TAG_UNLOAD_Constant:
    return new ConstantUnloadingRule();

  case DEG_TAG_UNLOAD_Takeda:
    return new TakedaUnloadingRule();

  case DEG_TAG_UNLOAD_Energy:
    return new EnergyUnloadingRule();

  case DEG_TAG_UNLOAD_Karsan:
    return new KarsanUnloadingRule();            
    
  default:
    opserr << "FEM_ObjectBrokerAllClasses::getNewUnloadingRule - ";
    opserr << " - no UnloadingRule type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

StiffnessDegradation *
FEM_ObjectBrokerAllClasses::getNewStiffnessDegradation(int classTag)
{
  switch (classTag) {
  case DEG_TAG_STIFF_Constant:
    return new ConstantStiffnessDegradation();

  case DEG_TAG_STIFF_Ductility:
    return new DuctilityStiffnessDegradation();

  case DEG_TAG_STIFF_Energy:
    return new EnergyStiffnessDegradation();

  case DEG_TAG_STIFF_Pincheira:
    return new PincheiraStiffnessDegradation();            
    
  default:
    opserr << "FEM_ObjectBrokerAllClasses::getStiffnessDegradation - ";
    opserr << " - no StiffnessDegradation type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

StrengthDegradation *
FEM_ObjectBrokerAllClasses::getNewStrengthDegradation(int classTag)
{
  switch (classTag) {
  case DEG_TAG_STRENGTH_Constant:
    return new ConstantStrengthDegradation();

  case DEG_TAG_STRENGTH_Ductility:
    return new DuctilityStrengthDegradation();

  case DEG_TAG_STRENGTH_Petrangeli:
    return new PetrangeliStrengthDegradation();

  case DEG_TAG_STRENGTH_Energy:
    return new EnergyStrengthDegradation();

  case DEG_TAG_STRENGTH_Section:
    return new SectionStrengthDegradation();

  case DEG_TAG_STRENGTH_ACI:
    return new ACIStrengthDegradation();    
    
  default:
    opserr << "FEM_ObjectBrokerAllClasses::getStrengthDegradation - ";
    opserr << " - no StrengthDegradation type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
}

UniaxialMaterial *
FEM_ObjectBrokerAllClasses::getNewUniaxialMaterial(int classTag)
{
    switch(classTag) {
    case MAT_TAG_ElasticBilin:
      return new ElasticBilin();
      
	case MAT_TAG_SPSW02:
		return new SPSW02(); // SAJalali
	case MAT_TAG_BoucWen:
		return new BoucWenMaterial();
	case MAT_TAG_BoucWenOriginal:
		return new BoucWenOriginal();
	case MAT_TAG_BWBN:
		return new BWBN();		
	case MAT_TAG_BoucWenInfill:
		return new BoucWenInfill();		
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

    case MAT_TAG_DamperMaterial:
      return new DamperMaterial();

    case MAT_TAG_Penalty:
      return new PenaltyMaterial();

    case MAT_TAG_Multiplier:
      return new MultiplierMaterial();

    case MAT_TAG_TensionOnly:
      return new TensionOnlyMaterial();      

	case MAT_TAG_ASD_SMA_3K:  
	     return new ASD_SMA_3K();

	case MAT_TAG_ASDConcrete1DMaterial:  
	     return new ASDConcrete1DMaterial();

	case MAT_TAG_Concrete01:  
	     return new Concrete01();

	case MAT_TAG_Concrete01WithSITC:  
	     return new Concrete01WithSITC();	     

	case MAT_TAG_Concrete02:  
	     return new Concrete02();

	case MAT_TAG_Concrete02IS:  
	     return new Concrete02IS();	     

	case MAT_TAG_Concrete04:  
	     return new Concrete04();

	case MAT_TAG_Concrete06:  
	     return new Concrete06();

	case MAT_TAG_Concrete07:  
	     return new Concrete07();

	case MAT_TAG_ConcretewBeta:  
	     return new ConcretewBeta();

    case MAT_TAG_CreepMaterial:
      return new CreepMaterial();

    case MAT_TAG_TDConcrete:
      return new TDConcrete();

    case MAT_TAG_TDConcreteNL:
      return new TDConcreteNL();

    case MAT_TAG_TDConcreteEXP:
      return new TDConcreteEXP();

    case MAT_TAG_TDConcreteMC10:
      return new TDConcreteMC10();

    case MAT_TAG_TDConcreteMC10NL:
      return new TDConcreteMC10NL();
      
	case MAT_TAG_Steel01:  
	     return new Steel01();

	case MAT_TAG_Steel02:  
	     return new Steel02();

	case MAT_TAG_SteelMPF:  
	     return new SteelMPF();	     

	case MAT_TAG_Steel2:  
	     return new Steel2();

	case MAT_TAG_Steel4:  
	     return new Steel4();

	case MAT_TAG_RambergOsgoodSteel:  
	     return new RambergOsgoodSteel();	     	     

	case MAT_TAG_SteelFractureDI:
		return new SteelFractureDI();

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

	case MAT_TAG_AxialSp:
	    return new AxialSp();

	case MAT_TAG_AxialSpHD:
	    return new AxialSpHD();

	case MAT_TAG_SelfCentering:
	    return new SelfCenteringMaterial();
	    
    case MAT_TAG_TzLiq1:
		return new TzLiq1();

	case MAT_TAG_QzSimple1:
		return new QzSimple1();

	case MAT_TAG_QzSimple2:
		return new QzSimple2();

    case MAT_TAG_QzLiq1:
		return new QzLiq1();

	case MAT_TAG_TzSandCPT:
		return new TzSandCPT();

	case MAT_TAG_QbSandCPT:
		return new QbSandCPT();

	case MAT_TAG_Hysteretic:
		return new HystereticMaterial();

	case MAT_TAG_HystereticSM:			
		return new HystereticSMMaterial();

	case MAT_TAG_OOHysteretic:			
		return new OOHystereticMaterial();		

	case MAT_TAG_ModIMKPeakOriented:
		return new ModIMKPeakOriented();

	case MAT_TAG_IMKBilin:
		return new IMKBilin();

	case MAT_TAG_IMKPeakOriented:
		return new IMKPeakOriented();

    case MAT_TAG_IMKPinching:
      return new IMKPinching();

	case MAT_TAG_UVCuniaxial:
	    return new UVCuniaxial();

	case MAT_TAG_SLModel:
	    return new SLModel();

	case MAT_TAG_SteelDRC:
	    return new SteelDRC();
	    
	case MAT_TAG_Cast:
	    return new Cast();

	case MAT_TAG_SnapClough:
		return new Clough();

	case MAT_TAG_LimitState:
		return new LimitStateMaterial();

	case MAT_TAG_EPPGap:
		return new EPPGapMaterial();

	case MAT_TAG_HookGap:
	    return new HookGap();
	    
	case MAT_TAG_Viscous:
		return new ViscousMaterial();

	case MAT_TAG_ViscousDamper:
		return new ViscousDamper();

        case MAT_TAG_APDVFD:
		return new APDVFD();
	    
	case MAT_TAG_APDMD:
		return new APDMD();
	    
	case MAT_TAG_APDFMD:
		return new APDFMD();

	case MAT_TAG_BilinearOilDamper:
	    return new BilinearOilDamper();

	case MAT_TAG_Maxwell:
	    return new Maxwell();
	    
	case MAT_TAG_ContinuumUniaxial:
		return new ContinuumUniaxial();

	case MAT_TAG_PathIndependent:
		return new PathIndependentMaterial();

	case MAT_TAG_Backbone:
		return new BackboneMaterial();		

	case MAT_TAG_SeriesMaterial:
		return new SeriesMaterial();

	case MAT_TAG_CableMaterial:
		return new CableMaterial();
	     
	case MAT_TAG_ENTMaterial:
		return new ENTMaterial();

	case MAT_TAG_GNG:
		return new GNGMaterial();
		
	case MAT_TAG_Ratchet:
		return new Ratchet();				
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

	case MAT_TAG_Hertzdamp:
		return new Hertzdamp();

	case MAT_TAG_JankowskiImpact:
		return new JankowskiImpact();

	case MAT_TAG_ViscoelasticGap:
		return new ViscoelasticGap();

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

        case MAT_TAG_FRCC:
	  return new FRCC();	  

        case MAT_TAG_SimpleFractureMaterial:
	  return new SimpleFractureMaterial();

        case MAT_TAG_ConfinedConcrete01:
            return new ConfinedConcrete01();

	case MAT_TAG_FRPConfinedConcrete:
	    return new FRPConfinedConcrete();

	case MAT_TAG_ConcreteCM:
	    return new ConcreteCM();
	    
	case MAT_TAG_BarSlip:
	    return new BarSlipMaterial();
		
	case MAT_TAG_HystereticPoly:			// Salvatore Sessa
	    return new HystereticPoly();
		    
	case MAT_TAG_HystereticSmooth:			// Salvatore Sessa
			return new HystereticSmooth();

	case MAT_TAG_HystereticAsym:			// Salvatore Sessa
		return new HystereticAsym();	    

	case MAT_TAG_DowelType:
		return new DowelType();

	case MAT_TAG_DuctileFracture:
		return new DuctileFracture();

	case MAT_TAG_GMG_CyclicReinforcedConcrete:
		return new GMG_CyclicReinforcedConcrete();

	case MAT_TAG_Pinching4:
		return new Pinching4Material();
	
	case MAT_TAG_CFSSSWP:
	    return new CFSSSWP();

	case MAT_TAG_CFSWSWP:
	    return new CFSWSWP();

	case MAT_TAG_SAWSMaterial:
		return new SAWSMaterial();

    case MAT_TAG_KikuchiAikenHDR:
      return new KikuchiAikenHDR();

    case MAT_TAG_KikuchiAikenLRB:
      return new KikuchiAikenLRB();      

	default:

	  UniaxialPackage *matCommands = theUniaxialPackage;
	  bool found = false;
	  while (matCommands != NULL && found == false) {
	    if ((matCommands->classTag == classTag) && (matCommands->funcPtr != 0)){
	      UniaxialMaterial *result = (*(matCommands->funcPtr))();
	      return result;
	    } 
	    matCommands = matCommands->next;
	  }	  

	  opserr << "FEM_ObjectBrokerAllClasses::getNewUniaxialMaterial - ";
	  opserr << " - no UniaxialMaterial type exists for class tag ";
	  opserr << classTag << endln;
	  return 0;
	  
    }        
}

SectionForceDeformation *
FEM_ObjectBrokerAllClasses::getNewSection(int classTag)
{
    switch(classTag) {
	case SEC_TAG_Elastic2d:
	     return new ElasticSection2d();
	     
	case SEC_TAG_Elastic3d:
	     return new ElasticSection3d();	     
	     
    case SEC_TAG_ElasticShear2d:
	     return new ElasticShearSection2d();

    case SEC_TAG_ElasticBDShear2d:
	     return new ElasticBDShearSection2d();
	     
	case SEC_TAG_ElasticShear3d:
	     return new ElasticShearSection3d();

	case SEC_TAG_ElasticTube3d:
	     return new ElasticTubeSection3d();

    case SEC_TAG_ElasticWarpingShear2d:
      return new ElasticWarpingShearSection2d();
	     
	case SEC_TAG_Generic1d:
	     return new GenericSection1d();
	     
	     //case SEC_TAG_GenericNd:
	     //return new GenericSectionNd();	     

	case SEC_TAG_Aggregator:
	     return new SectionAggregator();

	case SEC_TAG_Parallel:
	     return new ParallelSection();	     

	     //case SEC_TAG_Fiber:
	     //return new FiberSection();
	
	case SEC_TAG_FiberSection2d:
		return new FiberSection2d();
      
	case SEC_TAG_FiberSection3d:
		return new FiberSection3d();

	case SEC_TAG_FiberSection2dThermal:
		return new FiberSection2dThermal();
      
	case SEC_TAG_FiberSection3dThermal:
		return new FiberSection3dThermal();		

	case SEC_TAG_NDFiberSection2d:
		return new NDFiberSection2d();

	case SEC_TAG_NDFiberSectionWarping2d:
		return new NDFiberSectionWarping2d();		
      
	case SEC_TAG_NDFiberSection3d:
		return new NDFiberSection3d();		

	case SEC_TAG_FiberSectionAsym3d:
		return new FiberSectionAsym3d(); //Xinlong Du

	case SEC_TAG_FiberSectionWarping3d:
		return new FiberSectionWarping3d();
		
	case SEC_TAG_ElasticPlateSection:
		return new ElasticPlateSection();

	case SEC_TAG_ElasticMembranePlateSection:
		return new ElasticMembranePlateSection();

	case SEC_TAG_MembranePlateFiberSection:
		return new MembranePlateFiberSection();

	case SEC_TAG_DoubleMembranePlateFiberSection:
		return new DoubleMembranePlateFiberSection();

	//start Yuli Huang & Xinzheng Lu LayeredShellFiberSection
        case SEC_TAG_LayeredShellFiberSection:
	  return new LayeredShellFiberSection();
	//end Yuli Huang & Xinzheng Lu LayeredShellFiberSection

	case SEC_TAG_Bidirectional:
		return new Bidirectional();
	case SEC_TAG_Elliptical2:
		return new Elliptical2();
    case SEC_TAG_Isolator2spring:
      return new Isolator2spring();
	

	case SEC_TAG_ReinforcedConcreteLayeredMembraneSection:
		return new ReinforcedConcreteLayeredMembraneSection();

	case SEC_TAG_LayeredMembraneSection:
		return new LayeredMembraneSection();

	case SEC_TAG_ElasticMembraneSection:
		return new ElasticMembraneSection();

	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewSection - ";
	     opserr << " - no section type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}

NDMaterial*
FEM_ObjectBrokerAllClasses::getNewNDMaterial(int classTag)
{
  switch(classTag) {
  case ND_TAG_ElasticIsotropicPlaneStrain2d:
    return new ElasticIsotropicPlaneStrain2D();
    
  case ND_TAG_ElasticIsotropicPlaneStress2d:
    return new ElasticIsotropicPlaneStress2D();
		
  case ND_TAG_ElasticIsotropicAxiSymm:
    return new ElasticIsotropicAxiSymm();
    
  case ND_TAG_ElasticIsotropicPlateFiber:
    return new ElasticIsotropicPlateFiber();

  case ND_TAG_ElasticIsotropicBeamFiber:
    return new ElasticIsotropicBeamFiber();

  case ND_TAG_ElasticIsotropicBeamFiber2d:
    return new ElasticIsotropicBeamFiber2d();
    
  case ND_TAG_ElasticIsotropicThreeDimensional:
    return new ElasticIsotropicThreeDimensional();

  case ND_TAG_ElasticOrthotropicThreeDimensional:
    return new ElasticOrthotropicThreeDimensional();

  case ND_TAG_ElasticOrthotropicPlaneStress:
    return new ElasticOrthotropicPlaneStress();    
		  
  case ND_TAG_J2PlaneStrain:
    return new J2PlaneStrain();
    
  case ND_TAG_J2PlaneStress:
    return new J2PlaneStress();
    
  case ND_TAG_J2AxiSymm:
    return new J2AxiSymm();
    
  case ND_TAG_J2PlateFiber:
    return new J2PlateFiber();

  case ND_TAG_J2PlateFibre:
    return new J2PlateFibre();    

  case ND_TAG_J2BeamFiber2d:
    return new J2BeamFiber2d();

  case ND_TAG_J2BeamFiber3d:
    return new J2BeamFiber3d();
    
  case ND_TAG_J2ThreeDimensional:
    return new J2ThreeDimensional();

  case ND_TAG_SimplifiedJ2:
    return new SimplifiedJ2();

  case ND_TAG_PlaneStressSimplifiedJ2:
    return new PlaneStressSimplifiedJ2();
    
  case ND_TAG_PlaneStressMaterial:
    return new PlaneStressMaterial();

  case ND_TAG_PlaneStrainMaterial:
    return new PlaneStrainMaterial();

  case ND_TAG_OrthotropicMaterial:
    return new OrthotropicMaterial();

  case ND_TAG_Series3DMaterial:
    return new Series3DMaterial();

  case ND_TAG_Parallel3DMaterial:
    return new Parallel3DMaterial();

  case ND_TAG_PlaneStressRebarMaterial:
    return new PlaneStressRebarMaterial();

  case ND_TAG_PlaneStressLayeredMaterial:
    return new PlaneStressLayeredMaterial();

  //start Yuli Huang & Xinzheng 
  case ND_TAG_PlateRebarMaterial:
    return new PlateRebarMaterial();

  case ND_TAG_PlateFromPlaneStressMaterial:
    return new PlateFromPlaneStressMaterial();

  case ND_TAG_ConcreteS:
    return new ConcreteS();

  case ND_TAG_PlaneStressUserMaterial:
    return new PlaneStressUserMaterial();
  //end Yuli Huang & Xinzheng Lu 
		  
  case ND_TAG_PlateFiberMaterial:
    return new PlateFiberMaterial();
    
  case ND_TAG_FluidSolidPorousMaterial:
    return new FluidSolidPorousMaterial();

  case ND_TAG_PressureDependMultiYield:
    return new PressureDependMultiYield();

  case ND_TAG_PressureDependMultiYield02:
    return new PressureDependMultiYield02();

  case ND_TAG_PressureDependMultiYield03:
    return new PressureDependMultiYield03();
	
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

  case ND_TAG_InitStrainNDMaterial:
      return new InitStrainNDMaterial();

  case ND_TAG_MinMaxNDMaterial:
      return new MinMaxNDMaterial();      

  case ND_TAG_ASDConcrete3DMaterial:
      return new ASDConcrete3DMaterial();

  case ND_TAG_PlasticDamageConcrete3d:
      return new PlasticDamageConcrete3d();

  case ND_TAG_PlasticDamageConcretePlaneStress:
      return new PlasticDamageConcretePlaneStress();

  case ND_TAG_OrthotropicRotatingAngleConcreteT2DMaterial01:
	  return new OrthotropicRotatingAngleConcreteT2DMaterial01();

  case ND_TAG_SmearedSteelDoubleLayerT2DMaterial01:
	  return new SmearedSteelDoubleLayerT2DMaterial01();

  case ND_TAG_SAniSandMS: 
    return new SAniSandMS();

  case ND_TAG_SAniSandMSPlaneStrain: 
    return new SAniSandMSPlaneStrain();

  case ND_TAG_SAniSandMS3D: 
    return new SAniSandMS3D();

    
  default:
    opserr << "FEM_ObjectBrokerAllClasses::getNewNDMaterial - ";
    opserr << " - no NDMaterial type exists for class tag ";
    opserr << classTag << endln;
    return 0;   
  }
}

Fiber*
FEM_ObjectBrokerAllClasses::getNewFiber(int classTag)
{
	switch(classTag) {
	case FIBER_TAG_Uniaxial2d:
		return new UniaxialFiber2d();

	case FIBER_TAG_Uniaxial3d:
		return new UniaxialFiber3d();

	default:
		opserr << "FEM_ObjectBrokerAllClasses::getNewFiber - ";
		opserr << " - no Fiber type exists for class tag ";
		opserr << classTag << endln;
		return 0;
	}
}

FrictionModel *
FEM_ObjectBrokerAllClasses::getNewFrictionModel(int classTag)
{
    switch(classTag) {
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
	  opserr << "FEM_ObjectBrokerAllClasses::getNewFrictionModel - ";
	  opserr << " - no FrictionModel type exists for class tag ";
	  opserr << classTag << endln;
	  return 0;
    }        
}

ConvergenceTest *
FEM_ObjectBrokerAllClasses::getNewConvergenceTest(int classTag)
{
    switch(classTag) {
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
	     opserr << "FEM_ObjectBrokerAllClasses::getNewConvergenceTest - ";
	     opserr << " - no ConvergenceTest type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }
}


LoadPattern *
FEM_ObjectBrokerAllClasses::getNewLoadPattern(int classTag)
{
    switch(classTag) {
	case PATTERN_TAG_LoadPattern:
	     return new LoadPattern();

	case PATTERN_TAG_UniformExcitation:
	     return new UniformExcitation();

	case PATTERN_TAG_MultiSupportPattern:
	     return new MultiSupportPattern();

	case PATTERN_TAG_DRMLoadPattern:
	     return new DRMLoadPatternWrapper();

#ifdef _H5DRM
    case PATTERN_TAG_H5DRM:
         return new H5DRMLoadPattern();
#endif
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getPtrLoadPattern - ";
	     opserr << " - no Load type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}


GroundMotion *
FEM_ObjectBrokerAllClasses::getNewGroundMotion(int classTag)
{
    switch(classTag) {

        case GROUND_MOTION_TAG_GroundMotion:
	  return new GroundMotion(GROUND_MOTION_TAG_GroundMotion);

        case GROUND_MOTION_TAG_InterpolatedGroundMotion:
	  return new GroundMotion(GROUND_MOTION_TAG_InterpolatedGroundMotion);

	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getPtrGroundMotion - ";
	     opserr << " - no Load type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}

TimeSeries *
FEM_ObjectBrokerAllClasses::getNewTimeSeries(int classTag)
{
    switch(classTag) {
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
	  //Tang.S
		case TSERIES_TAG_MPAccSeries:
	  return new MPAccSeries;

	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getPtrTimeSeries - ";
	     opserr << " - no Load type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}

TimeSeriesIntegrator *
FEM_ObjectBrokerAllClasses::getNewTimeSeriesIntegrator(int classTag)
{
    switch(classTag) {
    case TIMESERIES_INTEGRATOR_TAG_Trapezoidal:
	  return new TrapezoidalTimeSeriesIntegrator();

	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getPtrTimeSeriesIntegrator - ";
	     opserr << " - no Load type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}


Matrix	  *
FEM_ObjectBrokerAllClasses::getPtrNewMatrix(int classTag, int noRows, int noCols)
{
    switch(classTag) {
	case MATRIX_TAG_Matrix:  
	     return new Matrix(noRows,noCols);
	     
	     
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getPtrNewMatrix - ";
	     opserr << " - no NodalLoad type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}


Vector	  *
FEM_ObjectBrokerAllClasses::getPtrNewVector(int classTag, int size)
{
    switch(classTag) {
	case VECTOR_TAG_Vector:  
	     return new Vector(size);
	     
	     
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getPtrNewVector - ";
	     opserr << " - no Vector type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}


ID	          *
FEM_ObjectBrokerAllClasses::getPtrNewID(int classTag, int size)
{
    switch(classTag) {
	case ID_TAG_ID:  
	     return new ID(size);
	     
	     
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getPtrNewID - ";
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
FEM_ObjectBrokerAllClasses::getPtrNewStream(int classTag)
{
    switch(classTag) {
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
      opserr << "FEM_ObjectBrokerAllClasses::getPtrNewStream - ";
      opserr << " - no DataOutputHandler type exists for class tag ";
      opserr << classTag << endln;
      return 0;
	     
	 }        
}

Recorder *
FEM_ObjectBrokerAllClasses::getPtrNewRecorder(int classTag)
{
    switch(classTag) {
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
#ifdef _HDF5
	case RECORDER_TAGS_MPCORecorder:
	  return new MPCORecorder();
#endif // _HDF5
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewRecordr - ";
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

ConstraintHandler   *
FEM_ObjectBrokerAllClasses::getNewConstraintHandler(int classTag)
{
    switch(classTag) {
	case HANDLER_TAG_PlainHandler:  
	     return new PlainHandler();
	     
	case HANDLER_TAG_PenaltyConstraintHandler:  
	     return new PenaltyConstraintHandler(1.0e12, 1.0e12);

	case HANDLER_TAG_LagrangeConstraintHandler:  
	     return new LagrangeConstraintHandler(1.0, 1.0);

	case HANDLER_TAG_TransformationConstraintHandler:  
	     return new TransformationConstraintHandler();

	case HANDLER_TAG_AutoConstraintHandler:  
	     return new AutoConstraintHandler();

	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewConstraintHandler - ";
	     opserr << " - no ConstraintHandler type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}


DOF_Numberer        *
FEM_ObjectBrokerAllClasses::getNewNumberer(int classTag)
{
    switch(classTag) {
	case NUMBERER_TAG_DOF_Numberer:  
	     return new DOF_Numberer();
	     
	     
	case NUMBERER_TAG_PlainNumberer:  
	     return new PlainNumberer();


#ifdef _PARALLEL_PROCESSING
	case NUMBERER_TAG_ParallelNumberer:  
	     return new ParallelNumberer();
#endif
	     
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewConstraintHandler - ";
	     opserr << " - no ConstraintHandler type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }
}


AnalysisModel       *
FEM_ObjectBrokerAllClasses::getNewAnalysisModel(int classTag)
{
    switch(classTag) {
	case AnaMODEL_TAGS_AnalysisModel:  
	     return new AnalysisModel();
	     
	     
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewAnalysisModel - ";
	     opserr << " - no AnalysisModel type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}


EquiSolnAlgo        *
FEM_ObjectBrokerAllClasses::getNewEquiSolnAlgo(int classTag)
{
    switch(classTag) {
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
	     opserr << "FEM_ObjectBrokerAllClasses::getNewEquiSolnAlgo - ";
	     opserr << " - no EquiSolnAlgo type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}

Accelerator        *
FEM_ObjectBrokerAllClasses::getAccelerator(int classTag)
{
    switch(classTag) {

    case ACCELERATOR_TAGS_Krylov:
      return new KrylovAccelerator;
    case ACCELERATOR_TAGS_Raphson:
      return new RaphsonAccelerator;

    default:
      opserr << "FEM_ObjectBrokerAllClasses::getAccelerator - ";
      opserr << " - no EquiSolnAlgo type exists for class tag ";
      opserr << classTag << endln;
      return 0;
      
    }        
}

LineSearch        *
FEM_ObjectBrokerAllClasses::getLineSearch(int classTag)
{
    switch(classTag) {

    case LINESEARCH_TAGS_BisectionLineSearch:
      return new BisectionLineSearch();

    case LINESEARCH_TAGS_InitialInterpolatedLineSearch:
      return new InitialInterpolatedLineSearch();

    case  LINESEARCH_TAGS_RegulaFalsiLineSearch:
      return new RegulaFalsiLineSearch();
    
    case  LINESEARCH_TAGS_SecantLineSearch:
      return new SecantLineSearch();
    default:
      opserr << "FEM_ObjectBrokerAllClasses::getNewEquiSolnAlgo - ";
      opserr << " - no EquiSolnAlgo type exists for class tag ";
      opserr << classTag << endln;
      return 0;
    }        
}


DomainDecompAlgo    *
FEM_ObjectBrokerAllClasses::getNewDomainDecompAlgo(int classTag)
{
    switch(classTag) {
	case DomDecompALGORITHM_TAGS_DomainDecompAlgo:  
	     return new DomainDecompAlgo();

	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewDomainDecompAlgo - ";
	     opserr << " - no DomainDecompAlgo type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }
}


StaticIntegrator    *
FEM_ObjectBrokerAllClasses::getNewStaticIntegrator(int classTag)
{
    switch(classTag) {
	case INTEGRATOR_TAGS_LoadControl:  
	     return new LoadControl(1.0,1,1.0,.10); // must recvSelf

    case INTEGRATOR_TAGS_StagedLoadControl:
        return new StagedLoadControl(1.0, 1, 1.0, .10); // must recvSelf

#ifdef _PARALLEL_PROCESSING
	case INTEGRATOR_TAGS_DistributedDisplacementControl:  
	     return new DistributedDisplacementControl(); // must recvSelf
#endif	     
	     
	case INTEGRATOR_TAGS_ArcLength:  
	     return new ArcLength(1.0);      // must recvSelf

	     
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewStaticIntegrator - ";
	     opserr << " - no StaticIntegrator type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }
}


TransientIntegrator *
FEM_ObjectBrokerAllClasses::getNewTransientIntegrator(int classTag)
{
    switch(classTag) {
	case INTEGRATOR_TAGS_AlphaOS:  
	     return new AlphaOS();

	case INTEGRATOR_TAGS_AlphaOS_TP:  
	     return new AlphaOS_TP();

	case INTEGRATOR_TAGS_AlphaOSGeneralized:  
	     return new AlphaOSGeneralized();

	case INTEGRATOR_TAGS_AlphaOSGeneralized_TP:  
	     return new AlphaOSGeneralized_TP();

	case INTEGRATOR_TAGS_CentralDifference:  
	     return new CentralDifference();      // must recvSelf

	case INTEGRATOR_TAGS_CentralDifferenceAlternative:  
	     return new CentralDifferenceAlternative();      // must recvSelf

    case INTEGRATOR_TAGS_CentralDifferenceNoDamping:  
	     return new CentralDifferenceNoDamping();      // must recvSelf

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

        case INTEGRATOR_TAGS_StagedNewmark:
        return new StagedNewmark();

    case INTEGRATOR_TAGS_NewmarkExplicit:  
	     return new NewmarkExplicit();

    case INTEGRATOR_TAGS_NewmarkHSFixedNumIter:  
	     return new NewmarkHSFixedNumIter();

    case INTEGRATOR_TAGS_NewmarkHSIncrLimit:  
	     return new NewmarkHSIncrLimit();

    case INTEGRATOR_TAGS_NewmarkHSIncrReduct:  
	     return new NewmarkHSIncrReduct();

    case INTEGRATOR_TAGS_PFEMIntegrator:
        return new PFEMIntegrator();

    case INTEGRATOR_TAGS_TRBDF2:  
	     return new TRBDF2();
            
    case INTEGRATOR_TAGS_TRBDF3:  
        return new TRBDF3();

    case INTEGRATOR_TAGS_WilsonTheta:  
        return new WilsonTheta();

	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewTransientIntegrator - ";
	     opserr << " - no TransientIntegrator type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }
}


IncrementalIntegrator *
FEM_ObjectBrokerAllClasses::getNewIncrementalIntegrator(int classTag)
{
    switch(classTag) {
	case INTEGRATOR_TAGS_LoadControl:  
	     return new LoadControl(1.0,1,1.0,1.0); // must recvSelf
	    
	     
	case INTEGRATOR_TAGS_ArcLength:  
	     return new ArcLength(1.0);      // must recvSelf
	     	     
	     
	case INTEGRATOR_TAGS_Newmark:  
	     return new Newmark();

#ifdef _PARALLEL_PROCESSING	     
	case INTEGRATOR_TAGS_DistributedDisplacementControl:  
	     return new DistributedDisplacementControl(); // must recvSelf
#endif
	     
	default:
	     opserr << "FEM_ObjectBrokerAllClasses::getNewIncrementalIntegrator - ";
	     opserr << " - no IncrementalIntegrator type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }
}

LinearSOE *
FEM_ObjectBrokerAllClasses::getNewLinearSOE(int classTagSOE)
{
    LinearSOE *theSOE =0;

    /*
      case LinSOE_TAGS_SlowLinearSOE:  
	if (classTagSolver == SOLVER_TAGS_SlowLinearSOESolver) {
	    theSlowSolver = new SlowLinearSOESolver();
	    theSOE = new SlowLinearSOE(*theSlowSolver);
	    lastLinearSolver = theSlowSolver;
	    return theSOE;
	} else {
	    opserr << "FEM_ObjectBrokerAllClasses::getNewLinearSOE - ";
	    opserr << " - no SlowLinearSOESolver type exists for class tag ";
	    opserr << classTagSolver << endln;
	    return 0;		 
	}
	
	*/

    
    switch(classTagSOE) {

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
        theSOE = new PetscSOE(*( new PetscSolver()));
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
	  opserr << "FEM_ObjectBrokerAllClasses::getNewLinearSOE - ";
	  opserr << " - no LinearSOE type exists for class tag ";
	  opserr << classTagSOE << endln;
	  return 0;
	  
      
    }
}


EigenSOE *
FEM_ObjectBrokerAllClasses::getNewEigenSOE(int classTagSOE)
{
    EigenSOE *theSOE =0;

    switch(classTagSOE) {

	case EigenSOE_TAGS_ArpackSOE:  
	  theSOE = new ArpackSOE();
	  return theSOE;

	default:
	  opserr << "FEM_ObjectBrokerAllClasses::getNewEigenSOE - ";
	  opserr << " - no EigenSOE type exists for class tag ";
	  opserr << classTagSOE << endln;
	  return 0;
	  
      
    }
}




DomainSolver *
FEM_ObjectBrokerAllClasses::getNewDomainSolver(void)
{
    return lastDomainSolver;
}
    
LinearSOE *
FEM_ObjectBrokerAllClasses::getPtrNewDDLinearSOE(int classTagSOE, 
				       int classTagDDSolver)
{
    ProfileSPDLinSubstrSolver *theProfileSPDSolver =0;    

    switch(classTagSOE) {
      case LinSOE_TAGS_ProfileSPDLinSOE:  

	if (classTagDDSolver == SOLVER_TAGS_ProfileSPDLinSubstrSolver) {
	    theProfileSPDSolver = new ProfileSPDLinSubstrSolver();
	    LinearSOE *theSOE = new ProfileSPDLinSOE(*theProfileSPDSolver);
	    lastDomainSolver = theProfileSPDSolver;
	    return theSOE;		 
	}
	else {
	    opserr << "FEM_ObjectBrokerAllClasses::getNewLinearSOE - ";
	    opserr << " - no ProfileSPD Domain Solver type exists for class tag ";
	    opserr << classTagDDSolver << endln;
	    return 0;		 
	}	     
	
					    
      default:
	opserr << "FEM_ObjectBrokerAllClasses::getNewLinearSOE - ";
	opserr << " - no LinearSOE type exists for class tag ";
	opserr << classTagSOE << endln;
	return 0;
	
    }
}


DomainDecompositionAnalysis *
FEM_ObjectBrokerAllClasses::getNewDomainDecompAnalysis(int classTag, 
						Subdomain &theSubdomain)
{
    switch(classTag) {
      case DomDecompANALYSIS_TAGS_DomainDecompositionAnalysis:  
	return new DomainDecompositionAnalysis(theSubdomain);

#ifdef _PARALLEL_PROCESSING
      case ANALYSIS_TAGS_StaticDomainDecompositionAnalysis:
	return new StaticDomainDecompositionAnalysis(theSubdomain);      

      case ANALYSIS_TAGS_TransientDomainDecompositionAnalysis:
	return new TransientDomainDecompositionAnalysis(theSubdomain);      
#endif
	
      default:
	opserr << "ObjectBrokerAllClasses::getNewDomainDecompAnalysis ";
	opserr << " - no DomainDecompAnalysis type exists for class tag " ;
	opserr << classTag << endln;
	return 0;
	
    }
}


Subdomain 	  *
FEM_ObjectBrokerAllClasses::getSubdomainPtr(int classTag)
{
    opserr << "FEM_ObjectBrokerAllClasses: NOT IMPLEMENTED YET";
    return 0;
}

/*
int 
FEM_ObjectBrokerAllClasses::addUniaxialMaterial(int classTag, 
				      const char *lib, 
				      const char *funcName, 
				      UniaxialMaterial *(*funcPtr)(void))
{
  // check to see if it's already added

  UniaxialPackage *matCommands = theUniaxialPackage;
  bool found = false;
  while (matCommands != NULL && found == false) {
    if ((strcmp(lib, matCommands->libName) == 0) && (strcmp(funcName, matCommands->funcName) == 0)) {
      return 0;
    }
  }

  //
  // if funPtr == 0; go get the handle
  //

  void *libHandle;
  if (funcPtr == 0) {
    if (getLibraryFunction(lib, funcName, &libHandle, (void **)&funcPtr) != 0) {
      opserr << "FEM_ObjectBrokerAllClasses::addUniaxialMaterial - could not find function\n";
      return -1;
    }
  } 
  
  //
  // add the new funcPtr
  //
  
  char *libNameCopy = new char[strlen(lib)+1];
  char *funcNameCopy = new char[strlen(funcName)+1];
  UniaxialPackage *theMat = new UniaxialPackage;
  if (libNameCopy == 0 || funcNameCopy == 0 || theMat == 0) {
      opserr << "FEM_ObjectBrokerAllClasses::addUniaxialMaterial - could not add lib, out of memory\n";
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
*/

Parameter *
FEM_ObjectBrokerAllClasses::getParameter(int classTag)
{
  Parameter *theRes = 0;

  switch(classTag) {
  case  PARAMETER_TAG_Parameter:
    theRes = new Parameter;
    break;

  case  PARAMETER_TAG_ElementParameter:
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

  default:
    ;
  }

  return theRes;
}

