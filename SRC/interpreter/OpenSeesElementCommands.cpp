/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */


// Written: Minjie

// Description: command to create element

#include <Element.h>
#include <elementAPI.h>
#include <map>
#include <Domain.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Block2D.h>
#include <Block3D.h>
#include <FourNodeQuad.h>
#include <SectionForceDeformation.h>
#include <ConstantPressureVolumeQuad.h>
#include <EnhancedQuad.h>
#include <SSPquad.h>
#include <SSPbrick.h>
#include <Brick.h>
#include <BbarBrick.h>
#include <ShellMITC4.h>
#include <FourNodeTetrahedron.h>

// no 'beamWithHinges', 'GenericClient', 'GenericCopy', 'flBrick', 'Adapter'

void* OPS_ZeroLengthND();
void* OPS_ZeroLengthSection();
void* OPS_ZeroLength();
void* OPS_TrussElement();
void* OPS_TrussSectionElement();
void* OPS_CorotTrussElement();
void* OPS_CorotTrussSectionElement();
void* OPS_ZeroLengthContactNTS2D();
void* OPS_ZeroLengthInterface2D();
void* OPS_ComponentElement2d();
void* OPS_ZeroLengthImpact3D();
void* OPS_ModElasticBeam2d();
void* OPS_ElasticTimoshenkoBeam2d();
void* OPS_ElasticTimoshenkoBeam3d();
extern "C" void* OPS_PY_Macro2D();
void* OPS_SimpleContact2D();
void* OPS_N4BiaxialTruss();
void* OPS_SimpleContact3D();
void* OPS_BeamContact3D();
void* OPS_BeamContact3Dp();
void* OPS_PileToe3D();
void* OPS_TFP_Bearing();
void* OPS_FPBearingPTV();
void* OPS_TripleFrictionPendulum();
void* OPS_HDR();
void* OPS_LeadRubberX();
void* OPS_ElastomericX();
void* OPS_MVLEM();
void* OPS_SFI_MVLEM();
void* OPS_MultiFP2d();
void* OPS_ShellMITC4();
void* OPS_ShellMITC9();
void* OPS_ShellDKGQ();
void* OPS_ShellNLDKGQ();
void* OPS_CoupledZeroLength();
void* OPS_BeamContact2D();
void* OPS_BeamContact2Dp();
void* OPS_BeamEndContact3D();
void* OPS_BeamEndContact3Dp();
void* OPS_Tri31(const ID& info);
void* OPS_SSPquad();
void* OPS_SSPquadUP();
void* OPS_SSPbrick();
void* OPS_SSPbrickUP();
void* OPS_SurfaceLoad();
void* OPS_TPB1D();
void* OPS_ElasticTubularJoint();
void* OPS_FourNodeQuad3d();
void* OPS_Quad4FiberOverlay();
void* OPS_Brick8FiberOverlay();
void* OPS_QuadBeamEmbedContact();
void* OPS_Truss2();
void* OPS_CorotTruss2();
void* OPS_AC3D8HexWithSensitivity();
void* OPS_AV3D4QuadWithSensitivity();
void* OPS_ElastomericBearingBoucWenMod3d();
void* OPS_VS3D4WuadWithSensitivity();
void* OPS_PFEMElement2DBubble(const ID& info);
void* OPS_PFEMElement3DBubble(const ID& info);
//void* OPS_TaylorHood2D();
void* OPS_PFEMElement2DCompressible(const ID& info);
void* OPS_PFEMElement2Dmini(const ID& info);
void* OPS_fElmt02();
void* OPS_ElasticBeam2d(const ID& info);
void* OPS_ElasticBeam3d();
void* OPS_DispBeamColumn2dInt();
void* OPS_ForceBeamColumn2d(const ID& info);
void* OPS_ForceBeamColumn3d();
void* OPS_DispBeamColumn2d(const ID& info);
void* OPS_DispBeamColumn3d();
void* OPS_ForceBeamColumnCBDI2d();
void* OPS_ForceBeamColumnCSBDI2d();
void* OPS_ForceBeamColumnWarping2d();
void* OPS_ElasticForceBeamColumnWarping2d();
void* OPS_DispBeamColumn2dThermal();
void* OPS_DispBeamColumn3dThermal();
void* OPS_ElasticForceBeamColumn2d();
void* OPS_ElasticForceBeamColumn3d();
void* OPS_DispBeamColumn3dWithSensitivity();
void* OPS_DispBeamColumn2dWithSensitivity();
void* OPS_FourNodeQuad();
void* OPS_FourNodeQuadWithSensitivity();
void* OPS_EnhancedQuad();
void* OPS_ConstantPressureVolumeQuad();
void* OPS_NineNodeMixedQuad();
void* OPS_FourNodeQuadUP();
void* OPS_BrickUP();
void* OPS_NineFourNodeQuadUP();
void* OPS_TwentyEightNodeBrickUP();
void* OPS_Twenty_Node_Brick();
void* OPS_BBarFourNodeQuadUP();
void* OPS_BBarBrickUP();
void* OPS_Brick();
void* OPS_BbarBrick();
void* OPS_BbarBrickWithSensitivity();
void* OPS_ZeroLengthRocking();
void* OPS_ZeroLengthContact2D();
void* OPS_ZeroLengthContact3D();
void* OPS_Joint2D();
void* OPS_Joint3D();
void* OPS_Inelastic2DYS01();
void* OPS_Inelastic2DYS02();
void* OPS_Inelastic2DYS03();
void* OPS_Elastic2DGNL();
void* OPS_BeamColumnJoint2d();
void* OPS_BeamColumnJoint3d();
void* OPS_Actuator();
void* OPS_ActuatorCorot();
void* OPS_FlatSliderSimple2d();
void* OPS_FlatSliderSimple3d();
void* OPS_SingleFPSimple2d();
void* OPS_SingleFPSimple3d();
void* OPS_RJWatsonEQS2d();
void* OPS_RJWatsonEQS3d();
void* OPS_ElastomericBearingPlasticity2d();
void* OPS_ElastomericBearingPlasticity3d();
void* OPS_ElastomericBearingBoucWen2d();
void* OPS_ElastomericBearingBoucWen3d();
void* OPS_ElastomericBearingUFRP2d();
void* OPS_TwoNodeLink();
void* OPS_MultipleShearSpring();
void* OPS_MultipleNormalSpring();
void* OPS_KikuchiBearing();
void* OPS_YamamotoBiaxialHDR();
void* OPS_FourNodeTetrahedron();

namespace {

    struct char_cmp {
	bool operator () (const char *a,const char *b) const
	    {
		return strcmp(a,b)<0;
	    }
    };

    typedef std::map<const char *, void *(*)(void), char_cmp> OPS_ParsingFunctionMap;


    static OPS_ParsingFunctionMap functionMap;

    static void* OPS_ForceBeamColumn()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    ID info;
	    return OPS_ForceBeamColumn2d(info);
	} else {
	    return OPS_ForceBeamColumn3d();
	}
    }

    static void* OPS_ElasticForceBeamColumn()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    return OPS_ElasticForceBeamColumn2d();
	} else {
	    return OPS_ElasticForceBeamColumn3d();
	}
    }

    static void* OPS_ElasticBeam()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    ID info;
	    return OPS_ElasticBeam2d(info);
	} else {
	    return OPS_ElasticBeam3d();
	}
    }

    static void* OPS_DispBeamColumn()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    ID info;
	    return OPS_DispBeamColumn2d(info);
	} else {
	    return OPS_DispBeamColumn3d();
	}
    }

    static void* OPS_DispBeamColumnThermal()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    return OPS_DispBeamColumn2dThermal();
	} else {
	    return OPS_DispBeamColumn3dThermal();
	}
    }

    static void* OPS_DispBeamColumnWithSensitivity()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    return OPS_DispBeamColumn2dWithSensitivity();
	} else {
	    return OPS_DispBeamColumn3dWithSensitivity();
	}
    }

    static void* OPS_ElasticTimoshenkoBeam()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    return OPS_ElasticTimoshenkoBeam2d();
	} else {
	    return OPS_ElasticTimoshenkoBeam3d();
	}
    }

    static void* OPS_BeamColumnJoint()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    return OPS_BeamColumnJoint2d();
	} else {
	    return OPS_BeamColumnJoint3d();
	}
    }

    static void* OPS_FlatSliderBearing()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    return OPS_FlatSliderSimple2d();
	} else {
	    return OPS_FlatSliderSimple3d();
	}
    }

    static void* OPS_SingleFPBearing()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    return OPS_SingleFPSimple2d();
	} else {
	    return OPS_SingleFPSimple3d();
	}
    }

    static void* OPS_RJWatsonEqsBearing()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    return OPS_RJWatsonEQS2d();
	} else {
	    return OPS_RJWatsonEQS3d();
	}
    }

    static void* OPS_ElastomericBearingPlasticity()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    return OPS_ElastomericBearingPlasticity2d();
	} else {
	    return OPS_ElastomericBearingPlasticity3d();
	}
    }

    static void* OPS_ElastomericBearingBoucWen()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    return OPS_ElastomericBearingBoucWen2d();
	} else {
	    return OPS_ElastomericBearingBoucWen3d();
	}
    }

    static void* OPS_ElastomericBearingUFRP()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    return OPS_ElastomericBearingUFRP2d();
	} else {
	    return 0;
	}
    }

    static void* OPS_PFEMElementBubble()
    {
	int ndm = OPS_GetNDM();
	ID info;
	if(ndm == 2) {
	    return OPS_PFEMElement2DBubble(info);
	} else {
	    return OPS_PFEMElement3DBubble(info);;
	}
    }

    static void* OPS_PFEMElementmini()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    ID info;
	    return OPS_PFEMElement2Dmini(info);
	} else {
	    return 0;
	}
    }

    static void* OPS_PFEMElementCompressible()
    {
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    ID info;
	    return OPS_PFEMElement2DCompressible(info);
	} else {
	    return 0;
	}
    }

    static void* OPS_Tri31NoInfo()
    {
	ID info;
	return OPS_Tri31(info);
    }

    static int setUpFunctions(void)
    {
	functionMap.insert(std::make_pair("KikuchiBearing", &OPS_KikuchiBearing));
	functionMap.insert(std::make_pair("YamamotoBiaxialHDR", &OPS_YamamotoBiaxialHDR));
	functionMap.insert(std::make_pair("MNS", &OPS_MultipleNormalSpring));
	functionMap.insert(std::make_pair("multipleNormalSpring", &OPS_MultipleNormalSpring));
	functionMap.insert(std::make_pair("MSS", &OPS_MultipleShearSpring));
	functionMap.insert(std::make_pair("multipleShearSpring", &OPS_MultipleShearSpring));
	functionMap.insert(std::make_pair("twoNodeLink", &OPS_TwoNodeLink));
	functionMap.insert(std::make_pair("elastomericBearingUFRP", &OPS_ElastomericBearingUFRP));
	functionMap.insert(std::make_pair("elastomericBearingPlasticity", &OPS_ElastomericBearingPlasticity));
	functionMap.insert(std::make_pair("ElastomericBearingBoucWen", &OPS_ElastomericBearingBoucWen));
	functionMap.insert(std::make_pair("elastomericBearing", &OPS_ElastomericBearingPlasticity));
	functionMap.insert(std::make_pair("RJWatsonEqsBearing", &OPS_RJWatsonEqsBearing));
	functionMap.insert(std::make_pair("singleFPBearing", &OPS_SingleFPBearing));
	functionMap.insert(std::make_pair("flatSliderBearing", &OPS_FlatSliderBearing));
	functionMap.insert(std::make_pair("corotActuator", &OPS_ActuatorCorot));
	functionMap.insert(std::make_pair("actuator", &OPS_Actuator));
	functionMap.insert(std::make_pair("beamColumnJoint", &OPS_BeamColumnJoint));
	functionMap.insert(std::make_pair("elastic2dGNL", &OPS_Elastic2DGNL));
	functionMap.insert(std::make_pair("element2dGNL", &OPS_Elastic2DGNL));
	functionMap.insert(std::make_pair("inelastic2dYS03", &OPS_Inelastic2DYS03));
	functionMap.insert(std::make_pair("inelastic2dYS02", &OPS_Inelastic2DYS02));
	functionMap.insert(std::make_pair("inelastic2dYS01", &OPS_Inelastic2DYS01));
	functionMap.insert(std::make_pair("Joint3D", &OPS_Joint3D));
	functionMap.insert(std::make_pair("Joint3d", &OPS_Joint3D));
	functionMap.insert(std::make_pair("Joint2D", &OPS_Joint2D));
	functionMap.insert(std::make_pair("Joint2d", &OPS_Joint2D));
	functionMap.insert(std::make_pair("zeroLengthContact2D", &OPS_ZeroLengthContact2D));
	functionMap.insert(std::make_pair("zeroLengthContact3D", &OPS_ZeroLengthContact3D));
	functionMap.insert(std::make_pair("zeroLengthRocking", &OPS_ZeroLengthRocking));
	functionMap.insert(std::make_pair("bbarBrickWithSensitivity", &OPS_BbarBrickWithSensitivity));
	functionMap.insert(std::make_pair("bbarBrick", &OPS_BbarBrick));
	functionMap.insert(std::make_pair("stdBrick", &OPS_Brick));
	functionMap.insert(std::make_pair("bbarBrickUP", &OPS_BBarBrickUP));
	functionMap.insert(std::make_pair("bbarQuadUP", &OPS_BBarFourNodeQuadUP));
	functionMap.insert(std::make_pair("20NodeBrick", &OPS_Twenty_Node_Brick));
	functionMap.insert(std::make_pair("20_8_BrickUP", &OPS_TwentyEightNodeBrickUP));
	functionMap.insert(std::make_pair("9_4_QuadUP", &OPS_NineFourNodeQuadUP));
	functionMap.insert(std::make_pair("brickUP", &OPS_BrickUP));
	functionMap.insert(std::make_pair("quadUP", &OPS_FourNodeQuadUP));
	functionMap.insert(std::make_pair("nineNodeMixedQuad", &OPS_NineNodeMixedQuad));
	functionMap.insert(std::make_pair("nineNodeQuad", &OPS_NineNodeMixedQuad));
	functionMap.insert(std::make_pair("bbarQuad", &OPS_ConstantPressureVolumeQuad));
	functionMap.insert(std::make_pair("mixedQuad", &OPS_ConstantPressureVolumeQuad));
	functionMap.insert(std::make_pair("enhancedQuad", &OPS_EnhancedQuad));
	functionMap.insert(std::make_pair("quadWithSensitivity", &OPS_FourNodeQuadWithSensitivity));
	functionMap.insert(std::make_pair("quad", &OPS_FourNodeQuad));
	functionMap.insert(std::make_pair("stdQuad", &OPS_FourNodeQuad));
	functionMap.insert(std::make_pair("dispBeamColumnWithSensitivity", &OPS_DispBeamColumnWithSensitivity));
	functionMap.insert(std::make_pair("elasticForceBeamColumn", &OPS_ElasticForceBeamColumn));
	functionMap.insert(std::make_pair("dispBeamColumnThermal", &OPS_DispBeamColumnThermal));
	functionMap.insert(std::make_pair("forceBeamColumnWarping", &OPS_ForceBeamColumnWarping2d));
	functionMap.insert(std::make_pair("elasticForceBeamColumnWarping", &OPS_ElasticForceBeamColumnWarping2d));
	functionMap.insert(std::make_pair("dispBeamColumnInt", &OPS_DispBeamColumn2dInt));
	functionMap.insert(std::make_pair("fTruss", &OPS_fElmt02));
	functionMap.insert(std::make_pair("PFEMElementCompressible", &OPS_PFEMElementCompressible));
	functionMap.insert(std::make_pair("PFEMElementBubble", &OPS_PFEMElementBubble));
	functionMap.insert(std::make_pair("MINI", &OPS_PFEMElementmini));
	//functionMap.insert(std::make_pair("TaylorHood2D", &OPS_TaylorHood2D));
	functionMap.insert(std::make_pair("VS3D4", &OPS_VS3D4WuadWithSensitivity));
	functionMap.insert(std::make_pair("elastomericBearingBoucWenMod", &OPS_ElastomericBearingBoucWenMod3d));
	functionMap.insert(std::make_pair("AV3D4", &OPS_AV3D4QuadWithSensitivity));
	functionMap.insert(std::make_pair("AC3D8", &OPS_AC3D8HexWithSensitivity));
	functionMap.insert(std::make_pair("CorotTruss2", &OPS_CorotTruss2));
	functionMap.insert(std::make_pair("Truss2", &OPS_Truss2));
	functionMap.insert(std::make_pair("QuadBeamEmbedContact", &OPS_QuadBeamEmbedContact));
	functionMap.insert(std::make_pair("Brick8FiberOverlay", &OPS_Brick8FiberOverlay));
	functionMap.insert(std::make_pair("Quad4FiberOverlay", &OPS_Quad4FiberOverlay));
	functionMap.insert(std::make_pair("quad3d", &OPS_FourNodeQuad3d));
	functionMap.insert(std::make_pair("Quad3d", &OPS_FourNodeQuad3d));
	functionMap.insert(std::make_pair("elasticTubularJoint", &OPS_ElasticTubularJoint));
	functionMap.insert(std::make_pair("ElasticTubularJoint", &OPS_ElasticTubularJoint));
	functionMap.insert(std::make_pair("TPB1D", &OPS_TPB1D));
	functionMap.insert(std::make_pair("Truss", &OPS_TrussElement));
	functionMap.insert(std::make_pair("truss", &OPS_TrussElement));
	functionMap.insert(std::make_pair("trussSection", &OPS_TrussSectionElement));
	functionMap.insert(std::make_pair("TrussSection", &OPS_TrussSectionElement));
	functionMap.insert(std::make_pair("corotTruss", &OPS_CorotTrussElement));
	functionMap.insert(std::make_pair("CorotTruss", &OPS_CorotTrussElement));
	functionMap.insert(std::make_pair("corotTrussSection", &OPS_CorotTrussSectionElement));
	functionMap.insert(std::make_pair("CorotTrussSection", &OPS_CorotTrussSectionElement));
	functionMap.insert(std::make_pair("zeroLengthContactNTS2D", &OPS_ZeroLengthContactNTS2D));
	functionMap.insert(std::make_pair("zeroLengthInterface2D", &OPS_ZeroLengthInterface2D));
	functionMap.insert(std::make_pair("componentElement2d", &OPS_ComponentElement2d));
	functionMap.insert(std::make_pair("zeroLengthImpact3D", &OPS_ZeroLengthImpact3D));
	functionMap.insert(std::make_pair("ModElasticBeam2d", &OPS_ModElasticBeam2d));
	functionMap.insert(std::make_pair("modElasticBeam2d", &OPS_ModElasticBeam2d));
	functionMap.insert(std::make_pair("ElasticTimoshenkoBeam", &OPS_ElasticTimoshenkoBeam));
	functionMap.insert(std::make_pair("elasticTimoshenkoBeam", &OPS_ElasticTimoshenkoBeam));
	functionMap.insert(std::make_pair("pyMacro2D", &OPS_PY_Macro2D));
	functionMap.insert(std::make_pair("PY_Macro2D", &OPS_PY_Macro2D));
	functionMap.insert(std::make_pair("SimpleContact2d", &OPS_SimpleContact2D));
	functionMap.insert(std::make_pair("SimpleContact2D", &OPS_SimpleContact2D));
	functionMap.insert(std::make_pair("N4BiaxialTruss", &OPS_N4BiaxialTruss));
	functionMap.insert(std::make_pair("SimpleContact3d", &OPS_SimpleContact3D));
	functionMap.insert(std::make_pair("SimpleContact3D", &OPS_SimpleContact3D));
	functionMap.insert(std::make_pair("BeamContact3d", &OPS_BeamContact3D));
	functionMap.insert(std::make_pair("BeamContact3D", &OPS_BeamContact3D));
	functionMap.insert(std::make_pair("BeamContact3dp", &OPS_BeamContact3Dp));
	functionMap.insert(std::make_pair("BeamContact3Dp", &OPS_BeamContact3Dp));
	functionMap.insert(std::make_pair("PileToe3d", &OPS_PileToe3D));
	functionMap.insert(std::make_pair("PileToe3D", &OPS_PileToe3D));
	functionMap.insert(std::make_pair("TFPbearing", &OPS_TFP_Bearing));
	functionMap.insert(std::make_pair("TFP", &OPS_TFP_Bearing));
	functionMap.insert(std::make_pair("FPBearingPTV", &OPS_FPBearingPTV));
	functionMap.insert(std::make_pair("TripleFrictionPendulum", &OPS_TripleFrictionPendulum));
	functionMap.insert(std::make_pair("HDR", &OPS_HDR));
	functionMap.insert(std::make_pair("LeadRubberX", &OPS_LeadRubberX));
	functionMap.insert(std::make_pair("ElastomericX", &OPS_ElastomericX));
	functionMap.insert(std::make_pair("MVLEM", &OPS_MVLEM));
	functionMap.insert(std::make_pair("SFI_MVLEM", &OPS_SFI_MVLEM));
	functionMap.insert(std::make_pair("MultiFP2d", &OPS_MultiFP2d));
	functionMap.insert(std::make_pair("shell", &OPS_ShellMITC4));
	functionMap.insert(std::make_pair("Shell", &OPS_ShellMITC4));
	functionMap.insert(std::make_pair("shellMITC4", &OPS_ShellMITC4));
	functionMap.insert(std::make_pair("ShellMITC4", &OPS_ShellMITC4));
	functionMap.insert(std::make_pair("shellNL", &OPS_ShellMITC9));
	functionMap.insert(std::make_pair("ShellNL", &OPS_ShellMITC9));
	functionMap.insert(std::make_pair("shellMITC9", &OPS_ShellMITC9));
	functionMap.insert(std::make_pair("ShellMITC9", &OPS_ShellMITC9));
	functionMap.insert(std::make_pair("ShellDKGQ", &OPS_ShellDKGQ));
	functionMap.insert(std::make_pair("shellDKGQ", &OPS_ShellDKGQ));
	functionMap.insert(std::make_pair("ShellNLDKGQ", &OPS_ShellNLDKGQ));
	functionMap.insert(std::make_pair("shellNLDKGQ", &OPS_ShellNLDKGQ));
	functionMap.insert(std::make_pair("CoupledZeroLength", &OPS_CoupledZeroLength));
	functionMap.insert(std::make_pair("ZeroLengthCoupled", &OPS_CoupledZeroLength));
	functionMap.insert(std::make_pair("BeamContact2d", &OPS_BeamContact2D));
	functionMap.insert(std::make_pair("BeamContact2D", &OPS_BeamContact2D));
	functionMap.insert(std::make_pair("BeamContact2dp", &OPS_BeamContact2Dp));
	functionMap.insert(std::make_pair("BeamContact2Dp", &OPS_BeamContact2Dp));
	functionMap.insert(std::make_pair("BeamEndContact3d", &OPS_BeamEndContact3D));
	functionMap.insert(std::make_pair("BeamEndContact3D", &OPS_BeamEndContact3D));
	functionMap.insert(std::make_pair("BeamEndContact3dp", &OPS_BeamEndContact3Dp));
	functionMap.insert(std::make_pair("BeamEndContact3Dp", &OPS_BeamEndContact3Dp));
	functionMap.insert(std::make_pair("Tri31", &OPS_Tri31NoInfo));
	functionMap.insert(std::make_pair("tri31", &OPS_Tri31NoInfo));
	functionMap.insert(std::make_pair("SSPquad", &OPS_SSPquad));
	functionMap.insert(std::make_pair("SSPQuad", &OPS_SSPquad));
	functionMap.insert(std::make_pair("SSPquadUP", &OPS_SSPquadUP));
	functionMap.insert(std::make_pair("SSPQuadUP", &OPS_SSPquadUP));
	functionMap.insert(std::make_pair("SSPbrick", &OPS_SSPbrick));
	functionMap.insert(std::make_pair("SSPBrick", &OPS_SSPbrick));
	functionMap.insert(std::make_pair("SSPbrickUP", &OPS_SSPbrickUP));
	functionMap.insert(std::make_pair("SSPBrickUP", &OPS_SSPbrickUP));
	functionMap.insert(std::make_pair("SurfaceLoad", &OPS_SurfaceLoad));
	functionMap.insert(std::make_pair("elasticBeamColumn", &OPS_ElasticBeam));
	functionMap.insert(std::make_pair("forceBeamColumn", &OPS_ForceBeamColumn));
	functionMap.insert(std::make_pair("nonlinearBeamColumn", &OPS_ForceBeamColumn));
	functionMap.insert(std::make_pair("dispBeamColumn", &OPS_DispBeamColumn));
	functionMap.insert(std::make_pair("forceBeamColumnCBDI", &OPS_ForceBeamColumnCBDI2d));
	functionMap.insert(std::make_pair("forceBeamColumnCSBDI", &OPS_ForceBeamColumnCSBDI2d));
	functionMap.insert(std::make_pair("zeroLength", &OPS_ZeroLength));
	functionMap.insert(std::make_pair("zeroLengthSection", &OPS_ZeroLengthSection));
	functionMap.insert(std::make_pair("zeroLengthND", &OPS_ZeroLengthND));
	functionMap.insert(std::make_pair("FourNodeTetrahedron", &OPS_FourNodeTetrahedron));

	return 0;
    }
}

int
OPS_Element()
{
    static bool initDone = false;
    if (initDone == false) {
	setUpFunctions();
	initDone = true;
    }

    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr<<"WARNING too few arguments: element type? tag? ...\n";
	return -1;
    }

    const char* type = OPS_GetString();

    OPS_ParsingFunctionMap::const_iterator iter = functionMap.find(type);
    if (iter == functionMap.end()) {
	opserr<<"WARNING element type " << type << " is unknown\n";
	return -1;
    }

    Element* theEle = (Element*) (*iter->second)();
    if (theEle == 0) {
	// for backward comatability
	if (strcmp(type, "truss")==0 || strcmp(type, "Truss")==0) {
	    theEle = (Element*) OPS_TrussSectionElement();
	    if (theEle == 0) return -1;
	} else {
	    return -1;
	}
    }

    // Now add the element to the domain
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    if (theDomain->addElement(theEle) == false) {
	opserr<<"ERROR could not add element to domain.\n";
	delete theEle;
	return -1;
    }

    return 0;

}

int OPS_doBlock2D()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    if (ndm < 2) {
	opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs? coords?";
	opserr << " : model dimension (ndm) must be at leat 2 \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() < 7) {
	opserr << "WARNING incorrect numer of args :block2D numX? numY? startNode? startEle? eleType? eleArgs? coords?";
	return -1;
    }

    // numX, numY, startNodeNum, startEleNum
    int idata[4];
    int numdata = 4;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid int inputs\n";
	return -1;
    }

    // element type
    const char* type = OPS_GetString();

    // get args
    const char* subtype = "";
    double thick = 1.0;
    int matTag=-1, secTag=-1;
    int cArg = 6;
    if (strcmp(type, "quad") == 0  || (strcmp(type,"stdQuad") == 0)) {
	if (OPS_GetNumRemainingInputArgs() < 3) {
	    opserr<<"WARNING: want - thick, type, matTag\n";
	    return -1;
	}
	int numdata = 1;
	if (OPS_GetDoubleInput(&numdata, &thick) < 0) {
	    opserr << "WARNING invalid thick\n";
	    return -1;
	}
	subtype = OPS_GetString();
	if (OPS_GetIntInput(&numdata, &matTag) < 0) {
	    opserr << "WARNING invalid matTag\n";
	    return -1;
	}
	cArg = 9;

    } else if (strcmp(type, "ShellMITC4") == 0 || strcmp(type, "shellMITC4") == 0 ||
	       strcmp(type, "shell") == 0 || strcmp(type, "Shell") == 0) {
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr<<"WARNING: want - secTag\n";
	    return -1;
	}
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &secTag) < 0) {
	    opserr << "WARNING invalid secTag\n";
	    return -1;
	}
	cArg = 7;

    } else if (strcmp(type, "bbarQuad") == 0 || strcmp(type,"mixedQuad") == 0) {
	if (OPS_GetNumRemainingInputArgs() < 2) {
	    opserr<<"WARNING: want - thick, matTag\n";
	    return -1;
	}
    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &thick) < 0) {
        opserr << "WARNING invalid thick\n";
        return -1;
    }
	if (OPS_GetIntInput(&numdata, &matTag) < 0) {
	    opserr << "WARNING invalid matTag\n";
	    return -1;
	}
	cArg = 8;

    } else if (strcmp(type, "enhancedQuad") == 0) {
	if (OPS_GetNumRemainingInputArgs() < 3) {
	    opserr<<"WARNING: want - thick, type, matTag\n";
	    return -1;
	}
    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &thick) < 0) {
        opserr << "WARNING invalid thick\n";
        return -1;
    }
	subtype = OPS_GetString();
	if (OPS_GetIntInput(&numdata, &matTag) < 0) {
	    opserr << "WARNING invalid matTag\n";
	    return -1;
	}
	cArg = 9;

    } else if (strcmp(type, "SSPquad") == 0 || strcmp(type, "SSPQuad") == 0) {
	if (OPS_GetNumRemainingInputArgs() < 3) {
	    opserr<<"WARNING: want - matTag, type, thick\n";
	    return -1;
	}
	int numdata = 1;
    if (OPS_GetIntInput(&numdata, &matTag) < 0) {
        opserr << "WARNING invalid matTag\n";
        return -1;
    }
	subtype = OPS_GetString();
    if (OPS_GetDoubleInput(&numdata, &thick) < 0) {
        opserr << "WARNING invalid thick\n";
        return -1;
    }
    cArg = 9;


    } else {
	opserr << "WARNING element type "<<type<<" is currently unknown by this command.\n";
	return -1;
    }

    // get num nodes
    int numEleNodes = 4;
    if (OPS_GetNumRemainingInputArgs() > 1) {
	const char* opt = OPS_GetString();
	if (strcmp(opt, "-numEleNodes") == 0) {
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &numEleNodes) < 0) {
		opserr<<"WARNING invalid numEleNodes\n";
		return -1;
	    }
	    if (numEleNodes != 4 && numEleNodes != 9) {
		opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs? ";
		opserr << "-numEleNodes numNodes?: invalid numNodes: " << numEleNodes << " 4 or 9 only\n";
		return -1;
	    }
	    if (numEleNodes == 9) {
		if (((idata[0] % 2) != 0) || ((idata[1] % 2) != 0)) {
		    opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs? ";
		    opserr << "-numEleNodes 9: numX and numY MUST BOTH BE EVEN\n";
		    return -1;
		}
	    }
	} else {
	    OPS_ResetCurrentInputArg(cArg);
	}
    }

    // get coords
    Matrix Coordinates(9,3);
    ID haveNode(9);
    Coordinates.Zero();
    for (int k=0; k<9; k++) {
	haveNode(k) = -1;
    }

    int numnodes = OPS_GetNumRemainingInputArgs() / (ndm+1);
    if (numnodes < 4) {
	opserr<<"WARNING four points (1-4) are required\n";
	return -1;
    }

    if (numnodes > 9) numnodes = 9;

    for (int i=0; i<numnodes; i++) {
	int tag;
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
	    opserr<<"WARNING failed to get node tag\n";
	    return -1;
	}
	haveNode(tag-1) = tag;

	Vector crds(ndm);
	if (OPS_GetDoubleInput(&ndm, &crds(0)) < 0) {
	    opserr<<"WARNING failed to get coordinates\n";
	    return -1;
	}
	for (int j=0; j<ndm; j++) {
	    Coordinates(tag-1,j) = crds(j);
	}
    }

    // create Block2D object
    Block2D theBlock(idata[0], idata[1], haveNode, Coordinates, numEleNodes);

    // create the nodes: (numX+1)*(numY+1) nodes to be created
    int nodeID = idata[2];
    Node* theNode = 0;
    for (int j=0; j<=idata[1]; j++) {
	for (int i=0; i<=idata[0]; i++) {
	    const Vector& nodeCoords = theBlock.getNodalCoords(i,j);
	    double xLoc = nodeCoords(0);
	    double yLoc = nodeCoords(1);

	    if (ndm == 2) {
		theNode = new Node(nodeID,ndf,xLoc, yLoc);
	    } else if (ndm == 3) {
		double zLoc = nodeCoords(2);
		theNode = new Node(nodeID,ndf,xLoc, yLoc, zLoc);
	    }

	    if (theNode == 0) {
		opserr << "WARNING ran out of memory creating node\n";
		opserr << "node: " << nodeID << "\n";
		return -1;
	    }

	    if (theDomain->addNode(theNode) == false) {
		opserr << "WARNING failed to add node to the domain\n";
		opserr << "node: " << nodeID << endln;
		delete theNode;
		return -1;
	    }
	    nodeID++;
	}
    }

    // create the elements: numX*numY elements to be created if 4 node elements
    //                      numX/2 * numY /2 nodes to be created if 9 node elements
    int eleID = idata[3];
    if (numnodes == 9) {
	idata[0] /= 2;
	idata[1] /= 2;
    }

    Element* theEle = 0;
    for (int j=0; j<idata[1]; j++) {
	for (int i=0; i<idata[0]; i++) {
	    const ID& nodeTags = theBlock.getElementNodes(i,j);

	    if (strcmp(type, "quad") == 0  || (strcmp(type,"stdQuad") == 0)) {

		if (numEleNodes != 4) {
		    opserr<<"WARNING quad element only needs four nodes\n";
		    return -1;
		}
		NDMaterial* mat = OPS_getNDMaterial(matTag);
		if (mat == 0) {
		    opserr << "WARNING material not found\n";
		    opserr << "Material: " << matTag;
		    opserr << "\nFourNodeQuad \n";
		    return -1;
		}
		int nd1 = nodeTags(0) + idata[2];
		int nd2 = nodeTags(1) + idata[2];
		int nd3 = nodeTags(2) + idata[2];
		int nd4 = nodeTags(3) + idata[2];
		theEle = new FourNodeQuad(eleID,nd1,nd2,nd3,nd4,*mat,subtype,thick);


	    } else if (strcmp(type, "ShellMITC4") == 0 || strcmp(type, "shellMITC4") == 0 ||
		       strcmp(type, "shell") == 0 || strcmp(type, "Shell") == 0) {

		if (numEleNodes != 4) {
		    opserr<<"WARNING ShellMITC4 element only needs four nodes\n";
		    return -1;
		}
		SectionForceDeformation *sec = OPS_getSectionForceDeformation(secTag);

		if (sec == 0) {
		    opserr << "WARNING:  section " << secTag << " not found\n";
		    return -1;
		}

		int nd1 = nodeTags(0) + idata[2];
		int nd2 = nodeTags(1) + idata[2];
		int nd3 = nodeTags(2) + idata[2];
		int nd4 = nodeTags(3) + idata[2];
		theEle = new ShellMITC4(eleID,nd1,nd2,nd3,nd4,*sec);


	    } else if (strcmp(type, "bbarQuad") == 0 || strcmp(type,"mixedQuad") == 0) {

		if (numEleNodes != 4) {
		    opserr<<"WARNING ConstantPressureVolumeQuad element only needs four nodes\n";
		    return -1;
		}
		NDMaterial* mat = OPS_getNDMaterial(matTag);
		if (mat == 0) {
		    opserr << "WARNING material not found\n";
		    opserr << "Material: " << matTag << "\n";
		    return -1;
		}

		int nd1 = nodeTags(0) + idata[2];
		int nd2 = nodeTags(1) + idata[2];
		int nd3 = nodeTags(2) + idata[2];
		int nd4 = nodeTags(3) + idata[2];
		theEle = new ConstantPressureVolumeQuad(eleID,nd1,nd2,nd3,nd4,*mat,thick);

	    } else if (strcmp(type, "enhancedQuad") == 0) {

		if (numEleNodes != 4) {
		    opserr<<"WARNING EnhancedQuad element only needs four nodes\n";
		    return -1;
		}
		NDMaterial* mat = OPS_getNDMaterial(matTag);
		if (mat == 0) {
		    opserr << "WARNING material not found\n";
		    opserr << "Material: " << matTag << "\n";
		    return -1;
		}

		int nd1 = nodeTags(0) + idata[2];
		int nd2 = nodeTags(1) + idata[2];
		int nd3 = nodeTags(2) + idata[2];
		int nd4 = nodeTags(3) + idata[2];
		theEle = new EnhancedQuad(eleID,nd1,nd2,nd3,nd4,*mat,subtype,thick);

	    } else if (strcmp(type, "SSPquad") == 0 || strcmp(type, "SSPQuad") == 0) {

		if (numEleNodes != 4) {
		    opserr<<"WARNING SSPquad element only needs four nodes\n";
		    return -1;
		}
		NDMaterial* mat = OPS_getNDMaterial(matTag);
		if (mat == 0) {
		    opserr << "WARNING material not found\n";
		    opserr << "Material: " << matTag;
		    opserr << "\nFourNodeQuad \n";
		    return -1;
		}
		int nd1 = nodeTags(0) + idata[2];
		int nd2 = nodeTags(1) + idata[2];
		int nd3 = nodeTags(2) + idata[2];
		int nd4 = nodeTags(3) + idata[2];
		theEle = new SSPquad(eleID,nd1,nd2,nd3,nd4,*mat,subtype,thick);
	    }

	    if (theDomain->addElement(theEle) == false) {
		opserr<<"WARNING failed to add element to domain\n";
		delete theEle;
		return -1;
	    }

	    eleID++;
	}
    }


    return 0;
}

int OPS_doBlock3D()
{
    int ndm = OPS_GetNDM();
    if (ndm < 3) {
	opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
	opserr << " : model dimension (ndm) must be at leat 3 \n";
	return -1;
    }

    int ndf = OPS_GetNDF();
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    if (OPS_GetNumRemainingInputArgs() < 8) {
	opserr << "WARNING incorrect numer of args :block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs? coords?";
	return -1;
    }

    // numX, numY, numZ, startNodeNum, startEleNum
    int idata[5];
    int numdata = 5;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid int inputs\n";
	return -1;
    }

    // element type
    const char* type = OPS_GetString();

    // get mattag
    int matTag=-1;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &matTag) < 0) {
	opserr << "WARNING invalid matTag\n";
	return -1;
    }

    // get coords
    Matrix Coordinates(27,3);
    ID haveNode(27);
    Coordinates.Zero();
    for (int k=0; k<27; k++) {
	haveNode(k) = -1;
    }

    int numnodes = OPS_GetNumRemainingInputArgs() / (ndm+1);
    if (numnodes < 8) {
	opserr<<"WARNING eight points (1-8) are required\n";
	return -1;
    }

    if (numnodes > 27) numnodes = 27;

    for (int i=0; i<numnodes; i++) {
	int tag;
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
	    opserr<<"WARNING failed to get node tag\n";
	    return -1;
	}
	if (tag < 1 || tag > 27) {
	    opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
	    opserr << " : node tag out of bounds [1, 27] \n";
	    return -1;
	}
	haveNode(tag-1) = tag;

	Vector crds(ndm);
	if (OPS_GetDoubleInput(&ndm, &crds(0)) < 0) {
	    opserr<<"WARNING failed to get coordinates\n";
	    return -1;
	}
	for (int j=0; j<ndm; j++) {
	    Coordinates(tag-1,j) = crds(j);
	}
    }

    // create Block3D object
    Block3D theBlock(idata[0], idata[1], idata[2], haveNode, Coordinates);

    // create the nodes: (numX+1)*(numY+1)*(numZ+1) nodes to be created
    int nodeID = idata[3];
    Node* theNode = 0;
    for (int k=0; k<=idata[2]; k++) {
	for (int j=0; j<=idata[1]; j++) {
	    for (int i=0; i<=idata[0]; i++) {
		const Vector& nodeCoords = theBlock.getNodalCoords(i,j,k);
		double xLoc = nodeCoords(0);
		double yLoc = nodeCoords(1);
		double zLoc = nodeCoords(2);

		theNode = new Node(nodeID, ndf, xLoc, yLoc, zLoc);

		if (theNode == 0) {
		    opserr << "WARNING ran out of memory creating node\n";
		    opserr << "node: " << nodeID << "\n";
		    return -1;
		}

		if (theDomain->addNode(theNode) == false) {
		    opserr << "WARNING failed to add node to the domain\n";
		    opserr << "node: " << nodeID << endln;
		    delete theNode;
		    return -1;
		}
		nodeID++;
	    }
	}
    }

    // create the elements: numX*numY*numZ elements to be created if 4 node elements
    //                      numX/2*numY/2*numZ/2 nodes to be created if 9 node elements
    int eleID = idata[4];
    Element* theEle = 0;

    NDMaterial* mat = OPS_getNDMaterial(matTag);
    if (mat == 0) {
	opserr << "WARNING material not found\n";
	opserr << "Material: " << matTag << "\n";
	return -1;
    }



    for (int k=0; k<idata[2]; k++) {
	for (int j=0; j<idata[1]; j++) {
	    for (int i=0; i<idata[0]; i++) {

		const ID& nodeTags = theBlock.getElementNodes(i,j,k);
		int nd1 = nodeTags(0) + idata[3];
		int nd2 = nodeTags(1) + idata[3];
		int nd3 = nodeTags(2) + idata[3];
		int nd4 = nodeTags(3) + idata[3];
		int nd5 = nodeTags(4) + idata[3];
		int nd6 = nodeTags(5) + idata[3];
		int nd7 = nodeTags(6) + idata[3];
		int nd8 = nodeTags(7) + idata[3];

		if (strcmp(type, "stdBrick") == 0) {

		    theEle = new Brick(eleID,nd1,nd2,nd3,nd4,nd5,nd6,nd7,nd8,
				       *mat,0.,0.,0.);


		} else if (strcmp(type, "bbarBrick") == 0) {
		    theEle = new BbarBrick(eleID,nd1,nd2,nd3,nd4,nd5,nd6,nd7,nd8,
					   *mat,0.,0.,0.);

		} else if (strcmp(type, "SSPbrick") == 0 || strcmp(type,"SSPBrick") == 0) {

		    theEle = new SSPbrick(eleID,nd1,nd2,nd3,nd4,nd5,nd6,nd7,nd8,
					   *mat,0.,0.,0.);

		} else {
		    opserr << "WARNING element type " << type << " is currently unknown by this command.\n";
		    return -1;
		}

		if (theDomain->addElement(theEle) == false) {
		    opserr<<"WARNING failed to add element to domain\n";
		    delete theEle;
		    return -1;
		}

		eleID++;
	    }
	}
    }

    return 0;
}
