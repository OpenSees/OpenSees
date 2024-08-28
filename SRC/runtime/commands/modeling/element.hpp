#include <string>
#include <unordered_map>

class G3_Runtime;
typedef void *OPS_Routine(G3_Runtime* , int, const char** const);

extern OPS_Routine OPS_ComponentElement2d;
extern OPS_Routine OPS_ComponentElement3d;
// extern  void *OPS_ComponentElementDamp2d(G3_Runtime*);
extern OPS_Routine OPS_TrussElement;
extern OPS_Routine OPS_TrussSectionElement;
extern OPS_Routine OPS_CorotTrussElement;
extern OPS_Routine OPS_CorotTrussSectionElement;
extern OPS_Routine OPS_ElasticTubularJoint;
extern OPS_Routine OPS_ZeroLength;
extern OPS_Routine OPS_ZeroLengthContactNTS2D;
extern OPS_Routine OPS_ZeroLengthVG_HG;
extern OPS_Routine OPS_ZeroLengthInterface2D;
extern OPS_Routine OPS_ZeroLengthImpact3D;
extern OPS_Routine OPS_ZeroLengthContactASDimplex; 
extern "C" OPS_Routine OPS_PY_Macro2D;
extern OPS_Routine OPS_SimpleContact2D;
extern OPS_Routine OPS_SimpleContact3D;

extern OPS_Routine OPS_SurfaceLoad;
extern OPS_Routine OPS_TriSurfaceLoad;

extern OPS_Routine OPS_ModElasticBeam2d;
extern OPS_Routine OPS_ModElasticBeam3d;
extern OPS_Routine OPS_ElasticTimoshenkoBeam2d;
extern OPS_Routine OPS_ElasticTimoshenkoBeam3d;
extern OPS_Routine OPS_AxEqDispBeamColumn2d;
extern OPS_Routine OPS_BeamGT;
// extern void* OPS_GradientInelasticBeamColumn2d();
// extern void* OPS_GradientInelasticBeamColumn3d();
extern OPS_Routine OPS_DispBeamColumnAsym3dTcl;  // Xinlong Du
extern OPS_Routine OPS_MixedBeamColumnAsym3dTcl; // Xinlong Du
#if defined(_HAVE_LHNMYS) || defined(OPSDEF_ELEMENT_LHNMYS)
  extern void *OPS_BeamColumn2DwLHNMYS(G3_Runtime*);
  extern void *OPS_Beam2dDamage(G3_Runtime*);
  extern void *OPS_BeamColumn2DwLHNMYS_Damage(G3_Runtime*);
  extern void *OPS_BeamColumn3DwLHNMYS(G3_Runtime*);
#endif

extern OPS_Routine OPS_FourNodeTetrahedron;
extern OPS_Routine OPS_TPB1D;
extern OPS_Routine OPS_TFP_Bearing;
extern OPS_Routine OPS_FPBearingPTV;
extern OPS_Routine OPS_MultiFP2d;
extern OPS_Routine OPS_CoupledZeroLength;
extern OPS_Routine OPS_FourNodeQuad3d;
extern OPS_Routine OPS_Quad4FiberOverlay;
extern OPS_Routine OPS_QuadBeamEmbedContact;
extern OPS_Routine OPS_ASID8QuadWithSensitivity;
extern OPS_Routine OPS_AV3D4QuadWithSensitivity;

extern OPS_Routine OPS_Brick8FiberOverlay;
extern OPS_Routine OPS_TripleFrictionPendulum;
extern OPS_Routine OPS_TripleFrictionPendulumX;
extern OPS_Routine OPS_Truss2;
extern OPS_Routine OPS_PML3D;
extern OPS_Routine OPS_PML2D;
extern OPS_Routine OPS_CorotTruss2;
extern OPS_Routine OPS_HDR;
extern OPS_Routine OPS_LeadRubberX;
extern OPS_Routine OPS_ElastomericX;
extern OPS_Routine OPS_N4BiaxialTruss;
extern OPS_Routine OPS_AC3D8HexWithSensitivity;
extern OPS_Routine OPS_VS3D4WuadWithSensitivity;
extern OPS_Routine OPS_MVLEM;        // Kristijan Kolozvari
extern OPS_Routine OPS_SFI_MVLEM;    // Kristijan Kolozvari
extern OPS_Routine OPS_MVLEM_3D;     // Kristijan Kolozvari
extern OPS_Routine OPS_SFI_MVLEM_3D; // Kristijan Kolozvari
extern OPS_Routine OPS_E_SFI_MVLEM_3D;
extern OPS_Routine OPS_E_SFI;
extern OPS_Routine OPS_MEFI;
extern OPS_Routine OPS_ElastomericBearingBoucWenMod3d;
extern OPS_Routine OPS_InertiaTrussElement; // Added by Xiaodong Ji, Yuhao Cheng, Yue Yu
extern OPS_Routine OPS_CatenaryCableElement;
extern OPS_Routine OPS_LysmerTriangle;
extern OPS_Routine OPS_ASDEmbeddedNodeElement;     // Massimo Petracca (ASDEA)
extern OPS_Routine OPS_ASDAbsorbingBoundary2D;     // Massimo Petracca (ASDEA)
extern OPS_Routine OPS_ASDAbsorbingBoundary3D;     // Massimo Petracca (ASDEA)
extern OPS_Routine OPS_FSIInterfaceElement2D;      // Massimo Petracca (ASDEA)
extern OPS_Routine OPS_FSIFluidBoundaryElement2D;  // Massimo Petracca (ASDEA)
extern OPS_Routine OPS_FSIFluidElement2D;          // Massimo Petracca (ASDEA)
extern OPS_Routine OPS_ASDShellT3;
extern OPS_Routine OPS_TwoNodeLink;
extern OPS_Routine OPS_LinearElasticSpring;
extern OPS_Routine OPS_Inerter;
extern OPS_Routine OPS_Inno3DPnPJoint;
extern OPS_Routine OPS_Adapter;
extern OPS_Routine OPS_Actuator;
extern OPS_Routine OPS_ActuatorCorot;
extern OPS_Routine OPS_ElastomericBearingPlasticity2d;
extern OPS_Routine OPS_ElastomericBearingPlasticity3d;
extern OPS_Routine OPS_ElastomericBearingBoucWen2d;
extern OPS_Routine OPS_ElastomericBearingBoucWen3d;
extern OPS_Routine OPS_ElastomericBearingUFRP2d;

extern OPS_Routine OPS_RJWatsonEQS2d;
extern OPS_Routine OPS_RJWatsonEQS3d;
extern OPS_Routine OPS_RockingBC;
extern OPS_Routine OPS_LehighJoint2d;
extern OPS_Routine OPS_MasonPan12;
extern OPS_Routine OPS_MasonPan3D;

#include <algorithm>
#include <string>
static
std::string toLower( const std::string & s )
{
    std::string copy = s;
    transform( copy.begin( ), copy.end( ), copy.begin( ), 
        [](unsigned char c) { return std::tolower(c); });
    return copy;
}

static bool 
equalsIgnoreCase( const std::string & lhs, const std::string & rhs )
{
    return toLower( lhs ) == toLower( rhs );
}

class CaseInsensitive
{
  public:
    size_t operator( ) ( const std::string & s ) const
    {  
        static std::hash<std::string> hf;
        return hf( toLower( s ) );
    }
    
    bool operator( ) ( const std::string & lhs, const std::string & rhs ) const
    {
        return equalsIgnoreCase( lhs, rhs );
    }
};

#if 0
// template <OPS_Routine fn2d, OPS_Routine fn3d> static int
class BasicModelBuilder;
struct Tcl_Interp;
class Domain;
template <OPS_Routine fn> static int
dispatch(BasicModelBuilder* builder, Tcl_Interp* interp, int argc, const char** const argv)
{
  int ndm = builder->getNDM();
  Domain* domain = builder->getDomain();
  G3_Runtime *rt = G3_getRuntime(interp);

  Element* theElement = (Element*)fn( rt, argc, argv );

  if (domain->addElement(*theElement) == false) {
    opserr << G3_ERROR_PROMPT << "Could not add element to the domain.\n";
    delete theElement;
    return TCL_ERROR;
  }
  return TCL_OK;
}
#endif
Tcl_CmdProc TclCommand_addTwoNodeLink;
// Plane
Tcl_CmdProc TclBasicBuilder_addFourNodeQuad;
Tcl_CmdProc TclBasicBuilder_addFourNodeQuadWithSensitivity;
Tcl_CmdProc TclBasicBuilder_addEnhancedQuad;
Tcl_CmdProc TclBasicBuilder_addConstantPressureVolumeQuad;
Tcl_CmdProc TclBasicBuilder_addNineNodeMixedQuad;
Tcl_CmdProc TclBasicBuilder_addNineNodeQuad;
Tcl_CmdProc TclBasicBuilder_addSixNodeTri;
Tcl_CmdProc TclBasicBuilder_addEightNodeQuad;
Tcl_CmdProc TclBasicBuilder_addFourNodeQuadUP;
Tcl_CmdProc TclBasicBuilder_addNineFourNodeQuadUP;
Tcl_CmdProc TclBasicBuilder_addBBarFourNodeQuadUP;
Tcl_CmdProc TclDispatch_newTri31;
// Brick
Tcl_CmdProc TclBasicBuilder_addBrickUP;
Tcl_CmdProc TclBasicBuilder_addBBarBrickUP;
Tcl_CmdProc TclBasicBuilder_addTwentyEightNodeBrickUP;
Tcl_CmdProc TclBasicBuilder_addTwentyNodeBrick;
Tcl_CmdProc TclBasicBuilder_addBrick;
Tcl_CmdProc TclCommand_SSP_Element;

Tcl_CmdProc TclCommand_addActuator;
Tcl_CmdProc TclCommand_addActuatorCorot;
Tcl_CmdProc TclCommand_addAdapter;
Tcl_CmdProc TclBasicBuilder_addRJWatsonEqsBearing;

static
std::unordered_map<std::string, Tcl_CmdProc *, CaseInsensitive, CaseInsensitive> 
element_dispatch_tcl = {
  {"twoNodeLink",               TclCommand_addTwoNodeLink},
//
// Plane
//
  {"Quad",                      TclBasicBuilder_addFourNodeQuad},
  {"stdQuad",                   TclBasicBuilder_addFourNodeQuad},

  {"quadWithSensitivity",       TclBasicBuilder_addFourNodeQuadWithSensitivity},

  {"enhancedQuad",              TclBasicBuilder_addEnhancedQuad},

  {"bbarQuad",                  TclBasicBuilder_addConstantPressureVolumeQuad},
  {"mixedQuad",                 TclBasicBuilder_addConstantPressureVolumeQuad},

  {"nineNodeMixedQuad",         TclBasicBuilder_addNineNodeMixedQuad},
  {"nineNodeQuad",              TclBasicBuilder_addNineNodeMixedQuad},

  {"quad9n",                    TclBasicBuilder_addNineNodeQuad},

  {"quad8n",                    TclBasicBuilder_addEightNodeQuad},

  {"tri6n",                     TclBasicBuilder_addSixNodeTri},
  {"tri31",                     TclDispatch_newTri31},
  {"quadUP",                    TclBasicBuilder_addFourNodeQuadUP},

  {"9_4_QuadUP",                TclBasicBuilder_addNineFourNodeQuadUP},

  {"bbarQuadUP",                TclBasicBuilder_addBBarFourNodeQuadUP},
//
// Brick
//
  {"BrickUP",                   TclBasicBuilder_addBrickUP},

  {"20_8_BrickUP",              TclBasicBuilder_addTwentyEightNodeBrickUP},

  {"20NodeBrick",               TclBasicBuilder_addTwentyNodeBrick},

  {"bbarBrickUP",               TclBasicBuilder_addBBarBrickUP},

  {"stdBrick",                  TclBasicBuilder_addBrick},
  {"bbarBrick",                 TclBasicBuilder_addBrick},
  {"bbarBrickWithSensitivity",  TclBasicBuilder_addBrick},
  {"flBrick",                   TclBasicBuilder_addBrick},
  
  {"SSPquad",                   TclCommand_SSP_Element},
  {"SSPquadUP",                 TclCommand_SSP_Element},
  {"SSPbrick",                  TclCommand_SSP_Element},

// Actuators
  {"actuator",                  TclCommand_addActuator},
  {"corotActuator",             TclCommand_addActuatorCorot},
  {"adapter",                   TclCommand_addAdapter},

// Bearing
  {"RJWatsonEqsBearing",        TclBasicBuilder_addRJWatsonEqsBearing},
  {"RJWatsonBearing",           TclBasicBuilder_addRJWatsonEqsBearing},
  {"EQSBearing",                TclBasicBuilder_addRJWatsonEqsBearing},
};

static
std::unordered_map<std::string, OPS_Routine *, CaseInsensitive, CaseInsensitive> 
element_dispatch = {
// Truss
  {"TrussSection",                 OPS_TrussSectionElement},
  {"CorotTrussSection",            OPS_CorotTrussSectionElement},
  {"N4BiaxialTruss",               OPS_N4BiaxialTruss},
  {"Truss2",                       OPS_Truss2},
  {"CorotTruss2",                  OPS_CorotTruss2},
  {"InertiaTruss",                 OPS_InertiaTrussElement},


// Shell
  {"ASDShellT3",                   OPS_ASDShellT3},

// Point
  {"zeroLengthContactNTS2D",       OPS_ZeroLengthContactNTS2D},
  {"zeroLengthInterface2D",        OPS_ZeroLengthInterface2D},
  {"zeroLengthImpact3D",           OPS_ZeroLengthImpact3D},

  {"componentElement2d",           OPS_ComponentElement2d},
  {"componentElement3d",           OPS_ComponentElement3d},

#if 0
  {"componentElementDamp2d", OPS_ComponentElementDamp2d},
#endif

  {"ModElasticBeam2d",             OPS_ModElasticBeam2d},
  {"ModElasticBeam3d",             OPS_ModElasticBeam3d},


  {"FPBearingPTV",                 OPS_FPBearingPTV},

  {"TripleFrictionPendulum",       OPS_TripleFrictionPendulum},
  {"TripleFrictionPendulumX",      OPS_TripleFrictionPendulumX},
  {"HDR",                          OPS_HDR},
  {"LeadRubberX",                  OPS_LeadRubberX},
  {"ElastomericX",                 OPS_ElastomericX},

  {"AxEqDispBeamColumn2d",         OPS_AxEqDispBeamColumn2d},

// MVLEM
  {"MVLEM",                        OPS_MVLEM},        // Kristijan Kolozvari
  {"SFI_MVLEM",                    OPS_SFI_MVLEM},    // Kristijan Kolozvari
  {"MVLEM_3D",                     OPS_MVLEM_3D},     // Kristijan Kolozvari
  {"SFI_MVLEM_3D",                 OPS_SFI_MVLEM_3D}, // Kristijan Kolozvari
  {"E_SFI_MVLEM_3D",               OPS_E_SFI_MVLEM_3D},
  {"E_SFI",                        OPS_E_SFI},
  {"MEFI",                         OPS_MEFI},

// Fluid
  {"FSIFluidElement2D",            OPS_FSIFluidElement2D },
  {"FSIInterfaceElement2D",        OPS_FSIInterfaceElement2D },
  {"FSIFluidBoundaryElement2D",    OPS_FSIFluidBoundaryElement2D },

// Joint
  {"ElasticTubularJoint",          OPS_ElasticTubularJoint},
  {"Inno3DPnPJoint",               OPS_Inno3DPnPJoint},

// Other
  {"MasonPan12",                   OPS_MasonPan12},
  {"MasonPan3D",                   OPS_MasonPan3D},
  {"BeamGT",                       OPS_BeamGT},
  {"ZeroLengthVG_HG",              OPS_ZeroLengthVG_HG},
  {"ZeroLengthContactASDimplex",   OPS_ZeroLengthContactASDimplex},
//{"twoNodeLink",                  OPS_TwoNodeLink},
  {"SurfaceLoad",                  OPS_SurfaceLoad},
  {"TriSurfaceLoad",               OPS_TriSurfaceLoad},
  {"TPB1D",                        OPS_TPB1D},
  {"quad3d",                       OPS_FourNodeQuad3d},
  {"AC3D8",                        OPS_AC3D8HexWithSensitivity},
  {"ASI3D8",                       OPS_ASID8QuadWithSensitivity},
  {"AV3D4",                        OPS_AV3D4QuadWithSensitivity},
  {"ElastomericBearingBoucWenMod", OPS_ElastomericBearingBoucWenMod3d},
  {"VS3D4",                        OPS_VS3D4WuadWithSensitivity},
  {"CatenaryCable",                OPS_CatenaryCableElement},
  {"ASDEmbeddedNodeElement",       OPS_ASDEmbeddedNodeElement},
  {"LysmerTriangle",               OPS_LysmerTriangle},
  {"ASDAbsorbingBoundary2D",       OPS_ASDAbsorbingBoundary2D},
  {"ASDAbsorbingBoundary3D",       OPS_ASDAbsorbingBoundary3D},
  {"FourNodeTetrahedron",          OPS_FourNodeTetrahedron},
  {"LinearElasticSpring",          OPS_LinearElasticSpring},
  {"Inerter",                      OPS_Inerter},
  {"Adapter",                      OPS_Adapter},
  {"Actuator",                     OPS_Actuator},
  {"CorotActuator",                OPS_ActuatorCorot},
  {"RockingBC",                    OPS_RockingBC},
  {"LehighJoint2D",                OPS_LehighJoint2d},
};

