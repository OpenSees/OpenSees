//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#include <tcl.h>

// modeling/model.cpp
extern Tcl_CmdProc  TclCommand_wipeModel;
extern Tcl_CmdProc  buildModel;

// modeling/nodes.cpp
extern Tcl_CmdProc  TclCommand_getNDM;
extern Tcl_CmdProc  TclCommand_getNDF;
extern Tcl_CmdProc  TclCommand_addNode;
extern Tcl_CmdProc  TclCommand_addNodalMass;
extern Tcl_CmdProc  TclCommand_addNodalLoad;
// 
extern Tcl_CmdProc  TclCommand_addSeries;
extern Tcl_CmdProc  TclCommand_addPattern;
extern Tcl_CmdProc  TclCommand_addTimeSeries;
extern Tcl_CmdProc  TclCommand_addGeomTransf;

// element.cpp
extern Tcl_CmdProc  TclCommand_addElement;

// blockND.cpp
extern Tcl_CmdProc  TclCommand_doBlock2D;
extern Tcl_CmdProc  TclCommand_doBlock3D;

// uniaxial.cpp
extern Tcl_CmdProc  TclCommand_addUniaxialMaterial;

// section.cpp
extern Tcl_CmdProc  TclCommand_addSection;
extern Tcl_CmdProc  TclCommand_addPatch;
extern Tcl_CmdProc  TclCommand_addReinfLayer;
// extern Tcl_CmdProc  TclCommand_addRemoFiber;
extern Tcl_CmdProc  TclCommand_addFiber;
extern Tcl_CmdProc  TclCommand_addHFiber;

// Constraints
extern Tcl_CmdProc TclCommand_addMP;
extern Tcl_CmdProc TclCommand_addSP;
extern Tcl_CmdProc TclCommand_addHomogeneousBC;
extern Tcl_CmdProc TclCommand_addHomogeneousBC_X;
extern Tcl_CmdProc TclCommand_addHomogeneousBC_Y; 
extern Tcl_CmdProc TclCommand_addHomogeneousBC_Z;
extern Tcl_CmdProc TclCommand_addEqualDOF_MP;
extern Tcl_CmdProc TclCommand_addEqualDOF_MP_Mixed;
extern Tcl_CmdProc TclCommand_RigidLink;
extern Tcl_CmdProc TclCommand_RigidDiaphragm;
extern Tcl_CmdProc TclCommand_addImposedMotionSP;
extern Tcl_CmdProc TclCommand_addGroundMotion;

// Loads
// extern Tcl_CmdProc  TclCommand_addNodalLoad;
Tcl_CmdProc TclCommand_addElementalLoad;

// Other
extern Tcl_CmdProc  TclCommand_addHystereticBackbone;
extern Tcl_CmdProc  TclCommand_updateMaterialStage;

// UpdatedLagrange
Tcl_CmdProc TclCommand_addCyclicModel;

Tcl_CmdProc TclCommand_addParameter;
Tcl_CmdProc TclCommand_mesh;
Tcl_CmdProc TclCommand_remesh;
Tcl_CmdProc TclCommand_backgroundMesh; 
Tcl_CmdProc TclCommand_addBeamIntegration;

//
Tcl_CmdProc TclCommand_addFrictionModel;
Tcl_CmdProc TclCommand_addLimitCurve;
Tcl_CmdProc TclCommand_addNDMaterial;

// invoking.cpp
Tcl_CmdProc TclCommand_invoke;

// printing.cpp
Tcl_CmdProc TclCommand_print;
Tcl_CmdProc TclCommand_classType;

struct char_cmd {
  const char* name;
  Tcl_CmdProc*  func;
}  const tcl_char_cmds[] =  {
  {"build",                buildModel},

  {"getNDM",               TclCommand_getNDM},
  {"getNDF",               TclCommand_getNDF},
  {"node",                 TclCommand_addNode},
  {"mass",                 TclCommand_addNodalMass},
  {"element",              TclCommand_addElement},

  {"print",                TclCommand_print},
  {"classType",            TclCommand_classType},
  {"printModel",           TclCommand_print},

  {"fix",                  TclCommand_addHomogeneousBC},
  {"fixX",                 TclCommand_addHomogeneousBC_X},
  {"fixY",                 TclCommand_addHomogeneousBC_Y},
  {"fixZ",                 TclCommand_addHomogeneousBC_Z},

// //
  {"with",                 TclCommand_invoke},
  {"invoke",               TclCommand_invoke},
// Materials & sections
  {"uniaxialMaterial",     TclCommand_addUniaxialMaterial},
  {"nDMaterial",           TclCommand_addNDMaterial},
  {"beamIntegration",      TclCommand_addBeamIntegration},

  {"section",              TclCommand_addSection},
  {"patch",                TclCommand_addPatch},
  {"fiber",                TclCommand_addFiber},
  {"layer",                TclCommand_addReinfLayer},
  {"Hfiber",               TclCommand_addHFiber},

  {"geomTransf",           TclCommand_addGeomTransf},
  {"transform",            TclCommand_addGeomTransf},

  {"pattern",              TclCommand_addPattern},
//   {"load",             TclCommand_addNodalLoad},
  {"nodalLoad",            TclCommand_addNodalLoad},
  {"timeSeries",           TclCommand_addTimeSeries},

  {"equalDOF",             TclCommand_addEqualDOF_MP},
  {"rigidLink",            TclCommand_RigidLink},
  
  {"sp",                   TclCommand_addSP},
  {"groundMotion",         TclCommand_addGroundMotion},
  {"imposedMotion",        TclCommand_addImposedMotionSP},
  {"imposedSupportMotion", TclCommand_addImposedMotionSP},

  {"eleLoad",              TclCommand_addElementalLoad},

  {"block2D",              TclCommand_doBlock2D},
  {"block3D",              TclCommand_doBlock3D},
  {"rigidDiaphragm",       &TclCommand_RigidDiaphragm},

/*
  {"mp",                   TclCommand_addMP},

  {"equalDOF_Mixed",       TclCommand_addEqualDOF_MP_Mixed},
  {"PySimple1Gen",         TclCommand_doPySimple1Gen},
  {"TzSimple1Gen",         TclCommand_doTzSimple1Gen},
  {"ShallowFoundationGen", BasicModelBuilder_doShallowFoundationGen},
*/

// OTHER OBJECT TYPES
  {"hystereticBackbone",   TclCommand_addHystereticBackbone},
  {          "backbone",   TclCommand_addHystereticBackbone},

  {"frictionModel",        TclCommand_addFrictionModel},

  {"cyclicModel",          TclCommand_addCyclicModel},
#if 0
  {"yieldSurface_BC",      TclCommand_addYieldSurface_BC},
  {"ysEvolutionModel",     TclCommand_addYS_EvolutionModel},
  {"plasticMaterial",      TclCommand_addYS_PlasticMaterial},
  {"limitCurve",           TclCommand_addLimitCurve},
  {"damageModel",          TclCommand_addDamageModel},
  {"stiffnessDegradation", TclCommand_addStiffnessDegradation},
  {"unloadingRule",        TclCommand_addUnloadingRule},
  {"strengthDegradation",  TclCommand_addStrengthDegradation},
  {"loadPackage",          TclCommand_Package},
#endif


// command for elast2plast in Multi-yield plasticity, by ZHY
  {"updateMaterialStage", TclCommand_updateMaterialStage},
#if 0
  {"updateMaterials",     TclCommand_UpdateMaterials},
#endif

#if 0
// command for updating properties of soil materials, by ZHY
   {"updateParameter", TclCommand_UpdateParameter},
#endif

};

Tcl_CmdProc TclCommand_Package;

// Added by Scott J. Brandenberg
Tcl_CmdProc TclCommand_doPySimple1Gen;
Tcl_CmdProc TclCommand_doTzSimple1Gen;

// Added by Prishati Raychowdhury (UCSD)
Tcl_CmdProc BasicModelBuilder_doShallowFoundationGen;
// End PRC
//Leo
Tcl_CmdProc BasicModelBuilder_addRemoHFiber;
Tcl_CmdProc TclCommand_addStiffnessDegradation;
Tcl_CmdProc TclCommand_addUnloadingRule;
Tcl_CmdProc TclCommand_addStrengthDegradation;
/// added by ZHY
Tcl_CmdProc TclCommand_UpdateMaterials;
Tcl_CmdProc TclCommand_UpdateParameter;
Tcl_CmdProc TclCommand_addElementRayleigh;

// Added by Alborz Ghofrani - U.Washington
Tcl_CmdProc TclCommand_GenerateInterfacePoints;

