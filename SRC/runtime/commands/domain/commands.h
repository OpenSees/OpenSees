//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//

// domain/node.cpp
Tcl_CmdProc nodeCoord;
Tcl_CmdProc nodeDOFs;
Tcl_CmdProc nodeMass;
Tcl_CmdProc nodePressure;
Tcl_CmdProc nodeDisp;
Tcl_CmdProc nodeReaction;
Tcl_CmdProc nodeUnbalance;
Tcl_CmdProc nodeEigenvector;
Tcl_CmdProc setNodeCoord;
Tcl_CmdProc nodeRotation;

// domain/region.cpp
Tcl_CmdProc TclCommand_addMeshRegion;

// domain/recorder.cpp
Tcl_CmdProc OPS_recorderValue;


// domain/section.cpp
Tcl_CmdProc sectionForce;
Tcl_CmdProc sectionDeformation;
Tcl_CmdProc sectionStiffness;
Tcl_CmdProc sectionFlexibility;
Tcl_CmdProc sectionLocation;
Tcl_CmdProc sectionWeight;
Tcl_CmdProc sectionTag;
Tcl_CmdProc sectionDisplacement;



namespace OpenSees {
namespace DomainCommands {
  // domain.cpp
  Tcl_ObjCmdProc removeObject;
  Tcl_ObjCmdProc fixedNodes;
  Tcl_ObjCmdProc constrainedNodes;
  Tcl_ObjCmdProc fixedDOFs;
  Tcl_ObjCmdProc constrainedDOFs;
  Tcl_ObjCmdProc domainChange;
  Tcl_CmdProc    retainedDOFs;
  Tcl_CmdProc    updateElementDomain;

  // domain/element.cpp
  Tcl_CmdProc addElementRayleigh;
  Tcl_CmdProc getEleTags;
  Tcl_CmdProc getNumElements;
  Tcl_CmdProc getEleClassTags;
  Tcl_CmdProc eleNodes;
  Tcl_CmdProc eleType;
  Tcl_CmdProc eleForce;
  Tcl_CmdProc localForce;
  Tcl_CmdProc eleDynamicalForce;
  Tcl_CmdProc eleResponse;

  // modal.cpp
  Tcl_CmdProc modalProperties;
}
}

// parameter.cpp
Tcl_CmdProc getParamTags;
Tcl_CmdProc TclCommand_parameter;
Tcl_CmdProc getParamValue;
Tcl_CmdProc TclCommand_setParameter;


//


// sensitivity.cpp
Tcl_CmdProc computeGradients;
Tcl_CmdProc sensNodeDisp;
Tcl_CmdProc sensLambda; // Abbas
Tcl_CmdProc sensNodeVel;
Tcl_CmdProc sensNodeAccel;
Tcl_CmdProc sensNodePressure;
Tcl_CmdProc sensSectionForce;
Tcl_CmdProc TclCommand_sensitivityAlgorithm;
// Tcl_CmdProc sensitivityIntegrator;


// Tcl_CmdProc startTimer;
// Tcl_CmdProc stopTimer;
Tcl_CmdProc TclCommand_getTime;
Tcl_CmdProc TclCommand_setTime;

Tcl_CmdProc rayleighDamping;

Tcl_CmdProc modalDamping;

Tcl_CmdProc modalDampingQ;

Tcl_CmdProc basicDeformation;

Tcl_CmdProc basicForce;

Tcl_CmdProc basicStiffness;

// added: Chris McGann, U.Washington for initial state analysis of nDMaterials
Tcl_CmdProc InitialStateAnalysis;



Tcl_CmdProc setLoadConst;

Tcl_CmdProc setCreep;

Tcl_CmdProc getLoadFactor;

Tcl_CmdProc printModelGID;

Tcl_CmdProc TclAddRecorder;

Tcl_CmdProc addAlgoRecorder;

Tcl_CmdProc addDatabase;

Tcl_CmdProc playbackRecorders;

Tcl_CmdProc playbackAlgorithmRecorders;

Tcl_CmdProc groundExcitation;


Tcl_CmdProc getEleLoadData;
Tcl_CmdProc getEleLoadClassTags;
Tcl_CmdProc getEleLoadTags;


Tcl_CmdProc findID;



//
Tcl_CmdProc nodeBounds;

Tcl_CmdProc nodeVel;

Tcl_CmdProc setNodeVel;

Tcl_CmdProc setNodeDisp;

Tcl_CmdProc setNodeAccel;

Tcl_CmdProc nodeAccel;

Tcl_CmdProc nodeResponse;

Tcl_CmdProc calculateNodalReactions;

Tcl_CmdProc getNodeTags;
Tcl_CmdProc retainedNodes;