//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#include <tcl.h>
#include "element.hpp"
#include <assert.h>
#include <stdlib.h>

#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#define strcmp strcasecmp

#include <runtimeAPI.h>
#include <BasicModelBuilder.h>

#include <OPS_Stream.h>
#include <G3_Logging.h>
#include <packages.h>
#include <Domain.h>
#include <Element.h>
#include <CrdTransf.h>
#include <NDMaterial.h>


#include <UniaxialMaterial.h>
#include <MultipleShearSpring.h>
#include <MultipleNormalSpring.h>
#include <KikuchiBearing.h>
#include <YamamotoBiaxialHDR.h>
#include <WheelRail.h>


typedef struct elementPackageCommand {
  char *funcName;
  void *(*funcPtr)();
  struct elementPackageCommand *next;
} ElementPackageCommand;

static ElementPackageCommand *theElementPackageCommands = nullptr;

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp, int cArg,
                          int mArg, TCL_Char ** const argv, Domain *);

//
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

#if 0 // cmp - commented out to eliminate use of TclBasicBuilder
extern int TclBasicBuilder_addFeapTruss(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv, Domain *, TclBasicBuilder *, int argStart);
extern int Tcl_addWrapperElement(eleObj *, ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv, Domain *, TclBuilder *);
// Added by Quan Gu and Yongdou Liu, et al. on 2018/10/31 (Xiamen University)
#endif
static Tcl_CmdProc TclBasicBuilder_addWheelRail;



extern OPS_Routine OPS_ElasticBeam3d;
extern void *OPS_ElasticBeam2d(G3_Runtime *, const ID &);


// Frame
Tcl_CmdProc TclBasicBuilder_addElasticBeam;
Tcl_CmdProc TclBasicBuilder_addGradientInelasticBeamColumn;
Tcl_CmdProc TclBasicBuilder_addForceBeamColumn;

// Zero-length
Tcl_CmdProc TclCommand_addZeroLength;
Tcl_CmdProc TclCommand_addZeroLengthSection;
Tcl_CmdProc TclCommand_addZeroLengthContact2D;
Tcl_CmdProc TclCommand_addZeroLengthContact3D;
Tcl_CmdProc TclCommand_addZeroLengthRocking;
Tcl_CmdProc TclCommand_addZeroLengthND;

Tcl_CmdProc TclBasicBuilder_addBeamWithHinges;
Tcl_CmdProc TclBasicBuilder_addDispBeamColumnInt;

// Joint
Tcl_CmdProc TclBasicBuilder_addJoint2D;
Tcl_CmdProc TclBasicBuilder_addJoint3D;
Tcl_CmdProc TclBasicBuilder_addBeamColumnJoint;

// Other
Tcl_CmdProc TclBasicBuilder_addElement2dYS;
Tcl_CmdProc TclBasicBuilder_addElastic2dGNL;
Tcl_CmdProc TclBasicBuilder_addKikuchiBearing;


Tcl_CmdProc TclBasicBuilder_addGenericCopy;
Tcl_CmdProc TclBasicBuilder_addGenericClient;

Tcl_CmdProc TclCommand_addFlatSliderBearing;
Tcl_CmdProc TclCommand_addSingleFPBearing;

class TclBasicBuilder;
typedef int (G3_TclElementCommand)(ClientData, Tcl_Interp*, int, const char** const, Domain*, TclBasicBuilder*);
G3_TclElementCommand TclBasicBuilder_addMultipleShearSpring;
G3_TclElementCommand TclBasicBuilder_addMultipleNormalSpring;
G3_TclElementCommand TclBasicBuilder_addYamamotoBiaxialHDR;
G3_TclElementCommand TclBasicBuilder_addMasonPan12;
G3_TclElementCommand TclBasicBuilder_addMasonPan3D;
G3_TclElementCommand TclBasicBuilder_addBeamGT;




// Shells
Element* TclDispatch_newASDShellQ4(ClientData, Tcl_Interp*, int, TCL_Char** const);
Element* TclDispatch_newShellANDeS(ClientData, Tcl_Interp*, int, TCL_Char** const);
Element* TclDispatch_newShellDKGQ(ClientData, Tcl_Interp*, int, TCL_Char** const);
Element* TclDispatch_newShellDKGT(ClientData, Tcl_Interp*, int, TCL_Char** const);
Element* TclDispatch_newShellMITC4(ClientData, Tcl_Interp*, int, TCL_Char** const);
Element* TclDispatch_newShellMITC4Thermal(ClientData, Tcl_Interp*, int, TCL_Char** const);
Element* TclDispatch_newShellMITC9(ClientData, Tcl_Interp*, int, TCL_Char** const);
Element* TclDispatch_newShellNLDKGQ(ClientData, Tcl_Interp*, int, TCL_Char** const);
Element* TclDispatch_newShellNLDKGQThermal(ClientData, Tcl_Interp*, int, TCL_Char** const);
Element* TclDispatch_newShellNLDKGT(ClientData, Tcl_Interp*, int, TCL_Char** const);



int
TclCommand_addElement(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  TclBasicBuilder *theTclBuilder = (TclBasicBuilder*)G3_getSafeBuilder(rt);

  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  Domain *theTclDomain = builder->getDomain();

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theTclDomain);

  // check at least two arguments so don't segemnt fault on strcmp
  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "insufficient arguments, expected:\n";
    opserr << "      element eleType <specific element args> .. \n";
    return TCL_ERROR;
  }

  void* theEle = nullptr;
  Element *theElement = nullptr;
  int ndm = builder->getNDM();

  auto cmd = element_dispatch.find(std::string(argv[1]));
  if (cmd != element_dispatch.end()) {
    theEle = (*cmd->second)(rt, argc, &argv[0]);
  }

  // Try Tcl element library
  auto tcl_cmd = element_dispatch_tcl.find(std::string(argv[1]));
  if (tcl_cmd != element_dispatch_tcl.end()) {
    return (*tcl_cmd->second)(clientData, interp, argc, &argv[0]);
  }

  if (strcasecmp(argv[1], "truss") == 0) {
    theEle = OPS_TrussElement(rt, argc, argv);

    // for backward compatibility
    if (theEle == nullptr)
      theEle = OPS_TrussSectionElement(rt, argc, argv);
  }

  else if ((strcasecmp(argv[1], "elasticBeamColumn") == 0) ||
           (strcasecmp(argv[1], "elasticBeam") == 0)  ||
           (strcasecmp(argv[1], "PrismFrame") == 0)) {

    return TclBasicBuilder_addElasticBeam(clientData, interp, argc, argv);
  }

  else if (strcasecmp(argv[1], "PML") == 0) {
    if (ndm == 2)
      theEle = OPS_PML2D(rt, argc, argv);
    else
      theEle = OPS_PML3D(rt, argc, argv);

#if 0
  } else if (strcmp(argv[1], "gradientInelasticBeamColumn") == 0) {

      Element *theEle = 0;
      if (ndm == 2)
        theEle = OPS_GradientInelasticBeamColumn2d(rt, argc, argv);
      else
        theEle = OPS_GradientInelasticBeamColumn3d(rt, argc, argv);

      if (theEle != 0)
        theElement = theEle;
      else {
        return TCL_ERROR;
      }
    }
#endif

#if defined(_HAVE_LHNMYS) || defined(OPSDEF_ELEMENT_LHNMYS)
  } else if (strcmp(argv[1], "beamColumn2DwLHNMYS") == 0) {
    theEle = OPS_BeamColumn2DwLHNMYS(rt, argc, argv);

  } else if (strcmp(argv[1], "beamColumn2dDamage") == 0) {
    theEle = OPS_Beam2dDamage(rt, argc, argv);

  } else if (strcmp(argv[1], "beamColumn2DwLHNMYS_Damage") == 0) {
    theEle = OPS_BeamColumn2DwLHNMYS_Damage(rt, argc, argv);

  } else if (strcmp(argv[1], "beamColumn3DwLHNMYS") == 0) {
    theEle = OPS_BeamColumn3DwLHNMYS(rt, argc, argv);
#endif

  } else if (strcmp(argv[1], "ElasticTimoshenkoBeam") == 0) {
    if (ndm == 2)
      theEle = OPS_ElasticTimoshenkoBeam2d(rt, argc, argv);
    else
      theEle = OPS_ElasticTimoshenkoBeam3d(rt, argc, argv);
  }

  else if ((strcmp(argv[1], "pyMacro2D") == 0) ||
             (strcmp(argv[1], "PY_Macro2D") == 0)) {

    theEle = OPS_PY_Macro2D(rt, argc, argv);
  }

  else if ((strcmp(argv[1], "TFPbearing") == 0) ||
             (strcmp(argv[1], "TFP") == 0) ||
             (strcmp(argv[1], "TPFbearing") == 0) ||
             (strcmp(argv[1], "TPF") == 0)) {

    theEle = OPS_TFP_Bearing(rt, argc, argv);
  }

  else if (strcasecmp(argv[1], "CorotTruss") == 0) {
    theEle = OPS_CorotTrussElement(rt, argc, argv);

    // for backward compatibility
    if (theEle == nullptr)
      theEle = OPS_CorotTrussSectionElement(rt, argc, argv);
  }

  else if ((strcmp(argv[1], "MultiFP2d") == 0) ||
            (strcmp(argv[1], "MultiFPB2d") == 0)) {

    theEle = OPS_MultiFP2d(rt, argc, argv);
  }

// Other
  else if ((strcmp(argv[1], "CoupledZeroLength") == 0) ||
             (strcmp(argv[1], "ZeroLengthCoupled") == 0)) {
    theEle = OPS_CoupledZeroLength(rt, argc, argv);
  }

  else if (strcmp(argv[1], "ElastomericBearing") == 0 ||
          (strcmp(argv[1], "ElastomericBearingPlasticity")) == 0) {

    if (ndm == 2)
      theEle = OPS_ElastomericBearingPlasticity2d(rt, argc, argv);
    else
      theEle = OPS_ElastomericBearingPlasticity3d(rt, argc, argv);
  }

  else if (strcmp(argv[1], "ElastomericBearingBoucWen") == 0 ||
          (strcmp(argv[1], "ElastomericBearingBW")) == 0) {
    if (ndm == 2)
      theEle = OPS_ElastomericBearingBoucWen2d(rt, argc, argv);
    else
      theEle = OPS_ElastomericBearingBoucWen3d(rt, argc, argv);
  }

  else if (strcmp(argv[1], "ElastomericBearingUFRP") == 0) {
    if (ndm == 2)
      theEle = OPS_ElastomericBearingUFRP2d(rt, argc, argv);
    else {;}
      // theEle = OPS_ElastomericBearingUFRP3d(rt, argc, argv);
  }

  else if (strcmp(argv[1], "FlatSliderBearing") == 0) {
    return TclCommand_addFlatSliderBearing(clientData,
                                           interp,
                                           argc,
                                           argv);
//  if (ndm == 2)
//    theEle = OPS_FlatSliderSimple2d(rt, argc, argv);
//  else
//    theEle = OPS_FlatSliderSimple3d(rt, argc, argv);
  }

  else if (strcmp(argv[1], "SingleFPBearing") == 0 ||
          (strcmp(argv[1], "SinglePFBearing")) == 0 ||
          (strcmp(argv[1], "SFPBearing")) == 0 ||
          (strcmp(argv[1], "SPFBearing")) == 0) {
    return TclCommand_addSingleFPBearing(clientData,
                                         interp,
                                         argc,
                                         argv);
//  if (ndm == 2)
//    theEle = OPS_SingleFPSimple2d(rt, argc, argv);
//  else
//    theEle = OPS_SingleFPSimple3d(rt, argc, argv);
  }

  // Xinlong Du
  else if ((strcmp(argv[1], "DispBeamColumnAsym") == 0) ||
           (strcmp(argv[1], "DispBeamAsym")) == 0) {
    if (ndm == 3)
      theEle = OPS_DispBeamColumnAsym3dTcl(rt, argc, argv);
  }

  else if ((strcmp(argv[1], "MixedBeamColumnAsym") == 0) ||
           (strcmp(argv[1], "MixedBeamAsym") == 0)) {

    if (ndm == 3)
      theEle = OPS_MixedBeamColumnAsym3dTcl(rt, argc, argv);
  }
  // Xinlong Du

//
// Shells
//
  else if ((strcmp(argv[1], "Shell") == 0) ||
           (strcmp(argv[1], "ShellMITC4") == 0)) {
    theEle = TclDispatch_newShellMITC4(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "ShellMITC4Thermal") == 0) {
    theEle = TclDispatch_newShellMITC4Thermal(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "ShellNLDKGQThermal") == 0) {
    theEle = TclDispatch_newShellNLDKGQThermal(clientData, interp, argc, argv);
  }

  else if ((strcmp(argv[1], "ShellNL") == 0) ||
             (strcmp(argv[1], "ShellMITC9") == 0)) {
    theEle = TclDispatch_newShellMITC9(clientData, interp, argc, argv);
  }

  else if ((strcmp(argv[1], "shellDKGQ") == 0) ||
           (strcmp(argv[1], "ShellDKGQ") == 0)) {
    theEle = TclDispatch_newShellDKGQ(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "ShellNLDKGQ") == 0) {
    theEle = TclDispatch_newShellNLDKGQ(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "ShellDKGT") == 0) {
    theEle = TclDispatch_newShellDKGT(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "ShellNLDKGT") == 0) {
    theEle = TclDispatch_newShellNLDKGT(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "ASDShellQ4") == 0) {
    theEle = TclDispatch_newASDShellQ4(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "ShellANDeS") == 0) {
    theEle = TclDispatch_newShellANDeS(clientData, interp, argc, argv);
  }


  // if one of the above worked
  theElement = (Element*)theEle;

  if (theElement != nullptr) {
    if (theTclDomain->addElement(theElement) == false) {
      opserr << "WARNING could not add element of with tag: "
             << theElement->getTag()
             << " and of type: " << theElement->getClassType()
             << " to the Domain\n";
      delete theElement;
      return TCL_ERROR;
    } else
      return TCL_OK;
  }

#if 0 && defined(OPSDEF_ELEMENT_FEAP)
  if (strcmp(argv[1], "fTruss") == 0) {
    int eleArgStart = 1;
    int result = TclBasicBuilder_addFeapTruss(
        clientData, interp, argc, argv, theTclDomain, theTclBuilder, eleArgStart);
    return result;

  }
#endif // _OPS_Element_FEAP

  //
  // Beams
  //
  if (strcmp(argv[1], "dispBeamColumnInt") == 0) {
    return TclBasicBuilder_addDispBeamColumnInt(clientData, interp, argc, argv);
  } 

  else if ((strcmp(argv[1], "WheelRail") == 0)) {
    return TclBasicBuilder_addWheelRail(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "DisplFrame") == 0 ||
           strcmp(argv[1], "CubicFrame") == 0 ||
           strcmp(argv[1], "ForceFrame") == 0 ||
           strcmp(argv[1], "ExactFrame") == 0 ||
           strcmp(argv[1], "ForceDeltaFrame") == 0 ||

           strcmp(argv[1], "ForceBeamColumn") == 0 ||
           strcmp(argv[1], "DispBeamColumn") == 0 ||
           strcmp(argv[1], "DispBeamColumn") == 0 ||
           strcmp(argv[1], "TimoshenkoBeamColumn") == 0 ||
           strcmp(argv[1], "ForceBeamColumnCBDI") == 0 ||
           strcmp(argv[1], "ForceBeamColumnCSBDI") == 0 ||
           strcmp(argv[1], "ForceBeamColumnWarping") == 0 ||
           strcmp(argv[1], "ForceBeamColumnThermal") == 0 ||
           strcmp(argv[1], "ElasticForceBeamColumnWarping") == 0 ||
           strcmp(argv[1], "DispBeamColumnNL") == 0 ||
           strcmp(argv[1], "DispBeamColumnThermal") == 0 ||
           strcmp(argv[1], "ElasticForceBeamColumn") == 0 ||
           strcmp(argv[1], "NonlinearBeamColumn") == 0 ||
           strcmp(argv[1], "DispBeamColumnWithSensitivity") == 0) {

    return TclBasicBuilder_addForceBeamColumn(clientData, interp, argc, argv);

  } else if ((strstr(argv[1], "BeamWithHinges") != 0) ||
             (strcmp(argv[1], "BeamWithHinges") == 0)) {
    return TclBasicBuilder_addBeamWithHinges(clientData, interp, argc, argv);

  //
  //
//
// Brick
//
//
// Zero-Length
//
  } else if (strcmp(argv[1], "zeroLength") == 0) {
    return TclCommand_addZeroLength(clientData, interp, argc, argv);

  } else if (strcmp(argv[1], "zeroLengthSection") == 0) {
    return TclCommand_addZeroLengthSection(clientData, interp, argc, argv);

  } else if (strcmp(argv[1], "zeroLengthRocking") == 0) {
    int result = TclCommand_addZeroLengthRocking(clientData, interp, argc, argv);
    return result;
  } else if (strcmp(argv[1], "zeroLengthContact2D") == 0) {
    int result = TclCommand_addZeroLengthContact2D(clientData, interp, argc, argv);
    return result;
  } else if (strcmp(argv[1], "zeroLengthContact3D") == 0) {
    return TclCommand_addZeroLengthContact3D(clientData, interp, argc, argv);

  } else if (strcmp(argv[1], "zeroLengthND") == 0) {
    return TclCommand_addZeroLengthND(clientData, interp, argc, argv);

  //
  // Joints
  //
  } else if ((strcmp(argv[1], "Joint2D") == 0) ||
             (strcmp(argv[1], "Joint2d") == 0)) {
    return
        TclBasicBuilder_addJoint2D(clientData, interp, argc, argv);


  } else if ((strcmp(argv[1], "Joint3D") == 0) ||
             (strcmp(argv[1], "Joint3d") == 0)) {
    return TclBasicBuilder_addJoint3D(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "genericClient") == 0) {
    return TclBasicBuilder_addGenericClient(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "genericCopy") == 0) {
    return TclBasicBuilder_addGenericCopy(clientData, interp, argc, argv);

  } else if ((strcmp(argv[1], "inelastic2dYS01") == 0) ||
             (strcmp(argv[1], "inelastic2dYS02") == 0) ||
             (strcmp(argv[1], "inelastic2dYS03") == 0) ||
             (strcmp(argv[1], "inelastic2dYS04") == 0) ||
             (strcmp(argv[1], "inelastic2dYS05") == 0)) {
    return TclBasicBuilder_addElement2dYS(clientData, interp, argc, argv);

  } else if ((strcmp(argv[1], "element2dGNL") == 0) ||
             (strcmp(argv[1], "elastic2dGNL") == 0)) {
    return TclBasicBuilder_addElastic2dGNL(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "beamColumnJoint") == 0) {
    return TclBasicBuilder_addBeamColumnJoint(clientData, interp, argc, argv);
  }

  // Kikuchi
  else if ((strcmp(argv[1], "multipleShearSpring") == 0) ||
           (strcmp(argv[1], "MSS") == 0)) {
    int result = TclBasicBuilder_addMultipleShearSpring(
        clientData, interp, argc, argv, theTclDomain, theTclBuilder);
    return result;
  }

  else if ((strcmp(argv[1], "multipleNormalSpring") == 0) ||
           (strcmp(argv[1], "MNS") == 0)) {
    int result = TclBasicBuilder_addMultipleNormalSpring(
        clientData, interp, argc, argv, theTclDomain, theTclBuilder);
    return result;
  }

  else if (strcmp(argv[1], "KikuchiBearing") == 0) {
    return TclBasicBuilder_addKikuchiBearing(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "YamamotoBiaxialHDR") == 0) {
    int result = TclBasicBuilder_addYamamotoBiaxialHDR(
        clientData, interp, argc, argv, theTclDomain, theTclBuilder);
    return result;
  }

  // MSN
  else if (strcmp(argv[1], "gradientInelasticBeamColumn") == 0) {
    int result = TclBasicBuilder_addGradientInelasticBeamColumn(clientData, interp, argc, argv);
    return result;
  }

  else {

    //
    // maybe element already loaded as c++ class from a package
    //

    // try existing loaded packages

    ElementPackageCommand *eleCommands = theElementPackageCommands;
    bool found = false;
    int result = TCL_ERROR;
    while (eleCommands != NULL && found == false) {
      if (strcmp(argv[1], eleCommands->funcName) == 0) {

        OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theTclDomain);
        void *theRes = (*(eleCommands->funcPtr))();
        if (theRes != 0) {
          Element *theEle = (Element *)theRes;
          result = theTclDomain->addElement(theEle);

          if (result >= 0)
            return TCL_OK;
          else
            return TCL_ERROR;
        }
        return TCL_ERROR;
        ;
      } else
        eleCommands = eleCommands->next;
    }
#if 0
    //
    // maybe element in a routine, check existing ones or try loading new ones
    //

    char *eleType = new char[strlen(argv[1]) + 1];
    strcpy(eleType, argv[1]);
    eleObj *eleObject = OPS_GetElementType(eleType, (int)strlen(eleType));

    delete[] eleType;

    if (eleObject != 0) {

      int result = Tcl_addWrapperElement(eleObject, clientData, interp, argc, argv,
                                         theTclDomain, theTclBuilder);

      if (result != 0)
        delete eleObject;
      else
        return result;
    }
#endif
    //
    // try loading new dynamic library containing a C++ class
    //

    void *libHandle;
    void *(*funcPtr)();
    int eleNameLength = (int)strlen(argv[1]);
    char *tclFuncName = new char[eleNameLength + 5];
    strcpy(tclFuncName, "OPS_");

    strcpy(&tclFuncName[4], argv[1]);

    opserr << "checking library: " << tclFuncName << endln;
    int res =
        getLibraryFunction(argv[1], tclFuncName, &libHandle, (void **)&funcPtr);

    delete[] tclFuncName;

    if (res == 0) {

      char *eleName = new char[eleNameLength + 1];
      strcpy(eleName, argv[1]);
      ElementPackageCommand *theEleCommand = new ElementPackageCommand;
      theEleCommand->funcPtr = funcPtr;
      theEleCommand->funcName = eleName;
      theEleCommand->next = theElementPackageCommands;
      theElementPackageCommands = theEleCommand;

      // OPS_ResetInput(clientData, interp, 2, argc, argv, theTclDomain,
      //                theTclBuilder);

      OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theTclDomain);
      void *theRes = (*funcPtr)();

      if (theRes != 0) {
        Element *theEle = (Element *)theRes;
        result = theTclDomain->addElement(theEle);
        if (result >= 0)
          return TCL_OK;
        else
          return TCL_ERROR;
      } else {
        return TCL_ERROR;
      }
    }
  }

  // If we get here, the element type is unknown
  opserr << "ERROR -- element of type " << argv[1] << " not known" << endln;
  return TCL_ERROR;
}

int
TclBasicBuilder_addMultipleShearSpring(ClientData clientData, Tcl_Interp *interp,
                                       int argc, TCL_Char ** const argv,
                                       Domain *theTclDomain, 
                                       [[maybe_unused]] TclBasicBuilder* unused)
{
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  if (builder == 0 || clientData == 0) {
    opserr << "WARNING builder has been destroyed - multipleShearSpring\n";
    return TCL_ERROR;
  }

  // 3-dim, 6-dof
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING multipleShearSpring command only works when ndm is 3 "
              "and ndf is 6"
           << endln;
    return TCL_ERROR;
  }

  // arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;
  int nSpring;
  int matTag;

  // material
  UniaxialMaterial *material = nullptr;
  UniaxialMaterial **theMaterials = nullptr;
  int recvMat = 0;

  // arguments (optional)
  double limDisp = 0.0;
  Vector oriX(0);
  Vector oriYp(3);
  oriYp(0) = 0.0;
  oriYp(1) = 1.0;
  oriYp(2) = 0.0;
  double mass = 0.0;

  //
  Element *theElement = nullptr;

  // error flag
  bool ifNoError = true;

  if (argc < 8) { // element multipleShearSpring eleTag? iNode? jNode? nSpring?
                  // -mat matTag?

    opserr << "WARNING insufficient arguments\n";
    ifNoError = false;

  } else {

    // argv[2~5]
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
      opserr << "WARNING invalid multipleShearSpring eleTag\n";
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
      opserr << "WARNING invalid iNode\n";
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode\n";
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[5], &nSpring) != TCL_OK || nSpring <= 0) {
      opserr << "WARNING invalid nSpring\n";
      ifNoError = false;
    }

    // argv[6~]
    for (int i = 6; i <= (argc - 1); ++i) {

      double value;

      if (strcmp(argv[i], "-mat") == 0 &&
          (i + 1) <= (argc - 1)) { // -mat matTag?

        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid matTag\n";
          ifNoError = false;
        }

        material = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (material == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "multipleShearSpring element: " << eleTag << endln;
          return TCL_ERROR;
        }

        recvMat++;
        i += 1;

      } else if (strcmp(argv[i], "-nMat") == 0 &&
                 (i + nSpring) <= (argc - 1)) { // -mat matTag?

        theMaterials = new UniaxialMaterial *[nSpring];
        for (int j = 0; j < nSpring; j++) {
          if (Tcl_GetInt(interp, argv[j + i + 1], &matTag) != TCL_OK) {
            opserr << "WARNING invalid matTag\n";
            ifNoError = false;
          }

          theMaterials[j] = builder->getTypedObject<UniaxialMaterial>(matTag);
          if (theMaterials[j] == 0) {
            opserr << "WARNING material model not found\n";
            opserr << "uniaxialMaterial: " << matTag << endln;
            opserr << "multipleShearSpring element: " << eleTag << endln;
            return TCL_ERROR;
          }
        }
        recvMat++;
        i += nSpring;

      } else if (strcmp(argv[i], "-orient") == 0 && (i + 6) <= (argc - 1) &&
                 Tcl_GetDouble(interp, argv[i + 4], &value) == TCL_OK) { 
        // <-orient x1? x2? x3? yp1? yp2? yp3?>

        oriX.resize(3);

        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            ifNoError = false;
          } else {
            oriX(j - 1) = value;
          }
        }

        i += 3;

        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            ifNoError = false;
          } else {
            oriYp(j - 1) = value;
          }
        }

        i += 3;

      } else if (strcmp(argv[i], "-orient") == 0 && (i + 3) <= (argc - 1)) { 
        // <-orient yp1? yp2? yp3?> の読み込み

        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            ifNoError = false;
          } else {
            oriYp(j - 1) = value;
          }
        }

        i += 3;

      } else if (strcmp(argv[i], "-mass") == 0 && (i + 1) <= (argc - 1)) { 
        // <-mass m?>

        if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK || mass <= 0) {
          opserr << "WARNING invalid mass\n";
          ifNoError = false;
        }

        i += 1;

      } else if (strcmp(argv[i], "-lim") == 0 && (i + 1) <= (argc - 1)) {
        // <-lim limDisp?>

        if (Tcl_GetDouble(interp, argv[i + 1], &limDisp) != TCL_OK || limDisp < 0) {
          opserr << "WARNING invalid limDisp\n";
          ifNoError = false;
        }

        i += 1;

      } else { // invalid option

        opserr << "WARNING invalid optional arguments \n";
        ifNoError = false;
        break;
      }
    }

  } // end input

  // confirm material
  if (recvMat != 1) {
    opserr << "WARNING wrong number of -mat inputs\n";
    opserr << "got " << recvMat << " inputs, but want 1 input\n";
    ifNoError = false;
  }

  // if error detected
  if (!ifNoError) {
    opserr << "Want: element multipleShearSpring eleTag? iNode? jNode? "
              "nSpring? -mat matTag? <-lim dsp> <-orient <x1? x2? x3?> yp1? "
              "yp2? yp3?> <-mass m?>\n";
    return TCL_ERROR;
  }

  // now create the multipleShearSpring
  if (theMaterials == 0) {
    theElement = new MultipleShearSpring(eleTag, iNode, jNode, nSpring,
                                         material, limDisp, oriYp, oriX, mass);
  } else {
    theElement = new MultipleShearSpring(eleTag, iNode, jNode, theMaterials,
                                         nSpring, limDisp, oriYp, oriX, mass);
    delete[] theMaterials;
  }

  // then add the multipleShearSpring to the domain
  if (theTclDomain->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "multipleShearSpring element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }

  // if get here we have successfully created the multipleShearSpring and added
  // it to the domain
  return TCL_OK;
}

static bool
errDetected(bool ifNoError, const char *msg)
{

 if (ifNoError) {
    opserr << "" << endln;
    opserr << "========================================" << endln;
    opserr << " element : input error detected" << endln;
    opserr << "------------------------------" << endln;
  }
  opserr << "  " << msg << endln;
  return false;
};

int
TclBasicBuilder_addMultipleNormalSpring(ClientData clientData, Tcl_Interp *interp,
                                        int argc, TCL_Char ** const argv,
                                        Domain *theTclDomain, TclBasicBuilder *theTclBuilder)
{

  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // 3-dim, 6-dof
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING multipleNormalSpring command only works when ndm is 3 "
              "and ndf is 6"
           << endln;
    return TCL_ERROR;
  }

  // arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;
  int nDivide;

  // arguments (necessary, input with -???)
  int matTag;
  UniaxialMaterial *material = nullptr;
  int shape = 0;
  double size;

  // arguments (optional, input with -???)
  double lambda = -1.0;
  Vector oriX(0);
  Vector oriYp(3);
  oriYp(0) = 0.0;
  oriYp(1) = 1.0;
  oriYp(2) = 0.0;
  double mass = 0.0;

  // input comfirmation
  int recvMat = 0;
  int recvShape = 0;
  int recvSize = 0;
  int recvLambda = 0;
  int recvOrient = 0;
  int recvMass = 0;

  //
  Element *theElement = nullptr;

  // error flag
  bool ifNoError = true;

  if (argc < 6) { // element multipleNormalSpring eleTag? iNode? jNode? nDivide?

    ifNoError = errDetected(ifNoError, "insufficient arguments");

  } else {

    // argv[2~5]
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
      ifNoError = errDetected(ifNoError, "invalid eleTag");
    }

    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
      ifNoError = errDetected(ifNoError, "invalid iNode");
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
      ifNoError = errDetected(ifNoError, "invalid jNode");
    }

    if (Tcl_GetInt(interp, argv[5], &nDivide) != TCL_OK || nDivide <= 0) {
      ifNoError = errDetected(ifNoError, "invalid nDivide");
    }

    // argv[6~]
    for (int i = 6; i <= (argc - 1); ++i) {

      double value;

      if (strcmp(argv[i], "-mat") == 0 &&
          (i + 1) <= (argc - 1)) { // -mat matTag?

        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          ifNoError = errDetected(ifNoError, "invalid matTag");
        }

        material = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (material == nullptr) {
          ifNoError = errDetected(ifNoError, "material model not found");
        }

        recvMat++;
        i += 1;

      } else if (strcmp(argv[i], "-shape") == 0 &&
                 (i + 1) <= (argc - 1)) { // -shape shape?

        if (strcmp(argv[i + 1], "round") == 0) {
          shape = 1; // round shape
        } else if (strcmp(argv[i + 1], "square") == 0) {
          shape = 2; // square
        } else {
          ifNoError = errDetected(
              ifNoError,
              "invalid shape (\"round\" or \"square\" are available)");
          goto error;
        }

        recvShape++;
        i += 1;

      } else if (strcmp(argv[i], "-size") == 0 &&
                 (i + 1) <= (argc - 1)) { // -size size?

        if (Tcl_GetDouble(interp, argv[i + 1], &size) != TCL_OK || size <= 0) {
          ifNoError = errDetected(ifNoError, "invalid size");
        }

        recvSize++;
        i += 1;

      } else if (strcmp(argv[i], "-lambda") == 0 &&
                 (i + 1) <= (argc - 1)) {
        // <-lambda lambda?>

        if (Tcl_GetDouble(interp, argv[i + 1], &lambda) != TCL_OK || lambda < 0) {
          ifNoError = errDetected(ifNoError, "invalid lambda");
        }

        recvLambda++;
        i += 1;

      } else if (strcmp(argv[i], "-orient") == 0 && (i + 6) <= (argc - 1) &&
                 Tcl_GetDouble(interp, argv[i + 4], &value) ==
                     TCL_OK) {
        // <-orient x1? x2? x3? yp1? yp2? yp3?>

        oriX.resize(3);

        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            ifNoError = errDetected(ifNoError, "invalid orient");
          } else {
            oriX(j - 1) = value;
          }
        }

        i += 3;

        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            ifNoError = errDetected(ifNoError, "invalid orient");
          } else {
            oriYp(j - 1) = value;
          }
        }

        recvOrient++;
        i += 3;

      } else if (strcmp(argv[i], "-orient") == 0 &&
                 (i + 3) <= (argc - 1)) {
        // <-orient yp1? yp2? yp3?>

        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            ifNoError = errDetected(ifNoError, "invalid orient");
          } else {
            oriYp(j - 1) = value;
          }
        }

        recvOrient++;
        i += 3;

      } else if (strcmp(argv[i], "-mass") == 0 &&
                 (i + 1) <= (argc - 1)) {
        // <-mass m?> の読み込み

        if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK || mass <= 0) {
          ifNoError = errDetected(ifNoError, "invalid mass");
        }

        recvMass++;
        i += 1;

      } else { // invalid option
        ifNoError = errDetected(ifNoError, "invalid optional arguments");
        break;
      }
    }

  } // end input

  // input cofirmation
  // necessary arguments
  if (recvMat != 1) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -mat inputs (got %d inputs, but want 1 input)",
            recvMat);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvShape != 1) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -shape inputs (got %d inputs, but want 1 input)",
            recvShape);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvSize != 1) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -size inputs (got %d inputs, but want 1 input)",
            recvSize);
    ifNoError = errDetected(ifNoError, buf);
  }

  // optional arguments
  if (recvLambda >= 2) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -lambda inputs (got %d inputs, but want 1 input)",
            recvLambda);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvOrient >= 2) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -ori inputs (got %d inputs, but want 1 input)",
            recvOrient);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvMass >= 2) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -mass inputs (got %d inputs, but want 1 input)",
            recvMass);
    ifNoError = errDetected(ifNoError, buf);
  }

  // if error detected
  if (!ifNoError) {
error:
    opserr << "Want: element multipleNormalSpring eleTag? iNode? jNode? "
              "\n    nDivide? -mat matTag? -shape shape? -size size? <-lambda "
              "\n    lambda?> <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>\n";
    opserr << "" << endln;
    return TCL_ERROR;
  }

  // now create the multipleNormalSpring
  theElement = new MultipleNormalSpring(eleTag, iNode, jNode, nDivide,
                    material, shape, size, lambda, oriYp, oriX, mass);

  // then add the multipleNormalSpring to the domain
  if (theTclDomain->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "multipleNormalSpring element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }

  // if get here we have successfully created the multipleNormalSpring and added
  // it to the domain
  return TCL_OK;
}

int
TclBasicBuilder_addKikuchiBearing(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // 3-dim, 6dof
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING KikuchiBearing command only works when ndm is 3 and ndf "
              "is 6"
           << endln;
    return TCL_ERROR;
  }

  // arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;

  // arguments (necessary, input with -???)
  int shape = 0;
  double size;
  double totalRubber;
  int nMSS;
  int matMSSTag;
  UniaxialMaterial *matMSS = nullptr;
  int nMNS;
  int matMNSTag;
  UniaxialMaterial *matMNS = nullptr;

  // arguments (optional, input with -???)
  double totalHeight = -1.0; // default: Norm(I->J)
  double limDisp = -1.0;     // default: INF
  double lambda = -1.0;      // default: INF
  Vector oriX(0);            // default: local-x Vec(I->J)
  Vector oriYp(3);
  oriYp(0) = 0.0;
  oriYp(1) = 1.0;
  oriYp(2) = 0.0; // default: global-Y
  double mass = 0.0;
  bool ifPDInput = true;
  bool ifTilt = true;
  double adjCi = 0.5;
  double adjCj = 0.5;
  bool ifBalance = false;
  double limFo = -1.0; // default: INF
  double limFi = -1.0; // default: INF
  int nIter = 1;

  // input comfirmation
  int recvShape = 0;
  int recvSize = 0;
  int recvHeight = 0;
  int recvNMSS = 0;
  int recvMatMSS = 0;
  int recvLimDisp = 0;
  int recvNMNS = 0;
  int recvMatMNS = 0;
  int recvLambda = 0;
  int recvOrient = 0;
  int recvMass = 0;
  int recvIfPD = 0;
  int recvIfTl = 0;
  int recvAdj = 0;
  int recvBal = 0;

  //
  Element *theElement = nullptr;

  // error flag
  bool ifNoError = true;

  if (argc < 5) { // element KikuchiBearing eleTag? iNode? jNode?
    ifNoError = errDetected(ifNoError, "insufficient arguments");

  } else {

    // argv[2~4]
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
      ifNoError = errDetected(ifNoError, "invalid eleTag");
    }

    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
      ifNoError = errDetected(ifNoError, "invalid iNode");
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
      ifNoError = errDetected(ifNoError, "invalid jNode");
    }

    // argv[5~]
    for (int i = 5; i <= (argc - 1); ++i) {

      double value;

      if (strcmp(argv[i], "-shape") == 0 &&
          (i + 1) <= (argc - 1)) { // -shape shape?

        if (strcmp(argv[i + 1], "round") == 0) {
          shape = 1; // round
        } else if (strcmp(argv[i + 1], "square") == 0) {
          shape = 2; // square
        } else {
          ifNoError = errDetected(
              ifNoError,
              "invalid shape (\"round\" or \"square\" are available)");
        }

        recvShape++;
        i += 1;

      } else if (strcmp(argv[i], "-size") == 0 && (i + 2) <= (argc - 1)) { 
        // -size size? totalRubber?

        if (Tcl_GetDouble(interp, argv[i + 1], &size) != TCL_OK || size <= 0) {
          ifNoError = errDetected(ifNoError, "invalid size");
        }

        if (Tcl_GetDouble(interp, argv[i + 2], &totalRubber) != TCL_OK ||
            totalRubber <= 0) {
          ifNoError = errDetected(ifNoError, "invalid totalRubber");
        }

        recvSize++;
        i += 2;

      } else if (strcmp(argv[i], "-totalHeight") == 0 && (i + 1) <= (argc - 1)) {
        // -totalHeight totalHeight?

        if (Tcl_GetDouble(interp, argv[i + 1], &totalHeight) != TCL_OK ||
            totalHeight <= 0) {
          ifNoError = errDetected(ifNoError, "invalid totalHeight");
        }

        recvHeight++;
        i += 1;

      } else if (strcmp(argv[i], "-nMSS") == 0 && (i + 1) <= (argc - 1)) {
        // -nMSS nMSS?

        if (Tcl_GetInt(interp, argv[i + 1], &nMSS) != TCL_OK || nMSS <= 0) {
          ifNoError = errDetected(ifNoError, "invalid nMSS");
        }

        recvNMSS++;
        i += 1;

      } else if (strcmp(argv[i], "-matMSS") == 0 && (i + 1) <= (argc - 1)) {
        // -matMSS matMSSTag?

        if (Tcl_GetInt(interp, argv[i + 1], &matMSSTag) != TCL_OK) {
          ifNoError = errDetected(ifNoError, "invalid matMSSTag");
        }

        matMSS = builder->getTypedObject<UniaxialMaterial>(matMSSTag);
        if (matMSS == 0) {
          ifNoError =
              errDetected(ifNoError, "material for MSS model not found");
        }

        recvMatMSS++;
        i += 1;

      } else if (strcmp(argv[i], "-limDisp") == 0 &&
                 (i + 1) <= (argc - 1)) {
        // <-limDisp limDisp?>

        if (Tcl_GetDouble(interp, argv[i + 1], &limDisp) != TCL_OK || limDisp < 0) {
          ifNoError = errDetected(ifNoError, "invalid limDisp");
        }

        recvLimDisp++;
        i += 1;

      } else if (strcmp(argv[i], "-nMNS") == 0 &&
                 (i + 1) <= (argc - 1)) { // -nMNS nMNS?

        if (Tcl_GetInt(interp, argv[i + 1], &nMNS) != TCL_OK || nMNS <= 0) {
          ifNoError = errDetected(ifNoError, "invalid nMNS");
        }

        recvNMNS++;
        i += 1;

      } else if (strcmp(argv[i], "-matMNS") == 0 &&
                 (i + 1) <= (argc - 1)) { // -matMNS matMNSTag?

        if (Tcl_GetInt(interp, argv[i + 1], &matMNSTag) != TCL_OK) {
          ifNoError = errDetected(ifNoError, "invalid matMNSTag");
        }

        matMNS = builder->getTypedObject<UniaxialMaterial>(matMNSTag);
        if (matMNS == 0) {
          ifNoError =
              errDetected(ifNoError, "material for MNS model not found");
        }

        recvMatMNS++;
        i += 1;

      } else if (strcmp(argv[i], "-lambda") == 0 &&
                 (i + 1) <= (argc - 1)) {
        // <-lambda lambda?>

        if (Tcl_GetDouble(interp, argv[i + 1], &lambda) != TCL_OK || lambda < 0) {
          ifNoError = errDetected(ifNoError, "invalid lambda");
        }

        recvLambda++;
        i += 1;

      } else if (strcmp(argv[i], "-orient") == 0 && (i + 6) <= (argc - 1) &&
                 Tcl_GetDouble(interp, argv[i + 4], &value) ==
                     TCL_OK) {
        // <-orient x1? x2? x3? yp1? yp2? yp3?>

        oriX.resize(3);

        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            ifNoError = errDetected(ifNoError, "invalid orient");
          } else {
            oriX(j - 1) = value;
          }
        }

        i += 3;

        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            ifNoError = errDetected(ifNoError, "invalid orient");
          } else {
            oriYp(j - 1) = value;
          }
        }

        recvOrient++;
        i += 3;

      } else if (strcmp(argv[i], "-orient") == 0 &&
                 (i + 3) <= (argc - 1)) {
        // <-orient yp1? yp2? yp3?>

        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            ifNoError = errDetected(ifNoError, "invalid orient");
          } else {
            oriYp(j - 1) = value;
          }
        }

        recvOrient++;
        i += 3;

      } else if (strcmp(argv[i], "-mass") == 0 &&
                 (i + 1) <= (argc - 1)) {
        // <-mass mass?>

        if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK || mass <= 0) {
          ifNoError = errDetected(ifNoError, "invalid mass");
        }

        recvMass++;
        i += 1;

      } else if (strcmp(argv[i], "-noPDInput") == 0) {
        // <-noPDInput>

        ifPDInput = false;

        recvIfPD++;
        i += 0;

      } else if (strcmp(argv[i], "-noTilt") == 0) {
        // <-noTilt>

        ifTilt = false;

        recvIfTl++;
        i += 0;

      } else if (strcmp(argv[i], "-adjustPDOutput") == 0 &&
                 (i + 2) <= (argc - 1)) { // -adjustPDOutput ci? cj?

        if (Tcl_GetDouble(interp, argv[i + 1], &adjCi) != TCL_OK) {
          ifNoError = errDetected(ifNoError, "invalid ci");
        }

        if (Tcl_GetDouble(interp, argv[i + 2], &adjCj) != TCL_OK) {
          ifNoError = errDetected(ifNoError, "invalid cj");
        }

        recvAdj++;
        i += 2;

      } else if (strcmp(argv[i], "-doBalance") == 0 &&
                 (i + 3) <= (argc - 1)) { // -doBalance limFo? limFi? nIter?

        if (Tcl_GetDouble(interp, argv[i + 1], &limFo) != TCL_OK || limFo <= 0) {
          ifNoError = errDetected(ifNoError, "invalid limFo");
        }

        if (Tcl_GetDouble(interp, argv[i + 2], &limFi) != TCL_OK || limFi <= 0) {
          ifNoError = errDetected(ifNoError, "invalid limFi");
        }

        if (Tcl_GetInt(interp, argv[i + 3], &nIter) != TCL_OK || nIter <= 0) {
          ifNoError = errDetected(ifNoError, "invalid nIter");
        }

        ifBalance = true;

        recvBal++;
        i += 3;

      } else { // invalid option

        ifNoError = errDetected(ifNoError, "invalid optional arguments");
        break;
      }
    }

  } // end input

  // input cofirmation
  // necessary arguments
  if (recvShape != 1) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -shape inputs (got %d inputs, but want 1 input)",
            recvShape);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvSize != 1) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -size inputs (got %d inputs, but want 1 input)",
            recvSize);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvNMSS != 1) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -NMSS inputs (got %d inputs, but want 1 input)",
            recvNMSS);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvMatMSS != 1) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -matMSS inputs (got %d inputs, but want 1 input)",
            recvMatMSS);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvNMNS != 1) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -NMNS inputs (got %d inputs, but want 1 input)",
            recvNMNS);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvMatMNS != 1) {
    char buf[100];
    snprintf(buf, 100,
            "wrong number of -matMNS inputs (got %d inputs, but want 1 input)",
            recvMatMNS);
    ifNoError = errDetected(ifNoError, buf);
  }

  // optional arguments
  if (recvHeight >= 2) {
    char buf[100];
    snprintf(buf,100,
        "wrong number of -totalHeight inputs (got %d inputs, but want 1 input)",
        recvHeight);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvLimDisp >= 2) {
    char buf[100];
    snprintf(buf,100,
            "wrong number of -limDisp inputs (got %d inputs, but want 1 input)",
            recvLimDisp);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvLambda >= 2) {
    char buf[100];
    snprintf(buf,100,
            "wrong number of -lambda inputs (got %d inputs, but want 1 input)",
            recvLambda);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvOrient >= 2) {
    char buf[100];
    snprintf(buf,100,
            "wrong number of -ori inputs (got %d inputs, but want 1 input)",
            recvOrient);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvMass >= 2) {
    char buf[100];
    snprintf(buf,100,
            "wrong number of -mass inputs (got %d inputs, but want 1 input)",
            recvMass);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvIfPD >= 2) {
    char buf[100];
    snprintf(buf, 100,
        "wrong number of -noPDInput inputs (got %d inputs, but want 1 input)",
        recvIfPD);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvIfTl >= 2) {
    char buf[100];
    snprintf(buf,100,
            "wrong number of -noTilt inputs (got %d inputs, but want 1 input)",
            recvIfTl);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvAdj >= 2) {
    char buf[100];
    snprintf(buf,100,
            "wrong number of -adjustPDOutput inputs (got %d inputs, but want 1 "
            "input)",
            recvAdj);
    ifNoError = errDetected(ifNoError, buf);
  }

  if (recvBal >= 2) {
    char buf[100];
    snprintf(buf, 100,
        "wrong number of -doBalance inputs (got %d inputs, but want 1 input)",
        recvBal);
    ifNoError = errDetected(ifNoError, buf);
  }

  // if error detected
  if (!ifNoError) {
    opserr << "Want: element KikuchiBearing eleTag? iNode? jNode?\n";
    opserr << "                             -shape shape? -size size? "
              "totalRubber? <-totalHeight totalHeight?>\n";
    opserr << "                             -nMSS nMSS? -matMSS matMSSTag? "
              "<-lim limDisp?>\n";
    opserr << "                             -nMNS nMNS? -matMNS matMNSTag? "
              "<-lambda lambda?>\n";
    opserr << "                             <-orient <x1? x2? x3?> yp1? yp2? "
              "yp3?> <-mass m?>\n";
    opserr << "                             <-noPDInput> <-noTilt> "
              "<-adjustPDOutput ci? cj?> <-doBalance limFo? limFi? nIter?>\n";
    opserr << "" << endln;
    return TCL_ERROR;
  }

  // now create the KikuchiBearing
  theElement = new KikuchiBearing(
      eleTag, iNode, jNode, shape, size, totalRubber, totalHeight, nMSS, matMSS,
      limDisp, nMNS, matMNS, lambda, oriYp, oriX, mass, ifPDInput, ifTilt,
      adjCi, adjCj, ifBalance, limFo, limFi, nIter);

  // then add the KikuchiBearing to the domain
  if (builder->getDomain()->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "KikuchiBearing element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclBasicBuilder_addYamamotoBiaxialHDR(ClientData clientData, Tcl_Interp *interp,
                                      int argc, TCL_Char ** const argv,
                                      [[maybe_unused]] Domain *theTclDomain_, 
                                      [[maybe_unused]] TclBasicBuilder *unused)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  Domain *theTclDomain = builder->getDomain();
  

  // 3-dim, 6-dof
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING YamamotoBiaxialHDR command only works when ndm is 3 and "
              "ndf is 6" << endln;
    return TCL_ERROR;
  }

  // arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;

  int Tp = 1;
  double DDo;
  double DDi;
  double Hr;

  // arguments (optional)
  double Cr = 1.0;
  double Cs = 1.0;
  Vector oriX(0);
  Vector oriYp(3);
  oriYp(0) = 0.0;
  oriYp(1) = 1.0;
  oriYp(2) = 0.0;
  double mass = 0.0;

  Element *theElement = nullptr;

  // error flag
  bool ifNoError = true;

  if (argc < 9) { 
    // element YamamotoBiaxialHDR eleTag? iNode? jNode? Tp? DDo? DDi? Hr?
    opserr << "WARNING insufficient arguments\n";
    ifNoError = false;

  } else {
    // argv[2~8]
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
      opserr << "WARNING invalid YamamotoBiaxialHDR eleTag\n";
      ifNoError = false;
    }

    // iNode
    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
      opserr << "WARNING invalid iNode\n";
      ifNoError = false;
    }

    // jNode
    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode\n";
      ifNoError = false;
    }

    // Tp
    if (strcmp(argv[5], "1") == 0) {
      Tp = 1; // Bridgestone X0.6R (EESD version)
    } else {
      opserr << "WARNING invalid YamamotoBiaxialHDR Tp" << endln;
      ifNoError = false;
    }

    // DDo
    if (Tcl_GetDouble(interp, argv[6], &DDo) != TCL_OK || DDo <= 0.0) {
      opserr << "WARNING invalid YamamotoBiaxialHDR DDo" << endln;
      ifNoError = false;
    }

    // DDi
    if (Tcl_GetDouble(interp, argv[7], &DDi) != TCL_OK || DDi < 0.0) {
      opserr << "WARNING invalid YamamotoBiaxialHDR DDi" << endln;
      ifNoError = false;
    }

    // Hr
    if (Tcl_GetDouble(interp, argv[8], &Hr) != TCL_OK || Hr <= 0.0) {
      opserr << "WARNING invalid YamamotoBiaxialHDR Hr" << endln;
      ifNoError = false;
    }

    // argv[9~]
    for (int i = 9; i <= (argc - 1); ++i) {
      double value;

      if (strcmp(argv[i], "-orient") == 0 && (i + 6) <= (argc - 1) &&
          Tcl_GetDouble(interp, argv[i + 4], &value) == TCL_OK) {
        // <-orient x1? x2? x3? yp1? yp2? yp3?>

        oriX.resize(3);

        // x1, x2, x3
        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            ifNoError = false;
          } else {
            oriX(j - 1) = value;
          }
        }

        i += 3;

        // yp1, yp2, yp3
        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            ifNoError = false;
          } else {
            oriYp(j - 1) = value;
          }
        }

        i += 3;

      } else if (strcmp(argv[i], "-orient") == 0 && (i + 3) <= (argc - 1)) {
        // <-orient yp1? yp2? yp3?>

        for (int j = 1; j <= 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            ifNoError = false;
          } else {
            oriYp(j - 1) = value;
          }
        }

        i += 3;

      } else if (strcmp(argv[i], "-mass") == 0 && (i + 1) <= (argc - 1)) {
        // <-mass m?>

        if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK || mass <= 0) {
          opserr << "WARNING invalid mass\n";
          ifNoError = false;
        }

        i += 1;

      } else if (strcmp(argv[i], "-coRS") == 0 && (i + 2) <= (argc - 1)) {
        // <-coRS cr? cs?>

        if (Tcl_GetDouble(interp, argv[i + 1], &Cr) != TCL_OK || Cr <= 0) {
          opserr << "WARNING invalid cr\n";
          ifNoError = false;
        }
        if (Tcl_GetDouble(interp, argv[i + 2], &Cs) != TCL_OK || Cs <= 0) {
          opserr << "WARNING invalid cs\n";
          ifNoError = false;
        }

        i += 2;

      } else {

        opserr << "WARNING invalid optional arguments \n";
        ifNoError = false;
        break;
      }
    }

  } // end input

  if (!ifNoError) {
    // want:
    opserr << "Want: element YamamotoBiaxialHDR eleTag? iNode? jNode? Tp? DDo? "
              "DDi? Hr?  <-coRS cr? cs?> <-orient <x1? x2? x3?> y1? y2? y3?> "
              "<-mass m?>\n";
    return TCL_ERROR;
  }

  // now create the YamamotoBiaxialHDR
  theElement = new YamamotoBiaxialHDR(eleTag, iNode, jNode, Tp, DDo, DDi, Hr,
                                      Cr, Cs, oriYp, oriX, mass);

  // then add the YamamotoBiaxialHDR to the domain
  if (theTclDomain->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "YamamotoBiaxialHDR element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }

  // if get here we have successfully created the YamamotoBiaxialHDR and added
  // it to the domain
  return TCL_OK;
}

int
TclBasicBuilder_addWheelRail(ClientData clientData, Tcl_Interp *interp, int argc,
                             TCL_Char ** const argv)
{
  constexpr static int eleArgStart = 1;
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  Element *theElement = nullptr;

  int pTag, pnLoad;
  //-------------Beginning of a 2D wheel-rail element(By Quan Gu, Yongdou Liu,
  // et al.) on 2018/10/29
  if (ndm == 2) {

    // check plane frame problem has 3 dof per node
    if (ndf != 3) {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for plane problem need 3 - elasticBeamColumn \n";
      return TCL_ERROR;
    }

    // check the number of arguments
    if ((argc - eleArgStart) < 8) {
      opserr << "WARNING bad command - want: elasticBeamColumn beamId iNode "
                "jNode A E I <alpha> <d> transTag <-mass m> <-cMass>\n";
      return TCL_ERROR;
    }

    // get the id, end nodes, and section properties
    int pNd1, transTag;

    double pDeltT, pVel, pInitLocation, pRWheel, pI, pE, pA;

    if (Tcl_GetInt(interp, argv[1 + eleArgStart], &pTag) != TCL_OK) {
      opserr << "WARNING invalid pTag: " << argv[1 + eleArgStart];
      opserr << " - WheelRail pTag iNode jNode";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[2 + eleArgStart], &pDeltT) != TCL_OK) {
      opserr << "WARNING invalid pDeltT - WheelRail " << pTag
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3 + eleArgStart], &pVel) != TCL_OK) {
      opserr << "WARNING invalid pVel - WheelRail " << pTag
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4 + eleArgStart], &pInitLocation) != TCL_OK) {
      opserr << "WARNING invalid pInitLocation - WheelRail " << pTag
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[5 + eleArgStart], &pNd1) != TCL_OK) {
      opserr << "WARNING invalid pNd1 - WheelRail " << pTag
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6 + eleArgStart], &pRWheel) != TCL_OK) {
      opserr << "WARNING invalid pRWheel - WheelRail " << pTag
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7 + eleArgStart], &pI) != TCL_OK) {
      opserr << "WARNING invalid pI - WheelRail " << pTag
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8 + eleArgStart], &pE) != TCL_OK) {
      opserr << "WARNING invalid pE - WheelRail " << pTag
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9 + eleArgStart], &pA) != TCL_OK) {
      opserr << "WARNING invalid pA - WheelRail " << pTag
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[10 + eleArgStart], &transTag) != TCL_OK) {
      opserr << "WARNING invalid transTag - WheelRail " << pTag
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    CrdTransf *theTransRWheel = builder->getTypedObject<CrdTransf>(transTag);

    if (Tcl_GetInt(interp, argv[11 + eleArgStart], &pnLoad) != TCL_OK) {
      opserr << "WARNING invalid I - WheelRail " << pTag
             << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    //----------------------------------
    Vector *pNodeList = 0;
    Vector *pDeltaYList = 0;
    Vector *pDeltaYLocationList = 0;

    if (strcmp(argv[12 + eleArgStart], "-NodeList") == 0) {
      int pathSize;
      TCL_Char **pathStrings;

      // int debug =
      //     Tcl_SplitList(interp, argv[13 + eleArgStart], &pathSize, &pathStrings);

      if (Tcl_SplitList(interp, argv[13 + eleArgStart], &pathSize, &pathStrings) !=
          TCL_OK) {
        opserr << "WARNING problem splitting path list "
               << argv[13 + eleArgStart] << " - ";
        opserr << " NodeList -values {path} ... \n";
        return TCL_OK;
      }
      pNodeList = new Vector(pathSize);
      for (int i = 0; i < pathSize; ++i) {
        double value;
        // int debug = Tcl_GetDouble(interp, pathStrings[i], &value);
        if (Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
          opserr << "WARNING problem reading path data value " << pathStrings[i]
                 << " - ";
          opserr << " -strain {path} ... \n";
          return 0;
        }
        (*pNodeList)(i) = value;
      } // for
    }
    if (strcmp(argv[14 + eleArgStart], "-DeltaYList") == 0) {
      int pathSize;
      TCL_Char **pathStrings;
      if (Tcl_SplitList(interp, argv[15 + eleArgStart], &pathSize, &pathStrings) !=
          TCL_OK) {
        opserr << "WARNING problem splitting path list "
               << argv[15 + eleArgStart] << " - ";
        opserr << " NodeList -values {path} ... \n";
        return TCL_OK;
      }
      pDeltaYList = new Vector(pathSize);
      for (int i = 0; i < pathSize; ++i) {
        double value;
        if (Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
          opserr << "WARNING problem reading path data value " << pathStrings[i]
                 << " - ";
          opserr << " -strain {path} ... \n";
          return 0;
        }
        (*pDeltaYList)(i) = value;
      } // for
    }
    if (strcmp(argv[16 + eleArgStart], "-LocationList") == 0) {
      int pathSize;
      TCL_Char **pathStrings;
      if (Tcl_SplitList(interp, argv[17 + eleArgStart], &pathSize, &pathStrings) !=
          TCL_OK) {
        opserr << "WARNING problem splitting path list "
               << argv[17 + eleArgStart] << " - ";
        opserr << " NodeList -values {path} ... \n";
        return TCL_OK;
      }
      pDeltaYLocationList = new Vector(pathSize);
      for (int i = 0; i < pathSize; ++i) {
        double value;
        if (Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
          opserr << "WARNING problem reading path data value " << pathStrings[i]
                 << " - ";
          opserr << " -strain {path} ... \n";
          return 0;
        }
        (*pDeltaYLocationList)(i) = value;
      }
    }
    theElement = new WheelRail(pTag, pDeltT, pVel, pInitLocation, pNd1, pRWheel,
                               pI, pE, pA, theTransRWheel, pnLoad, pNodeList,
                               pDeltaYList, pDeltaYLocationList);

  } 
  // -- End of a 2D wheel-rail element(By Quan Gu, Yongdou Liu, et al.) on 2018/10/29

  else if (ndm == 3) {
    opserr << G3_ERROR_PROMPT << "Unimplemented." << endln;
    return TCL_ERROR;
  }

  // add the WheelRail element to the Domain
  if (builder->getDomain()->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "YamamotoBiaxialHDR element: " << pTag << endln;
    delete theElement;
    return TCL_ERROR;
  }

  return 0;
}

