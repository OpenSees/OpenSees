//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#include <stdlib.h>
#include <string.h>

#include <tcl.h>
#include <Matrix.h>
#include <NDMaterial.h>
#include <BasicModelBuilder.h>
#include <G3_Logging.h>

typedef const char TCL_Char;

Tcl_CmdProc PlaneStress_usePlaneStressMaterial;
static Tcl_CmdProc PlaneStress_setStrainPlaneStressMaterial;
static Tcl_CmdProc PlaneStress_getStressPlaneStressMaterial;
static Tcl_CmdProc PlaneStress_getTangPlaneStressMaterial;


static int count;
static int countsTillCommit;

// constructor: the constructor will add certain commands to the interpreter
int TclCommand_usePlaneStress(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);


  // get the matID form command line
  int matID;
  if (Tcl_GetInt(interp, argv[1], &matID) != TCL_OK) {
    opserr << "WARNING could not read matID: plane stressTest matID?\n";
    return TCL_ERROR;
  }

  // get the material from the modelbuilder with matID
  // and set the testing material to point to a copy of it
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr) {
    opserr << "WARNING no material found with matID\n";
    return TCL_ERROR;
  } else {
    theMaterial = theMaterial->getCopy("PlaneStress");
  }


  Tcl_CreateCommand(interp, "setStrain",
                    PlaneStress_setStrainPlaneStressMaterial,
                    (ClientData)theMaterial, NULL);

  Tcl_CreateCommand(interp, "getStress",
                    PlaneStress_getStressPlaneStressMaterial,
                    (ClientData)theMaterial, NULL);

  Tcl_CreateCommand(interp, "getTangent",
                    PlaneStress_getTangPlaneStressMaterial,
                    (ClientData)theMaterial, NULL);







  Tcl_DeleteCommand(interp, "setMaterial");
  Tcl_DeleteCommand(interp, "setStrain");
  Tcl_DeleteCommand(interp, "getStress");
  Tcl_DeleteCommand(interp, "getTangent");

  return TCL_OK;
}

//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static int
PlaneStress_setStrainPlaneStressMaterial(ClientData clientData,
                                                          Tcl_Interp *interp,
                                                          int argc,
                                                          TCL_Char ** const argv)
{
  NDMaterial* theMaterial = (NDMaterial*)clientData;

  // check number of arguments in command line
  if (argc < 4) {
    opserr << "WARNING bad command - want: strainPlaneStressTest strain?\n";
    return TCL_ERROR;
  }

  // get the matID form command line
  static double strain[3];
  static Vector strainV(strain, 3);
  if (Tcl_GetDouble(interp, argv[1], &strain[0]) != TCL_OK) {
    opserr << "WARNING could not read strain: strainPlaneStressTest strain?\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &strain[1]) != TCL_OK) {
    opserr << "WARNING could not read strain: strainPlaneStressTest strain?\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &strain[2]) != TCL_OK) {
    opserr << "WARNING could not read strain: strainPlaneStressTest strain?\n";
    return TCL_ERROR;
  }

  theMaterial->setTrialStrain(strainV);
  if (count == countsTillCommit) {
    theMaterial->commitState();
    count = 1;
  } else
    count++;
  return TCL_OK;
}

static int
PlaneStress_getStressPlaneStressMaterial(ClientData clientData,
                                                          Tcl_Interp *interp,
                                                          int argc,
                                                          TCL_Char ** const argv)
{
  NDMaterial* theMaterial = (NDMaterial*)clientData;
  static Vector stress(3);

  // delete the old testing material
    stress = theMaterial->getStress();
    char buffer[60];
    sprintf(buffer, "%.10e %.10e %.10e", stress(0), stress(1), stress(2));
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
}

static int
PlaneStress_getTangPlaneStressMaterial(ClientData clientData,
                                                        Tcl_Interp *interp,
                                                        int argc,
                                                        TCL_Char ** const argv)
{
  static Matrix tangent(3, 3);
  NDMaterial* theMaterial = (NDMaterial*)clientData;

  tangent = theMaterial->getTangent();
  char buffer[180];
  sprintf(buffer, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e",
          tangent(0, 0), tangent(0, 1), tangent(0, 2), tangent(1, 0),
          tangent(1, 1), tangent(1, 2), tangent(2, 0), tangent(2, 1),
          tangent(2, 2));
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);
  return TCL_OK;
}
