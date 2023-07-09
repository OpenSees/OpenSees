#include <tcl.h>
#include <elementAPI.h>

#include <UniaxialMaterial.h>
#include <string.h>

#include <StandardStream.h>
StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

#include <SimulationInformation.h>
SimulationInformation simulationInfo;
SimulationInformation *theSimulationInfoPtr = 0;

#include <FE_Datastore.h>
FE_Datastore *theDatabase  =0;

#define TCL_Char const char

static UniaxialMaterial *theTestingUniaxialMaterial =0;
static double ops_strain = 0;

extern UniaxialMaterial *OPS_ParseUniaxialMaterialCommand(const char *matType);

extern int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp, int cArg, int mArg, TCL_Char **argv, Domain *domain);

const char * getInterpPWD(Tcl_Interp *interp) {
  static char *pwd = 0;

  if (pwd != 0)
    delete [] pwd;

#ifdef _TCL84
  Tcl_Obj *cwd = Tcl_FSGetCwd(interp);
  if (cwd != NULL) {
    int length;
    const char *objPWD = Tcl_GetStringFromObj(cwd, &length);
    pwd = new char[length+1];
    strcpy(pwd, objPWD);
    Tcl_DecrRefCount(cwd);	
  }
#else

  Tcl_DString buf;
  const char *objPWD = Tcl_GetCwd(interp, &buf);

  pwd = new char[strlen(objPWD)+1];
  strcpy(pwd, objPWD);

  Tcl_DStringFree(&buf);

#endif
  return pwd;
}


int  
OPS_UniaxialMaterialCommand(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv) {

  // check number of arguments
  if (argc < 3) {
    opserr << "WARNING insufficient number of uniaxial material arguments\n";
    opserr << "Want: uniaxialMaterial type? tag? <specific material args>" << endln;
    return TCL_ERROR;
  }

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, 0);

  const char *matType = argv[1];
  UniaxialMaterial *theMaterial = OPS_ParseUniaxialMaterialCommand(matType);


}

int OPS_testUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {

  int tag = 0;
  if (argc < 2) {
    Tcl_SetResult(interp, "WARNING bad command - want: testUniaxialMaterial tag?", TCL_STATIC);
    return TCL_ERROR;
  }    

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    Tcl_SetResult(interp, "WARNING could not read strain: setStrain strain?", TCL_STATIC);
    return TCL_ERROR;
  }

  if (theTestingUniaxialMaterial != 0) 
    delete theTestingUniaxialMaterial;

  theTestingUniaxialMaterial=OPS_getUniaxialMaterial(tag);

  if (theTestingUniaxialMaterial == 0) {
    Tcl_SetResult(interp, "WARNING no active UniaxialMaterial found", TCL_STATIC);    
    return TCL_ERROR;
  }

  return TCL_OK;
}

int OPS_setStrain(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {

  if (argc < 2) {
    Tcl_SetResult(interp, "WARNING bad command - want: strainUniaxialTest strain?", TCL_STATIC);
    return TCL_ERROR;
  }    

  if (Tcl_GetDouble(interp, argv[1], &ops_strain) != TCL_OK) {
    Tcl_SetResult(interp, "WARNING could not read strain: strainUniaxialTest strain?", TCL_STATIC);
    return TCL_ERROR;
  }

  if (theTestingUniaxialMaterial !=0 ) {
    theTestingUniaxialMaterial->setTrialStrain(ops_strain);
  } else {
    Tcl_SetResult(interp, "WARNING no active UniaxialMaterial - use testUniaxialMaterial command", TCL_STATIC);    
    return TCL_ERROR;
  }
  return TCL_OK;
}


int OPS_getStrain(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {

  double strain = 0.0;

  if (theTestingUniaxialMaterial !=0) {
    strain = theTestingUniaxialMaterial->getStrain();
    sprintf(interp->result,"%.10e",strain);
    return TCL_OK;
  } else {
    Tcl_SetResult(interp, "WARNING no active UniaxialMaterial - use testUniaxialMaterial command", TCL_STATIC);    
    return TCL_ERROR;
  }
}

int OPS_getStress(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {

  double stress = 0.0;

  if (theTestingUniaxialMaterial !=0) {
    stress = theTestingUniaxialMaterial->getStress();
    sprintf(interp->result,"%.10e",stress);
    return TCL_OK;
  } else {
    Tcl_SetResult(interp, "WARNING no active UniaxialMaterial - use testUniaxialMaterial command", TCL_STATIC);    
    return TCL_ERROR;
  }
}

int  OPS_getTangent(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
  double tangent = 0.0;

  // delete the old testing material
  if (theTestingUniaxialMaterial !=0) {
    tangent = theTestingUniaxialMaterial->getTangent();
    sprintf(interp->result,"%.10e",tangent);
    return TCL_OK;
  } else {
    Tcl_SetResult(interp, "WARNING no active UniaxialMaterial - use testUniaaxialMaterial command", TCL_STATIC);    
    return TCL_ERROR;
  }
}


int addOpenSeesCommands(Tcl_Interp *interp) {

  Tcl_CreateCommand(interp, "uniaxialMaterial", &OPS_UniaxialMaterialCommand,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 

  Tcl_CreateCommand(interp, "testUniaxialMaterial", &OPS_testUniaxialMaterial,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 

  Tcl_CreateCommand(interp, "setStrain", &OPS_setStrain,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
  
  Tcl_CreateCommand(interp, "getStrain", &OPS_getStrain,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
  
  Tcl_CreateCommand(interp, "getStress", &OPS_getStress,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
  
  Tcl_CreateCommand(interp, "getTangent", &OPS_getTangent,
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
  return 0;
}


