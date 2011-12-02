#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <Domain.h>
#include <TclModelBuilder.h>
#include "MyTruss.h"

#ifdef _USRDLL
#include <windows.h>
#define DllExport _declspec(dllexport)
#else
#define DllExport
#endif

extern "C" DllExport int
TclCommand_myTruss(ClientData clientData, Tcl_Interp *interp,  int argc, 
		   TCL_Char **argv, Domain*theTclDomain,
		   TclModelBuilder *theTclBuilder)
{
  // check the number of arguments is correct
  if (argc < 6) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element truss eleTag? iNode? jNode? A? matTag?\n";
    return TCL_ERROR;
  }    

  int ndm = theTclBuilder->getNDM();

  // get the id and end nodes 
  int trussId, iNode, jNode, matID;
  double A;
  if (Tcl_GetInt(interp, argv[2], &trussId) != TCL_OK) {
    opserr << "WARNING invalid truss eleTag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "truss element: " << trussId << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "truss element: " << trussId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK) {
    opserr << "WARNING invalid A\n";
    opserr << "truss element: " << trussId << endln;
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[6], &matID) != TCL_OK) {
    opserr << "WARNING invalid matTag\n";
    opserr << "truss element: " << trussId << endln;
    return TCL_ERROR;
  }
  
  UniaxialMaterial *theMaterial = theTclBuilder->getUniaxialMaterial(matID);
      
  if (theMaterial == 0) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << matID;
    opserr << "\ntruss element: " << trussId << endln;
    return TCL_ERROR;
  }
  
  // now create the truss and add it to the Domain
  Element *theTruss = new MyTruss(trussId, iNode, jNode, *theMaterial, A);

  if (theTruss == 0) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "truss element: " << trussId << endln;
    return TCL_ERROR;
  }

  if (theTclDomain->addElement(theTruss) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "truss element: " << trussId << endln;
    delete theTruss;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}
