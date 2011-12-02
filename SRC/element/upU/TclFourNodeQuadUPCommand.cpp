//
// Description: This file contains the implementation of the TclModelBuilder_addFourNodeQuadUP()
// command. 
//

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <FourNodeQuadUP.h>

#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

/*  *****************************************************************************
    
    Q U A D  U_P

    ***************************************************************************** */

int
TclModelBuilder_addFourNodeQuadUP(ClientData clientData, Tcl_Interp *interp,  
				int argc, 
				TCL_Char **argv, 
				Domain*theTclDomain,
				TclModelBuilder *theTclBuilder)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 3) {
		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with QuadUP element\n";
		return TCL_ERROR;
	}

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc-argStart) < 12) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element FourNodeQuadUP eleTag? iNode? jNode? kNode? lNode? thk? type? matTag? bulk? rho? perm_x? perm_y? <b1? b2? pressure? dM? dK?>\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int FourNodeQuadUPId, iNode, jNode, kNode, lNode, matID;
  double thickness, bk, r, perm1, perm2;
	double p = 0.0;		// uniform normal traction (pressure)
	double b1 = 0.0;
	double b2 = 0.0;
	double dM = 0.0;
	double dK = 0.0;

  TCL_Char *type;
  if (Tcl_GetInt(interp, argv[argStart], &FourNodeQuadUPId) != TCL_OK) {
    opserr << "WARNING invalid FourNodeQuadUP eleTag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[1+argStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2+argStart], &jNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[3+argStart], &kNode) != TCL_OK) {
     opserr << "WARNING invalid kNode\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }  
  
  if (Tcl_GetInt(interp, argv[4+argStart], &lNode) != TCL_OK) {
     opserr << "WARNING invalid lNode\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }  

  if (Tcl_GetDouble(interp, argv[5+argStart], &thickness) != TCL_OK) {
     opserr << "WARNING invalid thickness\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }  
  
  type = argv[6+argStart];
  
  if (Tcl_GetInt(interp, argv[7+argStart], &matID) != TCL_OK) {
     opserr << "WARNING invalid matID\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

	if (Tcl_GetDouble(interp, argv[8+argStart], &bk) != TCL_OK) {
     opserr << "WARNING invalid fluid bulk modulus\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }  

  if (Tcl_GetDouble(interp, argv[9+argStart], &r) != TCL_OK) {
     opserr << "WARNING invalid fluid mass density\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }  

  if (Tcl_GetDouble(interp, argv[10+argStart], &perm1) != TCL_OK) {
     opserr << "WARNING invalid lateral permeability\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }  

  if (Tcl_GetDouble(interp, argv[11+argStart], &perm2) != TCL_OK) {
     opserr << "WARNING invalid vertical permeability\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }  

	if ((argc-argStart) >= 13) {
		if (Tcl_GetDouble(interp, argv[12+argStart], &b1) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 14) {
		if (Tcl_GetDouble(interp, argv[13+argStart], &b2) != TCL_OK) {
			opserr << "WARNING invalid b2\n";
			opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 15) {
		if (Tcl_GetDouble(interp, argv[14+argStart], &p) != TCL_OK) {
			opserr << "WARNING invalid pressure\n";
			opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 16) {
		if (Tcl_GetDouble(interp, argv[15+argStart], &dM) != TCL_OK) {
			opserr << "WARNING invalid dM\n";
			opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 17) {
		if (Tcl_GetDouble(interp, argv[16+argStart], &dK) != TCL_OK) {
			opserr << "WARNING invalid dK\n";
			opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
			return TCL_ERROR;
		}
	}

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
      
  if (theMaterial == 0) {
      opserr << "WARNING material not found\n";
      opserr << "Material: " << matID;
      opserr << "\nFourNodeQuadUP element: " << FourNodeQuadUPId << endln;
      return TCL_ERROR;
  }
  
  // now create the FourNodeQuadUP and add it to the Domain
  FourNodeQuadUP *theFourNodeQuadUP = 
      new FourNodeQuadUP(FourNodeQuadUPId,iNode,jNode,kNode,lNode,
		       *theMaterial, type, thickness, bk, r, perm1, perm2, b1, b2, p, dM, dK);
  if (theFourNodeQuadUP == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "FourNodeQuad element: " << FourNodeQuadUPId << endln;
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theFourNodeQuadUP) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "FourNodeQuad element: " << FourNodeQuadUPId << endln;
      delete theFourNodeQuadUP;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}


