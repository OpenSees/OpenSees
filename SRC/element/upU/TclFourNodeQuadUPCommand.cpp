//
// Description: This file contains the implementation of
// TclModelBuilder_addFourNodeQuadUP() ,
// TclModelBuilder_addBrickUP() ,
// TclModelBuilder_addNineFourNodeQuadUP() , and
// TclModelBuilder_addTwentyEightNodeBrickUP()
//
// Zhaohui Yang and Jinchi Lu (May 2004)

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <FourNodeQuadUP.h>
#include <BrickUP.h>
#include <Nine_Four_Node_QuadUP.h>
#include <Twenty_Eight_Node_BrickUP.h>

#include <TclModelBuilder.h>

//#ifdef _DEBUG
//#define TCL_Char const char
//#endif

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

  if ((argc-argStart) < 11) {
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

  if (Tcl_GetInt(interp, argv[6+argStart], &matID) != TCL_OK) {
     opserr << "WARNING invalid matID\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

	if (Tcl_GetDouble(interp, argv[7+argStart], &bk) != TCL_OK) {
     opserr << "WARNING invalid fluid bulk modulus\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[8+argStart], &r) != TCL_OK) {
     opserr << "WARNING invalid fluid mass density\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[9+argStart], &perm1) != TCL_OK) {
     opserr << "WARNING invalid lateral permeability\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[10+argStart], &perm2) != TCL_OK) {
     opserr << "WARNING invalid vertical permeability\n";
     opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

	if ((argc-argStart) >= 12) {
		if (Tcl_GetDouble(interp, argv[11+argStart], &b1) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 13) {
		if (Tcl_GetDouble(interp, argv[12+argStart], &b2) != TCL_OK) {
			opserr << "WARNING invalid b2\n";
			opserr << "FourNodeQuadUP element: " << FourNodeQuadUPId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 14) {
		if (Tcl_GetDouble(interp, argv[13+argStart], &p) != TCL_OK) {
			opserr << "WARNING invalid pressure\n";
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
		       *theMaterial, "PlaneStrain", thickness, bk, r, perm1, perm2, b1, b2, p);
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


/*  *****************************************************************************

    BRICK  U_P

    ***************************************************************************** */

int
TclModelBuilder_addBrickUP(ClientData clientData, Tcl_Interp *interp,
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

	if (theTclBuilder->getNDM() != 3 || theTclBuilder->getNDF() != 4) {
		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with QuadUP element\n";
		return TCL_ERROR;
	}

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc-argStart) < 15) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element brickUP eleTag? N1? N2? N3? N4? N5? N6? N7? N8? matTag? bulk? rhof? perm_x? perm_y? perm_z? <b1? b2? b3?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int brickUPId, Nod[8], matID;
  double bk, r, perm1, perm2, perm3;
	double b1 = 0.0;
	double b2 = 0.0;
	double b3 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &brickUPId) != TCL_OK) {
    opserr << "WARNING invalid brickUP eleTag" << endln;
    return TCL_ERROR;
  }

  for (int i=0; i<8; i++)
    if (Tcl_GetInt(interp, argv[1+argStart+i], &(Nod[i])) != TCL_OK) {
      opserr << "WARNING invalid Node number\n";
      opserr << "brickUP element: " << brickUPId << endln;
      return TCL_ERROR;
	}

  if (Tcl_GetInt(interp, argv[9+argStart], &matID) != TCL_OK) {
     opserr << "WARNING invalid matID\n";
     opserr << "brickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

	if (Tcl_GetDouble(interp, argv[10+argStart], &bk) != TCL_OK) {
     opserr << "WARNING invalid fluid bulk modulus\n";
     opserr << "brickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[11+argStart], &r) != TCL_OK) {
     opserr << "WARNING invalid fluid mass density\n";
     opserr << "brickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[12+argStart], &perm1) != TCL_OK) {
     opserr << "WARNING invalid permeability_x\n";
     opserr << "brickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[13+argStart], &perm2) != TCL_OK) {
     opserr << "WARNING invalid permeability_y\n";
     opserr << "brickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[14+argStart], &perm3) != TCL_OK) {
     opserr << "WARNING invalid permeability_z\n";
     opserr << "brickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

	if ((argc-argStart) >= 16) {
		if (Tcl_GetDouble(interp, argv[15+argStart], &b1) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "brickUP element: " << brickUPId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 17) {
		if (Tcl_GetDouble(interp, argv[16+argStart], &b2) != TCL_OK) {
			opserr << "WARNING invalid b2\n";
			opserr << "brickUP element: " << brickUPId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 18) {
		if (Tcl_GetDouble(interp, argv[17+argStart], &b3) != TCL_OK) {
			opserr << "WARNING invalid b3\n";
			opserr << "brickUP element: " << brickUPId << endln;
			return TCL_ERROR;
		}
	}

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);

  if (theMaterial == 0) {
      opserr << "WARNING material not found\n";
      opserr << "Material: " << matID;
      opserr << "\nbrickUP element: " << brickUPId << endln;
      return TCL_ERROR;
  }

  // now create the brickUP and add it to the Domain
  BrickUP *theBrickUP =
      new BrickUP(brickUPId,Nod[0],Nod[1],Nod[2],Nod[3],Nod[4],Nod[5],Nod[6],Nod[7],
		       *theMaterial, bk, r, perm1, perm2, perm3, b1, b2, b3);
  if (theBrickUP == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "brickUP element: " << brickUPId << endln;
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theBrickUP) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "brickUP element: " << brickUPId << endln;
      delete theBrickUP;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}

/*  *****************************************************************************

    9-4-N O D E  Q U A D  U_P

    ***************************************************************************** */

int
TclModelBuilder_addNineFourNodeQuadUP(ClientData clientData, Tcl_Interp *interp,
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

	if (theTclBuilder->getNDM() != 2) {
		opserr << "WARNING -- model dimensions not compatible with 9-4-NodeQuadUP element\n";
		return TCL_ERROR;
	}

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc-argStart) < 16) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element FourNodeQuadUP eleTag? Node1? ... Node9? thk? type? matTag? bulk? rho? perm_x? perm_y? <b1? b2? pressure? dM? dK?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int NineFourNodeQuadUPId, Node[9], matID;
  double thickness, bk, r, perm1, perm2;
	double b1 = 0.0;
	double b2 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &NineFourNodeQuadUPId) != TCL_OK) {
    opserr << "WARNING invalid FourNodeQuadUP eleTag" << endln;
    return TCL_ERROR;
  }
  for (int i=1; i<=9; i++) {
     if (Tcl_GetInt(interp, argv[i+argStart], &Node[i-1]) != TCL_OK) {
       opserr << "WARNING invalid Node\n";
       opserr << "FourNodeQuadUP element: " << NineFourNodeQuadUPId << endln;
       return TCL_ERROR;
	 }
  }

  if (Tcl_GetDouble(interp, argv[10+argStart], &thickness) != TCL_OK) {
     opserr << "WARNING invalid thickness\n";
     opserr << "FourNodeQuadUP element: " << NineFourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[11+argStart], &matID) != TCL_OK) {
     opserr << "WARNING invalid matID\n";
     opserr << "FourNodeQuadUP element: " << NineFourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

	if (Tcl_GetDouble(interp, argv[12+argStart], &bk) != TCL_OK) {
     opserr << "WARNING invalid fluid bulk modulus\n";
     opserr << "FourNodeQuadUP element: " << NineFourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[13+argStart], &r) != TCL_OK) {
     opserr << "WARNING invalid fluid mass density\n";
     opserr << "FourNodeQuadUP element: " << NineFourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[14+argStart], &perm1) != TCL_OK) {
     opserr << "WARNING invalid lateral permeability\n";
     opserr << "FourNodeQuadUP element: " << NineFourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[15+argStart], &perm2) != TCL_OK) {
     opserr << "WARNING invalid vertical permeability\n";
     opserr << "FourNodeQuadUP element: " << NineFourNodeQuadUPId << endln;
     return TCL_ERROR;
  }

	if ((argc-argStart) >= 17) {
		if (Tcl_GetDouble(interp, argv[16+argStart], &b1) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "FourNodeQuadUP element: " << NineFourNodeQuadUPId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 18) {
		if (Tcl_GetDouble(interp, argv[17+argStart], &b2) != TCL_OK) {
			opserr << "WARNING invalid b2\n";
			opserr << "FourNodeQuadUP element: " << NineFourNodeQuadUPId << endln;
			return TCL_ERROR;
		}
	}

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);

  if (theMaterial == 0) {
      opserr << "WARNING material not found\n";
      opserr << "Material: " << matID;
      opserr << "\nFourNodeQuadUP element: " << NineFourNodeQuadUPId << endln;
      return TCL_ERROR;
  }

  // now create the FourNodeQuadUP and add it to the Domain
  NineFourNodeQuadUP *theNineFourNodeQuadUP =
      new NineFourNodeQuadUP(NineFourNodeQuadUPId,Node[0],Node[1],Node[2],
	                         Node[3],Node[4],Node[5],Node[6],Node[7],Node[8],
		                     *theMaterial, "PlaneStrain", thickness, bk, r, perm1, perm2, b1, b2);
  if (theNineFourNodeQuadUP == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "FourNodeQuad element: " << NineFourNodeQuadUPId << endln;
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theNineFourNodeQuadUP) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "FourNodeQuad element: " << NineFourNodeQuadUPId << endln;
      delete theNineFourNodeQuadUP;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}

/*  *****************************************************************************

    20-8 NODED BRICK  U_P

    ***************************************************************************** */

int
TclModelBuilder_addTwentyEightNodeBrickUP(ClientData clientData, Tcl_Interp *interp,
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

	if (theTclBuilder->getNDM() != 3 ) {
		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with 20_8_BrickUP element\n";
		return TCL_ERROR;
	}

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc-argStart) < 27) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element 20_8_BrickUP eleTag? N1? N2? N3? N4? N5? N6? N7? N8? N9? N10? N11? N12? N13? N14? N15? N16? N17? N18? N19? N20? matTag? bulk? rhof? perm_x? perm_y? perm_z? <b1? b2? b3?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int brickUPId, Nod[20], matID;
  double bk, r, perm1, perm2, perm3;
	double b1 = 0.0;
	double b2 = 0.0;
	double b3 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &brickUPId) != TCL_OK) {
    opserr << "WARNING invalid 20_8_BrickUP eleTag" << endln;
    return TCL_ERROR;
  }

  for (int i=0; i<20; i++)
    if (Tcl_GetInt(interp, argv[1+argStart+i], &(Nod[i])) != TCL_OK) {
      opserr << "WARNING invalid Node number\n";
      opserr << "20_8_BrickUP element: " << brickUPId << endln;
      return TCL_ERROR;
	}

  if (Tcl_GetInt(interp, argv[21+argStart], &matID) != TCL_OK) {
     opserr << "WARNING invalid matID\n";
     opserr << "20_8_BrickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

	if (Tcl_GetDouble(interp, argv[22+argStart], &bk) != TCL_OK) {
     opserr << "WARNING invalid fluid bulk modulus\n";
     opserr << "20_8_BrickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[23+argStart], &r) != TCL_OK) {
     opserr << "WARNING invalid fluid mass density\n";
     opserr << "20_8_BrickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[24+argStart], &perm1) != TCL_OK) {
     opserr << "WARNING invalid permeability_x\n";
     opserr << "20_8_BrickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[25+argStart], &perm2) != TCL_OK) {
     opserr << "WARNING invalid permeability_y\n";
     opserr << "20_8_BrickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[26+argStart], &perm3) != TCL_OK) {
     opserr << "WARNING invalid permeability_z\n";
     opserr << "20_8_BrickUP element: " << brickUPId << endln;
     return TCL_ERROR;
  }

	if ((argc-argStart) >= 28) {
		if (Tcl_GetDouble(interp, argv[27+argStart], &b1) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "20_8_BrickUP element: " << brickUPId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 29) {
		if (Tcl_GetDouble(interp, argv[28+argStart], &b2) != TCL_OK) {
			opserr << "WARNING invalid b2\n";
			opserr << "20_8_BrickUP element: " << brickUPId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 30) {
		if (Tcl_GetDouble(interp, argv[29+argStart], &b3) != TCL_OK) {
			opserr << "WARNING invalid b3\n";
			opserr << "20_8_BrickUP element: " << brickUPId << endln;
			return TCL_ERROR;
		}
	}

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);

  if (theMaterial == 0) {
      opserr << "WARNING material not found\n";
      opserr << "Material: " << matID;
      opserr << "\n20_8_BrickUP element: " << brickUPId << endln;
      return TCL_ERROR;
  }

  // now create the brickUP and add it to the Domain
  TwentyEightNodeBrickUP *theTwentyEightNodeBrickUP =
      new TwentyEightNodeBrickUP(brickUPId,Nod[0],Nod[1],Nod[2],Nod[3],Nod[4],Nod[5],Nod[6],Nod[7],
				Nod[8],Nod[9],Nod[10],Nod[11],Nod[12],Nod[13],Nod[14],
				Nod[15],Nod[16],Nod[17],Nod[18],Nod[19],
		       *theMaterial, bk, r, perm1, perm2, perm3, b1, b2, b3);
  if (theTwentyEightNodeBrickUP == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "20_8_BrickUP element: " << brickUPId << endln;
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theTwentyEightNodeBrickUP) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "20_8_BrickUP element: " << brickUPId << endln;
      delete theTwentyEightNodeBrickUP;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}

