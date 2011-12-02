//
// Description: This file contains the implementation of the
// 
//
// Jinchi Lu and Zhaohui Yang (May 2004)

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <Twenty_Node_Brick.h>

#include <TclModelBuilder.h>

//#ifdef _DEBUG
//#define TCL_Char const char
//#endif

extern void printCommand(int argc, TCL_Char **argv);

/*  *****************************************************************************

    20-NODE BRICK

    ***************************************************************************** */

int
TclModelBuilder_addTwentyNodeBrick(ClientData clientData, Tcl_Interp *interp,
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
		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with 20NodeBrick element\n";
		return TCL_ERROR;
	}

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc-argStart) < 22) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element 20NodeBrick eleTag? N1? N2? N3? N4? N5? N6? N7? N8? N9? N10? N11? N12? N13? N14? N15? N16? N17? N18? N19? N20? matTag? <b1? b2? b3?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int brickId, Nod[20], matID;
	double b1 = 0.0;
	double b2 = 0.0;
	double b3 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &brickId) != TCL_OK) {
    opserr << "WARNING invalid 20NodeBrick eleTag" << endln;
    return TCL_ERROR;
  }

  for (int i=0; i<20; i++)
    if (Tcl_GetInt(interp, argv[1+argStart+i], &(Nod[i])) != TCL_OK) {
      opserr << "WARNING invalid Node number\n";
      opserr << "20NodeBrick element: " << brickId << endln;
      return TCL_ERROR;
	}

  if (Tcl_GetInt(interp, argv[21+argStart], &matID) != TCL_OK) {
     opserr << "WARNING invalid matID\n";
     opserr << "20NodeBrick element: " << brickId << endln;
     return TCL_ERROR;
  }

	if ((argc-argStart) >= 23) {
		if (Tcl_GetDouble(interp, argv[22+argStart], &b1) != TCL_OK) {
			opserr << "WARNING invalid b1\n";
			opserr << "20NodeBrick element: " << brickId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 24) {
		if (Tcl_GetDouble(interp, argv[23+argStart], &b2) != TCL_OK) {
			opserr << "WARNING invalid b2\n";
			opserr << "20NodeBrick element: " << brickId << endln;
			return TCL_ERROR;
		}
	}
	if ((argc-argStart) >= 25) {
		if (Tcl_GetDouble(interp, argv[24+argStart], &b3) != TCL_OK) {
			opserr << "WARNING invalid b3\n";
			opserr << "20NodeBrick element: " << brickId << endln;
			return TCL_ERROR;
		}
	}

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);

  if (theMaterial == 0) {
      opserr << "WARNING material not found\n";
      opserr << "Material: " << matID;
      opserr << "\20NodeBrick element: " << brickId << endln;
      return TCL_ERROR;
  }

  // now create the brick and add it to the Domain
  Twenty_Node_Brick *theTwentyNodeBrick =
      new Twenty_Node_Brick(brickId,Nod[0],Nod[1],Nod[2],Nod[3],Nod[4],Nod[5],Nod[6],Nod[7],
				Nod[8],Nod[9],Nod[10],Nod[11],Nod[12],Nod[13],Nod[14],
				Nod[15],Nod[16],Nod[17],Nod[18],Nod[19],
		       *theMaterial, b1, b2, b3);
  if (theTwentyNodeBrick == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "20NodeBrick element: " << brickId << endln;
      return TCL_ERROR;
  }
  if (theTclDomain->addElement(theTwentyNodeBrick) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "20NodeBrick element: " << brickId << endln;
      delete theTwentyNodeBrick;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}

