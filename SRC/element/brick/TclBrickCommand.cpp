/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.8 $
// $Date: 2009-08-07 00:18:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/brick/TclBrickCommand.cpp,v $
                                                                        
// Written: fmk 
// Created: 03/01
//
// What: "@(#) TclBrickCommand.cpp, revA"

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <Brick.h>
#include <BbarBrick.h>
#include <BbarBrickWithSensitivity.h>

#ifdef _FLBrick
#include <FLBrick.h>
#endif

#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addBrick(ClientData clientData, Tcl_Interp *interp,  int argc, 
			 TCL_Char **argv, Domain*theTclDomain,
			 TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 11) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element Brick eleTag? Node1? Node2? Node3? Node4? Node5? Node6? Node7? Node 8? matTag?\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int BrickId, Node1, Node2, Node3, Node4, matID;
  int          Node5, Node6, Node7, Node8 ;

  if (Tcl_GetInt(interp, argv[1+eleArgStart], &BrickId) != TCL_OK) {
    opserr << "WARNING invalid Brick eleTag" << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2+eleArgStart], &Node1) != TCL_OK) {
    opserr << "WARNING invalid Node1\n";
    opserr << "Brick element: " << BrickId << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3+eleArgStart], &Node2) != TCL_OK) {
     opserr << "WARNING invalid Node2\n";
     opserr << "Brick element: " << BrickId << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4+eleArgStart], &Node3) != TCL_OK) {
     opserr << "WARNING invalid Node3\n";
     opserr << "Brick element: " << BrickId << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5+eleArgStart], &Node4) != TCL_OK) {
     opserr << "WARNING invalid Node4\n";
     opserr << "Brick element: " << BrickId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6+eleArgStart], &Node5) != TCL_OK) {
    opserr << "WARNING invalid Node5\n";
    opserr << "Brick element: " << BrickId << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[7+eleArgStart], &Node6) != TCL_OK) {
     opserr << "WARNING invalid Node6\n";
     opserr << "Brick element: " << BrickId << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[8+eleArgStart], &Node7) != TCL_OK) {
     opserr << "WARNING invalid Node7\n";
     opserr << "Brick element: " << BrickId << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[9+eleArgStart], &Node8) != TCL_OK) {
     opserr << "WARNING invalid Node8\n";
     opserr << "Brick element: " << BrickId << endln;
     return TCL_ERROR;
  }


  if (Tcl_GetInt(interp, argv[10+eleArgStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matTag\n";
    opserr << "Brick element: " << BrickId << endln;
    return TCL_ERROR;
  }

   NDMaterial *theMaterial = OPS_getNDMaterial(matID);      
      
  if (theMaterial== 0) {
    opserr << "WARNING material not found\n";
    opserr << "material tag: " << matID;
    opserr << "\nBrick element: " << BrickId << endln;
    return TCL_ERROR;
  }
  
  double b1=0.0;
  double b2=0.0; 
  double b3=0.0;
  if ((argc-eleArgStart) > 11) {
    if (Tcl_GetDouble(interp, argv[11+eleArgStart], &b1) != TCL_OK) {
       opserr << "WARNING invalid b1\n";
       opserr << "Brick element: " << BrickId << endln;
       return TCL_ERROR;
	}
  }

  if ((argc-eleArgStart) > 12) {
    if (Tcl_GetDouble(interp, argv[12+eleArgStart], &b2) != TCL_OK) {
       opserr << "WARNING invalid b2\n";
       opserr << "Brick element: " << BrickId << endln;
       return TCL_ERROR;
	}
  }

  if ((argc-eleArgStart) > 13) {
    if (Tcl_GetDouble(interp, argv[13+eleArgStart], &b3) != TCL_OK) {
       opserr << "WARNING invalid b3\n";
       opserr << "Brick element: " << BrickId << endln;
       return TCL_ERROR;
	}
  }

  // now create the Brick and add it to the Domain
  Element *theBrick = 0;
  if (strcmp(argv[1],"stdBrick") == 0) {
    theBrick = new Brick(BrickId,Node1,Node2,Node3,Node4,
			 Node5, Node6, Node7, Node8, *theMaterial,
			 b1, b2, b3);
  }
  else if (strcmp(argv[1],"bbarBrickWithSensitivity") == 0) {
    theBrick = new BbarBrickWithSensitivity(BrickId,Node1,Node2,Node3,Node4,
			 Node5, Node6, Node7, Node8, *theMaterial,
			 b1, b2, b3);
  }
  else if (strcmp(argv[1],"bbarBrick") == 0) {
    theBrick = new BbarBrick(BrickId,Node1,Node2,Node3,Node4,
			     Node5, Node6, Node7, Node8, *theMaterial, b1, b2, b3);
	#ifdef _FLBrick
	  } else if (strcmp(argv[1],"flBrick") == 0) {
		theBrick = new FLBrick(BrickId,Node1,Node2,Node3,Node4,
				   Node5, Node6, Node7, Node8, *theMaterial, b1, b2, b3);
	  }
	#endif

  } else {
    opserr << "WARNING element " << argv[1] << " type not recognized\n";
    return TCL_ERROR;
  }

  

  if (theBrick == 0) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "Brick element: " << BrickId << endln;
    return TCL_ERROR;
  }

  if (theTclDomain->addElement(theBrick) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "Brick element: " << BrickId << endln;
    delete theBrick;
    return TCL_ERROR;
  }

  // if get here we have successfully created the node and added it to the domain
  return TCL_OK;
}


