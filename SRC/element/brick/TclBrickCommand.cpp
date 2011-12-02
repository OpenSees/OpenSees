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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-25 23:32:49 $
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

   NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);      
      
  if (theMaterial== 0) {
    opserr << "WARNING material not found\n";
    opserr << "material tag: " << matID;
    opserr << "\nBrick element: " << BrickId << endln;
    return TCL_ERROR;
  }
  
  // now create the Brick and add it to the Domain
  Brick *theBrick = new Brick(BrickId,Node1,Node2,Node3,Node4,
			      Node5, Node6, Node7, Node8, *theMaterial);

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

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}


int
TclModelBuilder_addBBarBrick(ClientData clientData, Tcl_Interp *interp,  int argc, 
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
    opserr << "Want: element BbarBrick eleTag? Node1? Node2? Node3? Node4? Node5? Node6? Node7? Node 8? matTag?\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int BbarBrickId, Node1, Node2, Node3, Node4, matID;
  int          Node5, Node6, Node7, Node8 ;

  if (Tcl_GetInt(interp, argv[1+eleArgStart], &BbarBrickId) != TCL_OK) {
    opserr << "WARNING invalid BbarBrick eleTag" << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2+eleArgStart], &Node1) != TCL_OK) {
    opserr << "WARNING invalid Node1\n";
    opserr << "BbarBrick element: " << BbarBrickId << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3+eleArgStart], &Node2) != TCL_OK) {
     opserr << "WARNING invalid Node2\n";
     opserr << "BbarBrick element: " << BbarBrickId << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4+eleArgStart], &Node3) != TCL_OK) {
     opserr << "WARNING invalid Node3\n";
     opserr << "BbarBrick element: " << BbarBrickId << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5+eleArgStart], &Node4) != TCL_OK) {
     opserr << "WARNING invalid Node4\n";
     opserr << "BbarBrick element: " << BbarBrickId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6+eleArgStart], &Node5) != TCL_OK) {
    opserr << "WARNING invalid Node5\n";
    opserr << "BbarBrick element: " << BbarBrickId << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[7+eleArgStart], &Node6) != TCL_OK) {
     opserr << "WARNING invalid Node6\n";
     opserr << "BbarBrick element: " << BbarBrickId << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[8+eleArgStart], &Node7) != TCL_OK) {
     opserr << "WARNING invalid Node7\n";
     opserr << "BbarBrick element: " << BbarBrickId << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[9+eleArgStart], &Node8) != TCL_OK) {
     opserr << "WARNING invalid Node8\n";
     opserr << "BbarBrick element: " << BbarBrickId << endln;
     return TCL_ERROR;
  }


  if (Tcl_GetInt(interp, argv[10+eleArgStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matTag\n";
    opserr << "BbarBrick element: " << BbarBrickId << endln;
    return TCL_ERROR;
  }

   NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);      
      
  if (theMaterial== 0) {
    opserr << "WARNING material not found\n";
    opserr << "material tag: " << matID;
    opserr << "\nBbarBrick element: " << BbarBrickId << endln;
    return TCL_ERROR;
  }
  
  // now create the BbarBrick and add it to the Domain
  BbarBrick *theBbarBrick = new BbarBrick(BbarBrickId,Node1,Node2,Node3,Node4,
			      Node5, Node6, Node7, Node8, *theMaterial);

  if (theBbarBrick == 0) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "BbarBrick element: " << BbarBrickId << endln;
    return TCL_ERROR;
  }

  if (theTclDomain->addElement(theBbarBrick) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "BbarBrick element: " << BbarBrickId << endln;
    delete theBbarBrick;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

