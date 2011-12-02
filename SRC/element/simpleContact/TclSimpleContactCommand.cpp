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
                                                                        
// $Revision: 1.1 $
// $Date: 2004/03/25 23:27:22 $
// $Source: /home/cee/pmackenz/Software/OpenSees/SRC/element/SimpleContact/TclSimpleContactCommand.cpp,v $
                                                                        
// File: ~/element/TclSimpleContactCommand.C
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the implementation of the TclModelBuilder_addSimpleContact2D(),
// TclModelBuilder_addSimpleContact3D(), and TclModelBuilder_addBeamContact3D() commands. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include "SimpleContact2D.h"
#include "SimpleContact3D.h"
#include "BeamContact3D.h"

#include <CrdTransf.h>

#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

/*  *****************************************************************************
    
    S I M P L E   C O N T A C T   2 D

    ***************************************************************************** */

int
TclModelBuilder_addSimpleContact2D(ClientData clientData, Tcl_Interp *interp,  
				int argc, 
				TCL_Char **argv, 
				Domain*theTclDomain,
				TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }


  // check the number of arguments is correct

  if (argc-eleArgStart < 9) {
    opserr << "WARNING insufficient arguments here\n";
    printCommand(argc, argv);
    opserr << "Want: element SimpleContact2D eleTag? iNode? jNode? slaveNode? lambdaNode? matTag? tolGap? tolForce?\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int SimpleContact2DId, iNode, jNode, slaveNode, lambdaNode, matID;
  double tolG, tolF;

  if (Tcl_GetInt(interp, argv[1+eleArgStart], &SimpleContact2DId) != TCL_OK) {
    opserr << "WARNING invalid SimpleContact2D eleTag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "SimpleContact2D element: " << SimpleContact2DId << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "SimpleContact2D element: " << SimpleContact2DId << endln;
     return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[4+eleArgStart], &slaveNode) != TCL_OK) {
     opserr << "WARNING invalid slaveNode\n";
     opserr << "SimpleContact2D element: " << SimpleContact2DId << endln;
     return TCL_ERROR;
  }  
  
  if (Tcl_GetInt(interp, argv[5+eleArgStart], &lambdaNode) != TCL_OK) {
     opserr << "WARNING invalid lambdaNode\n";
     opserr << "SimpleContact2D element: " << SimpleContact2DId << endln;
     return TCL_ERROR;
  }  
 
  if (Tcl_GetInt(interp, argv[6+eleArgStart], &matID) != TCL_OK) {
     opserr << "WARNING invalid matID\n";
     opserr << "SimpleContact2D element: " << SimpleContact2DId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7+eleArgStart], &tolG) != TCL_OK) {
     opserr << "WARNING invalid tolGap\n";
     opserr << "SimpleContact2D element: " << SimpleContact2DId << endln;
     return TCL_ERROR;
  }
  
    if (Tcl_GetDouble(interp, argv[8+eleArgStart], &tolF) != TCL_OK) {
     opserr << "WARNING invalid tolForce\n";
     opserr << "SimpleContact2D element: " << SimpleContact2DId << endln;
     return TCL_ERROR;
  }
  
  
  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
      
  if (theMaterial == 0) {
      opserr << "WARNING material not found\n";
      opserr << "Material: " << matID;
      opserr << "\nSimpleContact2D element: " << SimpleContact2DId << endln;
      return TCL_ERROR;
  }
  
  // now create the SimpleContact2D and add it to the Domain
  SimpleContact2D *theSimpleContact2D = 
      new SimpleContact2D(SimpleContact2DId,iNode,jNode,slaveNode,
				lambdaNode, *theMaterial, tolG, tolF);
  if (theSimpleContact2D == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "SimpleContact2D element: " << SimpleContact2DId << endln;
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theSimpleContact2D) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "SimpleContact2D element: " << SimpleContact2DId << endln;
      delete theSimpleContact2D;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}

/*  *****************************************************************************
    
    S I M P L E   C O N T A C T   3 D

    ***************************************************************************** */

int
TclModelBuilder_addSimpleContact3D(ClientData clientData, Tcl_Interp *interp,  
				int argc, 
				TCL_Char **argv, 
				Domain*theTclDomain,
				TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }


  // check the number of arguments is correct

  if (argc-eleArgStart < 11) {
    opserr << "WARNING insufficient arguments here\n";
    printCommand(argc, argv);
    opserr << "Want: element SimpleContact3D eleTag?  iNode? jNode? kNode? lNode? slaveNode? lambdaNode? matTag? tolGap? tolF\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int SimpleContact3DId, iNode, jNode, kNode, lNode, slaveNode, lambdaNode, matID;
  double tolG, tolF;

  if (Tcl_GetInt(interp, argv[1+eleArgStart], &SimpleContact3DId) != TCL_OK) {
    opserr << "WARNING invalid SimpleContact3D eleTag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "SimpleContact3D element: " << SimpleContact3DId << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "SimpleContact3D element: " << SimpleContact3DId << endln;
     return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[4+eleArgStart], &kNode) != TCL_OK) {
     opserr << "WARNING invalid kNode\n";
     opserr << "SimpleContact3D element: " << SimpleContact3DId << endln;
     return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[5+eleArgStart], &lNode) != TCL_OK) {
     opserr << "WARNING invalid lNode\n";
     opserr << "SimpleContact3D element: " << SimpleContact3DId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6+eleArgStart], &slaveNode) != TCL_OK) {
     opserr << "WARNING invalid slaveNode\n";
     opserr << "SimpleContact3D element: " << SimpleContact3DId << endln;
     return TCL_ERROR;
  }  
  
  if (Tcl_GetInt(interp, argv[7+eleArgStart], &lambdaNode) != TCL_OK) {
     opserr << "WARNING invalid lambdaNode\n";
     opserr << "SimpleContact3D element: " << SimpleContact3DId << endln;
     return TCL_ERROR;
  } 
  
  if (Tcl_GetInt(interp, argv[8+eleArgStart], &matID) != TCL_OK) {
     opserr << "WARNING invalid matID\n";
     opserr << "SimpleContact3D element: " << SimpleContact3DId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[9+eleArgStart], &tolG) != TCL_OK) {
     opserr << "WARNING invalid tolGap\n";
     opserr << "SimpleContact3D element: " << SimpleContact3DId << endln;
     return TCL_ERROR;
  }  
  
  if (Tcl_GetDouble(interp, argv[10+eleArgStart], &tolF) != TCL_OK) {
     opserr << "WARNING invalid tolForce\n";
     opserr << "SimpleContact3D element: " << SimpleContact3DId << endln;
     return TCL_ERROR;
  } 
  
  
  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
      
  if (theMaterial == 0) {
      opserr << "WARNING material not found\n";
      opserr << "Material: " << matID;
      opserr << "\nSimpleContact3D element: " << SimpleContact3DId << endln;
      return TCL_ERROR;
  }
  
  // now create the SimpleContact3D and add it to the Domain
  SimpleContact3D *theSimpleContact3D = 
      new SimpleContact3D(SimpleContact3DId,iNode,jNode, kNode, lNode,
				 slaveNode,lambdaNode,*theMaterial, tolG, tolF);
  if (theSimpleContact3D == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "SimpleContact3D element: " << SimpleContact3DId << endln;
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theSimpleContact3D) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "SimpleContact3D element: " << SimpleContact3DId << endln;
      delete theSimpleContact3D;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}

/*  *****************************************************************************
    
    B E A M   C O N T A C T   3 D

    ***************************************************************************** */

int
TclModelBuilder_addBeamContact3D(ClientData clientData, Tcl_Interp *interp,  
				int argc, 
				TCL_Char **argv, 
				Domain*theTclDomain,
				TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }


  // check the number of arguments is correct

  if (argc-eleArgStart < 11) {
    opserr << "WARNING insufficient arguments here\n";
    printCommand(argc, argv);
    opserr << "Want: element BeamContact3D eleTag?  iNode? jNode? slaveNode? lambdaNode? radius? crdTransf? matTag? tolGap? tolF\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int BeamContact3DId, iNode, jNode, slaveNode, lambdaNode, transTag, matID;
  double rad, tolG, tolF;

  if (Tcl_GetInt(interp, argv[1+eleArgStart], &BeamContact3DId) != TCL_OK) {
    opserr << "WARNING invalid BeamContact3D eleTag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "BeamContact3D element: " << BeamContact3DId << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "BeamContact3D element: " << BeamContact3DId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4+eleArgStart], &slaveNode) != TCL_OK) {
     opserr << "WARNING invalid slaveNode\n";
     opserr << "BeamContact3D element: " << BeamContact3DId << endln;
     return TCL_ERROR;
  }  
  
  if (Tcl_GetInt(interp, argv[5+eleArgStart], &lambdaNode) != TCL_OK) {
     opserr << "WARNING invalid lambdaNode\n";
     opserr << "BeamContact3D element: " << BeamContact3DId << endln;
     return TCL_ERROR;
  } 

  if (Tcl_GetDouble(interp, argv[6+eleArgStart], &rad) != TCL_OK) {
     opserr << "WARNING invalid radius\n";
     opserr << "BeamContact3D element: " << BeamContact3DId << endln;
     return TCL_ERROR;
  }  

  if (Tcl_GetInt(interp, argv[7+eleArgStart], &transTag) != TCL_OK) {
     opserr << "WARNING invalid transTag\n";
     opserr << "BeamContact3D element: " << BeamContact3DId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[8+eleArgStart], &matID) != TCL_OK) {
     opserr << "WARNING invalid matID\n";
     opserr << "BeamContact3D element: " << BeamContact3DId << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[9+eleArgStart], &tolG) != TCL_OK) {
     opserr << "WARNING invalid tolGap\n";
     opserr << "BeamContact3D element: " << BeamContact3DId << endln;
     return TCL_ERROR;
  }  
  
  if (Tcl_GetDouble(interp, argv[10+eleArgStart], &tolF) != TCL_OK) {
     opserr << "WARNING invalid tolForce\n";
     opserr << "BeamContact3D element: " << BeamContact3DId << endln;
     return TCL_ERROR;
  } 
  

  CrdTransf *theTrans = OPS_GetCrdTransf(transTag);

  if (theTrans ==0) {
    opserr << "WARNING transformation object not found - beamContact3D" << BeamContact3DId;
    return TCL_ERROR;
  }
  
  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
      
  if (theMaterial == 0) {
      opserr << "WARNING material not found\n";
      opserr << "Material: " << matID;
      opserr << "\nBeamContact3D element: " << BeamContact3DId << endln;
      return TCL_ERROR;
  }
  
  // now create the BeamContact3D and add it to the Domain
  BeamContact3D *theBeamContact3D = 
      new BeamContact3D(BeamContact3DId,iNode,jNode, 
				 slaveNode,lambdaNode, rad, *theTrans,
				 *theMaterial, tolG, tolF);
  if (theBeamContact3D == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "BeamContact3D element: " << BeamContact3DId << endln;
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theBeamContact3D) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "BeamContact3D element: " << BeamContact3DId << endln;
      delete theBeamContact3D;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}
