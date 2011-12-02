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
// $Date: 2005-06-14 18:50:31 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewCommand/MyTrussPackage.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the myTruss_Init procedure that
// is invoked when the 'package' is added to the interpreter.


#include <tcl.h>
#include "MyTruss.h"
#include <Domain.h>
#include <TclModelBuilder.h>

extern Domain theDomain;
extern TclModelBuilder *theBuilder;

static int OpenSees_MyTruss(ClientData clientData, 
			    Tcl_Interp *interp, 
			    int argc, 
			    TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc != 6) {
      interp->result = "WARNING bad command - want: myTruss eleId iNode jNode Area matID";
      return TCL_ERROR;
  }    

  // get the id, x_loc and y_loc
  int trussId, iNode, jNode, matID;
  double A, M = 0.0;
  if (Tcl_GetInt(interp, argv[1], &trussId) != TCL_OK) {
     interp->result = "WARNING invalid eleId - myTruss eleId iNode jNode Area matID";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &iNode) != TCL_OK) {
     interp->result = "WARNING invalid iNode- myTruss eleId iNode jNode Area matID";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &jNode) != TCL_OK) {
     interp->result = "WARNING invalid jNode- myTruss eleId iNode jNode Area matID";
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &A) != TCL_OK) {
     interp->result = "WARNING invalid A- myTruss eleId iNode jNode Area matID";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5], &matID) != TCL_OK) {
     interp->result = "WARNING invalid matId- myTruss eleId iNode jNode Area matID";
      return TCL_ERROR;
  }
  if (argc == 8 && Tcl_GetDouble(interp, argv[7], &M) != TCL_OK) {
     interp->result = "WARNING invalid matId- myTruss eleId iNode jNode Area matID";
     return TCL_ERROR;
  }  
  
  UniaxialMaterial *theMaterial = theBuilder->getUniaxialMaterial(matID);

  if (theMaterial == 0) {
    opserr << "WARNING TclPlaneTruss - truss - no Material found with tag ";
    opserr << matID << endln;
    return TCL_ERROR;
  }

  // now create the truss and add it to the Domain
  MyTruss *theTruss = new MyTruss(trussId,iNode,jNode,*theMaterial,A);
  if (theTruss == 0) {
    opserr << "WARNING TclPlaneTruss - addMyTruss - ran out of memory for node ";
    opserr << trussId << endln;
    return TCL_ERROR;
  }
  if (theDomain.addElement(theTruss) == false) {
    delete theTruss;
    opserr << "WARNING TclPlaneTruss - addTruss - could not add Truss to domain ";
    opserr << trussId << endln;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}


/*
 *----------------------------------------------------------------------
 *
 * mytruss_Init --
 *
 *	This is a package initialization procedure, which is called
 *	by Tcl when this package is to be added to an interpreter.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

extern "C" int
Mytruss_Init(Tcl_Interp *interp)
{
    int code;

    if (Tcl_InitStubs(interp, TCL_VERSION, 1) == NULL) {
	return TCL_ERROR;
    }
    
    // add the package to list of available packages
    code = Tcl_PkgProvide(interp, "myTruss", "1.0");
    if (code != TCL_OK) {
	return code;
    }

    // add the myTruss command to the interpreter
    Tcl_CreateCommand(interp, "myTruss", OpenSees_MyTruss, (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    // done
    return TCL_OK;
}
