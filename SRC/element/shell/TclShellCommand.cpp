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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-25 23:33:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/TclShellCommand.cpp,v $
                                                                        
// Written: fmk 
// Created: 03/01
//
// What: "@(#) TclShellCommand.C, revA"

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <ShellMITC4.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addShellMITC4(ClientData clientData, Tcl_Interp *interp,  int argc, 
			 TCL_Char **argv, Domain*theTclDomain,
			 TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 7) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element ShellMITC4 eleTag? iNode? jNode? kNode? lNode? secTag?\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int ShellMITC4Id, iNode, jNode, kNode, lNode, matID;
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &ShellMITC4Id) != TCL_OK) {
    opserr << "WARNING invalid ShellMITC4 eleTag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "ShellMITC4 element: " << ShellMITC4Id << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "ShellMITC4 element: " << ShellMITC4Id << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4+eleArgStart], &kNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "ShellMITC4 element: " << ShellMITC4Id << endln;
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5+eleArgStart], &lNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "ShellMITC4 element: " << ShellMITC4Id << endln;
     return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6+eleArgStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matTag\n";
    opserr << "ShellMITC4 element: " << ShellMITC4Id << endln;
    return TCL_ERROR;
  }

  SectionForceDeformation *theSection = theTclBuilder->getSection(matID);      

  if (theSection == 0) {
    opserr << "WARNING section not found\n";
    opserr << "section tag: " << matID;
    opserr << "\nShellMITC4 element: " << ShellMITC4Id << endln;
    return TCL_ERROR;
  }
  
  // now create the ShellMITC4 and add it to the Domain
  ShellMITC4 *theShellMITC4 = new ShellMITC4(ShellMITC4Id,iNode,jNode,kNode,lNode,*theSection);
  if (theShellMITC4 == 0) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "ShellMITC4 element: " << ShellMITC4Id << endln;
    return TCL_ERROR;
  }

  if (theTclDomain->addElement(theShellMITC4) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "ShellMITC4 element: " << ShellMITC4Id << endln;
    delete theShellMITC4;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}

