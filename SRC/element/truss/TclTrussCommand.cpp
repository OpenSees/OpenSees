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
                                                                        
// $Revision: 1.6 $
// $Date: 2003-03-12 19:20:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/truss/TclTrussCommand.cpp,v $
                                                                        
                                                                        
// File: ~/element/TclTrussCommand.C
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the implementation of the TclModelBuilder_addTruss()
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <Truss.h>
#include <TrussSection.h>
#include <TclModelBuilder.h>
#include <CorotTruss.h>
#include <CorotTrussSection.h>


extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addTruss(ClientData clientData, Tcl_Interp *interp,  int argc, 
			 TCL_Char **argv, Domain*theTclDomain,
			 TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 5) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element truss eleTag? iNode? jNode? A? matTag?\n";
    opserr << " -or- element truss eleTag? iNode? jNode? sectTag?" << endln;
    return TCL_ERROR;
  }    

  int ndm = theTclBuilder->getNDM();

  // get the id and end nodes 
  int trussId, iNode, jNode, matID;
  double A;
  double rho = 0.0;
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &trussId) != TCL_OK) {
    opserr << "WARNING invalid truss eleTag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "truss element: " << trussId << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
     opserr << "WARNING invalid jNode\n";
     opserr << "truss element: " << trussId << endln;
     return TCL_ERROR;
  }

  for (int i = 4+eleArgStart; i < argc; i++) {
    if (i+1 < argc && strcmp(argv[i], "-rho") == 0) {
      if (Tcl_GetDouble(interp, argv[i+1], &rho) != TCL_OK) {
	opserr << "WARNING invalid rho\n";
	opserr << "truss element: " << trussId << endln;
	return TCL_ERROR;
      }
      argc -= 2;
      break;
    }
  }

  if ((argc-eleArgStart) == 6) {  
      if (Tcl_GetDouble(interp, argv[4+eleArgStart], &A) != TCL_OK) {
	  opserr << "WARNING invalid A\n";
	  opserr << "truss element: " << trussId << endln;
	  return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[5+eleArgStart], &matID) != TCL_OK) {
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
      Element *theTruss = 0;
      if (strcmp(argv[eleArgStart],"corotTruss") == 0)
          theTruss = new CorotTruss(trussId,ndm,iNode,jNode,*theMaterial,A,rho);
      else
          theTruss = new Truss(trussId,ndm,iNode,jNode,*theMaterial,A,rho);

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
  } else {
      int sectID;
      if (Tcl_GetInt(interp, argv[4+eleArgStart], &sectID) != TCL_OK) {
	  opserr << "WARNING invalid matTag\n";
	  opserr << "truss element: " << trussId << endln;
	  return TCL_ERROR;
      }
      
      SectionForceDeformation *theSection = theTclBuilder->getSection(sectID);

      if (theSection == 0) {
	      opserr << "WARNING section not found\n";
	      opserr << "section tag: " << sectID;
	      opserr << "\ntruss element: " << trussId << endln;
	      return TCL_ERROR;
      }
  
      // now create the truss and add it to the Domain
      Element *theTruss = 0;
      if (strcmp(argv[eleArgStart],"corotTruss") == 0)
          theTruss = new CorotTrussSection(trussId,ndm,iNode,jNode,*theSection,rho);
      else
          theTruss = new TrussSection(trussId,ndm,iNode,jNode,*theSection,rho);
      
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
  }
      

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



