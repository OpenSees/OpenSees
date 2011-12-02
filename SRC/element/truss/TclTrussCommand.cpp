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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
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
#include <iostream.h>
#include <Domain.h>

#include <Truss.h>
#include <TrussSection.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addTruss(ClientData clientData, Tcl_Interp *interp,  int argc, 
			 char **argv, Domain*theTclDomain,
			 TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 5) {
    cerr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    cerr << "Want: element truss eleTag? iNode? jNode? A? matTag?\n";
    cerr << " -or- element truss eleTag? iNode? jNode? sectTag?" << endl;
    return TCL_ERROR;
  }    

  int ndm = theTclBuilder->getNDM();

  // get the id and end nodes 
  int trussId, iNode, jNode, matID;
  double A;
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &trussId) != TCL_OK) {
    cerr << "WARNING invalid truss eleTag" << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
    cerr << "WARNING invalid iNode\n";
    cerr << "truss element: " << trussId << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
     cerr << "WARNING invalid jNode\n";
     cerr << "truss element: " << trussId << endl;
     return TCL_ERROR;
  }

  if ((argc-eleArgStart) == 6) {  
      if (Tcl_GetDouble(interp, argv[4+eleArgStart], &A) != TCL_OK) {
	  cerr << "WARNING invalid A\n";
	  cerr << "truss element: " << trussId << endl;
	  return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[5+eleArgStart], &matID) != TCL_OK) {
	  cerr << "WARNING invalid matTag\n";
	  cerr << "truss element: " << trussId << endl;
	  return TCL_ERROR;
      }
      
      UniaxialMaterial *theMaterial = theTclBuilder->getUniaxialMaterial(matID);
      
      if (theMaterial == 0) {
	  cerr << "WARNING material not found\n";
	  cerr << "Material: " << matID;
	  cerr << "\ntruss element: " << trussId << endl;
	  return TCL_ERROR;
      }
  
      // now create the truss and add it to the Domain
      Truss *theTruss = new Truss(trussId,ndm,iNode,jNode,*theMaterial,A);
      if (theTruss == 0) {
	  cerr << "WARNING ran out of memory creating element\n";
	  cerr << "truss element: " << trussId << endl;
	  return TCL_ERROR;
      }

      if (theTclDomain->addElement(theTruss) == false) {
	  cerr << "WARNING could not add element to the domain\n";
	  cerr << "truss element: " << trussId << endl;
	  delete theTruss;
	  return TCL_ERROR;
      }
  } else {
      int sectID;
      if (Tcl_GetInt(interp, argv[4+eleArgStart], &sectID) != TCL_OK) {
	  cerr << "WARNING invalid matTag\n";
	  cerr << "truss element: " << trussId << endl;
	  return TCL_ERROR;
      }
      
      SectionForceDeformation *theSection = theTclBuilder->getSection(sectID);

      if (theSection == 0) {
	      cerr << "WARNING section not found\n";
	      cerr << "section tag: " << sectID;
	      cerr << "\ntruss element: " << trussId << endl;
	      return TCL_ERROR;
      }
  
      // now create the truss and add it to the Domain
      TrussSection *theTruss = 0;       

	  theTruss = new TrussSection(trussId,ndm,iNode,jNode,*theSection);
      
      if (theTruss == 0) {
	  cerr << "WARNING ran out of memory creating element\n";
	  cerr << "truss element: " << trussId << endl;
	  return TCL_ERROR;
      }

      if (theTclDomain->addElement(theTruss) == false) {
	  cerr << "WARNING could not add element to the domain\n";
	  cerr << "truss element: " << trussId << endl;
	  delete theTruss;
	  return TCL_ERROR;
      }
  }
      

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



