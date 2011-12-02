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
// $Date: 2000-09-15 08:23:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/TclFourNodeQuadCommand.cpp,v $
                                                                        
// File: ~/element/TclFourNodeQuadCommand.C
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the implementation of the TclModelBuilder_addFourNodeQuad()
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <FourNodeQuad.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, char **argv);

int
TclModelBuilder_addFourNodeQuad(ClientData clientData, Tcl_Interp *interp,  
				int argc, 
				char **argv, 
				Domain*theTclDomain,
				TclModelBuilder *theTclBuilder)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 2) {
		cerr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
		return TCL_ERROR;
	}

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc-argStart) < 8) {
    cerr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    cerr << "Want: element FourNodeQuad eleTag? iNode? jNode? kNode? lNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int FourNodeQuadId, iNode, jNode, kNode, lNode, matID;
  double thickness;
	double p = 0.0;		// uniform normal traction (pressure)
	double r = 0.0;		// mass density
	double b1 = 0.0;
	double b2 = 0.0;
  char *type;
  if (Tcl_GetInt(interp, argv[argStart], &FourNodeQuadId) != TCL_OK) {
    cerr << "WARNING invalid FourNodeQuad eleTag" << endl;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[1+argStart], &iNode) != TCL_OK) {
    cerr << "WARNING invalid iNode\n";
    cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2+argStart], &jNode) != TCL_OK) {
     cerr << "WARNING invalid jNode\n";
     cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
     return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[3+argStart], &kNode) != TCL_OK) {
     cerr << "WARNING invalid kNode\n";
     cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
     return TCL_ERROR;
  }  
  
  if (Tcl_GetInt(interp, argv[4+argStart], &lNode) != TCL_OK) {
     cerr << "WARNING invalid lNode\n";
     cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
     return TCL_ERROR;
  }  

  if (Tcl_GetDouble(interp, argv[5+argStart], &thickness) != TCL_OK) {
     cerr << "WARNING invalid thickness\n";
     cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
     return TCL_ERROR;
  }  
  
  type = argv[6+argStart];
  
  if (Tcl_GetInt(interp, argv[7+argStart], &matID) != TCL_OK) {
     cerr << "WARNING invalid matID\n";
     cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
     return TCL_ERROR;
  }

	if ((argc-argStart) > 11) {
		if (Tcl_GetDouble(interp, argv[8+argStart], &p) != TCL_OK) {
			cerr << "WARNING invalid pressure\n";
			cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[9+argStart], &r) != TCL_OK) {
			cerr << "WARNING invalid rho\n";
			cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[10+argStart], &b1) != TCL_OK) {
			cerr << "WARNING invalid b1\n";
			cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[11+argStart], &b2) != TCL_OK) {
			cerr << "WARNING invalid b2\n";
			cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
			return TCL_ERROR;
		}

	}

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
      
  if (theMaterial == 0) {
      cerr << "WARNING material not found\n";
      cerr << "Material: " << matID;
      cerr << "\nFourNodeQuad element: " << FourNodeQuadId << endl;
      return TCL_ERROR;
  }
  
  // now create the FourNodeQuad and add it to the Domain
  FourNodeQuad *theFourNodeQuad = 
      new FourNodeQuad(FourNodeQuadId,iNode,jNode,kNode,lNode,
		       *theMaterial, type, thickness, p, r, b1, b2);
  if (theFourNodeQuad == 0) {
      cerr << "WARNING ran out of memory creating element\n";
      cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theFourNodeQuad) == false) {
      cerr << "WARNING could not add element to the domain\n";
      cerr << "FourNodeQuad element: " << FourNodeQuadId << endl;
      delete theFourNodeQuad;
      return TCL_ERROR;
  }

  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}



