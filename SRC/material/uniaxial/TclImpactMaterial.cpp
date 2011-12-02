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
                                                                        
// $Revision: 1.2 $
// $Date: 2009-03-27 19:19:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/TclImpactMaterial.cpp,v $

// File: ~/material/ImpactMaterial.C
//
// Written: md
// Created: 06/2008
//


#include <TclModelBuilder.h>
#include <ImpactMaterial.h>

#include <Vector.h>
#include <string.h>
#include <tcl.h>

int
TclCommand_ImpactMaterial(ClientData clientData, Tcl_Interp *interp, int argc, 
				 TCL_Char **argv, TclModelBuilder *theTclBuilder)
{

  int tag;
  double K1, K2, Delta_y, gap;
  UniaxialMaterial *theMaterial = 0;

  if (argc < 7) {
    opserr << "WARNING insufficient number of arguments\n";
    return 0;
  }

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    return 0;
  }

  if (Tcl_GetDouble(interp, argv[3], &K1) != TCL_OK) {
    opserr << "WARNING invalid K1\n";
    return 0;	
  }

  if (Tcl_GetDouble(interp, argv[4], &K2) != TCL_OK) {
    opserr << "WARNING invalid K2\n";
    return 0;	
  }
  
  if (Tcl_GetDouble(interp, argv[5], &Delta_y) != TCL_OK) {
    opserr << "WARNING invalid Delta_y\n";
    return 0;	
  }

  if (Tcl_GetDouble(interp, argv[6], &gap) != TCL_OK) {
    opserr << "WARNING invalid gap\n";
    return 0;	
  }

  theMaterial = new ImpactMaterial(tag, K1, K2, Delta_y, gap);
  
  if (theMaterial != 0) 
    return theTclBuilder->addUniaxialMaterial(*theMaterial);
  else
    return -1;
}

