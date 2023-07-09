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
// $Date: 2005-08-19 19:40:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/TclNewUniaxialMaterial.cpp,v $

// Written: fmk
// Created: Aug 2005
//
// Description: This file contains the implementation of the


#include <NewUniaxialMaterial.h>

#include <TclModelBuilder.h>
#include <string.h>

int
TclCommand_NewUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, int argc, 
				  TCL_Char **argv, TclModelBuilder *theTclBuilder)
{

  int tag;
  UniaxialMaterial *theMaterial = 0;

  if (argc < 3) {
    opserr << "WARNING insufficient number of arguments\n";
    return 0;
  }
  
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    return 0;
  }

  /*
  if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
    opserr << "WARNING invalid E\n";
    return 0;	
  }
  */
  
  theMaterial = new NewUniaxialMaterial(tag);

  if (theMaterial != 0) 
    return theTclBuilder->addUniaxialMaterial(*theMaterial);
  else
    return -1;
}
