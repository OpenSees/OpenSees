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
// $Date: 2008-04-15 18:29:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/TclHyperbolicGapMaterial.cpp,v $

// Written: MD
// Created: April 2008
//

#include <TclModelBuilder.h>
#include <HyperbolicGapMaterial.h>

#include <Vector.h>
#include <string.h>
#include <tcl.h>

int
TclCommand_HyperbolicGapMaterial(ClientData clientData, Tcl_Interp *interp, int argc, 
				 TCL_Char **argv, TclModelBuilder *theTclBuilder)
{

  int tag;
  double Kmax, Kur, Rf, Fult, gap;
  UniaxialMaterial *theMaterial = 0;

  if (argc < 8) {
    opserr << "WARNING insufficient number of arguments\n";
    return 0;
  }
  
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    return 0;
  }

  if (Tcl_GetDouble(interp, argv[3], &Kmax) != TCL_OK) {
    opserr << "WARNING invalid Kmax\n";
    return 0;	
  }

  if (Tcl_GetDouble(interp, argv[4], &Kur) != TCL_OK) {
    opserr << "WARNING invalid Kur\n";
    return 0;	
  }

  if (Tcl_GetDouble(interp, argv[5], &Rf) != TCL_OK) {
    opserr << "WARNING invalid Rf\n";
    return 0;	
  }

  if (Tcl_GetDouble(interp, argv[6], &Fult) != TCL_OK) {
    opserr << "WARNING invalid Fult\n";
    return 0;	
  }

  if (Tcl_GetDouble(interp, argv[7], &gap) != TCL_OK) {
    opserr << "WARNING invalid gap\n";
    return 0;	
  }
  
  theMaterial = new HyperbolicGapMaterial(tag, Kmax, Kur, Rf, Fult, gap);

  if (theMaterial != 0) 
    return theTclBuilder->addUniaxialMaterial(*theMaterial);
  else
    return -1;
}
