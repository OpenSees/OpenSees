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
// $Date: 2009-01-08 22:00:17 $
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
TclCommand_HyperbolicGapMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

  int tag;
  double Kmax, Kur, Rf, Fult, gap;
  UniaxialMaterial *theMaterial = 0;

  if (argc < 8) {
    opserr << "WARNING insufficient number of arguments\n";
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[3], &Kmax) != TCL_OK) {
    opserr << "WARNING invalid Kmax\n";
    return TCL_ERROR;	
  }

  if (Tcl_GetDouble(interp, argv[4], &Kur) != TCL_OK) {
    opserr << "WARNING invalid Kur\n";
    return TCL_ERROR;	
  }

  if (Tcl_GetDouble(interp, argv[5], &Rf) != TCL_OK) {
    opserr << "WARNING invalid Rf\n";
    return TCL_ERROR;	
  }

  if (Tcl_GetDouble(interp, argv[6], &Fult) != TCL_OK) {
    opserr << "WARNING invalid Fult\n";
    return TCL_ERROR;	
  }

  if (Tcl_GetDouble(interp, argv[7], &gap) != TCL_OK) {
    opserr << "WARNING invalid gap\n";
    return TCL_ERROR;	
  }
  
  theMaterial = new HyperbolicGapMaterial(tag, Kmax, Kur, Rf, Fult, gap);

  if (theMaterial != 0)  {
    if (OPS_addUniaxialMaterial(theMaterial) == true)
      return 0;
    else
      return -1;
  }

  return -1;
}
