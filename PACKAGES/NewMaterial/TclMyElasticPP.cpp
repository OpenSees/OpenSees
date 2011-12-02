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
// $Date: 2005-07-01 19:53:43 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/TclMyElasticPP.cpp,v $
                                                                        
// Written: fmk

// Description: This file contains the function invoked when the user invokes
// the MyElasticPP command in the interpreter. 
//
// What: "@(#) TclModelBuilderUniaxialMaterialCommand.C, revA"


#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <TclModelBuilder.h>
#include "MyElasticPP.h"

#ifdef _USRDLL
#define DllExport _declspec(dllexport)
#else
#define DllExport
#endif

extern "C" DllExport int
TclCommand_MyElasticPP(ClientData clientData, 
		       Tcl_Interp *interp,  
		       int argc, 
		       TCL_Char **argv, 
		       TclModelBuilder *theTclBuilder)
{
  
  // Pointer to a uniaxial material that will be added to the model builder
  UniaxialMaterial *theMaterial = 0;

  if (argc < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial ElasticPP tag? E? epsy?" << endln;
    return TCL_ERROR;
  }
  
  int tag;
  double E, ep;
  
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial ElasticPP tag" << endln;
    return TCL_ERROR;		
  }
  
  if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
    opserr << "WARNING invalid E\n";
    opserr << "uniaxialMaterial ElasticPP: " << tag << endln;
    return TCL_ERROR;	
  }
  
  if (Tcl_GetDouble(interp, argv[4], &ep) != TCL_OK) {
    opserr << "WARNING invalid epsy\n";
    opserr << "uniaxialMaterial ElasticPP: " << tag << endln;
    return TCL_ERROR;
  }
  
  theMaterial = new MyElasticPP(tag, E, ep);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial " << argv[1] << endln;
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (theTclBuilder->addUniaxialMaterial(*theMaterial) < 0) {
    opserr << "WARNING could not add uniaxialMaterial to the domain\n";
    opserr << *theMaterial << endln;
    delete theMaterial; // invoke the material objects destructor, otherwise mem leak
    return TCL_ERROR;
  }
    
  return TCL_OK;
}

