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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: August 2000
//
// Description: This file contains the parsing routines for the
// TCL stiffnessDegradation command.

#include <OPS_Globals.h>

#include <TclModelBuilder.h>

#include <DuctilityStiffnessDegradation.h>
#include <EnergyStiffnessDegradation.h>
#include <ConstantStiffnessDegradation.h>

#include <string.h>

extern void *OPS_DuctilityStiffnessDegradation(void);
extern void *OPS_EnergyStiffnessDegradation(void);
extern void *OPS_ConstantStiffnessDegradation(void);
extern void *OPS_PincheiraStiffnessDegradation(void);

#include <elementAPI.h>
#include <packages.h>

int
TclModelBuilderStiffnessDegradationCommand(ClientData clientData,
					   Tcl_Interp *interp,
					   int argc, TCL_Char **argv, Domain *theDomain)
{
  // Make sure there is a minimum number of arguments
  if (argc < 2) {
    opserr << "WARNING insufficient number of stiffnessDegradation arguments\n";
    opserr << "Want: stiffnessDegradation type? tag? <specific stiffnessDegradation args>" << endln;
    return TCL_ERROR;
  }
  
  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);	  

  // Pointer to a stiffnessDegradation that will be added to the model builder
  StiffnessDegradation *theState = 0;
  
  // Check argv[1] for stiffnessDegradation type	
  if (strcmp(argv[1],"Ductility") == 0) {
      void *theDegr = OPS_DuctilityStiffnessDegradation();
      if (theDegr != 0) 
	theState = (StiffnessDegradation *)theDegr;
      else 
	return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"Energy") == 0) {
      void *theDegr = OPS_EnergyStiffnessDegradation();
      if (theDegr != 0) 
	theState = (StiffnessDegradation *)theDegr;
      else 
	return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"Constant") == 0) {
      void *theDegr = OPS_ConstantStiffnessDegradation();
      if (theDegr != 0) 
	theState = (StiffnessDegradation *)theDegr;
      else 
	return TCL_ERROR;
  }

  else if (strcmp(argv[1],"Pincheira") == 0) {
      void *theDegr = OPS_PincheiraStiffnessDegradation();
      if (theDegr != 0) 
	theState = (StiffnessDegradation *)theDegr;
      else 
	return TCL_ERROR;
  }
  
  else  {
    opserr << "WARNING unknown type of stiffnessDegradation: " << argv[1];
    opserr << "\nValid types: Ductility, Energy, Constant\n";
    return TCL_ERROR;
  }
  
  // Ensure we have created the Degradation, out of memory if got here and no stiffnessDegradation
  if (theState == 0) {
    opserr << "WARNING ran out of memory creating stiffnessDegradation\n";
    opserr << argv[1] << endln;
    return TCL_ERROR;
  }
  
  // Now add the material to the modelBuilder
  if (OPS_addStiffnessDegradation(theState) == false) {
    opserr << "WARNING could not add stiffnessDegradation to the domain\n";
    opserr << *theState << endln;
    delete theState;	// Avoid memory leak
    return TCL_ERROR;
  }
  
  return TCL_OK;
}
