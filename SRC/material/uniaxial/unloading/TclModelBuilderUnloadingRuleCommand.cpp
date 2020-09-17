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
// TCL unloadingRule command.

#include <OPS_Globals.h>

#include <TakedaUnloadingRule.h>
#include <EnergyUnloadingRule.h>
#include <ConstantUnloadingRule.h>

#include <string.h>

extern void *OPS_TakedaUnloadingRule(void);
extern void *OPS_EnergyUnloadingRule(void);
extern void *OPS_ConstantUnloadingRule(void);
extern void *OPS_KarsanUnloadingRule(void);

#include <elementAPI.h>
#include <packages.h>

int
TclModelBuilderUnloadingRuleCommand(ClientData clientData,
				    Tcl_Interp *interp,
				    int argc, TCL_Char **argv, Domain *theDomain)
{
  // Make sure there is a minimum number of arguments
  if (argc < 2) {
    opserr << "WARNING insufficient number of unloadingRule arguments\n";
    opserr << "Want: unloadingRule type? tag? <specific unloadingRule args>" << endln;
    return TCL_ERROR;
  }
  
  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);	  

  // Pointer to a unloadingRule that will be added to the model builder
  UnloadingRule *theState = 0;
  
  // Check argv[1] for unloadingRule type	
  if (strcmp(argv[1],"Ductility") == 0 || strcmp(argv[1],"Takeda") == 0) {
    void *theDegr = OPS_TakedaUnloadingRule();
    if (theDegr != 0) 
      theState = (UnloadingRule *)theDegr;
    else 
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"Energy") == 0) {
    void *theDegr = OPS_EnergyUnloadingRule();
    if (theDegr != 0) 
      theState = (UnloadingRule *)theDegr;
    else 
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"Constant") == 0) {
    void *theDegr = OPS_ConstantUnloadingRule();
    if (theDegr != 0) 
      theState = (UnloadingRule *)theDegr;
    else 
      return TCL_ERROR;
  }

  else if (strcmp(argv[1],"Karsan") == 0) {
    void *theDegr = OPS_KarsanUnloadingRule();
    if (theDegr != 0) 
      theState = (UnloadingRule *)theDegr;
    else 
      return TCL_ERROR;
  }
  
  else  {
    opserr << "WARNING unknown type of unloadingRule: " << argv[1];
    opserr << "\nValid types: Ductility, Energy, Constant\n";
    return TCL_ERROR;
  }
  
  // Ensure we have created the Degradation, out of memory if got here and no unloadingRule
  if (theState == 0) {
    opserr << "WARNING ran out of memory creating unloadingRule\n";
    opserr << argv[1] << endln;
    return TCL_ERROR;
  }
  
  // Now add the material to the modelBuilder
  if (OPS_addUnloadingRule(theState) == false) {
    opserr << "WARNING could not add unloadingRule to the domain\n";
    opserr << *theState << endln;
    delete theState;	// Avoid memory leak
    return TCL_ERROR;
  }
  
  return TCL_OK;
}
