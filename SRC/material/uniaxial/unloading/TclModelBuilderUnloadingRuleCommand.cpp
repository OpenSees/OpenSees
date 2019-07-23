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

#include <TclModelBuilder.h>

#include <TakedaUnloadingRule.h>
#include <EnergyUnloadingRule.h>
#include <ConstantUnloadingRule.h>

#include <string.h>

static void printCommand(int argc, TCL_Char **argv)
{
  opserr << "Input command: ";
  for (int i=0; i<argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
} 

int
TclModelBuilderUnloadingRuleCommand(ClientData clienData,
				    Tcl_Interp *interp,
				    int argc, TCL_Char **argv)
{
  // Make sure there is a minimum number of arguments
  if (argc < 2) {
    opserr << "WARNING insufficient number of unloadingRule arguments\n";
    opserr << "Want: unloadingRule type? tag? <specific unloadingRule args>" << endln;
    return TCL_ERROR;
  }
  
  // Pointer to a unloadingRule that will be added to the model builder
  UnloadingRule *theState = 0;
  
  // Check argv[1] for unloadingRule type	
  if (strcmp(argv[1],"Ductility") == 0 || strcmp(argv[1],"Takeda") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: unloadingRule Takeda tag? C? beta?" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    double C, b;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid unloadingRule Takeda tag" << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &C) != TCL_OK) {
      opserr << "WARNING invalid C\n";
      opserr << "unloadingRule Takeda: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &b) != TCL_OK) {
      opserr << "WARNING invalid beta\n";
      opserr << "unloadingRule Takeda: " << tag << endln;
      return TCL_ERROR;	
    }
    
    theState = new TakedaUnloadingRule (tag, C, b);
  }
  
  else if (strcmp(argv[1],"Energy") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: unloadingRule Energy tag? Et? c?" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    double c, et;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid unloadingRule Energy tag" << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &et) != TCL_OK) {
      opserr << "WARNING invalid Et\n";
      opserr << "unloadingRule Energy: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &c) != TCL_OK) {
      opserr << "WARNING invalid c\n";
      opserr << "unloadingRule Energy: " << tag << endln;
      return TCL_ERROR;	
    }
    
    theState = new EnergyUnloadingRule (tag, et, c);
  }
  
  else if (strcmp(argv[1],"Constant") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: unloadingRule Constant tag? alpha? beta?" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    double a, b;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid unloadingRule Constant tag" << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &a) != TCL_OK) {
      opserr << "WARNING invalid alpha\n";
      opserr << "unloadingRule Constant: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &b) != TCL_OK) {
      opserr << "WARNING invalid beta\n";
      opserr << "unloadingRule Constant: " << tag << endln;
      return TCL_ERROR;	
    }
    
    theState = new ConstantUnloadingRule (tag, a, b);
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
