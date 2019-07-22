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

static void printCommand(int argc, TCL_Char **argv)
{
  opserr << "Input command: ";
  for (int i=0; i<argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
} 

int
TclModelBuilderStiffnessDegradationCommand(ClientData clienData,
					   Tcl_Interp *interp,
					   int argc, TCL_Char **argv)
{
  // Make sure there is a minimum number of arguments
  if (argc < 2) {
    opserr << "WARNING insufficient number of stiffnessDegradation arguments\n";
    opserr << "Want: stiffnessDegradation type? tag? <specific stiffnessDegradation args>" << endln;
    return TCL_ERROR;
  }
  
  // Pointer to a stiffnessDegradation that will be added to the model builder
  StiffnessDegradation *theState = 0;
  
  // Check argv[1] for stiffnessDegradation type	
  if (strcmp(argv[1],"Ductility") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: stiffnessDegradation Ductility tag? alpha? beta?" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    double a, b;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid stiffnessDegradation Ductility tag" << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &a) != TCL_OK) {
      opserr << "WARNING invalid alpha\n";
      opserr << "stiffnessDegradation Ductility: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &b) != TCL_OK) {
      opserr << "WARNING invalid beta\n";
      opserr << "stiffnessDegradation Ductility: " << tag << endln;
      return TCL_ERROR;	
    }
    
    theState = new DuctilityStiffnessDegradation (tag, a, b);
  }
  
  else if (strcmp(argv[1],"Energy") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: stiffnessDegradation Energy tag? Et? c?" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    double c, et;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid stiffnessDegradation Energy tag" << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &et) != TCL_OK) {
      opserr << "WARNING invalid Et\n";
      opserr << "stiffnessDegradation Energy: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &c) != TCL_OK) {
      opserr << "WARNING invalid c\n";
      opserr << "stiffnessDegradation Energy: " << tag << endln;
      return TCL_ERROR;	
    }
    
    theState = new EnergyStiffnessDegradation (tag, et, c);
  }
  
  else if (strcmp(argv[1],"Constant") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: stiffnessDegradation Constant tag? alpha? beta?" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    double a, b;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid stiffnessDegradation Constant tag" << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &a) != TCL_OK) {
      opserr << "WARNING invalid alpha\n";
      opserr << "stiffnessDegradation Constant: " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &b) != TCL_OK) {
      opserr << "WARNING invalid beta\n";
      opserr << "stiffnessDegradation Constant: " << tag << endln;
      return TCL_ERROR;	
    }
    
    theState = new ConstantStiffnessDegradation (tag, a, b);
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
