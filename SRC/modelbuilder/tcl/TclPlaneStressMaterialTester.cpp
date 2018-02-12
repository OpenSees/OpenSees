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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-25 23:34:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclUniaxialMaterialTester.cpp,v $
                                                                        
// File: ~/modelbuilder/tcl/TclUniaxialMaterialTester.C
// 
// Written: fmk 
// Created: 03/01
//
// Description: This file contains the implementaion of the TclUniaxialMaterialTester class.
//
// What: "@(#) TclUniaxialMaterialTester.C, revA"

#include <stdlib.h>
#include <string.h>

#include <ArrayOfTaggedObjects.h>
#include <NDMaterial.h>
#include <TclPlaneStressMaterialTester.h>
#include <Matrix.h>

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static TclPlaneStressMaterialTester *theTclBuilder =0;
static NDMaterial *theMaterial =0;

// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int  TclPlaneStressMaterialTester_setPlaneStressMaterial(ClientData clientData,
							 Tcl_Interp *interp, 
							 int argc,   
							 TCL_Char **argv);

int  TclPlaneStressMaterialTester_setStrainPlaneStressMaterial(ClientData clientData, 
							       Tcl_Interp *interp,
							       int argc,   
							       TCL_Char **argv);

int  TclPlaneStressMaterialTester_getStressPlaneStressMaterial(ClientData clientData,
							       Tcl_Interp *interp,
							       int argc,
							       TCL_Char **argv);


int  TclPlaneStressMaterialTester_getTangPlaneStressMaterial(ClientData clientData, 
							     Tcl_Interp *interp,
							     int argc,   
							     TCL_Char **argv);

//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

static int count;
static int countsTillCommit;
				    
// constructor: the constructor will add certain commands to the interpreter
TclPlaneStressMaterialTester::TclPlaneStressMaterialTester(Domain &theDomain, 
							   Tcl_Interp *interp, 
							   int cTC)
  :TclModelBuilder(theDomain, interp, 1, 1), theInterp(interp)
{
  countsTillCommit = cTC;
  Tcl_CreateCommand(interp, "setMaterial", TclPlaneStressMaterialTester_setPlaneStressMaterial,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "setStrain", TclPlaneStressMaterialTester_setStrainPlaneStressMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "getStress", TclPlaneStressMaterialTester_getStressPlaneStressMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "getTangent", TclPlaneStressMaterialTester_getTangPlaneStressMaterial,
		    (ClientData)NULL, NULL);
  
  
  // set the static pointers in this file
  theTclBuilder = this;
}

TclPlaneStressMaterialTester::~TclPlaneStressMaterialTester()
{

  theTclBuilder =0;
  opserr << "REMOVING COMMANDS";
  Tcl_DeleteCommand(theInterp, "setMaterial");
  Tcl_DeleteCommand(theInterp, "setStrain");
  Tcl_DeleteCommand(theInterp, "getStress");
  Tcl_DeleteCommand(theInterp, "getTangent");
}


//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int
TclPlaneStressMaterialTester_setPlaneStressMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
					      TCL_Char **argv)
{
  count = 1;
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }
  
  // check number of arguments in command line
  if (argc < 2) {
    opserr << "WARNING bad command - want: plane stressTest matID?\n";
    return TCL_ERROR;
  }    
  
  // get the matID form command line
  int matID;
  if (Tcl_GetInt(interp, argv[1], &matID) != TCL_OK) {
    opserr <<  "WARNING could not read matID: plane stressTest matID?\n";
    return TCL_ERROR;
  }
  
  // delete the old testing material
  if (theMaterial !=0) {
    delete theMaterial;
    theMaterial = 0;
  }

  // get the material from the modelbuilder with matID 
  // and set the testing material to point to a copy of it
  NDMaterial *theOrigMaterial = OPS_getNDMaterial(matID);
  if (theOrigMaterial == 0) {
    opserr << "WARNING no material found with matID\n";
    return TCL_ERROR;
  }  else {
    theMaterial = theOrigMaterial->getCopy("PlaneStress");
  }

  return TCL_OK;
}


int  
TclPlaneStressMaterialTester_setStrainPlaneStressMaterial(ClientData clientData, Tcl_Interp *interp,
						    int argc,   TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }
  
  // check number of arguments in command line
  if (argc < 4) {
    opserr <<  "WARNING bad command - want: strainPlaneStressTest strain?\n";
    return TCL_ERROR;
  }    

  // get the matID form command line
  static double strain[3];
  static Vector strainV(strain,3);
  if (Tcl_GetDouble(interp, argv[1], &strain[0]) != TCL_OK) {
    opserr <<  "WARNING could not read strain: strainPlaneStressTest strain?\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &strain[1]) != TCL_OK) {
    opserr << "WARNING could not read strain: strainPlaneStressTest strain?\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &strain[2]) != TCL_OK) {
    opserr << "WARNING could not read strain: strainPlaneStressTest strain?\n";
    return TCL_ERROR;
  }

  // delete the old testing material
  if (theMaterial !=0) {
    theMaterial->setTrialStrain(strainV);
    if (count == countsTillCommit) {
      theMaterial->commitState();    
      count = 1;
    } else count++;
  }
  return TCL_OK;
}



int  TclPlaneStressMaterialTester_getStressPlaneStressMaterial(ClientData clientData, Tcl_Interp *interp,
							 int argc,   TCL_Char **argv)
{
  static Vector stress(3);

  // delete the old testing material
  if (theMaterial !=0) {
    stress = theMaterial->getStress();
    char buffer[60];
    sprintf(buffer,"%.10e %.10e %.10e",stress(0),stress(1),stress(2));
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  } else {
    opserr << "WARNING no active Plane StressMaterial - use plane stressTest command\n";
    return TCL_ERROR;
  }
}

int  TclPlaneStressMaterialTester_getTangPlaneStressMaterial(ClientData clientData, Tcl_Interp *interp,
						       int argc,   TCL_Char **argv)
{
  static Matrix tangent(3,3);

  // delete the old testing material
  if (theMaterial !=0) {
    tangent = theMaterial->getTangent();
    char buffer[180];
    sprintf(buffer,"%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e",tangent(0,0), tangent(0,1), tangent(0,2), tangent(1,0), tangent(1,1), tangent(1,2), tangent(2,0),tangent(2,1),tangent(2,2));
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  } else {
    opserr << "WARNING no active PlaneStressMaterial - use plane stressTest command\n"; 
    return TCL_ERROR;
  }
}

  



