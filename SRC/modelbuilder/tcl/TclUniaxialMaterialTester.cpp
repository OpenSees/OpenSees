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
#include <UniaxialMaterial.h>
#include <TclUniaxialMaterialTester.h>

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static TclUniaxialMaterialTester *theTclBuilder =0;
static UniaxialMaterial *theTestingUniaxialMaterial =0;

// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int  TclUniaxialMaterialTester_setUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, 
						   int argc,   TCL_Char **argv);
				    
int  TclUniaxialMaterialTester_setStrainUniaxialMaterial(ClientData clientData, Tcl_Interp *interp,
							 int argc,   TCL_Char **argv);

int  TclUniaxialMaterialTester_getStressUniaxialMaterial(ClientData clientData, Tcl_Interp *interp,
							 int argc,   TCL_Char **argv);


int  TclUniaxialMaterialTester_getTangUniaxialMaterial(ClientData clientData, Tcl_Interp *interp,
						       int argc,   TCL_Char **argv);

//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

static int count;
static int countsTillCommit;
				    
// constructor: the constructor will add certain commands to the interpreter
TclUniaxialMaterialTester::TclUniaxialMaterialTester(Domain &theDomain, Tcl_Interp *interp, int cTC)
  :TclModelBuilder(theDomain, interp, 1, 1), theInterp(interp)
{
  countsTillCommit = cTC;
  Tcl_CreateCommand(interp, "uniaxialTest", TclUniaxialMaterialTester_setUniaxialMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "strainUniaxialTest", TclUniaxialMaterialTester_setStrainUniaxialMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "stressUniaxialTest", TclUniaxialMaterialTester_getStressUniaxialMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "tangUniaxialTest", TclUniaxialMaterialTester_getTangUniaxialMaterial,
		    (ClientData)NULL, NULL);
  
  
  // set the static pointers in this file
  theTclBuilder = this;
}

TclUniaxialMaterialTester::~TclUniaxialMaterialTester()
{

  theTclBuilder =0;

  Tcl_DeleteCommand(theInterp, "uniaxialTest");
  Tcl_DeleteCommand(theInterp, "strainUniaxialTest");
  Tcl_DeleteCommand(theInterp, "stressUniaxialTest");
  Tcl_DeleteCommand(theInterp, "tangUniaxialTest");
}


//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int
TclUniaxialMaterialTester_setUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
					      TCL_Char **argv)
{
  count = 1;
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    Tcl_SetResult(interp, "WARNING builder has been destroyed", TCL_STATIC);
    return TCL_ERROR;
  }

  // check number of arguments in command line
  if (argc < 2) {
    Tcl_SetResult(interp, "WARNING bad command - want: uniaxialTest matID?", TCL_STATIC);
    return TCL_ERROR;
  }    

  // get the matID form command line
  int matID;
  if (Tcl_GetInt(interp, argv[1], &matID) != TCL_OK) {
    Tcl_SetResult(interp, "WARNING could not read matID: uniaxialTest matID?", TCL_STATIC);
    return TCL_ERROR;
  }

  // delete the old testing material
  if (theTestingUniaxialMaterial !=0) {
    delete theTestingUniaxialMaterial;
    theTestingUniaxialMaterial = 0;
  }

  // get the material from the modelbuilder with matID 
  // and set the testing material to point to a copy of it
  UniaxialMaterial *theOrigMaterial = theTclBuilder->getUniaxialMaterial(matID);
  if (theOrigMaterial == 0) {
    Tcl_SetResult(interp, "WARNING no material found with matID", TCL_STATIC);
    return TCL_ERROR;
  }  else {
    theTestingUniaxialMaterial = theOrigMaterial->getCopy();
  }

  return TCL_OK;
}


int  
TclUniaxialMaterialTester_setStrainUniaxialMaterial(ClientData clientData, Tcl_Interp *interp,
						    int argc,   TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    Tcl_SetResult(interp, "WARNING builder has been destroyed", TCL_STATIC);
    return TCL_ERROR;
  }

  // check number of arguments in command line
  if (argc < 2) {
    Tcl_SetResult(interp, "WARNING bad command - want: strainUniaxialTest strain?", TCL_STATIC);
    return TCL_ERROR;
  }    

  // get the matID form command line
  double strain;
  if (Tcl_GetDouble(interp, argv[1], &strain) != TCL_OK) {
    Tcl_SetResult(interp, "WARNING could not read strain: strainUniaxialTest strain?", TCL_STATIC);
    return TCL_ERROR;
  }

  // delete the old testing material
  if (theTestingUniaxialMaterial !=0) {
    theTestingUniaxialMaterial->setTrialStrain(strain);
    if (count == countsTillCommit) {
      theTestingUniaxialMaterial->commitState();    
      count = 1;
    } else count++;
  }
  return TCL_OK;
}



int  TclUniaxialMaterialTester_getStressUniaxialMaterial(ClientData clientData, Tcl_Interp *interp,
							 int argc,   TCL_Char **argv)
{
  double stress = 0.0;

  // delete the old testing material
  if (theTestingUniaxialMaterial !=0) {
    stress = theTestingUniaxialMaterial->getStress();
    sprintf(interp->result,"%.10e",stress);
    return TCL_OK;
  } else {
    Tcl_SetResult(interp, "WARNING no active UniaxialMaterial - use uniaxialTest command", TCL_STATIC);    
    return TCL_ERROR;
  }
}

int  TclUniaxialMaterialTester_getTangUniaxialMaterial(ClientData clientData, Tcl_Interp *interp,
						       int argc,   TCL_Char **argv)
{
  double tangent = 0.0;

  // delete the old testing material
  if (theTestingUniaxialMaterial !=0) {
    tangent = theTestingUniaxialMaterial->getTangent();
    sprintf(interp->result,"%.10e",tangent);
    return TCL_OK;
  } else {
    Tcl_SetResult(interp, "WARNING no active UniaxialMaterial - use uniaxialTest command", TCL_STATIC);    
    return TCL_ERROR;
  }
}

  



