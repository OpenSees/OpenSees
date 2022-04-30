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
// $Date: 2008-12-09 21:24:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclSectionTester.cpp,v $
                                                                        
// Written: fmk 
// Created: 10/08
//
// Description: This file contains the implementation of the TclUniaxialMaterialTester class.
//
// What: "@(#) TclUniaxialMaterialTester.C, revA"

#include <stdlib.h>
#include <string.h>

#include <ArrayOfTaggedObjects.h>
#include <SectionForceDeformation.h>
#include <TclSectionTester.h>
#include <Vector.h>
#include <DummyStream.h>
#include <Response.h>

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static TclSectionTester *theTclBuilder =0;
static SectionForceDeformation *theTestingSection =0;

// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int  TclSectionTester_setSection(ClientData clientData, Tcl_Interp *interp, 
				 int argc,   TCL_Char **argv);

int  TclSectionTester_setStrainSection(ClientData clientData, Tcl_Interp *interp,
				       int argc,   TCL_Char **argv);

int  TclSectionTester_getStressSection(ClientData clientData, Tcl_Interp *interp,
				       int argc,   TCL_Char **argv);

int  TclSectionTester_getTangSection(ClientData clientData, Tcl_Interp *interp,
				     int argc,   TCL_Char **argv);

int  TclSectionTester_getResponseSection(ClientData clientData, Tcl_Interp* interp,
    int argc, TCL_Char** argv);

//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

static int count;
static int countsTillCommit;

				    
// constructor: the constructor will add certain commands to the interpreter
TclSectionTester::TclSectionTester(Domain &theDomain, Tcl_Interp *interp, int cTC)
  :TclModelBuilder(theDomain, interp, 3, 6), theInterp(interp)
{
  countsTillCommit = cTC;
  Tcl_CreateCommand(interp, "sectionTest", TclSectionTester_setSection,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "strainSectionTest", TclSectionTester_setStrainSection,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "stressSectionTest", TclSectionTester_getStressSection,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "tangSectionTest", TclSectionTester_getTangSection,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "responseSectionTest", TclSectionTester_getResponseSection,
      (ClientData)NULL, NULL);
  
  // set the static pointers in this file
  theTclBuilder = this;
}

TclSectionTester::~TclSectionTester()
{

  theTclBuilder =0;

  Tcl_DeleteCommand(theInterp, "sectionTest");
  Tcl_DeleteCommand(theInterp, "strainSectionTest");
  Tcl_DeleteCommand(theInterp, "stressSectionTest");
  Tcl_DeleteCommand(theInterp, "tangSectionTest");
  Tcl_DeleteCommand(theInterp, "responseSectionTest");
}


//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int
TclSectionTester_setSection(ClientData clientData, Tcl_Interp *interp, int argc,   
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
    opserr << "WARNING bad command - want: sectionTest secID?\n";
    return TCL_ERROR;
  }    

  // get the matID form command line
  int sectionID;
  if (Tcl_GetInt(interp, argv[1], &sectionID) != TCL_OK) {
    opserr << "WARNING could not read sectionID: sectionTest sectionID?\n";
    return TCL_ERROR;
  }

  // delete the old testing material
  if (theTestingSection !=0) {
    delete theTestingSection;
    theTestingSection = 0;
  }

  // get the material from the modelbuilder with sectionID 
  // and set the testing material to point to a copy of it
  SectionForceDeformation *theOrigMaterial = theTclBuilder->getSection(sectionID);
  if (theOrigMaterial == 0) {
    opserr << "WARNING no material found with sectionID\n";
    return TCL_ERROR;
  }  else {
    theTestingSection = theOrigMaterial->getCopy();
  }

  return TCL_OK;
}


int  
TclSectionTester_setStrainSection(ClientData clientData, Tcl_Interp *interp,
						    int argc,   TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }

  // check number of arguments in command line
  if (argc < 2) {
    opserr <<  "WARNING bad command - want: strainSectionTest strain?\n";
    return TCL_ERROR;
  }    

  // get the sectionID form command line
  // Need to set the data based on argc, otherwise it crashes when setting "data(i-1) = strain"
  static Vector data(argc-1);
  double strain;
  for (int i=1; i<argc; i++) {
    if (Tcl_GetDouble(interp, argv[i], &strain) != TCL_OK) {
      opserr << "WARNING could not read strain: strainSectionTest strain1? strain2? ... strainN?\n";
      return TCL_ERROR;
    }
    data(i-1) = strain;
  }

  // if the section exists, otherwise throw an error
  if (theTestingSection !=0) {
    theTestingSection->setTrialSectionDeformation(data);
    if (count == countsTillCommit) {
      theTestingSection->commitState();    
      count = 1;
    } else count++;
  }
  return TCL_OK;
}


int  TclSectionTester_getStressSection(ClientData clientData, Tcl_Interp *interp,
				       int argc,   TCL_Char **argv)
{
  // if the section exists, otherwise throw an error
  if (theTestingSection !=0) {
    const Vector &stress = theTestingSection->getStressResultant();
    for (int i=0; i<stress.Size(); i++) {
      char buffer[40];
      sprintf(buffer,"%.10e ",stress(i));
      Tcl_AppendResult(interp, buffer, NULL);
      //      sprintf(interp->result,"%.10e",stress(i));
    }
    return TCL_OK;
  } else {
    opserr <<  "WARNING no active Section - use sectionTest command\n";    
    return TCL_ERROR;
  }
}

int  TclSectionTester_getTangSection(ClientData clientData, Tcl_Interp *interp,
						       int argc,   TCL_Char **argv)
{

  // if the section exists, otherwise throw an error
  if (theTestingSection !=0) {
    const Matrix &tangent = theTestingSection->getSectionTangent();
    for (int i=0; i<tangent.noRows(); i++)
      for (int j=0; j<tangent.noCols(); j++) {
	char buffer[40];
	sprintf(buffer,"%.10e ",tangent(i,j));
	Tcl_AppendResult(interp, buffer, NULL);
	//	sprintf(interp->result,"%.10e",tangent(i,j));
      }
    return TCL_OK;
  } else {
    opserr << "WARNING no active Section - use sectionTest command\n";    
    return TCL_ERROR;
  }
}

int  TclSectionTester_getResponseSection(ClientData clientData, Tcl_Interp* interp,
    int argc, TCL_Char** argv)
{
    // if the section exists, otherwise throw an error
    if (theTestingSection != 0) {
        DummyStream dummy;
        Response* theResponse = theTestingSection->setResponse(argv+1, argc-1, dummy);
        if (theResponse == 0) {
            return TCL_ERROR;
        }
        if (theResponse->getResponse() < 0) {
            delete theResponse;
            return TCL_ERROR;
        } 
        Information &eleInfo = theResponse->getInformation();
        const Vector &data = eleInfo.getData();
        for (int i = 0; i < data.Size(); i++) {
            char buffer[40];
            sprintf(buffer, "%.10e ", data(i));
            Tcl_AppendResult(interp, buffer, NULL);
            //      sprintf(interp->result,"%.10e",stress(i));
        }
        // Now delete theResponse since I already retrieved the data
        delete theResponse;
        return TCL_OK;
    }
    else {
        opserr << "WARNING no active Section - use sectionTest command\n";
        return TCL_ERROR;
    }
}