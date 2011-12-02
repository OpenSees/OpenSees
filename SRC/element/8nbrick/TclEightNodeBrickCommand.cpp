///////////////////////////////////////////////////////////////////////////////
//
// COPYLEFT (C):     :-))
//``This  source code is Copyrighted in U.S., by the The Regents of the University
//of California, for an indefinite period, and anybody caught using it without our
//permission,  will  be  mighty  good friends of ourn, cause we don't give a darn.
//Hack  it.  Compile it. Debug it. Run it. Yodel it. Enjoy it. We wrote it, that's
//all we wanted to do.'' bj
// PROJECT:           Object Oriented Finite Element Program
// FILE:              EightNodeBrick.h
// CLASS:             EightNodeBrick
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class
// RETURN:
// VERSION:
// LANGUAGE:          C++.ver >= 3.0
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic, Zhaohui Yang and Xiaoyan Wu
// PROGRAMMER:        Boris Jeremic, Zhaohui Yang  and Xiaoyan Wu
// DATE:              Aug. 2000
// Description: This file contains the implementation of the TclModelBuilder_addEightNodeBrick() 
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <ErrorHandler.h>
#include <EightNodeBrick.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addEightNodeBrick(ClientData clientData, Tcl_Interp *interp,  int argc, 
				  TCL_Char **argv, Domain*theTclDomain,
				  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
      opserr << "command: element Brick8N - no modelbuilder\n";
      return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 15) {
    opserr << "command: element Brick8N - insufficient args - want: " <<
      "element Brick8N eleTag? node1? node2? .. node8? matTag? bforce1? bforce2? bforce3? massDensity?\n";
    return TCL_ERROR;
  }    

  // get the id and end nodes 
  int eleID, matID;
  int nodes[8];
  double bodyforces[3], massDensity;
  //char *type;
  
  // read the eleTag
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
    opserr << "command: element Brick8N - invalid integer tag " << argv[1+eleArgStart] << endln;
    return TCL_ERROR;
  }
  
  // read the 8 node tags
  int i;
  for (i=0; i<8; i++) {
      if (Tcl_GetInt(interp, argv[2+i+eleArgStart], &nodes[i]) != TCL_OK) {
	opserr << "command: element Brick8N " << eleID << " - invalid integer tag " << argv[2+i+eleArgStart] << endln;
	return TCL_ERROR;
      }
  }

  // read in material tag & check the material exists in the model builder
  if (Tcl_GetInt(interp, argv[10+eleArgStart], &matID) != TCL_OK) {
    opserr << "command: element Brick8N " << eleID << " - invalid matID tag " << argv[10+eleArgStart] << endln;
    return TCL_ERROR;
  }

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
  
  if (theMaterial == 0) {
    opserr << "command: element Brick8N " << eleID << " - no NDMaterial with tag " << argv[10+eleArgStart] << "exists\n";
    return TCL_ERROR;      
  }
  
  //type = argv[11+eleArgStart];

  // read the 3 bodyforce accel's 
  for (i=0; i<3; i++) {
      if (Tcl_GetDouble(interp, argv[11+i+eleArgStart], &bodyforces[i]) != TCL_OK) {
	opserr << "command: element Brick8N " << eleID << " - invalid bodyforces tag " << argv[11+i+eleArgStart] << endln;
	return TCL_ERROR;
      }
  }

  // now get the massDensity
  if (Tcl_GetDouble(interp, argv[14+eleArgStart], &massDensity) != TCL_OK) {
    opserr << "command: element Brick8N " << eleID << " - invalid massDensity " << argv[14+eleArgStart] << endln;
    return TCL_ERROR;
  }  
  
  // now create the EightNodeBrick and add it to the Domain
  EightNodeBrick *theEle = new EightNodeBrick(eleID,nodes[0], nodes[1], nodes[2], nodes[3], nodes[4],
                                              nodes[5],nodes[6], nodes[7], theMaterial, 
					      bodyforces[0], bodyforces[1], bodyforces[2], massDensity, 0.0);
					      
  if (theEle == 0) {
    opserr << "command: element Brick8N " << eleID << " - out of memory\n";      
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theEle) == false) {
    opserr << "command: element Brick8N " << eleID << "- could not add ele to domain\n"; 
    delete theEle;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



