///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              TwentyNodeBrick.cpp
// CLASS:             TwentyNodeBrick
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic, Zhaohui Yang and Xiaoyan Wu
// PROGRAMMER:        Boris Jeremic, Zhaohui Yang  and Xiaoyan Wu
// DATE:              Aug. 2001
// UPDATE HISTORY:
// Description: This file contains the implementation of the TclModelBuilder_addEightNodeBrick() 
// command. 
//

#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>
#include <Domain.h>

#include <ErrorHandler.h>
#include <TwentyNodeBrick.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addTwentyNodeBrick(ClientData clientData, Tcl_Interp *interp,  int argc, 
				  TCL_Char **argv, Domain*theTclDomain,
				  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
      opserr << "command: element Brick20N - no modelbuilder\n";
      return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 27) {
    opserr << "command: element Brick20N - insufficient args - want " << 
      "element Brick20N eleTag? node1? node2? .. node20?  matTag? bforce1? bforce2? bforce3? massDensity?\n";
      return TCL_ERROR;
  }    

  // get the id and end nodes 
  int eleID, matID;
  int nodes[20];
  double bodyforces[3], massDensity;
  //char *type;
  
  // read the eleTag
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
    opserr << "command: element Brick20N - invalid integer tag " << argv[1+eleArgStart] << endln;
    return TCL_ERROR;
  }
  
  // read the 20 node tags
  int i;
  for (i=0; i<20; i++) {
      if (Tcl_GetInt(interp, argv[2+i+eleArgStart], &nodes[i]) != TCL_OK) {
	opserr << "command: element Brick20N " << eleID << " - invalid integer tag " <<
	  argv[2+i+eleArgStart] << endln;
	  return TCL_ERROR;
      }
  }

  // read in material tag & check the material exists in the model builder
  if (Tcl_GetInt(interp, argv[22+eleArgStart], &matID) != TCL_OK) {
    opserr << "command: element Brick20N " << eleID << " - invalid matID tag " <<      
      argv[22+eleArgStart] << endln;      
    
    return TCL_ERROR;
  }

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
  
  if (theMaterial == 0) {
    opserr << "command: element Brick20N " << eleID << 
      " - no NDMaterial with tag " << argv[22+eleArgStart] << "exists\n";
    return TCL_ERROR;      
  }
  
  //type = argv[11+eleArgStart];

  // read the 3 bodyforce accel's 
  for (i=0; i<3; i++) {
      if (Tcl_GetDouble(interp, argv[23+i+eleArgStart], &bodyforces[i]) != TCL_OK) {
	opserr << "command: element Brick20N " << eleID << " - invalid bodyforces tag " <<
	  argv[23+i+eleArgStart] << endln;
	return TCL_ERROR;
      }
  }

  // now get the massDensity
  if (Tcl_GetDouble(interp, argv[26+eleArgStart], &massDensity) != TCL_OK) {
    opserr << "command: element Brick20N " << eleID << "- invalid massDensity " <<
      argv[26+eleArgStart] << endln;      
      return TCL_ERROR;
  }  
  
  // now create the EightNodeBrick and add it to the Domain
  TwentyNodeBrick *theEle = new TwentyNodeBrick(eleID, nodes[ 0], nodes [1], nodes[ 2], nodes[ 3], nodes[ 4],
                                                      nodes[ 5], nodes [6], nodes[ 7], nodes[ 8], nodes[ 9],
                                                      nodes[10], nodes[11], nodes[12], nodes[13], nodes[14],
                                                      nodes[15], nodes[16], nodes[17], nodes[18], nodes[19],
					       theMaterial, bodyforces[0], bodyforces[1], bodyforces[2], massDensity, 0.0);
					      
  if (theEle == 0) {
    opserr << "command: element Brick20N " << eleID << " - out of memory\n";      
      return TCL_ERROR;
  }

  if (theTclDomain->addElement(theEle) == false) {
    opserr << "command: element Brick20N  - could not add ele: " << eleID << " to domain\n";
    delete theEle;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



