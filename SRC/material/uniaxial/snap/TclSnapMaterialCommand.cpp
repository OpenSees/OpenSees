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

// $Revision: 1.3 $
// $Date: 2003-02-25 23:34:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/TclSnapMaterialCommand.cpp,v $

// Written: MHS
// Created: Feb 2002
//
// Description: This file contains the implementation of the
// TclModelBuilder_addSnapMaterial() function. 

//#include <SnapBilinearMaterial.h>
//#include <SnapCloughMaterial.h>
//#include <SnapPinchMaterial.h>

#include <Pinching.h>

#include <TclModelBuilder.h>
#include <Vector.h>
#include <string.h>

static void printCommand(int argc, TCL_Char **argv)
{
  opserr << "Input command: ";
  for (int i=0; i<argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
} 

UniaxialMaterial *
TclModelBuilder_addSnapMaterial(ClientData clientData, Tcl_Interp *interp, int argc, 
				TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
  if (argc < 3) {
    opserr << "WARNING insufficient number of arguments\n";
    printCommand(argc, argv);
    return 0;
  }
  
  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    printCommand(argc, argv);
    return 0;
  }
  
  UniaxialMaterial *theMaterial = 0;

  /*
  if (strcmp(argv[1],"Bilinear") == 0) {
    if (argc < 19) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: uniaxialMaterial Bilinear tag? ..." << endln;
      return 0;
    }
    
    Vector input(16);
    double temp;
    
    for (int i = 3, j = 0; j < 16; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
	opserr << "WARNING invalid input, data " << i << endln;
	printCommand(argc, argv);
	return 0;
      }
      input(j) = temp;
    }
    
    theMaterial = new SnapBilinearMaterial(tag, input);
  }

  if (strcmp(argv[1],"Clough") == 0) {
    if (argc < 19) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: uniaxialMaterial Clough tag? ..." << endln;
      return 0;
    }
    
    Vector input(16);
    double temp;
    
    for (int i = 3, j = 0; j < 16; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
	opserr << "WARNING invalid input, data " << i << endln;
	printCommand(argc, argv);
	return 0;
      }
      input(j) = temp;
    }
    
    theMaterial = new SnapCloughMaterial(tag, input);
  }

  if (strcmp(argv[1],"Pinch") == 0) {
    if (argc < 22) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: uniaxialMaterial Pinch tag? ..." << endln;
      return 0;
    }
    
    Vector input(19);
    double temp;
    
    for (int i = 3, j = 0; j < 19; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
	opserr << "WARNING invalid input, data " << i << endln;
	printCommand(argc, argv);
	return 0;
      }
      input(j) = temp;
    }
    
    theMaterial = new SnapPinchMaterial(tag, input);
  }
  */

  if (strcmp(argv[1],"Pinching") == 0) {
    if (argc < 22) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: uniaxialMaterial Pinching tag? ..." << endln;
      return 0;
    }
    
    Vector input(19);
    double temp;
    
    for (int i = 3, j = 0; j < 19; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
	opserr << "WARNING invalid input, data " << i << endln;
	printCommand(argc, argv);
	return 0;
      }
      input(j) = temp;
    }
    
    theMaterial = new Pinching(tag, input);
  }
  
  return theMaterial;
}
