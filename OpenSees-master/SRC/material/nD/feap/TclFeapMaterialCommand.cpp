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
// $Date: 2003-02-25 23:33:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/feap/TclFeapMaterialCommand.cpp,v $

#include <FeapMaterial01.h>
#include <FeapMaterial02.h>
#include <FeapMaterial03.h>

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

NDMaterial *
TclModelBuilder_addFeapMaterial(ClientData clientData, Tcl_Interp *interp,
				int argc, TCL_Char **argv,
				TclModelBuilder *theTclBuilder)
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
  
  NDMaterial *theMaterial = 0;

  if (strcmp(argv[1],"ElasticFeap") == 0 || strcmp(argv[1],"FeapElastic") == 0) {
    if (argc < 5) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc,argv);
      opserr << "Want: nDMaterial ElasticFeap tag? E? nu?" << endln;
      return 0;
    }
    
    double E, nu;

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      printCommand(argc, argv);
      return 0;       
    }
    if (Tcl_GetDouble(interp, argv[4], &nu) != TCL_OK) {
      opserr << "WARNING invalid nu\n";
      printCommand(argc, argv);
      return 0;       
    }

    theMaterial = new FeapMaterial01(tag, E, nu);
  }

  else if (strcmp(argv[1],"J2Feap") == 0 || strcmp(argv[1],"FeapJ2") == 0) {
    if (argc < 7) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc,argv);
      opserr << "Want: nDMaterial J2Feap tag? K? G? sigY? Hiso?" << endln;
      return 0;
    }

    double K, G, sigY, Hiso;

    if (Tcl_GetDouble(interp, argv[3], &K) != TCL_OK) {
      opserr << "WARNING invalid K\n";
      printCommand(argc, argv);
      return 0;       
    }
    if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      printCommand(argc, argv);
      return 0;       
    }
    if (Tcl_GetDouble(interp, argv[5], &sigY) != TCL_OK) {
      opserr << "WARNING invalid sigY\n";
      printCommand(argc, argv);
      return 0;       
    }
    if (Tcl_GetDouble(interp, argv[6], &Hiso) != TCL_OK) {
      opserr << "WARNING invalid Hiso\n";
      printCommand(argc, argv);
      return 0;       
    }

    theMaterial = new FeapMaterial03(tag, K, G, sigY, Hiso);
  }

  else if (strcmp(argv[1],"ViscousFeap") == 0 || strcmp(argv[1],"FeapViscous") == 0) {
    if (argc < 10) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc,argv);
      opserr << "Want: nDMaterial ViscousFeap tag? K? G? muK? muG? lamK? lamG? theta?" << endln;
      return 0;
    }

    double K, G, muK, muG, lamK, lamG, theta;

    if (Tcl_GetDouble(interp, argv[3], &K) != TCL_OK) {
      opserr << "WARNING invalid K\n";
      printCommand(argc, argv);
      return 0;       
    }
    if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      printCommand(argc, argv);
      return 0;
    }
    if (Tcl_GetDouble(interp, argv[5], &muK) != TCL_OK) {
      opserr << "WARNING invalid muK\n";
      printCommand(argc, argv);
      return 0;       
    }
    if (Tcl_GetDouble(interp, argv[6], &muG) != TCL_OK) {
      opserr << "WARNING invalid muG\n";
      printCommand(argc, argv);
      return 0;
    }
    if (Tcl_GetDouble(interp, argv[7], &lamK) != TCL_OK) {
      opserr << "WARNING invalid lamK\n";
      printCommand(argc, argv);
      return 0;       
    }
    if (Tcl_GetDouble(interp, argv[8], &lamG) != TCL_OK) {
      opserr << "WARNING invalid lamG\n";
      printCommand(argc, argv);
      return 0;
    }
    if (Tcl_GetDouble(interp, argv[9], &theta) != TCL_OK) {
      opserr << "WARNING invalid theta\n";
      printCommand(argc, argv);
      return 0;
    }

    theMaterial = new FeapMaterial02(tag, K, G, muK, muG, lamK, lamG, theta);
  }
  
  return theMaterial;
}
