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
// $Date: 2003-02-25 23:34:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PY/TclPyTzQzMaterialCommand.cpp,v $

//temp out BJ #include <Domain.h>     // RWB bringing in Domain for PyLiq
//temp out BJ // !!!!!!!!!!!!! IFF to be discussed and changed, this Domain, ROss needs
//temp out BJ // that to make PYliq working...

#include <TclModelBuilder.h>

//PY Springs: RWBoulanger and BJeremic
#include <PySimple1.h>// RWB
#include <TzSimple1.h>// RWB
#include <QzSimple1.h>// RWB
// out for the moment#include <PyLiq1.h>// RWB
// out for the moment#include <TzLiq1.h>// RWB

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
TclModelBuilder_addPyTzQzMaterial(ClientData clientData, Tcl_Interp *interp, int argc, 
				  TCL_Char **argv, TclModelBuilder *theTclBuilder)
  //temp out BJ					char **argv, Domain *theDomain, TclModelBuilder *theTclBuilder)
  //temp out BJ // !!!!!!!!!!!!! IFF to be discussed and changed, this theDomain argument, ROss needs
  //temp out BJ // that to make PYliq working...
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

//  INSERTING THE EXTRA LINES FOR PySimple1 //////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
	if (strcmp(argv[1],"PySimple1") == 0) {
	if (argc < 7) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial PySimple1 tag? soilType? pult? y50? drag? dashpot? " << endln;
	    return 0;
	}

	int tag, soilType;
	double pult, y50, drag, dashpot;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial PySimple1 tag" << endln;
	    return 0;		
	}

	if (Tcl_GetInt(interp, argv[3], &soilType) != TCL_OK) {
	    opserr << "WARNING invalid soilType\n";
	    opserr << "uniaxialMaterial PySimple1: " << tag << endln;
	    return 0;	
	}

	if (Tcl_GetDouble(interp, argv[4], &pult) != TCL_OK) {
	    opserr << "WARNING invalid pult\n";
	    opserr << "uniaxialMaterial PySimple1: " << tag << endln;
	    return 0;
	}

	if (Tcl_GetDouble(interp, argv[5], &y50) != TCL_OK) {
	    opserr << "WARNING invalid y50\n";
	    opserr << "uniaxialMaterial PySimple1: " << tag << endln;
	    return 0;	
	}

	if (Tcl_GetDouble(interp, argv[6], &drag) != TCL_OK) {
	    opserr << "WARNING invalid drag\n";
	    opserr << "uniaxialMaterial PySimple1: " << tag << endln;
	    return 0;	
	}

	if (argc == 7) dashpot = 0.0;

	if (argc > 7) {
		if (Tcl_GetDouble(interp, argv[7], &dashpot) != TCL_OK) {
			opserr << "WARNING invalid dashpot\n";
			opserr << "uniaxialMaterial PySimple1: " << tag << endln;
			return 0;	
		}
	}

	// Parsing was successful, allocate the material
	theMaterial = new PySimple1(tag,MAT_TAG_PySimple1,soilType, pult, y50, drag, dashpot);
    }

//  INSERTING THE EXTRA LINES FOR PyLiq1 //////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
//temp out BJ 	else if (strcmp(argv[1],"PyLiq1") == 0) {
//temp out BJ 	if (argc < 10) {
//temp out BJ 	    opserr << "WARNING insufficient arguments\n";
//temp out BJ 	    printCommand(argc,argv);
//temp out BJ 	    opserr << "Want: uniaxialMaterial PyLiq1 tag? soilType? pult? y50? drag? dashpot? solidElem1? solidElem2?" << endln;
//temp out BJ 	    return 0;
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	int tag, soilType, solidElem1, solidElem2;
//temp out BJ 	double pult, y50, drag, dashpot;
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid uniaxialMaterial PyLiq1 tag" << endln;
//temp out BJ 	    return 0;		
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[3], &soilType) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid soilType\n";
//temp out BJ 	    opserr << "uniaxialMaterial PyLiq1: " << tag << endln;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[4], &pult) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid pult\n";
//temp out BJ 	    opserr << "uniaxialMaterial PyLiq1: " << tag << endln;
//temp out BJ 	    return 0;
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[5], &y50) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid y50\n";
//temp out BJ 	    opserr << "uniaxialMaterial PyLiq1: " << tag << endln;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[6], &drag) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid drag\n";
//temp out BJ 	    opserr << "uniaxialMaterial PyLiq1: " << tag << endln;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[7], &dashpot) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid dashpot\n";
//temp out BJ 	    opserr << "uniaxialMaterial PyLiq1: " << tag << endln;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[8], &solidElem1) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid solidElem\n";
//temp out BJ 	    opserr << "uniaxialMaterial PyLiq1: " << tag << endln;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 	
//temp out BJ 	if (Tcl_GetInt(interp, argv[9], &solidElem2) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid solidElem\n";
//temp out BJ 	    opserr << "uniaxialMaterial PyLiq1: " << tag << endln;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	// Parsing was successful, allocate the material
//temp out BJ 	theMaterial = new PyLiq1(tag, MAT_TAG_PyLiq1,soilType, pult, y50, drag, dashpot,
//temp out BJ 							solidElem1, solidElem2, theDomain);
//temp out BJ     }

//  INSERTING THE EXTRA LINES FOR QzSimple1 //////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

	else if (strcmp(argv[1],"QzSimple1") == 0) {
	if (argc < 6) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial QzSimple1 tag? QzType? Qult? z50? suction? dashpot? " << endln;
	    return 0;
	}

	int tag, QzType;
	double Qult, z50, suction, dashpot;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial QzSimple1 tag" << endln;
	    return 0;		
	}

	if (Tcl_GetInt(interp, argv[3], &QzType) != TCL_OK) {
	    opserr << "WARNING invalid QzType\n";
	    opserr << "uniaxialMaterial QzSimple1: " << tag << endln;
	    return 0;	
	}

	if (Tcl_GetDouble(interp, argv[4], &Qult) != TCL_OK) {
	    opserr << "WARNING invalid Qult\n";
	    opserr << "uniaxialMaterial QzSimple1: " << tag << endln;
	    return 0;
	}

	if (Tcl_GetDouble(interp, argv[5], &z50) != TCL_OK) {
	    opserr << "WARNING invalid z50\n";
	    opserr << "uniaxialMaterial QzSimple1: " << tag << endln;
	    return 0;	
	}

	if (argc == 6) {
		suction = 0.0;
		dashpot = 0.0;
	}

	if (argc > 6) {
		if (Tcl_GetDouble(interp, argv[6], &suction) != TCL_OK) {
		    opserr << "WARNING invalid suction\n";
			opserr << "uniaxialMaterial QzSimple1: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[7], &dashpot) != TCL_OK) {
			opserr << "WARNING invalid dashpot\n";
			opserr << "uniaxialMaterial QzSimple1: " << tag << endln;
			return 0;	
		}
	}

	// Parsing was successful, allocate the material
	theMaterial = new QzSimple1(tag, QzType, Qult, z50, suction, dashpot);
    }


//  INSERTING THE EXTRA LINES FOR TzSimple1 //////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
	else if (strcmp(argv[1],"TzSimple1") == 0) {
	if (argc < 6) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial TzSimple1 tag? tzType? tult? z50? dashpot? " << endln;
	    return 0;
	}

	int tag, tzType;
	double tult, z50, dashpot;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial TzSimple1 tag" << endln;
	    return 0;		
	}

	if (Tcl_GetInt(interp, argv[3], &tzType) != TCL_OK) {
	    opserr << "WARNING invalid tzType\n";
	    opserr << "uniaxialMaterial TzSimple1: " << tag << endln;
	    return 0;	
	}

	if (Tcl_GetDouble(interp, argv[4], &tult) != TCL_OK) {
	    opserr << "WARNING invalid tult\n";
	    opserr << "uniaxialMaterial TzSimple1: " << tag << endln;
	    return 0;
	}

	if (Tcl_GetDouble(interp, argv[5], &z50) != TCL_OK) {
	    opserr << "WARNING invalid z50\n";
	    opserr << "uniaxialMaterial TzSimple1: " << tag << endln;
	    return 0;	
	}

	if (argc == 6) dashpot = 0.0;

	if (argc > 6) {
		if (Tcl_GetDouble(interp, argv[6], &dashpot) != TCL_OK) {
			opserr << "WARNING invalid dashpot\n";
			opserr << "uniaxialMaterial TzSimple1: " << tag << endln;
			return 0;	
		}
	}

	// Parsing was successful, allocate the material
	theMaterial = new TzSimple1(tag, MAT_TAG_TzSimple1, tzType, tult, z50, dashpot);
    }

	//  INSERTING THE EXTRA LINES FOR TzLiq1 //////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
//temp out BJ 	else if (strcmp(argv[1],"TzLiq1") == 0) {
//temp out BJ 	if (argc < 9) {
//temp out BJ 	    opserr << "WARNING insufficient arguments\n";
//temp out BJ 	    printCommand(argc,argv);
//temp out BJ 	    opserr << "Want: uniaxialMaterial TzLiq1 tag? tzType? tult? z50? dashpot? solidElem1? solidElem2?" << endln;
//temp out BJ 	    return 0;
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	int tag, tzType, solidElem1, solidElem2;
//temp out BJ 	double tult, z50, dashpot;
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid uniaxialMaterial TzLiq1 tag" << endln;
//temp out BJ 	    return 0;		
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[3], &tzType) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid tzType\n";
//temp out BJ 	    opserr << "uniaxialMaterial TzLiq1: " << tag << endln;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[4], &tult) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid tult\n";
//temp out BJ 	    opserr << "uniaxialMaterial TzLiq1: " << tag << endln;
//temp out BJ 	    return 0;
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[5], &z50) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid z50\n";
//temp out BJ 	    opserr << "uniaxialMaterial TzLiq1: " << tag << endln;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetDouble(interp, argv[6], &dashpot) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid dashpot\n";
//temp out BJ 	    opserr << "uniaxialMaterial TzLiq1: " << tag << endln;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	if (Tcl_GetInt(interp, argv[7], &solidElem1) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid solidElem\n";
//temp out BJ 	    opserr << "uniaxialMaterial TzLiq1: " << tag << endln;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 	
//temp out BJ 	if (Tcl_GetInt(interp, argv[8], &solidElem2) != TCL_OK) {
//temp out BJ 	    opserr << "WARNING invalid solidElem\n";
//temp out BJ 	    opserr << "uniaxialMaterial TzLiq1: " << tag << endln;
//temp out BJ 	    return 0;	
//temp out BJ 	}
//temp out BJ 
//temp out BJ 	// Parsing was successful, allocate the material
//temp out BJ 	theMaterial = new TzLiq1(tag, MAT_TAG_TzLiq1,tzType, tult, z50, dashpot,
//temp out BJ 							solidElem1, solidElem2, theDomain);
//temp out BJ     }
//temp out BJ 
//////////////////////////////////////////////////////////////////////

	return theMaterial;
}
