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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp,v $
                                                                        
                                                                        
// File: ~/material/nD/TclModelBuilderNDMaterialComand.C
// 
// Written: MHS 
// Created: June 2000
// Revision: A
//
// Description: This file contains the function invoked when the user invokes
// the nDMaterial command in the interpreter.
//
// What: "@(#) TclModelBuilderNDMaterialCommand.C, revA"

#include <TclModelBuilder.h>

#include <ElasticIsotropicMaterial.h>

#include <string.h>

static void printCommand(int argc, char **argv)
{
    cerr << "Input command: ";
    for (int i=0; i<argc; i++)
	cerr << argv[i] << " ";
    cerr << endl;
} 

int
TclModelBuilderNDMaterialCommand (ClientData clientData, Tcl_Interp *interp, int argc,
				  char **argv, TclModelBuilder *theTclBuilder)
{
    // Make sure there is a minimum number of arguments
    if (argc < 3) {
	cerr << "WARNING insufficient number of ND material arguments\n";
	cerr << "Want: nDMaterial type? tag? <specific material args>" << endl;
	return TCL_ERROR;
    }

    // Pointer to an ND material that will be added to the model builder
    NDMaterial *theMaterial = 0;

    // Check argv[1] for ND material type
    if (strcmp(argv[1],"ElasticIsotropic") == 0) {
	if (argc < 5) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: nDMaterial ElasticIsotropic tag? E? v? <type?>" << endl;
	    return TCL_ERROR;
	}    

	int tag;
	double E, v;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid ElasticIsotropic tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    cerr << "WARNING invalid E\n";
	    cerr << "nDMaterial ElasticIsotropic: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[4], &v) != TCL_OK) {
	    cerr << "WARNING invalid v\n";
	    cerr << "nDMaterial ElasticIsotropic: " << tag << endl;
	    return TCL_ERROR;	
	}

	NDMaterial *temp = new ElasticIsotropicMaterial (tag, E, v);
	
	// Obtain a specific copy if requested
	if (argc > 5)
	    theMaterial = temp->getCopy(argv[5]);
	else
	    theMaterial = temp;
	}	
    
    else {
	cerr << "WARNING unknown type of nDMaterial: " << argv[1];
	cerr << "\nValid types: ElasticIsotropic\n";
	return TCL_ERROR;
    }

    // Ensure we have created the Material, out of memory if got here and no material
    if (theMaterial == 0) {
	cerr << "WARNING ran out of memory creating nDMaterial\n";
	cerr << argv[1] << endl;
	return TCL_ERROR;
    }

    // Now add the material to the modelBuilder
    if (theTclBuilder->addNDMaterial(*theMaterial) < 0) {
	cerr << "WARNING could not add material to the domain\n";
	cerr << *theMaterial << endl;
	delete theMaterial; // invoke the material objects destructor, otherwise mem leak
	return TCL_ERROR;
    }
    
    return TCL_OK;
}
