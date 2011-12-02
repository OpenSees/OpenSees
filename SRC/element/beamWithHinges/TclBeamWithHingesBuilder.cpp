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
// $Date: 2000-09-15 08:23:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/beamWithHinges/TclBeamWithHingesBuilder.cpp,v $
                                                                        
                                                                        
// File: ~/tcl/TclElmtBuilder.C
// 
// Written: Remo M. de Souza
// Created: 08/99
// based on TclPlaneFrame.C by fmk and rms
//
// Description: This file contains the implementation of the commands used 
// to add sections and nonlinear frame elements to the model.

#include <stdlib.h>
#include <string.h>
#include <iostream.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <BeamWithHinges2d.h>
#include <BeamWithHinges3d.h>

#include <SectionForceDeformation.h>

#include <CrdTransf2d.h>
#include <CrdTransf3d.h>

#include <TclModelBuilder.h>

static void printCommand(int argc, char **argv)
{
    cerr << "Input command: ";
    for (int i=0; i<argc; i++)
	cerr << argv[i] << " ";
    cerr << endl;
} 

int
TclModelBuilder_addBeamWithHinges (ClientData clientData, Tcl_Interp *interp,
				   int argc, char **argv,
				   Domain *theDomain, TclModelBuilder *theBuilder)		
{
    int NDM = theBuilder->getNDM();
    int NDF = theBuilder->getNDF();
    
    // Plane frame element
    if (NDM == 2 && NDF == 3) {
	if (argc < 13) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: element beamWithHinges tag? ndI? ndJ? secTagI? lenI? secTagJ? lenJ? ";
	    cerr << "E? A? I? transfTag? <-shear shearLength?> <-mass massDens?> <-iter maxIters tolerance>" << endl;
	    return TCL_ERROR;
	}

	int tag, ndI, ndJ, secTagI, secTagJ, transfTag;
	double lenI, lenJ, E, A, I;
	double massDens = 0.0;
	int numIters = 1;
	double tol = 1.0e-10;
	double shearLength = 1.0;
    
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid beamWithHinges tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
	    cerr << "WARNING invalid ndI\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
	    cerr << "WARNING invalid ndJ\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetInt(interp, argv[5], &secTagI) != TCL_OK) {
	    cerr << "WARNING invalid secTagI\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble (interp, argv[6], &lenI) != TCL_OK) {
	    cerr << "WARNING invalid lenI\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[7], &secTagJ) != TCL_OK) {
	    cerr << "WARNING invalid ndJ\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble (interp, argv[8], &lenJ) != TCL_OK) {
	    cerr << "WARNING invalid lenJ\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble (interp, argv[9], &E) != TCL_OK) {
	    cerr << "WARNING invalid E\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble (interp, argv[10], &A) != TCL_OK) {
	    cerr << "WARNING invalid A\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble (interp, argv[11], &I) != TCL_OK) {
	    cerr << "WARNING invalid I\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetInt (interp, argv[12], &transfTag) != TCL_OK) {
	    cerr << "WARNING invalid transfTag\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}	

	if (argc > 13) {
	    for (int i = 13; i < argc; i++) {
		if (strcmp(argv[i],"-mass") == 0 && ++i < argc) {
		    if (Tcl_GetDouble(interp, argv[i], &massDens) != TCL_OK) {
			cerr << "WARNING invalid massDens\n";
			cerr << "BeamWithHinges: " << tag << endl;
			return TCL_ERROR;
		    }
		}
		
		if (strcmp(argv[i],"-shear") == 0 && ++i < argc) {
		    if (Tcl_GetDouble(interp, argv[i], &shearLength) != TCL_OK) {
			cerr << "WARNING invalid shearLength\n";
			cerr << "BeamWithHinges: " << tag << endl;
			return TCL_ERROR;
		    }
		}

		if (strcmp(argv[i],"-iter") == 0 && i+2 < argc) {
			if (Tcl_GetInt(interp, argv[++i], &numIters) != TCL_OK) {
				cerr << "WARNING invalid maxIters\n";
				cerr << "BeamWithHinges: " << tag << endl;
				return TCL_ERROR;
		    }
			if (Tcl_GetDouble(interp, argv[++i], &tol) != TCL_OK) {
				cerr << "WARNING invalid tolerance\n";
				cerr << "BeamWithHinges: " << tag << endl;
				return TCL_ERROR;
		    }
	    }
		}
	}	
	
	// Retrieve section I from the model builder	
	SectionForceDeformation *sectionI = theBuilder->getSection (secTagI);

	if (sectionI == 0) {
	    cerr << "WARNING section I does not exist\n";
	    cerr << "section: " << secTagI; 
	    cerr << "\nBeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	// Retrieve section J from the model builder	
	SectionForceDeformation *sectionJ = theBuilder->getSection (secTagJ);

	if (sectionJ == 0) {
	    cerr << "WARNING section J does not exist\n";
	    cerr << "section: " << secTagJ; 
	    cerr << "\nBeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	CrdTransf2d *theTransf = theBuilder->getCrdTransf2d (transfTag);

	if (theTransf == 0) {
	    cerr << "WARNING geometric transformation does not exist\n";
	    cerr << "geometric transformation: " << transfTag; 
	    cerr << "\nBeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}	
	
	Element *theElement = new BeamWithHinges2d (tag, ndI, ndJ, E, I, A, 1.0, 0.0,
		*sectionI, lenI, *sectionJ, lenJ, *theTransf, shearLength, massDens, numIters, tol);

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0) {
	    cerr << "WARNING ran out of memory creating element\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	if (theDomain->addElement(theElement) == false) {
	    cerr << "WARNING TclElmtBuilder - addBeamWithHinges - could not add element to domain ";
	    cerr << tag << endl;
	    return TCL_ERROR;
	}
    }
    
    else if (NDM == 3 && NDF == 6) {
	if (argc < 16) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: element beamWithHinges tag? ndI? ndJ? secTagI? lenI? secTagJ? lenJ? ";
	    cerr << "E? A? Iz? Iy? G? J? transfTag? <-shear shearLength?> <-mass massDens?> <-iter maxIters tolerance>" << endl;
	    return TCL_ERROR;
	}

	int tag, ndI, ndJ, secTagI, secTagJ, transfTag;
	double lenI, lenJ, E, A, Iz, Iy, G, J;
	double massDens = 0.0;
	int numIters = 1;
	double tol = 1.0e-10;
	double shearLength = 1.0;
    
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid beamWithHinges tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
	    cerr << "WARNING invalid ndI\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
	    cerr << "WARNING invalid ndJ\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetInt(interp, argv[5], &secTagI) != TCL_OK) {
	    cerr << "WARNING invalid secTagI\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble (interp, argv[6], &lenI) != TCL_OK) {
	    cerr << "WARNING invalid lenI\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[7], &secTagJ) != TCL_OK) {
	    cerr << "WARNING invalid ndJ\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble (interp, argv[8], &lenJ) != TCL_OK) {
	    cerr << "WARNING invalid lenJ\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble (interp, argv[9], &E) != TCL_OK) {
	    cerr << "WARNING invalid E\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble (interp, argv[10], &A) != TCL_OK) {
	    cerr << "WARNING invalid A\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetDouble (interp, argv[11], &Iz) != TCL_OK) {
	    cerr << "WARNING invalid Iz\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[12], &Iy) != TCL_OK) {
	    cerr << "WARNING invalid Iy\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble (interp, argv[13], &G) != TCL_OK) {
	    cerr << "WARNING invalid G\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetDouble (interp, argv[14], &J) != TCL_OK) {
	    cerr << "WARNING invalid J\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}	
	
	if (Tcl_GetInt (interp, argv[15], &transfTag) != TCL_OK) {
	    cerr << "WARNING invalid transfTag\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}	

	if (argc > 16) {
	    for (int i = 16; i < argc; i++) {
		if (strcmp(argv[i],"-mass") == 0 && ++i < argc) {
		    if (Tcl_GetDouble(interp, argv[i], &massDens) != TCL_OK) {
			cerr << "WARNING invalid massDens\n";
			cerr << "BeamWithHinges: " << tag << endl;
			return TCL_ERROR;
		    }
		}
		
		if (strcmp(argv[i],"-shear") == 0 && ++i < argc) {
		    if (Tcl_GetDouble(interp, argv[i], &shearLength) != TCL_OK) {
			cerr << "WARNING invalid shearLength\n";
			cerr << "BeamWithHinges: " << tag << endl;
			return TCL_ERROR;
		    }
		}

		if (strcmp(argv[i],"-iter") == 0 && i+2 < argc) {
			if (Tcl_GetInt(interp, argv[++i], &numIters) != TCL_OK) {
				cerr << "WARNING invalid maxIters\n";
				cerr << "BeamWithHinges: " << tag << endl;
				return TCL_ERROR;
		    }
			if (Tcl_GetDouble(interp, argv[++i], &tol) != TCL_OK) {
				cerr << "WARNING invalid tolerance\n";
				cerr << "BeamWithHinges: " << tag << endl;
				return TCL_ERROR;
		    }
	    }
	    }
	}	
	
	// Retrieve section I from the model builder	
	SectionForceDeformation *sectionI = theBuilder->getSection (secTagI);

	if (sectionI == 0) {
	    cerr << "WARNING section I does not exist\n";
	    cerr << "section: " << secTagI; 
	    cerr << "\nBeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	// Retrieve section J from the model builder	
	SectionForceDeformation *sectionJ = theBuilder->getSection (secTagJ);

	if (sectionJ == 0) {
	    cerr << "WARNING section J does not exist\n";
	    cerr << "section: " << secTagJ; 
	    cerr << "\nBeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	CrdTransf3d *theTransf = theBuilder->getCrdTransf3d (transfTag);

	if (theTransf == 0) {
	    cerr << "WARNING geometric transformation does not exist\n";
	    cerr << "geometric transformation: " << transfTag; 
	    cerr << "\nBeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}		

	Element *theElement = new BeamWithHinges3d (tag, ndI, ndJ, E, Iz, Iy, A, G, J, 0.0,
		*sectionI, lenI, *sectionJ, lenJ, *theTransf, shearLength, massDens, numIters, tol);

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0) {
	    cerr << "WARNING ran out of memory creating element\n";
	    cerr << "BeamWithHinges: " << tag << endl;
	    return TCL_ERROR;
	}

	if (theDomain->addElement(theElement) == false) {
	    cerr << "WARNING TclElmtBuilder - addBeamWithHinges - could not add element to domain ";
	    cerr << tag << endl;
	    return TCL_ERROR;
	}
    }
    
    else {
	cerr << "ERROR -- model dimension: " << NDM << " and nodal degrees of freedom: "
	    << NDF << " are incompatible for BeamWithHinges element" << endl;
	return TCL_ERROR;
    }
    
    return TCL_OK;
}
