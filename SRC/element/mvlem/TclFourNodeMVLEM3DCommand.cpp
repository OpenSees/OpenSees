// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								Kamiar Kalbasi
//								California State University, Fullerton 
//
// Created:
//
// Description: 
//
// References:
//
// Source: /usr/local/cvs/OpenSees/SRC/element/FourNodeMVLEM3D/FourNodeMVLEM3D.h
//
// Rev: 1

#include <TclModelBuilder.h>
#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>
#include <UniaxialMaterial.h>
#include "FourNodeMVLEM3D.h"

extern void printCommand(int argc, TCL_Char **argv);

int TclModelBuilder_addFourNodeMVLEM3D(ClientData clientData,
	Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
	TclModelBuilder *theTclBuilder, int eleArgStart)
{

	if (theTclBuilder == 0) {
		opserr << "WARNING builder has been destroyed - FourNodeMVLEM3D\n";
		return TCL_ERROR;
	}

	Element *theElement = 0;

	int ndm = theTclBuilder->getNDM();
	int ndf = theTclBuilder->getNDF();

	int ok = 0;
	if (ndm == 3 && ndf == 6)
		ok = 1;

	if (ok == 0) {
		opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
			<< " not compatible with FourNodeMVLEM3D element" << endln;
		return TCL_ERROR;
	}

	// !!! This doesnt see to be up to date
	// Check the number of arguments is correct
	if ((argc - eleArgStart) < 11) {
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: FourNodeMVLEM3D eleTag Dens iNode jNode m c nu Tfactor -thick {fiberThick} -width {fiberWidth} -rho {Rho} -matConcrete {matTagsConcrete} -matSteel {matTagsSteel} -matShear {matTagShear}\n";
		return TCL_ERROR;

	}      // get the id and end nodes 

	int tag, iNode, jNode, kNode, lNode, numMatC, numMatS, numMatSh, numThick, numWidth, numRho, matTag, i, m;
	double c, Dens, thick, width, Rho;
	double NUelastic, Tfactor;

	int argi = eleArgStart + 1;
	if (Tcl_GetInt(interp, argv[argi], &tag) != TCL_OK) {
		opserr << "WARNING invalid FourNodeMVLEM3D eleTag\n";
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetDouble(interp, argv[argi], &Dens) != TCL_OK) {
		opserr << "WARNING invalid Dens\n";
		opserr << "SFI_FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	// Order Nodes Appropriately for Element Interpolation ......................
	// Nodes should be assigned in counterclockwise. Model orders nodes as shown below:
	//    4........3 
	//    .        .
	//    .        .
	//    .        .
	//    1........2

	if (Tcl_GetInt(interp, argv[argi], &iNode) != TCL_OK) {
		opserr << "WARNING invalid iNode\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetInt(interp, argv[argi], &jNode) != TCL_OK) {
		opserr << "WARNING invalid jNode\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetInt(interp, argv[argi], &lNode) != TCL_OK) {
		opserr << "WARNING invalid lNode\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetInt(interp, argv[argi], &kNode) != TCL_OK) {
		opserr << "WARNING invalid kNode\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetInt(interp, argv[argi], &m) != TCL_OK) {
		opserr << "WARNING invalid m\n";
		opserr << "SFI_FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetDouble(interp, argv[argi], &c) != TCL_OK) {
		opserr << "WARNING invalid c\n";
		opserr << "SFI_FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetDouble(interp, argv[argi], &NUelastic) != TCL_OK) {
		opserr << "WARNING invalid NUelastic\n";
		opserr << "SFI_FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetDouble(interp, argv[argi], &Tfactor) != TCL_OK) {
		opserr << "WARNING invalid Tfactor\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	// Read fiber thicknesses ............................................................... 
	numThick = 0;
	if (strcmp(argv[argi], "-thick") != 0) {
		opserr << "WARNING expecting -thick value\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	i = argi;
	while (i < argc - 4 * (m + 1) - 2) {
		numThick++;
		i++;
	}

	if (numThick != m) {
		opserr << "WARNING :WRONG THICKNESS number\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}

	// create array of m fiber thickenss
	double *theThickness = new double[numThick];

	for (i = 0; i<numThick; i++) {
		theThickness[i] = 0;

		if (Tcl_GetDouble(interp, argv[argi], &thick) != TCL_OK) {
			opserr << "WARNING invalid thickness\n";
			opserr << "FourNodeMVLEM3D element: " << tag << endln;
			return TCL_ERROR;
		}

		theThickness[i] = thick;

		argi++;

	}

	// Read fiber widths ............................................................... 
	numWidth = 0;
	if (strcmp(argv[argi], "-width") != 0) {
		opserr << "WARNING expecting -width flag\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	i = argi;
	while (i < argc - 3 * (m + 1) - 2) {
		numWidth++;
		i++;
	}

	if (numWidth != m) {
		opserr << "WARNING :WRONG WIDTH number\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}

	// Create array of m fiber widths 
	double *theWidth = new double[numWidth];

	for (i = 0; i<numWidth; i++) {
		theWidth[i] = 0;

		if (Tcl_GetDouble(interp, argv[argi], &width) != TCL_OK) {
			opserr << "WARNING invalid width\n";
			opserr << "FourNodeMVLEM3D element: " << tag << endln;
			return TCL_ERROR;
		}

		theWidth[i] = width;

		argi++;
	}

	// Read fiber reinforcing ratios ............................................................... 
	numRho = 0;
	if (strcmp(argv[argi], "-rho") != 0) {
		opserr << "WARNING expecting -rho flag\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	i = argi;
	while (i < argc - 2 * (m + 1) - 2) {
		numRho++;
		i++;
	}

	if (numRho != m) {
		opserr << "WARNING :WRONG RHO number\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}

	// Create array of m fiber reinforcing ratios
	double *theRho = new double[numRho];

	for (i = 0; i<numRho; i++) {
		theRho[i] = 0;

		if (Tcl_GetDouble(interp, argv[argi], &Rho) != TCL_OK) {
			opserr << "WARNING invalid width\n";
			opserr << "FourNodeMVLEM3D element: " << tag << endln;
			return TCL_ERROR;
		}

		theRho[i] = Rho;

		argi++;
	}

	// Read fiber Concrete materials ................................
	numMatC = 0;
	if (strcmp(argv[argi], "-matConcrete") != 0) {
		opserr << "WARNING expecting -matConcrete flag\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	i = argi;
	while (i < argc - 1 * (m + 1) - 2) {
		numMatC++;
		i++;
	}

	if (numMatC != m) {
		opserr << "WARNING :WRONG Concrete material number\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;

	}

	// Create array of m uniaxial fiber Concrete materials
	UniaxialMaterial **theMaterialsConcrete = new UniaxialMaterial*[numMatC];

	for (i = 0; i<numMatC; i++) {
		theMaterialsConcrete[i] = 0;
		if (Tcl_GetInt(interp, argv[argi], &matTag) != TCL_OK) {
			opserr << "WARNING invalid matTag\n";
			opserr << "FourNodeMVLEM3D element: " << tag << endln;
			return TCL_ERROR;
		}

		theMaterialsConcrete[i] = OPS_getUniaxialMaterial(matTag);

		if (theMaterialsConcrete[i] == 0) {
			opserr << "WARNING Concrete material model not found\n";
			opserr << "uniaxialMaterial " << matTag << endln;
			opserr << "FourNodeMVLEM3D element: " << tag << endln;
			return TCL_ERROR;
		}
		argi++;
	}

	// Read fiber Steel materials ................................
	numMatS = 0;
	if (strcmp(argv[argi], "-matSteel") != 0) {
		opserr << "WARNING expecting -matSteel flag\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	i = argi;
	while (i < argc - 0 * (m + 1) - 2) {
		numMatS++;
		i++;
	}

	if (numMatS != m) {
		opserr << "WARNING :WRONG Steel material number\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}

	// Create array of m uniaxial fiber Steel materials
	UniaxialMaterial **theMaterialsSteel = new UniaxialMaterial*[numMatS];

	for (i = 0; i<numMatS; i++) {
		theMaterialsSteel[i] = 0;
		if (Tcl_GetInt(interp, argv[argi], &matTag) != TCL_OK) {
			opserr << "WARNING invalid matTag\n";
			opserr << "FourNodeMVLEM3D element: " << tag << endln;
			return TCL_ERROR;
		}

		theMaterialsSteel[i] = OPS_getUniaxialMaterial(matTag);

		if (theMaterialsSteel[i] == 0) {
			opserr << "WARNING Steel material model not found\n";
			opserr << "uniaxialMaterial " << matTag << endln;
			opserr << "FourNodeMVLEM3D element: " << tag << endln;
			return TCL_ERROR;
		}
		argi++;
	}

	// Read Shear material ................................
	numMatSh = 0;
	if (strcmp(argv[argi], "-matShear") != 0) {
		opserr << "WARNING expecting -matShear flag\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	i = argi;
	while (i < argc) {
		numMatSh++;
		i++;
	}

	if (numMatSh != 1) {
		opserr << "WARNING :WRONG Shear material number\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}

	// Create array of 1 uniaxial Shear material
	UniaxialMaterial **theMaterialsShear = new UniaxialMaterial*[numMatSh];

	for (i = 0; i<numMatSh; i++) {
		theMaterialsShear[i] = 0;
		if (Tcl_GetInt(interp, argv[argi], &matTag) != TCL_OK) {
			opserr << "WARNING invalid matTag\n";
			opserr << "FourNodeMVLEM3D element: " << tag << endln;
			return TCL_ERROR;
		}

		theMaterialsShear[i] = OPS_getUniaxialMaterial(matTag);

		if (theMaterialsShear[i] == 0) {
			opserr << "WARNING Shear material model not found\n";
			opserr << "uniaxialMaterial " << matTag << endln;
			opserr << "FourNodeMVLEM3D element: " << tag << endln;
			return TCL_ERROR;
		}
		argi++;
	}

	//.............................................................................

	// Get pointers to end nodes - passed to obtain element height
	Node * theNodei = theTclDomain->getNode(iNode);
	Node * theNodej = theTclDomain->getNode(jNode);
	Node * theNodek = theTclDomain->getNode(kNode);
	Node * theNodel = theTclDomain->getNode(lNode);

	// Create the FourNodeMVLEM3D element
	theElement = new FourNodeMVLEM3D(tag, Dens, iNode, jNode, kNode, lNode, theNodei, theNodej, theNodek, theNodel, theMaterialsConcrete, theMaterialsSteel, theMaterialsShear, theRho, theThickness, theWidth, m, c, NUelastic, Tfactor);

	// Cleanup dynamic memory
	if (theMaterialsConcrete != 0)
		delete[] theMaterialsConcrete;
	if (theMaterialsSteel != 0)
		delete[] theMaterialsSteel;
	if (theMaterialsShear != 0)
		delete[] theMaterialsShear;
	if (theElement == 0) {
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		return TCL_ERROR;

	}  	// then add the FourNodeMVLEM3D to the domain
	if (theTclDomain->addElement(theElement) == false) {
		opserr << "WARNING could not add element to the domain\n";
		opserr << "FourNodeMVLEM3D element: " << tag << endln;
		delete theElement;
		return TCL_ERROR;
	}
	return TCL_OK;
}