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

// !!! get rid of ####
#include <TclModelBuilder.h>
#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>
#include "FourNodeSFI_MVLEM3D.h"
#include "elementAPI.h" // CSUF

extern void printCommand(int argc, TCL_Char **argv);

int TclModelBuilder_addFourNodeSFI_MVLEM3D(ClientData clientData,
	Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
	TclModelBuilder *theTclBuilder, int eleArgStart)
{

	if (theTclBuilder == 0) {
		opserr << "WARNING builder has been destroyed - FourNodeSFI_MVLEM3D\n";
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
			<< " not compatible with FourNodeSFI_MVLEM3D element" << endln;
		return TCL_ERROR;
	}

	// FourNodeSFI_MVLEM3D eleTag iNode jNode kNode lNode m c nu Tfactor-thick fiberThick -width fiberWidth -mat matTags
	if ((argc - eleArgStart) < 9) { // 
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: FourNodeSFI_MVLEM3D eleTag iNode jNode kNode lNode m c nu Tfactor -thick fiberThick -width fiberWidth -mat matTags\n"; 	// !!! This doesnt see to be up to date @@@ DONE
		return TCL_ERROR;
	}

	// !!! create comments to make it clear what is going on with node ordering
	// !!! in the quad wall there is a number of checks to make sure that the user assigned the nodes correclyt. 
	// !!! I think you should implement them. Maybe they will be redundant based ont he alter checks. You should think it through so you dont miss anything and dont double count anything.
	// Order Nodes Appropriately for Element Interpolation ......................
	// Nodes can be assigned in any order. Model orders nodes for interpolation as shown below:
	//    4........3 
	//    .        .
	//    .        .
	//    .        .
	//    1........2

	// !!! screen messages below are not compatible with variable names - check and make consistent @@@ DONE
	// get the id and end nodes 
	int tag, iNode, jNode, kNode, lNode, numMat, numThick, numWidth, matTag, i, m;
	double c, thick, width;
	double NUelastic, Tfactor; //##### new elastic parameters

	int argi = eleArgStart + 1;
	if (Tcl_GetInt(interp, argv[argi], &tag) != TCL_OK) {
		opserr << "WARNING invalid FourNodeSFI_MVLEM3D eleTag\n";
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetInt(interp, argv[argi], &iNode) != TCL_OK) {
		opserr << "WARNING invalid iNode\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetInt(interp, argv[argi], &jNode) != TCL_OK) {
		opserr << "WARNING invalid jNode\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetInt(interp, argv[argi], &lNode) != TCL_OK) {
		opserr << "WARNING invalid lNode\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetInt(interp, argv[argi], &kNode) != TCL_OK) {
		opserr << "WARNING invalid kNode\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetInt(interp, argv[argi], &m) != TCL_OK) {
		opserr << "WARNING invalid m\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetDouble(interp, argv[argi], &c) != TCL_OK) {
		opserr << "WARNING invalid c\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetDouble(interp, argv[argi], &NUelastic) != TCL_OK) {
		opserr << "WARNING invalid nu\n";
		opserr << "SFI_MVLEM_3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	if (Tcl_GetDouble(interp, argv[argi], &Tfactor) != TCL_OK) {
		opserr << "WARNING invalid Tfactor\n";
		opserr << "SFI_MVLEM_3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	// read the number of cross-section thickness ............................................................... 
	numThick = 0;
	if (strcmp(argv[argi], "-thick") != 0) {
		opserr << "WARNING expecting -thick value\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	i = argi;
	while (i < argc - 2 * (m + 1)) { // why is this formula not the same as in MVLEM?
		numThick++;
		i++;
	}

	if (numThick != m) {
		opserr << "WARNING :WRONG THICKNESS number\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}

	// create array of thickenss
	double *theThickness = new double[numThick];

	for (i = 0; i<numThick; i++) {
		theThickness[i] = 0;

		if (Tcl_GetDouble(interp, argv[argi], &thick) != TCL_OK) {
			opserr << "WARNING invalid thickness\n";
			opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
			return TCL_ERROR;
		}

		theThickness[i] = thick;

		argi++;

	}

	// read the number of cross-section widths ............................................................... 
	numWidth = 0;
	if (strcmp(argv[argi], "-width") != 0) {
		opserr << "WARNING expecting -width flag\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	i = argi;
	while (i < argc - 1 * (m + 1)) { // why is this formula not the same as in MVLEM?
		numWidth++;
		i++;
	}

	if (numWidth != m) {
		opserr << "WARNING :WRONG WIDTH number\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}

	// create array of widths 
	double *theWidth = new double[numWidth];

	for (i = 0; i<numWidth; i++) {
		theWidth[i] = 0;

		if (Tcl_GetDouble(interp, argv[argi], &width) != TCL_OK) {
			opserr << "WARNING invalid width\n";
			opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
			return TCL_ERROR;
		}

		theWidth[i] = width;

		argi++;
	}

	// read the number of macro-fiber materials ................................
	numMat = 0;
	if (strcmp(argv[argi], "-mat") != 0) {
		opserr << "WARNING expecting -mat flag\n";
		opserr << "FourNodeSFI_MVLEM3D_4 element: " << tag << endln;
		return TCL_ERROR;
	}
	argi++;

	i = argi;
	while (i < argc) {
		numMat++;
		i++;
	}

	if (numMat != m) {
		opserr << "WARNING :WRONG MAT number\n";
		opserr << "FourNodeSFI_MVLEM3D_4 element: " << tag << endln;
		return TCL_ERROR;
	}

	// create array of NDMaterials
	NDMaterial **theMaterials = new NDMaterial*[numMat];

	for (i = 0; i<numMat; i++) {
		theMaterials[i] = 0;
		if (Tcl_GetInt(interp, argv[argi], &matTag) != TCL_OK) {
			opserr << "WARNING invalid matTag\n";
			opserr << "FourNodeSFI_MVLEM3D_4 element: " << tag << endln;
			return TCL_ERROR;
		}

		theMaterials[i] = OPS_GetNDMaterial(matTag); // CSUF
		if (theMaterials[i] == 0) {
			opserr << "WARNING material model not found\n";
			opserr << "NDMaterial " << matTag << endln;
			opserr << "FourNodeSFI_MVLEM3D_4 element: " << tag << endln;
			return TCL_ERROR;
		}
		argi++;
	}
	//.............................................................................

	// get pointers to end nodes - used to obtain their locations in the Element constructor needed for internal nodes
	Node * theNodei = theTclDomain->getNode(iNode);
	Node * theNodej = theTclDomain->getNode(jNode);
	Node * theNodek = theTclDomain->getNode(kNode);
	Node * theNodel = theTclDomain->getNode(lNode);

	theElement = new FourNodeSFI_MVLEM3D(tag, iNode, jNode, kNode, lNode, theNodei, theNodej, theNodek, theNodel, theMaterials, theThickness, theWidth, theTclDomain, m, c, NUelastic, Tfactor);

	// cleanup dynamic memory
	if (theMaterials != 0)
		delete[] theMaterials;
	if (theElement == 0) {
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		return TCL_ERROR;
	}

	// then add the FourNodeSFI_MVLEM3D to the domain
	if (theTclDomain->addElement(theElement) == false) {
		opserr << "WARNING could not add element to the domain\n";
		opserr << "FourNodeSFI_MVLEM3D element: " << tag << endln;
		delete theElement;
		return TCL_ERROR;
	}

	return TCL_OK;
}