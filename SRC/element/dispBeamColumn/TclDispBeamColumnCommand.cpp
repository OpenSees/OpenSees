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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-25 23:32:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumn/TclDispBeamColumnCommand.cpp,v $
                                                                        
// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the implementation of the 
// TclModelBuilder_addDispBeamColumn() command. 

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <DispBeamColumn2d.h>
#include <DispBeamColumn3d.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addDispBeamColumn(ClientData clientData, Tcl_Interp *interp,  
				int argc, 
				TCL_Char **argv, 
				Domain*theTclDomain,
				TclModelBuilder *theTclBuilder)
{
	// ensure the destructor has not been called - 
	if (theTclBuilder == 0) {
		opserr << "WARNING builder has been destroyed\n";    
		return TCL_ERROR;
	}

	int ndm = theTclBuilder->getNDM();
	int ndf = theTclBuilder->getNDF();

	int ok = 0;
	if (ndm == 2 && ndf == 3)
		ok = 1;
	if (ndm == 3 && ndf == 6)
		ok = 1;

	if (ok == 0) {
		opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
			<< " not compatible with dispBeamColumn element" << endln;
		return TCL_ERROR;
	}

	if (argc < 8) {
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
		return TCL_ERROR;
	}

	// get the id and end nodes 
	int eleTag, iNode, jNode, nIP, transfTag;
	int secTag[10]; // Max size of integration rule ... can change if needed
	int argi = 2;

	if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK) {
		opserr << "WARNING invalid dispBeamColumn eleTag" << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK) {
		opserr << "WARNING invalid iNode ";
		opserr << "dispBeamColumn element: " << eleTag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK) {
		opserr << "WARNING invalid jNode ";
		opserr << "dispBeamColumn element: " << eleTag << endln;
		return TCL_ERROR;
	}
  
	if (Tcl_GetInt(interp, argv[argi++], &nIP) != TCL_OK) {
		opserr << "WARNING invalid nIP ";
		opserr << "dispBeamColumn element: " << eleTag << endln;
		return TCL_ERROR;
	}  
  
	if (strcmp(argv[argi], "-sections") == 0) {
	  argi++;
	  if (argi+nIP > argc) {
	    opserr << "WARNING insufficient number of section tags - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	    return TCL_ERROR;
	  }
	  int section;
	  for (int i = 0; i < nIP; i++) {
	    if (Tcl_GetInt(interp, argv[argi+i], &section) != TCL_OK) {
	      opserr << "WARNING invalid secTag - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	      return TCL_ERROR;
	    }
	    secTag[i] = section;
	  }
	  argi += nIP;
	}
	
	else {
	  int section;
	  if (Tcl_GetInt(interp, argv[argi++], &section) != TCL_OK) {
	    opserr << "WARNING invalid secTag - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	    return TCL_ERROR;
	  }
	  for (int i = 0; i < nIP; i++)
	    secTag[i] = section;
	}
	
	if (argi >= argc || Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK) {
	  opserr << "WARNING invalid transfTag? - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	  return TCL_ERROR;
	}

	double massDens = 0.0;

	while (argi != argc) {
	  if (strcmp(argv[argi++],"-mass") == 0 && argi < argc) {
	    if (Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK) {
	      opserr << "WARNING invalid massDens - element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
	      return TCL_ERROR;
	    }
	  }
	}

	SectionForceDeformation **sections = new SectionForceDeformation* [nIP];
	
	if (!sections) {
	  opserr << "WARNING TclElmtBuilder - addFrameElement - Insufficient memory to create sections\n";
	  return TCL_ERROR;
	}
	
	for (int j=0; j<nIP; j++) {
	  SectionForceDeformation *theSection = theTclBuilder->getSection(secTag[j]);
	  
	  if (theSection == 0) {
	    opserr << "WARNING TclElmtBuilder - frameElement - no Section found with tag ";
	    opserr << secTag[j] << endln;
	    delete [] sections;
	    return TCL_ERROR;
	  }

	  sections[j] = theSection;
	}
	
	Element *theElement = 0;

	if (ndm == 2) {

		CrdTransf2d *theTransf = theTclBuilder->getCrdTransf2d(transfTag);
      
		if (theTransf == 0) {
			opserr << "WARNING transformation not found\n";
			opserr << "transformation: " << transfTag;
			opserr << "\ndispBeamColumn element: " << eleTag << endln;
			return TCL_ERROR;
		}

		// now create the DispBeamColumn and add it to the Domain
		theElement = new DispBeamColumn2d(eleTag,iNode,jNode,nIP,sections,*theTransf,massDens);

		delete [] sections;
	}

	if (ndm == 3) {

		CrdTransf3d *theTransf = theTclBuilder->getCrdTransf3d(transfTag);
      
		if (theTransf == 0) {
			opserr << "WARNING transformation not found\n";
			opserr << "transformation: " << transfTag;
			opserr << "\ndispBeamColumn element: " << eleTag << endln;
			return TCL_ERROR;
		}

		// now create the DispBeamColumn and add it to the Domain
		theElement = new DispBeamColumn3d(eleTag,iNode,jNode,nIP,sections,*theTransf,massDens);

		delete [] sections;
	}

	if (theElement == 0) {
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "dispBeamColumn element: " << eleTag << endln;
		return TCL_ERROR;
	}

	if (theTclDomain->addElement(theElement) == false) {
		opserr << "WARNING could not add element to the domain\n";
		opserr << "dispBeamColumn element: " << eleTag << endln;
		delete theElement;
		return TCL_ERROR;
	}

	// if get here we have sucessfully created the element and added it to the domain
	return TCL_OK;
}



