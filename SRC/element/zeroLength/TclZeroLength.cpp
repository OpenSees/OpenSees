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
// $Date: 2003-02-25 23:33:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/TclZeroLength.cpp,v $
                                                                        
                                                                        
// File: ~/element/zeroLength/TclZeroLength.C
// 
// Written: fmk
// Created: 01/00
//
// Description: This file contains the implementation of the command used 
// to add zero length elements to the model.


#include <stdlib.h>
#include <string.h>
#include <tcl.h>

#include <ZeroLength.h>
//#include <ZeroLengthND.h>
#include <ZeroLengthSection.h>
#include <TclModelBuilder.h>
#include <ID.h>
#include <Vector.h>
#include <Domain.h>


int
TclModelBuilder_addZeroLength(ClientData clientData, Tcl_Interp *interp,
			      int argc, TCL_Char **argv,
			      Domain *theDomain,
			      TclModelBuilder *theBuilder) {

    int ndm = theBuilder->getNDM(); // the spatial dimension of the problem

    //
    // first scan the command line to obtain eleID, iNode, jNode, material ID's
    // and their directions, and the orientation of ele xPrime and yPrime not
    // along the global x and y axis
    //
    
    int eleTag, iNode, jNode;
    
    // a quick check on number of args
    if (argc < 9) {
      opserr << "WARNING too few arguments " <<
	"want - element ZeroLength eleTag? iNode? jNode? " <<
	"-mat matID1? ... -dir dirMat1? .. " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";

	return TCL_ERROR;
    }

    // get the ele tag 
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
      opserr << "WARNING invalied eleTag " << argv[2] <<
	"- element ZeroLength eleTag? iNode? jNode? -mat matID1? ... -dir dirMat1? .. " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";
	return TCL_ERROR;
    }

    // get the two end nodes
    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
      opserr << "WARNING invalied iNode " << argv[3] << 
	"- element ZeroLength eleTag? iNode? jNode? " <<
	"-mat matID1? ... -dir dirMat1? .. " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";

	return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode " << argv[4] <<
	"- element ZeroLength eleTag? iNode? jNode? " <<
	"-mat matID1? ... -dir dirMat1? .. " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";
	return TCL_ERROR;
    }

    // create an array of material pointers, to do this first count
    // the materials to create the array then get matID's and from ModelBuilder
    // obtain pointers to the material objects

    // read the number of materials
    int numMat = 0;
    if (strcmp(argv[5],"-mat") != 0) {
      opserr << "WARNING expecting -mat flag %s %s %s %s\n" << argv[5] <<
	"- element ZeroLength eleTag? iNode? jNode? " <<
	"-mat matID1? ... -dir dirMat1? .. " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";

	return TCL_ERROR;	
    }

    int argi = 6; 
    while ((argi < argc) && (strcmp(argv[argi],"-dir") != 0)) {
	numMat++;
	argi++;
    }

    if (argi == argc) { // check we encounterd the -dirn flag
      opserr << "WARNING no -dirn flag encountered " <<
	"- element ZeroLength eleTag? iNode? jNode? " <<
	"-mat matID1? ... -dir dirMat1? .. " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";
	return TCL_ERROR;
    }	

    if (numMat == 0) {
      opserr << "WARNING no materials specified " << 
	"- element ZeroLength eleTag? iNode? jNode? " <<
	"-mat <matID1? ... -dir irMat1? .. " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";
	return TCL_ERROR;
    }

    // create the array
    UniaxialMaterial **theMats = new UniaxialMaterial *[numMat];
    if (theMats == 0) {
      opserr << "WARNING out of memory " <<
	"creating material array of size " << numMat <<
	"- element ZeroLength eleTag? iNode? jNode? " <<
	"-mat matID1? ... -dir dirMat1? .. " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";
	return TCL_ERROR;
    }

    // fill in the material array
    argi=6; 
    for (int i=0; i<numMat; i++) {

	int matID;	

	// read the material tag
	if (Tcl_GetInt(interp, argv[argi], &matID) != TCL_OK) {
	  opserr << "WARNING invalid matID " << argv[argi] <<
	    "- element ZeroLength eleTag? iNode? jNode? " <<
	    "-mat matID1? ... -dir dirMat1? .. " <<
	    "<-orient x1? x2? x3? y1? y2? y3?>\n";
	    delete [] theMats;	    
	    return TCL_ERROR;
	} else {

	    // get a pointer to the material from the modelbuilder	    
	    argi++;
	    UniaxialMaterial *theMat = theBuilder->getUniaxialMaterial(matID);
	    if (theMat == 0) {
	      opserr << "WARNING no material " << matID <<
		"exitsts - element ZeroLength eleTag? iNode? jNode? " <<
		"-mat matID1? ... -dir dirMat1? .. " <<
		"<-orient x1? x2? x3? y1? y2? y3?>\n"  ;
		delete [] theMats;		
		return TCL_ERROR;		
	    } else {
		
		// add the material to the array
		theMats[i] = theMat;
	    }
	}
    }
    

    // now read the dirn ID's for the materials added
    argi = 6 + numMat;
    if (strcmp(argv[argi],"-dir") != 0) {
      opserr << "WARNING expecting -dirn flag " << argv[argi] << 
	  "- element ZeroLength eleTag? iNode? jNode? " <<
	  "-mat matID1? ... -dir dirMat1? .. " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";

    delete [] theMats;
	return TCL_ERROR;	
    }    
    if ((argi + numMat) > argc) {
	opserr << "WARNING not enough directions provided for ele " << eleTag <<
	  "- element ZeroLength eleTag? iNode? jNode? " <<
	  "-mat matID1? ... -dir dirMat1? .. " <<
	  "<-orient x1? x2? x3? y1? y2? y3?>\n";
	
	delete [] theMats;	
	return TCL_ERROR;		
    }
    
    // create an ID to hold the directions
    ID theDirns(numMat);
    argi++; 
    int dirnID;	

    // read the dirn identifiers
    for (int j=0; j<numMat; j++) {
	if (Tcl_GetInt(interp, argv[argi], &dirnID) != TCL_OK) {
	  opserr << "WARNING invalid directiion " << argv[argi] <<
	    "- element ZeroLength eleTag? iNode? jNode? " <<
	    "-mat matID1? ... -dir dirMat1? .. " <<
	    "<-orient x1? x2? x3? y1? y2? y3?>\n";	
	  
	  delete [] theMats;	    
	  return TCL_ERROR;
	} else {
	  theDirns[j] = dirnID -1; // the minus g3 to C++
	  argi++;
	}
    }

    // create the vectors for the element orientation
    Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
    Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;

    // finally check the command line to see if user specified orientation
    if (argi < argc) {
	if (strcmp(argv[argi],"-orient") == 0) {
	    if (argc < (argi+7)) {
	      opserr << "WARNING not enough paramaters after -orient flag for ele " << eleTag <<
		"- element ZeroLength eleTag? iNode? jNode? " <<
		"-mat matID1? ... -dir dirMat1? .. " <<
		"<-orient x1? x2? x3? y1? y2? y3?>\n";	      
	      delete [] theMats;		
	      return TCL_ERROR;		
	    } else {
	      argi++;
	      double value;
	      // read the x values
	      for (int i=0; i<3; i++)  {
		if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
		  opserr << "WARNING invalid -orient value for ele  " << eleTag << argv[i] <<
		    "- element ZeroLength eleTag? iNode? jNode? " <<
		    "-mat matID1? ... -dir dirMat1? .. " <<
		    "<-orient x1? x2? x3? y1? y2? y3?>\n";	
		  delete [] theMats;			
		  return TCL_ERROR;
		    } else {
			argi++;
			x(i) = value;
		    }
		}
		// read the y values
		for (int j=0; j<3; j++)  {
		    if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
		      opserr << "WARNING invalid -orient value for ele  " <<
			eleTag << argv[argi] <<
			"- element ZeroLength eleTag? iNode? jNode? " <<
			"-mat matID1? ... -dir dirMat1? .. " <<
			"<-orient x1? x2? x3? y1? y2? y3?>\n";	
		      delete [] theMats;			
		      return TCL_ERROR;
		    } else {
		      argi++;
		      y(j) = value;		
		    }
		}
	    }
	}
    }
    
    //
    // now we create the element and add it to the domain
    //

    Element *theEle;
    theEle = new ZeroLength(eleTag, ndm, iNode, jNode, x, y, numMat, theMats, theDirns);
    if (theEle == 0) {
	delete [] theMats;	
	return TCL_ERROR;
    }
   
    if (theDomain->addElement(theEle) == false) {
	delete [] theMats;
	return TCL_ERROR;
    }

    // return the memory we stole and return OK
    delete [] theMats;    
    return TCL_OK;
}




int
TclModelBuilder_addZeroLengthSection(ClientData clientData, Tcl_Interp *interp,
			      int argc, TCL_Char **argv,
			      Domain *theDomain,
			      TclModelBuilder *theBuilder) {

    int ndm = theBuilder->getNDM(); // the spatial dimension of the problem

    //
    // first scan the command line to obtain eleID, iNode, jNode, material ID's
    // and their directions, and the orientation of ele xPrime and yPrime not
    // along the global x and y axis
    //
    
    int eleTag, iNode, jNode;
    
    // a quick check on number of args
    if (argc < 6) {
	opserr << "WARNING too few arguments " <<
	  "want - element zeroLengthSection eleTag? iNode? jNode? " <<
	  "secTag? " <<
	  "<-orient x1? x2? x3? y1? y2? y3?>\n";
	return TCL_ERROR;
    }

    // get the ele tag 
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
      opserr << "WARNING invalied eleTag " << argv[2] <<
	"- element zeroLengthSection eleTag? iNode? jNode? " <<
	"secTag? " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";	
	return TCL_ERROR;
    }

    // get the two end nodes
    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
      opserr << "WARNING invalied iNode " << argv[3] <<
	"- element zeroLengthSection eleTag? iNode? jNode? " <<
	"secTag? " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";	

      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode " << argv[4] <<
	"- element zeroLengthSection eleTag? iNode? jNode? " <<
	"secTag? " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";	
      return TCL_ERROR;
    }

	int secTag;

    if (Tcl_GetInt(interp, argv[5], &secTag) != TCL_OK) {
      opserr << "WARNING invalid secTag " << argv[5] <<
	"- element zeroLengthSection eleTag? iNode? jNode? " <<
	"secTag? " <<
	"<-orient x1? x2? x3? y1? y2? y3?>\n";	
      return TCL_ERROR;
    }

    // create the vectors for the element orientation
    Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
    Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;

	int argi = 6;

    // finally check the command line to see if user specified orientation
    if (argi < argc) {
      if (strcmp(argv[argi],"-orient") == 0) {
	if (argc < (argi+7)) {
	  opserr << "WARNING not enough paramaters after -orient flag for ele " <<
	    eleTag << "- element zeroLengthSection eleTag? iNode? jNode? secTag? " <<
	    "<-orient x1? x2? x3? y1? y2? y3?>\n";			
	  return TCL_ERROR;		
	} else {
	  argi++;
	  double value;
	  // read the x values
	  for (int i=0; i<3; i++)  {
	    if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
	      opserr << "WARNING invalid -orient value for ele  " <<
		eleTag << argv[argi] <<
		"- element zeroLengthSection eleTag? iNode? jNode secTag? " <<
		"<-orient x1? x2? x3? y1? y2? y3?>\n";	
	      return TCL_ERROR;
	    } else {
	      argi++;
	      x(i) = value;
	    }
	  }
	  // read the y values
	  for (int j=0; j<3; j++)  {
	    if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
	      opserr << "WARNING invalid -orient value for ele  " <<
		eleTag << argv[argi] <<
		"- element zeroLengthSection eleTag? iNode? jNode? secTag? " <<
		"<-orient x1? x2? x3? y1? y2? y3?>\n";	
	      return TCL_ERROR;
	    } else {
	      argi++;
	      y(j) = value;		
	    }
	  }
	}
      }
    }
    
    //
    // now we create the element and add it to the domain
    //

    SectionForceDeformation *theSection = theBuilder->getSection(secTag);
    
    if (theSection == 0) {
      opserr << "zeroLengthSection -- no section with tag " << secTag << " exists in Domain\n";
      return TCL_ERROR;		
    }
    
    Element *theEle = new ZeroLengthSection(eleTag, ndm, iNode, jNode, x, y, *theSection);
    
    if (theEle == 0)
      return TCL_ERROR;
    
    if (theDomain->addElement(theEle) == false)
      return TCL_ERROR;
    
    return TCL_OK;
}

/*
int
TclModelBuilder_addZeroLengthND(ClientData clientData, Tcl_Interp *interp,
			      int argc, TCL_Char **argv,
			      Domain *theDomain,
			      TclModelBuilder *theBuilder) {

    int ndm = theBuilder->getNDM(); // the spatial dimension of the problem

    //
    // first scan the command line to obtain eleID, iNode, jNode, material ID's
    // and their directions, and the orientation of ele xPrime and yPrime not
    // along the global x and y axis
    //
    
    int eleTag, iNode, jNode;
    
    // a quick check on number of args
    if (argc < 6) {
	opserr << "WARNING too few arguments %s %s %s\n" <<
				"want - element zeroLengthND eleTag? iNode? jNode? " <<
				"NDTag? <1DTag?>" <<
				"<-orient x1? x2? x3? y1? y2? y3?>\n";
	return TCL_ERROR;
    }

    // get the ele tag 
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
	opserr << "WARNING invalied eleTag %s %s %s %s\n" << 
				argv[2],
				"- element zeroLengthND eleTag? iNode? jNode? " <<
				"NDTag? <1DTag?>" <<
				"<-orient x1? x2? x3? y1? y2? y3?>\n";	
	return TCL_ERROR;
    }

    // get the two end nodes
    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
	opserr << "WARNING invalied iNode %s %s %s %s\n" << 
				argv[3],
				"- element zeroLengthND eleTag? iNode? jNode? " <<
				"NDTag? <1DTag?>" <<
				"<-orient x1? x2? x3? y1? y2? y3?>\n";	
	return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
	opserr << "WARNING invalid jNode %s %s %s %s\n" << 
				argv[4],
				"- element zeroLengthND eleTag? iNode? jNode? " <<
				"NDTag? <1DTag?>" <<
				"<-orient x1? x2? x3? y1? y2? y3?>\n";	
	return TCL_ERROR;
    }

	int NDTag;

    if (Tcl_GetInt(interp, argv[5], &NDTag) != TCL_OK) {
	opserr << "WARNING invalid NDTag %s %s %s %s\n" << 
				argv[5],
				"- element zeroLengthND eleTag? iNode? jNode? " <<
				"NDTag? <1DTag?>" <<
				"<-orient x1? x2? x3? y1? y2? y3?>\n";	
	return TCL_ERROR;
    }

	UniaxialMaterial *the1DMat = 0;

	int argi = 6;

	if (argc > 6 && strcmp(argv[6],"-orient") != 0) {

		int uniTag;

		if (Tcl_GetInt(interp, argv[6], &uniTag) != TCL_OK) {
			opserr << "WARNING invalid 1DTag %s %s %s %s\n" << 
					argv[5],
					"- element zeroLengthND eleTag? iNode? jNode? " <<
					"NDTag? <1DTag?>" <<
					"<-orient x1? x2? x3? y1? y2? y3?>\n";	
			return TCL_ERROR;
		}

		the1DMat = theBuilder->getUniaxialMaterial(uniTag);

		if (the1DMat == 0)
			opserr << "WARNING UniaxialMaterial %d not found in model, %s" <<
				uniTag, "proceeding without\n";

		argi = 7;
	}

    // create the vectors for the element orientation
    Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
    Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;

    // finally check the command line to see if user specified orientation
    if (argi < argc) {
	if (strcmp(argv[argi],"-orient") == 0) {
	    if (argc < (argi+7)) {
		opserr << "WARNING %s %d %s %s %s %s\n" << 
				"not enough paramaters after -orient flag for ele " <<
				eleTag,
				"- element zeroLengthND eleTag? iNode? jNode? " <<
				"NDTag? <1DTag?>" <<
				"<-orient x1? x2? x3? y1? y2? y3?>\n";	
		return TCL_ERROR;		
	    } else {
		argi++;
		double value;
		// read the x values
		for (int i=0; i<3; i++)  {
		    if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
			opserr << "WARNING %s %d %s %s %s %s\n" << 
						"invalid -orient value for ele  " <<
						eleTag,
						argv[argi],
						"- element zeroLengthND eleTag? iNode? jNode? " <<
						"NDTag? <1DTag?>" <<
						"<-orient x1? x2? x3? y1? y2? y3?>\n";	
			return TCL_ERROR;
		    } else {
			argi++;
			x(i) = value;
		    }
		}
		// read the y values
		for (int j=0; j<3; j++)  {
		    if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
			opserr << "WARNING %s %d %s %s %s %s\n" << 
						"invalid -orient value for ele  " <<
						eleTag,
						argv[argi],
						"- element zeroLengthND eleTag? iNode? jNode? " <<
						"NDTag? <1DTag?>" <<
						"<-orient x1? x2? x3? y1? y2? y3?>\n";	
			return TCL_ERROR;
		    } else {
			argi++;
			y(j) = value;		
		    }
		}
	    }
	}
    }
    
    //
    // now we create the element and add it to the domain
    //

	NDMaterial *theNDMat = theBuilder->getNDMaterial(NDTag);

	if (theNDMat == 0) {
		opserr << "%s -- no NDMaterial with tag %d exists in Domain\n" << 
			"zeroLengthND" << NDTag);	
		return TCL_ERROR;		
	}

	Element *theEle = 0;

	if (the1DMat == 0)
		theEle = new ZeroLengthND(eleTag, ndm, iNode, jNode, x, y, *theNDMat);
	else
		theEle = new ZeroLengthND(eleTag, ndm, iNode, jNode, x, y, *theNDMat, *the1DMat);
    
	if (theEle == 0)
		return TCL_ERROR;
   
    if (theDomain->addElement(theEle) == false)
		return TCL_ERROR;

    return TCL_OK;
}
*/
