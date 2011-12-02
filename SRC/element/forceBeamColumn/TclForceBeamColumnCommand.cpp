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
                                                                        
// $Revision: 1.8 $
// $Date: 2003-04-03 23:40:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/TclForceBeamColumnCommand.cpp,v $
                                                                        
// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the implementation of the 
// TclModelBuilder_addDispBeamColumn() command. 

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <TclModelBuilder.h>

#include <ForceBeamColumn2d.h>
#include <ForceBeamColumn3d.h>

#include <LobattoBeamIntegration.h>
#include <UserDefinedBeamIntegration.h>

#include <HingeMidpointBeamIntegration2d.h>
#include <HingeMidpointBeamIntegration3d.h>
#include <HingeRadauBeamIntegration2d.h>
#include <HingeRadauBeamIntegration3d.h>
#include <UserDefinedHingeIntegration2d.h>
#include <UserDefinedHingeIntegration3d.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addForceBeamColumn(ClientData clientData, Tcl_Interp *interp,  
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
	 << " not compatible with forceBeamColumn element" << endln;
    return TCL_ERROR;
  }
  

  if (argc < 6) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? ...\n";
    return TCL_ERROR;
  }

  int eleTag, iNode, jNode, transfTag;
  CrdTransf2d *theTransf2d = 0;
  CrdTransf3d *theTransf3d = 0;
  Element *theElement = 0;


  if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
    opserr << "WARNING invalid forceBeamColumn eleTag" << endln;
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "forceBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode\n";
    opserr << "forceBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }


  //
  // fmk UNDOCUMENTED FEATURE - 
  // all to take similar command to nonlinearBeamColumn & dispBeamColumn 
  // 

  if ((strcmp(argv[6],"Lobatto") != 0) && (strcmp(argv[6],"HingeMidpoint") != 0) &&
      (strcmp(argv[6],"HingeRadau") != 0) && (strcmp(argv[6],"UserDefined") != 0) &&
      (strcmp(argv[6],"UserHinge") != 0)) {

    int nIP, secTag;

    if (Tcl_GetInt(interp, argv[5], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[6], &secTag) != TCL_OK) {
      opserr << "WARNING invalid secTag\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);
    if (theSection == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTag;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &transfTag) != TCL_OK) {
      opserr << "WARNING invalid transfTag\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    int argi = 8;
    int numIter = 0;
    double tol = 0.0;
    if (argc > argi) {
      if (strcmp(argv[argi],"-iter") == 0) {
	if (argc < argi+3) {
	  opserr << "WARNING not enough -iter args need -iter numIter? tol?\n";
	  opserr << "forceBeamColumn element: " << eleTag << endln;
	  return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[argi+1], &numIter) != TCL_OK) {
	  opserr << "WARNING invalid numIter\n";
	  opserr << "forceBeamColumn element: " << eleTag << endln;
	  return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[argi+2], &tol) != TCL_OK) {
	  opserr << "WARNING invalid numIter\n";
	  opserr << "forceBeamColumn element: " << eleTag << endln;
	  return TCL_ERROR;
	}
      }
    }
	
      
    if (ndm == 2) {
      
      theTransf2d = theTclBuilder->getCrdTransf2d(transfTag);
      
      if (theTransf2d == 0) {
	opserr << "WARNING transformation not found\n";
	opserr << "transformation: " << transfTag;
	opserr << "\nforceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
    }
    
    if (ndm == 3) {
      
      theTransf3d = theTclBuilder->getCrdTransf3d(transfTag);
      
      if (theTransf3d == 0) {
	opserr << "WARNING transformation not found\n";
	opserr << "transformation: " << transfTag;
	opserr << "\nforceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
    }

    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (int i = 0; i < nIP; i++)
      sections[i] = theSection;

    LobattoBeamIntegration beamIntegr;


    if (ndm == 2) {
      if (tol == 0.0)
	theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					   beamIntegr, *theTransf2d);
      else
	theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					   beamIntegr, *theTransf2d, 0.0, numIter, tol);
    }
    else {
      if (tol == 0.0)      
	theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					   beamIntegr, *theTransf3d);
      else
	theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					   beamIntegr, *theTransf3d, 0.0, numIter, tol);
    }

    delete [] sections;    
    if (theElement == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    if (theTclDomain->addElement(theElement) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      delete theElement;
      return TCL_ERROR;
    }

    return TCL_OK;
  } 

  
  //
  // otherwise use correct format of command as found in current documentation
  //

  if (Tcl_GetInt(interp, argv[5], &transfTag) != TCL_OK) {
    opserr << "WARNING invalid transfTag\n";
    opserr << "forceBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  if (ndm == 2) {
    
    theTransf2d = theTclBuilder->getCrdTransf2d(transfTag);
    
    if (theTransf2d == 0) {
      opserr << "WARNING transformation not found\n";
      opserr << "transformation: " << transfTag;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
  }
  
  if (ndm == 3) {
    
    theTransf3d = theTclBuilder->getCrdTransf3d(transfTag);
    
    if (theTransf3d == 0) {
      opserr << "WARNING transformation not found\n";
      opserr << "transformation: " << transfTag;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
  }
  
  if (strcmp(argv[6],"Lobatto") == 0) {
    int secTag, nIP;
    
    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? Lobatto secTag? nIP?\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &secTag) != TCL_OK) {
      opserr << "WARNING invalid secTag\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[8], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    SectionForceDeformation *theSection = theTclBuilder->getSection(secTag);
    if (theSection == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTag;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (int i = 0; i < nIP; i++)
      sections[i] = theSection;
    
    LobattoBeamIntegration beamIntegr;

    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf3d);
    delete [] sections;
  }

  else if (strcmp(argv[6],"HingeMidpoint") == 0 ||
	   strcmp(argv[6],"HingeRadau") == 0) {
    
    if (argc < 14) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? type secTagI? lpI? secTagJ? lpJ? E? A? Iz? <Iy? G? J?>\n";
      return TCL_ERROR;
    }

    int secTagI, secTagJ;
    double lpI, lpJ;
    double E, A, Iz, Iy, G, J;
    
    if (Tcl_GetInt(interp, argv[7], &secTagI) != TCL_OK) {
      opserr << "WARNING invalid secTagI\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &lpI) != TCL_OK) {
      opserr << "WARNING invalid lpI\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[9], &secTagJ) != TCL_OK) {
      opserr << "WARNING invalid secTagJ\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &lpJ) != TCL_OK) {
      opserr << "WARNING invalid lpJ\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[11], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[12], &A) != TCL_OK) {
      opserr << "WARNING invalid A\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[13], &Iz) != TCL_OK) {
      opserr << "WARNING invalid I\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    if (ndm == 3 && argc > 16) {
      if (Tcl_GetDouble(interp, argv[14], &Iy) != TCL_OK) {
	opserr << "WARNING invalid Iy\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[15], &G) != TCL_OK) {
	opserr << "WARNING invalid G\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[16], &J) != TCL_OK) {
	opserr << "WARNING invalid J\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
    }

    SectionForceDeformation *sectionI = theTclBuilder->getSection(secTagI);
    if (sectionI == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTagI;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    SectionForceDeformation *sectionJ = theTclBuilder->getSection(secTagJ);
    if (sectionJ == 0) {
      opserr << "WARNING section not found\n";
      opserr << "Section: " << secTagJ;
      opserr << "\nforceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    SectionForceDeformation *sections[2];
    sections[0] = sectionI;
    sections[1] = sectionJ;

    BeamIntegration *beamIntegr = 0;

    if (ndm == 2) {
      if (strcmp(argv[6],"HingeMidpoint") == 0)
	beamIntegr = new HingeMidpointBeamIntegration2d(E, A, Iz, lpI, lpJ);
      else
	beamIntegr = new HingeRadauBeamIntegration2d(E, A, Iz, lpI, lpJ);
      
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, 2, sections,
					 *beamIntegr, *theTransf2d);
    }
    else {
      if (strcmp(argv[6],"HingeMidpoint") == 0)
	beamIntegr =
	  new HingeMidpointBeamIntegration3d(E, A, Iz, Iy, G, J, lpI, lpJ);
      else
	beamIntegr =
	  new HingeRadauBeamIntegration3d(E, A, Iz, Iy, G, J, lpI, lpJ);
      
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, 2, sections,
					 *beamIntegr, *theTransf3d);
    }

    delete beamIntegr;
  }

  else if (strcmp(argv[6],"UserDefined") == 0) {

    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? UserDefined nIP? secTag1? ... pt1? ... wt1? ...\n";
      return TCL_ERROR;
    }

    int nIP;
    
    if (Tcl_GetInt(interp, argv[7], &nIP) != TCL_OK) {
      opserr << "WARNING invalid nIP\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    
    ID secs(nIP);
    Vector pts(nIP);
    Vector wts(nIP);

    int i, j;
    for (i = 0, j = 8; i < nIP; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	opserr << "WARNING invalid sec\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+nIP], &pt) != TCL_OK) {
	opserr << "WARNING invalid pt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+2*nIP], &wt) != TCL_OK) {
	opserr << "WARNING invalid wt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      secs(i) = sec;
      pts(i)  = pt;
      wts(i)  = wt;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (i = 0; i < nIP; i++) {
      SectionForceDeformation *theSection = theTclBuilder->getSection(secs(i));
      if (theSection == 0) {
	opserr << "WARNING section not found\n";
	opserr << "Section: " << secs(i);
	opserr << "\nforceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      sections[i] = theSection;
    }
    
    UserDefinedBeamIntegration beamIntegr(nIP, pts, wts);

    if (ndm == 2)
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf2d);
    else
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf3d);
    
    delete [] sections;
  }

  else if (strcmp(argv[6],"UserHinge") == 0) {

    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: element forceBeamColumn eleTag? iNode? jNode? transfTag? UserHinge E? A? Iz? <Iy? G? J?> npL? secTagL1? ... ptL1? ... wtL1? ... npR? secTagR1? ... ptR1? ... wtR1? ...\n";
      return TCL_ERROR;
    }

    double E, A, Iz, Iy, G, J;
    
    if (Tcl_GetDouble(interp, argv[7], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &A) != TCL_OK) {
      opserr << "WARNING invalid A\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[9], &Iz) != TCL_OK) {
      opserr << "WARNING invalid I\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    int argStart = 10;
    if (ndm == 3) {
      if (Tcl_GetDouble(interp, argv[10], &Iy) != TCL_OK) {
	opserr << "WARNING invalid I\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[11], &G) != TCL_OK) {
	opserr << "WARNING invalid G\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[12], &J) != TCL_OK) {
	opserr << "WARNING invalid J\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      argStart = 13;
    }

    int npL, npR;
      
    if (Tcl_GetInt(interp, argv[argStart], &npL) != TCL_OK) {
      opserr << "WARNING invalid npL\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[argStart+3*npL+1], &npR) != TCL_OK) {
      opserr << "WARNING invalid npR\n";
      opserr << "forceBeamColumn element: " << eleTag << endln;
      return TCL_ERROR;
    }

    int nIP = npL+npR;

    ID secs(nIP);
    Vector ptsL(npL);
    Vector wtsL(npL);
    Vector ptsR(npR);
    Vector wtsR(npR);

    int i, j;
    for (i = 0, j = argStart+1; i < npL; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	opserr << "WARNING invalid sec\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+npL], &pt) != TCL_OK) {
	opserr << "WARNING invalid pt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+2*npL], &wt) != TCL_OK) {
	opserr << "WARNING invalid wt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      secs(i) = sec;
      ptsL(i) = pt;
      wtsL(i) = wt;
    }
    for (i = 0, j = 1+(argStart+1)+3*npL; i < npR; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
	opserr << "WARNING invalid sec\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+npR], &pt) != TCL_OK) {
	opserr << "WARNING invalid pt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j+2*npR], &wt) != TCL_OK) {
	opserr << "WARNING invalid wt\n";
	opserr << "forceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      secs(i+npL) = sec;
      ptsR(i)     = pt;
      wtsR(i)     = wt;
    }
    
    SectionForceDeformation **sections = new SectionForceDeformation *[nIP];
    for (i = 0; i < nIP; i++) {
      SectionForceDeformation *theSection = theTclBuilder->getSection(secs(i));
      if (theSection == 0) {
	opserr << "WARNING section not found\n";
	opserr << "Section: " << secs(i);
	opserr << "\nforceBeamColumn element: " << eleTag << endln;
	return TCL_ERROR;
      }
      sections[i] = theSection;
    }
    
    if (ndm == 2) {
      UserDefinedHingeIntegration2d beamIntegr(npL, ptsL, wtsL,
					       npR, ptsR, wtsR,
					       E, A, Iz);
    
      theElement = new ForceBeamColumn2d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf2d);
    }
    else {
      UserDefinedHingeIntegration3d beamIntegr(npL, ptsL, wtsL,
					       npR, ptsR, wtsR,
					       E, A, Iz, Iy, G, J);
    
      theElement = new ForceBeamColumn3d(eleTag, iNode, jNode, nIP, sections,
					 beamIntegr, *theTransf3d);
    }
    
    delete [] sections;
  }

  else {
    opserr << "Unknown integration type: " << argv[6] << endln;
    opserr << "forceBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "forceBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  if (theTclDomain->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "forceBeamColumn element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }
  
  // if get here we have sucessfully created the element and added it to the domain
  return TCL_OK;
}
