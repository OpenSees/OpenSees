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
                                                                        
// $Revision$
// $Date$
// $URL$
                                                                        
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the elasticBeamColumn element.

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <ElasticBeam2d.h>
#include <ElasticBeam3d.h>
#include <SectionForceDeformation.h>

#include <CrdTransf.h>
#include <CrdTransf.h>

#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addElasticBeam(ClientData clientData, Tcl_Interp *interp, int argc, 
			       TCL_Char **argv, Domain *theTclDomain, TclModelBuilder *theTclBuilder,
			       int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeamColumn \n";    
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  Element *theBeam = 0;

  if (ndm == 2) {
    // check plane frame problem has 3 dof per node
    if (ndf != 3) {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for plane problem need 3 - elasticBeamColumn \n";    
      return TCL_ERROR;
    } 

    // check the number of arguments
    if ((argc-eleArgStart) < 8) {
      opserr << "WARNING bad command - want: elasticBeamColumn beamId iNode jNode A E I <alpha> <d> transTag <-mass m> <-cMass>\n";
      printCommand(argc, argv);
      return TCL_ERROR;
    }    

    // get the id, end nodes, and section properties
    int beamId, iNode, jNode, transTag;
    double A,E,I;
    if (Tcl_GetInt(interp, argv[1+eleArgStart], &beamId) != TCL_OK) {
      opserr << "WARNING invalid beamId: " << argv[1+eleArgStart];
      opserr << " - elasticBeamColumn beamId iNode jNode A E I\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
      opserr << "WARNING invalid iNode - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4+eleArgStart], &A) != TCL_OK) {
      opserr << "WARNING invalid A - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5+eleArgStart], &E) != TCL_OK) {
      opserr << "WARNING invalid E - elasticBeam " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[6+eleArgStart], &I) != TCL_OK) {
      opserr << "WARNING invalid I - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }  

    double alpha = 0.0;
    double d = 0.0;
    double mass = 0.0;
    int cMass = 0;
    int argi = 0;

    if ( ((argc-eleArgStart) == 10) && (strcmp(argv[eleArgStart+8],"-mass") != 0)){    

      if (Tcl_GetDouble(interp, argv[7+eleArgStart], &alpha) != TCL_OK) {
	opserr << "WARNING invalid alpha - elasticBeamColumn " << beamId << " iNode jNode A E I alpha d \n";
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[8+eleArgStart], &d) != TCL_OK) {
	opserr << "WARNING invalid d - elasticBeamColumn " << beamId << " iNode jNode A E I alpha d \n";
	return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[9+eleArgStart], &transTag) != TCL_OK) {
	opserr << "WARNING invalid transTag - elasticBeamColumn " << beamId << " iNode jNode A E I alpha d transTag\n";
	return TCL_ERROR;
      }
      argi = 10;
    } 

    else {
      if (Tcl_GetInt(interp, argv[7+eleArgStart], &transTag) != TCL_OK) {
	opserr << "WARNING invalid transTag - elasticBeamColumn " << beamId << " iNode jNode A E I transTag\n";
	return TCL_ERROR;
      }
      argi = 8;
    }

    CrdTransf *theTrans = OPS_getCrdTransf(transTag);
    
    if (theTrans == 0) {
	opserr << "WARNING transformation object not found - elasticBeamColumn " << beamId;
	return TCL_ERROR;
    }

    while (argi < argc) {
      if (strcmp(argv[argi],"-mass") == 0) {
	if (argc < argi+2) {
	  opserr << "WARNING not enough -mass args need -mass mass??\n";
	  opserr << argv[1] << " element: " << beamId << endln;
	  return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[argi+1], &mass) != TCL_OK) {
	  opserr << "WARNING invalid mass\n";
	  opserr << argv[1] << " element: " << beamId << endln;
	  return TCL_ERROR;
	}
	argi += 2;
      } else if ((strcmp(argv[argi],"-lMass") == 0) || (strcmp(argv[argi],"lMass") == 0)) {
	cMass = 0;  // lumped mass matrix (default)
    argi++;
      } else if ((strcmp(argv[argi],"-cMass") == 0) || (strcmp(argv[argi],"cMass") == 0)) {
	cMass = 1;  // consistent mass matrix
    argi++;
      } else
	argi++;
    }

    
    // now create the beam and add it to the Domain
    theBeam = new ElasticBeam2d (beamId,A,E,I,iNode,jNode, *theTrans, alpha, d, mass, cMass);
    
    if (theBeam == 0) {
      opserr << "WARNING ran out of memory creating beam - elasticBeamColumn ";	
      opserr << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
  }

  else if (ndm == 3) {
    // check space frame problem has 6 dof per node
    if (ndf != 6) {
      opserr << "WARNING invalid ndof: " << ndf;
      opserr << ", for 3d problem  need 6 - elasticBeamColumn \n";    
      return TCL_ERROR;
    } 

    // check the number of arguments
    if (((argc-eleArgStart) < 11) && ((argc-eleArgStart) != 6)) {
      opserr << "WARNING bad command - want: elasticBeamColumn beamId iNode jNode";
      opserr << " A E G Jx Iy Iz transTag <-mass m> <-cMass>" << endln;
      printCommand(argc, argv);
      return TCL_ERROR;
    }    

    // get the id, end nodes, and section properties
    int beamId, iNode, jNode, transTag;
    double A,E,G,Jx,Iy,Iz;
    if (Tcl_GetInt(interp, argv[1+eleArgStart], &beamId) != TCL_OK) {
      opserr << "WARNING invalid beamId: " << argv[1+eleArgStart];
      opserr << " - elasticBeamColumn beamId iNode jNode A E G Jx Iy Iz\n ";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
      opserr << "WARNING invalid iNode - elasticBeamColumn " << beamId;
      opserr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode - elasticBeamColumn " << beamId;
      opserr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }

    if ((argc-eleArgStart) == 6) {
      int section;
      if (Tcl_GetInt(interp, argv[4+eleArgStart], &section) != TCL_OK) {
	opserr << "WARNING invalid secTag - elasticBeamColumn iNode jNode sectionTag? transfTag?\n";
	return TCL_ERROR;
      }
      
      SectionForceDeformation *theSection = theTclBuilder->getSection(section);

      if (Tcl_GetInt(interp, argv[5+eleArgStart], &transTag) != TCL_OK) {
	opserr << "WARNING invalid transTag - elasticBeamColumn " << beamId;
	opserr << " iNode jNode sectionTag? transfTag?\n";
	return TCL_ERROR;
      }      
      
      CrdTransf *theTrans = OPS_getCrdTransf(transTag);
    
      if (theTrans == 0) {
	opserr << "WARNING transformation object not found - elasticBeamColumn " << beamId;
	return TCL_ERROR;
      }

      double mass = 0.0;
      int cMass = 0;
      int argi = 6+eleArgStart;

      while (argi < argc) {
	if (strcmp(argv[argi],"-mass") == 0) {
	  if (argc < argi+2) {
	    opserr << "WARNING not enough -mass args need -mass mass??\n";
	    opserr << argv[1] << " element: " << beamId << endln;
	    return TCL_ERROR;
	  }
	  if (Tcl_GetDouble(interp, argv[argi+1], &mass) != TCL_OK) {
	    opserr << "WARNING invalid mass\n";
	    opserr << argv[1] << " element: " << beamId << endln;
	    return TCL_ERROR;
	  }
	  argi += 2;
    } else if ((strcmp(argv[argi],"-lMass") == 0) || (strcmp(argv[argi],"lMass") == 0)) {
	  cMass = 0;  // lumped mass matrix (default)
      argi++;
    } else if ((strcmp(argv[argi],"-cMass") == 0) || (strcmp(argv[argi],"cMass") == 0)) {
	  cMass = 1;  // consistent mass matrix
      argi++;
	} else
	  argi++;
      }
      
      // now create the beam and add it to the Domain
      theBeam = new ElasticBeam3d (beamId, iNode, jNode, theSection, *theTrans, mass, cMass);      

    } else {

      if (Tcl_GetDouble(interp, argv[4+eleArgStart], &A) != TCL_OK) {
	opserr << "WARNING invalid A - elasticBeamColumn " << beamId;
	opserr << " iNode jNode A E G Jx Iy Iz\n";
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[5+eleArgStart], &E) != TCL_OK) {
	opserr << "WARNING invalid E - elasticBeamColumn " << beamId;
	opserr << " iNode jNode A E G Jx Iy Iz\n";
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[6+eleArgStart], &G) != TCL_OK) {
	opserr << "WARNING invalid G - elasticBeamColumn " << beamId;
	opserr << " iNode jNode A E G Jx Iy Iz\n";
	return TCL_ERROR;
      }  
      if (Tcl_GetDouble(interp, argv[7+eleArgStart], &Jx) != TCL_OK) {
	opserr << "WARNING invalid Jx - elasticBeamColumn " << beamId;
	opserr << " iNode jNode A E G Jx Iy Iz\n";
	return TCL_ERROR;
      }  
      if (Tcl_GetDouble(interp, argv[8+eleArgStart], &Iy) != TCL_OK) {
	opserr << "WARNING invalid Iy - elasticBeamColumn " << beamId;
	opserr << " iNode jNode A E G Jx Iy Iz\n";
	return TCL_ERROR;
      }  
      if (Tcl_GetDouble(interp, argv[9+eleArgStart], &Iz) != TCL_OK) {
	opserr << "WARNING invalid Iz - elasticBeamColumn " << beamId;
	opserr << " iNode jNode A E G Jx Iy Iz\n";
	return TCL_ERROR;
      }  
      if (Tcl_GetInt(interp, argv[10+eleArgStart], &transTag) != TCL_OK) {
	opserr << "WARNING invalid transTag - elasticBeamColumn " << beamId;
	opserr << " iNode jNode A E G Jx Iy Iz\n";
	return TCL_ERROR;
      }      
      
      CrdTransf *theTrans = OPS_getCrdTransf(transTag);
    
      if (theTrans == 0) {
	opserr << "WARNING transformation object not found - elasticBeamColumn " << beamId;
	return TCL_ERROR;
      }



      double mass = 0.0;
      int cMass = 0;
      int argi = 11+eleArgStart;
      
      while (argi < argc) {
	if (strcmp(argv[argi],"-mass") == 0) {
	  if (argc < argi+2) {
	    opserr << "WARNING not enough -mass args need -mass mass??\n";
	    opserr << argv[1] << " element: " << beamId << endln;
	    return TCL_ERROR;
	  }
	  if (Tcl_GetDouble(interp, argv[argi+1], &mass) != TCL_OK) {
	    opserr << "WARNING invalid mass\n";
	    opserr << argv[1] << " element: " << beamId << endln;
	    return TCL_ERROR;
	  }
	  argi += 2;
    } else if ((strcmp(argv[argi],"-lMass") == 0) || (strcmp(argv[argi],"lMass") == 0)) {
	  cMass = 0;  // lumped mass matrix (default)
      argi++;
    } else if ((strcmp(argv[argi],"-cMass") == 0) || (strcmp(argv[argi],"cMass") == 0)) {
	  cMass = 1;  // consistent mass matrix
      argi++;
	} else
	  argi++;
      }

      
      // now create the beam and add it to the Domain
      theBeam = new ElasticBeam3d (beamId,A,E,G,Jx,Iy,Iz,iNode,jNode, *theTrans, mass, cMass);
    }

    if (theBeam == 0) {
      opserr << "WARNING ran out of memory creating beam - elasticBeamColumn ";	
      opserr << beamId << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }
  }

  else {
    opserr << "WARNING elasticBeamColumn command only works when ndm is 2 or 3, ndm: ";
    opserr << ndm << endln;
    return TCL_ERROR;
  }

  // now add the beam to the domain
  if (theTclDomain->addElement(theBeam) == false) {
    opserr << "WARNING TclModelBuilder - addBeam - could not add beam to domain ";
    opserr << *theBeam;
    delete theBeam; // clean up the memory to avoid leaks
    return TCL_ERROR;
  }

  // if get here we have successfully created the node and added it to the domain
  return TCL_OK;
}
