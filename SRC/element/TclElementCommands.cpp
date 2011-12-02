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
                                                                        
// $Revision: 1.2 $
// $Date: 2000-10-28 05:50:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/TclElementCommands.cpp,v $
                                                                        
                                                                        
// File: ~/element/TclElementCommands.C
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the implementation of the TclElementCommands.
// The file contains the routine TclElementCommands which is invoked by the
// TclModelBuilder.
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <Domain.h>

#include <fElmt02.h>
#include <Truss.h>
#include <TrussSection.h>
#include <ElasticBeam2d.h>
#include <beam2d02.h>
#include <ElasticBeam3d.h>

#include <CrdTransf2d.h>
#include <CrdTransf3d.h>

#include <TclModelBuilder.h>

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

extern void printCommand(int argc, char **argv);

// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

extern int
TclModelBuilder_addFeapTruss(ClientData clientData, Tcl_Interp *interp,  int argc, 
			     char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addTruss(ClientData clientData, Tcl_Interp *interp,  int argc, 
			 char **argv, Domain*, TclModelBuilder *, int argStart);		      
		       
int
TclModelBuilder_addElasticBeam(ClientData clientData, Tcl_Interp *interp,  int argc, 
			       char **argv, Domain*, TclModelBuilder *, int argStart);	   

// GLF			       
extern int 
TclModelBuilder_addZeroLength(ClientData, Tcl_Interp *, int, char **,
			      Domain*, TclModelBuilder *);

// MHS			       
extern int 
TclModelBuilder_addZeroLengthSection(ClientData, Tcl_Interp *, int, char **,
			      Domain*, TclModelBuilder *);

// MHS
extern int 
TclModelBuilder_addZeroLengthND(ClientData, Tcl_Interp *, int, char **,
			      Domain*, TclModelBuilder *);


// REMO
extern int 
TclModelBuilder_addFrameElement(ClientData, Tcl_Interp *, int, char **,
				Domain*, TclModelBuilder *);
			
// MHS
extern int 
TclModelBuilder_addBeamWithHinges(ClientData, Tcl_Interp *, int, char **,
				  Domain*, TclModelBuilder *);
extern int 
TclModelBuilder_addFourNodeQuad(ClientData, Tcl_Interp *, int, char **,
				Domain*, TclModelBuilder *);
		   
int
TclModelBuilderElementCommand(ClientData clientData, Tcl_Interp *interp,
			      int argc, char **argv, 
			      Domain *theTclDomain, TclModelBuilder *theTclBuilder)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

  // check at least two arguments so don't segemnt fault on strcmp  
  if (argc < 2) {
    cerr << "WARNING need to specify an element type\n";
    cerr << "Want: element eleType <specific element args>\n";
    cerr << "Valid types: truss, elasticBeamColumn, nonlinearBeamColumn\n";
    return TCL_ERROR;
  }

  if (strcmp(argv[1],"fTruss") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addFeapTruss(clientData, interp, argc, argv,
					      theTclDomain, theTclBuilder, eleArgStart);
    return result;
  } else if (strcmp(argv[1],"truss") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addTruss(clientData, interp, argc, argv,
					      theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }else if (strcmp(argv[1],"elasticBeamColumn") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addElasticBeam(clientData, interp, argc, argv,
					      theTclDomain, theTclBuilder, eleArgStart);    
    return result;
  } else if (strcmp(argv[1],"nonlinearBeamColumn") == 0) {
    int result = TclModelBuilder_addFrameElement(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"beamWithHinges") == 0) {
	  int result = TclModelBuilder_addBeamWithHinges(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"quad") == 0) {
	  int result = TclModelBuilder_addFourNodeQuad(clientData, interp, argc, argv,
						       theTclDomain, theTclBuilder);
	  return result;
  } else if (strcmp(argv[1],"zeroLength") == 0) {
    int result = TclModelBuilder_addZeroLength(clientData, interp, argc, argv,
					       theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"zeroLengthSection") == 0) {
    int result = TclModelBuilder_addZeroLengthSection(clientData, interp, argc, argv,
					       theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"zeroLengthND") == 0) {
    int result = TclModelBuilder_addZeroLengthND(clientData, interp, argc, argv,
					       theTclDomain, theTclBuilder);
    return result;
  } else {
    cerr << "WARNING unknown element type: " <<  argv[1];
    cerr << "Valid types: truss, elasticBeamColumn, nonlinearBeamColumn, beamWithHinges, zeroLength\n";
    return TCL_ERROR;
  }
}



//
// command for the beam2d and beam3d elements
//


int
TclModelBuilder_addElasticBeam(ClientData clientData, Tcl_Interp *interp, int argc, 
			       char **argv, Domain *theTclDomain, TclModelBuilder *theTclBuilder,
			       int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - elasticBeamColumn \n";    
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  Element *theBeam = 0;

  if (ndm == 2) {
    // check plane frame problem has 3 dof per node
    if (ndf != 3) {
      cerr << "WARNING invalid ndf: " << ndf;
      cerr << ", for plane problem need 3 - elasticBeamColumn \n";    
      return TCL_ERROR;
    } 

    // check the number of arguments
    if ((argc-eleArgStart) < 8) {
      cerr << "WARNING bad command - want: elasticBeamColumn beamId iNode jNode A E I transTag\n";
      printCommand(argc, argv);
      return TCL_ERROR;
    }    

    // get the id, end nodes, and section properties
    int beamId, iNode, jNode, transTag;
    double A,E,I;
    if (Tcl_GetInt(interp, argv[1+eleArgStart], &beamId) != TCL_OK) {
      cerr << "WARNING invalid beamId: " << argv[1+eleArgStart];
      cerr << " - elasticBeamColumn beamId iNode jNode A E I\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
      cerr << "WARNING invalid iNode - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
      cerr << "WARNING invalid jNode - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4+eleArgStart], &A) != TCL_OK) {
      cerr << "WARNING invalid A - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5+eleArgStart], &E) != TCL_OK) {
      cerr << "WARNING invalid E - elasticBeam " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[6+eleArgStart], &I) != TCL_OK) {
      cerr << "WARNING invalid I - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }  
    if (Tcl_GetInt(interp, argv[7+eleArgStart], &transTag) != TCL_OK) {
      cerr << "WARNING invalid transTag - elasticBeamColumn " << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
    
    CrdTransf2d *theTrans = theTclBuilder->getCrdTransf2d(transTag);
    
    if (theTrans == 0) {
	cerr << "WARNING transformation object not found - elasticBeamColumn " << beamId;
	return TCL_ERROR;
    }
    
    // now create the beam and add it to the Domain
    theBeam = new ElasticBeam2d (beamId,A,E,I,iNode,jNode, *theTrans);
    
    if (theBeam == 0) {
      cerr << "WARNING ran out of memory creating beam - elasticBeamColumn ";	
      cerr << beamId << " iNode jNode A E I\n";
      return TCL_ERROR;
    }
  }

  else if (ndm == 3) {
    // check space frame problem has 6 dof per node
    if (ndf != 6) {
      cerr << "WARNING invalid ndof: " << ndf;
      cerr << ", for 3d problem  need 6 - elasticBeamColumn \n";    
      return TCL_ERROR;
    } 

    // check the number of arguments
    if ((argc-eleArgStart) < 11) {
      cerr << "WARNING bad command - want: elasticBeamColumn beamId iNode jNode";
      cerr << " A E G Jx Iy Iz transTag" << endl;
      printCommand(argc, argv);
      return TCL_ERROR;
    }    

    // get the id, end nodes, and section properties
    int beamId, iNode, jNode, transTag;
    double A,E,G,Jx,Iy,Iz;
    if (Tcl_GetInt(interp, argv[1+eleArgStart], &beamId) != TCL_OK) {
      cerr << "WARNING invalid beamId: " << argv[1+eleArgStart];
      cerr << " - elasticBeamColumn beamId iNode jNode A E G Jx Iy Iz\n ";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK) {
      cerr << "WARNING invalid iNode - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK) {
      cerr << "WARNING invalid jNode - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4+eleArgStart], &A) != TCL_OK) {
      cerr << "WARNING invalid A - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5+eleArgStart], &E) != TCL_OK) {
      cerr << "WARNING invalid E - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6+eleArgStart], &G) != TCL_OK) {
      cerr << "WARNING invalid G - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }  
    if (Tcl_GetDouble(interp, argv[7+eleArgStart], &Jx) != TCL_OK) {
      cerr << "WARNING invalid Jx - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }  
    if (Tcl_GetDouble(interp, argv[8+eleArgStart], &Iy) != TCL_OK) {
      cerr << "WARNING invalid Iy - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }  
    if (Tcl_GetDouble(interp, argv[9+eleArgStart], &Iz) != TCL_OK) {
      cerr << "WARNING invalid Iz - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }  
    if (Tcl_GetInt(interp, argv[10+eleArgStart], &transTag) != TCL_OK) {
      cerr << "WARNING invalid transTag - elasticBeamColumn " << beamId;
      cerr << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }      

    CrdTransf3d *theTrans = theTclBuilder->getCrdTransf3d(transTag);
    
    if (theTrans == 0) {
	cerr << "WARNING transformation object not found - elasticBeamColumn " << beamId;
	return TCL_ERROR;
    }
    
    // now create the beam and add it to the Domain
    theBeam = new ElasticBeam3d (beamId,A,E,G,Jx,Iy,Iz,iNode,jNode, *theTrans);
    
    if (theBeam == 0) {
      cerr << "WARNING ran out of memory creating beam - elasticBeamColumn ";	
      cerr << beamId << " iNode jNode A E G Jx Iy Iz\n";
      return TCL_ERROR;
    }
    
  }

  else {
    cerr << "WARNING elasticBeamColumn command only works when ndm is 2 or 3, ndm: ";
    cerr << ndm << endl;
    return TCL_ERROR;
  }

  // now add the beam to the domain
  if (theTclDomain->addElement(theBeam) == false) {
    cerr << "WARNING TclModelBuilder - addBeam - could not add beam to domain ";
    cerr << *theBeam;
    delete theBeam; // clean up the memory to avoid leaks
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}
