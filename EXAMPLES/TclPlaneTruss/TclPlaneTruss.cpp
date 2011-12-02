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
// $Date: 2003-02-19 21:49:00 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/TclPlaneTruss/TclPlaneTruss.cpp,v $
                                                                        
// Written: fmk 
// Created: 04/98
// Revision: A
//
// Description: This file contains the class implementation for TclPlaneTruss
// TclPlaneTruss is a class for building a Plane Frame model in an interpreted
// enviroment. The constructor is used to add new commands to the interpreter,
// these commands are also defined in this file.
//
// What: "@(#) TclPlaneTruss.C, revA"

#include <stdlib.h>
#include <string.h>

#include <Domain.h>
#include <Node.h>
#include <Truss.h>
#include "MyTruss.h"
#include <fElmt02.h>

#include <ElasticMaterial.h>
#include <ElasticPPMaterial.h>
#include <ParallelMaterial.h>
#include <SP_Constraint.h>
#include <NodalLoad.h>
#include <ArrayOfTaggedObjects.h>
#include <LoadPattern.h>

#include "TclPlaneTruss.h"

//
// some static variables used in the functions
//

static Domain *theTclPlaneTrussDomain = 0;
static TclPlaneTruss *theTclPlaneTruss =0;
extern LoadPattern *theTclLoadPattern;
static int numSPs = 0;
static int loadTag = 0;
// 
// the functions that will be invoked by the interpreter while building the model
//

int
TclPlaneTruss_addNode(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv);

int
TclPlaneTruss_addTruss(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv);
		      
int
TclPlaneTruss_addMyTruss(ClientData clientData, Tcl_Interp *interp, int argc, 
			 char **argv);		      
			 
int
TclPlaneTruss_addfTruss(ClientData clientData, Tcl_Interp *interp, int argc, 
			char **argv);		      			 
		      
int
TclModelBuilder_addMyPattern(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv);

int
TclPlaneTruss_addMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
		    char **argv);

int
TclPlaneTruss_addSP(ClientData clientData, Tcl_Interp *interp, int argc,   
		    char **argv);

int
TclPlaneTruss_addNodalLd(ClientData clientData, Tcl_Interp *interp, int argc,   
			 char **argv);

int
TclPlaneTruss_done(ClientData clientData, Tcl_Interp *interp, int argc,   
		   char **argv);


//
// the class constructor, destructor and methods
//

TclPlaneTruss::TclPlaneTruss(Domain &theDomain, Tcl_Interp *interp)
  :ModelBuilder(theDomain)
{
  theMaterials = new ArrayOfTaggedObjects(32);

  // set the static pointer used in the class
  theTclPlaneTrussDomain = &theDomain;
  
  // call Tcl_CreateCommand for class specific commands
  Tcl_CreateCommand(interp, "node", TclPlaneTruss_addNode,
		    (ClientData) NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "truss", TclPlaneTruss_addTruss,
		    (ClientData) NULL, (Tcl_CmdDeleteProc *)NULL);
  
  Tcl_CreateCommand(interp, "myTruss", TclPlaneTruss_addMyTruss,
		    (ClientData) NULL, (Tcl_CmdDeleteProc *)NULL);  
  
  Tcl_CreateCommand(interp, "fTruss", TclPlaneTruss_addfTruss,
		    (ClientData) NULL, (Tcl_CmdDeleteProc *)NULL);    
  
  Tcl_CreateCommand(interp, "material", TclPlaneTruss_addMaterial,
		    (ClientData) NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "fix", TclPlaneTruss_addSP,
		    (ClientData) NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "pattern", TclModelBuilder_addMyPattern,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "load", TclPlaneTruss_addNodalLd,
		    (ClientData) NULL, (Tcl_CmdDeleteProc *)NULL);

  Tcl_CreateCommand(interp, "doneModel", TclPlaneTruss_done,
		    (ClientData) NULL, (Tcl_CmdDeleteProc *)NULL);

  // set the static pointer in this file
  theTclPlaneTruss = this;
  numSPs = 0;
  theTclLoadPattern = 0;  
}

TclPlaneTruss::~TclPlaneTruss()
{
  // delete the material objects
  theMaterials->clearAll();

  delete theMaterials;
  theTclLoadPattern =0;
  
  // may possibly invoke Tcl_DeleteCommand() later  
}
    
int 
TclPlaneTruss::buildFE_Model(void)
{
  return 0;
}


int 
TclPlaneTruss::addMaterial(UniaxialMaterial &theMaterial)
{
  bool result = theMaterials->addComponent(&theMaterial);
  if (result == true)
    return 0;
  else
    return -1;
}

UniaxialMaterial *
TclPlaneTruss::getMaterial(int tag)
{
  TaggedObject *mc = theMaterials->getComponentPtr(tag);
  if (mc == 0) return 0;
  UniaxialMaterial *result = (UniaxialMaterial *)mc;
  return result;
}



//
// the functions for adding nodes, beams, constraints and loads to the model
// that are called by the interpreter
//


int
TclPlaneTruss_addNode(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 4) {
      interp->result = "WARNING bad command - want node node_id x_loc y_loc <m_x m_y>";
      return TCL_ERROR;
  }    

  // get the id, x_loc and y_loc
  int nodeId;
  double xLoc, yLoc;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
      interp->result = "WARNING invalid node_id - node node_id x_loc y_loc ";
      return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      interp->result = "WARNING invalid x_loc - node node_id x_loc y_loc ";
      return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      interp->result = "WARNING invalid x_loc - node node_id x_loc y_loc ";
      return TCL_ERROR;
  }

  
  // now create the node and add it to the Domain
  Node *theNode = new Node(nodeId,2,xLoc,yLoc);
  if (theNode == 0) {
    opserr << "WARNING TclPlaneTruss - addNode - ran out of memory for node ";
    opserr << nodeId << endln;
    return TCL_ERROR;
  }
  
  // the node may have some mass associated with it
  if (argc == 6) {
      double mx, my;
      if (Tcl_GetDouble(interp, argv[4], &mx) != TCL_OK) {
	  interp->result = "WARNING invalid x_mass - node node_id x_loc y_loc x_mass y_mass";
	  return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[5], &my) != TCL_OK) {
	  interp->result = "WARNING invalid x_loc - node node_id x_loc y_loc ";
	  return TCL_ERROR;
      }      
      
      Matrix mass(2,2);
      mass(0,0) = mx;
      mass(1,1) = my;
      theNode->setMass(mass);
  }  
  

  if (theTclPlaneTrussDomain->addNode(theNode) == false) {
    opserr << "WARNING TclPlaneTruss - addNode - could not add node to domain ";
    opserr << nodeId << endln;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}


// 
// to create a beam of type beam2d02 and add to the domain
//
int
TclPlaneTruss_addTruss(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc != 6) {
      interp->result = "WARNING bad command - want: truss eleId iNode jNode A matID";
      return TCL_ERROR;
  }    

  // get the id, x_loc and y_loc
  int trussId, iNode, jNode, matID;
  double A;
  if (Tcl_GetInt(interp, argv[1], &trussId) != TCL_OK) {
     interp->result = "WARNING invalid trussId- truss trussId iNode jNode A matID";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &iNode) != TCL_OK) {
     interp->result = "WARNING invalid iNode- truss trussId iNode jNode A matID";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &jNode) != TCL_OK) {
     interp->result = "WARNING invalid jNode- truss trussId iNode jNode A matID";
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &A) != TCL_OK) {
     interp->result = "WARNING invalid A- truss trussId iNode jNode A matID";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5], &matID) != TCL_OK) {
     interp->result = "WARNING invalid matId- truss trussId iNode jNode A matID";
      return TCL_ERROR;
  }
  
  UniaxialMaterial *theMaterial = theTclPlaneTruss->getMaterial(matID);

  if (theMaterial == 0) {
    opserr << "WARNING TclPlaneTruss - truss - no Material found with tag ";
    opserr << matID << endln;
    return TCL_ERROR;
  }

  // now create the truss and add it to the Domain
  Truss *theTruss = new Truss(trussId,2,iNode,jNode,*theMaterial,A);
  if (theTruss == 0) {
    opserr << "WARNING TclPlaneTruss - addTruss - ran out of memory for node ";
    opserr << trussId << endln;
    return TCL_ERROR;
  }
  if (theTclPlaneTrussDomain->addElement(theTruss) == false) {
    opserr << "WARNING TclPlaneTruss - addTruss - could not add Truss to domain ";
    opserr << trussId << endln;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}


// 
// to create an element of type MyTruss
//
int
TclPlaneTruss_addMyTruss(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc != 6 && argc != 7) {
      interp->result = "WARNING bad command - want: myTruss eleId iNode jNode Area matID";
      return TCL_ERROR;
  }    

  // get the id, x_loc and y_loc
  int trussId, iNode, jNode, matID;
  double A, M = 0.0;
  if (Tcl_GetInt(interp, argv[1], &trussId) != TCL_OK) {
     interp->result = "WARNING invalid eleId - myTruss eleId iNode jNode Area matID";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &iNode) != TCL_OK) {
     interp->result = "WARNING invalid iNode- myTruss eleId iNode jNode Area matID";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &jNode) != TCL_OK) {
     interp->result = "WARNING invalid jNode- myTruss eleId iNode jNode Area matID";
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &A) != TCL_OK) {
     interp->result = "WARNING invalid A- myTruss eleId iNode jNode Area matID";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5], &matID) != TCL_OK) {
     interp->result = "WARNING invalid matId- myTruss eleId iNode jNode Area matID";
      return TCL_ERROR;
  }
  if (argc == 7 && Tcl_GetDouble(interp, argv[6], &M) != TCL_OK) {
     interp->result = "WARNING invalid matId- myTruss eleId iNode jNode Area matID";
     return TCL_ERROR;
  }  
  
  UniaxialMaterial *theMaterial = theTclPlaneTruss->getMaterial(matID);

  if (theMaterial == 0) {
    opserr << "WARNING TclPlaneTruss - truss - no Material found with tag ";
    opserr << matID << endln;
    return TCL_ERROR;
  }

  // now create the truss and add it to the Domain
  MyTruss *theTruss = new MyTruss(trussId,iNode,jNode,*theMaterial,A,M);
  if (theTruss == 0) {
    opserr << "WARNING TclPlaneTruss - addMyTruss - ran out of memory for node ";
    opserr << trussId << endln;
    return TCL_ERROR;
  }
  if (theTclPlaneTrussDomain->addElement(theTruss) == false) {
    delete theTruss;
    opserr << "WARNING TclPlaneTruss - addTruss - could not add Truss to domain ";
    opserr << trussId << endln;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}





// 
// to create an element of type fElmt02 - fTruss
//
int
TclPlaneTruss_addfTruss(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc != 6 && argc != 7) {
      interp->result = "WARNING bad command - want: myTruss eleId iNode jNode Area E <rho>";
      return TCL_ERROR;
  }    

  // get the id, x_loc and y_loc
  int trussId, iNode, jNode;
  double A, E, M = 0.0;
  if (Tcl_GetInt(interp, argv[1], &trussId) != TCL_OK) {
     interp->result = "WARNING invalid eleId - myTruss eleId iNode jNode Area E";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &iNode) != TCL_OK) {
     interp->result = "WARNING invalid iNode- myTruss eleId iNode jNode Area E";
     return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &jNode) != TCL_OK) {
     interp->result = "WARNING invalid jNode- myTruss eleId iNode jNode Area E";
     return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &A) != TCL_OK) {
     interp->result = "WARNING invalid A- myTruss eleId iNode jNode Area E";
     return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5], &E) != TCL_OK) {
     interp->result = "WARNING invalid E - myTruss eleId iNode jNode Area E";
      return TCL_ERROR;
  }
  if (argc == 7 && Tcl_GetDouble(interp, argv[6], &M) != TCL_OK) {
     interp->result = "WARNING invalid rho - myTruss eleId iNode jNode Area E rho";
     return TCL_ERROR;
  }  
  
  // now create the truss and add it to the Domain
  fElmt02 *theTruss = new fElmt02(trussId,iNode,jNode,A,E,M);
  if (theTruss == 0) {
    opserr << "WARNING TclPlaneTruss - addMyTruss - ran out of memory for node ";
    opserr << trussId << endln;
    return TCL_ERROR;
  }
  if (theTclPlaneTrussDomain->addElement(theTruss) == false) {
    delete theTruss;
    opserr << "WARNING TclPlaneTruss - addTruss - could not add Truss to domain ";
    opserr << trussId << endln;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclPlaneTruss_addMaterial(ClientData clientData, 
			  Tcl_Interp *interp, int argc,    
			  char **argv)
{
  // make sure at least one other argument to contain integrator
  if (argc < 2) {
      interp->result = "WARNING need to specify a Material type ";
      return TCL_ERROR;
  }    

  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1],"Elastic") == 0) {
    if (argc != 4) {
      interp->result = "WARNING want: material Elastic tag E ";
      return TCL_ERROR;
    }    
    int tag;
    double E;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)	
      return TCL_ERROR;		
    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK)	
      return TCL_ERROR;	


    UniaxialMaterial *theMaterial = new ElasticMaterial(tag,E);       
    if (theTclPlaneTruss->addMaterial(*theMaterial) < 0)
      return TCL_ERROR;
    else
      return TCL_OK;
  }
  else if (strcmp(argv[1],"ElasticPP") == 0) {
    if (argc != 5) {
      interp->result = "WARNING want: material Elastic tag E ep";
      return TCL_ERROR;
    }    
    int tag;
    double E, ep;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)	
      return TCL_ERROR;		
    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK)	
      return TCL_ERROR;	
    if (Tcl_GetDouble(interp, argv[4], &ep) != TCL_OK)	
      return TCL_ERROR;	


    UniaxialMaterial *theMaterial = new ElasticPPMaterial(tag,E,ep);       
    if (theTclPlaneTruss->addMaterial(*theMaterial) < 0)
      return TCL_ERROR;
    else
      return TCL_OK;
  }


  else if (strcmp(argv[1],"Parallel") == 0) {
    if (argc < 3) {
      interp->result = "WARNING want: material Parallel tag matTags";
      return TCL_ERROR;
    }    
    // get material tag
    int tag;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)	
      return TCL_ERROR;		

    // get tags of material in parallel model
    int numMaterials = argc-3;

    UniaxialMaterial **theMats = new UniaxialMaterial *[numMaterials];

    for (int i=0; i<numMaterials; i++) {
      int matTag;
      if (Tcl_GetInt(interp, argv[3+i], &matTag) != TCL_OK)	
	return TCL_ERROR;		

      UniaxialMaterial *theMaterial = 
	theTclPlaneTruss->getMaterial(matTag);

      if (theMaterial == 0) {
	interp->result = "WARNING material Parallel tags - unknown tag ";
	return TCL_ERROR;
      }    

      theMats[i] = theMaterial;

    }

    // create the parallel model
    UniaxialMaterial *theMaterial = 
      new ParallelMaterial(tag, numMaterials, theMats);       

    delete [] theMats;

    if (theTclPlaneTruss->addMaterial(*theMaterial) < 0)
      return TCL_ERROR;
    else
      return TCL_OK;
  }
  else {
    interp->result = "WARNING No Material type exists ";
    return TCL_ERROR;
  }    
  return TCL_OK;
}



extern int TclPatternCommand(ClientData clientData, 
			     Tcl_Interp *interp, int argc,    
			     char **argv, Domain *theDomain);
int
TclModelBuilder_addMyPattern(ClientData clientData, Tcl_Interp *interp, 
			   int argc, char **argv)
			  
{
  return TclPatternCommand(clientData, interp, 
			   argc, argv, theTclPlaneTrussDomain);
}




int
TclPlaneTruss_addSP(ClientData clientData, Tcl_Interp *interp, int argc,   
		    char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 4) {
      interp->result = "WARNING bad command - want: fix nodeId fixX fixY";
      return TCL_ERROR;
  }    

  // get the id, x_loc and y_loc
  int nodeId, xFix, yFix;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
      interp->result = "WARNING invalid nodeId - fix nodeId fixX fixY";
      return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &xFix) != TCL_OK) {
      interp->result = "WARNING invalid fixX - fix nodeId fixX fixY";
      return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &yFix) != TCL_OK) {
      interp->result = "WARNING invalid fixY - fix nodeId fixX fixY";
      return TCL_ERROR;
  }

  // now for every constrained dof create a SP_Constraint
  SP_Constraint *theSP;
  if (xFix != 0) {
    theSP = new SP_Constraint(numSPs++, nodeId, 0, 0.0);
    if (theSP == 0) {
      opserr << "WARNING TclPlaneTruss - fix - ran out of memory for SP_Constraint ";
      opserr << nodeId << endln;
      return TCL_ERROR;
    }
    if (theTclPlaneTrussDomain->addSP_Constraint(theSP) == false) {
      opserr << "WARNING TclPlaneTruss - fix - could not add SP_Constraint to domain ";
      opserr << nodeId << endln;
      return TCL_ERROR;
    }
  }

  if (yFix != 0) {
    theSP = new SP_Constraint(numSPs++, nodeId, 1, 0.0);
    if (theSP == 0) {
      opserr << "WARNING TclPlaneTruss - fix - ran out of memory for SP_Constraint ";
      opserr << nodeId << endln;
      return TCL_ERROR;
    }
    if (theTclPlaneTrussDomain->addSP_Constraint(theSP) == false) {
      opserr << "WARNING TclPlaneTruss - fix - could not add SP_Constraint to domain ";
      opserr << nodeId << endln;
      return TCL_ERROR;
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}







int
TclPlaneTruss_addNodalLd(ClientData clientData, Tcl_Interp *interp, int argc,   
			 char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 4) {
      interp->result = "WARNING bad command - want: load nodeId forceX forceY <tStart duration>";
      return TCL_ERROR;
  }    

  // get the id, x_loc and y_loc
  int nodeId;
  double xForce, yForce;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
      interp->result = "WARNING invalid nodeId - load nodeId forceX forceY ";
      return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &xForce) != TCL_OK) {
      interp->result = "WARNING invalid forceX - load nodeId forceX forceY ";
      return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &yForce) != TCL_OK) {
      interp->result = "WARNING invalid forceY - load nodeId forceX forceY ";
      return TCL_ERROR;
  }

  // now create the load and add it to the Domain
  Vector forces(2);
  forces(0)=xForce; forces(1)=yForce;

  bool isLoadConst = false;
  int loadPatternTag = 123456789; // some pattern that will never be used!

  // allow some additional options at end of command
  int endMarker = 4;
  while (endMarker != argc) {
    if (strcmp(argv[endMarker],"-const") == 0) {
      // allow user to specify const load
      isLoadConst = true;
    } else if (strcmp(argv[endMarker],"-pattern") == 0) {
      // allow user to specify load pattern other than current
      endMarker++;
      if (endMarker == argc || 
	  Tcl_GetInt(interp, argv[endMarker], &loadPatternTag) != TCL_OK) {

	opserr << "WARNING invalid patternTag - load " << nodeId << " ";
	opserr << "fX fY -pattern patternTag\n";
	return TCL_ERROR;
      }
    }
    endMarker++;
  }

  // get the current pattern tag if no tag given in i/p
  if (loadPatternTag == 123456789)
    if (theTclLoadPattern == 0) {
	opserr << "WARNING no current load pattern - load " << nodeId;
	opserr << " fX fY\n";
	return TCL_ERROR;
    } else	
      loadPatternTag = theTclLoadPattern->getTag();


  // get a tag for the load - use the number of nodal loads in the domain
  loadTag++;

  // create the load
  NodalLoad *theLoad = new NodalLoad(loadTag, nodeId, forces, isLoadConst);
  if (theLoad == 0) {
    opserr << "WARNING ran out of memory for load  - load " << nodeId;
    opserr << " fX fY\n";
    return TCL_ERROR;
  }

  // add the load to the domain
  if (theTclPlaneTrussDomain->addNodalLoad(theLoad, loadPatternTag) 
      == false) { 
      
    opserr << "WARNING TclModelBuilder - could not add load to domain ";
    delete theLoad;
    return TCL_ERROR;
  }  
  
  // if get here we have sucessfully created the node and added it 
  return TCL_OK;
}


int
TclPlaneTruss_done(ClientData clientData, Tcl_Interp *interp, int argc,   
		   char **argv)
{
  return TCL_OK;
}


