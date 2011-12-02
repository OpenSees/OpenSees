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
                                                                        
// $Revision: 1.3 $
// $Date: 2000-12-19 04:03:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclModelBuilder.cpp,v $
                                                                        
                                                                        
// File: ~/modelbuilder/tcl/TclModelBuilder.C
// 
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the class definition for TclModelBuilder.
// A TclModelBuilder adds the commands to create the model for the standard
// models that can be generated using the elements released with the g3 
// framework. currently these elements include:
//	1) linear-elastic 2 and 3d beam-column elements
//	2) non-linear material truss
//	3) non-linear 2 and 3d fiber-beam-column elements

//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <iostream.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <ArrayOfTaggedObjects.h>

#include <Domain.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <NodalLoad.h>
#include <LoadPattern.h>

#include <SectionForceDeformation.h>
#include <SectionRepres.h>

#include <CrdTransf2d.h>
#include <CrdTransf3d.h>

#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <TclModelBuilder.h>
#include <ImposedMotionSP.h>
#include <ImposedMotionSP1.h>
#include <MultiSupportPattern.h>

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static Domain *theTclDomain =0;
static TclModelBuilder *theTclBuilder =0;
extern LoadPattern *theTclLoadPattern;
extern MultiSupportPattern *theTclMultiSupportPattern;
static int eleArgStart = 0;
static int nodeLoadTag = 0;
// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int
TclModelBuilder_addNode(ClientData clientData, Tcl_Interp *interp, int argc, 
			char **argv);

int
TclModelBuilder_addElement(ClientData clientData, Tcl_Interp *interp,  int argc, 
			   char **argv);

int
TclModelBuilder_addUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
				    char **argv);

int
TclModelBuilder_addNDMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
			    char **argv);

int
TclModelBuilder_addSection(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv);
			    
int
TclModelBuilder_addPattern(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv);

int
TclModelBuilder_addSeries(ClientData clientData, Tcl_Interp *interp, int argc,   
			  char **argv);

int
TclModelBuilder_addHomogeneousBC(ClientData clientData, Tcl_Interp *interp, int argc,   
				 char **argv);

int
TclModelBuilder_addEqualDOF_MP (ClientData clientData, Tcl_Interp *interp,
                                int argc, char **argv);

int
TclModelBuilder_addMP(ClientData clientData, Tcl_Interp *interp, int argc,   
		      char **argv);

int
TclModelBuilder_addNodalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
			     char **argv);
int
TclModelBuilder_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc,   
			     char **argv);
int
TclModelBuilder_addSP(ClientData clientData, Tcl_Interp *interp, int argc,   
		      char **argv);

int
TclModelBuilder_addImposedMotionSP(ClientData clientData, 
				   Tcl_Interp *interp, 
				   int argc,    
				   char **argv);	

int
TclModelBuilder_addRemoPatch(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc,   
			     char **argv);  

int
TclModelBuilder_addRemoLayer(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc,   
			     char **argv);   
			       
int
TclModelBuilder_addRemoFiber(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc,    
			     char **argv);   

int
TclModelBuilder_addRemoGeomTransf(ClientData clientData, 
				  Tcl_Interp *interp, 
				  int argc,   
				  char **argv); 

				  

int
TclModelBuilder_addGroundMotion(ClientData clientData, 
				Tcl_Interp *interp, 
				int argc,    
				char **argv);

/// added by ZHY
int
TclModelBuilder_UpdateMaterialStage(ClientData clientData, 
				    Tcl_Interp *interp,  
				    int argc,
				    char **argv);
			   
// REMO
extern int
TclModelBuilder_addPatch (ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv,
			  TclModelBuilder *theTclBuilder);

			  
extern int
TclModelBuilder_addFiber (ClientData clientData, Tcl_Interp *interp,
			  int argc, char **argv,
			  TclModelBuilder *theTclBuilder);
			  

extern int
TclModelBuilder_addReinfLayer (ClientData clientData, Tcl_Interp *interp,
			       int argc, char **argv,
			       TclModelBuilder *theTclBuilder);


extern int
TclModelBuilder_addGeomTransf(ClientData, Tcl_Interp *, int, char **,
			      Domain*, TclModelBuilder *);   
	
					 
					 
//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

// constructor: the constructor will add certain commands to the interpreter
TclModelBuilder::TclModelBuilder(Domain &theDomain, Tcl_Interp *interp, int NDM, int NDF)
  :ModelBuilder(theDomain), ndm(NDM), ndf(NDF)
{
  theUniaxialMaterials = new ArrayOfTaggedObjects(32);
  theNDMaterials = new ArrayOfTaggedObjects(32);
  theSections  = new ArrayOfTaggedObjects(32);
  theSectionRepresents = new ArrayOfTaggedObjects(32);  
  the2dGeomTransfs = new ArrayOfTaggedObjects(32);  
  the3dGeomTransfs = new ArrayOfTaggedObjects(32);  

  // call Tcl_CreateCommand for class specific commands
  Tcl_CreateCommand(interp, "node", TclModelBuilder_addNode,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "element", TclModelBuilder_addElement,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "uniaxialMaterial", TclModelBuilder_addUniaxialMaterial,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "nDMaterial", TclModelBuilder_addNDMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "section", TclModelBuilder_addSection,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "pattern", TclModelBuilder_addPattern,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "load", TclModelBuilder_addNodalLoad,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "mass", TclModelBuilder_addNodalMass,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fix", TclModelBuilder_addHomogeneousBC,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "sp", TclModelBuilder_addSP,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "imposedSupportMotion", 
		    TclModelBuilder_addImposedMotionSP,
		    (ClientData)NULL, NULL);  
  
  Tcl_CreateCommand(interp, "groundMotion", 
		    TclModelBuilder_addGroundMotion,
		    (ClientData)NULL, NULL);    

  Tcl_CreateCommand(interp, "equalDOF", TclModelBuilder_addEqualDOF_MP,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "mp", TclModelBuilder_addMP,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "patch", TclModelBuilder_addRemoPatch,
		    (ClientData)NULL, NULL);  

  Tcl_CreateCommand(interp, "layer", TclModelBuilder_addRemoLayer,
		    (ClientData)NULL, NULL);    
  
  Tcl_CreateCommand(interp, "fiber", TclModelBuilder_addRemoFiber,
		    (ClientData)NULL, NULL);    

  Tcl_CreateCommand(interp, "geomTransf", TclModelBuilder_addRemoGeomTransf,
		    (ClientData)NULL, NULL);    

  
  ///new command for elast2plast in Multi-yield plasticity, by ZHY
  Tcl_CreateCommand(interp, "updateMaterialStage", 
		    TclModelBuilder_UpdateMaterialStage,
		    (ClientData)NULL, NULL);
  
  // set the static pointers in this file
  theTclBuilder = this;
  theTclDomain = &theDomain;
  theTclLoadPattern = 0;
  theTclMultiSupportPattern = 0;  

  nodeLoadTag = 0;
  eleArgStart = 0;
}

TclModelBuilder::~TclModelBuilder()
{
  theUniaxialMaterials->clearAll();
  theNDMaterials->clearAll();
  theSections->clearAll(); 
  theSectionRepresents->clearAll();
  the2dGeomTransfs->clearAll();
  the3dGeomTransfs->clearAll();

  // free up memory allocated in the constructor
  delete theUniaxialMaterials;
  delete theNDMaterials;
  delete theSections;
  delete theSectionRepresents;
  delete the2dGeomTransfs;
  delete the3dGeomTransfs;

  // set the pointers to 0 
  theTclDomain =0;
  theTclBuilder =0;
  theTclLoadPattern =0;
  theTclMultiSupportPattern = 0;  
  
  // may possibly invoke Tcl_DeleteCommand() later
}


//
// CLASS METHODS
//

int 
TclModelBuilder::buildFE_Model(void)
{
  // does nothing
  return 0;
}

int 
TclModelBuilder::getNDM(void) const
{
  return ndm;
}

int 
TclModelBuilder::getNDF(void) const
{
  return ndf;
}

int 
TclModelBuilder::addUniaxialMaterial(UniaxialMaterial &theMaterial)
{
  bool result = theUniaxialMaterials->addComponent(&theMaterial);
  if (result == true)
    return 0;
  else {
    cerr << "TclModelBuilder::addUniaxialMaterial() - failed to add material: " << theMaterial;
    return -1;
  }
}


UniaxialMaterial *
TclModelBuilder::getUniaxialMaterial(int tag)
{
  TaggedObject *mc = theUniaxialMaterials->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  UniaxialMaterial *result = (UniaxialMaterial *)mc;
  return result;
}

int 
TclModelBuilder::addNDMaterial(NDMaterial &theMaterial)
{
  bool result = theNDMaterials->addComponent(&theMaterial);
  if (result == true)
    return 0;
  else {
    cerr << "TclModelBuilder::addNDMaterial() - failed to add material: " << theMaterial;
    return -1;
  }
}


NDMaterial *
TclModelBuilder::getNDMaterial(int tag)
{
  TaggedObject *mc = theNDMaterials->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // otherweise we do a cast and return
  NDMaterial *result = (NDMaterial *)mc;
  return result;
}

int 
TclModelBuilder::addSection(SectionForceDeformation &theSection)
{
  bool result = theSections->addComponent(&theSection);
  if (result == true)
    return 0;
  else {
    cerr << "TclModelBuilder::addSection() - failed to add section: " << theSection;
    return -1;
  }
}



SectionForceDeformation *
TclModelBuilder::getSection(int tag)
{
  TaggedObject *mc = theSections->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // do a cast and return
  SectionForceDeformation *result = (SectionForceDeformation *)mc;
  return result;
}



int 
TclModelBuilder::addSectionRepres(SectionRepres &theSectionRepres)
{
  bool result = theSectionRepresents->addComponent(&theSectionRepres);

  if (result == true)
    return 0;
  else {
      cerr << "TclModelBuilder::addSectionRepres() - failed to add SectionRepres\n";
      return -1;
  }
}


SectionRepres *
TclModelBuilder::getSectionRepres(int tag)
{
  TaggedObject *mc = theSectionRepresents->getComponentPtr(tag);
  if (mc == 0) return 0;
  SectionRepres *result = (SectionRepres *)mc;
  return result;
}



int 
TclModelBuilder::addCrdTransf2d(CrdTransf2d &theCrdTransf)
{
  bool result = the2dGeomTransfs->addComponent(&theCrdTransf);
  if (result == true)
    return 0;
  else {
    cerr << "TclModelBuilder::addCrdTransf() - failed to add crdTransf: " << theCrdTransf;
    return -1;
  }
}


int 
TclModelBuilder::addCrdTransf3d(CrdTransf3d &theCrdTransf)
{
  bool result = the3dGeomTransfs->addComponent(&theCrdTransf);
  if (result == true)
    return 0;
  else {
    cerr << "TclModelBuilder::addCrdTransf() - failed to add crdTransf: " << theCrdTransf;
    return -1;
  }
}



CrdTransf2d *
TclModelBuilder::getCrdTransf2d(int tag)
{
  TaggedObject *mc = the2dGeomTransfs->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // do a cast and return
  CrdTransf2d *result = (CrdTransf2d *)mc;
  return result;
}


CrdTransf3d *
TclModelBuilder::getCrdTransf3d(int tag)
{
  TaggedObject *mc = the3dGeomTransfs->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // do a cast and return
  CrdTransf3d *result = (CrdTransf3d *)mc;
  return result;
}

//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//

void printCommand(int argc, char **argv)
{
  cerr << "Input command: ";
  for (int i=0; i<argc; i++)
    cerr << argv[i] << " ";
  cerr << endl;
} 

int
TclModelBuilder_addNode(ClientData clientData, Tcl_Interp *interp, int argc, 
                        char **argv)
{

  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed" << endl;
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  // make sure corect number of arguments on command line
  if (argc < 2+ndm) {
    cerr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    cerr << "Want: node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }    

  Node *theNode = 0;

  // get the nodal id
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    cerr << "WARNING invalid nodeTag\n";
    cerr << "Want: node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }

  // read in the coordinates and create the node
  double xLoc, yLoc, zLoc;
  if (ndm == 1) { 
    // create a node in 1d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      cerr << "WARNING invalid XCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    theNode = new Node(nodeId,ndf,xLoc);
  } 

  else if (ndm == 2) { 
    // create a node in 2d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      cerr << "WARNING invalid XCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      cerr << "WARNING invalid YCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    theNode = new Node(nodeId,ndf,xLoc,yLoc);
  } 

  else if (ndm == 3) { 
    // create a node in 3d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      cerr << "WARNING invalid XCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      cerr << "WARNING invalid YCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &zLoc) != TCL_OK) {
      cerr << "WARNING invalid ZCoordinate\n";
      cerr << "node: " << nodeId << endl;
      return TCL_ERROR;
    }
    theNode = new Node(nodeId,ndf,xLoc,yLoc,zLoc);
  } else {
      cerr << "WARNING invalid ndm\n";
      cerr << "node: " << nodeId << endl;;
      return TCL_ERROR;
  }

  if (theNode == 0) {
    cerr << "WARNING ran out of memory creating node\n";
    cerr << "node: " << nodeId << endl;
    return TCL_ERROR;
  }

  if (theTclDomain->addNode(theNode) == false) {
    cerr << "WARNING failed to add node to the domain\n";
    cerr << "node: " << nodeId << endl;
    delete theNode; // otherwise memory leak
    return TCL_ERROR;
  }

  // check for mass terms
  if (argc > 2+ndm) {
    if (strcmp(argv[2+ndm],"-mass") == 0) {
      if (argc < 3+ndm+ndf) {
        cerr << "WARNING incorrect number of nodal mass terms\n";
        cerr << "node: " << nodeId << endl;
        return TCL_ERROR;      
      }	
      Matrix mass(ndf,ndf);
      double theMass;
      for (int i=0; i<ndf; i++) {
	if (Tcl_GetDouble(interp, argv[i+3+ndm], &theMass) != TCL_OK) {
	  cerr << "WARNING invalid nodal mass term\n";
	  cerr << "node: " << nodeId << ", dof: " << i+1 << endl;
	  return TCL_ERROR;
	}
	mass(i,i) = theMass;
      }
      theNode->setMass(mass);      
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}





// the function for creating ne material objects and patterns is in a seperate file.
// this allows new material and patternobjects to be added without touching this file.
// does so at the expense of an extra procedure call.


extern int 
TclModelBuilderElementCommand(ClientData clientData, 
			      Tcl_Interp *interp, int argc,    
			      char **argv, 
			      Domain *theDomain, TclModelBuilder *theTclBuilder);
int
TclModelBuilder_addElement(ClientData clientData, Tcl_Interp *interp, 
			   int argc,    char **argv)
                          
{
  return TclModelBuilderElementCommand(clientData, interp, 
				       argc, argv, theTclDomain, theTclBuilder);
}


extern int
TclModelBuilderUniaxialMaterialCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				 char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, 
				    int argc, char **argv)
                          
{
  return TclModelBuilderUniaxialMaterialCommand(clientData, interp, 
						argc, argv, theTclBuilder);
}

extern int
TclModelBuilderNDMaterialCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				  char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addNDMaterial(ClientData clientData, Tcl_Interp *interp, 
			    int argc,    char **argv)
                          
{
  return TclModelBuilderNDMaterialCommand(clientData, interp, 
					  argc, argv, theTclBuilder);
}

extern int
TclModelBuilderSectionCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				  char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addSection(ClientData clientData, Tcl_Interp *interp, 
			    int argc,    char **argv)
                          
{
  return TclModelBuilderSectionCommand(clientData, interp, 
				       argc, argv, theTclBuilder);
}



extern int
TclPatternCommand(ClientData clientData, Tcl_Interp *interp, 
			   int argc, char **argv, Domain *theDomain);
			   
int
TclModelBuilder_addPattern(ClientData clientData, Tcl_Interp *interp, 
			   int argc, char **argv)
			  
{
  return TclPatternCommand(clientData, interp, argc, argv, theTclDomain);
}




extern int
TclGroundMotionCommand(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc,    
		       char **argv,
		       MultiSupportPattern *thePattern);

int
TclModelBuilder_addGroundMotion(ClientData clientData, Tcl_Interp *interp, 
			   int argc, char **argv)
			  
{
  return TclGroundMotionCommand(clientData, interp, argc, argv, 
				theTclMultiSupportPattern);
}


int
TclModelBuilder_addNodalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
			 char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - load \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  NodalLoad *theLoad = 0;
  
  // make sure at least one other argument to contain type of system
  if (argc < (2 + ndf)) {
    cerr << "WARNING bad command - want: load nodeId " << ndf << " forces\n";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    cerr << "WARNING invalid nodeId: " << argv[1];
    cerr << " - load nodeId " << ndf << " forces\n";
    return TCL_ERROR;
  }

  // get the load vector
  Vector forces(ndf);
  for (int i=0; i<ndf; i++) {
    double theForce;
    if (Tcl_GetDouble(interp, argv[2+i], &theForce) != TCL_OK) {
      cerr << "WARNING invalid force " << i+1 << " - load " << nodeId;
      cerr << " " << ndf << " forces\n";
      return TCL_ERROR;
    } else
      forces(i) = theForce;
  }

  bool isLoadConst = false;
  int loadPatternTag = 123456789; // some pattern that will never be used!

  // allow some additional options at end of command
  int endMarker = 2+ndf;
  while (endMarker != argc) {
    if (strcmp(argv[endMarker],"-const") == 0) {
      // allow user to specify const load
      isLoadConst = true;
    } else if (strcmp(argv[endMarker],"-pattern") == 0) {
      // allow user to specify load pattern other than current
      endMarker++;
      if (endMarker == argc || 
	  Tcl_GetInt(interp, argv[endMarker], &loadPatternTag) != TCL_OK) {

	cerr << "WARNING invalid patternTag - load " << nodeId << " ";
	cerr << ndf << " forces pattern patterntag\n";
	return TCL_ERROR;
      }
    }
    endMarker++;
  }

  // get the current pattern tag if no tag given in i/p
  if (loadPatternTag == 123456789)
    if (theTclLoadPattern == 0) {
	cerr << "WARNING no current load pattern - load " << nodeId;
	cerr << " " << ndf << " forces\n";
	return TCL_ERROR;
    } else 
	loadPatternTag = theTclLoadPattern->getTag();

  // create the load
  theLoad = new NodalLoad(nodeLoadTag, nodeId, forces, isLoadConst);
  if (theLoad == 0) {
    cerr << "WARNING ran out of memory for load  - load " << nodeId;
    cerr << " " << ndf << " forces\n";
    return TCL_ERROR;
  }

  // add the load to the domain
  if (theTclDomain->addNodalLoad(theLoad, loadPatternTag) == false) {
    cerr << "WARNING TclModelBuilder - could not add load to domain ";
    printCommand(argc, argv);
    delete theLoad;
    return TCL_ERROR;
  }
  nodeLoadTag++;

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclModelBuilder_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc, 
                        char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - load \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // make sure at least one other argument to contain type of system
  if (argc < (2 + ndf)) {
    cerr << "WARNING bad command - want: mass nodeId " << ndf << " mass values\n";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    cerr << "WARNING invalid nodeId: " << argv[1];
    cerr << " - mass nodeId " << ndf << " forces\n";
    return TCL_ERROR;
  }

  Node *theNode = 0;
  theNode = theTclDomain->getNode(nodeId);

  if (theNode == 0)
  {
    cerr << "WARNING failed to get node pointer from the domain\n";
    cerr << "node: " << nodeId << endl;
    return TCL_ERROR;
  }

  // check for mass terms
  Matrix mass(ndf,ndf);
  double theMass;
  for (int i=0; i<ndf; i++) 
  {
     if (Tcl_GetDouble(interp, argv[i+2], &theMass) != TCL_OK) 
     {
	  cerr << "WARNING invalid nodal mass term\n";
	  cerr << "node: " << nodeId << ", dof: " << i+1 << endl;
	  return TCL_ERROR;
      }
      mass(i,i) = theMass;
  }
      
  theNode->setMass(mass);      

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}





int
TclModelBuilder_addHomogeneousBC(ClientData clientData, Tcl_Interp *interp, int argc,   
				 char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();
  int numSPs = theTclDomain->getNumSPs();

  // check number of arguments
  if (argc < (2 + ndf)) {
    cerr << "WARNING bad command - want: fix nodeId " << ndf << " [0,1] conditions";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
      cerr << "WARNING invalid nodeId - fix nodeId " << ndf << " [0,1] conditions\n";
      return TCL_ERROR;
  }

  // get the fixity condition and add the constraint if fixed
  for (int i=0; i<ndf; i++) {
    int theFixity;
    if (Tcl_GetInt(interp, argv[2+i], &theFixity) != TCL_OK) {
      cerr << "WARNING invalid fixity " << i+1 << " - load " << nodeId;
      cerr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    } else {
      if (theFixity != 0) {
	// create a homogeneous constraint
	SP_Constraint *theSP = new SP_Constraint(numSPs, nodeId, i, 0.0);
	if (theSP == 0) {
	  cerr << "WARNING ran out of memory for SP_Constraint ";
	  cerr << "fix " << nodeId << " " << ndf << " [0,1] conditions\n";
	  return TCL_ERROR;
	}
	if (theTclDomain->addSP_Constraint(theSP) == false) {
	  cerr << "WARNING could not add SP_Constraint to domain - fix";
	  cerr << nodeId << " " << ndf << " [0,1] conditions\n";
	  delete theSP;
	  return TCL_ERROR;
	}
	numSPs++;      }
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}


int
TclModelBuilder_addSP(ClientData clientData, Tcl_Interp *interp, int argc,   
		      char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - sp \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // check number of arguments
  if (argc < 4) {
    cerr << "WARNING bad command - want: sp nodeId dofID value";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the nodeID, dofId and value of the constraint
  int nodeId, dofId;
  double value;

  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    cerr << "WARNING invalid nodeId: " << argv[1] << " -  sp nodeId dofID value\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dofId) != TCL_OK) {
    cerr << "WARNING invalid dofId: " << argv[2] << " -  sp ";
    cerr << nodeId << " dofID value\n";
      return TCL_ERROR;
  }
  dofId--; // DECREMENT THE DOF VALUE BY 1 TO GO TO OUR C++ INDEXING

  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    cerr << "WARNING invalid value: " << argv[3] << " -  sp ";
    cerr << nodeId << " dofID value\n";
      return TCL_ERROR;
  }

  bool isSpConst = false;
  int loadPatternTag = 123456789; // some pattern that will never be used!

  // allow some additional options at end of command
  theTclLoadPattern->getTag();
  int endMarker = 4;
  while (endMarker != argc) {
    if (strcmp(argv[endMarker],"-const") == 0) {
      // allow user to specify const load
      isSpConst = true;
    } else if (strcmp(argv[endMarker],"-pattern") == 0) {
      // allow user to specify load pattern other than current
      endMarker++;
      if (endMarker == argc || 
	  Tcl_GetInt(interp, argv[endMarker], &loadPatternTag) != TCL_OK) {

	cerr << "WARNING invalid patternTag - load " << nodeId << " ";
	cerr << ndf << " forces pattern patterntag\n";
	return TCL_ERROR;
      }
    }  
    endMarker++;
  }

  // if load pattern tag has not changed - get the pattern tag from current one
  if (loadPatternTag == 123456789) {
    if (theTclLoadPattern == 0) {
      cerr << "WARNING no current pattern - sp " << nodeId << " dofID value\n";
      return TCL_ERROR;
    } else	
      loadPatternTag = theTclLoadPattern->getTag();
  }
  
  LoadPattern *thePattern = theTclDomain->getLoadPattern(loadPatternTag);
  SP_ConstraintIter &theSPs = thePattern->getSPs();
  int numSPs = 0;
  SP_Constraint *theSP2;
  while ((theSP2 = theSPs()) != 0)
      numSPs++;
  
  
  // create a homogeneous constraint
  SP_Constraint *theSP = new SP_Constraint(numSPs, nodeId, dofId, value, isSpConst);

  if (theSP == 0) {
    cerr << "WARNING ran out of memory for SP_Constraint ";
    cerr << " - sp " << nodeId << " dofID value\n";
    return TCL_ERROR;
  }
  if (theTclDomain->addSP_Constraint(theSP, loadPatternTag) == false) {
    cerr << "WARNING could not add SP_Constraint to domain ";
    printCommand(argc, argv);
    delete theSP;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclModelBuilder_addImposedMotionSP(ClientData clientData, 
				   Tcl_Interp *interp, 
				   int argc,   
				   char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    cerr << "WARNING builder has been destroyed - sp \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // check number of arguments
  if (argc < 4) {
    cerr << "WARNING bad command - want: imposedSupportMotion nodeId dofID gMotionID\n";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the nodeID, dofId and value of the constraint
  int nodeId, dofId, gMotionID;

  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    cerr << "WARNING invalid nodeId: " << argv[1];
    cerr << " - imposedSupportMotion nodeId dofID gMotionID\n";    
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dofId) != TCL_OK) {
    cerr << "WARNING invalid dofId: " << argv[2] << " -  imposedSupportMotion ";
    cerr << nodeId << " dofID gMotionID\n";    
      return TCL_ERROR;
  }
  dofId--; // DECREMENT THE DOF VALUE BY 1 TO GO TO OUR C++ INDEXING

  if (Tcl_GetInt(interp, argv[3], &gMotionID) != TCL_OK) {
    cerr << "WARNING invalid gMotionID: " << argv[3] << " -  imposedSupportMotion ";
    cerr << nodeId << " dofID gMotionID\n";
    return TCL_ERROR;
  }

  bool alt = false;
  if (argc == 5) {
    if (strcmp(argv[4],"-other") == 0) 
      alt = true;
  }

  MultiSupportPattern *thePattern = theTclMultiSupportPattern;
  int loadPatternTag = thePattern->getTag();
  
  GroundMotion *theGMotion = thePattern->getMotion(gMotionID);
  if (theGMotion == 0) {
    cerr << "WARNING no GroundMotion with tag: " << argv[3];
    cerr << " in current MultipleSupportPattern\n";
    return TCL_ERROR;      
  }
  
  SP_ConstraintIter &theSPs = thePattern->getSPs();
  int numSPs = 0;
  SP_Constraint *theSP2;
  while ((theSP2 = theSPs()) != 0)
      numSPs++;
  
  // create a new ImposedMotionSP
  SP_Constraint *theSP;
  if (alt == true) {
    theSP = new ImposedMotionSP1(numSPs, nodeId, dofId, 
				*theGMotion, false);
  }
  else {
    theSP = new ImposedMotionSP(numSPs, nodeId, dofId, 
				*theGMotion, false);
  }
  if (theSP == 0) {
    cerr << "WARNING ran out of memory for ImposedMotionSP ";
    cerr << " -  imposedSupportMotion ";
    cerr << nodeId << " " << dofId++ << " " << gMotionID << endl;
    return TCL_ERROR;
  }
  if (thePattern->addSP_Constraint(theSP) == false) {
    cerr << "WARNING could not add SP_Constraint to pattern ";
    printCommand(argc, argv);
    delete theSP;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}






int
TclModelBuilder_addEqualDOF_MP (ClientData clientData, Tcl_Interp *interp,
                                int argc, char **argv)
{
        // Ensure the destructor has not been called
        if (theTclBuilder == 0)
        {
                cerr << "WARNING builder has been destroyed - equalDOF \n";
                return TCL_ERROR;
        }

        // Check number of arguments
        if (argc < 4)
        {
                cerr << "WARNING bad command - want: equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
                printCommand (argc, argv);
                return TCL_ERROR;
        }

        // Read in the node IDs and the DOF
        int RnodeID, CnodeID, dofID;

        if (Tcl_GetInt (interp, argv[1], &RnodeID) != TCL_OK)
        {
                cerr << "WARNING invalid RnodeID: " << argv[1]
                         << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
                return TCL_ERROR;
        }
        if (Tcl_GetInt (interp, argv[2], &CnodeID) != TCL_OK)
        {
                cerr << "WARNING invalid CnodeID: " << argv[2]
                         << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
                return TCL_ERROR;
        }

        // The number of DOF to be coupled
        int numDOF = argc - 3;

        // The constraint matrix ... U_c = C_cr * U_r
        Matrix Ccr (numDOF, numDOF);
        Ccr.Zero();

        // The vector containing the retained and constrained DOFs
        ID rcDOF (numDOF);

        int i, j;
        // Read the degrees of freedom which are to be coupled
        for (i = 3, j = 0; i < argc; i++, j++)
        {
                if (Tcl_GetInt (interp, argv[i], &dofID) != TCL_OK)
                {
                        cerr << "WARNING invalid dofID: " << argv[3]
                                 << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
                        return TCL_ERROR;
                }

                rcDOF (j) = --dofID;    // Decrement for C++ indexing
                Ccr (j,j) = 1;
        }

        // Use this for the MP tag
        int numMPs = theTclDomain->getNumMPs();

        // Create the multi-point constraint
        MP_Constraint *theMP = new MP_Constraint (numMPs, RnodeID, CnodeID, Ccr, rcDOF, rcDOF);
        if (theMP == 0)
        {
                cerr << "WARNING ran out of memory for equalDOF MP_Constraint ";
                printCommand (argc, argv);
                return TCL_ERROR;
        }

        // Add the multi-point constraint to the domain
        if (theTclDomain->addMP_Constraint (theMP) == false)
        {
                cerr << "WARNING could not add equalDOF MP_Constraint to domain ";
                printCommand(argc, argv);
                delete theMP;
                return TCL_ERROR;
        }

        return TCL_OK;
}


int
TclModelBuilder_addMP(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv)
{
  cerr << "WARNING - TclModelBuilder_addMP() not yet implemented\n";
  return TCL_OK;
}


int
TclModelBuilder_addRemoPatch(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv)
{
  return TclModelBuilder_addPatch(clientData, interp, argc,argv,
				    theTclBuilder);
}

int
TclModelBuilder_addRemoFiber(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv)
{
  return TclModelBuilder_addFiber(clientData, interp, argc,argv,
				  theTclBuilder);
}

int
TclModelBuilder_addRemoLayer(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv)
{
  return TclModelBuilder_addReinfLayer(clientData, interp, argc,argv,
				       theTclBuilder);
}


					 
int
TclModelBuilder_addRemoGeomTransf(ClientData clientData, Tcl_Interp *interp, int argc,   
			   char **argv)
{
  return TclModelBuilder_addGeomTransf(clientData, interp, argc,argv,
				       theTclDomain,
				       theTclBuilder);
}


/// added by ZHY
extern int 
TclModelBuilderUpdateMaterialStageCommand(ClientData clientData, 
					  Tcl_Interp *interp, 
					  int argc, 
					  char **argv, 
					  TclModelBuilder *theTclBuilder);
int
TclModelBuilder_UpdateMaterialStage(ClientData clientData, 
				    Tcl_Interp *interp,  
				    int argc, 
				    char **argv)
{
  return TclModelBuilderUpdateMaterialStageCommand(clientData, interp, 
				       argc, argv, theTclBuilder);
}
