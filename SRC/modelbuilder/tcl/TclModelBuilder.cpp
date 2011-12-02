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
                                                                        
// $Revision: 1.27 $
// $Date: 2005-01-03 22:32:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/TclModelBuilder.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 07/99
//
// Description: This file contains the class definition for TclModelBuilder.
// A TclModelBuilder adds the commands to create the model for the standard
// models that can be generated using the elements released with the g3 
// framework. currently these elements include:
//	1) linear-elastic 2 and 3d beam-column elements
//	2) non-linear material truss
//	3) non-linear 2 and 3d fiber-beam-column elements

//
// What: "@(#) TclModelBuilder.cpp, revA"

#include <stdlib.h>
#include <string.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <ArrayOfTaggedObjects.h>

#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <NodalLoad.h>
#include <Beam2dPointLoad.h>
#include <Beam2dUniformLoad.h>
#include <Beam2dTempLoad.h>
#include <Beam3dPointLoad.h>
#include <Beam3dUniformLoad.h>
#include <BrickSelfWeight.h>
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

#include <Block2D.h>
#include <Block3D.h>
// Added by Scott J. Brandenberg (sjbrandenberg@ucdavis.edu)
#include <PySimple1Gen.h>
#include <TzSimple1Gen.h>
// End added by SJB

#include <YieldSurface_BC.h>
#include <YS_Evolution.h>
#include <PlasticHardeningMaterial.h>
#include <CyclicModel.h> //!!
#include <DamageModel.h> //!!

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

static Domain *theTclDomain =0;
static TclModelBuilder *theTclBuilder =0;
extern LoadPattern *theTclLoadPattern;
extern MultiSupportPattern *theTclMultiSupportPattern;
static int eleArgStart = 0;
static int nodeLoadTag = 0;
static int eleLoadTag = 0;

static int currentSpTag = 0;

// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int
TclModelBuilder_addNode(ClientData clientData, Tcl_Interp *interp, int argc, 
			TCL_Char **argv);

int
TclModelBuilder_addElement(ClientData clientData, Tcl_Interp *interp,  int argc, 
			   TCL_Char **argv);

int
TclModelBuilder_addUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
				    TCL_Char **argv);

int
TclModelBuilder_addNDMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
			    TCL_Char **argv);

int
TclModelBuilder_addSection(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv);

int
TclModelBuilder_addYieldSurface_BC(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv);

int
TclModelBuilder_addYS_EvolutionModel(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv);

int
TclModelBuilder_addYS_PlasticMaterial(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv);

int //!!
TclModelBuilder_addCyclicModel(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv);			    
int //!!
TclModelBuilder_addDamageModel(ClientData clientData, Tcl_Interp *interp,
			       int argc, TCL_Char **argv);	

int
TclModelBuilder_addPattern(ClientData clientData, Tcl_Interp *interp, int argc,
			   TCL_Char **argv);

int
TclModelBuilder_addPattern(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv);

int
TclModelBuilder_addSeries(ClientData clientData, Tcl_Interp *interp, int argc,   
			  TCL_Char **argv);

int
TclModelBuilder_addHomogeneousBC(ClientData clientData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv);
int
TclModelBuilder_addHomogeneousBC_X(ClientData clientData, Tcl_Interp *interp, int argc,
				   TCL_Char **argv);
int
TclModelBuilder_addHomogeneousBC_Y(ClientData clientData, Tcl_Interp *interp, int argc,
				   TCL_Char **argv);
int
TclModelBuilder_addHomogeneousBC_Z(ClientData clientData, Tcl_Interp *interp, int argc,
				   TCL_Char **argv);
int
TclModelBuilder_addEqualDOF_MP (ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char **argv);

int
TclModelBuilder_addMP(ClientData clientData, Tcl_Interp *interp, int argc,   
		      TCL_Char **argv);

int
TclModelBuilder_addNodalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
			     TCL_Char **argv);

int
TclModelBuilder_addElementalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
				 TCL_Char **argv);

int
TclModelBuilder_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc,   
			     TCL_Char **argv);
int
TclModelBuilder_addSP(ClientData clientData, Tcl_Interp *interp, int argc,   
		      TCL_Char **argv);

int
TclModelBuilder_addImposedMotionSP(ClientData clientData, 
				   Tcl_Interp *interp, 
				   int argc,    
				   TCL_Char **argv);	
// Added by Scott J. Brandenberg
int
TclModelBuilder_doPySimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char **argv);

int
TclModelBuilder_doTzSimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char **argv);
// End added by SJB

int
TclModelBuilder_doBlock2D(ClientData clientData, Tcl_Interp *interp, int argc, 
			  TCL_Char **argv);

int
TclModelBuilder_doBlock3D(ClientData clientData, Tcl_Interp *interp, int argc, 
			  TCL_Char **argv);

int
TclModelBuilder_addRemoPatch(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc,   
			     TCL_Char **argv);  

int
TclModelBuilder_addRemoLayer(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc,   
			     TCL_Char **argv);   
			       
int
TclModelBuilder_addRemoFiber(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc,    
			     TCL_Char **argv);   

int
TclModelBuilder_addRemoGeomTransf(ClientData clientData, 
				  Tcl_Interp *interp, 
				  int argc,   
				  TCL_Char **argv); 

int
TclModelBuilder_addGroundMotion(ClientData clientData, 
				Tcl_Interp *interp, 
				int argc,    
				TCL_Char **argv);

/// added by ZHY
int
TclModelBuilder_UpdateMaterialStage(ClientData clientData, 
				    Tcl_Interp *interp,  
				    int argc,
				    TCL_Char **argv);
			   

/// added by ZHY			   
int
TclModelBuilder_UpdateParameter(ClientData clientData, 
				    Tcl_Interp *interp,  
				    int argc, 
				    TCL_Char **argv);

// REMO
extern int
TclModelBuilder_addPatch (ClientData clientData, Tcl_Interp *interp,
			  int argc, TCL_Char **argv,
			  TclModelBuilder *theTclBuilder);

			  
extern int
TclModelBuilder_addFiber (ClientData clientData, Tcl_Interp *interp,
			  int argc, TCL_Char **argv,
			  TclModelBuilder *theTclBuilder);
			  

extern int
TclModelBuilder_addReinfLayer (ClientData clientData, Tcl_Interp *interp,
			       int argc, TCL_Char **argv,
			       TclModelBuilder *theTclBuilder);


extern int
TclModelBuilder_addGeomTransf(ClientData, Tcl_Interp *, int, TCL_Char **,
			      Domain*, TclModelBuilder *);   
	
//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

// constructor: the constructor will add certain commands to the interpreter
TclModelBuilder::TclModelBuilder(Domain &theDomain, Tcl_Interp *interp, int NDM, int NDF)
  :ModelBuilder(theDomain), ndm(NDM), ndf(NDF), theInterp(interp)
{
  theUniaxialMaterials = new ArrayOfTaggedObjects(32);
  theNDMaterials = new ArrayOfTaggedObjects(32);
  theSections  = new ArrayOfTaggedObjects(32);
  theSectionRepresents = new ArrayOfTaggedObjects(32);  
  the2dGeomTransfs = new ArrayOfTaggedObjects(32);  
  the3dGeomTransfs = new ArrayOfTaggedObjects(32);  
  theYieldSurface_BCs = new ArrayOfTaggedObjects(32);
  theCycModels = new ArrayOfTaggedObjects(32); //!!
  theDamageModels = new ArrayOfTaggedObjects(32); //!!
  theYS_EvolutionModels = new ArrayOfTaggedObjects(32);
  thePlasticMaterials = new ArrayOfTaggedObjects(32);

  SP_ConstraintIter &theSPs = theDomain.getSPs();
  SP_Constraint *theSP;
  while ((theSP = theSPs()) != 0) {
    int spTag = theSP->getTag();
    if (spTag >= currentSpTag)
      currentSpTag = spTag+1;
  }    

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

  Tcl_CreateCommand(interp, "yieldSurface_BC", TclModelBuilder_addYieldSurface_BC,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "ysEvolutionModel", TclModelBuilder_addYS_EvolutionModel,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "plasticMaterial", TclModelBuilder_addYS_PlasticMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "cyclicModel", TclModelBuilder_addCyclicModel,
		    (ClientData)NULL, NULL); //!!

  Tcl_CreateCommand(interp, "damageModel", TclModelBuilder_addDamageModel,
		    (ClientData)NULL, NULL); //!!

  Tcl_CreateCommand(interp, "pattern", TclModelBuilder_addPattern,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "load", TclModelBuilder_addNodalLoad,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "eleLoad", TclModelBuilder_addElementalLoad,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "mass", TclModelBuilder_addNodalMass,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fix", TclModelBuilder_addHomogeneousBC,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fixX", TclModelBuilder_addHomogeneousBC_X,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fixY", TclModelBuilder_addHomogeneousBC_Y,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fixZ", TclModelBuilder_addHomogeneousBC_Z,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "sp", TclModelBuilder_addSP,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "imposedMotion", 
		    TclModelBuilder_addImposedMotionSP,
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

  // Added by Scott J. Brandenberg
  Tcl_CreateCommand(interp, "PySimple1Gen", TclModelBuilder_doPySimple1Gen,
                        (ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "TzSimple1Gen", TclModelBuilder_doTzSimple1Gen,
                        (ClientData)NULL, NULL);
  // End added by SJB

  Tcl_CreateCommand(interp, "block2D", TclModelBuilder_doBlock2D,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "block3D", TclModelBuilder_doBlock3D,
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
  
  ///new command for updating properties of soil materials, by ZHY
  Tcl_CreateCommand(interp, "updateParameter", 
		    TclModelBuilder_UpdateParameter,
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
  currentSpTag = 0;

  theUniaxialMaterials->clearAll();
  theNDMaterials->clearAll();
  theSections->clearAll(); 
  theSectionRepresents->clearAll();
  the2dGeomTransfs->clearAll();
  the3dGeomTransfs->clearAll();
  theYieldSurface_BCs->clearAll();
  theYS_EvolutionModels->clearAll();
  thePlasticMaterials->clearAll();
  theCycModels->clearAll();//!!
  theDamageModels->clearAll();//!!

  // free up memory allocated in the constructor
  delete theUniaxialMaterials;
  delete theNDMaterials;
  delete theSections;
  delete theSectionRepresents;
  delete the2dGeomTransfs;
  delete the3dGeomTransfs;
  delete theYieldSurface_BCs;
  delete theYS_EvolutionModels;
  delete thePlasticMaterials;
  delete theCycModels;//!!
  delete theDamageModels;//!!

  // set the pointers to 0 
  theTclDomain =0;
  theTclBuilder =0;
  theTclLoadPattern =0;
  theTclMultiSupportPattern = 0;  
  
  // may possibly invoke Tcl_DeleteCommand() later
  Tcl_DeleteCommand(theInterp, "node");
  Tcl_DeleteCommand(theInterp, "element");
  Tcl_DeleteCommand(theInterp, "uniaxialMaterial");
  Tcl_DeleteCommand(theInterp, "nDMaterial");
  Tcl_DeleteCommand(theInterp, "section");
  Tcl_DeleteCommand(theInterp, "pattern");
  Tcl_DeleteCommand(theInterp, "load");
  Tcl_DeleteCommand(theInterp, "mass");
  Tcl_DeleteCommand(theInterp, "fix");
  Tcl_DeleteCommand(theInterp, "fixX");
  Tcl_DeleteCommand(theInterp, "fixY");
  Tcl_DeleteCommand(theInterp, "fixZ");
  Tcl_DeleteCommand(theInterp, "sp");
  Tcl_DeleteCommand(theInterp, "imposedSupportMotion");
  Tcl_DeleteCommand(theInterp, "groundMotion");
  Tcl_DeleteCommand(theInterp, "equalDOF");
  Tcl_DeleteCommand(theInterp, "mp");
  Tcl_DeleteCommand(theInterp, "PySimple1Gen");  // Added by Scott J. Brandenberg
  Tcl_DeleteCommand(theInterp, "TzSimple1Gen");  // Added by Scott J. Brandenberg
  Tcl_DeleteCommand(theInterp, "block2D");
  Tcl_DeleteCommand(theInterp, "block3D");
  Tcl_DeleteCommand(theInterp, "patch");
  Tcl_DeleteCommand(theInterp, "layer");
  Tcl_DeleteCommand(theInterp, "fiber");
  Tcl_DeleteCommand(theInterp, "geomTransf");
  Tcl_DeleteCommand(theInterp, "updateMaterialStage");
  Tcl_DeleteCommand(theInterp, "updateParameter");
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
    opserr << "TclModelBuilder::addUniaxialMaterial() - failed to add material: " << theMaterial;
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
    opserr << "TclModelBuilder::addNDMaterial() - failed to add material: " << theMaterial;
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
    opserr << "TclModelBuilder::addSection() - failed to add section: " << theSection;
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
TclModelBuilder::addYS_EvolutionModel(YS_Evolution &theModel)
{
  bool result = theYS_EvolutionModels->addComponent(&theModel);
  if (result == true)
    return 0;
  else {
    opserr << "TclModelBuilder::addYS_EvolutionModel() - failed to add model " << theModel;
    return -1;
  }
}


YS_Evolution *
TclModelBuilder::getYS_EvolutionModel(int tag)
{
  TaggedObject *mc = theYS_EvolutionModels->getComponentPtr(tag);
  if (mc == 0)
    return 0;

  // otherweise we do a cast and return
  YS_Evolution *result = (YS_Evolution *)mc;
  return result;
}

int
TclModelBuilder::addYieldSurface_BC(YieldSurface_BC &theYS)
{
//	TaggedObject *mc = &theYS;

  bool result = theYieldSurface_BCs->addComponent(&theYS);
  if (result == true)
    return 0;
  else {
    opserr << "TclModelBuilder::addYieldSurfaceBC() - failed to add YS: " << theYS;
    return -1;
  }
}

int //!!
TclModelBuilder::addCyclicModel(CyclicModel &theCM)
{
//	TaggedObject *mc = &theYS;

  bool result = theCycModels->addComponent(&theCM);
  if (result == true)
    return 0;
  else {
    opserr << "TclModelBuilder::addCyclicModel() - failed to add : " << theCM;
    return -1;
  }
}

YieldSurface_BC *
TclModelBuilder::getYieldSurface_BC(int tag)
{
  TaggedObject *mc = theYieldSurface_BCs->getComponentPtr(tag);
  if (mc == 0)
    return 0;

  // otherweise we do a cast and return
  YieldSurface_BC *result = (YieldSurface_BC *)mc;
  return result;
}

CyclicModel * //!!
TclModelBuilder::getCyclicModel(int tag)
{
  TaggedObject *mc = theCycModels->getComponentPtr(tag);
  if (mc == 0)
    return 0;

  // otherweise we do a cast and return
  CyclicModel *result = (CyclicModel *)mc;
  return result;
}

int
TclModelBuilder::addPlasticMaterial(PlasticHardeningMaterial &theMat)
{
//	TaggedObject *mc = &theYS;

  bool result = thePlasticMaterials->addComponent(&theMat);
  if (result == true)
    return 0;
  else {
    opserr << "TclModelBuilder::addPlasticMaterial() - failed to add Material: " << theMat;
    return -1;
  }
}

PlasticHardeningMaterial *
TclModelBuilder::getPlasticMaterial(int tag)
{
  TaggedObject *mc = thePlasticMaterials->getComponentPtr(tag);
  if (mc == 0)
    return 0;

  // otherweise we do a cast and return
  PlasticHardeningMaterial *result = (PlasticHardeningMaterial *)mc;
  return result;
}

int 
TclModelBuilder::addSectionRepres(SectionRepres &theSectionRepres)
{
  bool result = theSectionRepresents->addComponent(&theSectionRepres);

  if (result == true)
    return 0;
  else {
      opserr << "TclModelBuilder::addSectionRepres() - failed to add SectionRepres\n";
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
    opserr << "TclModelBuilder::addCrdTransf() - failed to add crdTransf: " << theCrdTransf;
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
    opserr << "TclModelBuilder::addCrdTransf() - failed to add crdTransf: " << theCrdTransf;
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

void printCommand(int argc, TCL_Char **argv)
{
  opserr << "Input command: ";
  for (int i=0; i<argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
} 

int
TclModelBuilder_addNode(ClientData clientData, Tcl_Interp *interp, int argc, 
                        TCL_Char **argv)
{

  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed" << endln;
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  // make sure corect number of arguments on command line
  if (argc < 2+ndm) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }    

  Node *theNode = 0;

  // get the nodal id
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeTag\n";
    opserr << "Want: node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }

  // read in the coordinates and create the node
  double xLoc, yLoc, zLoc;
  if (ndm == 1) { 
    // create a node in 1d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << "WARNING invalid XCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    theNode = new Node(nodeId,ndf,xLoc);
  } 

  else if (ndm == 2) { 
    // create a node in 2d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << "WARNING invalid XCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      opserr << "WARNING invalid YCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    theNode = new Node(nodeId,ndf,xLoc,yLoc);
  } 

  else if (ndm == 3) { 
    // create a node in 3d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << "WARNING invalid XCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      opserr << "WARNING invalid YCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &zLoc) != TCL_OK) {
      opserr << "WARNING invalid ZCoordinate\n";
      opserr << "node: " << nodeId << endln;
      return TCL_ERROR;
    }
    theNode = new Node(nodeId,ndf,xLoc,yLoc,zLoc);
  } else {
      opserr << "WARNING invalid ndm\n";
      opserr << "node: " << nodeId << endln;;
      return TCL_ERROR;
  }

  if (theNode == 0) {
    opserr << "WARNING ran out of memory creating node\n";
    opserr << "node: " << nodeId << endln;
    return TCL_ERROR;
  }

  if (theTclDomain->addNode(theNode) == false) {
    opserr << "WARNING failed to add node to the domain\n";
    opserr << "node: " << nodeId << endln;
    delete theNode; // otherwise memory leak
    return TCL_ERROR;
  }

  // check for mass terms
  if (argc > 2+ndm) {
    if (strcmp(argv[2+ndm],"-mass") == 0) {
      if (argc < 3+ndm+ndf) {
        opserr << "WARNING incorrect number of nodal mass terms\n";
        opserr << "node: " << nodeId << endln;
        return TCL_ERROR;      
      }	
      Matrix mass(ndf,ndf);
      double theMass;
      for (int i=0; i<ndf; i++) {
	if (Tcl_GetDouble(interp, argv[i+3+ndm], &theMass) != TCL_OK) {
	  opserr << "WARNING invalid nodal mass term\n";
	  opserr << "node: " << nodeId << ", dof: " << i+1 << endln;
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
			      TCL_Char **argv, 
			      Domain *theDomain, TclModelBuilder *theTclBuilder);
int
TclModelBuilder_addElement(ClientData clientData, Tcl_Interp *interp, 
			   int argc,    TCL_Char **argv)
                          
{
  return TclModelBuilderElementCommand(clientData, interp, 
				       argc, argv, theTclDomain, theTclBuilder);
}


extern int
TclModelBuilderUniaxialMaterialCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclBuilder, Domain *theDomain);

int
TclModelBuilder_addUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, 
				    int argc, TCL_Char **argv)
                          
{
  return TclModelBuilderUniaxialMaterialCommand(clientData, interp, 
						argc, argv, theTclBuilder, theTclDomain);
}

extern int
TclModelBuilderNDMaterialCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				  TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addNDMaterial(ClientData clientData, Tcl_Interp *interp, 
			    int argc,    TCL_Char **argv)
                          
{
  return TclModelBuilderNDMaterialCommand(clientData, interp, 
					  argc, argv, theTclBuilder);
}

extern int
TclModelBuilderSectionCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				  TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addSection(ClientData clientData, Tcl_Interp *interp, 
			    int argc,    TCL_Char **argv)
                          
{
  return TclModelBuilderSectionCommand(clientData, interp, 
				       argc, argv, theTclBuilder);
}



extern int
TclModelBuilderYieldSurface_BCCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addYieldSurface_BC(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv)

{
  return TclModelBuilderYieldSurface_BCCommand(clientData, interp,
						argc, argv, theTclBuilder);
}

extern int
TclModelBuilderYS_EvolutionModelCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addYS_EvolutionModel(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv)

{
  return TclModelBuilderYS_EvolutionModelCommand(clientData, interp,
						argc, argv, theTclBuilder);
}

extern int
TclModelBuilderPlasticMaterialCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addYS_PlasticMaterial(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv)

{
  return TclModelBuilderPlasticMaterialCommand(clientData, interp,
						argc, argv, theTclBuilder);
}

//!!
extern int TclModelBuilderCyclicModelCommand(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclBuilder);
int
TclModelBuilder_addCyclicModel(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv)

{
  return TclModelBuilderCyclicModelCommand(clientData, interp,
						argc, argv, theTclBuilder);
}

//!!
extern int TclModelBuilderDamageModelCommand(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilder_addDamageModel(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv)

{
  return TclModelBuilderDamageModelCommand(clientData, interp,
						argc, argv, theTclBuilder);
}

extern int
TclPatternCommand(ClientData clientData, Tcl_Interp *interp, 
			   int argc, TCL_Char **argv, Domain *theDomain);
			   
int
TclModelBuilder_addPattern(ClientData clientData, Tcl_Interp *interp, 
			   int argc, TCL_Char **argv)
			  
{
  return TclPatternCommand(clientData, interp, argc, argv, theTclDomain);
}




extern int
TclGroundMotionCommand(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc,    
		       TCL_Char **argv,
		       MultiSupportPattern *thePattern);

int
TclModelBuilder_addGroundMotion(ClientData clientData, Tcl_Interp *interp, 
			   int argc, TCL_Char **argv)
			  
{
  return TclGroundMotionCommand(clientData, interp, argc, argv, 
				theTclMultiSupportPattern);
}


int
TclModelBuilder_addNodalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - load \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  NodalLoad *theLoad = 0;
  
  // make sure at least one other argument to contain type of system
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: load nodeId " << ndf << " forces\n";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeId: " << argv[1];
    opserr << " - load nodeId " << ndf << " forces\n";
    return TCL_ERROR;
  }

  // get the load vector
  Vector forces(ndf);
  for (int i=0; i<ndf; i++) {
    double theForce;
    if (Tcl_GetDouble(interp, argv[2+i], &theForce) != TCL_OK) {
      opserr << "WARNING invalid force " << i+1 << " - load " << nodeId;
      opserr << " " << ndf << " forces\n";
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

	opserr << "WARNING invalid patternTag - load " << nodeId << " ";
	opserr << ndf << " forces pattern patterntag\n";
	return TCL_ERROR;
      }
    }
    endMarker++;
  }

  // get the current pattern tag if no tag given in i/p
  if (loadPatternTag == 123456789)
    if (theTclLoadPattern == 0) {
	opserr << "WARNING no current load pattern - load " << nodeId;
	opserr << " " << ndf << " forces\n";
	return TCL_ERROR;
    } else 
	loadPatternTag = theTclLoadPattern->getTag();

  // create the load
  theLoad = new NodalLoad(nodeLoadTag, nodeId, forces, isLoadConst);
  if (theLoad == 0) {
    opserr << "WARNING ran out of memory for load  - load " << nodeId;
    opserr << " " << ndf << " forces\n";
    return TCL_ERROR;
  }

  // add the load to the domain
  if (theTclDomain->addNodalLoad(theLoad, loadPatternTag) == false) {
    opserr << "WARNING TclModelBuilder - could not add load to domain ";
    printCommand(argc, argv);
    delete theLoad;
    return TCL_ERROR;
  }
  nodeLoadTag++;

  // if get here we have sucessfully created the load and added it to the domain
  return TCL_OK;
}




int
TclModelBuilder_addElementalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING current builder has been destroyed - eleLoad\n";    
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  ElementalLoad *theLoad = 0;

  ID theEleTags(0,16);

  // we first create an ID containing the ele tags of all elements
  // for which the load applies.
  int count = 1;
  int doneEle = 0;
  int eleCount = 0;
  while (doneEle == 0 && count < argc) {
    if (strcmp(argv[count],"-ele") == 0) {
      count ++;
      int eleStart = count;
      int eleEnd = 0;
      int eleID;
      while (count < argc && eleEnd == 0) {
	if (Tcl_GetInt(interp, argv[count], &eleID) != TCL_OK)
	  eleEnd = count;
	else
	  count++;
      }
      if (eleStart != eleEnd) {
	for (int i=eleStart; i<eleEnd; i++) {
	  Tcl_GetInt(interp, argv[i], &eleID);
	  theEleTags[eleCount++] = eleID;
	}
      }
    }
    else if (strcmp(argv[count],"-range") == 0) {
      count ++;
      int eleStart, eleEnd;
      if (Tcl_GetInt(interp, argv[count], &eleStart) != TCL_OK) {
	opserr << "WARNING eleLoad -range invalid eleStart " << argv[count] << "\n";
	return TCL_ERROR;
      }
      count++;
      if (Tcl_GetInt(interp, argv[count], &eleEnd) != TCL_OK) {
	opserr << "WARNING eleLoad -range invalid eleEnd " << argv[count] << "\n";	
	return TCL_ERROR;
      }
      count++;
      for (int i=eleStart; i<=eleEnd; i++)
	theEleTags[eleCount++] = i;	
    } else
      doneEle = 1;
  }


  // we then create the load
  if (strcmp(argv[count],"-type") != 0) {
    opserr << "WARNING eleLoad - expecting -type option but got "
	 << argv[count] << endln;
    return TCL_ERROR;
  } 
  count++;
  if (strcmp(argv[count],"-beamUniform") == 0 ||
      strcmp(argv[count],"beamUniform") == 0){
    count++;
    if (ndm == 2) {
      double wt;
      double wa = 0.0;
      if (Tcl_GetDouble(interp, argv[count], &wt) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wt " << argv[count]
	     << " for beamUniform \n";
	return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wa) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wa " << argv[count]
	     << " for beamUniform \n";
	return TCL_ERROR;
      }
      theLoad = new Beam2dUniformLoad(eleLoadTag, wt, wa, theEleTags);    
    }
    else if (ndm == 3) {
      double wy, wz;
      double wx = 0.0;
      if (Tcl_GetDouble(interp, argv[count], &wy) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wy " << argv[count]
	     << " for beamUniform \n";
	return TCL_ERROR;
      }
      count++;
      if (Tcl_GetDouble(interp, argv[count], &wz) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wz " << argv[count]
	     << " for beamUniform \n";
	return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wx) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wx " << argv[count]
	     << " for beamUniform \n";
	return TCL_ERROR;
      }
      theLoad = new Beam3dUniformLoad(eleLoadTag, wy, wz, wx, theEleTags);    
    }
    else { 
      opserr << "WARNING eleLoad beamUniform currently only valid only for ndm=2 or 3\n";     
      return TCL_ERROR;
    }
  } else if (strcmp(argv[count],"-beamPoint") == 0 ||
	     strcmp(argv[count],"beamPoint") == 0 ) {
    count++;
    if (ndm == 2) {
      double P, x;
      double N = 0.0;
      if (Tcl_GetDouble(interp, argv[count], &P) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid P " << argv[count] << " for beamPoint\n";		
	return TCL_ERROR;
      } 
      if (Tcl_GetDouble(interp, argv[count+1], &x) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid xDivL " << argv[count+1] << " for beamPoint\n";	
	return TCL_ERROR;
      } 
      if (count+2 < argc && Tcl_GetDouble(interp, argv[count+2], &N) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid N " << argv[count+2] << " for beamPoint\n";		
	return TCL_ERROR;
      } 

      if (x < 0.0 || x > 1.0) {
	opserr << "WARNING eleLoad - invalid xDivL of " << x;
	opserr << " for beamPoint (valid range [0.0, 1.0]\n";
	return TCL_ERROR;
      }

      theLoad = new Beam2dPointLoad(eleLoadTag, P, x, theEleTags, N);    
    }
    else if (ndm == 3) {
      double Py, Pz, x;
      double N = 0.0;
      if (Tcl_GetDouble(interp, argv[count], &Py) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid Py " << argv[count] << " for beamPoint\n";		
	return TCL_ERROR;
      } 
      if (Tcl_GetDouble(interp, argv[count+1], &Pz) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid Pz " << argv[count] << " for beamPoint\n";		
	return TCL_ERROR;
      } 
      if (Tcl_GetDouble(interp, argv[count+2], &x) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid xDivL " << argv[count+2] << " for beamPoint\n";	
	return TCL_ERROR;
      } 
      if (count+3 < argc && Tcl_GetDouble(interp, argv[count+3], &N) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid N " << argv[count+3] << " for beamPoint\n";		
	return TCL_ERROR;
      } 

      if (x < 0.0 || x > 1.0) {
	opserr << "WARNING eleLoad - invalid xDivL of " << x;
	opserr << " for beamPoint (valid range [0.0, 1.0]\n";
	return TCL_ERROR;
      }

      theLoad = new Beam3dPointLoad(eleLoadTag, Py, Pz, x, theEleTags, N);    
    }
    else {
      opserr << "WARNING eleLoad beamPoint type currently only valid only for ndm=2 or 3\n";
      return TCL_ERROR;
    }  
  }
  // Added Joey Yang UC Davis
  else if (strcmp(argv[count],"-BrickW") == 0) {
      theLoad = new BrickSelfWeight(eleLoadTag, theEleTags);
  }

  // Added by Scott R. Hamilton   - Stanford
  else if (strcmp(argv[count],"-beamTemp") == 0) {
    count++;
    if (ndm == 2) {
      double temp1, temp2, temp3, temp4;
      
      // Four temps given, Temp change at top node 1, bottom node 1, top node 2, bottom node 2.
      if (argc-count == 4){
	if (Tcl_GetDouble(interp, argv[count], &temp1) != TCL_OK) {
	  opserr << "WARNING eleLoad - invalid Ttop1 " << argv[count] << " for -beamTemp\n";		
	  return TCL_ERROR;
	} 
      
	if (Tcl_GetDouble(interp, argv[count+1],&temp2 ) != TCL_OK) {
	  opserr << "WARNING eleLoad - invalid Tbot1 " << argv[count+1] << " for -beamTemp\n";	
	  return TCL_ERROR;
	} 
	if (Tcl_GetDouble(interp, argv[count+2], &temp3) != TCL_OK) {
	  opserr << "WARNING eleLoad - invalid Ttop2 " << argv[count+1] << " for -beamTemp\n";	
	  return TCL_ERROR;
	} 
	if (Tcl_GetDouble(interp, argv[count+3], &temp4) != TCL_OK) {
	  opserr << "WARNING eleLoad - invalid Tbot2 " << argv[count+1] << " for -beamTemp\n";	
	  return TCL_ERROR;
	} 
	
	theLoad=0;
	theLoad = new Beam2dTempLoad(eleLoadTag, temp1, temp2, temp3, temp4, theEleTags);
      }
      // Two temps given, temp change at top, temp at bottom of element
      else if (argc-count == 2) {
	if (Tcl_GetDouble(interp, argv[count], &temp1) != TCL_OK) {
	  opserr << "WARNING eleLoad - invalid Ttop " << argv[count] << " for -beamTemp\n";		
	  return TCL_ERROR;
	} 
	
	if (Tcl_GetDouble(interp, argv[count+1],&temp2 ) != TCL_OK) {
	  opserr << "WARNING eleLoad - invalid Tbot " << argv[count+1] << " for -beamTemp\n";	
	  return TCL_ERROR;
	}
	theLoad=0;
	theLoad = new Beam2dTempLoad(eleLoadTag, temp1, temp2, theEleTags);
      }
      // One twmp change give, uniform temp change in element
      else if (argc-count == 1) {
	if (Tcl_GetDouble(interp, argv[count],&temp1 ) != TCL_OK) {
	  opserr << "WARNING eleLoad - invalid Tbot " << argv[count+1] << " for -beamTemp\n";	
	  return TCL_ERROR;
	}
	theLoad=0;
	theLoad = new Beam2dTempLoad(eleLoadTag, temp1, theEleTags);
      }
      // No temps, no change in temp in element--not a case likely to be used
      else if (argc-count == 0){
	theLoad=0;
	theLoad = new Beam2dTempLoad(eleLoadTag, theEleTags);
      }
      else {
	opserr << "WARNING eleLoad -beamTempLoad invalid number of temperature aguments,/n looking for 0, 1, 2 or 4 arguments.\n";
      }
    } else {
      opserr << "WARNING eleLoad -beamTempLoad type currently only valid only for ndm=2\n";
      return TCL_ERROR;
    }  
  }

  if (theLoad == 0) {
    opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count] ;
    return TCL_ERROR;
  }

  // get the current pattern tag if no tag given in i/p
  int loadPatternTag = theTclLoadPattern->getTag();

  // add the load to the domain
  if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
    opserr << "WARNING eleLoad - could not add following load to domain:\n ";
    opserr << theLoad;
    delete theLoad;
    return TCL_ERROR;
  }
  eleLoadTag++;

  // if get here we have sucessfully created the load and added it to the domain
  return TCL_OK;
}



int
TclModelBuilder_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc, 
                        TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - load \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // make sure at least one other argument to contain type of system
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: mass nodeId " << ndf << " mass values\n";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeId: " << argv[1];
    opserr << " - mass nodeId " << ndf << " forces\n";
    return TCL_ERROR;
  }

  Node *theNode = 0;
  theNode = theTclDomain->getNode(nodeId);

  if (theNode == 0)
  {
    opserr << "WARNING failed to get node pointer from the domain\n";
    opserr << "node: " << nodeId << endln;
    return TCL_ERROR;
  }

  // check for mass terms
  Matrix mass(ndf,ndf);
  double theMass;
  for (int i=0; i<ndf; i++) 
  {
     if (Tcl_GetDouble(interp, argv[i+2], &theMass) != TCL_OK) 
     {
	  opserr << "WARNING invalid nodal mass term\n";
	  opserr << "node: " << nodeId << ", dof: " << i+1 << endln;
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
				 TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: fix nodeId " << ndf << " [0,1] conditions";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
      opserr << "WARNING invalid nodeId - fix nodeId " << ndf << " [0,1] conditions\n";
      return TCL_ERROR;
  }

  // get the fixity condition and add the constraint if fixed
  for (int i=0; i<ndf; i++) {
    int theFixity;
    if (Tcl_GetInt(interp, argv[2+i], &theFixity) != TCL_OK) {
      opserr << "WARNING invalid fixity " << i+1 << " - load " << nodeId;
      opserr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    } else {
      if (theFixity != 0) {
	// create a homogeneous constraint
	SP_Constraint *theSP = new SP_Constraint(currentSpTag, nodeId, i, 0.0);
	if (theSP == 0) {
	  opserr << "WARNING ran out of memory for SP_Constraint ";
	  opserr << "fix " << nodeId << " " << ndf << " [0,1] conditions\n";
	  return TCL_ERROR;
	}
	if (theTclDomain->addSP_Constraint(theSP) == false) {
	  opserr << "WARNING could not add SP_Constraint to domain - fix";
	  opserr << nodeId << " " << ndf << " [0,1] conditions\n";
	  delete theSP;
	  return TCL_ERROR;
	}
	currentSpTag++;      }
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}


int
TclModelBuilder_addHomogeneousBC_X(ClientData clientData, Tcl_Interp *interp, 
				   int argc, TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: fixX xLoc " << ndf << " [0,1] conditions";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the xCrd of nodes to be constrained
  double xLoc;
  if (Tcl_GetDouble(interp, argv[1], &xLoc) != TCL_OK) {
      opserr << "WARNING invalid xCrd - fixX xLoc " << ndf << " [0,1] conditions\n";
      return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i=0; i<ndf; i++) {
    if (Tcl_GetInt(interp, argv[2+i], &fixity(i)) != TCL_OK) {
      opserr << "WARNING invalid fixity " << i+1 << " - fixX " << xLoc;
      opserr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    } 
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2+ndf],"-tol") == 0)
    if (Tcl_GetDouble(interp, argv[3+ndf], &tol) != TCL_OK) {
      opserr << "WARNING invalid tol specified - fixX " << xLoc << endln;
      return TCL_ERROR;
    }       
  }

  NodeIter &theNodes = theTclDomain->getNodes();
  Node *theNode;

  // loop over all the nodes
  while ((theNode = theNodes()) != 0) {
    const Vector &theCrd = theNode->getCrds();
    double nodeX = theCrd(0);

    // add a single point constraint if Xcrd of node is within tol of xLoc
    if (fabs(nodeX - xLoc) < tol) {
      int nodeId = theNode->getTag();
      int theFixity = 0;

      // loop over all the ndf values valid for the node
      int numDOF = theNode->getNumberDOF();
      if (numDOF  < ndf) numDOF = ndf;

      for (int i=0; i<numDOF; i++) {
	theFixity = fixity(i);
	if (theFixity != 0) {
	  // create a homogeneous constraint
	  SP_Constraint *theSP = new SP_Constraint(currentSpTag, nodeId, i, 0.0);
	  if (theSP == 0) {
	    opserr << "WARNING ran out of memory for SP_Constraint at node " << nodeId;
	    opserr << " - fixX " << xLoc << " " << ndf << " [0,1] conditions\n";
	    return TCL_ERROR;
	  }
	  if (theTclDomain->addSP_Constraint(theSP) == false) {
	    opserr << "WARNING could not add SP_Constraint to domain for node " << nodeId;
	    opserr << " - fixX " << xLoc << " " << ndf << " [0,1] conditions\n";
	    delete theSP;
	    return TCL_ERROR;
	  }
	  currentSpTag++;      
	}
      }
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}




int
TclModelBuilder_addHomogeneousBC_Y(ClientData clientData, Tcl_Interp *interp, 
				   int argc, TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: fixY yLoc " << ndf << " [0,1] conditions";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the yCrd of nodes to be constrained
  double yLoc;
  if (Tcl_GetDouble(interp, argv[1], &yLoc) != TCL_OK) {
      opserr << "WARNING invalid yCrd - fixY yLoc " << ndf << " [0,1] conditions\n";
      return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i=0; i<ndf; i++) {
    if (Tcl_GetInt(interp, argv[2+i], &fixity(i)) != TCL_OK) {
      opserr << "WARNING invalid fixity " << i+1 << " - fixY " << yLoc;
      opserr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    } 
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2+ndf],"-tol") == 0)
    if (Tcl_GetDouble(interp, argv[3+ndf], &tol) != TCL_OK) {
      opserr << "WARNING invalid tol specified - fixY " << yLoc << endln;
      return TCL_ERROR;
    }       
  }

  NodeIter &theNodes = theTclDomain->getNodes();
  Node *theNode;

  // loop over all the nodes
  while ((theNode = theNodes()) != 0) {
    const Vector &theCrd = theNode->getCrds();
    if (theCrd.Size() > 1) {
      double nodeY = theCrd(1);

      // add a single point constraint if Xcrd of node is within tol of yLoc
      if (fabs(nodeY - yLoc) < tol) {

	int nodeId = theNode->getTag();
	int theFixity = 0;

	// loop over all the ndf values valid for the node
	int numDOF = theNode->getNumberDOF();
	if (numDOF  < ndf) numDOF = ndf;

	for (int i=0; i<numDOF; i++) {
	  theFixity = fixity(i);
	  if (theFixity != 0) {
	    // create a homogeneous constraint
	    SP_Constraint *theSP = new SP_Constraint(currentSpTag, nodeId, i, 0.0);
	    if (theSP == 0) {
	      opserr << "WARNING ran out of memory for SP_Constraint at node " << nodeId;
	      opserr << " - fixY " << yLoc << " " << ndf << " [0,1] conditions\n";
	      return TCL_ERROR;
	    }
	    if (theTclDomain->addSP_Constraint(theSP) == false) {
	      opserr << "WARNING could not add SP_Constraint to domain for node " << nodeId;
	      opserr << " - fixY " << yLoc << " " << ndf << " [0,1] conditions\n";
	      delete theSP;
	      return TCL_ERROR;
	    }
	    currentSpTag++;      
	  }
	}
      }
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclModelBuilder_addHomogeneousBC_Z(ClientData clientData, Tcl_Interp *interp, 
				   int argc, TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << "WARNING bad command - want: fixZ zLoc " << ndf << " [0,1] conditions";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the yCrd of nodes to be constrained
  double zLoc;
  if (Tcl_GetDouble(interp, argv[1], &zLoc) != TCL_OK) {
      opserr << "WARNING invalid zCrd - fixZ zLoc " << ndf << " [0,1] conditions\n";
      return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i=0; i<ndf; i++) {
    if (Tcl_GetInt(interp, argv[2+i], &fixity(i)) != TCL_OK) {
      opserr << "WARNING invalid fixity " << i+1 << " - fixZ " << zLoc;
      opserr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    } 
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2+ndf],"-tol") == 0)
    if (Tcl_GetDouble(interp, argv[3+ndf], &tol) != TCL_OK) {
      opserr << "WARNING invalid tol specified - fixZ " << zLoc << endln;
      return TCL_ERROR;
    }       
  }

  NodeIter &theNodes = theTclDomain->getNodes();
  Node *theNode;

  // loop over all the nodes
  while ((theNode = theNodes()) != 0) {
    const Vector &theCrd = theNode->getCrds();
    if (theCrd.Size() > 2) {
      double nodeZ = theCrd(2);

      // add a single point constraint if Xcrd of node is within tol of zLoc
      if (fabs(nodeZ - zLoc) < tol) {

	int nodeId = theNode->getTag();
	int theFixity = 0;

	// loop over all the ndf values valid for the node
	int numDOF = theNode->getNumberDOF();
	if (numDOF  < ndf) numDOF = ndf;

	for (int i=0; i<numDOF; i++) {
	  theFixity = fixity(i);
	  if (theFixity != 0) {
	    // create a homogeneous constraint
	    SP_Constraint *theSP = new SP_Constraint(currentSpTag, nodeId, i, 0.0);
	    if (theSP == 0) {
	      opserr << "WARNING ran out of memory for SP_Constraint at node " << nodeId;
	      opserr << " - fixZ " << zLoc << " " << ndf << " [0,1] conditions\n";
	      return TCL_ERROR;
	    }
	    if (theTclDomain->addSP_Constraint(theSP) == false) {
	      opserr << "WARNING could not add SP_Constraint to domain for node " << nodeId;
	      opserr << " - fixZ " << zLoc << " " << ndf << " [0,1] conditions\n";
	      delete theSP;
	      return TCL_ERROR;
	    }
	    currentSpTag++;      
	  }
	}
      }
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}





int
TclModelBuilder_addSP(ClientData clientData, Tcl_Interp *interp, int argc,   
		      TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - sp \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // check number of arguments
  if (argc < 4) {
    opserr << "WARNING bad command - want: sp nodeId dofID value";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the nodeID, dofId and value of the constraint
  int nodeId, dofId;
  double value;

  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeId: " << argv[1] << " -  sp nodeId dofID value\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dofId) != TCL_OK) {
    opserr << "WARNING invalid dofId: " << argv[2] << " -  sp ";
    opserr << nodeId << " dofID value\n";
      return TCL_ERROR;
  }
  dofId--; // DECREMENT THE DOF VALUE BY 1 TO GO TO OUR C++ INDEXING

  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr << "WARNING invalid value: " << argv[3] << " -  sp ";
    opserr << nodeId << " dofID value\n";
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

	opserr << "WARNING invalid patternTag - load " << nodeId << " ";
	opserr << ndf << " forces pattern patterntag\n";
	return TCL_ERROR;
      }
    }  
    endMarker++;
  }

  // if load pattern tag has not changed - get the pattern tag from current one
  if (loadPatternTag == 123456789) {
    if (theTclLoadPattern == 0) {
      opserr << "WARNING no current pattern - sp " << nodeId << " dofID value\n";
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
    opserr << "WARNING ran out of memory for SP_Constraint ";
    opserr << " - sp " << nodeId << " dofID value\n";
    return TCL_ERROR;
  }
  if (theTclDomain->addSP_Constraint(theSP, loadPatternTag) == false) {
    opserr << "WARNING could not add SP_Constraint to domain ";
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
				   TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - sp \n";    
    return TCL_ERROR;
  }

  int ndf = theTclBuilder->getNDF();

  // check number of arguments
  if (argc < 4) {
    opserr << "WARNING bad command - want: imposedMotion nodeId dofID gMotionID\n";
    printCommand(argc, argv);
    return TCL_ERROR;
  }    

  // get the nodeID, dofId and value of the constraint
  int nodeId, dofId, gMotionID;

  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << "WARNING invalid nodeId: " << argv[1];
    opserr << " - imposedMotion nodeId dofID gMotionID\n";    
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dofId) != TCL_OK) {
    opserr << "WARNING invalid dofId: " << argv[2] << " -  imposedMotion ";
    opserr << nodeId << " dofID gMotionID\n";    
      return TCL_ERROR;
  }
  dofId--; // DECREMENT THE DOF VALUE BY 1 TO GO TO OUR C++ INDEXING

  if (Tcl_GetInt(interp, argv[3], &gMotionID) != TCL_OK) {
    opserr << "WARNING invalid gMotionID: " << argv[3] << " -  imposedMotion ";
    opserr << nodeId << " dofID gMotionID\n";
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
    opserr << "WARNING no GroundMotion with tag: " << argv[3];
    opserr << " in current MultipleSupportPattern\n";
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
    opserr << "WARNING ran out of memory for ImposedMotionSP ";
    opserr << " -  imposedMotion ";
    opserr << nodeId << " " << dofId++ << " " << gMotionID << endln;
    return TCL_ERROR;
  }
  if (thePattern->addSP_Constraint(theSP) == false) {
    opserr << "WARNING could not add SP_Constraint to pattern ";
    printCommand(argc, argv);
    delete theSP;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}






int
TclModelBuilder_addEqualDOF_MP (ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char **argv)
{
        // Ensure the destructor has not been called
        if (theTclBuilder == 0) {
	  opserr << "WARNING builder has been destroyed - equalDOF \n";
	  return TCL_ERROR;
        }

        // Check number of arguments
        if (argc < 4) {
	  opserr << "WARNING bad command - want: equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
	  printCommand (argc, argv);
	  return TCL_ERROR;
        }

        // Read in the node IDs and the DOF
        int RnodeID, CnodeID, dofID;

        if (Tcl_GetInt (interp, argv[1], &RnodeID) != TCL_OK) {
	  opserr << "WARNING invalid RnodeID: " << argv[1]
	       << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
	  return TCL_ERROR;
        }
        if (Tcl_GetInt (interp, argv[2], &CnodeID) != TCL_OK) {
	  opserr << "WARNING invalid CnodeID: " << argv[2]
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
        for (i = 3, j = 0; i < argc; i++, j++) {
	  if (Tcl_GetInt (interp, argv[i], &dofID) != TCL_OK) {
	    opserr << "WARNING invalid dofID: " << argv[3]
		 << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
	    return TCL_ERROR;
	  }

	  dofID -= 1; // Decrement for C++ indexing
	  if (dofID < 0) {
	    opserr << "WARNING invalid dofID: " << argv[i]
		   << " must be >= 1";
	    return TCL_ERROR;
	  }
	  rcDOF (j) = dofID;    
	  Ccr (j,j) = 1.0;
        }

        // Use this for the MP tag
        int numMPs = theTclDomain->getNumMPs();

        // Create the multi-point constraint
        MP_Constraint *theMP = new MP_Constraint (numMPs, RnodeID, CnodeID, Ccr, rcDOF, rcDOF);
        if (theMP == 0) {
	  opserr << "WARNING ran out of memory for equalDOF MP_Constraint ";
	  printCommand (argc, argv);
	  return TCL_ERROR;
        }

        // Add the multi-point constraint to the domain
        if (theTclDomain->addMP_Constraint (theMP) == false) {
	  opserr << "WARNING could not add equalDOF MP_Constraint to domain ";
	  printCommand(argc, argv);
	  delete theMP;
	  return TCL_ERROR;
        }

        return TCL_OK;
}



int
TclModelBuilder_addMP(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv)
{
  opserr << "WARNING - TclModelBuilder_addMP() not yet implemented\n";
  return TCL_OK;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Added by Scott J. Brandenberg, UC Davis, sjbrandenberg@ucdavis.edu
int
TclModelBuilder_doPySimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char **argv)
{
    	if(argc < 6 || argc > 7){
		opserr << "WARNING PySimple1Gen file1? file2? file3? file4? file5? <file6?>";
		opserr << "Must have either 5 or 6 arguments." << endln;
	}
	
	PySimple1Gen *thePySimple1Gen;
	thePySimple1Gen = new PySimple1Gen;

	if(argc==6)
		thePySimple1Gen->WritePySimple1(argv[1], argv[2], argv[3], argv[4], argv[5]);
	if(argc==7)
		thePySimple1Gen->WritePySimple1(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);

	delete thePySimple1Gen;

	return TCL_OK;
}

int
TclModelBuilder_doTzSimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char **argv)
{
	if(argc < 6 || argc > 7){
		opserr << "WARNING TzSimple1Gen file1? file2? file3? file4? file5? <file6?>";
		opserr << "Must have either 5 or 6 arguments." << endln;
	}

	TzSimple1Gen *theTzSimple1Gen;
	theTzSimple1Gen = new TzSimple1Gen;

	if(argc==6)
		theTzSimple1Gen->WriteTzSimple1(argv[1], argv[2], argv[3], argv[4], argv[5]);
	if(argc==7)
		theTzSimple1Gen->WriteTzSimple1(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);

	delete theTzSimple1Gen;

	return TCL_OK;
}
// End Added by Scott J. Brandenberg
///////////////////////////////////////////////////////////////////////////////////////////////////	

int
TclModelBuilder_doBlock2D(ClientData clientData, Tcl_Interp *interp, int argc,   
			  TCL_Char **argv)
{

  int ndm = theTclBuilder->getNDM();
  if (ndm < 2) {
    opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : model dimension (ndm) must be at leat 2 " << endln;
    return TCL_ERROR;
  }

  if (argc < 8) {
    opserr << "WARNING incorrect numer of args :block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    return TCL_ERROR;
  }
  int numX, numY, startNodeNum, startEleNum;
  if (Tcl_GetInt (interp, argv[1], &numX) != TCL_OK) {
    opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numX: " << argv[1] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[2], &numY) != TCL_OK) {
    opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numY: " << argv[2] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[3], &startNodeNum) != TCL_OK) {
    opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startNode: " << argv[3] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[4], &startEleNum) != TCL_OK) {
    opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startEle: " << argv[4] << endln;
    return TCL_ERROR;
  }


  static Matrix Coordinates(9,3);
  static ID     haveNode(9);
  Coordinates.Zero();
  for (int k=0; k<9; k++) haveNode(k) = -1;

  int numNodes = 4;
  if (argc == 10) {
    if (strcmp(argv[7],"-numEleNodes") == 0) 
      if (Tcl_GetInt (interp, argv[8], &numNodes) != TCL_OK) {
	opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
	opserr << " -numEleNodes numNodes?: invalid numNodes: " << argv[8] << endln;
	return TCL_ERROR;
      }
    if (numNodes != 4 && numNodes != 9) {
      opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs? ";
      opserr << "-numEleNodes numNodes?: invalid numNodes: " << argv[8] << " 4 or 9 only\n";
      return TCL_ERROR;
    }

    if (numNodes == 9) {
      if (((numX % 2) != 0) || ((numY % 2) != 0)) {
	opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs? ";
	opserr << "-numEleNodes 9: numX and numY MUST BOTH BE EVEN\n";
	return TCL_ERROR;
      }
    }
  }


  //  char *nodalInfo;
  int nodalInfo = 9;
  if (numNodes == 4)
    //nodalInfo = argv[7];
    nodalInfo=7;
  //  else
  //nodalInfo = argv[9];

  TCL_Char **argvNodes;
  int  argcNodes;
  
  Tcl_SplitList(interp, argv[nodalInfo], &argcNodes, &argvNodes);

  int ndf = theTclBuilder->getNDF();
  
  int count = 0;
  while (count < argcNodes) {
    if ((count + ndm + 1) >  argcNodes) {
      opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid number of node args: " << argv[7] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR; 
    }
    int nodeTag;
    double value;
    if (Tcl_GetInt (interp, argvNodes[count], &nodeTag) != TCL_OK) {
      opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid node tag: " << argvNodes[count] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR; 
    }
    if (nodeTag < 1 || nodeTag > 9) {
      opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid node tag out of bounds [1,9]: " << argvNodes[count] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    for (int i=0; i<ndm; i++) {
      if (Tcl_GetDouble(interp, argvNodes[count+1+i], &value) != TCL_OK) {
	opserr << "WARNING block2D numX? numY? startNode? startEle? eleType? eleArgs?";
	opserr << " : invalid node coordinate for node: " << argvNodes[count] << endln;
	Tcl_Free((char *)argvNodes);
	return TCL_ERROR;
      }
      Coordinates(nodeTag-1,i) = value;
      haveNode(nodeTag-1) = nodeTag;
    }      
    count += 1 + ndm;
  }

  Tcl_Free((char *)argvNodes);

  Block2D  theBlock(numX, numY, haveNode, Coordinates, numNodes);

  // create the nodes: (numX+1)*(numY+1) nodes to be created
  int nodeID = startNodeNum;
  int jj;
  for (jj=0; jj<=numY; jj++) {
    for (int ii=0; ii<=numX; ii++) {
      const Vector &nodeCoords = theBlock.getNodalCoords(ii,jj);
      double xLoc = nodeCoords(0);
      double yLoc = nodeCoords(1);
      double zLoc = nodeCoords(2);
      Node *theNode = 0;
      if (ndm == 2) {
	theNode = new Node(nodeID,ndf,xLoc, yLoc);
      } else if (ndm == 3) {
	theNode = new Node(nodeID,ndf,xLoc, yLoc, zLoc);
      } 

      if (theNode == 0) {
	opserr << "WARNING ran out of memory creating node\n";
	opserr << "node: " << nodeID << endln;
	return TCL_ERROR;
      }

      if (theTclDomain->addNode(theNode) == false) {
	opserr << "WARNING failed to add node to the domain\n";
	opserr << "node: " << nodeID << endln;
	delete theNode; // otherwise memory leak
	return TCL_ERROR;
      }

      nodeID++;
    }
  }
    
  // create the elements: numX*numY elements to be created if 4 node elements
  //                      numX/2 * numY /2 nodes to be v=created if 9 node elements
  TCL_Char *eleType = argv[5];
  TCL_Char *additionalEleArgs = argv[6];
  //  const ID &nodeTags = theBlock.getElementNodes(0,0);  
  //  int numNodes = nodeTags.Size();

  // assumes 15 is largest string for individual nodeTags
  count = 10 + strlen(eleType) + strlen(additionalEleArgs) + 15 * (numNodes+1);
  char *eleCommand = new char[count];
  int initialCount = 8 + strlen(eleType);

  int  eleID = startEleNum; 
  if (numNodes == 9) {
    numX /= 2;
    numY /= 2;
  }
    

  for (jj=0; jj<numY; jj++) {
    for (int ii=0; ii<numX; ii++) {
      count = initialCount;

      const ID &nodeTags = theBlock.getElementNodes(ii,jj);
      
      // create the string to be evaluated
      strcpy(eleCommand, "element ");
      strcpy(&eleCommand[8], eleType);
      count += sprintf(&eleCommand[count], " %d ", eleID);
      for (int i=0; i<numNodes; i++) {
	int nodeTag = nodeTags(i)+startNodeNum;
	count += sprintf(&eleCommand[count], " %d ", nodeTag);
      }
      strcat(eleCommand, additionalEleArgs);      

      // now to create the element we get the string eveluated
      if (Tcl_Eval(interp, eleCommand) != TCL_OK) {
          delete [] eleCommand;
	return TCL_ERROR;
      }
      eleID++;
    }
  }

  delete [] eleCommand;
  return TCL_OK;
}


int
TclModelBuilder_doBlock3D(ClientData clientData, Tcl_Interp *interp, int argc,   
			  TCL_Char **argv)
{

  int ndm = theTclBuilder->getNDM();
  if (ndm < 3) {
    opserr << "WARNING block3D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : model dimension (ndm) must be at leat 2 " << endln;
    return TCL_ERROR;
  }

  int numX, numY, numZ, startNodeNum, startEleNum;
  if (Tcl_GetInt (interp, argv[1], &numX) != TCL_OK) {
    opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numX: " << argv[1] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[2], &numY) != TCL_OK) {
    opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numY: " << argv[2] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[3], &numZ) != TCL_OK) {
    opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numZ: " << argv[3] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[4], &startNodeNum) != TCL_OK) {
    opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startNode: " << argv[4] << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt (interp, argv[5], &startEleNum) != TCL_OK) {
    opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startEle: " << argv[5] << endln;
    return TCL_ERROR;
  }

  static Matrix Coordinates(27,3);
  static ID     haveNode(27);
  Coordinates.Zero();
  for (int k=0; k<27; k++) haveNode(k) = -1;

  TCL_Char *nodalInfo = argv[8];
  TCL_Char **argvNodes;
  int  argcNodes;
  
  Tcl_SplitList(interp, nodalInfo, &argcNodes, &argvNodes);

  int ndf = theTclBuilder->getNDF();
  
  int count = 0;
  while (count < argcNodes) {
    if ((count + ndm + 1) > argcNodes) {
      opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid number of node args: " << argv[8] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR; 
    }
    int nodeTag;
    double value;
    if (Tcl_GetInt (interp, argvNodes[count], &nodeTag) != TCL_OK) {
      opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid node id in node args: " << argvNodes[count] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR; 
    }
    if (nodeTag < 1 || nodeTag > 27) {
      opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
      opserr << " : node tag out of bounds [1, 27]: " << argvNodes[count] << endln;
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    for (int i=0; i<ndm; i++) {
      if (Tcl_GetDouble(interp, argvNodes[count+1+i], &value) != TCL_OK) {
	opserr << "WARNING block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
	opserr << " : invalid coordinate in node args: " << argvNodes[count] << endln;
	Tcl_Free((char *)argvNodes);
	return TCL_ERROR;
      }
      Coordinates(nodeTag-1,i) = value;
      haveNode(nodeTag-1) = nodeTag;
    }      
    count += 1 + ndm;
  }

  Tcl_Free((char *)argvNodes);

  Block3D  theBlock(numX, numY, numZ, haveNode, Coordinates);

  // create the nodes: (numX+1)*(numY+1) nodes to be created
  int nodeID = startNodeNum;
  int kk;
  for (kk=0; kk<=numZ; kk++) {
    for (int jj=0; jj<=numY; jj++) {
      for (int ii=0; ii<=numX; ii++) {
	const Vector &nodeCoords = theBlock.getNodalCoords(ii,jj,kk);
	double xLoc = nodeCoords(0);
	double yLoc = nodeCoords(1);
	double zLoc = nodeCoords(2);
	Node *theNode = 0;
	theNode = new Node(nodeID,ndf,xLoc, yLoc, zLoc);
	
	if (theNode == 0) {
	  opserr << "WARNING ran out of memory creating node\n";
	  opserr << "node: " << nodeID << endln;
	  return TCL_ERROR;
	}

	if (theTclDomain->addNode(theNode) == false) {
	  opserr << "WARNING failed to add node to the domain\n";
	  opserr << "node: " << nodeID << endln;
	  delete theNode; // otherwise memory leak
	  return TCL_ERROR;
	}
	
	nodeID++;
      }
    }
  }
    
  // create the elements: numX*numY elements to be created
  TCL_Char *eleType = argv[6];
  TCL_Char *additionalEleArgs = argv[7];
  const ID &nodeTags = theBlock.getElementNodes(0,0,0);  
  int numNodes = nodeTags.Size();

  // assumes 15 is largest string for individual nodeTags
  count = 10 + strlen(eleType) + strlen(additionalEleArgs) + 15 * (numNodes+1);
  char *eleCommand = new char[count];
  int initialCount = 8 + strlen(eleType);

  int  eleID = startEleNum;  
  for (kk=0; kk<numZ; kk++) {
    for (int jj=0; jj<numY; jj++) {
      for (int ii=0; ii<numX; ii++) {
	count = initialCount;

	const ID &nodeTags = theBlock.getElementNodes(ii,jj,kk);
      
	// create the string to be evaluated
	strcpy(eleCommand, "element ");
	strcpy(&eleCommand[8], eleType);
	count += sprintf(&eleCommand[count], " %d ", eleID);
	for (int i=0; i<numNodes; i++) {
	  int nodeTag = nodeTags(i)+startNodeNum;
	  count += sprintf(&eleCommand[count], " %d ", nodeTag);
	}
	strcat(eleCommand, additionalEleArgs);      
	
	// now to create the element we get the string eveluated
	if (Tcl_Eval(interp, eleCommand) != TCL_OK) {
        delete [] eleCommand;
	  return TCL_ERROR;
	}
	eleID++;
      }
    }
  }

  delete [] eleCommand;
  return TCL_OK;
}




int
TclModelBuilder_addRemoPatch(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv)
{
  return TclModelBuilder_addPatch(clientData, interp, argc,argv,
				    theTclBuilder);
}

int
TclModelBuilder_addRemoFiber(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv)
{
  return TclModelBuilder_addFiber(clientData, interp, argc,argv,
				  theTclBuilder);
}

int
TclModelBuilder_addRemoLayer(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv)
{
  return TclModelBuilder_addReinfLayer(clientData, interp, argc,argv,
				       theTclBuilder);
}


					 
int
TclModelBuilder_addRemoGeomTransf(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv)
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
					  TCL_Char **argv, 
					  TclModelBuilder *theTclBuilder);
int
TclModelBuilder_UpdateMaterialStage(ClientData clientData, 
				    Tcl_Interp *interp,  
				    int argc, 
				    TCL_Char **argv)
{
  return TclModelBuilderUpdateMaterialStageCommand(clientData, interp, 
				       argc, argv, theTclBuilder);
}

/// added by ZHY
extern int 
TclModelBuilderUpdateParameterCommand(ClientData clientData, 
					  Tcl_Interp *interp, 
					  int argc, 
					  TCL_Char **argv, 
					  TclModelBuilder *theTclBuilder);
int
TclModelBuilder_UpdateParameter(ClientData clientData, 
				    Tcl_Interp *interp,  
				    int argc, 
				    TCL_Char **argv)
{
  return TclModelBuilderUpdateParameterCommand(clientData, interp, 
				       argc, argv, theTclBuilder);
}

int //!!
TclModelBuilder::addDamageModel(DamageModel &theDM)
{
  bool result = theDamageModels->addComponent(&theDM);
  if (result == true)
    return 0;
  else {
    opserr << "TclModelBuilder::addDamageModel() - failed to add : " << theDM;
    return -1;
  }
}

DamageModel * //!!
TclModelBuilder::getDamageModel(int tag)
{
  TaggedObject *mc = theDamageModels->getComponentPtr(tag);
  if (mc == 0)
    return 0;

  // otherweise we do a cast and return
  DamageModel *result = (DamageModel *)mc;
  return result;
}
