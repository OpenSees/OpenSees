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
                                                                        
// $Revision: 1.52 $
// $Date: 2010-05-12 20:17:50 $
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

#include <RigidRod.h>
#include <RigidBeam.h>
#include <RigidDiaphragm.h>

#include <CrdTransf.h>

#include <NodalLoad.h>
#include <Beam2dPointLoad.h>
#include <Beam2dUniformLoad.h>
#include <Beam2dPartialUniformLoad.h>
#include <Beam2dTempLoad.h>
#include <Beam2dThermalAction.h> //L.Jiang [SIF]
#include <Beam3dThermalAction.h>  //L.Jiang [SIF]
#include <ShellThermalAction.h>   //L.Jiang [SIF]
#include <ThermalActionWrapper.h>  //L.Jiang [SIF]
#include <NodalThermalAction.h>   //L.Jiang [SIF]

#include <Beam3dPointLoad.h>
#include <Beam3dUniformLoad.h>
#include <Beam3dPartialUniformLoad.h>
#include <BrickSelfWeight.h>
#include <SurfaceLoader.h>
#include <SelfWeight.h>
#include <LoadPattern.h>



#include <SectionForceDeformation.h>
#include <SectionRepres.h>

#include <UniaxialMaterial.h>
#include <LimitCurve.h>
#include <NDMaterial.h>
#include <TclModelBuilder.h>
#include <ImposedMotionSP.h>
#include <ImposedMotionSP1.h>
#include <MultiSupportPattern.h>

#include <TimeSeries.h>
#include <PathTimeSeriesThermal.h>       //L.Jiang [SIF]
#include <vector>						//L.Jiang [SIF] 
using std::vector;							//L.Jiang [SIF]
#include <SimulationInformation.h>				//L.Jiang [SIF]
extern SimulationInformation simulationInfo;		//L.Jiang [SIF]
extern const char * getInterpPWD(Tcl_Interp *interp);  //L.Jiang [SIF]

#include <Block2D.h>
#include <Block3D.h>
// Added by Scott J. Brandenberg (sjbrandenberg@ucdavis.edu)
#include <PySimple1Gen.h>
#include <TzSimple1Gen.h>
// End added by SJB

// Added by Prishati Raychowdhury  (PRC)
#include <ShallowFoundationGen.h>
//end PRC

#include <YieldSurface_BC.h>
#include <YS_Evolution.h>
#include <PlasticHardeningMaterial.h>
#include <CyclicModel.h> //!!
#include <DamageModel.h> //!!

#include <FrictionModel.h>

#ifdef OO_HYSTERETIC
#include <StiffnessDegradation.h>
#include <UnloadingRule.h>
#include <StrengthDegradation.h>
#endif
#include <HystereticBackbone.h>
#include <BeamIntegration.h>

////////////////////// gnp adding damping 
#include <Element.h>
////////////////////////////////////////////

extern void TCL_OPS_setModelBuilder(TclModelBuilder *theNewBuilder);
extern int OPS_ResetInput(ClientData clientData, 
			  Tcl_Interp *interp,  
			  int cArg, 
			  int mArg, 
			  TCL_Char **argv, 
			  Domain *domain,
			  TclModelBuilder *builder);
#include <packages.h>

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

// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

int
TclCommand_addParameter(ClientData clientData, Tcl_Interp *interp, int argc, 
			TCL_Char **argv);

int
TclCommand_addNode(ClientData clientData, Tcl_Interp *interp, int argc, 
		   TCL_Char **argv);

int
TclCommand_addElement(ClientData clientData, Tcl_Interp *interp,  int argc, 
		      TCL_Char **argv);

int
TclCommand_mesh(ClientData clientData, Tcl_Interp *interp,  int argc, 
		TCL_Char **argv);
int
TclCommand_remesh(ClientData clientData, Tcl_Interp *interp,  int argc, 
		  TCL_Char **argv);

int 
TclCommand_backgroundMesh(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
TclCommand_addUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
TclCommand_addBeamIntegration(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
TclCommand_addLimitCurve(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
TclCommand_addNDMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv);

int
TclCommand_addSection(ClientData clientData, Tcl_Interp *interp, int argc,   
		      TCL_Char **argv);

int
TclCommand_addYieldSurface_BC(ClientData clientData, Tcl_Interp *interp,
			      int argc, TCL_Char **argv);

int
TclCommand_addYS_EvolutionModel(ClientData clientData, Tcl_Interp *interp,
				int argc, TCL_Char **argv);

int
TclCommand_addYS_PlasticMaterial(ClientData clientData, Tcl_Interp *interp,
				 int argc, TCL_Char **argv);

int //!!
TclCommand_addCyclicModel(ClientData clientData, Tcl_Interp *interp,
			  int argc, TCL_Char **argv);			    
int //!!
TclCommand_addDamageModel(ClientData clientData, Tcl_Interp *interp,
			  int argc, TCL_Char **argv);	

int
TclCommand_addTimeSeries(ClientData clientData, Tcl_Interp *interp, int argc,
			 TCL_Char **argv);

int
TclCommand_addPattern(ClientData clientData, Tcl_Interp *interp, int argc,   
		      TCL_Char **argv);

int
TclCommand_addSeries(ClientData clientData, Tcl_Interp *interp, int argc,   
		     TCL_Char **argv);

int
TclCommand_addHomogeneousBC(ClientData clientData, Tcl_Interp *interp, int argc,
			    TCL_Char **argv);
int
TclCommand_addHomogeneousBC_X(ClientData clientData, Tcl_Interp *interp, int argc,
			      TCL_Char **argv);
int
TclCommand_addHomogeneousBC_Y(ClientData clientData, Tcl_Interp *interp, int argc,
			      TCL_Char **argv);
int
TclCommand_addHomogeneousBC_Z(ClientData clientData, Tcl_Interp *interp, int argc,
			      TCL_Char **argv);
int
TclCommand_addEqualDOF_MP (ClientData clientData, Tcl_Interp *interp,
			   int argc, TCL_Char **argv);

int
TclCommand_addEqualDOF_MP_Mixed (ClientData clientData, Tcl_Interp *interp,
			   int argc, TCL_Char **argv);

int 
TclCommand_RigidLink(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
TclCommand_RigidDiaphragm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);


int
TclCommand_addMP(ClientData clientData, Tcl_Interp *interp, int argc,   
		 TCL_Char **argv);

int
TclCommand_addNodalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
			TCL_Char **argv);

int
TclCommand_addElementalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
			    TCL_Char **argv);

int
TclCommand_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc,   
			TCL_Char **argv);
int
TclCommand_addSP(ClientData clientData, Tcl_Interp *interp, int argc,   
		      TCL_Char **argv);

int
TclCommand_addImposedMotionSP(ClientData clientData, 
			      Tcl_Interp *interp, 
			      int argc,    
			      TCL_Char **argv);	
// Added by Scott J. Brandenberg
int
TclCommand_doPySimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char **argv);

int
TclCommand_doTzSimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char **argv);
// End added by SJB

// Added by Prishati Raychowdhury (UCSD)
int
TclModelBuilder_doShallowFoundationGen(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char **argv);
// End PRC

int
TclCommand_doBlock2D(ClientData clientData, Tcl_Interp *interp, int argc, 
		     TCL_Char **argv);

int
TclCommand_doBlock3D(ClientData clientData, Tcl_Interp *interp, int argc, 
		     TCL_Char **argv);

int
TclCommand_addRemoPatch(ClientData clientData, 
			Tcl_Interp *interp, 
			int argc,   
			TCL_Char **argv);  

int
TclCommand_addRemoLayer(ClientData clientData, 
			Tcl_Interp *interp, 
			int argc,   
			TCL_Char **argv);   

int
TclCommand_addRemoFiber(ClientData clientData, 
			Tcl_Interp *interp, 
			int argc,    
			TCL_Char **argv);   


//Leo
int
TclModelBuilder_addRemoHFiber(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc,    
			     TCL_Char **argv);  

int
TclCommand_addRemoGeomTransf(ClientData clientData, 
			     Tcl_Interp *interp, 
			     int argc,   
			     TCL_Char **argv); 

int
TclCommand_addFrictionModel(ClientData clientData,
				 Tcl_Interp *interp,
				 int argc,   
				 TCL_Char **argv);

#ifdef OO_HYSTERETIC
int
TclCommand_addStiffnessDegradation(ClientData clientData,
				   Tcl_Interp *interp,
				   int argc, TCL_Char **argv);

int
TclCommand_addUnloadingRule(ClientData clientData,
			    Tcl_Interp *interp,
			    int argc, TCL_Char **argv);

int
TclCommand_addStrengthDegradation(ClientData clientData,
				  Tcl_Interp *interp,
				  int argc, TCL_Char **argv);
#endif

int
TclCommand_addHystereticBackbone(ClientData clientData,
				 Tcl_Interp *interp,
				 int argc, TCL_Char **argv);

int
TclCommand_addGroundMotion(ClientData clientData, 
			   Tcl_Interp *interp, 
			   int argc,    
			   TCL_Char **argv);

/// added by ZHY
int
TclCommand_UpdateMaterialStage(ClientData clientData, 
			       Tcl_Interp *interp,  
			       int argc,
			       TCL_Char **argv);
int
TclCommand_UpdateMaterials(ClientData clientData, 
			   Tcl_Interp *interp,  
			   int argc,
			   TCL_Char **argv);


/// added by ZHY			   
int
TclCommand_UpdateParameter(ClientData clientData, 
			   Tcl_Interp *interp,  
			   int argc, 
			   TCL_Char **argv);

////////////////gnp adding rayleigh //////////////////////////
int 
TclCommand_addElementRayleigh(ClientData clientData, 
			      Tcl_Interp *interp,  
			      int argc, 
			      TCL_Char **argv);
///////////////////////////////////////////////////////////////



// REMO
extern int
TclCommand_addPatch (ClientData clientData, Tcl_Interp *interp,
		     int argc, TCL_Char **argv,
		     TclModelBuilder *theTclBuilder);


extern int
TclCommand_addFiber (ClientData clientData, Tcl_Interp *interp,
		     int argc, TCL_Char **argv,
		     TclModelBuilder *theTclBuilder);


extern int
TclCommand_addHFiber (ClientData clientData, Tcl_Interp *interp,
			  int argc, TCL_Char **argv,
			  TclModelBuilder *theTclBuilder);


extern int
TclCommand_addReinfLayer (ClientData clientData, Tcl_Interp *interp,
			  int argc, TCL_Char **argv,
			  TclModelBuilder *theTclBuilder);


extern int
TclCommand_addGeomTransf(ClientData, Tcl_Interp *, int, TCL_Char **,
			 Domain*, TclModelBuilder *);   


int
TclCommand_Package(ClientData clientData, Tcl_Interp *interp, int argc, 
		   TCL_Char **argv);
		   
// Added by Alborz Ghofrani - U.Washington
int
TclCommand_GenerateInterfacePoints(ClientData clientData, Tcl_Interp *interp, int argc,
  		                   TCL_Char **argv);
// End Added by Alborz				  

//
// CLASS CONSTRUCTOR & DESTRUCTOR
//

// constructor: the constructor will add certain commands to the interpreter
TclModelBuilder::TclModelBuilder(Domain &theDomain, Tcl_Interp *interp, int NDM, int NDF)
  :ModelBuilder(theDomain), ndm(NDM), ndf(NDF), theInterp(interp)
{
  theSections  = new ArrayOfTaggedObjects(32);
  theSectionRepresents = new ArrayOfTaggedObjects(32);  
#ifdef OO_HYSTERETIC
  theStiffnessDegradations = new ArrayOfTaggedObjects(32);
  theUnloadingRules = new ArrayOfTaggedObjects(32);
  theStrengthDegradations = new ArrayOfTaggedObjects(32);
#endif
  //theHystereticBackbones= new ArrayOfTaggedObjects(32);
  theYieldSurface_BCs = new ArrayOfTaggedObjects(32);
  theCycModels = new ArrayOfTaggedObjects(32); //!!
  theYS_EvolutionModels = new ArrayOfTaggedObjects(32);
  thePlasticMaterials = new ArrayOfTaggedObjects(32);

  // call Tcl_CreateCommand for class specific commands
  Tcl_CreateCommand(interp, "parameter", TclCommand_addParameter,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "addToParameter", TclCommand_addParameter,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "updateParameter", TclCommand_addParameter,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "node", TclCommand_addNode,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "element", TclCommand_addElement,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "mesh", TclCommand_mesh,
		    (ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "remesh", TclCommand_remesh,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "background", &TclCommand_backgroundMesh, 
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);


  Tcl_CreateCommand(interp, "uniaxialMaterial", TclCommand_addUniaxialMaterial,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "beamIntegration", TclCommand_addBeamIntegration,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "limitCurve", TclCommand_addLimitCurve,
		    (ClientData)NULL, NULL);
  
  Tcl_CreateCommand(interp, "nDMaterial", TclCommand_addNDMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "section", TclCommand_addSection,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "yieldSurface_BC", TclCommand_addYieldSurface_BC,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "ysEvolutionModel", TclCommand_addYS_EvolutionModel,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "plasticMaterial", TclCommand_addYS_PlasticMaterial,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "cyclicModel", TclCommand_addCyclicModel,
		    (ClientData)NULL, NULL); //!!

  Tcl_CreateCommand(interp, "damageModel", TclCommand_addDamageModel,
		    (ClientData)NULL, NULL); //!!

  Tcl_CreateCommand(interp, "pattern", TclCommand_addPattern,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "timeSeries", TclCommand_addTimeSeries,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "load", TclCommand_addNodalLoad,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "eleLoad", TclCommand_addElementalLoad,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "mass", TclCommand_addNodalMass,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fix", TclCommand_addHomogeneousBC,
  		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fixX", TclCommand_addHomogeneousBC_X,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fixY", TclCommand_addHomogeneousBC_Y,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "fixZ", TclCommand_addHomogeneousBC_Z,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "sp", TclCommand_addSP,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "imposedMotion", 
		    TclCommand_addImposedMotionSP,
		    (ClientData)NULL, NULL);  

  Tcl_CreateCommand(interp, "imposedSupportMotion", 
		    TclCommand_addImposedMotionSP,
		    (ClientData)NULL, NULL);  
  
  Tcl_CreateCommand(interp, "groundMotion", 
		    TclCommand_addGroundMotion,
		    (ClientData)NULL, NULL);    

  Tcl_CreateCommand(interp, "equalDOF", TclCommand_addEqualDOF_MP,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "equalDOF_Mixed", TclCommand_addEqualDOF_MP_Mixed,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "rigidLink", &TclCommand_RigidLink, 
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);                

  Tcl_CreateCommand(interp, "rigidDiaphragm", &TclCommand_RigidDiaphragm, 
		    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);   

  Tcl_CreateCommand(interp, "mp", TclCommand_addMP,
		    (ClientData)NULL, NULL);

  // Added by Scott J. Brandenberg
  Tcl_CreateCommand(interp, "PySimple1Gen", TclCommand_doPySimple1Gen,
                        (ClientData)NULL, NULL);
  Tcl_CreateCommand(interp, "TzSimple1Gen", TclCommand_doTzSimple1Gen,
                        (ClientData)NULL, NULL);
  // End added by SJB

// Added by Prishati Raychowdhury (UCSD)
  Tcl_CreateCommand(interp, "ShallowFoundationGen", TclModelBuilder_doShallowFoundationGen,
                        (ClientData)NULL, NULL);
// End PRC

  Tcl_CreateCommand(interp, "block2D", TclCommand_doBlock2D,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "block3D", TclCommand_doBlock3D,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "patch", TclCommand_addRemoPatch,
		    (ClientData)NULL, NULL);  

  Tcl_CreateCommand(interp, "layer", TclCommand_addRemoLayer,
		    (ClientData)NULL, NULL);    
  
  Tcl_CreateCommand(interp, "fiber", TclCommand_addRemoFiber,
		    (ClientData)NULL, NULL);    

  //LEO
  Tcl_CreateCommand(interp, "Hfiber", TclModelBuilder_addRemoHFiber,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "geomTransf", TclCommand_addRemoGeomTransf,
		    (ClientData)NULL, NULL);    

  Tcl_CreateCommand(interp, "frictionModel",
		    TclCommand_addFrictionModel,
		    (ClientData)NULL, NULL);

#ifdef OO_HYSTERETIC
  Tcl_CreateCommand(interp, "stiffnessDegradation",
		    TclCommand_addStiffnessDegradation,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "unloadingRule",
		    TclCommand_addUnloadingRule,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "strengthDegradation",
		    TclCommand_addStrengthDegradation,
		    (ClientData)NULL, NULL);
#endif
  Tcl_CreateCommand(interp, "hystereticBackbone",
		    TclCommand_addHystereticBackbone,
		    (ClientData)NULL, NULL);

  ///new command for elast2plast in Multi-yield plasticity, by ZHY
  Tcl_CreateCommand(interp, "updateMaterialStage", 
		    TclCommand_UpdateMaterialStage,
		    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "updateMaterials", 
		    TclCommand_UpdateMaterials,
		    (ClientData)NULL, NULL);
  
  ///new command for updating properties of soil materials, by ZHY
  //Tcl_CreateCommand(interp, "updateParameter", 
  //		    TclCommand_UpdateParameter,
  //	    (ClientData)NULL, NULL);

  Tcl_CreateCommand(interp, "loadPackage", TclCommand_Package,
		    (ClientData)NULL, NULL);


  ////// gnp adding per element damping ///////////////////////////////
  Tcl_CreateCommand(interp, "setElementRayleighFactors",
		    TclCommand_addElementRayleigh,
		    (ClientData)NULL, NULL);
  /////////////////////////////////////////////////////////////////////

  // Added by Alborz Ghofrani - U.Washington
  Tcl_CreateCommand(interp, "generateInterfacePoints",
                    TclCommand_GenerateInterfacePoints,
                    (ClientData)NULL, NULL);
  // End Added by Alborz
  // set the static pointers in this file
  theTclBuilder = this;
  theTclDomain = &theDomain;
  theTclLoadPattern = 0;
  theTclMultiSupportPattern = 0;  

  TCL_OPS_setModelBuilder(this);

  nodeLoadTag = 0;
  eleArgStart = 0;
}

TclModelBuilder::~TclModelBuilder()
{
  OPS_clearAllTimeSeries();
  //  OPS_clearAllUniaxialMaterial();
  //  OPS_clearAllNDMaterial();
  //  OPS_clearAllSectionForceDeformation();
  OPS_clearAllCrdTransf();
  OPS_clearAllFrictionModel();
  OPS_clearAllLimitCurve();
  OPS_clearAllDamageModel();
  OPS_clearAllFrictionModel();
  OPS_clearAllHystereticBackbone();
  //  OPS_clearAllNDMaterial();

  theSections->clearAll(); 
  theSectionRepresents->clearAll();
  theYieldSurface_BCs->clearAll();
  theYS_EvolutionModels->clearAll();
  thePlasticMaterials->clearAll();
  theCycModels->clearAll();//!!

#ifdef OO_HYSTERETIC
  theStiffnessDegradations->clearAll();
  theUnloadingRules->clearAll();
  theStrengthDegradations->clearAll();
#endif
  // theHystereticBackbones->clearAll();

  // free up memory allocated in the constructor
  delete theSections;
  delete theSectionRepresents;
  delete theYieldSurface_BCs;
  delete theYS_EvolutionModels;
  delete thePlasticMaterials;
  delete theCycModels;//!!

#ifdef OO_HYSTERETIC
  delete theStiffnessDegradations;
  delete theUnloadingRules;
  delete theStrengthDegradations;
#endif
  // delete theHystereticBackbones;

  // set the pointers to 0 
  theTclDomain =0;
  theTclBuilder =0;
  theTclLoadPattern =0;
  theTclMultiSupportPattern = 0;  
  TCL_OPS_setModelBuilder(0);
  
  // may possibly invoke Tcl_DeleteCommand() later
  Tcl_DeleteCommand(theInterp, "parameter");
  Tcl_DeleteCommand(theInterp, "addToParameter");
  Tcl_DeleteCommand(theInterp, "updateParameter");
  Tcl_DeleteCommand(theInterp, "node");
  Tcl_DeleteCommand(theInterp, "element");
  Tcl_DeleteCommand(theInterp, "mesh");
  Tcl_DeleteCommand(theInterp, "remesh");
  Tcl_DeleteCommand(theInterp, "background");
  Tcl_DeleteCommand(theInterp, "uniaxialMaterial");
  Tcl_DeleteCommand(theInterp, "nDMaterial");
  Tcl_DeleteCommand(theInterp, "section");
  Tcl_DeleteCommand(theInterp, "pattern");
  Tcl_DeleteCommand(theInterp, "timeSeries");
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
  Tcl_DeleteCommand(theInterp, "Hfiber"); //LEO
  Tcl_DeleteCommand(theInterp, "geomTransf");
  Tcl_DeleteCommand(theInterp, "updateMaterialStage");
  Tcl_DeleteCommand(theInterp, "updateMaterials");

  Tcl_DeleteCommand(theInterp, "frictionModel");

#ifdef OO_HYSTERETIC
  Tcl_DeleteCommand(theInterp, "unloadingRule");
  Tcl_DeleteCommand(theInterp, "stiffnessDegradation");
  Tcl_DeleteCommand(theInterp, "strengthDegradation");
#endif
  Tcl_DeleteCommand(theInterp, "hystereticBackbone");

  Tcl_DeleteCommand(theInterp, "yieldSurface_BC");
  Tcl_DeleteCommand(theInterp, "ysEvolutionModel");
  Tcl_DeleteCommand(theInterp, "plasticMaterial");
  Tcl_DeleteCommand(theInterp, "cyclicModel");
  Tcl_DeleteCommand(theInterp, "damageModel");

  Tcl_DeleteCommand(theInterp, "loadPackage");
    Tcl_DeleteCommand(theInterp, "generateInterfacePoints"); // Added by Alborz Ghofrani - U.Washington
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

#ifdef OO_HYSTERETIC
int
TclModelBuilder::addStiffnessDegradation(StiffnessDegradation &theDegr)
{
  bool result = theStiffnessDegradations->addComponent(&theDegr);
  if (result == true)
    return 0;
  else {
    opserr << "TclModelBuilder::addStiffnessDegradation() - failed to add StiffnessDegradation: " << theDegr;
    return -1;
  }
}

StiffnessDegradation*
TclModelBuilder::getStiffnessDegradation(int tag)
{
  TaggedObject *mc = theStiffnessDegradations->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // do a cast and return
  StiffnessDegradation *result = (StiffnessDegradation *)mc;
  return result;
}

int
TclModelBuilder::addUnloadingRule(UnloadingRule &theDegr)
{
  bool result = theUnloadingRules->addComponent(&theDegr);
  if (result == true)
    return 0;
  else {
    opserr << "TclModelBuilder::addUnloadingRule() - failed to add UnloadingRule: " << theDegr;
    return -1;
  }
}

UnloadingRule*
TclModelBuilder::getUnloadingRule(int tag)
{
  TaggedObject *mc = theUnloadingRules->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // do a cast and return
  UnloadingRule *result = (UnloadingRule *)mc;
  return result;
}

int
TclModelBuilder::addStrengthDegradation(StrengthDegradation &theDegr)
{
  bool result = theStrengthDegradations->addComponent(&theDegr);
  if (result == true)
    return 0;
  else {
    opserr << "TclModelBuilder::addStrengthDegradation() - failed to add StrengthDegradation: " << theDegr;
    return -1;
  }
}

StrengthDegradation*
TclModelBuilder::getStrengthDegradation(int tag)
{
  TaggedObject *mc = theStrengthDegradations->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // do a cast and return
  StrengthDegradation *result = (StrengthDegradation *)mc;
  return result;
}
#endif

/*
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
*/

int 
TclModelBuilder::addSection(SectionForceDeformation &theSection)
{
  //  bool result = theSections->addComponent(&theSection);
  bool result = OPS_addSectionForceDeformation(&theSection);
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
  return OPS_getSectionForceDeformation(tag);
  /*
  TaggedObject *mc = theSections->getComponentPtr(tag);
  if (mc == 0) 
    return 0;

  // do a cast and return
  SectionForceDeformation *result = (SectionForceDeformation *)mc;
  return result;
  */
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
TclCommand_addNode(ClientData clientData, Tcl_Interp *interp, int argc, 
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
    //    theNode = new Node(nodeId,ndf,xLoc);
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
    //    theNode = new Node(nodeId,ndf,xLoc,yLoc);
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
    //    theNode = new Node(nodeId,ndf,xLoc,yLoc,zLoc);
  } else {
      opserr << "WARNING invalid ndm\n";
      opserr << "node: " << nodeId << endln;;
      return TCL_ERROR;
  }

  // check for -ndf override option
  int currentArg = 2+ndm;  
  if (currentArg < argc && strcmp(argv[currentArg],"-ndf") == 0) {
    if (Tcl_GetInt(interp, argv[currentArg+1], &ndf) != TCL_OK) {
      opserr << "WARNING invalid nodal ndf given for node " << nodeId << endln;
      return TCL_ERROR;
    }
    currentArg += 2;
  }

  //
  // create the node
  //

  if (ndm == 1) 
    theNode = new Node(nodeId,ndf,xLoc);        
  else if (ndm == 2)
    theNode = new Node(nodeId,ndf,xLoc,yLoc);    
  else
    theNode = new Node(nodeId,ndf,xLoc,yLoc,zLoc);    

  //
  // add the node to the domain
  //

  if (theTclDomain->addNode(theNode) == false) {
    opserr << "WARNING failed to add node to the domain\n";
    opserr << "node: " << nodeId << endln;
    delete theNode; // otherwise memory leak                                    
    return TCL_ERROR;
  }

  while (currentArg < argc) {
    if (strcmp(argv[currentArg],"-mass") == 0) {
      currentArg++;
      if (argc < currentArg+ndf) {
	opserr << "WARNING incorrect number of nodal mass terms\n";
	opserr << "node: " << nodeId << endln;
	return TCL_ERROR;      
      }	
      Matrix mass(ndf,ndf);
      double theMass;
      for (int i=0; i<ndf; i++) {
	if (Tcl_GetDouble(interp, argv[currentArg++], &theMass) != TCL_OK) {
	  opserr << "WARNING invalid nodal mass term\n";
	  opserr << "node: " << nodeId << ", dof: " << i+1 << endln;
	  return TCL_ERROR;
	}
	mass(i,i) = theMass;
      }
      theNode->setMass(mass);      
    } else if (strcmp(argv[currentArg],"-dispLoc") == 0) {
      currentArg++;
      if (argc < currentArg+ndm) {
	opserr << "WARNING incorrect number of nodal display location terms, need ndm\n";
	opserr << "node: " << nodeId << endln;
	return TCL_ERROR;      
      }	
      Vector displayLoc(ndm);
      double theCrd;
      for (int i=0; i<ndm; i++) {
	if (Tcl_GetDouble(interp, argv[currentArg++], &theCrd) != TCL_OK) {
	  opserr << "WARNING invalid nodal mass term\n";
	  opserr << "node: " << nodeId << ", dof: " << i+1 << endln;
	  return TCL_ERROR;
	}
	displayLoc(i) = theCrd;
      }
      theNode->setDisplayCrds(displayLoc);

    } else if (strcmp(argv[currentArg],"-disp") == 0) {
      currentArg++;
      if (argc < currentArg+ndf) {
	opserr << "WARNING incorrect number of nodal disp terms\n";
	opserr << "node: " << nodeId << endln;
	return TCL_ERROR;      
      }	
      Vector disp(ndf);
      double theDisp;
      for (int i=0; i<ndf; i++) {
	if (Tcl_GetDouble(interp, argv[currentArg++], &theDisp) != TCL_OK) {
	  opserr << "WARNING invalid nodal disp term\n";
	  opserr << "node: " << nodeId << ", dof: " << i+1 << endln;
	  return TCL_ERROR;
	}
	disp(i) = theDisp;
      }
      theNode->setTrialDisp(disp);      
      theNode->commitState();

    } else if (strcmp(argv[currentArg],"-vel") == 0) {
      currentArg++;
      if (argc < currentArg+ndf) {
	opserr << "WARNING incorrect number of nodal vel terms\n";
	opserr << "node: " << nodeId << endln;
	return TCL_ERROR;      
      }	
      Vector disp(ndf);
      double theDisp;
      for (int i=0; i<ndf; i++) {
	if (Tcl_GetDouble(interp, argv[currentArg++], &theDisp) != TCL_OK) {
	  opserr << "WARNING invalid nodal vel term\n";
	  opserr << "node: " << nodeId << ", dof: " << i+1 << endln;
	  return TCL_ERROR;
	}
	disp(i) = theDisp;
      }
      theNode->setTrialVel(disp); 
      theNode->commitState();

    } else
      currentArg++;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



/////////////////////////////   gnp adding element damping 
int 
TclCommand_addElementRayleigh(ClientData clientData, 
			      Tcl_Interp *interp,  
			      int argc, 
			      TCL_Char **argv) 
{
  
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed" << endln;
    return TCL_ERROR;
  }

  // make sure corect number of arguments on command line
  if (argc < 6) {
    opserr << "WARNING insufficient arguments\n";
    printCommand(argc, argv);
    opserr << "Want: setElementRayleighFactors elementTag?  alphaM? $betaK? $betaKinit? $betaKcomm? \n";
    return TCL_ERROR;
  }    
  
  int eleTag =0;
  
  if (Tcl_GetInt(interp, argv[1], &eleTag) != TCL_OK) {
    opserr << "WARNING: setElementRayleighFactors invalid eleTag: " << argv[1];
    opserr << " \n";
    return TCL_ERROR;
  }
  
  double alphaM,betaK,betaKinit,betaKcomm;
  
  if (Tcl_GetDouble(interp, argv[2], &alphaM) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "alphaM: " << argv[2] << endln;
    return TCL_ERROR;
  }
  
  if (Tcl_GetDouble(interp, argv[3], &betaK) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "betaK: " << argv[3] << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &betaKinit) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "betaKinit: " << argv[4] << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &betaKcomm) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "betaKcomm: " << argv[5] << endln;
    return TCL_ERROR;
  }
  
  Element* elePtr = theTclDomain->getElement(eleTag);
  
  if (elePtr == 0) 
    opserr << "WARNING : setElementRayleighFactors invalid eleTag: " << eleTag << " the element does not exist in the domain \n";
  
  
  if ( elePtr->setRayleighDampingFactors(alphaM, betaK, betaKinit, betaKcomm) != 0 ) {
    opserr << "ERROR : setElementRayleighFactors: FAILED to add damping factors for element " << eleTag << "\n";
    
  }
  
  return TCL_OK;
}
/////////////////////////////   gnp adding element damping 


// the function for creating ne material objects and patterns is in a seperate file.
// this allows new material and patternobjects to be added without touching this file.
// does so at the expense of an extra procedure call.

extern int 
TclModelBuilderParameterCommand(ClientData clientData, 
				Tcl_Interp *interp, int argc,    
				TCL_Char **argv, 
				Domain *theDomain, TclModelBuilder *theTclBuilder);
int
TclCommand_addParameter(ClientData clientData, Tcl_Interp *interp, 
			     int argc, TCL_Char **argv)
                          
{
  return TclModelBuilderParameterCommand(clientData, interp, 
					 argc, argv, theTclDomain, theTclBuilder);
}


extern int 
TclModelBuilderElementCommand(ClientData clientData, 
			      Tcl_Interp *interp, int argc,    
			      TCL_Char **argv, 
			      Domain *theDomain, TclModelBuilder *theTclBuilder);
int
TclCommand_addElement(ClientData clientData, Tcl_Interp *interp, 
			   int argc,    TCL_Char **argv)
                          
{
  return TclModelBuilderElementCommand(clientData, interp, 
				       argc, argv, theTclDomain, theTclBuilder);
}

// extern int OPS_LineMesh(Domain& domain, int ndm);
// extern int OPS_TriMesh(Domain& domain);
// extern int OPS_TriReMesh(Domain& domain, int ndf);
int
TclCommand_mesh(ClientData clientData, Tcl_Interp *interp,  int argc, 
		TCL_Char **argv) 
{
    // ensure the destructor has not been called - 
    if (theTclBuilder == 0) {
	opserr << "WARNING builder has been destroyed" << endln;
	return TCL_ERROR;
    }

    int ndm = theTclBuilder->getNDM();

    // make sure corect number of arguments on command line
    if (argc < 2) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: mesh type? ...>\n";
	return TCL_ERROR;
    }

    OPS_ResetInput(clientData,interp,2,argc,argv,theTclDomain,theTclBuilder);

    // mesh type
    int res = 0;
    // if (strcmp(argv[1], "line") == 0) {
    // 	res = OPS_LineMesh(*theTclDomain,ndm);
    // } else if (strcmp(argv[1], "tri") == 0) {
    // 	res = OPS_TriMesh(*theTclDomain);
    // } else {
    // 	opserr<<"WARNING: mesh type "<<argv[1]<<" is unknown\n";
    // 	return TCL_ERROR;
    // }

    if (res < 0) {
	return TCL_ERROR;
    }

    return 0;

}

int
TclCommand_remesh(ClientData clientData, Tcl_Interp *interp,  int argc, 
		TCL_Char **argv) 
{
    // ensure the destructor has not been called - 
    if (theTclBuilder == 0) {
	opserr << "WARNING builder has been destroyed" << endln;
	return TCL_ERROR;
    }

    int ndf = theTclBuilder->getNDF();

    // make sure corect number of arguments on command line
    if (argc < 2) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: mesh type? ...>\n";
	return TCL_ERROR;
    }

    OPS_ResetInput(clientData,interp,2,argc,argv,theTclDomain,theTclBuilder);

    // mesh type
    int res = 0;
    // if (strcmp(argv[1], "line") == 0) {
    // 	//res = OPS_LineMesh(*theTclDomain,ndm);
    // } else if (strcmp(argv[1], "tri") == 0) {
    // 	res = OPS_TriReMesh(*theTclDomain,ndf);
    // } else {
    // 	opserr<<"WARNING: remesh type "<<argv[1]<<" is unknown\n";
    // 	return TCL_ERROR;
    // }

    if (res < 0) {
	return TCL_ERROR;
    }
    
    return 0;

}

extern int OPS_BgMesh();

int 
TclCommand_backgroundMesh(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // ensure the destructor has not been called - 
    if (theTclBuilder == 0) {
	opserr << "WARNING builder has been destroyed" << endln;
	return TCL_ERROR;
    }
    
    OPS_ResetInput(clientData, interp, 1, argc, argv, theTclDomain, theTclBuilder);

    if(OPS_BgMesh() >= 0) return TCL_OK;
    else return TCL_ERROR;
    return TCL_OK;
}

extern void* OPS_LobattoBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_LegendreBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewtonCotesBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_RadauBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_TrapezoidalBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_CompositeSimpsonBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_UserDefinedBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_FixedLocationBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_LowOrderBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_MidDistanceBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_UserHingeBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_HingeMidpointBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_HingeRadauBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_HingeRadauTwoBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_HingeEndpointBeamIntegration(int& integrationTag, ID& secTags);

int
TclCommand_addBeamIntegration(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    if (argc < 2) {
	opserr << "WARNING: want beamIntegration type itag...\n";
	return TCL_ERROR;
    }

    OPS_ResetInput(clientData,interp,2,argc,argv,theTclDomain,theTclBuilder);
    
    int iTag;
    ID secTags;
    BeamIntegration* bi = 0;
    if (strcmp(argv[1],"Lobatto") == 0) {
	bi = (BeamIntegration*)OPS_LobattoBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"Legendre") == 0) {
	bi = (BeamIntegration*)OPS_LegendreBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"NewtoCotes") == 0) {
	bi = (BeamIntegration*)OPS_NewtonCotesBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"Radau") == 0) {
	bi = (BeamIntegration*)OPS_RadauBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"Trapezoidal") == 0) {
	bi = (BeamIntegration*)OPS_TrapezoidalBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"CompositeSimpson") == 0) {
	bi = (BeamIntegration*)OPS_CompositeSimpsonBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"UserDefined") == 0) {
	bi = (BeamIntegration*)OPS_UserDefinedBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"FixedLocation") == 0) {
	bi = (BeamIntegration*)OPS_FixedLocationBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"LowOrder") == 0) {
	bi = (BeamIntegration*)OPS_LowOrderBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"MidDistance") == 0) {
	bi = (BeamIntegration*)OPS_MidDistanceBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"UserHinge") == 0) {
	bi = (BeamIntegration*)OPS_UserHingeBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"HingeMidpoint") == 0) {
	bi = (BeamIntegration*)OPS_HingeMidpointBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"HingeRadau") == 0) {
	bi = (BeamIntegration*)OPS_HingeRadauBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"HingeRadauTwo") == 0) {
	bi = (BeamIntegration*)OPS_HingeRadauTwoBeamIntegration(iTag,secTags);
    } else if (strcmp(argv[1],"HingeEndpoint") == 0) {
	bi = (BeamIntegration*)OPS_HingeEndpointBeamIntegration(iTag,secTags);
    } else {
	opserr<<"WARNING: integration type "<<argv[1]<<" is unknown\n";
	return TCL_ERROR;
    }

    if (bi == 0) {
	opserr<<"WARNING: failed to create beam integration\n";
	return TCL_ERROR;
    }

    BeamIntegrationRule* rule = new BeamIntegrationRule(iTag,bi,secTags);
    if (rule == 0) {
	opserr<<"WARNING: failed to create beam integration\n";
	delete bi;
	return TCL_ERROR;
    }

    // Now add the 
    if(OPS_addBeamIntegrationRule(rule) == false) {
	opserr<<"WARNING: could not add BeamIntegrationRule.";
	delete rule; // invoke the destructor, otherwise mem leak
	return TCL_ERROR;;
    }

    return TCL_OK;
}

extern int
TclModelBuilderUniaxialMaterialCommand (ClientData clienData, Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theDomain);
				 
int
TclCommand_addUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  return TclModelBuilderUniaxialMaterialCommand(clientData, interp, argc, argv, theTclDomain);
}


extern int
Tcl_AddLimitCurveCommand (ClientData clienData, Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theDomain);

int
TclCommand_addLimitCurve(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  return Tcl_AddLimitCurveCommand(clientData, interp, argc, argv, theTclDomain);
}

extern int
TclModelBuilderNDMaterialCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				  TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclCommand_addNDMaterial(ClientData clientData, Tcl_Interp *interp, 
			    int argc,    TCL_Char **argv)
                          
{
  return TclModelBuilderNDMaterialCommand(clientData, interp, 
					  argc, argv, theTclBuilder);
}

extern int
TclModelBuilderSectionCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				  TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclCommand_addSection(ClientData clientData, Tcl_Interp *interp, 
			    int argc,    TCL_Char **argv)
                          
{
  return TclModelBuilderSectionCommand(clientData, interp, 
				       argc, argv, theTclBuilder);
}



extern int
TclModelBuilderYieldSurface_BCCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclCommand_addYieldSurface_BC(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv)

{
  return TclModelBuilderYieldSurface_BCCommand(clientData, interp,
						argc, argv, theTclBuilder);
}

extern int
TclModelBuilderYS_EvolutionModelCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclCommand_addYS_EvolutionModel(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv)

{
  return TclModelBuilderYS_EvolutionModelCommand(clientData, interp,
						argc, argv, theTclBuilder);
}

extern int
TclModelBuilderPlasticMaterialCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclCommand_addYS_PlasticMaterial(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv)

{
  return TclModelBuilderPlasticMaterialCommand(clientData, interp,
						argc, argv, theTclBuilder);
}

//!!
extern int TclModelBuilderCyclicModelCommand(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclBuilder);
int
TclCommand_addCyclicModel(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv)

{
  return TclModelBuilderCyclicModelCommand(clientData, interp,
						argc, argv, theTclBuilder);
}

extern int TclModelBuilderDamageModelCommand(ClientData clienData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
TclCommand_addDamageModel(ClientData clientData, Tcl_Interp *interp,
				    int argc, TCL_Char **argv)

{
  return TclModelBuilderDamageModelCommand(clientData, interp, argc, argv);
						
}

extern int
TclPatternCommand(ClientData clientData, Tcl_Interp *interp, 
			   int argc, TCL_Char **argv, Domain *theDomain);
			   
int
TclCommand_addPattern(ClientData clientData, Tcl_Interp *interp, 
			   int argc, TCL_Char **argv)
{
  return TclPatternCommand(clientData, interp, argc, argv, theTclDomain);
}


extern TimeSeries *
TclTimeSeriesCommand(ClientData clientData, Tcl_Interp *interp, 
		     int argc, TCL_Char **argv, Domain *theDomain);

int
TclCommand_addTimeSeries(ClientData clientData, Tcl_Interp *interp, 
			 int argc, TCL_Char **argv)
{
  TimeSeries *theSeries = TclTimeSeriesCommand(clientData, interp, argc-1, &argv[1], 0);

  if (theSeries != 0) {
    if (OPS_addTimeSeries(theSeries) == true)
      return TCL_OK;
    else
      return TCL_ERROR;
  }
  return TCL_ERROR;
}




extern int
TclGroundMotionCommand(ClientData clientData, 
		       Tcl_Interp *interp, 
		       int argc,    
		       TCL_Char **argv,
		       MultiSupportPattern *thePattern);

int
TclCommand_addGroundMotion(ClientData clientData, Tcl_Interp *interp, 
			   int argc, TCL_Char **argv)
			  
{
  return TclGroundMotionCommand(clientData, interp, argc, argv, 
				theTclMultiSupportPattern);
}


int
TclCommand_addNodalLoad(ClientData clientData, Tcl_Interp *interp, int argc,
	TCL_Char **argv)
{
	// ensure the destructor has not been called - 
	if (theTclBuilder == 0) {
		opserr << "WARNING builder has been destroyed - load \n";
		return TCL_ERROR;
	}

	//  int ndf = theTclBuilder->getNDF();

	int ndf = argc - 2;
	NodalLoad *theLoad = 0;

	bool isLoadConst = false;
	bool userSpecifiedPattern = false;
	int loadPatternTag = 0;
	//The above definition are moved forward for the use in both cases

	//-------------Adding Proc for NodalThermalAction, By Liming Jiang, [SIF] 2017
	if ((strcmp(argv[2], "-NodalThermal") == 0) || (strcmp(argv[2], "-nodalThermal") == 0)) {

		int nodeId;
		if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
			opserr << "WARNING invalid nodeId: " << argv[1] << endln;
			return TCL_ERROR;
		}

		Vector* thecrds = new Vector();
		Node* theNode = theTclDomain->getNode(nodeId);
		if (theNode == 0) {
			opserr << "WARNING invalid nodeID: " << argv[1] << endln;
			return TCL_ERROR;
		}
		(*thecrds) = theNode->getCrds();

		int count = 3;
		if (strcmp(argv[count], "-source") == 0) {
			count++;
			const char *pwd = getInterpPWD(interp);
			simulationInfo.addInputFile(argv[count], pwd);
			TimeSeries* theSeries;

			int dataLen = 9;//default num of temperature input for nodal ThermalAction;

			if (argc - count == 5) {
				//which indicates the nodal thermal action is applied to 3D I section Beam;
				dataLen = 15;
				theSeries = new PathTimeSeriesThermal(nodeId, argv[count], dataLen);
				count++;
				double RcvLoc1, RcvLoc2, RcvLoc3, RcvLoc4;
				if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
					opserr << "WARNING NodalLoad - invalid loc1  " << argv[count] << " for NodalThermalAction\n";
					return TCL_ERROR;
				}
				if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
					opserr << "WARNING NodalLoad - invalid loc2  " << argv[count + 1] << " for NodalThermalAction\n";
					return TCL_ERROR;
				}
				if (Tcl_GetDouble(interp, argv[count + 2], &RcvLoc3) != TCL_OK) {
					opserr << "WARNING NodalLoad - invalid loc3  " << argv[count + 2] << " for NodalThermalAction\n";
					return TCL_ERROR;
				}
				if (Tcl_GetDouble(interp, argv[count + 3], &RcvLoc4) != TCL_OK) {
					opserr << "WARNING NodalLoad - invalid loc4  " << argv[count + 3] << " for NodalThermalAction\n";
					return TCL_ERROR;
				}
				//end of recieving data;
				theLoad = new NodalThermalAction(nodeLoadTag, nodeId, RcvLoc1, RcvLoc2, RcvLoc3, RcvLoc4, theSeries, thecrds);
			}
			//end of for 15 data input;
			else if (argc - count == 3 || argc - count == 10) {

				theSeries = new PathTimeSeriesThermal(nodeId, argv[count]);
				count++;
				Vector locy;
				if (argc - count == 2) {
					double RcvLoc1, RcvLoc2;
					if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
						opserr << "WARNING NodalLoad - invalid loc1  " << argv[count] << " for NodalThermalAction\n";
						return TCL_ERROR;
					}
					if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
						opserr << "WARNING NodalLoad - invalid loc2  " << argv[count + 1] << " for NodalThermalAction\n";
						return TCL_ERROR;
					}
					locy = Vector(9);
					locy(0) = RcvLoc1; locy(1) = (7 * RcvLoc1 + 1 * RcvLoc2) / 8; locy(2) = (6 * RcvLoc1 + 2 * RcvLoc2) / 8;
					locy(3) = (5 * RcvLoc1 + 3 * RcvLoc2) / 8; locy(4) = (4 * RcvLoc1 + 4 * RcvLoc2) / 8;
					locy(5) = (3 * RcvLoc1 + 5 * RcvLoc2) / 8; locy(6) = (2 * RcvLoc1 + 6 * RcvLoc2) / 8;
					locy(7) = (1 * RcvLoc1 + 7 * RcvLoc2) / 8;  locy(8) = RcvLoc2;

				}//end of if only recieving one loc data;
				else if (argc - count == 9) {
					double indata[9];
					double BufferData;

					for (int i = 0; i<9; i++) {
						if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
							opserr << "WARNING eleLoad - invalid data " << argv[count] << " for -beamThermal 3D\n";
							return TCL_ERROR;
						}
						indata[i] = BufferData;
						count++;
					}
					locy = Vector(indata, 9);
					//temp1,loc1,temp2,loc2...temp9,loc9
				}//end of if only recieving 9 loc data;

				theLoad = new NodalThermalAction(nodeLoadTag, nodeId, locy, theSeries, thecrds);
				delete thecrds;
			}
			//end of recieving 9 temp data in external file;
			else {
				opserr << "WARNING NodalThermalAction - invalid dataLen\n";
			}
			//end of definition for different data input length(9 or15)
		}
		//end for detecting source 
		else {
			if (argc - count == 4) {
				double t1, t2, locY1, locY2;
				if (Tcl_GetDouble(interp, argv[count], &t1) != TCL_OK) {
					opserr << "WARNING eleLoad - invalid T1 " << argv[count] << " for NodalThermalAction\n";
					return TCL_ERROR;
				}
				if (Tcl_GetDouble(interp, argv[count + 1], &locY1) != TCL_OK) {
					opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1] << " for NodalThermalAction\n";
					return TCL_ERROR;
				}
				if (Tcl_GetDouble(interp, argv[count + 2], &t2) != TCL_OK) {
					opserr << "WARNING eleLoad - invalid T1 " << argv[count] << " for NodalThermalAction\n";
					return TCL_ERROR;
				}
				if (Tcl_GetDouble(interp, argv[count + 3], &locY2) != TCL_OK) {
					opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1] << " for NodalThermalAction\n";
					return TCL_ERROR;
				}

				theLoad = new NodalThermalAction(nodeLoadTag, nodeId, t1, locY1, t2, locY2, thecrds);
			}
			//for defining a uniform gradient thermal action
		}
		//end for source or no source 
		if (theLoad == 0) {
			opserr << "WARNING NodalLoad - out of memory creating load " << argv[1];
			return TCL_ERROR;
		}
		// get the current pattern tag if no tag given in i/p
		if (userSpecifiedPattern == false) {
			if (theTclLoadPattern == 0) {
				opserr << "WARNING no current load pattern - NodalThermalAction " << nodeId;
				return TCL_ERROR;
			}
			else
				loadPatternTag = theTclLoadPattern->getTag();

		}

	}
	//end of adding NodalThermalAction -------------end---------Liming,[SIF] 2017

	//start of else block, Liming [SIF]
	else{
	

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
	for (int i = 0; i < ndf; i++) {
		double theForce;
		if (Tcl_GetDouble(interp, argv[2 + i], &theForce) != TCL_OK) {
			opserr << "WARNING invalid force " << i + 1 << " - load " << nodeId;
			opserr << " " << ndf << " forces\n";
			return TCL_ERROR;
		}
		else
			forces(i) = theForce;
	}

	// allow some additional options at end of command
	int endMarker = 2 + ndf;
	while (endMarker != argc) {
		if (strcmp(argv[endMarker], "-const") == 0) {
			// allow user to specify const load
			isLoadConst = true;
		}
		else if (strcmp(argv[endMarker], "-pattern") == 0) {
			// allow user to specify load pattern other than current
			endMarker++;
			userSpecifiedPattern = true;
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
	if (userSpecifiedPattern == false) {
		if (theTclLoadPattern == 0) {
			opserr << "WARNING no current load pattern - load " << nodeId;
			opserr << " " << ndf << " forces\n";
			return TCL_ERROR;
		}
		else
			loadPatternTag = theTclLoadPattern->getTag();
	}

	// create the load
	theLoad = new NodalLoad(nodeLoadTag, nodeId, forces, isLoadConst);
	if (theLoad == 0) {
		opserr << "WARNING ran out of memory for load  - load " << nodeId;
		opserr << " " << ndf << " forces\n";
		return TCL_ERROR;
	}

} // end of Liming change for nodal thermal action , putting the above block into else{ }

  // add the load to the domain
  if (theTclDomain->addNodalLoad(theLoad, loadPatternTag) == false) {
    opserr << "WARNING TclModelBuilder - could not add load to domain\n";
    printCommand(argc, argv);
    delete theLoad;
    return TCL_ERROR;
  }
  nodeLoadTag++;

  // if get here we have sucessfully created the load and added it to the domain
  return TCL_OK;
}




int
TclCommand_addElementalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   
			 TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING current builder has been destroyed - eleLoad\n";    
    return TCL_ERROR;
  }

  if (theTclLoadPattern == 0) {
    opserr << "WARNING no active load pattern - eleLoad\n";    
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
      double wta;
      double waa = 0.0;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &wta) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wt for beamUniform \n";
	return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &waa) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wa for beamUniform \n";
	return TCL_ERROR;
      }
      double aL = 0.0;
      double bL = 1.0;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &aL) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid aOverL for beamUniform \n";
	return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &bL) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid bOverL for beamUniform \n";
	return TCL_ERROR;
      }
      double wab = waa;
      double wtb = wta;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wtb) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wt for beamUniform \n";
	return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wab) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wa for beamUniform \n";
	return TCL_ERROR;
      }      
      for (int i=0; i<theEleTags.Size(); i++) {
	if (aL > 0.0 || bL < 1.0 || wta != wtb || waa != wab)
	  theLoad = new Beam2dPartialUniformLoad(eleLoadTag, wta, wtb, waa, wab, aL, bL, theEleTags(i));
	else 
	  theLoad = new Beam2dUniformLoad(eleLoadTag, wta, waa, theEleTags(i));    

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
      }
      
      return 0;
    }
    else if (ndm == 3) {
      double wy, wz;
      double wx = 0.0;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &wy) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wy for beamUniform \n";
	return TCL_ERROR;
      }
      count++;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &wz) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wz for beamUniform \n";
	return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wx) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid wx for beamUniform \n";
	return TCL_ERROR;
      }
      double aL = 0.0;
      double bL = 1.0;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &aL) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid aOverL for beamUniform \n";
	return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &bL) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid bOverL for beamUniform \n";
	return TCL_ERROR;
      }

      for (int i=0; i<theEleTags.Size(); i++) {
	if (aL > 0.0 || bL < 1.0)
	  theLoad = new Beam3dPartialUniformLoad(eleLoadTag, wy, wz, wx, aL, bL, theEleTags(i));
	else 
	  theLoad = new Beam3dUniformLoad(eleLoadTag, wy, wz, wx, theEleTags(i));    	

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
      }
      
      return 0;

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
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &P) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid P for beamPoint\n";		
	return TCL_ERROR;
      } 
      if (count+1 >= argc || Tcl_GetDouble(interp, argv[count+1], &x) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid xDivL for beamPoint\n";	
	return TCL_ERROR;
      } 
      if (count+2 < argc && Tcl_GetDouble(interp, argv[count+2], &N) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid N for beamPoint\n";		
	return TCL_ERROR;
      } 

      if (x < 0.0 || x > 1.0) {
	opserr << "WARNING eleLoad - invalid xDivL of " << x;
	opserr << " for beamPoint (valid range [0.0, 1.0]\n";
	return TCL_ERROR;
      }


      for (int i=0; i<theEleTags.Size(); i++) {
	theLoad = new Beam2dPointLoad(eleLoadTag, P, x, theEleTags(i), N);    

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
      }
      
      return 0;

    }
    else if (ndm == 3) {
      double Py, Pz, x;
      double N = 0.0;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &Py) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid Py for beamPoint\n";		
	return TCL_ERROR;
      } 
      if (count+1 >= argc || Tcl_GetDouble(interp, argv[count+1], &Pz) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid Pz  for beamPoint\n";		
	return TCL_ERROR;
      } 
      if (count+2 >= argc || Tcl_GetDouble(interp, argv[count+2], &x) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid xDivL for beamPoint\n";	
	return TCL_ERROR;
      } 
      if (count+3 < argc && Tcl_GetDouble(interp, argv[count+3], &N) != TCL_OK) {
	opserr << "WARNING eleLoad - invalid N for beamPoint\n";		
	return TCL_ERROR;
      } 

      if (x < 0.0 || x > 1.0) {
	opserr << "WARNING eleLoad - invalid xDivL of " << x;
	opserr << " for beamPoint (valid range [0.0, 1.0]\n";
	return TCL_ERROR;
      }

      for (int i=0; i<theEleTags.Size(); i++) {
	theLoad = new Beam3dPointLoad(eleLoadTag, Py, Pz, x, theEleTags(i), N);    

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
      }
      return 0;
    }
    else {
      opserr << "WARNING eleLoad beamPoint type currently only valid only for ndm=2 or 3\n";
      return TCL_ERROR;
    }  
  }
  // Added Joey Yang UC Davis
  else if (strcmp(argv[count],"-BrickW") == 0) {

      for (int i=0; i<theEleTags.Size(); i++) {
	theLoad = new BrickSelfWeight(eleLoadTag, theEleTags(i));

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
      }
      return 0;
  }
  // Added: C.McGann, U.Washington
  else if ((strcmp(argv[count],"-surfaceLoad") == 0) || (strcmp(argv[count],"-SurfaceLoad") == 0)) {
	  count++;
  	  for (int i=0; i<theEleTags.Size(); i++) {
		  theLoad = new SurfaceLoader(eleLoadTag, theEleTags(i));

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
      }
      return 0;
  }
  // Added: C.McGann, U.Washington
  else if ((strcmp(argv[count],"-selfWeight") == 0) || (strcmp(argv[count],"-SelfWeight") == 0)) {
    count++;
    
    double xf = 0.0, yf = 0.0, zf = 0.0;
    if (Tcl_GetDouble(interp, argv[count], &xf) != TCL_OK) {
	  opserr << "WARNING eleLoad - invalid xFactor " << argv[count] << " for -selfWeight\n";
	  return TCL_ERROR;
	}
    if (Tcl_GetDouble(interp, argv[count+1], &yf) != TCL_OK) {
	  opserr << "WARNING eleLoad - invalid yFactor " << argv[count+1] << " for -selfWeight\n";
	  return TCL_ERROR;
	}
    if (count+2 < argc) { // adding to stop seg faults
      if (Tcl_GetDouble(interp, argv[count+2], &zf) != TCL_OK) {
	    opserr << "WARNING eleLoad - invalid zFactor " << argv[count+2] << " for -selfWeight\n";
	    return TCL_ERROR;
      }
    }
    
    for (int i=0; i<theEleTags.Size(); i++) {
      theLoad = new SelfWeight(eleLoadTag, xf, yf, zf, theEleTags(i));
      
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
    }
    return 0;
  }

 ///---------------------- Adding identifier for ThermalAction : [END] [SIF]----------------------------------------------//
  //-----------------Adding tcl command for shell thermal action, 2013..[Begin]---------------------
  else if (strcmp(argv[count], "-shellThermal") == 0) {
	  count++;
	  //so far three kinds of temperature distribution
	  //(1) 9 temperature points, i.e. 8 layers
	  //(2) 5 temperature points, i.e. 4 layers
	  //(3) 2 temperature points, i.e. 1 layers: linear or uniform

	  double t1, locY1, t2, locY2, t3, locY3, t4, locY4, t5, locY5,
		  t6, locY6, t7, locY7, t8, locY8, t9, locY9;
	  // 9 temperature points are given,i.e. 8 layers are defined; Also the 9 corresponding vertical coordinate is given.
	  // the temperature at each fiber is obtained by interpolating of temperatures at the nearby temperature points.
	  //Start to add source file
	  if (strcmp(argv[count], "-source") == 0) {
		  if (strcmp(argv[count + 1], "-node") != 0) {
			  count++;

			  const char *pwd = getInterpPWD(interp);
			  simulationInfo.addInputFile(argv[count], pwd);
			  TimeSeries* theSeries = new PathTimeSeriesThermal(eleLoadTag, argv[count]);

			  count++;

			  double RcvLoc1, RcvLoc2;
			  if (argc - count == 2) {

				  if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid single loc  " << argv[count] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }
				  if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid single loc  " << argv[count + 1] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }

			  }
			  else {
				  opserr << "WARNING eleLoad - invalid input for -shellThermal\n";
			  }

			  for (int i = 0; i<theEleTags.Size(); i++) {
				  theLoad = new ShellThermalAction(eleLoadTag, RcvLoc1, RcvLoc2,
					  theSeries, theEleTags(i));
				  if (theLoad == 0) {
					  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
			  }
		  }
		  //if not using nodal thermal action input
		  else {
			  for (int i = 0; i<theEleTags.Size(); i++) {
				  theLoad = new ShellThermalAction(eleLoadTag, theEleTags(i));
				  if (theLoad == 0) {
					  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
			  }//end of for loop
			  return 0;
		  }//end of <if(strcmp(argv[count+1],"-node") = 0)>
	  }
	  //end of the interface for importing temperature data from external file
	  else
	  {
		  if (argc - count == 18) {
			  double indata[18];
			  double BufferData;

			  for (int i = 0; i<18; i++) {
				  if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid data " << argv[count] << " for -beamThermal 3D\n";
					  return TCL_ERROR;
				  }
				  indata[i] = BufferData;
				  count++;
			  }

			  //temp1,loc1,temp2,loc2...temp9,loc9

			  for (int i = 0; i<theEleTags.Size(); i++) {
				  theLoad = new ShellThermalAction(eleLoadTag,
					  indata[0], indata[1], indata[2], indata[3],
					  indata[4], indata[5], indata[6], indata[7],
					  indata[8], indata[9], indata[10], indata[11],
					  indata[12], indata[13], indata[14], indata[15],
					  indata[16], indata[17], theEleTags(i));


				  if (theLoad == 0) {
					  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
			  }
			  return 0;
		  }
		  // 5 temperatures are given, i.e. 4 layers are defined.
		  else if (argc - count == 10) {
			  double indata[10];
			  double BufferData;

			  for (int i = 0; i<10; i++) {
				  if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid data " << argv[count] << " for -beamThermal 3D\n";
					  return TCL_ERROR;
				  }
				  indata[i] = BufferData;
				  count++;
			  }

			  //temp1,loc1,temp2,loc2...temp5,loc5

			  for (int i = 0; i<theEleTags.Size(); i++) {
				  theLoad = new ShellThermalAction(eleLoadTag,
					  indata[0], indata[1], indata[2], indata[3],
					  indata[4], indata[5], indata[6], indata[7],
					  indata[8], indata[9], theEleTags(i));


				  if (theLoad == 0) {
					  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
			  }
			  return 0;
		  }
		  // two temperature is given, 
		  //if the two temperatures are equal,i.e. uniform Temperature change in element
		  //if the two temperatures are different,i.e. linear Temperature change in element
		  else if (argc - count == 4) {
			  if (Tcl_GetDouble(interp, argv[count], &t1) != TCL_OK) {
				  opserr << "WARNING eleLoad - invalid T1 " << argv[count] << " for -shellThermal\n";
				  return TCL_ERROR;
			  }

			  if (Tcl_GetDouble(interp, argv[count + 1], &locY1) != TCL_OK) {
				  opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1] << " for -shellThermal\n";
				  return TCL_ERROR;
			  }
			  if (Tcl_GetDouble(interp, argv[count + 2], &t2) != TCL_OK) {
				  opserr << "WARNING eleLoad - invalid T2 " << argv[count] << " for -shellThermal\n";
				  return TCL_ERROR;
			  }

			  if (Tcl_GetDouble(interp, argv[count + 3], &locY2) != TCL_OK) {
				  opserr << "WARNING eleLoad - invalid LocY2 " << argv[count + 1] << " for -shellThermal\n";
				  return TCL_ERROR;
			  }

			  for (int i = 0; i<theEleTags.Size(); i++) {
				  theLoad = new ShellThermalAction(eleLoadTag,
					  t1, locY1, t2, locY2, theEleTags(i));


				  if (theLoad == 0) {
					  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
			  }
			  return 0;
		  }
		  //finish the temperature arguments
		  else {
			  opserr << "WARNING eleLoad -shellThermalLoad invalid number of temperature aguments,/n looking for 0, 1, 2 or 4 arguments.\n";
		  }
	  }
	  //end of if(recieved argument is not "source" or direct temperature input)//Liming,2014
  }
  //-----------------Adding tcl command for shell thermal action, 2013..[End]-----------------------

  else if (strcmp(argv[count], "-ThermalWrapper") == 0 || strcmp(argv[count], "-thermalWrapper") == 0) {

	  count++;
	  Vector loc = 0;;
	  ID NodalThermal = 0;
	  int numNodal = 0;

	  if (strcmp(argv[count], "-nodeLoc") == 0 || strcmp(argv[count], "nodeLoc") == 0)
	  {
		  count++;
		  numNodal = (argc - count) / 2;
		  loc.resize(numNodal);
		  NodalThermal.resize(numNodal);

		  for (int i = 0; i<numNodal; i++) {

			  double Dblloc;
			  int NodalTtag;

			  if (Tcl_GetDouble(interp, argv[count + i * 2 + 1], &Dblloc) != TCL_OK) {
				  opserr << "WARNING NodalLoad - invalid loc  " << argv[count] << " for NodalThermalAction\n";
				  return TCL_ERROR;
			  }

			  if (Tcl_GetInt(interp, argv[count + 2 * i], &NodalTtag) != TCL_OK) {
				  opserr << "WARNING invalid nodeId: " << argv[1];
				  return TCL_ERROR;
			  }

			  loc(i) = Dblloc; NodalThermal(i) = NodalTtag;
		  }
		  //end of for loop over numNodal
	  }
	  //end of nodelLoc
	  else if (strcmp(argv[count], "-node") == 0 || strcmp(argv[count], "node") == 0)
	  {
		  count++;
		  numNodal = argc - count;
		  NodalThermal.resize(numNodal);

		  for (int i = 0; i<numNodal; i++) {

			  int NodalTtag;

			  if (Tcl_GetInt(interp, argv[count + i], &NodalTtag) != TCL_OK) {
				  opserr << "WARNING invalid nodeId: " << argv[1];
				  return TCL_ERROR;
			  }

			  NodalThermal(i) = NodalTtag;
		  }
		  //end of for loop over numNodal
	  }
	  //end of node tag

	  //Obtain Pointers to NodalThermalAction;
	  Node* theNode = 0;
	  NodalThermalAction** theNodalThermals = 0;
	  theNodalThermals = new NodalThermalAction*[numNodal];
	  
	  for (int i = 0; i<numNodal; i++) {
	    theNode = theTclDomain->getNode(NodalThermal(i));
	    theNodalThermals[i] = theNode->getNodalThermalActionPtr();
	    if (theNodalThermals[i] == 0) {
	      opserr << "WARNING:: An empty nodalThermalAction detected for ThermalActionWrapper" << endln;
	      return TCL_ERROR;
	    }
	  }
	  
	  for (int i = 0; i<theEleTags.Size(); i++) {
	    if (numNodal == 2) {
	      theLoad = new ThermalActionWrapper(eleLoadTag, theEleTags(i), theNodalThermals[0], theNodalThermals[1]);
	    }
	    else if (numNodal == 3) {
	      theLoad = new ThermalActionWrapper(eleLoadTag, theEleTags(i), theNodalThermals[0], theNodalThermals[1], theNodalThermals[2]);
	    }
	    else if (numNodal == 4) {
	      theLoad = new ThermalActionWrapper(eleLoadTag, theEleTags(i), theNodalThermals[0], theNodalThermals[1], theNodalThermals[2], theNodalThermals[3]);
	    }
	    else if (numNodal == 5) {
	      theLoad = new ThermalActionWrapper(eleLoadTag, theEleTags(i), theNodalThermals[0], theNodalThermals[1],
						 theNodalThermals[2], theNodalThermals[3], theNodalThermals[4]);
	    }
	    else if (numNodal == 6) {
	      theLoad = new ThermalActionWrapper(eleLoadTag, theEleTags(i), theNodalThermals[0], theNodalThermals[1],
						 theNodalThermals[2], theNodalThermals[3], theNodalThermals[4], theNodalThermals[5]);
	    }
	    
	    
	    if (theLoad == 0) {
	      opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
			  return TCL_ERROR;
		  }

		  if (loc != 0)
			  ((ThermalActionWrapper*)theLoad)->setRatios(loc);

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
	  }//end of for loop
	  return 0;
  }
  //------------------------end  of using ThermalActionWrapper--------------------------
  //-----------------Adding tcl command for beam thermal action(2D&3D), 2013..[Begin]---------------

  else if (strcmp(argv[count], "-beamThermal") == 0) {
	  count++;
	  //For two dimensional model
	  if (ndm == 2) {
		  // The thermal action can be defined with external data file, which is identified with "-source"
		  // The external file itself can either be elemental data or nodal data, the latter is identified with "-node"
		  if (strcmp(argv[count], "-source") == 0) {

			  if (strcmp(argv[count + 1], "-node") == 0) {
				  for (int i = 0; i<theEleTags.Size(); i++) {
					  theLoad = new Beam2dThermalAction(eleLoadTag, theEleTags(i));
					  if (theLoad == 0) {
						  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
				  }//end of for loop
				  return 0;

			  }
			  //end of <if(strcmp(argv[count+1],"-node") != 0)>
			  else {
				  count++;
				  const char *pwd = getInterpPWD(interp);
				  simulationInfo.addInputFile(argv[count], pwd);
				  TimeSeries* theSeries = new PathTimeSeriesThermal(eleLoadTag, argv[count]);

				  count++;
				  Vector locs(9);
				  //---------------for recieving 2 arguments
				  if (argc - count == 2) {
					  double RcvLoc1, RcvLoc2;
					  if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
						  opserr << "WARNING eleLoad - invalid single loc  " << argv[count] << " for -beamThermal\n";
						  return TCL_ERROR;
					  }
					  if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
						  opserr << "WARNING eleLoad - invalid single loc  " << argv[count + 1] << " for -beamThermal\n";
						  return TCL_ERROR;
					  }

					  locs(0) = RcvLoc1; locs(8) = RcvLoc2;
					  for (int i = 1; i<8; i++) {
						  locs(i) = locs(0) - i*(locs(0) - locs(8)) / 8;
					  }
				  }
				  //----------------for recieving 9 arguments
				  else if (argc - count == 9) {

					  int ArgStart = count;
					  int ArgEnd = argc;
					  double data;

					  if (ArgStart != ArgEnd) {
						  for (int i = ArgStart; i<ArgEnd; i++) {
							  Tcl_GetDouble(interp, argv[i], &data);
							  locs(i - ArgStart) = data;
						  }
					  }

				  }
				  //end of recieving 9 arguments
				  else {
					  opserr << "WARNING eleLoad - invalid input for -beamThermal\n";
				  }
#ifdef _DEBUG 
				  opserr << "TclModelBuilder:: locs" << locs << endln;
#endif
				  for (int i = 0; i<theEleTags.Size(); i++) {
					  theLoad = new Beam2dThermalAction(eleLoadTag, locs,
						  theSeries, theEleTags(i));
					  if (theLoad == 0) {
						  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
				  }//end of for loop
				  return 0;
			  }//end of <if(strcmp(argv[count+1],"-node") = 0)>
			   //--------------------------end for beam2DThermalAction with time series ----------------------------------------	
		  }
		  else
		  {
			  //(1) 9 temperature points, i.e. 8 layers
			  //(2) 5 temperature points, i.e. 4 layers
			  //(3) 2 temperature points, i.e. 1 layers: linear or uniform
			  //double t1, locY1, t2, locY2, t3, locY3, t4, locY4, t5, locY5, t6, locY6, t7, locY7, t8, locY8, t9, locY9;
			  double Temp[9]; double Loc[9];
			  // 9 temperature points are given,i.e. 8 layers are defined; Also the 9 corresponding vertical coordinate is given.
			  // the temperature at each fiber is obtained by interpolating of temperatures at the nearby temperature points.
			  if (argc - count == 18) {
				  double indata[18];
				  double BufferData;

				  for (int i = 0; i<18; i++) {
					  if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
						  opserr << "WARNING eleLoad - invalid data " << argv[count] << " for -beamThermal 3D\n";
						  return TCL_ERROR;
					  }
					  indata[i] = BufferData;
					  count++;
				  }

				  for (int i = 0; i<9; i++) {
					  Temp[i] = indata[2 * i];
					  Loc[i]  = indata[2 * i + 1];
				  }

			  }

			  // 5 temperatures are given, i.e. 4 layers are defined.
			  else if (argc - count == 10) {

				  double indata[10];
				  double BufferData;

				  for (int i = 0; i<10; i++) {
					  if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
						  opserr << "WARNING eleLoad - invalid data " << argv[count] << " for -beamThermal 3D\n";
						  return TCL_ERROR;
					  }
					  indata[i] = BufferData;
					  count++;
				  }

				  Temp[0] = indata[0]; Temp[2] = indata[2]; Temp[4] = indata[4]; Temp[6] = indata[6]; Temp[8] = indata[8];
				  Loc[0]  = indata[1]; Loc[2]  = indata[3]; Loc[4]  = indata[5]; Loc[6]  = indata[7]; Loc[8]  = indata[9];

				  for (int i = 1; i<5; i++) {
					  Temp[2 * i - 1] = (Temp[2 * i - 2] + Temp[2 * i]) / 2;
					  Loc[2 * i - 1]  = (Loc[2 * i - 2] + Loc[2 * i]) / 2;
				  }
			  }
			  //End for 5 inputs
			  // two temperature is given, 
			  //if the two temperatures are equal,i.e. uniform Temperature change in element
			  //if the two temperatures are different,i.e. linear Temperature change in element
			  else if (argc - count == 4) {
				  double indata[4];
				  double BufferData;

				  for (int i = 0; i<4; i++) {
					  if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
						  opserr << "WARNING eleLoad - invalid data " << argv[count] << " for -beamThermal 3D\n";
						  return TCL_ERROR;
					  }
					  indata[i] = BufferData;
					  count++;
				  }

				  Temp[0] = indata[0]; Temp[8] = indata[2];
				  Loc[0]  = indata[1]; Loc[8]  = indata[3];
				  for (int i = 1; i<8; i++) {
					  Temp[i] = Temp[0] - i*(Temp[0] - Temp[8]) / 8;
					  Loc[i] = Loc[0] - i*(Loc[0] - Loc[8]) / 8;
				  }

			  }
			  //end for 2 inputs

			  for (int i = 0; i<theEleTags.Size(); i++) {
				  theLoad = new Beam2dThermalAction(eleLoadTag,
					  Temp[0], Loc[0], Temp[1], Loc[1],
					  Temp[2], Loc[2], Temp[3], Loc[3],
					  Temp[4], Loc[4], Temp[5], Loc[5],
					  Temp[6], Loc[6], Temp[7], Loc[7],
                      Temp[8], Loc[8], theEleTags(i));

				  if (theLoad == 0) {
					  opserr << "WARNING eleLoad - out of memory creating load of type beamThermal";
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
			  }

			  return 0;
			  //-End of Adding BeamThermalAction defined with direct input.
		  }
		  //end for no sourced pattern
	  }
	  // End of for the if (ndm==2)
	  else if (ndm == 3)
	  {
		  //so far three kinds of temperature distribution
		  double t1, locY1, t2, locY2, t3, locY3, t4, locY4, t5, locY5,
			  t6, t7, locZ1, t8, t9, locZ2, t10, t11, locZ3, t12, t13, locZ4, t14, t15, locZ5;

		  // the temperature at each fiber is obtained by interpolating of temperatures at the nearby temperature points.
		  if (strcmp(argv[count], "-source") == 0) {
			  count++;
			  if (strcmp(argv[count], "-node") == 0) {
				  for (int i = 0; i<theEleTags.Size(); i++) {
					  theLoad = new Beam3dThermalAction(eleLoadTag, theEleTags(i));
					  if (theLoad == 0) {
						  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
				  }	//end of loop tf all elements defined
				  return 0;
			  }//end for defing thermal action with nodal input
			  else {
				  const char *pwd = getInterpPWD(interp);
				  simulationInfo.addInputFile(argv[count], pwd);
				  count++;
				  bool using2Ddata = false;

				  double RcvLoc1, RcvLoc2, RcvLoc3, RcvLoc4;
				  TimeSeries* theSeries;

				  if (argc - count == 4) {
					  theSeries = new PathTimeSeriesThermal(eleLoadTag, argv[count - 1], 15);
					  using2Ddata = false;

					  if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
						  opserr << "WARNING eleLoad - invalid single loc  " << argv[count] << " for -beamThermal\n";
						  return TCL_ERROR;
					  }
					  if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
						  opserr << "WARNING eleLoad - invalid single loc  " << argv[count + 1] << " for -beamThermal\n";
						  return TCL_ERROR;
					  }
					  if (Tcl_GetDouble(interp, argv[count + 2], &RcvLoc3) != TCL_OK) {
						  opserr << "WARNING eleLoad - invalid single loc  " << argv[count + 2] << " for -beamThermal\n";
						  return TCL_ERROR;
					  }
					  if (Tcl_GetDouble(interp, argv[count + 3], &RcvLoc4) != TCL_OK) {
						  opserr << "WARNING eleLoad - invalid single loc  " << argv[count + 3] << " for -beamThermal\n";
						  return TCL_ERROR;
					  }

					  //end for recieving input
					  for (int i = 0; i<theEleTags.Size(); i++) {

						  theLoad = new Beam3dThermalAction(eleLoadTag, RcvLoc1, RcvLoc2, RcvLoc3, RcvLoc4,
							  theSeries, theEleTags(i));

						  if (theLoad == 0) {
							  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
					  }
					  return 0;

				  }
				  //end of defining 15 data points with external file
				  else if (argc - count == 2 || argc - count == 9) {
					  // for receiving data which has the similiar structure as 2D beam section
					  Vector locs(9);
					  using2Ddata = true;
					  TimeSeries* theSeries = new PathTimeSeriesThermal(eleLoadTag, argv[count - 1], 9);
					  if (argc - count == 2) {

						  double RcvLoc1, RcvLoc2;
						  if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
							  opserr << "WARNING eleLoad - invalid single loc  " << argv[count] << " for -beamThermal\n";
							  return TCL_ERROR;
						  }
						  if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
							  opserr << "WARNING eleLoad - invalid single loc  " << argv[count + 1] << " for -beamThermal\n";
							  return TCL_ERROR;
						  }

						  locs(0) = RcvLoc1; locs(8) = RcvLoc2;
						  for (int i = 1; i<8; i++) {
							  locs(i) = locs(0) - i*(locs(0) - locs(8)) / 8;
						  }
					  }
					  //end of receiving 2 data points
					  else {
						  double indata[9];
						  double BufferData;
						  for (int i = 0; i<9; i++) {
							  if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
								  opserr << "WARNING eleLoad - invalid data " << argv[count] << " for -beamThermal 3D\n";
								  return TCL_ERROR;
							  }
							  indata[i] = BufferData;
							  count++;
						  }

						  for (int i = 0; i<9; i++) {
							  locs(i) = indata[i];
						  }

					  }
					  //end of receiving 9data points

					  for (int i = 0; i<theEleTags.Size(); i++) {
						  theLoad = new Beam3dThermalAction(eleLoadTag, locs,
							  theSeries, theEleTags(i));
						  if (theLoad == 0) {
							  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
					  }//end of for loop
					  return 0;
				  }

				  else {
					  opserr << "WARNING eleLoad - invalid input for -beamThermal\n";
				  }

			  }//end for source beam element temperature data  -source -filePath		

		  }
		  //-------------------------end for importing temp data from external files --source -------------------------
		  else {

			  //double t1, locY1, t2, locY2, t3, locY3, t4, locY4, t5, locY5, 
			  //t6, t7, locZ1, t8, t9, locZ2, t10,t11, locZ3, t12, t13, locZ4, t14,t15, locZ5;

			  if (argc - count == 25) {
				  double indata[25];
				  double BufferData;

				  for (int i = 0; i<25; i++) {
					  if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
						  opserr << "WARNING eleLoad - invalid data " << argv[count] << " for -beamThermal 3D\n";
						  return TCL_ERROR;
					  }
					  indata[i] = BufferData;
					  count++;
				  }

				  for (int i = 0; i<theEleTags.Size(); i++) {
					  theLoad = new Beam3dThermalAction(eleLoadTag,
						  indata[0], indata[1], indata[2], indata[3], indata[4], indata[5], indata[6], indata[7],
						  indata[8], indata[9], indata[10], indata[11], indata[12], indata[13], indata[14], indata[15],
						  indata[16], indata[17], indata[18], indata[19], indata[20], indata[21], indata[22], indata[23],
						  indata[24], theEleTags(i));
					  if (theLoad == 0) {
						  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
				  }
				  return 0;
			  }
			  //end of  if (argc-count == 25){
			  else if (argc - count == 4) {

				  if (Tcl_GetDouble(interp, argv[count], &t1) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid T1 " << argv[count] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }
				  if (Tcl_GetDouble(interp, argv[count + 1], &locY1) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }
				  if (Tcl_GetDouble(interp, argv[count + 2], &t5) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid T1 " << argv[count] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }
				  if (Tcl_GetDouble(interp, argv[count + 3], &locY5) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }

				  locY2 = locY1 + (locY5 - locY1) / 4;
				  locY3 = locY1 + (locY5 - locY1) / 2;
				  locY4 = locY1 + 3 * (locY5 - locY1) / 4;
				  t2 = t1 + (t5 - t1) / 4;
				  t3 = t1 + (t5 - t1) / 2;
				  t4 = t1 + 3 * (t5 - t1) / 4;
				  locZ1 = locZ2 = locZ3 = locZ4 = locZ5 = 0;
				  t6 = t7 = t8 = t9 = t10 = 0;
				  t11 = t12 = t13 = t14 = t15 = 0;

				  for (int i = 0; i<theEleTags.Size(); i++) {
					  theLoad = new Beam3dThermalAction(eleLoadTag,
						  t1, locY1, t2, locY2, t3, locY3, t4, locY4,
						  t5, locY5, t6, t7, locZ1, t8, t9, locZ2, t10, t11, locZ3,
						  t12, t13, locZ4, t14, t15, locZ5, theEleTags(i));
					  if (theLoad == 0) {
						  opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count];
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
				  }
				  return 0;
			  }
			  //end of  if (argc-count == 4){
			  else {
				  opserr << "WARNING eleLoad Beam3dThermalAction: invalid number of temperature aguments,/n looking for arguments for Temperatures and cordinates.\n";
			  }
		  } //end for not sourced pattern
	  }//else if ndm=3  
  }//else if '-beamThermal'

   //--Adding identifier for ThermalAction:[END] [SIF]--//

  
  
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
	
	for (int i=0; i<theEleTags.Size(); i++) {
	  theLoad = new Beam2dTempLoad(eleLoadTag, temp1, temp2, temp3, temp4, theEleTags(i));
	  
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
	}
	
	return 0;

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

	for (int i=0; i<theEleTags.Size(); i++) {
	  theLoad = new Beam2dTempLoad(eleLoadTag, temp1, temp2, theEleTags(i));
	  
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
	}
      }
      // One twmp change give, uniform temp change in element
      else if (argc-count == 1) {
	if (Tcl_GetDouble(interp, argv[count],&temp1 ) != TCL_OK) {
	  opserr << "WARNING eleLoad - invalid Tbot " << argv[count+1] << " for -beamTemp\n";	
	  return TCL_ERROR;
	}
	theLoad=0;

	for (int i=0; i<theEleTags.Size(); i++) {
	  theLoad = new Beam2dTempLoad(eleLoadTag, temp1, theEleTags(i));

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
	}
	
	return 0;

      }

      else {
	opserr << "WARNING eleLoad -beamTempLoad invalid number of temperature aguments,/n looking for 0, 1, 2 or 4 arguments.\n";
      }
    } else {
      opserr << "WARNING eleLoad -beamTempLoad type currently only valid only for ndm=2\n";
      return TCL_ERROR;
    }  
  }

  // if get here we have sucessfully created the load and added it to the domain
  return TCL_OK;
}



int
TclCommand_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc, 
                        TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - load \n";    
    return TCL_ERROR;
  }

  //  int ndf = theTclBuilder->getNDF();
  int ndf = argc - 2;

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

  if (theTclDomain->setMass(mass, nodeId) != 0) {
    opserr << "WARNING failed to set mass at node " << nodeId << endln;
    return TCL_ERROR;
  }
    
  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclCommand_addHomogeneousBC(ClientData clientData, Tcl_Interp *interp, int argc,   
				 TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  //  int ndf = theTclBuilder->getNDF();
  int ndf = argc - 2;

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

  char buffer[80];
  strcpy(buffer,"");
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
	SP_Constraint *theSP = new SP_Constraint(nodeId, i, 0.0, true);
	if (theSP == 0) {
	  opserr << "WARNING ran out of memory for SP_Constraint ";
	  opserr << "fix " << nodeId << " " << ndf << " [0,1] conditions\n";
	  return TCL_ERROR;
	}

	// add it to the domain
	if (theTclDomain->addSP_Constraint(theSP) == false) {
	  opserr << "WARNING could not add SP_Constraint to domain using fix command - node may already be constrained\n";
	  sprintf(buffer, "%d ", 0);
	  delete theSP;
	} else {
	  sprintf(buffer, "%d ", theSP->getTag());
	  Tcl_AppendResult(interp, buffer, NULL);
	}
      }
    }
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}


int
TclCommand_addHomogeneousBC_X(ClientData clientData, Tcl_Interp *interp, 
				   int argc, TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  //  int ndf = theTclBuilder->getNDF();
  int ndf = argc - 2;
  if (strcmp(argv[argc-2],"-tol") == 0)
    ndf -= 2;

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

  theTclDomain->addSP_Constraint(0, xLoc, fixity, tol);

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}




int
TclCommand_addHomogeneousBC_Y(ClientData clientData, Tcl_Interp *interp, 
				   int argc, TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  //  int ndf = theTclBuilder->getNDF();
  int ndf = argc - 2;
  if (strcmp(argv[argc-2],"-tol") == 0)
    ndf -= 2;

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


  theTclDomain->addSP_Constraint(1, yLoc, fixity, tol);

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclCommand_addHomogeneousBC_Z(ClientData clientData, Tcl_Interp *interp, 
				   int argc, TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - elasticBeam \n";    
    return TCL_ERROR;
  }

  //int ndf = theTclBuilder->getNDF();
  int ndf = argc - 2;
  if (strcmp(argv[argc-2],"-tol") == 0)
    ndf -= 2;

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

  theTclDomain->addSP_Constraint(2, zLoc, fixity, tol);

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclCommand_addSP(ClientData clientData, Tcl_Interp *interp, int argc,   
		      TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - sp \n";    
    return TCL_ERROR;
  }

  //int ndf = theTclBuilder->getNDF();

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
  bool userSpecifiedPattern = false;
  int loadPatternTag = 0; // some pattern that will never be used!

  int endMarker = 4;
  while (endMarker != argc) {
    if (strcmp(argv[endMarker],"-const") == 0) {
      // allow user to specify const load
      isSpConst = true;
    } else if (strcmp(argv[endMarker],"-pattern") == 0) {
      // allow user to specify load pattern other than current
      endMarker++;
      userSpecifiedPattern = true;
      if (endMarker == argc || 
	  Tcl_GetInt(interp, argv[endMarker], &loadPatternTag) != TCL_OK) {

	opserr << "WARNING invalid patternTag - load " << nodeId << "\n";
	return TCL_ERROR;
      }
    }  
    endMarker++;
  }

  // if load pattern tag has not changed - get the pattern tag from current one
  if (userSpecifiedPattern == false) {
    if (theTclLoadPattern == 0) {
      opserr << "WARNING no current pattern - sp " << nodeId << " dofID value\n";
      return TCL_ERROR;
    } else	
      loadPatternTag = theTclLoadPattern->getTag();
  }
  
  LoadPattern *thePattern = theTclDomain->getLoadPattern(loadPatternTag);
  
  // create a homogeneous constraint
  SP_Constraint *theSP = new SP_Constraint(nodeId, dofId, value, isSpConst);

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
TclCommand_addImposedMotionSP(ClientData clientData, 
				   Tcl_Interp *interp, 
				   int argc,   
				   TCL_Char **argv)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed - sp \n";    
    return TCL_ERROR;
  }

  //  int ndf = theTclBuilder->getNDF();

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

  //
  // check valid node & dof
  //

  Node *theNode = theTclDomain->getNode(nodeId);
  if (theNode == 0) {
    opserr << "WARNING invalid node " << argv[2] << " node not found\n ";
    return -1;
  }
  int nDof = theNode->getNumberDOF();
  if (dofId < 0 || dofId >= nDof) {
    opserr << "WARNING invalid dofId: " << argv[2] << " dof specified cannot be <= 0 or greater than num dof at nod\n ";
    return -2;
  }


  MultiSupportPattern *thePattern = theTclMultiSupportPattern;
  int loadPatternTag = thePattern->getTag();
  
  // create a new ImposedMotionSP
  SP_Constraint *theSP;
  if (alt == true) {
    theSP = new ImposedMotionSP1(nodeId, dofId, loadPatternTag, gMotionID);
  }
  else {
    theSP = new ImposedMotionSP(nodeId, dofId, loadPatternTag, gMotionID);
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
TclCommand_addEqualDOF_MP (ClientData clientData, Tcl_Interp *interp,
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

        // Create the multi-point constraint
        MP_Constraint *theMP = new MP_Constraint (RnodeID, CnodeID, Ccr, rcDOF, rcDOF);
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

	char buffer[80];
	sprintf(buffer, "%d", theMP->getTag());
	Tcl_SetResult(interp, buffer, TCL_VOLATILE);

        return TCL_OK;
}


int
TclCommand_addEqualDOF_MP_Mixed(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char **argv)
{
        // Ensure the destructor has not been called
        if (theTclBuilder == 0) {
	  opserr << "WARNING builder has been destroyed - equalDOF \n";
	  return TCL_ERROR;
        }

        // Check number of arguments
        if (argc < 4) {
	  opserr << "WARNING bad command - want: equalDOFmixed RnodeID? CnodeID? numDOF? RDOF1? CDOF1? ... ...";
	  printCommand (argc, argv);
	  return TCL_ERROR;
        }

        // Read in the node IDs and the DOF
        int RnodeID, CnodeID, dofIDR, dofIDC, numDOF;

        if (Tcl_GetInt (interp, argv[1], &RnodeID) != TCL_OK) {
	  opserr << "WARNING invalid RnodeID: " << argv[1]
	       << " equalDOF RnodeID? CnodeID? numDOF? RDOF1? CDOF1? ...";
	  return TCL_ERROR;
        }
        if (Tcl_GetInt (interp, argv[2], &CnodeID) != TCL_OK) {
	  opserr << "WARNING invalid CnodeID: " << argv[2]
	       << " equalDOF RnodeID? CnodeID? numDOF? RDOF1? CDOF1? ...";
	  return TCL_ERROR;
        }

        if (Tcl_GetInt (interp, argv[3], &numDOF) != TCL_OK) {
	  opserr << "WARNING invalid numDOF: " << argv[2]
	       << " equalDOF RnodeID? CnodeID? numDOF? RDOF1? CDOF1? ...";
	  return TCL_ERROR;
        }

        // The number of DOF to be coupled
	//        int numDOF = argc - 3;

        // The constraint matrix ... U_c = C_cr * U_r
        Matrix Ccr (numDOF, numDOF);
        Ccr.Zero();

        // The vector containing the retained and constrained DOFs
        ID rDOF (numDOF);
        ID cDOF (numDOF);

        int i, j, k;
        // Read the degrees of freedom which are to be coupled
        for (i = 4, j = 5, k = 0; k < numDOF; i+=2, j+=2, k++) {
	  if (Tcl_GetInt (interp, argv[i], &dofIDR) != TCL_OK) {
	    opserr << "WARNING invalid dofID: " << argv[3]
		 << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
	    return TCL_ERROR;
	  }
	  if (Tcl_GetInt (interp, argv[j], &dofIDC) != TCL_OK) {
	    opserr << "WARNING invalid dofID: " << argv[3]
		 << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
	    return TCL_ERROR;
	  }

	  dofIDR -= 1; // Decrement for C++ indexing
	  dofIDC -= 1;
	  if (dofIDC < 0 || dofIDR < 0) {
	    opserr << "WARNING invalid dofID: " << argv[i]
		   << " must be >= 1";
	    return TCL_ERROR;
	  }
	  rDOF(k) = dofIDR;    
	  cDOF(k) = dofIDC;    
	  Ccr(k,k) = 1.0;
        }

        // Create the multi-point constraint
        MP_Constraint *theMP = new MP_Constraint (RnodeID, CnodeID, Ccr, cDOF, rDOF);
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

	char buffer[80];
	sprintf(buffer, "%d", theMP->getTag());
	Tcl_SetResult(interp, buffer, TCL_VOLATILE);

        return TCL_OK;
}


int 
TclCommand_RigidLink(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 4) {
      opserr << "WARNING rigidLink linkType? rNode? cNode?\n";
      return TCL_ERROR;
  }    

  int rNode, cNode;
  if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
      opserr << "WARNING rigidLink linkType? rNode? cNode? - could not read rNode \n";
      return TCL_ERROR;	        
  }
  if (Tcl_GetInt(interp, argv[3], &cNode) != TCL_OK) {
      opserr << "WARNING rigidLink linkType? rNode? cNode? - could not read CNode \n";
      return TCL_ERROR;	        
  }

  // construct a rigid rod or beam depending on 1st arg
  if ((strcmp(argv[1],"-bar") == 0) || (strcmp(argv[1],"bar") == 0)) {
    RigidRod theLink(*theTclDomain, rNode, cNode);
  } else if ((strcmp(argv[1],"-beam") == 0) || (strcmp(argv[1],"beam") == 0)) {
    RigidBeam theLink(*theTclDomain, rNode, cNode);
  } else {
      opserr << "WARNING rigidLink linkType? rNode? cNode? - unrecognised link type (-bar, -beam) \n";
      return TCL_ERROR;	        
  }

  return TCL_OK;
}

int 
TclCommand_RigidDiaphragm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 3) {
      opserr << "WARNING rigidLink perpDirn? rNode? <cNodes?>\n";
      return TCL_ERROR;
  }    

  int rNode, perpDirn;
  if (Tcl_GetInt(interp, argv[1], &perpDirn) != TCL_OK) {
      opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read perpDirn? \n";
      return TCL_ERROR;	        
  }

  if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
      opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read rNode \n";
      return TCL_ERROR;	        
  }
  
  // read in the constrained Nodes
  int numConstrainedNodes = argc - 3;
  ID constrainedNodes(numConstrainedNodes);
  for (int i=0; i<numConstrainedNodes; i++) {
      int cNode;
      if (Tcl_GetInt(interp, argv[3+i], &cNode) != TCL_OK) {
	  opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read a cNode\n";
	  return TCL_ERROR;	        
      }
      constrainedNodes(i) = cNode;
  }

  RigidDiaphragm theLink(*theTclDomain, rNode, constrainedNodes, perpDirn-1);
	

  return TCL_OK;
}




int
TclCommand_addMP(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv)
{
  opserr << "WARNING - TclCommand_addMP() not yet implemented\n";
  return TCL_OK;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Added by Scott J. Brandenberg, UC Davis, sjbrandenberg@ucdavis.edu
int
TclCommand_doPySimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
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
TclCommand_doTzSimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
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

// Added by Prishati Raychowdhury (UCSD)
int
TclModelBuilder_doShallowFoundationGen(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char **argv)
{
    	if(argc != 5){
		opserr << "WARNING ShallowFoundationGen FoundationID? ConnectingNode? InputDataFile? FoundationMatType?";
		opserr << "Must have 4 arguments." << endln;
	}
	
	ShallowFoundationGen *theShallowFoundationGen;
	theShallowFoundationGen = new ShallowFoundationGen;

	
      // Checking for error
        int FoundationID; int ConnectingNode; int FoundationMatType;

		if (Tcl_GetInt (interp, argv[1], &FoundationID) != TCL_OK) {
	  opserr << "WARNING invalid FoundationID: " << argv[1]
	       << ". ShallowFoundationGen FoundationID? ConnectingNode? InputDataFile? FoundationMatType? ";
	  return TCL_ERROR;
        }
        if (Tcl_GetInt (interp, argv[2], &ConnectingNode) != TCL_OK) {
	  opserr << "WARNING invalid ConnectingNode: " << argv[2]
	       << ". ShallowFoundationGen FoundationID? ConnectingNode? InputDataFile? FoundationMatType? ";
	  return TCL_ERROR;
        }
        if (Tcl_GetInt (interp, argv[4], &FoundationMatType) != TCL_OK) {
	  opserr << "WARNING invalid FoundationMatType: " << argv[4]
	       << ". ShallowFoundationGen FoundationID? ConnectingNode? InputDataFile? FoundationMatType? ";
	  return TCL_ERROR;
        }
     
	theShallowFoundationGen->GetShallowFoundation(argv[1], argv[2], argv[3], argv[4]);
	delete theShallowFoundationGen;

	return TCL_OK;
}
// End PRC

int
TclCommand_doBlock2D(ClientData clientData, Tcl_Interp *interp, int argc,   
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
TclCommand_doBlock3D(ClientData clientData, Tcl_Interp *interp, int argc,   
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
TclCommand_addRemoPatch(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv)
{
  return TclCommand_addPatch(clientData, interp, argc,argv,
				    theTclBuilder);
}

int
TclCommand_addRemoFiber(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv)
{
  return TclCommand_addFiber(clientData, interp, argc,argv,
				  theTclBuilder);
}

int				
TclModelBuilder_addRemoHFiber(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv)
{
  return TclCommand_addHFiber(clientData, interp, argc,argv,theTclBuilder);
				  
}

int
TclCommand_addRemoLayer(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv)
{
  return TclCommand_addReinfLayer(clientData, interp, argc,argv,
				       theTclBuilder);
}

					 
int
TclCommand_addRemoGeomTransf(ClientData clientData, Tcl_Interp *interp, int argc,   
			   TCL_Char **argv)
{
  return TclCommand_addGeomTransf(clientData, interp, argc,argv,
				       theTclDomain,
				       theTclBuilder);
}

#ifdef OO_HYSTERETIC
extern int
TclModelBuilderStiffnessDegradationCommand(ClientData clientData,
					   Tcl_Interp *interp,
					   int argc, TCL_Char **argv,
					   TclModelBuilder *theTclBuilder);

int
TclCommand_addStiffnessDegradation(ClientData clientData,
					Tcl_Interp *interp,
					int argc, TCL_Char **argv)
{
  return TclModelBuilderStiffnessDegradationCommand(clientData, interp, 
						    argc, argv, theTclBuilder);
}

extern int
TclModelBuilderUnloadingRuleCommand(ClientData clientData,
				    Tcl_Interp *interp,
				    int argc, TCL_Char **argv,
				    TclModelBuilder *theTclBuilder);

int
TclCommand_addUnloadingRule(ClientData clientData,
				 Tcl_Interp *interp,
				 int argc, TCL_Char **argv)
{
  return TclModelBuilderUnloadingRuleCommand(clientData, interp, 
					     argc, argv, theTclBuilder);
}

extern int
TclModelBuilderStrengthDegradationCommand(ClientData clientData,
					  Tcl_Interp *interp,
					  int argc, TCL_Char **argv,
					  TclModelBuilder *theTclBuilder);

int
TclCommand_addStrengthDegradation(ClientData clientData,
				       Tcl_Interp *interp,
				       int argc, TCL_Char **argv)
{
  return TclModelBuilderStrengthDegradationCommand(clientData, interp, 
						   argc, argv, theTclBuilder);
}
#endif

extern int
TclModelBuilderHystereticBackboneCommand(ClientData clientData,
					 Tcl_Interp *interp,
					 int argc, TCL_Char **argv);

int
TclCommand_addHystereticBackbone(ClientData clientData,
				      Tcl_Interp *interp,
				      int argc,	TCL_Char **argv)
{
  return TclModelBuilderHystereticBackboneCommand(clientData, interp, argc, argv);
}

/// added by ZHY
extern int 
TclModelBuilderUpdateMaterialStageCommand(ClientData clientData, 
					  Tcl_Interp *interp, 
					  int argc, 
					  TCL_Char **argv, 
					  TclModelBuilder *theTclBuilder,
					  Domain *theDomain);
int
TclCommand_UpdateMaterialStage(ClientData clientData, 
				    Tcl_Interp *interp,  
				    int argc, 
				    TCL_Char **argv)
{
  return TclModelBuilderUpdateMaterialStageCommand(clientData, interp, 
						   argc, argv, theTclBuilder, theTclDomain);
}

/// added by ZHY
extern int 
TclCommand_UpdateMaterialsCommand(ClientData clientData, 
				  Tcl_Interp *interp, 
				  int argc, 
				  TCL_Char **argv, 
				  TclModelBuilder *theTclBuilder,
				  Domain *theDomain);
int
TclCommand_UpdateMaterials(ClientData clientData, 
			   Tcl_Interp *interp,  
			   int argc, 
			   TCL_Char **argv)
{
  return TclCommand_UpdateMaterialsCommand(clientData, interp, 
					   argc, argv, theTclBuilder, theTclDomain);
}

/// added by ZHY
extern int 
TclModelBuilderUpdateParameterCommand(ClientData clientData, 
					  Tcl_Interp *interp, 
					  int argc, 
					  TCL_Char **argv, 
					  TclModelBuilder *theTclBuilder);
int
TclCommand_UpdateParameter(ClientData clientData, 
				    Tcl_Interp *interp,  
				    int argc, 
				    TCL_Char **argv)
{
  return TclModelBuilderUpdateParameterCommand(clientData, interp, 
				       argc, argv, theTclBuilder);
}

extern int
TclModelBuilderFrictionModelCommand (ClientData clienData,
				     Tcl_Interp *interp, int argc, TCL_Char **argv,
				     Domain *theDomain);

int
TclCommand_addFrictionModel(ClientData clientData,
                    Tcl_Interp *interp, int argc, TCL_Char **argv)                      
{
  return TclModelBuilderFrictionModelCommand(clientData, interp, argc, argv, theTclDomain);
}

int 
TclCommand_Package(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  
  void *libHandle;
  int (*funcPtr)(ClientData clientData, Tcl_Interp *interp,  int argc, 
		 TCL_Char **argv, Domain*, TclModelBuilder*);       
  
  const char *funcName = 0;
  int res = -1;
  
  if (argc == 2) {
    res = getLibraryFunction(argv[1], argv[1], &libHandle, (void **)&funcPtr);
  } else if (argc == 3) {
    res = getLibraryFunction(argv[1], argv[2], &libHandle, (void **)&funcPtr);
  }

  if (res == 0) {
    int result = (*funcPtr)(clientData, interp,
			    argc, 
			    argv,
			    theTclDomain,
			    theTclBuilder);	
  } else {
    opserr << "Error: Could not find function: " << argv[1] << endln;
    return -1;
  }

  return res;
}



// Added by Alborz Ghofrani - U.Washington
extern int
TclCommand_GenerateInterfacePoints(ClientData clientData, Tcl_Interp *interp, int argc,
 		TCL_Char **argv);
// End Added by Alborz
