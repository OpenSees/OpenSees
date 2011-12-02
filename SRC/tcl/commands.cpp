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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/commands.cpp,v $
                                                                        
                                                                        
// File: ~/tcl/commands.C
// 
// Written: fmk 
// Created: 04/98
// Revision: A
//
// Description: This file contains the functions that will be called by
// the interpreter when the appropriate command name is specified,
// see tkAppInit.C for command names.
//
// What: "@(#) commands.C, revA"

extern "C" {
#include <tcl.h>
#include <tk.h>
}

#include <G3Globals.h>

#include <RigidRod.h>
#include <RigidBeam.h>
#include <RigidDiaphragm.h>

#include <iostream.h>
#include <fstream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Timer.h>

#include <ModelBuilder.h>
#include "commands.h"

// domain
#include <Domain.h>
#include <Element.h>
#include <Node.h>
#include <ElementIter.h>
#include <NodeIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <ElementalLoad.h>
#include <ElementalLoadIter.h>


// analysis model
#include <AnalysisModel.h>

// convergence tests
#include <CTestNormUnbalance.h>
#include <CTestNormDispIncr.h>
#include <CTestEnergyIncr.h>

// soln algorithms
#include <Linear.h>
#include <NewtonRaphson.h>
#include <ModifiedNewton.h>
#include <FrequencyAlgo.h>

// constraint handlers
#include <PlainHandler.h>
#include <PenaltyConstraintHandler.h>
#include <LagrangeConstraintHandler.h>
#include <TransformationConstraintHandler.h>

// numberers
#include <PlainNumberer.h>
#include <DOF_Numberer.h>

// integrators
#include <LoadControl.h>
#include <ArcLength.h>
#include <ArcLength1.h>
#include <MinUnbalDispNorm.h>
#include <DisplacementControl.h>
#include <Newmark.h>
#include <WilsonTheta.h>
#include <HHT.h>
#include <HHT1.h>
#include <Newmark1.h> 
#include <EigenIntegrator.h>

// analysis
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <EigenAnalysis.h>

// system of eqn and solvers
// #include <SlowLinearSOE.h>
// #include <SlowLinearSOESolver.h>
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>



#include <FullGenLinSOE.h>
#include <FullGenLinLapackSolver.h>
#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>
// #include <ProfileSPDLinDirectBlockSolver.h>
// #include <ProfileSPDLinDirectThreadSolver.h>
// #include <ProfileSPDLinDirectSkypackSolver.h>
// #include <BandSPDLinThreadSolver.h>
#include <SparseGenColLinSOE.h>
#include <SuperLU.h>
#include <SymSparseLinSOE.h>
#include <SymSparseLinSolver.h>

#include <SparseGenColLinSOE.h>
#include <SuperLU.h>
#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>

#include <EigenSOE.h>
#include <EigenSolver.h>
#include <SymArpackSOE.h>
#include <SymArpackSolver.h>
#include <BandArpackSOE.h>
#include <BandArpackSolver.h>

// graph
#include <RCM.h>

#include <ErrorHandler.h>
#include <ConsoleErrorHandler.h>

#include <TclVideoPlayer.h>

ModelBuilder *theBuilder =0;
Domain theDomain;
static AnalysisModel *theAnalysisModel =0;
static EquiSolnAlgo *theAlgorithm =0;
static ConstraintHandler *theHandler =0;
static DOF_Numberer *theNumberer =0;
static LinearSOE *theSOE =0;
static StaticAnalysis *theStaticAnalysis = 0;
static DirectIntegrationAnalysis *theTransientAnalysis = 0;
static StaticIntegrator *theStaticIntegrator =0;
static TransientIntegrator *theTransientIntegrator =0;
static ConvergenceTest *theTest =0;
static bool builtModel = false;

static EigenAnalysis *theEigenAnalysis = 0;

ErrorHandler *g3ErrorHandler =0;
TclVideoPlayer *theTclVideoPlayer =0;

// g3AppInit() is the method called by tkAppInit() when the
// interpreter is being set up .. this is where all the
// commands defined in this file are registered with the interpreter.

extern int myCommands(Tcl_Interp *interp);

int g3AppInit(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "wipe", wipeModel,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
    Tcl_CreateCommand(interp, "wipeAnalysis", wipeAnalysis,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
    Tcl_CreateCommand(interp, "reset", resetModel,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);        
    Tcl_CreateCommand(interp, "initialize", initializeAnalysis,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);        
    Tcl_CreateCommand(interp, "loadConst", setLoadConst,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);            
    Tcl_CreateCommand(interp, "time", setTime,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);     
    Tcl_CreateCommand(interp, "build", buildModel,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "analyze", analyzeModel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "print", printModel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "analysis", specifyAnalysis, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "system", specifySOE, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "numberer", specifyNumberer, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "constraints", specifyConstraintHandler, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "algorithm", specifyAlgorithm, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "test", specifyCTest, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
    Tcl_CreateCommand(interp, "integrator", specifyIntegrator, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "recorder", addRecorder, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "algorithmRecorder", addAlgoRecorder, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "database", addDatabase, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "playback", playbackRecorders, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);        
    Tcl_CreateCommand(interp, "rigidLink", rigidLink, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);                
    Tcl_CreateCommand(interp, "rigidDiaphragm", rigidDiaphragm, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);   
    Tcl_CreateCommand(interp, "eigen", eigenAnalysis, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "video", videoPlayer, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "remove", removeObject, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       

    theAlgorithm =0;
    theHandler =0;
    theNumberer =0;
    theAnalysisModel =0;  
    theSOE =0;
    theStaticIntegrator =0;
    theTransientIntegrator =0;
    theStaticAnalysis =0;
    theTransientAnalysis =0;    
    theTest = 0;

    // create an error handler
    g3ErrorHandler = new ConsoleErrorHandler();
    theTclVideoPlayer = 0;

    return myCommands(interp);
}


int 
wipeModel(ClientData clientData, Tcl_Interp *interp, int argc, 
	   char **argv)
{
  // to build the model make sure the ModelBuilder has been constructed
  // and that the model has not already been constructed
  if (theBuilder != 0) {
    delete theBuilder;
    builtModel = false;
    theBuilder = 0;
  }

  if (theStaticAnalysis != 0) {
      theStaticAnalysis->clearAll();
      delete theStaticAnalysis;
  }
  
  if (theTransientAnalysis != 0) {
      theTransientAnalysis->clearAll();
      delete theTransientAnalysis;  
  }

  /*
  if (theEigenAnalysis != 0) {
    delete theEigenAnalysis;
    theEigenAnalysis = 0;
  }
  */

  theDomain.clearAll();

  if (theTest != 0)
      delete theTest;
  
  if (theTclVideoPlayer != 0) {
	  delete theTclVideoPlayer;
	  theTclVideoPlayer = 0;
  }

  theAlgorithm =0;
  theHandler =0;
  theNumberer =0;
  theAnalysisModel =0;  
  theSOE =0;
  theStaticIntegrator =0;
  theTransientIntegrator =0;
  theStaticAnalysis =0;
  theTransientAnalysis =0;    
  theTest = 0;
  
  // the domain deletes the record objects, 
  // just have to delete the private array
  return TCL_OK;  
}

int 
wipeAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, 
	   char **argv)
{


  if (theStaticAnalysis != 0) {
      theStaticAnalysis->clearAll();
      delete theStaticAnalysis;
  }
  
  if (theTransientAnalysis != 0) {
      theTransientAnalysis->clearAll();
      delete theTransientAnalysis;  
  }

  if (theTest != 0)
      delete theTest;
  
  theAlgorithm =0;
  theHandler =0;
  theNumberer =0;
  theAnalysisModel =0;  
  theSOE =0;
  theStaticIntegrator =0;
  theTransientIntegrator =0;
  theStaticAnalysis =0;
  theTransientAnalysis =0;    
  theTest = 0;
  
  // the domain deletes the record objects, 
  // just have to delete the private array
  return TCL_OK;  
}


int 
resetModel(ClientData clientData, Tcl_Interp *interp, int argc, 
	   char **argv)
{
  theDomain.revertToStart();
  return TCL_OK;
}

int
initializeAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, 
	   char **argv)
{
  if (theTransientAnalysis != 0)
    theTransientAnalysis->initialize();

  return TCL_OK;
}


int 
setLoadConst(ClientData clientData, Tcl_Interp *interp, int argc, 
	   char **argv)
{
  theDomain.setLoadConstant();
  if (argc == 3) {
      if( strcmp(argv[1],"-time")) {
	  double newTime;
	  if (Tcl_GetDouble(interp, argv[2], &newTime) != TCL_OK) {
	      cerr << "WARNING readingvalue - loadConst -time value ";
	      return TCL_ERROR;
	  } else
	      theDomain.setCurrentTime(newTime);
      }    	  
  }
	  
  return TCL_OK;
}

int 
setTime(ClientData clientData, Tcl_Interp *interp, int argc, 
	   char **argv)
{
  theDomain.setLoadConstant();
  if (argc < 2) {
      cerr << "WARNING illegal command - time pseudoTime? ";
      return TCL_ERROR;
  }
  double newTime;
  if (Tcl_GetDouble(interp, argv[1], &newTime) != TCL_OK) {
      cerr << "WARNING reading time value - time pseudoTime? ";
      return TCL_ERROR;
  } else
      theDomain.setCurrentTime(newTime);
  return TCL_OK;
}


// command invoked to build the model, i.e. to invoke buildFE_Model() 
// on the ModelBuilder
int 
buildModel(ClientData clientData, Tcl_Interp *interp, int argc, 
	   char **argv)
{
    // to build the model make sure the ModelBuilder has been constructed
  // and that the model has not already been constructed
  if (theBuilder != 0 && builtModel == false) {
    builtModel = true;
    return theBuilder->buildFE_Model();
  }  else if (theBuilder != 0 && builtModel == true) {
      interp->result = "WARNING Model has already been built - not built again ";
      return TCL_ERROR;
  }
  else {
      interp->result = "WARNING No ModelBuilder type has been specified ";
      return TCL_ERROR;
  }    
}


//
// command invoked to build the model, i.e. to invoke analyze() 
// on the Analysis object
//
int 
analyzeModel(ClientData clientData, Tcl_Interp *interp, int argc, 
	     char **argv)
{

  if (theStaticAnalysis != 0) {
    if (argc < 2) {
      interp->result = "WARNING static analysis: analysis numIncr?";
      return TCL_ERROR;
    }
    int numIncr;
    if (Tcl_GetInt(interp, argv[1], &numIncr) != TCL_OK)	
      return TCL_ERROR;	      
    return theStaticAnalysis->analyze(numIncr);
  } else if (theTransientAnalysis != 0) {
    if (argc < 3) {
      interp->result = "WARNING transient analysis: analysis numIncr? deltaT?";
      return TCL_ERROR;
    }
    int numIncr;
    if (Tcl_GetInt(interp, argv[1], &numIncr) != TCL_OK)	
      return TCL_ERROR;
    double dT;
    if (Tcl_GetDouble(interp, argv[2], &dT) != TCL_OK)	
      return TCL_ERROR;
    return theTransientAnalysis->analyze(numIncr, dT);

  } else {
    interp->result = "WARNING No Analysis type has been specified ";
    return TCL_ERROR;
  }    
}


int 
printElement(ClientData clientData, Tcl_Interp *interp, int argc, 
	  char **argv, int nodeArg, ostream &output);


int 
printNode(ClientData clientData, Tcl_Interp *interp, int argc, 
	  char **argv, int nodeArg, ostream &output);
	  
int 
printIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, 
		char **argv, int nodeArg, ostream &output);	  
		
int 
printAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, 
		char **argv, int nodeArg, ostream &output);	  		

int 
printModel(ClientData clientData, Tcl_Interp *interp, int argc, 
			 char **argv)
{
  // if just 'print' then print out the entire domain
  if (argc == 1) {
    cerr << theDomain;
    return TCL_OK;
  }    

  // if 'print ele i j k..' print out some elements
  if ((strcmp(argv[1],"-ele") == 0) || (strcmp(argv[1],"ele") == 0))
    return printElement(clientData, interp, argc, argv, 3, cerr);    

  // if 'print node i j k ..' print out some nodes
  else if ((strcmp(argv[1],"-node") == 0) || (strcmp(argv[1],"node") == 0)) 
      return printNode(clientData, interp, argc, argv, 3, cerr);
  
  // if 'print integrator flag' print out the integrator
  else if ((strcmp(argv[1],"integrator") == 0) || 
	   (strcmp(argv[1],"-integrator") == 0)) 
    return printIntegrator(clientData, interp, argc, argv, 3, cerr);  
  
  // if 'print algorithm flag' print out the algorithm
  else if ((strcmp(argv[1],"algorithm") == 0) || 
	   (strcmp(argv[1],"-algorithm") == 0))
    return printAlgorithm(clientData, interp, argc, argv, 3, cerr);    

  else { // it must be a file we are going to print to
    ofstream output(argv[1],ios::app); // open for appending to
    if (!output) {
      cerr << "print <filename> .. - failed to open file: " << argv[1] << endl;
      return TCL_ERROR;
    }
    // if just 'print <filename>' then print out the entire domain to eof
    if (argc == 2) {
      output << theDomain;
      return TCL_OK;
    }    

    int pos = 2;
    if ((strcmp(argv[pos],"string") == 0) || 
	(strcmp(argv[pos],"-string") == 0)) {
	output << argv[3] << endl;
	pos +=2;
    }
    int res = TCL_OK;    

    // if 'print <filename> ele i j k..' print out some elements
    if ((strcmp(argv[pos],"ele") == 0) || 
	(strcmp(argv[pos],"-ele") == 0))
      res = printElement(clientData, interp, argc, argv, pos+2, output);    

    // if 'print <filename> node i j k ..' print out some nodes
    else if ((strcmp(argv[pos],"node") == 0) || (strcmp(argv[pos],"-node") == 0))
      res = printNode(clientData, interp, argc, argv, pos+2, output);
    
    // if 'print integrator flag' print out the integrator
    else if ((strcmp(argv[pos],"integrator") == 0) 
	     || (strcmp(argv[pos],"-integrator") == 0))
	return printIntegrator(clientData, interp, argc, argv, pos+2, cerr);  
  
    // if 'print algorithm flag' print out the algorithm
    else if ((strcmp(argv[pos],"-algorithm") == 0)|| 
	     (strcmp(argv[pos],"algorithm") == 0))
	return printAlgorithm(clientData, interp, argc, argv, pos+2, cerr);    


    // close the output file
    output.close();
    return res;
  }
  
}


// printNode():
// function to print out the nodal information conatined in line
//     print <filename> node <flag int> <int int int>
// input: nodeArg: integer equal to arg count to node plus 1
//        output: output stream to which the results are sent
// 
int 
printNode(ClientData clientData, Tcl_Interp *interp, int argc, 
	  char **argv, int nodeArg, ostream &output)
{
  int flag = 0; // default flag sent to a nodes Print() method

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc < nodeArg) { 
      NodeIter &theNodes = theDomain.getNodes();
      Node *theNode;
      while ((theNode = theNodes()) != 0)
	theNode->Print(output);
      return TCL_OK;
  }    

  // if 'print <filename> node flag int <int int ..>' get the flag
  if ((strcmp(argv[nodeArg-1],"flag") == 0) ||
      (strcmp(argv[nodeArg-1],"-flag") == 0)) { 
      // get the specified flag
    if (argc <= nodeArg) {
      cerr << "WARNING print <filename> node <flag int> no int specified \n";
      return TCL_ERROR;
    }    
    if (Tcl_GetInt(interp, argv[nodeArg], &flag) != TCL_OK) {
      cerr << "WARNING print node failed to get integer flag: ";
      cerr << argv[nodeArg] << endl; 
      return TCL_ERROR;
    }    
    nodeArg += 2;
  }

  // now print the nodes with the specified flag, 0 by default

  // if 'print <filename> node flag' 
  //     print out all the nodes in the domain with flag
  if (argc < nodeArg) { 
    NodeIter &theNodes = theDomain.getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0)
      theNode->Print(output, flag);
    return TCL_OK;
  } else { 
    // otherwise print out the specified nodes i j k .. with flag
    for (int i= nodeArg-1; i<argc; i++) {
      int nodeTag;
      if (Tcl_GetInt(interp, argv[i], &nodeTag) != TCL_OK) {
	cerr << "WARNING print node failed to get integer: " << argv[i] << endl;
	return TCL_ERROR;
      }
      Node *theNode = theDomain.getNode(nodeTag);
      if (theNode != 0)
	theNode->Print(output,flag);
    }
    return TCL_OK;
  }
}


int 
printElement(ClientData clientData, Tcl_Interp *interp, int argc, 
	  char **argv, int eleArg, ostream &output)
{
  int flag = 0; // default flag sent to a nodes Print() method

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc < eleArg) { 
      ElementIter &theElements = theDomain.getElements();
      Element *theElement;
      while ((theElement = theElements()) != 0)
	theElement->Print(output);
      return TCL_OK;
  }    

  // if 'print <filename> Element flag int <int int ..>' get the flag
  if ((strcmp(argv[eleArg-1],"flag") == 0) ||
      (strcmp(argv[eleArg-1],"-flag")) == 0) { // get the specified flag
    if (argc <= eleArg) {
      cerr << "WARNING print <filename> ele <flag int> no int specified \n";
      return TCL_ERROR;
    }    
    if (Tcl_GetInt(interp, argv[eleArg], &flag) != TCL_OK) {
      cerr << "WARNING print ele failed to get integer flag: ";
      cerr << argv[eleArg] << endl; 
      return TCL_ERROR;
    }    
    eleArg += 2;
  }

  // now print the Elements with the specified flag, 0 by default

  // if 'print <filename> Element flag' 
  //     print out all the Elements in the domain with flag
  if (argc < eleArg) { 
    ElementIter &theElements = theDomain.getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0)      
      theElement->Print(output, flag);
    return TCL_OK;
  } else { 
    // otherwise print out the specified Elements i j k .. with flag
    for (int i= eleArg-1; i<argc; i++) {
      int ElementTag;
      if (Tcl_GetInt(interp, argv[i], &ElementTag) != TCL_OK) {
	cerr << "WARNING print ele failed to get integer: " << argv[i] << endl;
	return TCL_ERROR;
      }
      Element *theElement = theDomain.getElement(ElementTag);
      if (theElement != 0)
	theElement->Print(output,flag);
    }
    return TCL_OK;
  }
}


int 
printAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, 
	       char **argv, int eleArg, ostream &output)
{

  if (theAlgorithm == 0)
      return TCL_OK;

  // if just 'print <filename> algorithm'- no flag
  if (argc < eleArg) { 
      theAlgorithm->Print(output);
      return TCL_OK;
  }    

  // if 'print <filename> Algorithm flag' get the flag
  int flag;  
  if (Tcl_GetInt(interp, argv[eleArg-1], &flag) != TCL_OK) {  
      cerr << "WARNING print algorithm failed to get integer flag: ";
      cerr << argv[eleArg] << endl; 
      return TCL_ERROR;
  }    
  theAlgorithm->Print(output,flag);
  return TCL_OK;  
}


int 
printIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, 
	       char **argv, int eleArg, ostream &output)
{

  if (theStaticIntegrator == 0 && theTransientIntegrator == 0)
      return TCL_OK;
  
  IncrementalIntegrator *theIntegrator;
  if (theStaticIntegrator != 0)
      theIntegrator = theStaticIntegrator;
  else
      theIntegrator = theTransientIntegrator;

  // if just 'print <filename> algorithm'- no flag
  if (argc < eleArg) { 
      theIntegrator->Print(output);
      return TCL_OK;
  }    

  // if 'print <filename> Algorithm flag' get the flag
  int flag;  
  if (Tcl_GetInt(interp, argv[eleArg-1], &flag) != TCL_OK) {  
      cerr << "WARNING print algorithm failed to get integer flag: ";
      cerr << argv[eleArg] << endl; 
      return TCL_ERROR;
  }    
  theIntegrator->Print(output,flag);
  return TCL_OK;  
}


//
// command invoked to allow the Analysis object to be built
//
int 
specifyAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, 
		char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	interp->result = "WARNING need to specify an analysis type (Static, Transient)";
	return TCL_ERROR;
    }    

    // delete the old analysis
    if (theStaticAnalysis != 0)
	delete theStaticAnalysis;
    if (theTransientAnalysis != 0)
	delete theTransientAnalysis;
    
    // check argv[1] for type of SOE and create it
    if (strcmp(argv[1],"Static") == 0) {
	// make sure all the components have been built,
	// otherwise print a warning and use some defaults
	if (theAnalysisModel == 0) 
	    theAnalysisModel = new AnalysisModel();
	
	if (theAlgorithm == 0) {
	    cerr << "WARNING analysis Static - no Algorithm yet specified, ";
	    cerr << " NewtonRaphson default will be used\n";	    

	    if (theTest == 0) 
		theTest = new CTestNormUnbalance(1.0e-6,25,0);       
	    theAlgorithm = new NewtonRaphson(*theTest); 
	}
	if (theHandler == 0) {
	    cerr << "WARNING analysis Static - no ConstraintHandler yet specified, ";
	    cerr << " PlainHandler default will be used\n";
	    theHandler = new PlainHandler();       
	}
	if (theNumberer == 0) {
	    cerr << "WARNING analysis Static - no Numberer specified, ";
	    cerr << " RCM default will be used\n";
	    RCM *theRCM = new RCM();	
	    theNumberer = new DOF_Numberer(*theRCM);    	
	}
	if (theStaticIntegrator == 0) {
	    cerr << "WARNING analysis Static - no Integrator specified, ";
	    cerr << " StaticIntegrator default will be used\n";
	    theStaticIntegrator = new LoadControl(1, 1, 1, 1);       
	}
	if (theSOE == 0) {
	    cerr << "WARNING analysis Static - no LinearSOE specified, ";
	    cerr << " ProfileSPDLinSOE default will be used\n";
	    ProfileSPDLinSolver *theSolver;
	    theSolver = new ProfileSPDLinDirectSolver(); 	
	    theSOE = new ProfileSPDLinSOE(*theSolver);      
	}
    
	theStaticAnalysis = new StaticAnalysis(theDomain,
					       *theHandler,
					       *theNumberer,
					       *theAnalysisModel,
					       *theAlgorithm,
					       *theSOE,
					       *theStaticIntegrator);

    } else if (strcmp(argv[1],"Transient") == 0) {
	// make sure all the components have been built,
	// otherwise print a warning and use some defaults
	if (theAnalysisModel == 0) 
	    theAnalysisModel = new AnalysisModel();
	
	if (theAlgorithm == 0) {
	    cerr << "WARNING analysis Transient - no Algorithm yet specified, ";
	    cerr << " NewtonRaphson default will be used\n";	    

	    if (theTest == 0) 
		theTest = new CTestNormUnbalance(1.0e-6,25,0);       
	    theAlgorithm = new NewtonRaphson(*theTest); 
	}
	if (theHandler == 0) {
	    cerr << "WARNING analysis Transient dt tFinal - no ConstraintHandler";
	    cerr << " yet specified, PlainHandler default will be used\n";
	    theHandler = new PlainHandler();       
	}
	if (theNumberer == 0) {
	    cerr << "WARNING analysis Transient dt tFinal - no Numberer specified, ";
	    cerr << " RCM default will be used\n";
	    RCM *theRCM = new RCM();	
	    theNumberer = new DOF_Numberer(*theRCM);    	
	}
	if (theTransientIntegrator == 0) {
	    cerr << "WARNING analysis Transient dt tFinal - no Integrator specified, ";
	    cerr << " Newmark(.5,.25) default will be used\n";
	    theTransientIntegrator = new Newmark(0.5,0.25);       
	}
	if (theSOE == 0) {
	    cerr << "WARNING analysis Transient dt tFinal - no LinearSOE specified, ";
	    cerr << " ProfileSPDLinSOE default will be used\n";
	    ProfileSPDLinSolver *theSolver;
	    theSolver = new ProfileSPDLinDirectSolver(); 	
	    theSOE = new ProfileSPDLinSOE(*theSolver);      
	}
    
	theTransientAnalysis = new DirectIntegrationAnalysis(theDomain,
							     *theHandler,
							     *theNumberer,
							     *theAnalysisModel,
							     *theAlgorithm,
							     *theSOE,
							     *theTransientIntegrator);
	
    } else {
	interp->result = "WARNING No Analysis type exists (Static Transient only) ";
	return TCL_ERROR;
    }
    return TCL_OK;
}


//
// command invoked to allow the SystemOfEqn and Solver objects to be built
//
int 
specifySOE(ClientData clientData, Tcl_Interp *interp, int argc, 
		    char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 2) {
      interp->result = "WARNING need to specify a model type ";
      return TCL_ERROR;
  }    

  // check argv[1] for type of SOE and create it
  // BAND GENERAL SOE & SOLVER
  if (strcmp(argv[1],"BandGeneral") == 0) {
      BandGenLinSolver    *theSolver = new BandGenLinLapackSolver();
      theSOE = new BandGenLinSOE(*theSolver);      
  } 

  // BAND SPD SOE & SOLVER
  else if (strcmp(argv[1],"BandSPD") == 0) {
      BandSPDLinSolver    *theSolver = new BandSPDLinLapackSolver();   
      theSOE = new BandSPDLinSOE(*theSolver);        
  } 

  // PROFILE SPD SOE * SOLVER
  else if (strcmp(argv[1],"ProfileSPD") == 0) {
    // now must determine the type of solver to create from rest of args
    ProfileSPDLinSolver *theSolver = new ProfileSPDLinDirectSolver(); 	

    /* *********** Some misc solvers i play with ******************
    else if (strcmp(argv[2],"Normal") == 0) {
      theSolver = new ProfileSPDLinDirectSolver(); 	
    } 

    else if (strcmp(argv[2],"Block") == 0) {  
      int blockSize = 4;
      if (argc == 4) {
	if (Tcl_GetInt(interp, argv[3], &blockSize) != TCL_OK)
	  return TCL_ERROR;
      }
      theSolver = theSolver = new ProfileSPDLinDirectBlockSolver(1.0e-12,blockSize); 
    }

    
      int blockSize = 4;
      int numThreads = 1;
      if (argc == 5) {
	if (Tcl_GetInt(interp, argv[3], &blockSize) != TCL_OK)
	  return TCL_ERROR;
	if (Tcl_GetInt(interp, argv[4], &numThreads) != TCL_OK)
	  return TCL_ERROR;
      }
      theSolver = new ProfileSPDLinDirectThreadSolver(numThreads,blockSize,1.0e-12); 
      } else if (strcmp(argv[2],"Thread") == 0) {  
      int blockSize = 4;
      int numThreads = 1;
      if (argc == 5) {
	if (Tcl_GetInt(interp, argv[3], &blockSize) != TCL_OK)
	  return TCL_ERROR;
	if (Tcl_GetInt(interp, argv[4], &numThreads) != TCL_OK)
	  return TCL_ERROR;
      }
      theSolver = new ProfileSPDLinDirectThreadSolver(numThreads,blockSize,1.0e-12); 
    } 
    else if (strcmp(argv[2],"Skypack") == 0) {  
      if (argc == 5) {
	int mCols, mRows;
	if (Tcl_GetInt(interp, argv[3], &mCols) != TCL_OK)
	  return TCL_ERROR;
	if (Tcl_GetInt(interp, argv[4], &mRows) != TCL_OK)
	  return TCL_ERROR;
	theSolver = new ProfileSPDLinDirectSkypackSolver(mCols, mRows); 
      } else 
	theSolver = new ProfileSPDLinDirectSkypackSolver(); 	
    }
    else 
      theSolver = new ProfileSPDLinDirectSolver(); 	
    ***************************************************************  */
    
      theSOE = new ProfileSPDLinSOE(*theSolver);      
   }


  // SPARSE GENERAL SOE * SOLVER
  else if (strcmp(argv[1],"SparseGeneral") == 0) {
    // now must determine the type of solver to create from rest of args
    SparseGenColLinSolver *theSolver;
    if (argc == 3) {
	if (strcmp(argv[2],"p") == 0 || strcmp(argv[2],"piv"))
	    theSolver = new SuperLU(0,1.0); 	
        else
	    theSolver = new SuperLU(); 		    
    } else 
	theSolver = new SuperLU(); 		        

    theSOE = new SparseGenColLinSOE(*theSolver);      
  }	

  else if (strcmp(argv[1],"SparseSPD") == 0) {
    // now must determine the type of solver to create from rest of args
    SymSparseLinSolver *theSolver = new SymSparseLinSolver();
    theSOE = new SymSparseLinSOE(*theSolver);      
  }	
  
  else if (strcmp(argv[1],"UmfPack") == 0) {
    // now must determine the type of solver to create from rest of args
      UmfpackGenLinSolver *theSolver = new UmfpackGenLinSolver();
      theSOE = new UmfpackGenLinSOE(*theSolver);      
  }	  


  else {
    interp->result = "WARNING No SystemOfEqn type exists) ";
    return TCL_ERROR;
  }

  // if the analysis exists - we want to change the SOE
  if (theStaticAnalysis != 0)
    theStaticAnalysis->setLinearSOE(*theSOE);

  return TCL_OK;
}




//
// command invoked to allow the Numberer objects to be built
//
int 
specifyNumberer(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv)
{
  // make sure at least one other argument to contain numberer
  if (argc < 2) {
      interp->result = "WARNING need to specify a Nemberer type ";
      return TCL_ERROR;
  }    

  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1],"Plain") == 0) 
    theNumberer = new PlainNumberer();       
  else if (strcmp(argv[1],"RCM") == 0) {
    RCM *theRCM = new RCM();	
    theNumberer = new DOF_Numberer(*theRCM);    	
  }
  else {
    interp->result = "WARNING No Numberer type exists (Plain, RCM only) ";
    return TCL_ERROR;
  }    
  return TCL_OK;
}




//
// command invoked to allow the ConstraintHandler object to be built
//
int 
specifyConstraintHandler(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv)
{
  // make sure at least one other argument to contain numberer
  if (argc < 2) {
      interp->result = "WARNING need to specify a Nemberer type ";
      return TCL_ERROR;
  }    

  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1],"Plain") == 0) 
    theHandler = new PlainHandler();       

  else if (strcmp(argv[1],"Penalty") == 0) {
    if (argc < 4) {
      interp->result = "WARNING: need to specify alpha: handler Penalty alpha ";
      return TCL_ERROR;
    }    
    double alpha1, alpha2;
    if (Tcl_GetDouble(interp, argv[2], &alpha1) != TCL_OK)	
      return TCL_ERROR;	
    if (Tcl_GetDouble(interp, argv[3], &alpha2) != TCL_OK)	
      return TCL_ERROR;	
    theHandler = new PenaltyConstraintHandler(alpha1, alpha2);
  }
  
  else if (strcmp(argv[1],"Lagrange") == 0) {
    double alpha1 = 1.0;
    double alpha2 = 1.0;
    if (argc == 4) {
      if (Tcl_GetDouble(interp, argv[2], &alpha1) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[3], &alpha2) != TCL_OK)	
	return TCL_ERROR;	
    }
    theHandler = new LagrangeConstraintHandler(alpha1, alpha2);
  }  
  
  else if (strcmp(argv[1],"Transformation") == 0) {
    theHandler = new TransformationConstraintHandler();
  }    

  else {
    cerr << "WARNING No ConstraintHandler type exists (Plain, Penalty,";
    cerr << " Lagrange, Transformation) only\n";
    return TCL_ERROR;
  }    
  return TCL_OK;
}



//
// command invoked to allow the SolnAlgorithm object to be built
//
int
specifyAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv)
{
  // make sure at least one other argument to contain numberer
  if (argc < 2) {
      interp->result = "WARNING need to specify an Algorithm type ";
      return TCL_ERROR;
  }    

  // check argv[1] for type of Algorithm and create the object
  if (strcmp(argv[1],"Linear") == 0) 
    theAlgorithm = new Linear();       

  else if (strcmp(argv[1],"Newton") == 0) {
      if (theTest == 0) {
	  interp->result = "ERROR: No ConvergenceTest yet specified\n";
	  return TCL_ERROR;	  
      }
      theAlgorithm = new NewtonRaphson(*theTest); 
  }
  
  else if (strcmp(argv[1],"ModifiedNewton") == 0) {
      if (theTest == 0) {
	  interp->result = "ERROR: No ConvergenceTest yet specified\n";
	  return TCL_ERROR;	  
      }
      theAlgorithm = new ModifiedNewton(*theTest); 
  }  

  else {
      interp->result = "WARNING No EquiSolnAlgo type exists (Linear, Newton only) ";
      return TCL_ERROR;
  }    

  // if the analysis exists - we want to change the SOE
  if (theStaticAnalysis != 0)
    theStaticAnalysis->setAlgorithm(*theAlgorithm);
  else if (theTransientAnalysis != 0)
    theTransientAnalysis->setAlgorithm(*theAlgorithm);  
  
  return TCL_OK;
}


//
// command invoked to allow the SolnAlgorithm object to be built
//
int
specifyCTest(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv)
{
  // make sure at least one other argument to contain numberer
  if (argc < 2) {
      interp->result = "WARNING need to specify a ConvergenceTest Type type ";
      return TCL_ERROR;
  }    

  // get the tolerence first
  double tol;
  int numIter;
  int printIt = 0;
  
  if (argc == 3) {
      if (Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
	  return TCL_ERROR;		
      numIter = 25;      
  } else if (argc == 4) {
      if (Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
	  return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)	
	  return TCL_ERROR;			  
  } else if (argc == 5) {
      if (Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
	  return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)	
	  return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[4], &printIt) != TCL_OK)	
	  return TCL_ERROR;			  
  }  else {
     tol = 1.0e-6;
     numIter = 25;      
  }
      
  if (strcmp(argv[1],"NormUnbalance") == 0) 
    theTest = new CTestNormUnbalance(tol,numIter,printIt);       
  else if (strcmp(argv[1],"NormDispIncr") == 0) 
    theTest = new CTestNormDispIncr(tol,numIter,printIt);             
  else if (strcmp(argv[1],"EnergyIncr") == 0) 
    theTest = new CTestEnergyIncr(tol,numIter,printIt);             
  else {
    cerr << "WARNING No ConvergenceTest type (NormUnbalance, NormDispIncr, EnergyIncr)"; 
    return TCL_ERROR;
  }    

  // if the algorithm exists - we want to change the test
//  if (theAlgorithm != 0)
//    theAlgorithm->setTest(*theTest);  
  return TCL_OK;
}





//
// command invoked to allow the Integrator object to be built
//
int 
specifyIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, 
		      char **argv)
{
  // make sure at least one other argument to contain integrator
  if (argc < 2) {
      interp->result = "WARNING need to specify an Integrator type ";
      return TCL_ERROR;
  }    

  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1],"LoadControl") == 0) {
      double dLambda;
      double minIncr, maxIncr;
      int numIter;
      if (argc < 6) {
	cerr << argc; for (int j=0; j<argc; j++) cerr << argv[j] << endl;
	interp->result = "WARNING incorrect # args - integrator LoadControl alpha Jd minAlpha maxAlpha";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &dLambda) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[4], &minIncr) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[5], &maxIncr) != TCL_OK)	
	return TCL_ERROR;	      
      theStaticIntegrator = new LoadControl(dLambda, numIter, minIncr, maxIncr);       
  }

  
  else if (strcmp(argv[1],"ArcLength") == 0) {
      double arcLength;
      double alpha;
      if (argc != 4) {
	interp->result = "WARNING integrator ArcLength arcLength alpha ";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK)	
	return TCL_ERROR;	
      theStaticIntegrator = new ArcLength(arcLength,alpha);       
  }

  else if (strcmp(argv[1],"ArcLength1") == 0) {
      double arcLength;
      double alpha;
      if (argc != 4) {
	interp->result = "WARNING integrator ArcLength1 arcLength alpha ";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK)	
	return TCL_ERROR;	
      theStaticIntegrator = new ArcLength1(arcLength,alpha);       
  }

  else if (strcmp(argv[1],"MinUnbalDispNorm") == 0) {
      double lambda11, minlambda, maxlambda;
      int numIter;
      if (argc != 6) {
	cerr << "WARNING integrator MinUnbalDispNorm lambda11 Jd minLambda1j maxLambda1j\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &lambda11) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[4], &minlambda) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[5], &maxlambda) != TCL_OK)	
	return TCL_ERROR;	
      theStaticIntegrator = new MinUnbalDispNorm(lambda11,numIter,minlambda,maxlambda);
  }
  
  else if (strcmp(argv[1],"DisplacementControl") == 0) {
      int node;
      int dof;
      double increment, minIncr, maxIncr;
      int numIter;
      if (argc != 8) {
	cerr << "WARNING integrator DisplacementControl node dof dU ";
	cerr << "Jd minIncrement maxIncrement\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetInt(interp, argv[2], &node) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[4], &increment) != TCL_OK)	
	return TCL_ERROR;	      
      if (Tcl_GetInt(interp, argv[5], &numIter) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[6], &minIncr) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[7], &maxIncr) != TCL_OK)	
	return TCL_ERROR;	      


      theStaticIntegrator = new DisplacementControl(node,dof-1,increment, &theDomain,
						    numIter, minIncr, maxIncr);
  }  
  
  else if (strcmp(argv[1],"Newmark") == 0) {
      double gamma;
      double beta;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 4 && argc != 8) {
	interp->result = "WARNING integrator Newmark gamma beta <alphaM> <betaKcurrent> <betaKi> <betaKlastCommitted>";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &gamma) != TCL_OK) {
	  interp->result = "WARNING integrator Newmark gamma beta - undefined gamma";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
	  interp->result = "WARNING integrator Newmark gamma beta - undefined beta";
	  return TCL_ERROR;	
      }
      if (argc == 8) {
	  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
	      cerr << "WARNING integrator Newmark gamma beta alphaM betaK betaKi betaKc - alphaM";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
	      cerr << "WARNING integrator Newmark gamma beta alphaM betaK betaKi betaKc - betaK";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
	      cerr << "WARNING integrator Newmark gamma beta alphaM betaK betaKi betaKc - betaKi";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
	      cerr << "WARNING integrator Newmark gamma beta alphaM betaK betaKi betaKc - betaKc";
	      return TCL_ERROR;	
	  }
      }
      if (argc == 4)
	  theTransientIntegrator = new Newmark(gamma,beta);       
      else
	  theTransientIntegrator = new Newmark(gamma,beta,alphaM,betaK,betaKi,betaKc);
  }  
  
  else if (strcmp(argv[1],"HHT") == 0) {
      double alpha;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 7) {
	interp->result = "WARNING integrator HHT alpha <alphaM> <betaKcurrent> <betaKi> <betaKlastCommitted>";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &alpha) != TCL_OK) {
	  interp->result = "WARNING integrator HHT alpha - undefined alpha";	  
	  return TCL_ERROR;	
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      cerr << "WARNING integrator HHT alpha alphaM betaK betaKi betaKc - alphaM";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      cerr << "WARNING integrator HHT alpha alphaM betaK betaKi betaKc - betaK";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      cerr << "WARNING integrator HHT alpha alphaM betaK betaKi betaKc - betaKi";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	      cerr << "WARNING integrator HHT alpha alphaM betaK betaKi betaKc - betaKc";
	      return TCL_ERROR;	
	  }
      }      
      if (argc == 3)
	  theTransientIntegrator = new HHT(alpha);       
      else
	  theTransientIntegrator = new HHT(alpha,alphaM,betaK, betaKi, betaKc);       
  }    


  else if (strcmp(argv[1],"Newmark1") == 0) {
      double gamma;
      double beta;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 4 && argc != 8) {
	interp->result = "WARNING integrator Newmark1 gamma beta <alphaM> <betaKcurrent> <betaKi> <betaKlastCommitted>";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &gamma) != TCL_OK) {
	  interp->result = "WARNING integrator Newmark1 gamma beta - undefined gamma";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
	  interp->result = "WARNING integrator Newmark1 gamma beta - undefined beta";
	  return TCL_ERROR;	
      }
      if (argc == 8) {
	  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
	      cerr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi betaKc - alphaM";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
	      cerr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi betaKc - betaK";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
	      cerr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi betaKc - betaKi";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
	      cerr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi betaKc - betaKc";
	      return TCL_ERROR;	
	  }
      }
      if (argc == 4)
	  theTransientIntegrator = new Newmark1(gamma,beta);       
      else
	  theTransientIntegrator = new Newmark1(gamma,beta,alphaM,betaK,betaKi,betaKc);
  }

  else if (strcmp(argv[1],"HHT1") == 0) {
      double alpha;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 7) {
	interp->result = "WARNING integrator HHT1 alpha <alphaM> <betaKcurrent> <betaKi> <betaKlastCommitted>";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &alpha) != TCL_OK) {
	  interp->result = "WARNING integrator HHT alpha - undefined alpha";	  
	  return TCL_ERROR;	
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      cerr << "WARNING integrator HHT1 gamma beta alphaM betaK betaKi betaKc - alphaM";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      cerr << "WARNING integrator HHT1 gamma beta alphaM betaK betaKi betaKc - betaK";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      cerr << "WARNING integrator HHT1 gamma beta alphaM betaK betaKi betaKc - betaKi";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	      cerr << "WARNING integrator HHT1 gamma beta alphaM betaK betaKi betaKc - betaKc";
	      return TCL_ERROR;	
	  }
      }      
      if (argc == 3)
	  theTransientIntegrator = new HHT1(alpha);       
      else
	  theTransientIntegrator = new HHT1(alpha,alphaM,betaK, betaKi, betaKc);       
  }    

  
  else if (strcmp(argv[1],"WilsonTheta") == 0) {
      double theta, alphaM,betaK;
      if (argc != 3 && argc != 5) {
	interp->result = "WARNING integrator WilsonTheta theta <alphaM> <betaK>";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &theta) != TCL_OK) {
	  interp->result = "WARNING integrator WilsonTheta theta - undefined theta";	  
	  return TCL_ERROR;	
      }
      if (argc == 5) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      cerr << "WARNING integrator WilsonTheta gamma beta alphaM betaK - alphaM";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      cerr << "WARNING integrator WilsonTheta gamma beta alphaM betaK  - betaK";
	      return TCL_ERROR;	
	  }
      }            
      if (argc == 3)
	  theTransientIntegrator = new WilsonTheta(theta);       
      else
	  theTransientIntegrator = new WilsonTheta(theta,alphaM,betaK);       	  
  }      

  else {
    interp->result = "WARNING No Integrator type exists ";
    return TCL_ERROR;
  }    

  // if the analysis exists - we want to change the Integrator
  if (theStaticAnalysis != 0)
    theStaticAnalysis->setIntegrator(*theStaticIntegrator);
  else if (theTransientAnalysis != 0)
    theTransientAnalysis->setIntegrator(*theTransientIntegrator);

  return TCL_OK;
}


extern int
TclAddRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
	       char **argv, Domain &theDomain);

int 
addRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
	    char **argv)
{
    return TclAddRecorder(clientData, interp, argc, argv, theDomain);
}

extern int
TclAddAlgorithmRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
			char **argv, Domain &theDomain, EquiSolnAlgo *theAlgorithm);

int 
addAlgoRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
	    char **argv)
{
	if (theAlgorithm != 0)
		return TclAddAlgorithmRecorder(clientData, interp, argc, argv,
			theDomain, theAlgorithm);

	else
		return 0;
}

extern int
TclAddDatabase(ClientData clientData, Tcl_Interp *interp, int argc, 
	       char **argv, Domain &theDomain);

int 
addDatabase(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)
{
  return TclAddDatabase(clientData, interp, argc, argv, theDomain);
}


int 
playbackRecorders(ClientData clientData, Tcl_Interp *interp, int argc, 
		  char **argv)
{
  // make sure at least one other argument to contain integrator
  if (argc < 2) {
      interp->result = "WARNING need to specify the commitTag ";
      return TCL_ERROR;
  }    

  // check argv[1] for type of Numberer and create the object
  int cTag;
  if (Tcl_GetInt(interp, argv[1], &cTag) != TCL_OK)	
      return TCL_ERROR;	

  theDomain.playback(cTag);
  return TCL_OK;

}

/*
int 
groundExcitation(ClientData clientData, Tcl_Interp *interp, int argc, 
		  char **argv)
{
  // make sure at least one other argument to contain integrator
  if (argc < 2) {
      interp->result = "WARNING need to specify the commitTag ";
      return TCL_ERROR;
  }    

  if (strcmp(argv[1],"Single") == 0) {
      if (argc < 4) {
	interp->result = "WARNING quake single dof motion";
	return TCL_ERROR;
      }    

      int dof;
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK)	
	  return TCL_ERROR;	      
      
      // read in the ground motion
      GroundMotion *theMotion;
      if (strcmp(argv[3],"ElCentro") == 0) {
	  double fact = 1.0;
	  if (argc == 5) {
	      if (Tcl_GetDouble(interp, argv[4], &fact) != TCL_OK)	
		  return TCL_ERROR;	
	  }
	  theMotion = new ElCentroGroundMotion(fact);
      } else {
	  interp->result = "WARNING quake Single motion - no motion type exists ";
	  return TCL_ERROR;      
      }

      Load *theLoad = new SingleExcitation(*theMotion, dof, nextTag++);
      theDomain.addOtherLoad(theLoad);
      return TCL_OK;
  }  
  
  else {
    interp->result = "WARNING No quake type exists ";
    return TCL_ERROR;
  }    
}
*/

int 
rigidLink(ClientData clientData, Tcl_Interp *interp, int argc, 
	  char **argv)
{
  if (argc < 4) {
      interp->result = "WARNING rigidLink linkType? rNode? cNode?\n";
      return TCL_ERROR;
  }    

  int numMPs = theDomain.getNumMPs();
  int rNode, cNode;
  if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
      cerr << "WARNING rigidLink linkType? rNode? cNode? - could not read rNode ";
      return TCL_ERROR;	        
  }
  if (Tcl_GetInt(interp, argv[3], &cNode) != TCL_OK) {
      cerr << "WARNING rigidLink linkType? rNode? cNode? - could not read CNode ";
      return TCL_ERROR;	        
  }

  // construct a rigid rod or beam depending on 1st arg
  if ((strcmp(argv[1],"-bar") == 0) || (strcmp(argv[1],"bar") == 0)) {
    RigidRod theLink(theDomain, rNode, cNode, numMPs);
  } else if ((strcmp(argv[1],"-beam") == 0) || (strcmp(argv[1],"beam") == 0)) {
    RigidBeam theLink(theDomain, rNode, cNode, numMPs);
  } else {
      cerr << "WARNING rigidLink linkType? rNode? cNode? - unrecognised link type (-bar, -beam) ";
      return TCL_ERROR;	        
  }

  return TCL_OK;
}



int 
rigidDiaphragm(ClientData clientData, Tcl_Interp *interp, int argc, 
	   char **argv)
{
  if (argc < 3) {
      interp->result = "WARNING rigidLink perpDirn? rNode? <cNodes?>";
      return TCL_ERROR;
  }    

  int rNode, perpDirn;
  if (Tcl_GetInt(interp, argv[1], &perpDirn) != TCL_OK) {
      cerr << "WARNING rigidLink perpDirn rNode cNodes - could not read perpDirn? ";
      return TCL_ERROR;	        
  }

  if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
      cerr << "WARNING rigidLink perpDirn rNode cNodes - could not read rNode ";
      return TCL_ERROR;	        
  }
  
  // read in the constrained Nodes
  int numConstrainedNodes = argc - 3;
  ID constrainedNodes(numConstrainedNodes);
  for (int i=0; i<numConstrainedNodes; i++) {
      int cNode;
      if (Tcl_GetInt(interp, argv[3+i], &cNode) != TCL_OK) {
	  cerr << "WARNING rigidLink perpDirn rNode cNodes - could not read a cNode";
	  return TCL_ERROR;	        
      }
      constrainedNodes(i) = cNode;
  }
  int numMPs = theDomain.getNumMPs();
  RigidDiaphragm theLink(theDomain, rNode, constrainedNodes, 
	perpDirn-1, numMPs);

  return TCL_OK;
}


int 
eigenAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, 
		char **argv)
{
     // make sure at least one other argument to contain type of system
    if (argc < 2) {
	interp->result = "WARNING want - eigen #eigenValues?";
	return TCL_ERROR;
    }    

    // check argv[1] for type of SOE and create it
    int numEigen;
    if ((Tcl_GetInt(interp, argv[1], &numEigen) != TCL_OK) || numEigen < 0) {
	interp->result = "WARNING analysis Static numIncr  - illegal numIncr";    
	return TCL_ERROR;	
    }	

    if (theEigenAnalysis == 0) {
      EigenAlgorithm *theEigenAlgo = new FrequencyAlgo();
      EigenIntegrator  *theEigenIntegrator = new EigenIntegrator();    
      AnalysisModel *theEigenModel = new AnalysisModel();

      /*
      SymArpackSOE *theEigenSOE;
      SymArpackSolver *theEigenSolver;
      theEigenSolver = new SymArpackSolver(numEigen); 
      theEigenSOE = new SymArpackSOE(*theEigenSolver, *theEigenModel);    
      */

      BandArpackSOE *theEigenSOE;
      BandArpackSolver *theEigenSolver;    
      theEigenSolver = new BandArpackSolver(numEigen); 
      theEigenSOE = new BandArpackSOE(*theEigenSolver, *theEigenModel);    

      RCM *theRCM = new RCM();	
      DOF_Numberer *theEigenNumberer = new DOF_Numberer(*theRCM);    	

      // ConstraintHandler *theEigenHandler = new PlainHandler();
      ConstraintHandler *theEigenHandler = new TransformationConstraintHandler();
      theEigenAnalysis = new EigenAnalysis(theDomain,
					   *theEigenHandler,
					   *theEigenNumberer,
      					   *theEigenModel,
					   *theEigenAlgo,
					   *theEigenSOE,
					   *theEigenIntegrator);

    }    
    
    if (theEigenAnalysis->analyze(numEigen) == 0) {
      const Vector &eigenvalues = theDomain.getEigenvalues();
      cerr << "EigenValues: " << eigenvalues;
    }

    delete theEigenAnalysis;
    theEigenAnalysis = 0;

    return TCL_OK;

}



int 
videoPlayer(ClientData clientData, Tcl_Interp *interp, int argc, 
	    char **argv)
{
    
    // make sure at least one other argument to contain type of system
    if (argc < 5) {
	interp->result = "WARNING want - video -window windowTitle? -file fileName?\n";
	return TCL_ERROR;
    }    

    char *wTitle =0;
    char *fName = 0;
	char *imageName = 0;
    int endMarker = 1;
    while (endMarker < (argc-1)) {
      if (strcmp(argv[endMarker],"-window") == 0) {
	wTitle = argv[endMarker+1];
	endMarker+=2;
      } else if (strcmp(argv[endMarker],"-file") == 0) {
	fName = argv[endMarker+1];
	endMarker+=2;
      } else if (strcmp(argv[endMarker],"-image") == 0) {
	imageName = argv[endMarker+1];
	endMarker += 2;
	  }
      else {
	g3ErrorHandler->warning("WARNING unknown %s want - video -window windowTitle? -file fileName?\n", 
				argv[endMarker]);
	return TCL_ERROR;
      }
    }
    
    if (wTitle != 0 && fName != 0) {
      // delete the old video player if one exists
      if (theTclVideoPlayer != 0)
	delete theTclVideoPlayer;

      // create a new player
      theTclVideoPlayer = new TclVideoPlayer(wTitle, fName, imageName, interp);
    }
    else
      return TCL_ERROR;

    return TCL_OK;
}



int 
removeObject(ClientData clientData, Tcl_Interp *interp, int argc, 
	     char **argv)
{
    
    // make sure at least one other argument to contain type of system
    if (argc < 3) {
	interp->result = "WARNING want - remove objectType? objectTag?\n";
	return TCL_ERROR;
   }    

    int tag;
    if (strcmp(argv[1],"element") == 0) {
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	cerr << "WARNING remove element tag? failed to read tag: " << argv[2] << endl;
	return TCL_ERROR;
      }      
      Element *theEle = theDomain.removeElement(tag);
      if (theEle != 0) {
	// we also have to remove any elemental loads from the domain
	LoadPatternIter &theLoadPatterns = theDomain.getLoadPatterns();
	LoadPattern *thePattern;
	
	// go through all load patterns
	while ((thePattern = theLoadPatterns()) != 0) {
	  ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
	  ElementalLoad *theLoad;

	  // go through all elemental loads in the pattern
	  while ((theLoad = theEleLoads()) != 0) {

	    // remove & destroy elemental load if tag corresponds to element being removed
	    if (theLoad->getElementTag() == tag) {
	      thePattern->removeElementalLoad(theLoad->getTag());
	      delete theLoad;
	    }
	  }
	}

	// finally invoke the destructor on the element
	delete theEle;
      }
    }

    else
      cerr << "WARNING remove element tag? - only command available at the moment: " << endl;

    return TCL_OK;
}





