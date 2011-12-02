///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              EightNodeBrick.cpp
// CLASS:             EightNodeBrick
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic, Zhaohui Yang and Xiaoyan Wu
// PROGRAMMER:        Boris Jeremic, Zhaohui Yang  and Xiaoyan Wu
// DATE:              Aug. 2000
// UPDATE HISTORY:			 Modified from Brick3D and FourNodeQuad.hh  07/06/00
//																			 Sept. - Oct 2000 connected to OpenSees by Zhaohui
//
//  Purpose: This file contains a main procedure which is used to test the       !
//           EightNodeBrick element with ElasticIsotropic3D 			 !
//
// CONTACT:           jeremic@ucdavis.edu
///////////////////////////////////////////////////////////////////////////////


// standard C++ includes
#include <stdlib.h>
#include <iostream.h>

#include <G3Globals.h>
#include <ConsoleErrorHandler.h>

#include <ArrayOfTaggedObjects.h>
#include <MapOfTaggedObjects.h>
#include <MapOfTaggedObjectsIter.h>

// includes for the domain classes
#include <Domain.h>
#include <Node.h>
//#include <Truss.h>
#include <EightNodeBrick.h>

//#include <ElasticMaterial.h>
#include <ElasticIsotropic3D.h>

#include <SP_Constraint.h>
#include <LoadPattern.h>
#include <LinearSeries.h>
#include <NodalLoad.h>

// includes for the analysis classes
#include <StaticAnalysis.h>
#include <AnalysisModel.h>
#include <Linear.h>
#include <PenaltyConstraintHandler.h>
#include <DOF_Numberer.h>
#include <RCM.h>
#include <LoadControl.h>
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>

ErrorHandler *g3ErrorHandler;

int main(int argc, char **argv)
{
    //  first build our error handler
    g3ErrorHandler = new ConsoleErrorHandler();

    //
    //	now create a domain and a modelbuilder
    //  and build the model
    //

    //  Domain *theDomain = new Domain();

    MapOfTaggedObjects theStorage;
    //  ArrayOfTaggedObjects theStorage(32);    
    Domain *theDomain = new Domain(theStorage);
    
    // create the nodes using constructor: 
    //		Node(tag, ndof, crd1, crd2, crd3 )
    // and then add them to the domain
    
    Node *node1 = new Node(1, 3, 1.0, 1.0, 0.0 );
    Node *node2 = new Node(2, 3, 0.0, 1.0, 0.0 );
    Node *node3 = new Node(3, 3, 0.0, 0.0, 0.0 );    
    Node *node4 = new Node(4, 3, 1.0, 0.0, 0.0 );        
    Node *node5 = new Node(5, 3, 1.0, 1.0, 1.0 );
    Node *node6 = new Node(6, 3, 0.0, 1.0, 1.0 );
    Node *node7 = new Node(7, 3, 0.0, 0.0, 1.0 );    
    Node *node8 = new Node(8, 3, 1.0, 0.0, 1.0 );        
   
    theDomain->addNode(node1);
    theDomain->addNode(node2);
    theDomain->addNode(node3);
    theDomain->addNode(node4);
    theDomain->addNode(node5);
    theDomain->addNode(node6);
    theDomain->addNode(node7);
    theDomain->addNode(node8);
    
    // create an elastic isotropic 3D material using constriuctor:  
    //		ElasticIsotropic3D(tag, E, v)

    NDMaterial *theMaterial = new ElasticIsotropic3D(1, 3000, 0.3);
    
    // create the EightNodeBrick elements using constructor:
    //		EightNodeBrick(tag, nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd8, 
    //                           Material &, const *char type, double bodyforce1, double bodyforce2,
    //			         double pressure, double rho, EPState *InitEPS)
    // and then add them to the domain
        
    // For elastic material, no EPState, just supply null.
    EPState *eps = 0; 
    EightNodeBrick *brick = new EightNodeBrick(1, 5, 6, 7, 8, 1, 2, 3, 4,
                              theMaterial, "ElasticIsotropic3D", 
			      0.0, 0.0, 0.0, 1.8, eps);        
    
    theDomain->addElement( brick );    
    
    // create the single-point constraint objects using constructor:
    //		SP_Constraint(tag, nodeTag, dofID, value)
    // and then add them to the domain
    
    SP_Constraint *sp1  = new SP_Constraint(1,  1, 0, 0.0);
    SP_Constraint *sp2  = new SP_Constraint(2,  1, 1, 0.0);    
    SP_Constraint *sp3  = new SP_Constraint(3,  1, 2, 0.0);
    SP_Constraint *sp4  = new SP_Constraint(4,  2, 0, 0.0);    
    SP_Constraint *sp5  = new SP_Constraint(5,  2, 1, 0.0);
    SP_Constraint *sp6  = new SP_Constraint(6,  2, 2, 0.0);        
    SP_Constraint *sp7  = new SP_Constraint(7,  3, 0, 0.0);
    SP_Constraint *sp8  = new SP_Constraint(8,  3, 1, 0.0);    
    SP_Constraint *sp9  = new SP_Constraint(9,  3, 2, 0.0);
    SP_Constraint *sp10 = new SP_Constraint(10, 4, 0, 0.0);    
    SP_Constraint *sp11 = new SP_Constraint(11, 4, 1, 0.0);
    SP_Constraint *sp12 = new SP_Constraint(12, 4, 2, 0.0);        

    theDomain->addSP_Constraint(sp1 );
    theDomain->addSP_Constraint(sp2 );
    theDomain->addSP_Constraint(sp3 );
    theDomain->addSP_Constraint(sp4 );    
    theDomain->addSP_Constraint(sp5 );    
    theDomain->addSP_Constraint(sp6 );    
    theDomain->addSP_Constraint(sp7 );
    theDomain->addSP_Constraint(sp8 );
    theDomain->addSP_Constraint(sp9 );
    theDomain->addSP_Constraint(sp10);    
    theDomain->addSP_Constraint(sp11);    
    theDomain->addSP_Constraint(sp12);    

    // construct a linear time series object using constructor:
    //		LinearSeries()
    
    TimeSeries *theSeries = new LinearSeries();
    
    // construct a load pattren using constructor:
    //		LoadPattern(tag)
    // and then set it's TimeSeries and add it to the domain
    
    LoadPattern *theLoadPattern = new LoadPattern(1);
    theLoadPattern->setTimeSeries(theSeries);
    theDomain->addLoadPattern(theLoadPattern);
    
    // construct a nodal load using constructor:
    //		NodalLoad(tag, nodeID, Vector &)
    // first construct a Vector of size 3 and set the values NOTE C INDEXING
    // then construct the load and add it to the domain
    
    Vector theLoadValues(3);
    theLoadValues(0) = 0.0;
    theLoadValues(1) = 0.0;
    theLoadValues(2) = -100.0;
    
    NodalLoad *theLoad = new NodalLoad(1, 5, theLoadValues);
    theDomain->addNodalLoad(theLoad, 1);
    
    theLoad = new NodalLoad(2, 6, theLoadValues);
    theDomain->addNodalLoad(theLoad, 1);
    
    theLoad = new NodalLoad(3, 7, theLoadValues);
    theDomain->addNodalLoad(theLoad, 1);
    
    theLoad = new NodalLoad(4, 8, theLoadValues);
    theDomain->addNodalLoad(theLoad, 1);

    // create an Analysis object to perform a static analysis of the model
    //  - constructs:
    //    AnalysisModel of type AnalysisModel,
    //	  EquiSolnAlgo of type Linear
    //	  StaticIntegrator of type LoadControl
    //	  ConstraintHandler of type Penalty
    //    DOF_Numberer which uses RCM
    //    LinearSOE of type Band SPD
    // and then the StaticAnalysis object
    
    AnalysisModel     *theModel = new AnalysisModel();
    EquiSolnAlgo      *theSolnAlgo = new Linear();
    StaticIntegrator  *theIntegrator = new LoadControl(1.0, 1.0, 1.0, 1.0);
    ConstraintHandler *theHandler = new PenaltyConstraintHandler(1.0e8,1.0e8);
    RCM               *theRCM = new RCM();
    DOF_Numberer      *theNumberer = new DOF_Numberer(*theRCM);    
    BandSPDLinSolver  *theSolver = new BandSPDLinLapackSolver();       
    LinearSOE         *theSOE = new BandSPDLinSOE(*theSolver);        

    StaticAnalysis    theAnalysis(*theDomain,
				  *theHandler,
				  *theNumberer,
				  *theModel,
				  *theSolnAlgo,
				  *theSOE,
				  *theIntegrator);

    // perform the analysis & print out the results for the domain
    int numSteps = 1;
    theAnalysis.analyze(numSteps);
    cerr << *theDomain;

    //brick->displaySelf();

    exit(0);
}	
	
