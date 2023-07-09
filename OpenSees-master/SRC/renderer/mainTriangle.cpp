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
// $Date: 2003-02-14 23:01:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/mainTriangle.cpp,v $
                                                                        
                                                                        
// File: ~/model/main.C
//
// Written: fmk 12/95
// Revised:
//
// Purpose: This file is a driver to create a 2d plane-frame
// model and to perform a linear static analysis on that model.
// 
//

#include <X11Renderer.h>
#include <DofColorMap.h>

#define numP 8

#define BLAS_EXISTS 1
// #define PETSC_EXISTS 1

#include <stdlib.h>
#include <iOPS_Stream.h>

#include <Timer.h>
#include <TrianglePlaneStress.h>

#include <Node.h>
#include <Domain.h>

#include <AnalysisModel.h>

#include <Linear.h>
#include <NewtonRaphson.h>

#include <PlainHandler.h>
#include <PenaltyConstraintHandler.h>

#include <PlainNumberer.h>
#include <DOF_Numberer.h>

#include <StaticIntegrator.h>
#include <LoadControl.h>

#include <CTestNormUnbalance.h>
#include <CTestNormDispIncr.h>
#include <CTestEnergyIncr.h>

#include <StaticAnalysis.h>

#include <SlowLinearSOE.h>
#include <SlowLinearSOESolver.h>
#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>
#include <ProfileSPDLinDirectBlockSolver.h>

#ifdef BLAS_EXISTS
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>
#include <FullGenLinSOE.h>
#include <FullGenLinLapackSolver.h>
#include <ProfileSPDLinDirectThreadSolver.h>
#include <ProfileSPDLinDirectSkypackSolver.h>
#include <SparseGenLinSOE.h>
#include <SparseGenSuperLuSolver.h>
#include <SparseGenSuperLuThreadSolver.h>
#endif

#ifdef PETSC_EXISTS
#include <PetscSOE.h>
#include <PetscSolver.h>
#include <mpi.h>
#endif

#include <Graph.h>
#include <RCM.h>
#include <Metis.h>

int main(int argc, char **argv)
{
    //
    // read in the number of analysis iterations
    //    - the first iteartion will be expensive due to
    //      construction of FE_Element, DOF_Groups, node
    //      and dof graphs
    // 

    int numReps;
    opserr << "Enter Number of Iterations: ";
    cin >> numReps;
    
#ifdef PETSC_EXISTS
    int numP; // number of processes started
    int me;  // the integer id of the process .. 0 through numP-1
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &numP);
    static char help[] = "My daft little program\n";
    PetscInitialize(&argc, &argv, (char *)0, help);
#endif	
    
    //
    //	now create a domain and a modelbuilder
    //  and build the model
    //

    Domain *theDomain = new Domain();
    TrianglePlaneStress  *theModelBuilder 
      = new TrianglePlaneStress(*theDomain);
    theModelBuilder->buildFE_Model();


    theModelBuilder->fixNodeXdirnXloc(100.0,0.0);
    theModelBuilder->fixNodeYdirnYloc(100.0,0.0);

    theModelBuilder->addPointLoad(0.0,100.0, -50.0, 0.0);

    //
    // create an Analysis object to perform a static analysis of the model
    //  - constructs:
    //    AnalysisModel of type AnalysisModel,
    //	  EquiSolnAlgo of type Linear
    //	  StaticIntegrator of type StaticIntegrator
    //	  ConstraintHandler of type PlainHandler (only homo SP constraints)
    //    DOF_Numberer which does either RCM or takes order provided
    //    LinearSOE of types full and band general, band and profile SPD

    AnalysisModel 	*theModel = new AnalysisModel();
    EquiSolnAlgo 	*theSolnAlgo = new Linear();
    ConvergenceTest *theTest = new CTestNormUnbalance(1.0e-8);
    // EquiSolnAlgo 	*theSolnAlgo = new NewtonRaphson(10,*theTest); 
    StaticIntegrator	*theIntegrator = new LoadControl(1);

    ConstraintHandler 	*theHandler;
    cout << "Enter ConstraintHandler type: [1]-Plain, [2]-Penalty: ";
    int handlerType;
    cin >> handlerType;
    if (handlerType == 1)
	theHandler = new PlainHandler();
    else
	theHandler = new PenaltyConstraintHandler(1.0e16,1.0e16);
    
    // construct the DOF_Numberer
    cout << "[1]-RCM, [2]-Metis : ";
    int numType;
    cin >> numType;    
    
    DOF_Numberer *theNumberer;
    GraphNumberer *theGN;
    if (numType == 1) {
      theGN = new RCM;
      theNumberer = new DOF_Numberer(*theGN);       
    } else {
      theGN = new Metis(numP);
      theNumberer = new DOF_Numberer(*theGN);       
    }

    // construct the LinearSOE
    LinearSOE *theSOE;

#ifdef PETSC_EXISTS        
    cout << "Enter SystemOfEquation type: [1]-FullGen, [2]-BandGen, ";
    cout << "[3]-BandSPD, [4]-ProfileSPD, [5]-SparseGEN, [6]-Petsc : ";
    int soeType;
    cin >> soeType;    
    if (soeType == 1) {
      FullGenLinSolver    *theSolver = new FullGenLinLapackSolver();
      theSOE = new FullGenLinSOE(*theSolver);            
    } else if (soeType == 2) {
      BandGenLinSolver    *theSolver = new BandGenLinLapackSolver();
      theSOE = new BandGenLinSOE(*theSolver);      
    } else if (soeType == 3) {
      BandSPDLinSolver    *theSolver = new BandSPDLinLapackSolver();   
      theSOE = new BandSPDLinSOE(*theSolver);        
    } else if (soeType == 4) { // a ProfileSPD
      // can choose from a variety of solver for a profile SPD
      ProfileSPDLinSolver *theSolver;
      cout << "Enter Solver Type: [1]-Normal, [2]-Block, [3]-Skypack: ";
      int solverType;
      cin >> solverType;
      if (solverType == 1) 
	theSolver = new ProfileSPDLinDirectSolver(); 
      else if (solverType == 2) {
	cout << "Enter blockSize: ";
	int blockSize;
	cin >> blockSize;
	theSolver = new ProfileSPDLinDirectBlockSolver(1.0e-12,blockSize); 
      } else {
	cout << "Enter mRows, mCols [e.g. 64 128] ";
	int mRows, mCols;
	cin >> mRows >> mCols;
	theSolver = new ProfileSPDLinDirectSkypackSolver(mCols, mRows); 
      }
      theSOE = new ProfileSPDLinSOE(*theSolver);      
    }
    else if (soeType == 5) {
      // cout << "Enter Num Threads: ";
      //      int NP;
      //      cin >> NP;
      cout << "Enter Ordering [0]: natural, [1]-min Degree A^tA ";
      cout << "[2]-min Degree A^t+A: ";
      int perm;
      cin >> perm;
      SparseGenLinSolver  *
	theSolver = new SparseGenSuperLuSolver(perm,3,3,0.0);   
      //      SparseGenLinSolver  
      //	*theSolver = new SparseGenSuperLuThreadSolver(NP,perm,3,3,0.0);   
      theSOE = new SparseGenLinSOE(*theSolver);        
    } 
    else {
      PetscSolver *theSolver = new PetscSolver(KSPCG, PCJACOBI);
      theSOE = new PetscSOE(*theSolver);

  }
#endif
    
#ifdef BLAS_EXISTS        
    cout << "Enter SystemOfEquation type: [1]-FullGen, [2]-BandGen, ";
    cout << "[3]-BandSPD, [4]-ProfileSPD, [5]-SparseGEN : ";
    int soeType;
    cin >> soeType;    
    if (soeType == 1) {
      FullGenLinSolver    *theSolver = new FullGenLinLapackSolver();
      theSOE = new FullGenLinSOE(*theSolver);            
    } else if (soeType == 2) {
      BandGenLinSolver    *theSolver = new BandGenLinLapackSolver();
      theSOE = new BandGenLinSOE(*theSolver);      
    } else if (soeType == 3) {
      BandSPDLinSolver    *theSolver = new BandSPDLinLapackSolver();   
      theSOE = new BandSPDLinSOE(*theSolver);        
    } else if (soeType == 4) { // a ProfileSPD
      // can choose from a variety of solver for a profile SPD
      ProfileSPDLinSolver *theSolver;
      cout << "Enter Solver Type: [1]-Normal, [2]-Block, [3]-Skypack: ";
      int solverType;
      cin >> solverType;
      if (solverType == 1) 
	theSolver = new ProfileSPDLinDirectSolver(); 
      else if (solverType == 2) {
	cout << "Enter blockSize: ";
	int blockSize;
	cin >> blockSize;
	theSolver = new ProfileSPDLinDirectBlockSolver(1.0e-12,blockSize); 
      } else {
	cout << "Enter mRows, mCols [e.g. 64 128] ";
	int mRows, mCols;
	cin >> mRows >> mCols;
	theSolver = new ProfileSPDLinDirectSkypackSolver(mCols, mRows); 
      }
      theSOE = new ProfileSPDLinSOE(*theSolver);      
    }
    else {
      // cout << "Enter Num Threads: ";
      //      int NP;
      //      cin >> NP;
      cout << "Enter Ordering [0]: natural, [1]-min Degree A^tA ";
      cout << "[2]-min Degree A^t+A: ";
      int perm;
      cin >> perm;
      SparseGenLinSolver  *
	theSolver = new SparseGenSuperLuSolver(perm,3,3,0.0);   
      //      SparseGenLinSolver  
      //	*theSolver = new SparseGenSuperLuThreadSolver(NP,perm,3,3,0.0);   
      theSOE = new SparseGenLinSOE(*theSolver);        
  }
#endif    
    


    StaticAnalysis 	theAnalysis(*theDomain,
				    *theHandler,
				    *theNumberer,
				    *theModel,
				    *theSolnAlgo,
				    *theSOE,
				    *theIntegrator);
    
    //
    // start the timing
    // 

    Timer timer;
    Timer timer1;

    timer.start();
    timer1.start();


    for (int j=0; j<numReps; j++) {
	timer.start();
	theAnalysis.analyze();    
	timer.pause();
	theSolnAlgo->Print(opserr);
        timer1.pause();
	cout << "Main::Time To Analyse Rep: " << 
	  j+1 << "  real:  " << timer.getReal();
        cout << "  cpu: " << timer.getCPU() << 
	  "  page: " << timer.getNumPageFaults() << endln; 
	cout << "Main::Time Since Start: " 
	     << j+1 << "  real:  " << timer1.getReal();
        cout << "  cpu: " << timer1.getCPU() 
	  << "  page: " << timer1.getNumPageFaults() << endln; 

    }

    // done
    //    opserr << *theDomain;
    
    Node *theNode = theDomain->getNode(2);
    if (theNode != 0) {
	opserr << *theNode;
    }

    //const Vector &theResult = theSOE->getX();
    //opserr << "x(0): " << theResult(0) << endln;

    DofColorMap theMap(theSOE->getNumEqn(),numP);
    WindowRenderer *theViewer = new X11Renderer(600,600,*theDomain, theMap);
    theViewer->setPRP(50, 50, 100);
    theViewer->setVRP(50, 50, 0);
    theViewer->setVUP(0,1,0);
    theViewer->setVPN(0,0,1);
    theViewer->setViewWindow(-5,105,-5, 105);
    theViewer->setPlaneDist(0,150);
    theViewer->setPortWindow(-1,1,-1, 1);
    theViewer->setProjectionMode(0); // 0 parallel mode, 1 perspective
    theViewer->setFillMode(1); // 0 fill, 1 wire
    
    theViewer->displayModel(0,0);

    int data;
    cin >> data;

#ifdef PETSC_EXISTS    
    MPI_Finalize();    
    PetscFinalize();
#endif
    exit(0);
}	
	
