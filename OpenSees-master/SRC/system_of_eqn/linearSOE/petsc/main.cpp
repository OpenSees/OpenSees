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
// $Date: 2006-01-12 23:59:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/petsc/main.cpp,v $
                                                                        
                                                                        

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
//#include <PetscSOE.h>
//#include <PetscSolver.h>
#include <PetscSOE.h>
#include <PetscSolver.h>
#include <Graph.h>
#include <Vertex.h>
#include <SuperLU.h>
#include <SparseGenColLinSOE.h>
#include <PetscSparseSeqSolver.h>
#include <SparseGenRowLinSOE.h>

// init the global variabled defined in OPS_Globals.h
#include <OPS_Globals.h>
#include <StandardStream.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;
 
double        ops_Dt = 0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;

main() 
{
    Matrix a(2,2);
    a(0,0) = 2.0;
    a(0,1) = 1.0;
    a(1,0) = 1.0;
    a(1,1) = 2.0;
    
    Matrix b(1,1);
    b(0,0) = 10;
    
    Matrix c(3,3);
    c(0,0) = 9.0;
    c(0,1) = 2.0;
    c(0,2) = 2.0;
    c(2,0) = 2.0;
    c(2,1) = 2.0;
    c(2,2) = 9.0;
    c(1,0) = 2.0;
    c(1,1) = 9.0;
    c(1,2) = 2.0;		      
    
    ID d(1);    d(0) =  1;
    ID e(1);    e(0) =  2;
    ID f(2);    f(0) =  1; f(1) = 0;
    ID g(2);    g(0) = -1; g(1) = 2;
    ID h(3);    h(0) =  2; h(1) = 1; h(2) = 3;
    ID m(3);    m(0) =  5; m(1) = 4; m(2) = 3;    

    Graph theGraph(6);  
    Vertex v0(0,0);  theGraph.addVertex(&v0);
    Vertex v1(1,1);  theGraph.addVertex(&v1);
    Vertex v2(2,2);  theGraph.addVertex(&v2);
    Vertex v3(3,3);  theGraph.addVertex(&v3);
    Vertex v4(4,4);  theGraph.addVertex(&v4);
    Vertex v5(5,5);  theGraph.addVertex(&v5);
    theGraph.addEdge(0,1);
    theGraph.addEdge(2,1);
    theGraph.addEdge(1,3);
    theGraph.addEdge(5,4);
    theGraph.addEdge(4,3);

    Vector x(3); x(0) = 2; x(1) = 3; x(2) = 4.5;

    SuperLU *theSolver1 = new SuperLU;
    SparseGenColLinSOE *theSOE1 = new SparseGenColLinSOE(*theSolver1);

    opserr << "Test One with SuperLU\n";

    theSOE1->setSize(theGraph);
    theSOE1->addA(a,f);
    theSOE1->addA(a,g);
    theSOE1->addA(c,h);
    theSOE1->addA(c,m,2.0);    

    theSOE1->addB(x,h);
    theSOE1->addB(x,m,-1);

    theSOE1->solve();
    opserr << "B: " << theSOE1->getB();
    opserr << "X: " << theSOE1->getX();

    delete theSOE1;
    
    KSPType method = KSPCG;            // KSPCG KSPGMRES
    PCType preconditioner = PCJACOBI; // PCJACOBI PCILU PCBJACOBI
    
    PetscSparseSeqSolver *theSolver2 = new     PetscSparseSeqSolver(method, preconditioner);
    SparseGenRowLinSOE *theSOE2 = new SparseGenRowLinSOE(*theSolver2);

    opserr << "Test Two with Petsc\n";

    theSOE2->setSize(theGraph);

    theSOE2->addA(a,f);
    theSOE2->addA(a,g);
    theSOE2->addA(c,h);
    theSOE2->addA(c,m,2.0);    

    theSOE2->addB(x,h);
    theSOE2->addB(x,m,-1);

    theSOE2->solve();
    opserr << "B: " << theSOE2->getB();
    opserr << "X: " << theSOE2->getX();

    // solve again
    theSOE2->solve();
    opserr << "  results on a second solve with same B\n";
    opserr << "B: " << theSOE2->getB();
    opserr << "X: " << theSOE2->getX();

    delete theSOE2;

    PetscSparseSeqSolver *theSolver3 = new PetscSparseSeqSolver(method, preconditioner);
    SparseGenRowLinSOE *theSOE3 = new SparseGenRowLinSOE(*theSolver3);

    opserr << "Test Three with Petsc\n";

    theSOE3->setSize(theGraph);

    theSOE3->addA(a,f);
    theSOE3->addA(a,g);
    theSOE3->addA(c,h);
    theSOE3->addA(c,m,2.0);    

    theSOE3->addB(x,h);
    theSOE3->addB(x,m,-1);

    theSOE2->solve();
    opserr << "B: " << theSOE2->getB();
    opserr << "X: " << theSOE2->getX();
    
    // perturb A a bit
    theSOE3->addA(c,m,0.001);    
    theSOE3->solve();
    opserr << "  results on a second solve with A perturbed\n";
    opserr << "B: " << theSOE2->getB();
    opserr << "X: " << theSOE2->getX();

    // perturb B a bit
    theSOE3->addB(x,m,0.001);
    theSOE3->solve();
    opserr << "  results on a second solve with B perturbed\n";
    opserr << "B: " << theSOE2->getB();
    opserr << "X: " << theSOE2->getX();

    delete theSOE3;
    exit(0);
}
