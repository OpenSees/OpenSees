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
// $Date: 2007-06-01 17:06:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/example.cpp,v $

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <MumpsParallelSOE.h>
#include <MumpsParallelSolver.h>
#include <MPI_Channel.h>
#include <FEM_ObjectBroker.h>
#include <Graph.h>
#include <Vertex.h>

// init the global variabled defined in OPS_Globals.h
#include <OPS_Globals.h>
#include <StandardStream.h>
#include <mpi.h>
#include <SimulationInformation.h>


StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;
SimulationInformation simulationInfo;
 
double        ops_Dt = 0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;

int main(int argc, char ** argv)
{
  int ierr, rank, np;
  int soeType =0;  // 0 - unsymmetric, 1 - symmetric positive definite, 2 - symmetric
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &np);

  ID *myID = 0;
  FEM_ObjectBroker theBroker;

  MumpsParallelSolver *theSolver = new MumpsParallelSolver();

  MumpsParallelSOE *theSOE = new MumpsParallelSOE(*theSolver, soeType);

    double numP = np * 1.0;

  if (rank == 0) {
    for (int i=1; i<np; i++) {
      MPI_Channel *theChannel = new MPI_Channel(i);
      theSOE->sendSelf(0, *theChannel);
    }
    myID = new ID(6);
    Graph graph;
    for (int i=0; i<6; i++) {
      Vertex *vertex = new Vertex(i, i);
      graph.addVertex(vertex);
      (*myID)(i) = i;
    }


    Matrix a(6,6);
    a.Zero();
    a(0,0) = 4169.95/numP;
    a(1,0) = 0.0;
    a(2,0) = 10075/numP;
    a(3,0) = -4030.0;
    a(4,0) = 0.0;
    a(5,0) = 0.0;
    a(1,1) = 10084/numP;
    a(2,1) = -1612/numP;
    a(3,1) = 0.0;
    a(4,1) = -8.9556;
    a(5,1) = -1612;
    a(2,2) = 1354080.0/numP;
    a(3,2) = 0.0;
    a(4,2) = 1612;
    a(5,2) = 193440.0;
    a(3,3) = 4169.93;
    a(4,3) = 0.0;
    a(5,3) = 10075;
    a(4,4) = 10084;
    a(5,4) = 1612;
    a(5,5) = 1354000.0;

    Vector b(6);
    b += 2.0;

    for (int i=0; i<6; i++)
      for (int j=i+1; j<6; j++) {
	graph.addEdge(i,j);
	a(i,j) = a(j,i);
      }

    theSOE->setSize(graph);

    theSOE->zeroA();
    theSOE->zeroB();
    theSOE->addA(a, *myID);
    theSOE->addB(b, *myID);

  } else {

    MPI_Channel *theChannel = new MPI_Channel(0);
    theSOE->recvSelf(0, *theChannel, theBroker);

    myID = new ID(3);
    Graph graph;
    for (int i=0; i<3; i++) {
      Vertex *vertex = new Vertex(i, i);
      graph.addVertex(vertex);
      (*myID)(i) = i;
    }

    Matrix a(3,3);
    Vector b(3);
    a(0,0) = 4169.95/numP;
    a(1,0) = 0.0;
    a(2,0) = 10075.0/numP;
    a(1,1) = 10084.0/numP;
    a(2,1) = -1612.0/numP;
    a(2,2) = 1354080.0/numP;


    for (int i=0; i<3; i++)
      for (int j=i+1; j<3; j++) {
	graph.addEdge(i,j);
	a(i,j) = a(j,i);
      }


    theSOE->setSize(graph);

    theSOE->zeroA();
    theSOE->zeroB();
    theSOE->addA(a, *myID);
    theSOE->addB(b, *myID);

  }

  Matrix a(6,6);
  a.Zero();
  a(0,0) = 4169.95;
  a(1,0) = 0.0;
  a(2,0) = 10075;
  a(3,0) = -4030.0;
  a(4,0) = 0.0;
  a(5,0) = 0.0;
  a(1,1) = 10084;
  a(2,1) = -1612;
  a(3,1) = 0.0;
  a(4,1) = -8.9556;
  a(5,1) = -1612;
  a(2,2) = 1354080.0;
  a(3,2) = 0.0;
  a(4,2) = 1612;
  a(5,2) = 193440.0;
  a(3,3) = 4169.93;
  a(4,3) = 0.0;
  a(5,3) = 10075;
  a(4,4) = 10084;
  a(5,4) = 1612;
  a(5,5) = 1354000.0;
  for (int i=0; i<6; i++)
    for (int j=i+1; j<6; j++) 
      a(i,j) = a(j,i);
  
  

  Vector b(6),x(6);

  b= theSOE->getB();  

  theSOE->solve();
  x = theSOE->getX();
  
  if (rank == np-1) {
    opserr << "A:\n";
    opserr << "B:\n" << b;
    opserr << "X:\n " << x;
    opserr << "A*X:\n" << a*x;
  }
  
  theSOE->zeroB();
  if (rank == 0) {
    Vector b(6);
    b += 2.0/numP;
    theSOE->addB(b, *myID);
  } else {
    Vector b(3);
    b += 2.0/numP;
    theSOE->addB(b, *myID);
  }

  b= theSOE->getB();  
  theSOE->solve();
  x = theSOE->getX();  

  if (rank == np-1) {
    opserr << "A:\n";
    opserr << "B:\n" << b;
    opserr << "X:\n " << x;
    opserr << "A*X:\n" << a*x;
  }

  delete theSOE;

  MPI_Finalize();
  return 0;
}
