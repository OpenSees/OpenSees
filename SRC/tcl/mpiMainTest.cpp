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
                                                                        
// $Revision: 1.1 $
// $Date: 2006-01-13 19:35:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/mpiMainTest.cpp,v $

extern "C" {
#include <tcl.h>
#include <tk.h>
}

#include <stdio.h>
#include <string.h>

#include <PartitionedDomain.h>
#include <MPI_MachineBroker.h>
#include <ShadowSubdomain.h>
#include <ActorSubdomain.h>
#include <FEM_ObjectBroker.h>

extern PartitionedDomain theDomain;

extern int OPS_PARALLEL_PROCESSING;
extern int OPS_NUM_SUBDOMAINS;
extern bool OPS_PARTITIONED;
extern FEM_ObjectBroker *OPS_OBJECT_BROKER;
extern MachineBroker    *OPS_MACHINE;

static MPI_MachineBroker *theMachineBroker = 0;

#include <Node.h>
#include <Matrix.h>

int
main(int argc, char **argv)
{
  FEM_ObjectBroker theBroker;
  MPI_MachineBroker theMachine(&theBroker, argc, argv);
  theMachineBroker = &theMachine;

  int rank = theMachine.getPID();
  int np = theMachine.getNP();

  //
  // depending on rank we do something
  //
  if (rank != 0) {

    //
    // on secondary processes we spin waiting to create & run actors
    //
    fprintf(stderr, "Secondary Process Running\n");
    theMachine.runActors();

  } else {

    //
    // on process 0 we create some ShadowSubdomains & then start the OpenSees interpreter
    //
    fprintf(stderr, "Primary Process Running OpenSees Interpreter\n");   

    //
    // set some global parameters
    //
    OPS_OBJECT_BROKER = &theBroker;
    OPS_MACHINE = &theMachine;
    OPS_PARALLEL_PROCESSING = np;
    OPS_NUM_SUBDOMAINS = np - 1;
    OPS_PARTITIONED    = false;
    Channel  **OPS_theChannels = 0;
    OPS_theChannels = new Channel *[OPS_NUM_SUBDOMAINS];
    Subdomain **OPS_theSubdomains = 0;
    OPS_theSubdomains = new Subdomain *[OPS_NUM_SUBDOMAINS];

    // start the remote actors
    for (int i=1; i<=OPS_NUM_SUBDOMAINS; i++) {
      ShadowSubdomain *theSubdomain = new ShadowSubdomain(i, 
							  *OPS_MACHINE, 
							  *OPS_OBJECT_BROKER);
      theDomain.addSubdomain(theSubdomain);
      OPS_theChannels[i-1] = theSubdomain->getChannelPtr();
      OPS_theSubdomains[i-1] = theSubdomain;
    }

    Subdomain *theSubdomain = OPS_theSubdomains[0];
    Node *node1 = new Node(1, 2,   0.0,  0.0);
    Node *node2 = new Node(2, 2, 144.0,  0.0);
    Node *node3 = new Node(3, 2, 168.0,  0.0);    
    Node *node4 = new Node(4, 2,  72.0, 96.0);        
    Matrix mass(2,2); mass(0,0)=1.0; mass(1,1)=2.0;
    node1->setMass(mass);
    node2->setMass(mass);

    theSubdomain->addNode(node1);
    theSubdomain->addExternalNode(node2);
    //    theDomain.addNode(node2);
    theDomain.addNode(node3);
    theDomain.addNode(node4);
    theDomain.Print(opserr);

    // some clean up to shut the remotes down if still running
    theDomain.clearAll();

    // test that everything is shutting down correctly
  }

  //
  // mpi clean up
  //

  theMachine.shutdown();

  fprintf(stderr, "Process Terminating %d\n", rank);

  return 0;
}

int OpenSeesExit(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  //
  // mpi clean up
  //

  if (theMachineBroker != 0) {
    theMachineBroker->shutdown();
    fprintf(stderr, "Process Terminating\n");
  }

  return 0;
}
