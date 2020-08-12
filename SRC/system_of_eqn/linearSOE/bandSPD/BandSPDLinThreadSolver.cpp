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
// $Date: 2003-02-14 23:02:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandSPD/BandSPDLinThreadSolver.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/bandSPD/BandSPDLinThreadSolver.h
//
// Written: fmk 
// Created: Mar, 1998
// Revision: A
//
// Description: This file contains the class definition for 
// BandSPDLinThreadSolver. It solves the BandSPDLinSOE object by calling
// Thread routines.
//
// What: "@(#) BandSPDLinThreadSolver.h, revA"

#include <BandSPDLinThreadSolver.h>
#include <BandSPDLinSOE.h>
#include <f2c.h>
#include <math.h>
#include <thread.h>

// global data that will be needed by the threads
struct thread_control_block {
  mutex_t start_mutex;
  mutex_t end_mutex;
  mutex_t startBlock_mutex;
  mutex_t endBlock_mutex;

  cond_t start_cond;
  cond_t end_cond;
  cond_t startBlock_cond;
  cond_t endBlock_cond;

  double *A,*X,*B;
  int n, kd;
  int blockSize, currentBlock;
  int workToDo;
  int info;
  int threadID;
  int numThreads;
  int threadsDone;
  int threadsDoneBlock;
} TCB_BandSPDLinThreadSolver;


BandSPDLinThreadSolver::BandSPDLinThreadSolver()
:BandSPDLinSolver(SOLVER_TAGS_BandSPDLinThreadSolver), NP(1), 
 running(false), blockSize(1)
{
  
}

BandSPDLinThreadSolver::BandSPDLinThreadSolver(int numProcessors, int blckSize)
:BandSPDLinSolver(SOLVER_TAGS_BandSPDLinThreadSolver), NP(numProcessors),
 running(false), blockSize(blckSize)
{

}

BandSPDLinThreadSolver::~BandSPDLinThreadSolver()
{
    
}


extern "C" int dpbsv_(char *UPLO, int *N, int *KD, int *NRHS, 
		      double *A, int *LDA, double *B, int *LDB, 
		      int *INFO);

extern "C" int dpbtrf_(char *UPLO, int *N, int *KD, double *A, 
		       int *LDA, int *INFO);

extern "C" int dpbtrs_(char *UPLO, int *N, int *KD, int *NRHS, 
		       double *A, int *LDA, double *B, int *LDB, 
		       int *INFO);
		       

void *BandSPDLinThreadSolver_Worker(void *arg);


int
BandSPDLinThreadSolver::solve(void)
{
    if (theSOE == 0) {
	opserr << "BandSPDLinThreadSolver::solve(void)- ";
	opserr << " No LinearSOE object has been set\n";
	return -1;
    }

    int n = theSOE->size;
    int kd = theSOE->half_band -1;
    int ldA = kd +1;
    int nrhs = 1;
    int ldB = n;
    int info;
    double *Aptr = theSOE->A;
    double *Xptr = theSOE->X;
    double *Bptr = theSOE->B;

    // first copy B into X
    for (int i=0; i<n; i++)
	*(Xptr++) = *(Bptr++);
    Xptr = theSOE->X;

    // now solve AX = Y
    if (theSOE->factored == false) {

      // start threads running
      if (running == false) {
	mutex_init(&TCB_BandSPDLinThreadSolver.start_mutex, USYNC_THREAD, 0);
	mutex_init(&TCB_BandSPDLinThreadSolver.end_mutex, USYNC_THREAD, 0);
	mutex_init(&TCB_BandSPDLinThreadSolver.startBlock_mutex, USYNC_THREAD, 0);
	mutex_init(&TCB_BandSPDLinThreadSolver.endBlock_mutex, USYNC_THREAD, 0);
	cond_init(&TCB_BandSPDLinThreadSolver.start_cond, USYNC_THREAD, 0);
	cond_init(&TCB_BandSPDLinThreadSolver.end_cond, USYNC_THREAD, 0);
	cond_init(&TCB_BandSPDLinThreadSolver.startBlock_cond, USYNC_THREAD, 0);
	cond_init(&TCB_BandSPDLinThreadSolver.endBlock_cond, USYNC_THREAD, 0);
	TCB_BandSPDLinThreadSolver.workToDo = 0;
        
	// Create the threads - Bound daemon threads
	for (int j = 0; j < NP; j++) 
	  thr_create(NULL,0, BandSPDLinThreadSolver_Worker, NULL, 
		     THR_BOUND|THR_DAEMON, NULL);

        running = true;
      }
      
  
      mutex_lock(&TCB_BandSPDLinThreadSolver.start_mutex); 
        // assign the global variables the threads will need
        TCB_BandSPDLinThreadSolver.A = Aptr;
	TCB_BandSPDLinThreadSolver.X = Xptr;
	TCB_BandSPDLinThreadSolver.B = Bptr;
	TCB_BandSPDLinThreadSolver.n = n;
	TCB_BandSPDLinThreadSolver.kd = kd;
	TCB_BandSPDLinThreadSolver.info = 0;
        TCB_BandSPDLinThreadSolver.workToDo = NP;
        TCB_BandSPDLinThreadSolver.blockSize = blockSize;
        TCB_BandSPDLinThreadSolver.currentBlock = -1;
	TCB_BandSPDLinThreadSolver.threadID = 0;
	TCB_BandSPDLinThreadSolver.numThreads = NP;
        TCB_BandSPDLinThreadSolver.threadsDone = 0;
        TCB_BandSPDLinThreadSolver.threadsDoneBlock = NP;

        // wake up the threads
        cond_broadcast(&TCB_BandSPDLinThreadSolver.start_cond);
      mutex_unlock(&TCB_BandSPDLinThreadSolver.start_mutex); 

      // yield the LWP
      thr_yield();
 
      // wait for all the threads to finish
      mutex_lock(&TCB_BandSPDLinThreadSolver.end_mutex);

         while (TCB_BandSPDLinThreadSolver.threadsDone < NP) 
           cond_wait(&TCB_BandSPDLinThreadSolver.end_cond, &TCB_BandSPDLinThreadSolver.end_mutex);

      mutex_unlock(&TCB_BandSPDLinThreadSolver.end_mutex);
    }

    // solve using factored matrix
    dpbtrs_("U",&n,&kd,&nrhs,Aptr,&ldA,Xptr,&ldB,&info);
    

    // check if successful
    if (info != 0)
	return info;

    theSOE->factored = true;
    return 0;
}
    


int
BandSPDLinThreadSolver::setSize()
{
    // nothing to do    
    return 0;
}


int
BandSPDLinThreadSolver::sendSelf(Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}

int
BandSPDLinThreadSolver::recvSelf(Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}




// Thread routine called from thr_create() as a Bound Daemon Thread
void *BandSPDLinThreadSolver_Worker(void *arg)
{
  // Do this loop forever - or until all the Non-Daemon threads have exited
  while(true)
    {

      // Wait for some work to do
      mutex_lock(&TCB_BandSPDLinThreadSolver.start_mutex);

        while (TCB_BandSPDLinThreadSolver.workToDo == 0)
          cond_wait(&TCB_BandSPDLinThreadSolver.start_cond, &TCB_BandSPDLinThreadSolver.start_mutex);

        TCB_BandSPDLinThreadSolver.workToDo--;

        // get an id;
	int myID = TCB_BandSPDLinThreadSolver.threadID++;

      mutex_unlock(&TCB_BandSPDLinThreadSolver.start_mutex);

      // do some initialisation for the factorization
      int currentBlock =0;
      int n = TCB_BandSPDLinThreadSolver.n;
      int kd = TCB_BandSPDLinThreadSolver.kd;
      int blckSize = TCB_BandSPDLinThreadSolver.blockSize;
      double *Aptr = TCB_BandSPDLinThreadSolver.A;
      double *Bptr = TCB_BandSPDLinThreadSolver.B;
      double *Xptr = TCB_BandSPDLinThreadSolver.X;
      int numThreads = TCB_BandSPDLinThreadSolver.numThreads;
    
      int nBlocks = n/blckSize;

      for (int i=0; i<nBlocks; i++) {
        int currentDiagThread = i%numThreads;
        if (myID == currentDiagThread) {

	  // wait till last block row is done
          mutex_lock(&TCB_BandSPDLinThreadSolver.endBlock_mutex);
            while(TCB_BandSPDLinThreadSolver.threadsDoneBlock < numThreads)
	      cond_wait(&TCB_BandSPDLinThreadSolver.endBlock_cond, &TCB_BandSPDLinThreadSolver.endBlock_mutex);

	    TCB_BandSPDLinThreadSolver.threadsDoneBlock = 0;
          mutex_unlock(&TCB_BandSPDLinThreadSolver.endBlock_mutex);        

          // FACTOR DIAG BLOCK
       
          // allow other threads to now proceed
          mutex_lock(&TCB_BandSPDLinThreadSolver.startBlock_mutex);
	    TCB_BandSPDLinThreadSolver.currentBlock = i;

            cond_broadcast(&TCB_BandSPDLinThreadSolver.startBlock_cond);
          mutex_unlock(&TCB_BandSPDLinThreadSolver.startBlock_mutex);        

          
          // FACTOR REST OF BLOCK ROW BELONGING TO THREAD

          // mark that the thread has finished with row
          mutex_lock(&TCB_BandSPDLinThreadSolver.endBlock_mutex);
	    TCB_BandSPDLinThreadSolver.threadsDoneBlock++;

            cond_signal(&TCB_BandSPDLinThreadSolver.endBlock_cond);
          mutex_unlock(&TCB_BandSPDLinThreadSolver.endBlock_mutex);                   
	}
        else {

          // wait till diag i is done 
          mutex_lock(&TCB_BandSPDLinThreadSolver.startBlock_mutex);
             while (TCB_BandSPDLinThreadSolver.currentBlock < i)
                cond_wait(&TCB_BandSPDLinThreadSolver.startBlock_cond, 
			  &TCB_BandSPDLinThreadSolver.startBlock_mutex); 
          mutex_unlock(&TCB_BandSPDLinThreadSolver.startBlock_mutex);

          // FACTOR REST OF BLOCK ROW BELONGING TO THREAD

          // mark that the thread has finished with row
          mutex_lock(&TCB_BandSPDLinThreadSolver.endBlock_mutex);
	    TCB_BandSPDLinThreadSolver.threadsDoneBlock++;

            cond_signal(&TCB_BandSPDLinThreadSolver.endBlock_cond);
          mutex_unlock(&TCB_BandSPDLinThreadSolver.endBlock_mutex);                             
	}
      }

      // signal the main thread that this thread is done with the work
      mutex_lock(&TCB_BandSPDLinThreadSolver.end_mutex);
         TCB_BandSPDLinThreadSolver.threadsDone++;

         cond_signal(&TCB_BandSPDLinThreadSolver.end_cond);
      mutex_unlock(&TCB_BandSPDLinThreadSolver.end_mutex);
    } 
  return NULL;
}







