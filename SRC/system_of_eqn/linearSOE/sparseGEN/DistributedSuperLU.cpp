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
// $Date: 2005-12-06 22:21:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/DistributedSuperLU.cpp,v $
                                                                        
                                                                        
// Written: fmk 
//
// Description: This file contains the implementation of DistributedSuperLU
//
// What: "@(#) DistributedSuperLU.h, revA"

#include <DistributedSuperLU.h>
#include <SparseGenColLinSOE.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>

#include <superlu_ddefs.h>

#ifdef SUPERLU_DIST_MAJOR_VERSION 
  // SuperLU_Dist 6.2.0 brought a major change that separated several structures tobe "precision-dependent"
  // which changes the name (prefixes with 'd' or 'z')
  #if (SUPERLU_DIST_MAJOR_VERSION >= 6 && SUPERLU_DIST_MINOR_VERSION >= 2)
    #include <superlu_FCnames.h>
    #ifndef _SUPERLU_DIST_6 
      #define _SUPERLU_DIST_6
    #endif
  #endif

  // SuperLU_DIST 'options' was redefined starting in Version 5.X.X
  #if SUPERLU_DIST_MAJOR_VERSION >= 5
      superlu_dist_options_t options;
  #else
      superlu_options_t options;
  #endif

#else

  #if defined(SUPERLU_DIST_MAJOR_VERSION) && SUPERLU_DIST_MAJOR_VERSION >= 5
    superlu_dist_options_t options;
  #else
    superlu_options_t options;
  #endif

#endif

SuperLUStat_t stat;
SuperMatrix A;
gridinfo_t grid;
MPI_Comm comm_SuperLU;

#ifdef _SUPERLU_DIST_6
  dScalePermstruct_t ScalePermstruct;
  dLUstruct_t LUstruct;
#else
  ScalePermstruct_t ScalePermstruct;
  LUstruct_t LUstruct;
#endif


DistributedSuperLU::DistributedSuperLU(int npR, int npC)
  :SparseGenColLinSolver(SOLVER_TAGS_DistributedSuperLU), 
   gridInit(false), npRow(npR), npCol(npC),
   processID(0), numChannels(0), theChannels(0)
{
  opserr << "DistributedSuperLU::DistributedSuperLU(int npR, int npC)\n";
}


DistributedSuperLU::DistributedSuperLU()
  :SparseGenColLinSolver(SOLVER_TAGS_DistributedSuperLU), 
   gridInit(false), npRow(0), npCol(0),
   processID(0), numChannels(0), theChannels(0)
{
  opserr << "DistributedSuperLU::DistributedSuperLU()\n";
}


DistributedSuperLU::~DistributedSuperLU()
{
  //Destroy_LU(theSOE->size, &grid, &LUstruct); 
  #ifdef _SUPERLU_DIST_6
    dScalePermstructFree(&ScalePermstruct);
    dLUstructFree(&LUstruct); 
  #else
    ScalePermstructFree(&ScalePermstruct);
    LUstructFree(&LUstruct); 
  #endif

  //superlu_gridexit(&grid);

  if (theChannels != 0)
    delete [] theChannels;

  //  MPI_Comm_free(&comm_SuperLU);
}

int
DistributedSuperLU::solve(void)
{
  if (theSOE == 0) {
    opserr << "WARNING DistributedSuperLU::solve(void)- ";
    opserr << " No LinearSOE object has been set\n";
    return -1;
  }

  // check for quick return
  if (theSOE->size == 0)
    return 0;

  // if subprocess recv A & B from p0
  if (processID != 0) {
    Channel *theChannel = theChannels[0];
    theChannel->recvVector(0, 0, (*theSOE->vectB));
    Vector vectA(theSOE->A, theSOE->nnz);    
    theChannel->recvVector(0, 0, vectA);
  } 

  //
  // if main process, send B & A to all, solve and send back X, B & result

  else {

    Vector vectA(theSOE->A, theSOE->nnz);    

    // send B & A to p1 through n-1
    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      theChannel->sendVector(0, 0, *(theSOE->vectB));
      theChannel->sendVector(0, 0, vectA);
    }
  }

  int info;

  int iam = grid.iam;

  if (iam < (npRow * npCol)) {
    int n = theSOE->size;
    int nnz = theSOE->nnz;
    int ldb = n;
    int nrhs = 1;
    static double berr[1];

    // first copy B into X
    double *Xptr = theSOE->X;
    double *Bptr = theSOE->B;

    for (int i=0; i<n; i++)
      *(Xptr++) = *(Bptr++);
    Xptr = theSOE->X;

    //
    // set the Fact options:
    //   leave DOFACT if never been factored
    //   otherwise use SamePattern if factored at least once
    //

    if ((options.Fact == FACTORED) && (theSOE->factored == false)) {
      options.Fact = SamePattern;
      for (int i=0; i<nnz; i++) rowA[i] = theSOE->rowA[i];
    }

    //
    // solve & mark as factored
    //

    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, Xptr, ldb, nrhs, &grid,
		     &LUstruct, berr, &stat, &info);

    if (theSOE->factored == false) {
      options.Fact = FACTORED;      
      theSOE->factored = true;
    }
  }

  /*
  // synch processes again
  if (processID != 0) {
    static ID idData(1); 
    idData(0) = info;
    Channel *theChannel = theChannels[0];
    theChannel->sendID(0, 0, idData);
  } 
  else {
    for (int j=0; j<numChannels; j++) {
      static ID idData(1);
      Channel *theChannel = theChannels[j];
      theChannel->recvID(0, 0, idData);
    }
  }
  */

  return 0;
}

int
DistributedSuperLU::setSize()
{
  int n = theSOE->size;

  //
  // init the super lu process grid
  //
  if (gridInit == false) {

    MPI_Comm  comm_world;  //
    MPI_Group group_world; // group_worker;

    comm_world = MPI_COMM_WORLD;
    MPI_Comm_group(comm_world, &group_world);
    //    MPI_Group_excl(group_world, 1, 0, &group_worker);  /* process 0 not member */

    MPI_Comm_create(comm_world, group_world, &comm_SuperLU);

    superlu_gridinit(comm_SuperLU, npRow, npCol, &grid);

  // free old structures if resize already called
  } else {
    
    #ifdef _SUPERLU_DIST_6
      dDestroy_LU(theSOE->size, &grid, &LUstruct); 
      dScalePermstructFree(&ScalePermstruct);
      dLUstructFree(&LUstruct); 
    #else
      Destroy_LU(theSOE->size, &grid, &LUstruct); 
      ScalePermstructFree(&ScalePermstruct);
      LUstructFree(&LUstruct); 
    #endif
  }
  
  //
  // Initialize the statistics variables.
  //
  PStatInit(&stat);
  
  //
  // Create compressed column matrix for A. 
  //
			      
  if (n > 0) {
    
    // create the SuperMatrix A	
    int nnz = theSOE->nnz;
    rowA = new int[nnz];
    for (int i=0; i<nnz; i++) rowA[i] = theSOE->rowA[i];

    dCreate_CompCol_Matrix_dist(&A, n, n, nnz, theSOE->A, 
				rowA, theSOE->colStartA, 
				SLU_NC, SLU_D, SLU_GE);


    //
    // Initialize ScalePermstruct and LUstruct.
    //
    #ifdef _SUPERLU_DIST_6
      dScalePermstructInit(n, n, &ScalePermstruct);
      dLUstructInit(n, &LUstruct);
    #else
      ScalePermstructInit(n, n, &ScalePermstruct);
      LUstructInit(n, &LUstruct);
    #endif
  }  
			      
			      
  //
  //  set default options
  //
  set_default_options_dist(&options);    
  options.PrintStat=NO;

  //
  // Initialize the statistics variables. 
  //
  PStatInit(&stat);

  return 0;
}


int
DistributedSuperLU::setProcessID(int dTag) 
{
  processID = dTag;
  return 0;
}

int
DistributedSuperLU::setChannels(int nChannels, Channel **theC)
{
  numChannels = nChannels;

  if (theChannels != 0)
    delete [] theChannels;

  theChannels = new Channel *[numChannels];
  if (theChannels == 0) {
    opserr << "DistributedSuperLU::sendSelf() - failed to allocate channel array of size: " << 
      numChannels << endln;
    return -1;
  }

  for (int i=0; i<numChannels; i++)
    theChannels[i] = theC[i];

  return 0;
}

int
DistributedSuperLU::sendSelf(int cTag, Channel &theChannel)
{
  int sendID =0;

  // if P0 check if already sent. If already sent use old processID; if not allocate a new process 
  // id for remote part of object, enlarge channel * to hold a channel * for this remote object.

  // if not P0, send current processID

  if (processID == 0) {

    // check if already using this object
    bool found = false;
    for (int i=0; i<numChannels; i++)
      if (theChannels[i] == &theChannel) {
	sendID = i+1;
	found = true;
      }

    // if new object, enlarge Channel pointers to hold new channel * & allocate new ID
    if (found == false) {
      int nextNumChannels = numChannels + 1;
      Channel **nextChannels = new Channel *[nextNumChannels];
      if (nextNumChannels == 0) {
	opserr << "DistributedSuperLU::sendSelf() - failed to allocate channel array of size: " << 
	  nextNumChannels << endln;
	return -1;
      }
      for (int i=0; i<numChannels; i++)
	nextChannels[i] = theChannels[i];
      nextChannels[numChannels] = &theChannel;

      numChannels = nextNumChannels;
      
      if (theChannels != 0)
	delete [] theChannels;

      theChannels = nextChannels;
      
      // allocate new processID for remote object
      sendID = numChannels;
    }

  } else 
    sendID = processID;
    opserr << "DistributedSuperLU::sendSelf(int cTag, Channel &theChannel) - 5\n";
  static ID idData(3);
  idData(0) = sendID;
  idData(1) = npRow;
  idData(2) = npCol;
  int res = theChannel.sendID(0, cTag, idData);
  if (res < 0) {
    opserr <<"WARNING DistributedSuperLU::sendSelf() - failed to send data\n";
    return -1;
  }	      

  return 0;
}

int
DistributedSuperLU::recvSelf(int cTag, 
			     Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  opserr << "DistributedSuperLU::recvSelf(int cTag, Channel &theChannel) - START\n";
  static ID idData(3);

  int res = theChannel.recvID(0, cTag, idData);
  if (res < 0) {
    opserr <<"WARNING DistributedSuperLU::recvSelf() - failed to receive data\n";
    return -1;
  }	      

  processID = idData(0);
  npRow = idData(1);
  npCol = idData(2);

  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  opserr << "DistributedSuperLU::recvSelf(int cTag, Channel &theChannel) - END\n";
  return 0;
}











