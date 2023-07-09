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
// $Date: 2009-05-11 21:14:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/machineBroker/MPI_MachineBroker.cpp,v $
                                                                        
                                                                        
// Written: fmk
// Revision: A


#include <FEM_ObjectBroker.h>
#include <MPI_MachineBroker.h>
#include <MPI_Channel.h>
#include <ID.h>

#include <mpi.h>

MPI_MachineBroker::MPI_MachineBroker(FEM_ObjectBroker *theBroker, int argc, char **argv)
  :MachineBroker(theBroker)
{
  int flag = 0;
  MPI_Initialized(&flag);
  if (!flag) {
      MPI_Init(&argc, &argv);
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  theChannels = new MPI_Channel *[size];
  for (int i=0; i<size; i++) {
    theChannels[i] = new MPI_Channel(i);
  }
  usedChannels = new ID(size);
  usedChannels->Zero();
}


MPI_MachineBroker::~MPI_MachineBroker()
{
  for (int i=0; i<size; i++) {
      delete theChannels[i]; 
  }

  delete [] theChannels;
  delete usedChannels;

  MPI_Finalize();
}


int 
MPI_MachineBroker::getPID(void)
{
  return rank;
}


int 
MPI_MachineBroker::getNP(void)
{
  return size;
}



Channel *
MPI_MachineBroker::getMyChannel(void)
{
  return theChannels[0];
}



Channel *
MPI_MachineBroker::getRemoteProcess(void)
{
  if (rank != 0) {
    opserr << "MPI_MachineBroker::getRemoteProcess() - child process cannot not yet allocate processes\n";
    return 0;
  }
      
  for (int i=0; i<size; i++)
    if (i != rank)
      if ((*usedChannels)(i) == 0) {
	(*usedChannels)(i) = 1;
	return theChannels[i];
      }
  
  // no processes available
  return 0;
}


int 
MPI_MachineBroker::freeProcess
(Channel *theChannel)
{
  for (int i=0; i<size; i++)
    if (i != rank)
      if (theChannels[i] == theChannel) {
	(*usedChannels)(i) = 0;
	return 0;
      }
  
  // channel not found!
  return -1;
}

