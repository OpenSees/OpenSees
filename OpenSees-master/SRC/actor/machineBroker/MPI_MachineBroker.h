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
// $Date: 2003-08-29 07:17:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/machineBroker/MPI_MachineBroker.h,v $
                                                                        
                                                                        
// Written: fmk
// Revision: A
//
// Purpose: This file contains the class definition for MPI_MachineBroker.
// MPI_MachineBroker is the broker responsible for monitoring the usage of
// the processes in an mpi run.
//
// What: "@(#) MPI_MachineBroker.h, revA"

#ifndef MPI_MachineBroker_h
#define MPI_MachineBroker_h

#include <MachineBroker.h>
class ID;
class MPI_Channel;
class FEM_ObjectBroker;

class MPI_MachineBroker : public MachineBroker
{
  public:
    MPI_MachineBroker(FEM_ObjectBroker *theBroker, int argc, char **argv);
    ~MPI_MachineBroker();

    // methods to return info about local process id and num processes
    int getPID(void);
    int getNP(void);

    // methods to get and free Channels (processes)
    Channel *getMyChannel(void);
    Channel *getRemoteProcess(void);
    int freeProcess(Channel *);

  protected:
    
  private:
    int rank;
    int size;
    ID *usedChannels;
    MPI_Channel **theChannels;
};

#endif
