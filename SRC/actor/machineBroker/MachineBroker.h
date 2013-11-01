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
// $Date: 2003-08-29 07:16:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/machineBroker/MachineBroker.h,v $
                                                                        
                                                                        
// Written: fmk
// Revision: A
//
// Purpose: This file contains the class definition for MachineBroker.
// MachineBroker is an abstract base class, a subclass of which must
// be written for each parallel machine. A MachineBroker is responsible
// for getting an actor process running on the parallel machine.
//
// What: "@(#) MachineBroker.h, revA"

#ifndef MachineBroker_h
#define MachineBroker_h

class Channel;
class FEM_ObjectBroker;
class ID;

class MachineBroker
{
  public:
    MachineBroker(FEM_ObjectBroker *theObjectBroker);
    virtual ~MachineBroker();

    // methods to return info about local process id and num processes
    virtual int getPID(void) = 0;
    virtual int getNP(void)  = 0;

    // methods to get and free Actors
    virtual int shutdown(void);    
    virtual int runActors(void);
    virtual Channel *startActor(int actorType, int compDemand = 0);
    virtual int finishedWithActor(Channel *);

    // methods to get and free Channels (processes)
    virtual Channel *getMyChannel(void)        =0;
    virtual Channel *getRemoteProcess(void)    =0;
    virtual int freeProcess(Channel *)         =0;
    
    void setObjectBroker(FEM_ObjectBroker *theBroker);

    /* ************ THE OLD INTERFACE ***************
    virtual int startActor(char *actorProgram, 
			   Channel &theChannel,
			   int compDemand =0) =0;
    ********************************************** */

  protected:
    
  private:
    FEM_ObjectBroker *theObjectBroker;

    Channel **actorChannels; // channels owned with running actor processes
    int numActorChannels;
    int numActiveChannels;
    ID *activeChannels;
};

#endif
