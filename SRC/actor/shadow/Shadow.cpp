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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/shadow/Shadow.cpp,v $
                                                                        
                                                                        
// File: ~/actor/Shadow.C
//
// Written: fmk
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementation of Shadow.
//
// What: "@(#) Shadow.C, revA"

#include <Shadow.h>
#include <stdlib.h>

#include <Channel.h>
#include <MachineBroker.h>
#include <Message.h>
#include <MovableObject.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>


Shadow::Shadow(Channel &theChan, 
	       FEM_ObjectBroker &myBroker,
	       ChannelAddress &theAddress)
:theChannel(&theChan),theBroker(&myBroker),theRemoteActorsAddress(&theAddress)
{

}

Shadow::Shadow(char *program,
	       Channel &theChan, 
	       FEM_ObjectBroker &myBroker,	       
	       MachineBroker &theMachineBroker,
	       int compDemand,
	       bool startShadow)
:theChannel(&theChan),theBroker(&myBroker),theRemoteActorsAddress(0)
{
    // start the remote actor process running
  if (startShadow == true) {
    int res = theMachineBroker.startActor(program,theChan,compDemand);
    if (res < 0) {
	cerr << "Shadow::Shadow - could not start remote actor\n";
	cerr << " using program " << *program << endl;
	exit(-1);
    }
  }

  // now call setUpShadow on the channel
  theChan.setUpShadow();
  theRemoteActorsAddress = theChan.getLastSendersAddress();
}

Shadow::~Shadow()
{
    
}    

int
Shadow::sendObject(MovableObject &theObject)
{
    return theChannel->sendObj(0, theObject,theRemoteActorsAddress);
}

int
Shadow::recvObject(MovableObject &theObject)
{
    return theChannel->recvObj(0, theObject,*theBroker,theRemoteActorsAddress);
}


int
Shadow::recvMessage(Message &theMessage)
{
    return theChannel->recvMsg(0, 0, theMessage,theRemoteActorsAddress);
}

int
Shadow::sendMessage(const Message &theMessage)
{
    return theChannel->sendMsg(0, 0, theMessage,theRemoteActorsAddress);
}



int
Shadow::sendMatrix(const Matrix &theMatrix)
{
    return theChannel->sendMatrix(0, 0, theMatrix,theRemoteActorsAddress);
}

int
Shadow::recvMatrix(Matrix &theMatrix)
{
    return theChannel->recvMatrix(0, 0, theMatrix,theRemoteActorsAddress);
}

int
Shadow::sendVector(const Vector &theVector)
{
    return theChannel->sendVector(0, 0, theVector,theRemoteActorsAddress);
}

int
Shadow::recvVector(Vector &theVector)
{
    return theChannel->recvVector(0, 0, theVector,theRemoteActorsAddress);
}

int
Shadow::sendID(const ID &theID)
{
    return theChannel->sendID(0, 0, theID,theRemoteActorsAddress);
}

int
Shadow::recvID(ID &theID)
{
    return theChannel->recvID(0, 0, theID,theRemoteActorsAddress);
}


Channel *
Shadow::getChannelPtr(void) const
{
    return theChannel;
}

FEM_ObjectBroker *
Shadow::getObjectBrokerPtr(void) const
{
    return theBroker;
}

ChannelAddress *
Shadow::getActorAddressPtr(void) const
{
    return theRemoteActorsAddress;
}


