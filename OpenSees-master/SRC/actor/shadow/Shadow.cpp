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
                                                                        
// $Revision: 1.7 $
// $Date: 2009-05-11 21:28:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/shadow/Shadow.cpp,v $
                                                                        

// Written: fmk
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
	       FEM_ObjectBroker &myBroker)
  :theChannel(&theChan), theObjectBroker(&myBroker), theMachineBroker(0),
   theRemoteActorsAddress(0), commitTag(0)
{
  if (theChannel->setUpConnection() != 0)  {
    opserr << "Shadow::Shadow() "
	   << "- failed to setup connection\n";
    exit(-1);
  }
}


Shadow::Shadow(Channel &theChan, 
	       FEM_ObjectBroker &myBroker,
	       ChannelAddress &theAddress)
  :theChannel(&theChan), theObjectBroker(&myBroker), theMachineBroker(0),
   theRemoteActorsAddress(&theAddress), commitTag(0)
{
  if (theChannel->setUpConnection() != 0)  {
    opserr << "Shadow::Shadow() "
	   << "- failed to setup connection\n";
    exit(-1);
  }
}

Shadow::Shadow(int actorType,
	       FEM_ObjectBroker &myBroker,	       
	       MachineBroker &theMachineBrokr,
	       int compDemand)
  :theObjectBroker(&myBroker), theMachineBroker(&theMachineBrokr), 
   theRemoteActorsAddress(0), commitTag(0)
{
    // start the remote actor process running
  theChannel = theMachineBroker->startActor(actorType, compDemand);
  if (theChannel == 0) {
    opserr << "Shadow::Shadow - could not start remote actor\n";
    opserr << " using program " << actorType << endln;
    exit(-1);
  }
  
  // now call setUpShadow on the channel
  if (theChannel->setUpConnection() != 0)  {
    opserr << "Shadow::Shadow() "
	   << "- failed to setup connection\n";
    exit(-1);
  }
  theRemoteActorsAddress = theChannel->getLastSendersAddress();
}

Shadow::~Shadow()
{
  if (theMachineBroker != 0)
    theMachineBroker->finishedWithActor(theChannel);    
}    

int
Shadow::sendObject(MovableObject &theObject)
{
    return theChannel->sendObj(commitTag, theObject, theRemoteActorsAddress);
}

int
Shadow::recvObject(MovableObject &theObject)
{
    return theChannel->recvObj(commitTag, theObject,*theObjectBroker, theRemoteActorsAddress);
}


int
Shadow::recvMessage(Message &theMessage)
{
    return theChannel->recvMsg(0, commitTag, theMessage, theRemoteActorsAddress);
}

int
Shadow::sendMessage(const Message &theMessage)
{
    return theChannel->sendMsg(0, commitTag, theMessage, theRemoteActorsAddress);
}

int
Shadow::sendMatrix(const Matrix &theMatrix)
{
    return theChannel->sendMatrix(0, commitTag, theMatrix, theRemoteActorsAddress);
}

int
Shadow::recvMatrix(Matrix &theMatrix)
{
    return theChannel->recvMatrix(0, commitTag, theMatrix, theRemoteActorsAddress);
}

int
Shadow::sendVector(const Vector &theVector)
{
    return theChannel->sendVector(0, commitTag, theVector, theRemoteActorsAddress);
}

int
Shadow::recvVector(Vector &theVector)
{
    return theChannel->recvVector(0, commitTag, theVector, theRemoteActorsAddress);
}

int
Shadow::sendID(const ID &theID)
{
    return theChannel->sendID(0, commitTag, theID, theRemoteActorsAddress);
}

int
Shadow::recvID(ID &theID)
{
    return theChannel->recvID(0, commitTag, theID, theRemoteActorsAddress);
}


void
Shadow::setCommitTag(int tag)
{
  commitTag = tag;
}


Channel *
Shadow::getChannelPtr(void) const
{
    return theChannel;
}

FEM_ObjectBroker *
Shadow::getObjectBrokerPtr(void) const
{
    return theObjectBroker;
}

ChannelAddress *
Shadow::getActorAddressPtr(void) const
{
    return theRemoteActorsAddress;
}


