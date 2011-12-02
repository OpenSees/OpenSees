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
                                                                        
// $Revision: 1.5 $
// $Date: 2008-09-23 22:47:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/actor/Actor.cpp,v $
                                                                        
                                                                        
// Written: fmk
// Revision: A
//
// Purpose: This file contains the implementation of Actor.
//
// What: "@(#) Actor.C, revA"

#include <Actor.h>
#include <Channel.h>
#include <ChannelAddress.h>
#include <Message.h>
#include <MovableObject.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <stdlib.h>

// Actor
//	constructor to init the list.


Actor::Actor(Channel &theChan,
	     FEM_ObjectBroker &myBroker,
	     int numActorMethods)
:theBroker(&myBroker), theChannel(&theChan),
 numMethods(0), maxNumMethods(numActorMethods), actorMethods(0), 
 theRemoteShadowsAddress(0), commitTag(0)
{
    // call setUpActor on the channel and get shadows address
    if (theChannel->setUpConnection() != 0)  {
        opserr << "Actor::Actor() "
            << "- failed to setup connection\n";
        exit(-1);
    }
    theRemoteShadowsAddress = theChan.getLastSendersAddress();

    if (numActorMethods != 0)
        actorMethods = new ActorMethod *[numActorMethods];

    if (actorMethods == 0)
        maxNumMethods = 0;

    for (int i=0; i<numMethods; i++)
        actorMethods[i] = 0;
}


Actor::~Actor()
{
  // delete the array of actorMethods if constructed one
  if (actorMethods != 0)
    delete actorMethods;
}


// void AddMethod(int tag, int (*fp)()):
//	Method to add a function to the list of avaiable actor methods.
//	The function will be identified as tag, it is a function with
//	no args that returns an int.

int 
Actor::addMethod(int tag, int (*fp)())
{
    // check we are not over our limit
    if (numMethods >= maxNumMethods)
	return -2;

    // check no other with the same tag exists
    ActorMethod *methodPtr;
    for (int i=0; i<numMethods; i++) {
	methodPtr = actorMethods[i];
	if (methodPtr->tag == tag)
	    return -1;
    }

    // add the new method
    ActorMethod *newMethod = new ActorMethod;
    if (newMethod == 0)
	return -3;
    
    newMethod->tag = tag;
    newMethod->theMethod = fp;

    actorMethods[numMethods] = newMethod;
    numMethods++;
    return 0;
}





// int GetMethod():
//	Method to return the integer tag of the next method the actor
//	has been asked to invoke.

int 
Actor::getMethod()
{
    int method = -1;
    Message msg(&method,1);
    this->recvMessage(msg);
    return method;
}




// int ProcMethod(int tag):
//	Method to process the function whose id is tag.

int 
Actor::processMethod(int tag)
{
    ActorMethod *current =0;

    for (int i=0; i<numMethods; i++)
	if (actorMethods[i]->tag == tag) {
	    current = actorMethods[tag];
	}

    if (current == 0)
	return -1;
    
    return (*current).theMethod();
}



int
Actor::sendObject(MovableObject &theObject, 
		  ChannelAddress *theAddress )
{
    if (theAddress == 0)
	return theChannel->sendObj(commitTag, theObject,theRemoteShadowsAddress);
    else
	return theChannel->sendObj(commitTag, theObject,theAddress);	
}

int
Actor::recvObject(MovableObject &theObject, 
		  ChannelAddress *theAddress )
{
    if (theAddress == 0)
	return theChannel->recvObj(commitTag, theObject,*theBroker,
				   theRemoteShadowsAddress); 
    else
	return theChannel->recvObj(commitTag, theObject,*theBroker,theAddress);	
}


int
Actor::recvMessage(Message &theMessage, ChannelAddress *theAddress )
{
    if (theAddress == 0)
	return theChannel->recvMsg(0, commitTag, theMessage,theRemoteShadowsAddress);
    else
	return theChannel->recvMsg(0, commitTag, theMessage,theAddress);	
}

int
Actor::sendMessage(const Message &theMessage, ChannelAddress *theAddress )
{
    if (theAddress == 0)
	return theChannel->sendMsg(0, commitTag, theMessage,theRemoteShadowsAddress);
    else
	return theChannel->sendMsg(0, commitTag, theMessage,theAddress);	
}



int
Actor::sendMatrix(const Matrix &theMatrix, ChannelAddress *theAddress )
{
    if (theAddress == 0)    
	return theChannel->sendMatrix(0, commitTag, theMatrix,theRemoteShadowsAddress);
    else
	return theChannel->sendMatrix(0, commitTag, theMatrix,theAddress);	
}

int
Actor::recvMatrix(Matrix &theMatrix, ChannelAddress *theAddress )
{
    if (theAddress == 0)
	return theChannel->recvMatrix(0, commitTag, theMatrix,theRemoteShadowsAddress);
    else
	return theChannel->recvMatrix(0, commitTag, theMatrix,theAddress);	
}

int
Actor::sendVector(const Vector &theVector, ChannelAddress *theAddress )
{
    if (theAddress == 0)
	return theChannel->sendVector(0, commitTag, theVector,theRemoteShadowsAddress);
    else
	return theChannel->sendVector(0, commitTag, theVector,theAddress);	
}

int
Actor::recvVector(Vector &theVector, ChannelAddress *theAddress )
{
    if (theAddress == 0)
	return theChannel->recvVector(0, commitTag, theVector,theRemoteShadowsAddress);
    else
	return theChannel->recvVector(0, commitTag, theVector,theAddress);	
}

int
Actor::sendID(const ID &theID, ChannelAddress *theAddress )
{
    if (theAddress == 0)
	return theChannel->sendID(0, commitTag, theID,theRemoteShadowsAddress);
    else
	return theChannel->sendID(0, commitTag, theID,theAddress);	
}

int
Actor::recvID(ID &theID, ChannelAddress *theAddress )
{
    if (theAddress == 0)
	return theChannel->recvID(0, commitTag, theID,theRemoteShadowsAddress);
    else
	return theChannel->recvID(0, commitTag, theID,theAddress);	
}


void
Actor::setCommitTag(int tag)
{
  commitTag = tag;
}


Channel *
Actor::getChannelPtr(void) const
{
    return theChannel;
}

FEM_ObjectBroker *
Actor::getObjectBrokerPtr(void) const
{
    return theBroker;
}



ChannelAddress *
Actor::getShadowsAddressPtr(void) const
{
    return theRemoteShadowsAddress;
}


// barrier check:
//
int
Actor::barrierCheck(int myResult = 0)
{
  int result;
  static ID data(1);
  data(0) = myResult; 
  theChannel->sendID(0, commitTag, data); 
  theChannel->recvID(0, commitTag, data);
  result = data(0);
  return result;
}

