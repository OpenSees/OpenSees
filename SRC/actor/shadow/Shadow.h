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
// $Source: /usr/local/cvs/OpenSees/SRC/actor/shadow/Shadow.h,v $
                                                                        
                                                                        
// File: ~/actor/Shadow.h
//
// Written: fmk
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class definition for Shadow.
// Shadow is meant to be an abstract base class and thus no objects of it's type
// can be instantiated. 
//
// What: "@(#) Shadow.h, revA"

#ifndef Shadow_h
#define Shadow_h

#include <bool.h>

class MachineBroker;
class Message;
class Channel;
class ChannelAddress;
class MovableObject;
class Matrix;
class Vector;
class ID;

class FEM_ObjectBroker;

class Shadow
{
  public:
    Shadow(Channel &theChannel, 
	   FEM_ObjectBroker &theBroker,
	   ChannelAddress &theAddress); // if actor process up and running
    
    Shadow(char *program,
	   Channel &theChannel, 
	   FEM_ObjectBroker &theBroker,
	   MachineBroker &theMachineBroker,
	   int compDemand,
           bool startShadow);   // to get an actor process up and running

    virtual ~Shadow();

    virtual int sendObject(MovableObject &theObject);  
    virtual int recvObject(MovableObject &theObject);      
    virtual int sendMessage(const Message &theMessage);  
    virtual int recvMessage(Message &theMessage);  
    virtual int sendMatrix(const Matrix &theMatrix);  
    virtual int recvMatrix(Matrix &theMatrix);      
    virtual int sendVector(const Vector &theVector);  
    virtual int recvVector(Vector &theVector);      
    virtual int sendID(const ID &theID);  
    virtual int recvID(ID &theID);      

    Channel 		  *getChannelPtr(void) const;
    FEM_ObjectBroker 	*getObjectBrokerPtr(void) const;        
    ChannelAddress 	  *getActorAddressPtr(void) const;            

  protected:
    FEM_ObjectBroker      *theBroker;
    Channel	 	  *theChannel;    

  private:
    ChannelAddress 	  *theRemoteActorsAddress;    
};

#endif
