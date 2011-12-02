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
                                                                        
// $Revision: 1.4 $
// $Date: 2005-11-23 18:24:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/actor/Actor.h,v $
                                                                        
                                                                        
#ifndef Actor_h
#define Actor_h

// Written: fmk
// Revision: A
//
// Purpose: This file contains the class definition for Actor.
// Actor is meant to be an abstract base class and thus no objects of it's type
// can be instantiated. 
//
// What: "@(#) Actor.h, revA"

class ObjectBroker;
class Message;
class Channel;
class ChannelAddress;
class MovableObject;
class Matrix;
class Vector;
class ID;
class FEM_ObjectBroker;

class ActorMethod
{
  public:
    int tag;
    int (*theMethod)();
};

class Actor
{
  public:
    Actor(Channel &theChannel, 
	  FEM_ObjectBroker &theBroker,
	  int numActorMethods =0);
    
    virtual ~Actor();
    
    virtual int  run(void) = 0;

    virtual int  addMethod(int tag, int (*fp)());
    virtual int  getMethod();
    virtual int  processMethod(int tag);

    virtual int sendObject(MovableObject &theObject, 
			   ChannelAddress *theAddress =0);  
    virtual int recvObject(MovableObject &theObject, 
			   ChannelAddress *theAddress =0);

    virtual int sendMessage(const Message &theMessage, 
			    ChannelAddress *theAddress =0);   
    virtual int recvMessage(Message &theMessage, 
			    ChannelAddress *theAddress =0);  
    
    virtual int sendMatrix(const Matrix &theMatrix, 
			   ChannelAddress *theAddress =0);   
    virtual int recvMatrix(Matrix &theMatrix, 
			   ChannelAddress *theAddress =0);  
    
    virtual int sendVector(const Vector &theVector, 
			   ChannelAddress *theAddress =0);   
    virtual int recvVector(Vector &theVector, 
			   ChannelAddress *theAddress =0);  
    
    virtual int sendID(const ID &theID, 
		       ChannelAddress *theAddress =0);   
    virtual int recvID(ID &theID, 
		       ChannelAddress *theAddress =0);  

    Channel 		*getChannelPtr(void) const;
    FEM_ObjectBroker 	*getObjectBrokerPtr(void) const;    
    ChannelAddress  	*getShadowsAddressPtr(void) const;            

    virtual int barrierCheck(int result);
    void setCommitTag(int commitTag);

  protected:
    FEM_ObjectBroker *theBroker; 
    Channel *theChannel;    
    
  private:	
    int numMethods, maxNumMethods;
    ActorMethod **actorMethods;
    ChannelAddress *theRemoteShadowsAddress;

    int commitTag;
};

#endif

