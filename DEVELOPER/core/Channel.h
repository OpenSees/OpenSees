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
// $Date: 2009-05-14 22:51:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/Channel.h,v $
                                                                        
                                                                        
#ifndef Channel_h
#define Channel_h

// Written: fmk
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class definition for Channel.
// Channel is an abstract base class which defines the channel interface.
// A channel is a point of communication in a program, a mailbox to/from
// which data enters/leaves a program.
//
// What: "@(#) Channel.h, revA"



class ChannelAddress;
class Message;
class MovableObject;
class Matrix;
class Vector;
class ID;
class FEM_ObjectBroker;

class Channel
{
  public:
    Channel ();    
    virtual ~Channel();

    // methods to set up the channel in an actor space
    virtual char *addToProgram(void) =0;
    virtual int setUpConnection(void) =0;
    virtual int setNextAddress(const ChannelAddress &theAddress) =0;
    virtual ChannelAddress *getLastSendersAddress(void) =0;

    virtual int isDatastore(void);
    virtual int getDbTag(void);
    int getTag(void);
    
    // methods to send/receive messages and objects on channels.
    virtual int sendObj(int commitTag,
			MovableObject &theObject, 
		    ChannelAddress *theAddress =0) =0;

    virtual int recvObj(int commitTag,
			MovableObject &theObject, 
			FEM_ObjectBroker &theBroker, 
		    ChannelAddress *theAddress =0) =0; 

    virtual int sendMsg(int dbTag, int commitTag, 
			const Message &theMessage, 
		    ChannelAddress *theAddress =0) =0;  

    virtual int recvMsg(int dbTag, int commitTag, 
			Message &theMessage, 
			ChannelAddress *theAddress =0) =0;  
    
    virtual int recvMsgUnknownSize(int dbTag, int commitTag, 
			Message &theMessage, 
			ChannelAddress *theAddress =0) =0;  

    virtual int sendMatrix(int dbTag, int commitTag, 
			const Matrix &theMatrix, 
			ChannelAddress *theAddress =0) =0;  

    virtual int recvMatrix(int dbTag, int commitTag, 
			Matrix &theMatrix, 
			ChannelAddress *theAddress =0) =0;  
    
    virtual int sendVector(int dbTag, int commitTag, 
			const Vector &theVector, 
			ChannelAddress *theAddress =0) =0;  

    virtual int recvVector(int dbTag, int commitTag, 
			Vector &theVector, 
			ChannelAddress *theAddress =0) =0;  
    
    virtual int sendID(int dbTag, int commitTag, 
		    const ID &theID, 
		    ChannelAddress *theAddress =0) =0;  

    virtual int recvID(int dbTag, int commitTag, 
		    ID &theID, 
		    ChannelAddress *theAddress =0) =0;      

  protected:
    
  private:
    static int numChannel;
    int tag;
};

#endif
