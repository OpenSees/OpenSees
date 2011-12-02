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
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/MPI_Channel.cpp,v $
                                                                        
                                                                        
// File: ~/actor/channel/MPI_Channel.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the implementation of the methods needed
// to define the MPI_Channel class interface.

#include <MPI_Channel.h>
#include <Matrix.h>
#include <ID.h>
#include <Vector.h>
#include <Message.h>
#include <MPI_ChannelAddress.h>
#include <MovableObject.h>

// MPI_Channel(unsigned int other_Port, char *other_InetAddr): 
// 	constructor to open a socket with my inet_addr and with a port number 
//	given by the OS. 

MPI_Channel::MPI_Channel(int other)
 :otherTag(other), otherComm(MPI_COMM_WORLD)
{
  
}    

// ~MPI_Channel():
//	destructor

MPI_Channel::~MPI_Channel()
{

}


int 
MPI_Channel::setUpActor(void) 
{
  return 0;    
}


int 
MPI_Channel::setUpShadow(void) 
{
    return 0;
}    

ChannelAddress *
MPI_Channel::getLastSendersAddress(void) 
{
    cerr << "MPI_Channel::getLastSendersAddress(void) - ";
    cerr << " this should not be called - need MPI-2.0 to use\n";

    return 0;
}    


int
MPI_Channel::setNextAddress(const ChannelAddress &theAddress)
{	
    MPI_ChannelAddress *theMPI_ChannelAddress = 0;
    if (theAddress.getType() == MPI_TYPE) {
      theMPI_ChannelAddress = (MPI_ChannelAddress *)(&theAddress);    
      otherTag = theMPI_ChannelAddress->otherTag;
      otherComm= theMPI_ChannelAddress->otherComm;
    }
    else {
	cerr << "MPI_Channel::setNextAddress() - an MPI_Channel ";
	cerr << "can only communicate with an MPI_Channel";
	cerr << " address given is not of type MPI_ChannelAddress\n"; 
	return -1;	    
    }		    	
	
    return 0;
}



int 
MPI_Channel::sendObj(int commitTag,
		     MovableObject &theObject, 
		    ChannelAddress *theAddress) 
{
    // first check address is the only address a MPI_Channel can send to
    MPI_ChannelAddress *theMPI_ChannelAddress = 0;
    if (theAddress != 0) {
      if (theAddress->getType() == MPI_TYPE) {
	theMPI_ChannelAddress = (MPI_ChannelAddress *)theAddress;
	otherTag = theMPI_ChannelAddress->otherTag;
	otherComm= theMPI_ChannelAddress->otherComm;
      }
	else {
	    cerr << "MPI_Channel::sendObj() - a MPI_Channel ";
	    cerr << "can only communicate with a MPI_Channel";
	    cerr << " address given is not of type MPI_ChannelAddress\n"; 
	    return -1;	    
	}		    
    }    
    return theObject.sendSelf(commitTag, *this);
}

int 
MPI_Channel::recvObj(int commitTag,
		     MovableObject &theObject, 
		     FEM_ObjectBroker &theBroker, 
		     ChannelAddress *theAddress)
{
    // first check address is the only address a MPI_Channel can send to
    MPI_ChannelAddress *theMPI_ChannelAddress = 0;
    if (theAddress != 0) {
      if (theAddress->getType() == MPI_TYPE) {
	theMPI_ChannelAddress = (MPI_ChannelAddress *)theAddress;
	otherTag = theMPI_ChannelAddress->otherTag;
	otherComm= theMPI_ChannelAddress->otherComm;
      } else {
	cerr << "MPI_Channel::recvObj() - a MPI_Channel ";
	cerr << "can only communicate with a MPI_Channel";
	cerr << " address given is not of type MPI_ChannelAddress\n"; 
	return -1;	    
      }		    
    }
    return theObject.recvSelf(commitTag, *this, theBroker);
}


// void Recv(Message &):
// 	Method to receive a message, also sets other_Addr to that of sender

int 
MPI_Channel::recvMsg(int dbTag, int commitTag, Message &msg, ChannelAddress *theAddress)
{	
    // first check address is the only address a MPI_Channel can send to
    MPI_ChannelAddress *theMPI_ChannelAddress = 0;
    if (theAddress != 0) {
      if (theAddress->getType() == MPI_TYPE) {
	theMPI_ChannelAddress = (MPI_ChannelAddress *)theAddress;
	otherTag = theMPI_ChannelAddress->otherTag;
	otherComm= theMPI_ChannelAddress->otherComm;
      } else {
	cerr << "MPI_Channel::recvMesg() - a MPI_Channel ";
	cerr << "can only communicate with a MPI_Channel";
	cerr << " address given is not of type MPI_ChannelAddress\n"; 
	return -1;	    
      }		    
    }

    // if o.k. get a ponter to the data in the message and 
    // place the incoming data there
    int nleft,nread;
    char *gMsg;
    gMsg = msg.data;
    nleft = msg.length;

    MPI_Status status;
    MPI_Recv((void *)gMsg, nleft, MPI_CHAR, otherTag, 0, otherComm, &status);
    int count =0;
    MPI_Get_count(&status, MPI_CHAR, &count);
    if (count != nleft) {
      cerr << "MPI_Channel::recvMesg() -";
      cerr << " incorrect size of Message received ";
      return -1;
    }
    else
      return 0;
}


// void Send(Message &):
// 	Method to send a message to an address given by other_Addr.

int 
MPI_Channel::sendMsg(int dbTag, int commitTag, const Message &msg, ChannelAddress *theAddress)
{	
    // first check address is the only address a MPI_Channel can send to
    MPI_ChannelAddress *theMPI_ChannelAddress = 0;
    if (theAddress != 0) {
      if (theAddress->getType() == MPI_TYPE) {
	theMPI_ChannelAddress = (MPI_ChannelAddress *)theAddress;
	otherTag = theMPI_ChannelAddress->otherTag;
	otherComm= theMPI_ChannelAddress->otherComm;
      } else {
	cerr << "MPI_Channel::sendMsg() - a MPI_Channel ";
	cerr << "can only communicate with a MPI_Channel";
	cerr << " address given is not of type MPI_ChannelAddress\n"; 
	return -1;	    
      }		    
    }

    // if o.k. get a ponter to the data in the message and 
    // place the incoming data there
    int nwrite, nleft;    
    char *gMsg;
    gMsg = msg.data;
    nleft = msg.length;

    MPI_Send((void *)gMsg, nleft, MPI_CHAR, otherTag, 0, otherComm);
    return 0;
}

int 
MPI_Channel::recvMatrix(int dbTag, int commitTag, Matrix &theMatrix, ChannelAddress *theAddress)
		    
{	
    // first check address is the only address a MPI_Channel can send to
    MPI_ChannelAddress *theMPI_ChannelAddress = 0;
    if (theAddress != 0) {
      if (theAddress->getType() == MPI_TYPE) {
	theMPI_ChannelAddress = (MPI_ChannelAddress *)theAddress;
	otherTag = theMPI_ChannelAddress->otherTag;
	otherComm= theMPI_ChannelAddress->otherComm;
      } else {
	cerr << "MPI_Channel::recvMatrix() - a MPI_Channel ";
	cerr << "can only communicate with a MPI_Channel";
	cerr << " address given is not of type MPI_ChannelAddress\n"; 
	return -1;	    
      }		    
    }

    // if o.k. get a ponter to the data in the Matrix and 
    // place the incoming data there
    int nleft,nread;
    double *data = theMatrix.data;
    char *gMsg = (char *)data;;
    nleft =  theMatrix.dataSize;

    MPI_Status status;
    MPI_Recv((void *)gMsg, nleft, MPI_DOUBLE, otherTag, 0, 
	     otherComm, &status);
    int count = 0;
    MPI_Get_count(&status, MPI_DOUBLE, &count);
    if (count != nleft) {
      cerr << "MPI_Channel::recvMatrix() -";
      cerr << " incorrect number of entries for Matrix received ";
      return -1;
    }
    else
      return 0;
}


// void Send(Matrix &):
// 	Method to send a Matrix to an address given by other_Addr.

int 
MPI_Channel::sendMatrix(int dbTag, int commitTag, const Matrix &theMatrix, ChannelAddress *theAddress)
{	
    // first check address is the only address a MPI_Channel can send to
    MPI_ChannelAddress *theMPI_ChannelAddress = 0;
    if (theAddress != 0) {
      if (theAddress->getType() == MPI_TYPE) {
	theMPI_ChannelAddress = (MPI_ChannelAddress *)theAddress;
	otherTag = theMPI_ChannelAddress->otherTag;
	otherComm= theMPI_ChannelAddress->otherComm;
      } else {
	cerr << "MPI_Channel::sendMatrix() - a MPI_Channel ";
	cerr << "can only communicate with a MPI_Channel";
	cerr << " address given is not of type MPI_ChannelAddress\n"; 
	return -1;	    
      }		    
    }

    // if o.k. get a ponter to the data in the Matrix and 
    // place the incoming data there
    int nwrite, nleft;    
    double *data = theMatrix.data;
    char *gMsg = (char *)data;
    nleft =  theMatrix.dataSize;

    MPI_Send((void *)gMsg, nleft, MPI_DOUBLE, otherTag, 0, otherComm);

    return 0;
}








int 
MPI_Channel::recvVector(int dbTag, int commitTag, Vector &theVector, ChannelAddress *theAddress)
		    
{	
    // first check address is the only address a MPI_Channel can send to
    MPI_ChannelAddress *theMPI_ChannelAddress = 0;
    if (theAddress != 0) {
      if (theAddress->getType() == MPI_TYPE) {
	theMPI_ChannelAddress = (MPI_ChannelAddress *)theAddress;
	otherTag = theMPI_ChannelAddress->otherTag;
	otherComm= theMPI_ChannelAddress->otherComm;
      } else {
	cerr << "MPI_Channel::recvVector() - a MPI_Channel ";
	cerr << "can only communicate with a MPI_Channel";
	cerr << " address given is not of type MPI_ChannelAddress\n"; 
	return -1;	    
      }		    
    }

    // if o.k. get a ponter to the data in the Vector and 
    // place the incoming data there
    int nleft,nread;
    double *data = theVector.theData;
    char *gMsg = (char *)data;;
    nleft =  theVector.sz;

    MPI_Status status;
    MPI_Recv((void *)gMsg, nleft, MPI_DOUBLE, otherTag, 0, otherComm, &status);
    int count =0;
    MPI_Get_count(&status, MPI_DOUBLE, &count);
    if (count != nleft) {
      cerr << "MPI_Channel::recvVector() -";
      cerr << " incorrect number of entries for Vector received ";
      return -1;
    }
    else
      return 0;
}


// void Send(Vector &):
// 	Method to send a Vector to an address given by other_Addr.

int 
MPI_Channel::sendVector(int dbTag, int commitTag, const Vector &theVector, ChannelAddress *theAddress)
{	
    // first check address is the only address a MPI_Channel can send to
    MPI_ChannelAddress *theMPI_ChannelAddress = 0;
    if (theAddress != 0) {
      if (theAddress->getType() == MPI_TYPE) {
	theMPI_ChannelAddress = (MPI_ChannelAddress *)theAddress;
	otherTag = theMPI_ChannelAddress->otherTag;
	otherComm= theMPI_ChannelAddress->otherComm;
      } else {
	cerr << "MPI_Channel::sendVector() - a MPI_Channel ";
	cerr << "can only communicate with a MPI_Channel";
	cerr << " address given is not of type MPI_ChannelAddress\n"; 
	return -1;	    
      }		    
    }

    // if o.k. get a ponter to the data in the Vector and 
    // place the incoming data there
    int nwrite, nleft;    
    double *data = theVector.theData;
    char *gMsg = (char *)data;
    nleft =  theVector.sz;

    MPI_Send((void *)gMsg, nleft, MPI_DOUBLE, otherTag, 0, otherComm);
    
    return 0;
}




int 
MPI_Channel::recvID(int dbTag, int commitTag, ID &theID, ChannelAddress *theAddress)
{	
    // first check address is the only address a MPI_Channel can send to
    MPI_ChannelAddress *theMPI_ChannelAddress = 0;
    if (theAddress != 0) {
      if (theAddress->getType() == MPI_TYPE) {
	theMPI_ChannelAddress = (MPI_ChannelAddress *)theAddress;
	otherTag = theMPI_ChannelAddress->otherTag;
	otherComm= theMPI_ChannelAddress->otherComm;
      } else {
	cerr << "MPI_Channel::recvID() - a MPI_Channel ";
	cerr << "can only communicate with a MPI_Channel";
	cerr << " address given is not of type MPI_ChannelAddress\n"; 
	return -1;	    
      }		    
    }

    // if o.k. get a ponter to the data in the ID and 
    // place the incoming data there
    int nleft,nread;
    int *data = theID.data;
    char *gMsg = (char *)data;;
    nleft =  theID.sz;

    MPI_Status status;
    MPI_Recv((void *)gMsg, nleft, MPI_INT, otherTag, 0, otherComm, &status);
    int count =0;
    MPI_Get_count(&status, MPI_INT, &count);
    if (count != nleft) {
      cerr << "MPI_Channel::recvID() -";
      cerr << " incorrect number of entries for ID received ";
      return -1;
    }
    else
      return 0;
}


// void Send(ID &):
// 	Method to send a ID to an address given by other_Addr.

int 
MPI_Channel::sendID(int dbTag, int commitTag, const ID &theID, ChannelAddress *theAddress)
{	
    // first check address is the only address a MPI_Channel can send to
    MPI_ChannelAddress *theMPI_ChannelAddress = 0;
    if (theAddress != 0) {
      if (theAddress->getType() == MPI_TYPE) {
	theMPI_ChannelAddress = (MPI_ChannelAddress *)theAddress;
	otherTag = theMPI_ChannelAddress->otherTag;
	otherComm= theMPI_ChannelAddress->otherComm;
      } else {
	cerr << "MPI_Channel::sendID() - a MPI_Channel ";
	cerr << "can only communicate with a MPI_Channel";
	cerr << " address given is not of type MPI_ChannelAddress\n"; 
	return -1;	    
      }		    
    }

    // if o.k. get a ponter to the data in the ID and 
    // place the incoming data there
    int nwrite, nleft;    
    int *data = theID.data;
    char *gMsg = (char *)data;
    nleft =  theID.sz;

    MPI_Send((void *)gMsg, nleft, MPI_INT, otherTag, 0, otherComm);
    
    return 0;
}


/*
int 
MPI_Channel::getPortNumber(void) const
{
    return otherTag;
}
*/

char *
MPI_Channel::addToProgram(void)
{
    cerr << "MPI_Channel::addToProgram(void) - ";
    cerr << " this should not be called - need MPI-2.0\n";
    char *newStuff =(char *)malloc(10*sizeof(char));
    for (int i=0; i<10; i++) 
	newStuff[i] = ' ';

    return newStuff;
}


