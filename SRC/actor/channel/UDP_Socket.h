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
// $Date: 2010-08-06 21:51:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/UDP_Socket.h,v $

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class definition for UDP_Socket.
// UDP_Socket is a sub-class of channel. It is implemented with
// Berkeley datagram sockets using the UDP protocol. Messages delivery
// is thus unreliable.

#ifndef UDP_Socket_h
#define UDP_Socket_h

#include <bool.h>
#include <Socket.h>
#include <Channel.h>

class UDP_Socket : public Channel
{
  public:
    UDP_Socket();
    UDP_Socket(unsigned int port, bool checkEndianness = false);   
    UDP_Socket(unsigned int other_Port, char *other_InetAddr,
        bool checkEndianness = false);  
    ~UDP_Socket();
    
    char *addToProgram();
    
    virtual int setUpConnection();

    int setNextAddress(const ChannelAddress &otherChannelAddress);
    virtual ChannelAddress *getLastSendersAddress(){ return 0;};

    int sendObj(int commitTag,
		MovableObject &theObject, 
		ChannelAddress *theAddress =0);
    int recvObj(int commitTag,
		MovableObject &theObject, 
		FEM_ObjectBroker &theBroker,
		ChannelAddress *theAddress =0);
		
    int sendMsg(int dbTag, int commitTag, 
		const Message &, 
		ChannelAddress *theAddress =0);    
    int recvMsg(int dbTag, int commitTag, 
		Message &, 
		ChannelAddress *theAddress =0);        
    int recvMsgUnknownSize(int dbTag, int commitTag, 
		Message &, 
		ChannelAddress *theAddress =0);        

    int sendMatrix(int dbTag, int commitTag, 
		   const Matrix &theMatrix, 
		   ChannelAddress *theAddress =0);
    int recvMatrix(int dbTag, int commitTag, 
		   Matrix &theMatrix, 
		   ChannelAddress *theAddress =0);
    
    int sendVector(int dbTag, int commitTag, 
		   const Vector &theVector, ChannelAddress *theAddress =0);
    int recvVector(int dbTag, int commitTag, 
		   Vector &theVector, 
		   ChannelAddress *theAddress =0);
    
    int sendID(int dbTag, int commitTag, 
	       const ID &theID, 
	       ChannelAddress *theAddress =0);
    int recvID(int dbTag, int commitTag, 
	       ID &theID, 
	       ChannelAddress *theAddress =0);    
    
  protected:
    unsigned int getPortNumber() const;
    unsigned int getBytesAvailable();
    
  private:
    socket_type sockfd;
    
    union {
	  struct sockaddr    addr;
	  struct sockaddr_in addr_in;
    } my_Addr;
    union {
      struct sockaddr    addr;
      struct sockaddr_in addr_in;
    } other_Addr;
    
    socklen_type addrLength;
    
    unsigned int myPort;
    int connectType;
    bool checkEndianness;
    bool endiannessProblem;
};

#endif 
