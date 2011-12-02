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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-10-15 00:31:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/TCP_Socket.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class definition for TCP_Socket.
// TCP_Socket is a sub-class of channel. It is implemented with Berkeley
// stream sockets using the TCP protocol. Messages delivery is garaunteed. 
// Communication is full-duplex between a pair of connected sockets.
//
// What: "@(#) TCP_Socket.h, revA"

#ifndef TCP_Socket_h
#define TCP_Socket_h

#include <bool.h>
#include <Socket.h>
#include <Channel.h>
#include <SocketAddress.h>

class TCP_Socket : public Channel
{
  public:
    TCP_Socket();        
    TCP_Socket(unsigned int);    
    TCP_Socket(unsigned int other_Port, char *other_InetAddr); 
    ~TCP_Socket();

    char *addToProgram(void);
    
    virtual int setUpConnection(void);

    int setNextAddress(const ChannelAddress &otherChannelAddress);
    virtual ChannelAddress *getLastSendersAddress(void){ return 0;};

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
    unsigned int getPortNumber(void) const;
    
  private:
    int sockfd;

    union {
	  struct sockaddr    addr;
	  struct sockaddr_in addr_in;
    } my_Addr;
    union {
      struct sockaddr    addr;
      struct sockaddr_in addr_in;
    } other_Addr;

    socklen_t addrLength;

    unsigned int myPort;
    int connectType;

    char add[40];
};


#endif 






 

