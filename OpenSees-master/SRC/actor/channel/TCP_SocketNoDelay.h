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
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/TCP_SocketNoDelay.h,v $
                                                                        
                                                                        
// File: ~/actor/TCP_SocketNoDelay.h
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for TCP_SocketNoDelay.
// TCP_SocketNoDelay is a sub-class of channel. It is implemented with Berkeley
// stream sockets using the TCP protocol. Messages delivery is garaunteed. 
// Communication is full-duplex between a pair of connected sockets.

#ifndef TCP_SocketNoDelay_h
#define TCP_SocketNoDelay_h

#include <Socket.h>
#include <Channel.h>
#include <SocketAddress.h>

class TCP_SocketNoDelay : public Channel
{
  public:
    TCP_SocketNoDelay();        
    TCP_SocketNoDelay(unsigned int);    
    TCP_SocketNoDelay(unsigned int other_Port, char *other_InetAddr); 
    ~TCP_SocketNoDelay();

    virtual int setUpShadow(void);
    virtual int setUpActor(void);

    int setNextAddress(const ChannelAddress &otherChannelAddress);
    virtual ChannelAddress *getLastSendersAddress(void){ return 0;};

    int sendObj(MovableObject &theObject, 
		FEM_ObjectBroker &theBroker,
		ChannelAddress *theAddress =0);
    
    int recvObj(MovableObject &theObject, 
		FEM_ObjectBroker &theBroker,
		ChannelAddress *theAddress =0);
		
    int recvMsg(Message &, ChannelAddress *theAddress =0);       
    int sendMsg(const Message &, ChannelAddress *theAddress =0);    

    int sendMatrix(const Matrix &theMatrix, ChannelAddress *theAddress =0);
    int recvMatrix(Matrix &theMatrix, ChannelAddress *theAddress =0);
    
    int sendVector(const Vector &theVector, ChannelAddress *theAddress =0);
    int recvVector(Vector &theVector, ChannelAddress *theAddress =0);
    
    int sendID(const ID &theID, ChannelAddress *theAddress =0);
    int recvID(ID &theID, ChannelAddress *theAddress =0);    
    
    unsigned int getPortNumber(void) const;
    char *addToProgram(void);
    
  protected:
	
  private:
    int sockfd;
    struct sockaddr_in  my_Addr, other_Addr;
    int addrLength;
    unsigned int myPort;
    char add[40];
};


#endif 






 

