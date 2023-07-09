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
// $Date: 2007-07-16 22:56:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/address/SocketAddress.h,v $
                                                                        
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for SocketAddress.
// It is used to encapsulate the addresses used to send/recv messages
// using Berkeley sockets. SocketAddress is needed as a friend by 
// UDP_Socket & TCP_Socket.


#ifndef SocketAddress_h
#define SocketAddress_h

#include "Socket.h"
#include <ChannelAddress.h>

class SocketAddress: public ChannelAddress
{
  public:
    SocketAddress();
    SocketAddress(char *, unsigned int);
    virtual ~SocketAddress();

  private:
    union {
      struct sockaddr_in  addr_in;
      struct sockaddr  addr;      
    } address;
    int addrLength;    

    friend class UDP_Socket;
    friend class TCP_Socket;
    friend class TCP_SocketSSL;
    friend class TCP_SocketNoDelay;
};



#endif
