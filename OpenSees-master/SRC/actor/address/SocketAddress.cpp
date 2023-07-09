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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-10-15 00:32:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/address/SocketAddress.cpp,v $
                                                                        
                                                                        
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for SocketAddress.
// It is used to encapsulate the addresses used to send/recv messages
// using Berkeley sockets. SocketAddress is needed as a friend by 
// UDP_Socket & TCP_Socket.


#include <SocketAddress.h>



SocketAddress::SocketAddress()
:ChannelAddress(SOCKET_TYPE)
{

}

SocketAddress::SocketAddress(char *, unsigned int)
:ChannelAddress(SOCKET_TYPE)
{

}

SocketAddress::~SocketAddress()
{

}

