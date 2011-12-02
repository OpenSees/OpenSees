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
// $Date: 2000-09-15 08:23:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/address/MPI_ChannelAddress.cpp,v $
                                                                        
                                                                        
// File: ~/actor/MPI_ChannelAddress.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for MPI_ChannelAddress.
// It is used to encapsulate the addresses used to send/recv messages
// using Berkeley sockets. MPI_ChannelAddress is needed as a friend by 
// UDP_Socket & TCP_Socket.

#include <MPI_ChannelAddress.h>

MPI_ChannelAddress::MPI_ChannelAddress(int other)
:ChannelAddress(MPI_TYPE), otherTag(other), otherComm(MPI_COMM_WORLD)
{

}

MPI_ChannelAddress::MPI_ChannelAddress(int other, MPI_Comm otherC)
:ChannelAddress(MPI_TYPE), otherTag(other), otherComm(otherC)
{

}

MPI_ChannelAddress::~MPI_ChannelAddress()
{

}

