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
// $Source: /usr/local/cvs/OpenSees/SRC/actor/address/ChannelAddress.cpp,v $
                                                                        
                                                                        
// File: ~/actor/ChannelAddress.h
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ChannelAddress.
// It is used to encapsulate the addresses used to send/recv messages
// using Channels. 


#include <ChannelAddress.h>




ChannelAddress::ChannelAddress(int type)
:theType(type)
{
}

ChannelAddress::~ChannelAddress()
{
}

int
ChannelAddress::getType(void) const
{
    return theType;
}

