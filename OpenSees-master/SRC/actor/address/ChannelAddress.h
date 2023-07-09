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
// $Date: 2008-04-14 17:28:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/address/ChannelAddress.h,v $
                                                                        
                                                                        
// File: ~/actor/ChannelAddress.h
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ChannelAddress.
// It is used to encapsulate the addresses used to send/recv messages
// using Channels. 


#ifndef ChannelAddress_h
#define ChannelAddress_h

#define SOCKET_TYPE 1

class ChannelAddress
{
  public:
    ChannelAddress(int type);
    virtual ~ChannelAddress();

    virtual int getType(void) const;

  private:
    int theType;
};


#endif
