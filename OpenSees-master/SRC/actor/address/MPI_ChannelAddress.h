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
// $Source: /usr/local/cvs/OpenSees/SRC/actor/address/MPI_ChannelAddress.h,v $
                                                                        
                                                                        
// File: ~/actor/MPI_ChannelAddress.h
//
// Written: fmk 06/98
// Revised:
//
// Purpose: This file contains the class definition for MPI_ChannelAddress.
// It is used to encapsulate the addresses used to send/recv messages
// using the MPI library.


#ifndef MPI_ChannelAddress_h
#define MPI_ChannelAddress_h

#ifdef _PARALLEL_PROCESSING
#define MPI_TYPE 2
#endif

#ifdef _PARALLEL_INTERPRETERS
#define MPI_TYPE 2
#endif



#include <ChannelAddress.h>
#include <mpi.h>

class MPI_ChannelAddress: public ChannelAddress
{
  public:
    MPI_ChannelAddress(int otherProcessTag);
    MPI_ChannelAddress(int otherProcessTag, MPI_Comm otherComm);
    virtual ~MPI_ChannelAddress();

    friend class MPI_Channel;

  private:
    int otherTag;
    MPI_Comm otherComm;

};

#endif

