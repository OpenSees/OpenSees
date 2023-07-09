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
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/MPI_Channel.h,v $
                                                                        
                                                                        
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for MPI_Channel.
// MPI_Channel is a sub-class of channel. It is implemented with Berkeley
// stream sockets using the TCP protocol. Messages delivery is garaunteed. 
// Communication is full-duplex between a pair of connected sockets.

#ifndef MPI_Channel_h
#define MPI_Channel_h

#include <mpi.h>
#include <Channel.h>

class MPI_Channel : public Channel
{
  public:
    MPI_Channel(int otherProcess);        
    ~MPI_Channel();

    char *addToProgram(void);
    
    virtual int setUpConnection(void);

    int setNextAddress(const ChannelAddress &otherChannelAddress);
    virtual ChannelAddress *getLastSendersAddress(void);

    int sendObj(int commitTag,
		MovableObject &theObject, 
		ChannelAddress *theAddress =0);
    
    int recvObj(int commitTag,
		MovableObject &theObject, 
		FEM_ObjectBroker &theBroker,
		ChannelAddress *theAddress =0);

    int recvMsg(int dbTag, int commitTag, Message &, ChannelAddress *theAddress =0);       
    int sendMsg(int dbTag, int commitTag, const Message &, ChannelAddress *theAddress =0);    
    int recvMsgUnknownSize(int dbTag, int commitTag, 
			   Message &, 
			   ChannelAddress *theAddress =0);        

    int sendMatrix(int dbTag, int commitTag, const Matrix &theMatrix, ChannelAddress *theAddress =0);
    int recvMatrix(int dbTag, int commitTag, Matrix &theMatrix, ChannelAddress *theAddress =0);
    
    int sendVector(int dbTag, int commitTag, const Vector &theVector, ChannelAddress *theAddress =0);
    int recvVector(int dbTag, int commitTag, Vector &theVector, ChannelAddress *theAddress =0);
    
    int sendID(int dbTag, int commitTag, const ID &theID, ChannelAddress *theAddress =0);
    int recvID(int dbTag, int commitTag, ID &theID, ChannelAddress *theAddress =0);    
    
    
  protected:
	
  private:
    int otherTag;
    MPI_Comm otherComm;    
};


#endif 






 

