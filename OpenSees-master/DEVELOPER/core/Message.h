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
// $Date: 2007-07-16 22:56:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/message/Message.h,v $
                                                                        
                                                                        
// File: ~/actor/Message.h
//
// Written: fmk 11/96
// Revised:
//
// Purpose: This file contains the class definition for Message.

#ifndef Message_h
#define Message_h

class Message
{
  public:
    Message();
    Message(double *, int);
    Message(int *, int);
    Message(char *, int);
    virtual ~Message();

	virtual void setData(char *theData, int length);
    virtual int putData(char *theData, int startLoc, int endLoc);    
    virtual const char *getData(void);
    virtual int getSize(void);

    friend class UDP_Socket;
    friend class TCP_Socket;
    friend class TCP_SocketSSL;
    friend class TCP_SocketNoDelay;
    friend class MPI_Channel;
    
  private:
    int length;
    char *data;    
};

#endif
