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
                                                                        
// $Revision: 1.1 $
// $Date: 2010-02-25 20:46:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/TCP_Stream.cpp,v $

// Written: fmk

#include <TCP_Stream.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <ID.h>
#include <Channel.h>
#include <Message.h>

#include <TCP_Socket.h>

TCP_Stream::TCP_Stream()
  :OPS_Stream(OPS_STREAM_TAGS_TCP_Stream), sendSize(0), data(1), theChannel(0) 
{
  theChannel = new TCP_Socket();
}


TCP_Stream::TCP_Stream(unsigned int other_Port, 
		       const char *other_InetAddr,
		       bool checkEndianness)
  :OPS_Stream(OPS_STREAM_TAGS_TCP_Stream), sendSize(0), data(1), theChannel(0) 
{
  theChannel = new TCP_Socket(other_Port, other_InetAddr, checkEndianness);
  if (theChannel->setUpConnection() < 0) {
    opserr << "TCP_Stream - Failed to set up connection\n";
    delete theChannel;
    theChannel = 0;
  }
}


TCP_Stream::~TCP_Stream()
{
  data(0) = -1;
  if (theChannel != 0 && theChannel->sendVector(0,0, data) < 0)
    opserr << "TCP_Stream - failed to send close signal\n";

  if (theChannel != 0) 
    delete theChannel;
}

int 
TCP_Stream::setFile(const char *name, openMode mode)
{
  return 0;
}

int 
TCP_Stream::open(void)
{
  return 0;
}

int 
TCP_Stream::close(void)
{
  return 0;
}

int 
TCP_Stream::tag(const char *tagName)
{
  return 0;
}

int 
TCP_Stream::tag(const char *tagName, const char *value)
{
  return 0;
}


int 
TCP_Stream::endTag()
{
  return 0;
}

int 
TCP_Stream::attr(const char *name, int value)
{
  return 0;
}

int 
TCP_Stream::attr(const char *name, double value)
{
  return 0;
}

int 
TCP_Stream::attr(const char *name, const char *value)
{
  return 0;
}

int 
TCP_Stream::write(Vector &dataToSend)
{
  int sizeToSend = dataToSend.Size();
  if (sizeToSend == 0 || theChannel == 0)
    return 0;

  if (sizeToSend != sendSize) {
    data(0) = sizeToSend;
    if (theChannel->sendVector(0, 0, data) < 0) {
      opserr << "TCP_Stream - failed to send new size\n";
      return -1;
    }

    data.resize(sizeToSend+1);
    sendSize = sizeToSend;
    data(0) = sendSize;
  }

  for (int i=0, j=1; i<sendSize; i++, j++)
    data(j) = dataToSend(i);

  if (theChannel->sendVector(0,0, data) < 0) {
      opserr << "TCP_Stream - failed to send data\n";
      return -1;
  }
    
  return 0;
}

OPS_Stream& 
TCP_Stream::write(const char *s,int n)
{
  return *this;
}

OPS_Stream& 
TCP_Stream::write(const unsigned char*s,int n)
{
  return *this;
}
OPS_Stream& 
TCP_Stream::write(const signed char*s,int n)
{
  return *this;
}

OPS_Stream& 
TCP_Stream::write(const void *s, int n)
{
  return *this;
}

OPS_Stream& 
TCP_Stream::operator<<(char c)
{
  return *this;
}

OPS_Stream& 
TCP_Stream::operator<<(unsigned char c)
{  
  return *this;
}

OPS_Stream& 
TCP_Stream::operator<<(signed char c)
{
  return *this;
}

OPS_Stream& 
TCP_Stream::operator<<(const char *s)
{
  return *this;
}

OPS_Stream& 
TCP_Stream::operator<<(const unsigned char *s)
{
  return *this;
}

OPS_Stream& 
TCP_Stream::operator<<(const signed char *s)
{
  return *this;
}

OPS_Stream& 
TCP_Stream::operator<<(const void *p)
{
  return *this;
}

OPS_Stream& 
TCP_Stream::operator<<(int n)
{
  return *this;
}
OPS_Stream& 
TCP_Stream::operator<<(unsigned int n)
{  
  return *this;
}
OPS_Stream& 
TCP_Stream::operator<<(long n)
{
  return *this;
}
OPS_Stream& 
TCP_Stream::operator<<(unsigned long n)
{
  return *this;
}
OPS_Stream& 
TCP_Stream::operator<<(short n)
{
  return *this;
}
OPS_Stream& 
TCP_Stream::operator<<(unsigned short n)
{
  return *this;
}
OPS_Stream& 
TCP_Stream::operator<<(bool b)
{
  return *this;
}
OPS_Stream& 
TCP_Stream::operator<<(double n)
{
  return *this;
}
OPS_Stream& 
TCP_Stream::operator<<(float n)
{
  return *this;
}


int 
TCP_Stream::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
TCP_Stream::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}

