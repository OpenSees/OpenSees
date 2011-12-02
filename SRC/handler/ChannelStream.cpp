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
// $Source: /usr/local/cvs/OpenSees/SRC/handler/ChannelStream.cpp,v $

// Written:fmk

#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <ID.h>
#include <Channel.h>
#include <Message.h>

#include <ChannelStream.h>

ChannelStream::ChannelStream(Channel *theC)
  :OPS_Stream(OPS_STREAM_TAGS_ChannelStream), sendSize(0), data(1), theChannel(theC) 
{
  if (theChannel->setUpConnection() < 0) {
    opserr << "ChannelStream - Failed to set up connection\n";
    delete theChannel;
    theChannel = 0;
  }
}

ChannelStream::~ChannelStream()
{
  data(0) = -1;
  if (theChannel != 0 && theChannel->sendVector(0,0, data) < 0)
    opserr << "ChannelStream - failed to send close signal\n";

  if (theChannel != 0) 
    delete theChannel;
}

int 
ChannelStream::setFile(const char *name, openMode mode)
{
  return 0;
}

int 
ChannelStream::open(void)
{
  return 0;
}

int 
ChannelStream::close(void)
{
  return 0;
}

int 
ChannelStream::tag(const char *tagName)
{
  return 0;
}

int 
ChannelStream::tag(const char *tagName, const char *value)
{
  return 0;
}


int 
ChannelStream::endTag()
{
  return 0;
}

int 
ChannelStream::attr(const char *name, int value)
{
  return 0;
}

int 
ChannelStream::attr(const char *name, double value)
{
  return 0;
}

int 
ChannelStream::attr(const char *name, const char *value)
{
  return 0;
}

int 
ChannelStream::write(Vector &dataToSend)
{
  int sizeToSend = dataToSend.Size();
  if (sizeToSend == 0 || theChannel == 0)
    return 0;

  if (sizeToSend != sendSize) {
    data(0) = sizeToSend;
    if (theChannel->sendVector(0, 0, data) < 0) {
      opserr << "ChannelStream - failed to send new size\n";
      return -1;
    }

    data.resize(sizeToSend+1);
    sendSize = sizeToSend;
    data(0) = sendSize;
  }

  for (int i=0, j=1; i<sendSize; i++, j++)
    data(j) = dataToSend(i);

  if (theChannel->sendVector(0,0, data) < 0) {
      opserr << "ChannelStream - failed to send data\n";
      return -1;
  }
    
  return 0;
}

OPS_Stream& 
ChannelStream::write(const char *s,int n)
{
  return *this;
}

OPS_Stream& 
ChannelStream::write(const unsigned char*s,int n)
{
  return *this;
}
OPS_Stream& 
ChannelStream::write(const signed char*s,int n)
{
  return *this;
}

OPS_Stream& 
ChannelStream::write(const void *s, int n)
{
  return *this;
}

OPS_Stream& 
ChannelStream::operator<<(char c)
{
  return *this;
}

OPS_Stream& 
ChannelStream::operator<<(unsigned char c)
{  
  return *this;
}

OPS_Stream& 
ChannelStream::operator<<(signed char c)
{
  return *this;
}

OPS_Stream& 
ChannelStream::operator<<(const char *s)
{
  return *this;
}

OPS_Stream& 
ChannelStream::operator<<(const unsigned char *s)
{
  return *this;
}

OPS_Stream& 
ChannelStream::operator<<(const signed char *s)
{
  return *this;
}

OPS_Stream& 
ChannelStream::operator<<(const void *p)
{
  return *this;
}

OPS_Stream& 
ChannelStream::operator<<(int n)
{
  return *this;
}
OPS_Stream& 
ChannelStream::operator<<(unsigned int n)
{  
  return *this;
}
OPS_Stream& 
ChannelStream::operator<<(long n)
{
  return *this;
}
OPS_Stream& 
ChannelStream::operator<<(unsigned long n)
{
  return *this;
}
OPS_Stream& 
ChannelStream::operator<<(short n)
{
  return *this;
}
OPS_Stream& 
ChannelStream::operator<<(unsigned short n)
{
  return *this;
}
OPS_Stream& 
ChannelStream::operator<<(bool b)
{
  return *this;
}
OPS_Stream& 
ChannelStream::operator<<(double n)
{
  return *this;
}
OPS_Stream& 
ChannelStream::operator<<(float n)
{
  return *this;
}


int 
ChannelStream::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
ChannelStream::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}

