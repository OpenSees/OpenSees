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

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 03/23
// Revision: A
//
// Description: This file contains the implementation of the GenericClient class.

#include <SocketStream.h>

#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <ID.h>
#include <Channel.h>
#include <Message.h>
#include <TCP_Socket.h>
#include <UDP_Socket.h>


SocketStream::SocketStream()
    : OPS_Stream(OPS_STREAM_TAGS_SocketStream),
    sendSize(0), data(1), theChannel(0)
{
    theChannel = new TCP_Socket();
}


SocketStream::SocketStream(
    unsigned int other_Port,
    const char* other_InetAddr,
    bool udp,
    bool checkEndianness)
    : OPS_Stream(OPS_STREAM_TAGS_SocketStream), 
    sendSize(0), data(1), theChannel(0)
{
    // setup the connection
    if (udp) {
        if (other_InetAddr == 0)
            theChannel = new UDP_Socket(other_Port, "127.0.0.1", checkEndianness);
        else
            theChannel = new UDP_Socket(other_Port, other_InetAddr, checkEndianness);
    }
    else {
        if (other_InetAddr == 0)
            theChannel = new TCP_Socket(other_Port, "127.0.0.1", checkEndianness);
        else
            theChannel = new TCP_Socket(other_Port, other_InetAddr, checkEndianness);
    }
    if (!theChannel) {
        opserr << "SocketStream::SocketStream() "
            << "- failed to create channel\n";
        exit(-1);
    }
    if (theChannel->setUpConnection() != 0) {
        opserr << "SocketStream::SocketStream() "
            << "- failed to setup connection\n";
        delete theChannel;
        theChannel = 0;
    }
}


SocketStream::~SocketStream()
{
    data(0) = -1;
    if (theChannel != 0 && theChannel->sendVector(0, 0, data) < 0)
        opserr << "SocketStream::~SocketStream() "
        << "- failed to send close signal\n";
    
    if (theChannel != 0)
        delete theChannel;
}


/*int SocketStream::setFile(const char* name, openMode mode)
{
    return 0;
}*/


int SocketStream::open(void)
{
    return 0;
}


int SocketStream::close(void)
{
    return 0;
}


int SocketStream::tag(const char* tagName)
{
    return 0;
}


int SocketStream::tag(const char* tagName, const char* value)
{
    return 0;
}


int SocketStream::endTag()
{
    return 0;
}


int SocketStream::attr(const char* name, int value)
{
    return 0;
}


int SocketStream::attr(const char* name, double value)
{
    return 0;
}


int SocketStream::attr(const char* name, const char* value)
{
    return 0;
}


int SocketStream::write(Vector& dataToSend)
{
    int sizeToSend = dataToSend.Size();
    if (sizeToSend == 0 || theChannel == 0)
        return 0;
    
    if (sizeToSend != sendSize) {
        data(0) = sizeToSend;
        if (theChannel->sendVector(0, 0, data) < 0) {
            opserr << "SocketStream::write() "
                << "- failed to send new size\n";
            return -1;
        }
        data.resize(sizeToSend + 1);
        sendSize = sizeToSend;
        data(0) = sendSize;
    }
    
    for (int i = 0, j = 1; i < sendSize; i++, j++)
        data(j) = dataToSend(i);
    
    if (theChannel->sendVector(0, 0, data) < 0) {
        opserr << "SocketStream::write() "
            << "- failed to send data\n";
        return -2;
    }
    
    return 0;
}


OPS_Stream& SocketStream::write(const char* s, int n)
{
    return *this;
}


OPS_Stream& SocketStream::write(const unsigned char* s, int n)
{
    return *this;
}


OPS_Stream& SocketStream::write(const signed char* s, int n)
{
    return *this;
}


OPS_Stream& SocketStream::write(const void* s, int n)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(char c)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(unsigned char c)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(signed char c)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(const char* s)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(const unsigned char* s)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(const signed char* s)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(const void* p)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(int n)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(unsigned int n)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(long n)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(unsigned long n)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(short n)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(unsigned short n)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(bool b)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(double n)
{
    return *this;
}


OPS_Stream& SocketStream::operator<<(float n)
{
    return *this;
}


int SocketStream::sendSelf(int commitTag, Channel& theChannel)
{
    return -1;
}

int SocketStream::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    return -1;
}
