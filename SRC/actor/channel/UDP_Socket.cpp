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

// $Revision: 1.7 $
// $Date: 2010-08-06 21:51:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/UDP_Socket.cpp,v $

// Written: fmk
// Created: 11/95
// Revision: A
//
// Purpose: This file contains the implementation of the methods needed
// to define the UDP_Socket class interface.

#include "UDP_Socket.h"
#include <string.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Message.h>
#include <ChannelAddress.h>
#include <MovableObject.h>
#include <SocketAddress.h>

static int GetHostAddr(char *host, char *IntAddr);
static void inttoa(unsigned int no, char *string, int *cnt);

#ifndef _WIN32
static void byte_swap(void *array, long long nArray,int size);
#endif


// UDP_Socket(): 
//	constructor to open a socket with my inet_addr and with a port number 
//	assigned by the OS from the available free port numbers.
UDP_Socket::UDP_Socket() 
    : myPort(0), connectType(0),
    checkEndianness(false), endiannessProblem(false)
{
    // initialize sockets
    startup_sockets();
    
    // set up my_Addr 
    bzero((char *) &my_Addr, sizeof(my_Addr));    
    my_Addr.addr_in.sin_family = AF_INET;
    my_Addr.addr_in.sin_port = htons(0);
    
#ifdef _WIN32
    my_Addr.addr_in.sin_addr.S_un.S_addr = htonl(INADDR_ANY);
#else
    my_Addr.addr_in.sin_addr.s_addr = htonl(INADDR_ANY);
#endif
    
    // open a socket
    if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
        opserr << "UDP_Socket::UDP_Socket() - could not open socket\n";
    }
    
    // bind local address to it
    if (bind(sockfd, &my_Addr.addr, sizeof(my_Addr.addr)) < 0) {
        opserr << "UDP_Socket::UDP_Socket() - could not bind local address\n";
    }
    
    // get my_address info
    addrLength = sizeof(my_Addr.addr);
    myPort = ntohs(my_Addr.addr_in.sin_port);
}


// UDP_Socket(unsigned int port): 
//	constructor to open a socket with my inet_addr and with a port number port.
UDP_Socket::UDP_Socket(unsigned int port, bool checkendianness) 
    : myPort(0), connectType(0),
    checkEndianness(checkendianness), endiannessProblem(false)
{
    // initialize sockets
    startup_sockets();
    
    // set up my_Addr.addr_in with address given by port and internet address
    // of machine on which the process that uses this routine is running.
    
    // set up my_Addr
    bzero((char *) &my_Addr, sizeof(my_Addr));
    my_Addr.addr_in.sin_family = AF_INET;
    my_Addr.addr_in.sin_port = htons(port);
    
#ifdef _WIN32
    my_Addr.addr_in.sin_addr.S_un.S_addr = htonl(INADDR_ANY);
#else
    my_Addr.addr_in.sin_addr.s_addr = htonl(INADDR_ANY);
#endif
    
    // open a socket
    if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
        opserr << "UDP_Socket::UDP_Socket() - could not open socket\n";
    }
    
    // bind local address to it
    if (bind(sockfd, &my_Addr.addr, sizeof(my_Addr.addr)) < 0) {
        opserr << "UDP_Socket::UDP_Socket() - could not bind local address\n";
    }
    
    // get my_address info
    addrLength = sizeof(my_Addr.addr);
    myPort = ntohs(my_Addr.addr_in.sin_port);
}


// UDP_Socket(unsigned int other_Port, char *other_InetAddr): 
// 	constructor to open a socket with my inet_addr and with a port number 
//	given by the OS. It then sends a messege comprising of a single character
// 	to an address given by other_InetAddr and other_Port. This is to allow
// 	a shadow object to find address of the actor it initiatites.
UDP_Socket::UDP_Socket(unsigned int other_Port,
    char *other_InetAddr, bool checkendianness) 
    : myPort(0), connectType(1),
    checkEndianness(checkendianness), endiannessProblem(false)
{
    // initialize sockets
    startup_sockets();
    
    // set up remote address
    bzero((char *) &other_Addr, sizeof(other_Addr));
    other_Addr.addr_in.sin_family = AF_INET;
    other_Addr.addr_in.sin_port = htons(other_Port);
    
#ifdef _WIN32
    other_Addr.addr_in.sin_addr.S_un.S_addr = inet_addr(other_InetAddr);
#else
    other_Addr.addr_in.sin_addr.s_addr = inet_addr(other_InetAddr);
#endif
    
    // set up my_Addr.addr_in 
    bzero((char *) &my_Addr, sizeof(my_Addr));
    my_Addr.addr_in.sin_family = AF_INET;
    my_Addr.addr_in.sin_port = htons(0);
    
#ifdef _WIN32
    my_Addr.addr_in.sin_addr.S_un.S_addr = htonl(INADDR_ANY);
#else
    my_Addr.addr_in.sin_addr.s_addr = htonl(INADDR_ANY);
#endif
    
    // open a socket
    if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
        opserr << "UDP_Socket::UDP_Socket() - could not open socket\n";
    }
    
    // bind local address to it
    if (bind(sockfd, &my_Addr.addr, sizeof(my_Addr.addr)) < 0) {
        opserr << "UDP_Socket::UDP_Socket() - could not bind local address\n";
    }
    
    addrLength = sizeof(my_Addr.addr);
    getsockname(sockfd, &my_Addr.addr, &addrLength);
    myPort = ntohs(my_Addr.addr_in.sin_port);    
}


// ~UDP_Socket():
//	destructor
UDP_Socket::~UDP_Socket()
{
#ifdef _WIN32
    closesocket(sockfd);
#else
    close(sockfd);
#endif
    
    // cleanup sockets
    cleanup_sockets();
}


int 
UDP_Socket::setUpConnection()
{
    if (connectType == 1) {
        
        // send a message to address
        char data = 'a';
        sendto(sockfd, &data, 1, 0, &other_Addr.addr, addrLength);
        
        // receive a message from other
        recvfrom(sockfd, &data, 1, 0, &other_Addr.addr, &addrLength);    
        
        // check for endianness problem if requested
        if (checkEndianness) {
            int i = 1;
            int j;
            
            int *data = &i;
            char *gMsg = (char *)data;
            sendto(sockfd, gMsg, 1, 0, &other_Addr.addr, addrLength);
            
            data = &j;
            gMsg = (char *)data;
            recvfrom(sockfd, gMsg, 1, 0, &other_Addr.addr, &addrLength);
            
            if (i != j) {
                int k = 0x41424344;
                char *c = (char *)&k;
                
                if (*c == 0x41) {
                    endiannessProblem = true;
                }
            }
        }
        
    } else {
        
        // wait for remote process to send message
        char data;
        recvfrom(sockfd, &data, 1, 0, &other_Addr.addr, &addrLength);    
        
        // then send a message back
        data = 'b';
        sendto(sockfd, &data, 1, 0, &other_Addr.addr, addrLength);        
        
        // check for endianness problem if requested
        if (checkEndianness) {
            int i;
            int j = 1;
            
            int *data = &i;
            char *gMsg = (char *)data;
            recvfrom(sockfd, gMsg, 1, 0, &other_Addr.addr, &addrLength);
            
            data = &j;
            gMsg = (char *)data;
            sendto(sockfd, gMsg, 1, 0, &other_Addr.addr, addrLength);
            
            if (i != j) {
                int k = 0x41424344;
                char *c = (char *)&k;
                if (*c == 0x41)
                    endiannessProblem = true;
            }
        }
    }

    return 0;    
}    


int
UDP_Socket::setNextAddress(const ChannelAddress &theAddress)
{	
    SocketAddress *theSocketAddress = 0;
    if (theAddress.getType() == SOCKET_TYPE) {
        theSocketAddress = (SocketAddress *)(&theAddress);    
        // set up the address of the Socket to which data will be sent
        bcopy((char *) &theSocketAddress->address.addr, (char *) &other_Addr, 
            theSocketAddress->addrLength);
        addrLength = theSocketAddress->addrLength;	

        return  0;	    
    }
    else {
        opserr << "UDP_Socket::setNextAddress() - a UDP_Socket ";
        opserr << "can only communicate with a UDP_Socket";
        opserr << " address given is not of type SocketAddress\n"; 
        return -1;	    
    }		    	
}


int 
UDP_Socket::sendObj(int commitTag,
    MovableObject &theObject, ChannelAddress *theAddress)
{
    // set up the address of the Socket to which data will be sent
    if (theAddress != 0) {
        SocketAddress *theSocketAddress = 0;
        
        if (theAddress->getType() == SOCKET_TYPE) {
            theSocketAddress = (SocketAddress *)theAddress;
            
            bcopy((char *) &theSocketAddress->address.addr, (char *) &other_Addr.addr, 
                theSocketAddress->addrLength);
            addrLength = theSocketAddress->addrLength;
        }
        else {
            opserr << "UDP_Socket::sendObj() - a UDP_Socket ";
            opserr << "can only communicate with a UDP_Socket";
            opserr << " address given is not of type SocketAddress\n";
            return -1;
        }
    }
    
    return theObject.sendSelf(commitTag, *this);
}


int 
UDP_Socket::recvObj(int commitTag,
    MovableObject &theObject, FEM_ObjectBroker &theBroker, 
    ChannelAddress *theAddress)
{
    int res = theObject.recvSelf(commitTag, *this, theBroker);
    if (res < 0) 
        return res;
    
    // check the address that message came from was correct
    if (theAddress != 0) {
        SocketAddress *theSocketAddress = 0;
        
        if (theAddress->getType() == SOCKET_TYPE) {
            theSocketAddress = (SocketAddress *)theAddress;
            
            if (memcmp((char *) &theSocketAddress->address.addr, (char *) &other_Addr.addr, 
                theSocketAddress->addrLength) != 0) {
                    opserr << "UDP_Socket::recvObj() - a UDP_Socket ";
                    opserr << "can only look at first incoming message\n"; 
                    opserr << "The last message did not come from write scource\n";
                    return -1;
            }
        }
        else {
            opserr << "UDP_Socket::recvObj() - a UDP_Socket ";
            opserr << "can only communicate with a UDP_Socket";
            opserr << " address given is not of type SocketAddress\n";
            return -1;
        }
    }
    
    return 0;
}


int
UDP_Socket::recvMsg(int dbTag, int commitTag,
    Message &msg, ChannelAddress *theAddress)
{	
    // if o.k. get a pointer to the data in the message and 
    // place the incoming data there
    int size;
    char *gMsg;
    gMsg = msg.data;
    size = msg.length;
    
    while (size > 0) {
        if (size <= MAX_UDP_DATAGRAM) {
            recvfrom(sockfd, gMsg, size, 0, &other_Addr.addr, &addrLength);
            size = 0;
        }
        else {
            recvfrom(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &other_Addr.addr, &addrLength);
            gMsg += MAX_UDP_DATAGRAM;
            size -= MAX_UDP_DATAGRAM;
        }
    }
    
    // check the address that message came from was correct
    if (theAddress != 0) {
        SocketAddress *theSocketAddress = 0;
        
        if (theAddress->getType() == SOCKET_TYPE) {
            theSocketAddress = (SocketAddress *)theAddress;
            
            if (memcmp((char *) &theSocketAddress->address.addr, (char *) &other_Addr.addr, 
                theSocketAddress->addrLength) != 0) {
                    opserr << "UDP_Socket::recvMsg() - a UDP_Socket ";
                    opserr << "can only look at first incoming message\n"; 
                    opserr << "The last message did not come from write scource\n";
                    return -1;
            }
        }
        else {
            opserr << "UDP_Socket::sendObj() - a UDP_Socket ";
            opserr << "can only communicate with a UDP_Socket";
            opserr << " address given is not of type SocketAddress\n";
            return -1;
        }
    }
    
    return 0;
}


int
UDP_Socket::recvMsgUnknownSize(int dbTag, int commitTag,
    Message &msg, ChannelAddress *theAddress)
{	
    opserr << "UDP_Socket::recvMsgUnknownSize() - not implemented yet.";
    return -1;	    
}

    
int
UDP_Socket::sendMsg(int dbTag, int commitTag,
    const Message &msg, ChannelAddress *theAddress)
{	
    // set up the address of the Socket to which data will be sent
    if (theAddress != 0) {
        SocketAddress *theSocketAddress = 0;
        
        if (theAddress->getType() == SOCKET_TYPE) {
            theSocketAddress = (SocketAddress *)theAddress;
            
            bcopy((char *) &theSocketAddress->address.addr, (char *) &other_Addr.addr, 
                theSocketAddress->addrLength);
            addrLength = theSocketAddress->addrLength;
        }
        else {
            opserr << "UDP_Socket::sendMsg() - a UDP_Socket ";
            opserr << "can only communicate with a UDP_Socket";
            opserr << " address given is not of type SocketAddress\n";
            return -1;
        }
    }
    
    // send the data
    int size; 
    char *gMsg;
    gMsg = msg.data;
    size = msg.length;
    
    while (size > 0) {
        if (size <= MAX_UDP_DATAGRAM) {
            sendto(sockfd, gMsg, size, 0, &other_Addr.addr, addrLength);
            size = 0;
        }
        else {	
            sendto(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &other_Addr.addr, addrLength);
            gMsg += MAX_UDP_DATAGRAM;
            size -= MAX_UDP_DATAGRAM;
        }
    }
    
    return 0;
}


int
UDP_Socket::recvMatrix(int dbTag, int commitTag,
    Matrix &theMatrix, ChannelAddress *theAddress)
{	
    // if o.k. get a pointer to the data in the message and 
    // place the incoming data there
    int size;
    double *data = theMatrix.data;
    char *gMsg = (char *)data;;
    size = theMatrix.dataSize * sizeof(double);
    
    while (size > 0) {
        if (size <= MAX_UDP_DATAGRAM) {
            recvfrom(sockfd, gMsg, size, 0, &other_Addr.addr, &addrLength);
            size = 0;
        }
        else {
            recvfrom(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &other_Addr.addr, &addrLength);
            gMsg += MAX_UDP_DATAGRAM;
            size -= MAX_UDP_DATAGRAM;
        }
    }
    
#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theMatrix.dataSize, sizeof(double));
    }
#endif
    
    // check the address that message came from was correct
    if (theAddress != 0) {
        SocketAddress *theSocketAddress = 0;
        
        if (theAddress->getType() == SOCKET_TYPE) {
            theSocketAddress = (SocketAddress *)theAddress;	
            
            if (memcmp((char *) &theSocketAddress->address.addr, (char *) &other_Addr.addr, 
                theSocketAddress->addrLength) != 0) {
                    opserr << "UDP_Socket::recvMsg() - a UDP_Socket ";
                    opserr << "can only look at first incoming message\n"; 
                    opserr << "The last message did not come from write scource\n";
                    return -1;
            }
        }
        else {
            opserr << "UDP_Socket::recvMatrix() - a UDP_Socket ";
            opserr << "can only communicate with a UDP_Socket";
            opserr << " address given is not of type SocketAddress\n";
            return -1;
        }
    }
    
    return 0;
}


int
UDP_Socket::sendMatrix(int dbTag, int commitTag,
    const Matrix &theMatrix, ChannelAddress *theAddress)
{	
    // set up the address of the Socket to which data will be sent
    if (theAddress != 0) {
        SocketAddress *theSocketAddress = 0;
        
        if (theAddress->getType() == SOCKET_TYPE) {
            theSocketAddress = (SocketAddress *)theAddress;
            
            bcopy((char *) &theSocketAddress->address.addr, (char *) &other_Addr.addr, 
                theSocketAddress->addrLength);
            addrLength = theSocketAddress->addrLength;	
        }
        else {
            opserr << "UDP_Socket::sendMatrix() - a UDP_Socket ";
            opserr << "can only communicate with a UDP_Socket";
            opserr << " address given is not of type SocketAddress\n";
            return -1;
        }
    }
    
    // send the data
    int size; 
    double *data = theMatrix.data;
    char *gMsg = (char *)data;;
    size = theMatrix.dataSize * sizeof(double);
    
#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theMatrix.dataSize,  sizeof(double));
    }
#endif
    
    while (size > 0) {
        if (size <= MAX_UDP_DATAGRAM) {
            sendto(sockfd, gMsg, size, 0, &other_Addr.addr, addrLength);
            size = 0;
        }
        else {
            sendto(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &other_Addr.addr, addrLength);
            gMsg += MAX_UDP_DATAGRAM;
            size -= MAX_UDP_DATAGRAM;
        }
    }
    
#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theMatrix.dataSize,  sizeof(double));
    }
#endif
    
    return 0;
}


int
UDP_Socket::recvVector(int dbTag, int commitTag,
    Vector &theVector, ChannelAddress *theAddress)
{	
    // if o.k. get a pointer to the data in the message and 
    // place the incoming data there
    int size;
    double *data = theVector.theData;
    char *gMsg = (char *)data;;
    size = theVector.sz * sizeof(double);
    
    while (size > 0) {
        if (size <= MAX_UDP_DATAGRAM) {
            recvfrom(sockfd, gMsg, size, 0, &other_Addr.addr, &addrLength);
            size = 0;
        }
        else {
            recvfrom(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &other_Addr.addr, &addrLength);
            gMsg += MAX_UDP_DATAGRAM;
            size -= MAX_UDP_DATAGRAM;
        }
    }
    
#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theVector.sz,  sizeof(double));
    }
#endif
    
    // check the address that message came from was correct
    if (theAddress != 0) {
        SocketAddress *theSocketAddress = 0;
        
        if (theAddress->getType() == SOCKET_TYPE) {
            theSocketAddress = (SocketAddress *)theAddress;
            
            if (memcmp((char *) &theSocketAddress->address.addr, (char *) &other_Addr.addr, 
                theSocketAddress->addrLength) != 0) {
                    opserr << "UDP_Socket::recvMsg() - a UDP_Socket ";
                    opserr << "can only look at first incoming message\n";
                    opserr << "The last message did not come from write scource\n";
                    return -1;
            }
        }
        else {
            opserr << "UDP_Socket::recvVector() - a UDP_Socket ";
            opserr << "can only communicate with a UDP_Socket";
            opserr << " address given is not of type SocketAddress\n";
            return -1;
        }
    }
    
    return 0;
}


int
UDP_Socket::sendVector(int dbTag, int commitTag,
    const Vector &theVector, ChannelAddress *theAddress)
{	
    // set up the address of the Socket to which data will be sent
    if (theAddress != 0) {
        SocketAddress *theSocketAddress = 0;
        
        if (theAddress->getType() == SOCKET_TYPE) {
            theSocketAddress = (SocketAddress *)theAddress;
            
            bcopy((char *) &theSocketAddress->address.addr, (char *) &other_Addr.addr, 
                theSocketAddress->addrLength);
            addrLength = theSocketAddress->addrLength;	
        }
        else {
            opserr << "UDP_Socket::sendVector() - a UDP_Socket ";
            opserr << "can only communicate with a UDP_Socket";
            opserr << " address given is not of type SocketAddress\n";
            return -1;
        }
    }
    
    // send the data
    int size; 
    double *data = theVector.theData;
    char *gMsg = (char *)data;;
    size = theVector.sz * sizeof(double);
    
#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theVector.sz,  sizeof(double));
    }
#endif
    
    while (size > 0) {
        if (size <= MAX_UDP_DATAGRAM) {
            sendto(sockfd, gMsg, size, 0, &other_Addr.addr, addrLength);
            size = 0;
        }
        else {
            sendto(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &other_Addr.addr, addrLength);
            gMsg += MAX_UDP_DATAGRAM;
            size -= MAX_UDP_DATAGRAM;
        }
    }
    
#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theVector.sz,  sizeof(double));
    }
#endif
    
    return 0;
}


int
UDP_Socket::recvID(int dbTag, int commitTag,
    ID &theID, ChannelAddress *theAddress)
{	
    // if o.k. get a pointer to the data in the message and 
    // place the incoming data there
    int size;
    int *data = theID.data;
    char *gMsg = (char *)data;;
    size = theID.sz * sizeof(int);
    
    while (size > 0) {
        if (size <= MAX_UDP_DATAGRAM) {
            recvfrom(sockfd, gMsg, size, 0, &other_Addr.addr, &addrLength);
            size = 0;
        }
        else {
            recvfrom(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &other_Addr.addr, &addrLength);
            gMsg += MAX_UDP_DATAGRAM;
            size -= MAX_UDP_DATAGRAM;
        }
    }
    
#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theID.sz, sizeof(int));
    }
#endif
    
    // check the address that message came from was correct
    if (theAddress != 0) {
        SocketAddress *theSocketAddress = 0;
        
        if (theAddress->getType() == SOCKET_TYPE) {
            theSocketAddress = (SocketAddress *)theAddress;
            
            if (memcmp((char *) &theSocketAddress->address.addr, (char *) &other_Addr.addr, 
                theSocketAddress->addrLength) != 0) {
                    opserr << "UDP_Socket::recvMsg() - a UDP_Socket ";
                    opserr << "can only look at first incoming message\n";
                    opserr << "The last message did not come from write scource\n";
                    return -1;
            }
        }
        else {
            opserr << "UDP_Socket::recvID() - a UDP_Socket ";
            opserr << "can only communicate with a UDP_Socket";
            opserr << " address given is not of type SocketAddress\n";
            return -1;
        }
    }
    
    return 0;
}


int
UDP_Socket::sendID(int dbTag, int commitTag,
    const ID &theID, ChannelAddress *theAddress)
{	
    // set up the address of the Socket to which data will be sent
    if (theAddress != 0) {
        SocketAddress *theSocketAddress = 0;
        
        if (theAddress->getType() == SOCKET_TYPE) {
            theSocketAddress = (SocketAddress *)theAddress;
            
            bcopy((char *) &theSocketAddress->address.addr, (char *) &other_Addr.addr, 
                theSocketAddress->addrLength);
            addrLength = theSocketAddress->addrLength;
        }
        else {
            opserr << "UDP_Socket::sendID() - a UDP_Socket ";
            opserr << "can only communicate with a UDP_Socket";
            opserr << " address given is not of type SocketAddress\n";
            return -1;
        }
    }
    
    // send the data
    int size; 
    int *data = theID.data;
    char *gMsg = (char *)data;;
    size = theID.sz * sizeof(int);
    
#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theID.sz,  sizeof(int));
    }
#endif
    
    while (size > 0) {
        if (size <= MAX_UDP_DATAGRAM) {
            sendto(sockfd, gMsg, size, 0, &other_Addr.addr, addrLength);
            size = 0;
        }
        else {
            sendto(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &other_Addr.addr, addrLength);
            gMsg += MAX_UDP_DATAGRAM;
            size -= MAX_UDP_DATAGRAM;
        }
    }
    
#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theID.sz,  sizeof(int));
    }
#endif
    
    return 0;
}


unsigned int 
UDP_Socket::getPortNumber() const
{
    return myPort;
}


unsigned int
UDP_Socket::getBytesAvailable()
{
    unsigned long bytesAvailable;

#ifdef _WIN32
    ioctlsocket(sockfd,FIONREAD,&bytesAvailable);
#else
    ioctl(sockfd,FIONREAD,&bytesAvailable);
#endif

    return bytesAvailable;
}


char *
UDP_Socket::addToProgram()
{
    const char *tcp = " 2 ";
    
    char  my_InetAddr[MAX_INET_ADDR];
    char  myPortNum[8];
    char  me[30];
    unsigned int thePort = this->getPortNumber();
    
    int start = 0;
    inttoa(thePort,myPortNum,&start);
    gethostname(me,MAX_INET_ADDR);
    GetHostAddr(me,my_InetAddr);
    
    char *newStuff =(char *)malloc(100*sizeof(char));
    for (int i=0; i<100; i++) 
        newStuff[i] = ' ';
    
    strcpy(newStuff,tcp);
    strcat(newStuff," ");          
    strcat(newStuff,my_InetAddr);
    strcat(newStuff," ");
    strcat(newStuff,myPortNum);
    strcat(newStuff," ");    
    
    return newStuff;
}


static int
GetHostAddr(char *host, char *IntAddr)
{
    register struct hostent *hostptr;

    if ( (hostptr = gethostbyname(host)) == NULL) 
        return (-1);

    switch(hostptr->h_addrtype) {
      case AF_INET:
          strcpy(IntAddr,inet_ntoa(*(struct in_addr *)*hostptr->h_addr_list));
          return (0);

      default:
          return (-2);
    }
}


static void
inttoa(unsigned int no, char *string, int *cnt) {
    if (no /10) {
        inttoa(no/10, string, cnt);
        *cnt = *cnt+1;
    }
    string[*cnt] = no % 10 + '0';
}


#ifndef _WIN32
static void
byte_swap(void *array, long long nArray,int size)
{

    long long i;
    int j;
    unsigned char *p= (unsigned char *) array;
    int half = size/2;
    unsigned char temp;
    unsigned char *out;

    if(size < 2)
        return;

    for(i=0; i < nArray;i++) {
        out = p + size -1;
        for (j = 0; j < half; ++j)
        {
            temp = *out;
            *out = *p;
            *p = temp;
            /*process next byte*/
            ++p;
            --out;
        }
        p += half;
    }
}
#endif
