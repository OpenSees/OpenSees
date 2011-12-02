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
// $Date: 2005-11-23 23:43:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/TCP_Socket.cpp,v $
                                                                        
                                                                        
// Written: fmk
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementation of the methods needed
// to define the TCP_Socket class interface.
//
// What: "@(#) TCP_Socket.C, revA"

#include "TCP_Socket.h"
#include <string.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Message.h>
#include <ChannelAddress.h>
#include <MovableObject.h>

static int GetHostAddr(char *host, char *IntAddr);
static void inttoa(unsigned int no, char *string, int *cnt);

// TCP_Socket(unsigned int other_Port, char *other_InetAddr): 
// 	constructor to open a socket with my inet_addr and with a port number 
//	given by the OS. 

TCP_Socket::TCP_Socket()
  :myPort(0), connectType(0)
{
    // set up my_Addr 
    bzero((char *) &my_Addr, sizeof(my_Addr));    
    my_Addr.addr_in.sin_family = AF_INET;
    my_Addr.addr_in.sin_port = htons(0);

    #ifdef _WIN32
      my_Addr.addr_in.sin_addr.S_un.S_addr = htonl(INADDR_ANY);
    #else
      my_Addr.addr_in.sin_addr.s_addr = htonl(INADDR_ANY);
    #endif

    addrLength = sizeof(my_Addr.addr_in);
    
    // open a socket
    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
 	opserr << "TCP_Socket::TCP_Socket - could not open socket\n";
    }

    // bind local address to it
    if (bind(sockfd, (struct sockaddr *) &my_Addr.addr_in, sizeof(my_Addr.addr_in)) < 0) {
	opserr << "TCP_Socket::TCP_Socket - could not bind local address\n";
    }
    
    // get my_address info
    getsockname(sockfd, &my_Addr.addr, &addrLength);
    myPort = ntohs(my_Addr.addr_in.sin_port);
}    





// TCP_Socket(unsigned int port): 
//	constructor to open a socket with my inet_addr and with a port number port.

TCP_Socket::TCP_Socket(unsigned int port) 
  :myPort(0), connectType(0)
{
    // set up my_Addr.addr_in with address given by port and internet address of
    // machine on which the process that uses this routine is running.

    // set up my_Addr 
    bzero((char *) &my_Addr.addr_in, sizeof(my_Addr.addr_in));
    my_Addr.addr_in.sin_family = AF_INET;
    my_Addr.addr_in.sin_port = htons(port);

    #ifdef _WIN32
      my_Addr.addr_in.sin_addr.S_un.S_addr = htonl(INADDR_ANY);
    #else
      my_Addr.addr_in.sin_addr.s_addr = htonl(INADDR_ANY);
    #endif

    addrLength = sizeof(my_Addr.addr_in);

    // open a socket
    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
	opserr << "TCP_Socket::TCP_Socket - could not open socket\n";
    }
    
    // bind local address to it
    if (bind(sockfd,&my_Addr.addr, sizeof(my_Addr.addr_in)) < 0) {
	opserr << "TCP_Socket::TCP_Socket - could not bind local address\n";
    }    

    // get my_address info
    getsockname(sockfd, &my_Addr.addr, &addrLength);
    myPort = ntohs(my_Addr.addr_in.sin_port);    
}



// TCP_Socket(unsigned int other_Port, char *other_InetAddr): 
// 	constructor to open a socket with my inet_addr and with a port number 
//	given by the OS. Then to connect with a TCP_Socket whose address is
//	given by other_Port and other_InetAddr. 

TCP_Socket::TCP_Socket(unsigned int other_Port, char *other_InetAddr)
  :myPort(0), connectType(1)
{
    // set up remote address
    bzero((char *) &other_Addr.addr_in, sizeof(other_Addr.addr_in));
    other_Addr.addr_in.sin_family      = AF_INET;
    other_Addr.addr_in.sin_port        = htons(other_Port);

    #ifdef _WIN32
      other_Addr.addr_in.sin_addr.S_un.S_addr = inet_addr(other_InetAddr);
    #else
      other_Addr.addr_in.sin_addr.s_addr = inet_addr(other_InetAddr);
    #endif

    // set up my_Addr.addr_in 
    bzero((char *) &my_Addr.addr_in, sizeof(my_Addr.addr_in));    
    my_Addr.addr_in.sin_family = AF_INET;
    my_Addr.addr_in.sin_port = htons(0);

    #ifdef _WIN32
      my_Addr.addr_in.sin_addr.S_un.S_addr = htonl(INADDR_ANY);
    #else
      my_Addr.addr_in.sin_addr.s_addr = htonl(INADDR_ANY);
    #endif

    addrLength = sizeof(my_Addr.addr_in);
    
    // open a socket
    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
 	opserr << "TCP_Socket::TCP_Socket - could not open socket\n";
    }

    // bind local address to it
    if (bind(sockfd, (struct sockaddr *) &my_Addr.addr_in,sizeof(my_Addr.addr_in)) < 0) {
	opserr << "TCP_Socket::TCP_Socket - could not bind local address\n";
    }
    myPort = ntohs(my_Addr.addr_in.sin_port);    
}    


// ~TCP_Socket():
//	destructor

TCP_Socket::~TCP_Socket()
{
  #ifdef _WIN32
    closesocket(sockfd);
  #else
    close(sockfd);
  #endif
}



int 
TCP_Socket::setUpConnection(void)
{
  if (connectType == 1) {

    // now try to connect to socket with remote address.
    if (connect(sockfd, (struct sockaddr *) &other_Addr.addr_in, 
		sizeof(other_Addr.addr_in))< 0) {
      
	opserr << "TCP_Socket::TCP_Socket - could not connect\n";
	return -1;
    }
    // get my_address info
    getsockname(sockfd, &my_Addr.addr, &addrLength);
    

  } else {
    
    // wait for other process to contact me & set up connection
    socket_type newsockfd;
    listen(sockfd, 1);    
    newsockfd = accept(sockfd, (struct sockaddr *) &other_Addr.addr_in, &addrLength);
    
    if (newsockfd < 0) {
      opserr << "TCP_Socket::TCP_Socket - could not accept connection\n";
      return -1;
    }    
    
    // close old socket & reset sockfd
    // we can close as we are not going to wait for others to connect
    #ifdef _WIN32
      closesocket(sockfd);
    #else
      close(sockfd);
    #endif

    sockfd = newsockfd;
    
    // get my_address info
    getsockname(sockfd, &my_Addr.addr, &addrLength);
    myPort = ntohs(my_Addr.addr_in.sin_port);    
  }    

  // set socket so no delay    
  /*
    int optlen;
    optlen = 1;
    if ((setsockopt(sockfd,IPPROTO_TCP, TCP_NODELAY, 
    (char *) &optlen, sizeof(int))) < 0) { 
    opserr << "TCP_Socket::TCP_Socket - could not set TCP_NODELAY\n";
    }	  
  */
  /*
    int flag=sizeof(int);
    if ((getsockopt(sockfd,IPPROTO_TCP, TCP_NODELAY, 
    (char *) &optlen, &flag)) < 0) { 
    opserr << "TCP_Socket::TCP_Socket - could not set TCP_NODELAY\n";
    }	        
    opserr << "TCP_Socket::TCP_Socket - " << optlen << " flag " << flag <<  endln;
  */
    

  return 0;
}    

int
TCP_Socket::setNextAddress(const ChannelAddress &theAddress)
{	
    SocketAddress *theSocketAddress = 0;
    if (theAddress.getType() == SOCKET_TYPE) {
	    theSocketAddress = (SocketAddress *)(&theAddress);    
	    // check address is the only address a TCP_socket can send to
	    if (bcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr, 
		     theSocketAddress->addrLength) != 0) {
		
		opserr << "TCP_Socket::recvMsg() - a TCP_Socket ";
		opserr << "can only communicate with one other TCP_Socket\n"; 
		return -1;
	    }
	}
    else {
	opserr << "TCP_Socket::setNextAddress() - a TCP_Socket ";
	opserr << "can only communicate with a TCP_Socket";
	opserr << " address given is not of type SocketAddress\n"; 
	return -1;	    
    }		    	
	
    return 0;
}



int 
TCP_Socket::sendObj(int commitTag,
		    MovableObject &theObject, 
		    ChannelAddress *theAddress) 
{
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_Socket::sendObj() - a TCP_Socket ";
	    opserr << "can only communicate with a TCP_Socket";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		    
	    
	if (bcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr, 
		 theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_Socket::sendObj() - a TCP_Socket ";
	    opserr << "can only communicate with one other TCP_Socket";
	    opserr << " address given is not that address\n"; 
	    return -1;	    
	}	
    }    
    return theObject.sendSelf(commitTag, *this);
}

int 
TCP_Socket::recvObj(int commitTag,
		    MovableObject &theObject, 
		    FEM_ObjectBroker &theBroker, 
		    ChannelAddress *theAddress)
{
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_Socket::sendObj() - a TCP_Socket ";
	    opserr << "can only communicate with a TCP_Socket";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		    
	if (bcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_Socket::recvMsg() - a TCP_Socket ";
	    opserr << "can only communicate with one other TCP_Socket\n"; 
	    return -1;
	}
    }
    return theObject.recvSelf(commitTag, *this, theBroker);
}


// void Recv(Message &):
// 	Method to receive a message, also sets other_Addr.addr_in to that of sender

int 
TCP_Socket::recvMsg(int dbTag, int commitTag,
		    Message &msg, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_Socket::sendObj() - a TCP_Socket ";
	    opserr << "can only communicate with a TCP_Socket";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		    
	if (bcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_Socket::recvMsg() - a TCP_Socket ";
	    opserr << "can only communicate with one other TCP_Socket\n"; 
	    return -1;
	}
    }

    // if o.k. get a ponter to the data in the message and 
    // place the incoming data there
    int nleft,nread;
    char *gMsg;
    gMsg = msg.data;
    nleft = msg.length;

    while (nleft > 0) {
	nread = recv(sockfd,gMsg,nleft,0);
	nleft -= nread;
	gMsg +=  nread;
    }
    return 0;
}


// void Send(Message &):
// 	Method to send a message to an address given by other_Addr.addr_in.

int 
TCP_Socket::sendMsg(int dbTag, int commitTag,
		    const Message &msg, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_Socket::sendObj() - a TCP_Socket ";
	    opserr << "can only communicate with a TCP_Socket";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		    
	if (bcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_Socket::recvMsg() - a TCP_Socket ";
	    opserr << "can only communicate with one other TCP_Socket\n"; 
	    return -1;
	}
    }

    // if o.k. get a ponter to the data in the message and 
    // place the incoming data there
    int nwrite, nleft;    
    char *gMsg;
    gMsg = msg.data;
    nleft = msg.length;

    while (nleft > 0) {
	nwrite = send(sockfd,gMsg,nleft,0);
	nleft -= nwrite;

	gMsg +=  nwrite;
    }
    return 0;
}



int 
TCP_Socket::recvMatrix(int dbTag, int commitTag,
		       Matrix &theMatrix, ChannelAddress *theAddress)
		    
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_Socket::sendObj() - a TCP_Socket ";
	    opserr << "can only communicate with a TCP_Socket";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		    
	if (bcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_Socket::recvMatrix() - a TCP_Socket ";
	    opserr << "can only communicate with one other TCP_Socket\n"; 
	    return -1;
	}
    }

    // if o.k. get a ponter to the data in the Matrix and 
    // place the incoming data there
    int nleft,nread;
    double *data = theMatrix.data;
    char *gMsg = (char *)data;;
    nleft =  theMatrix.dataSize * sizeof(double);

    while (nleft > 0) {
	nread = recv(sockfd,gMsg,nleft,0);
	nleft -= nread;
	gMsg +=  nread;
    }
    return 0;
}


// void Send(Matrix &):
// 	Method to send a Matrix to an address given by other_Addr.addr_in.

int 
TCP_Socket::sendMatrix(int dbTag, int commitTag,
		       const Matrix &theMatrix, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_Socket::sendObj() - a TCP_Socket ";
	    opserr << "can only communicate with a TCP_Socket";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		        SocketAddress *theSocketAddress = 0;

	if (bcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_Socket::recvMatrix() - a TCP_Socket ";
	    opserr << "can only communicate with one other TCP_Socket\n"; 
	    return -1;
	}
    }

    // if o.k. get a ponter to the data in the Matrix and 
    // place the incoming data there
    int nwrite, nleft;    
    double *data = theMatrix.data;
    char *gMsg = (char *)data;
    nleft =  theMatrix.dataSize * sizeof(double);

    while (nleft > 0) {
	nwrite = send(sockfd,gMsg,nleft,0);
	nleft -= nwrite;
	
	gMsg +=  nwrite;
    }
    return 0;
}








int 
TCP_Socket::recvVector(int dbTag, int commitTag,
		       Vector &theVector, ChannelAddress *theAddress)
		    
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_Socket::sendObj() - a TCP_Socket ";
	    opserr << "can only communicate with a TCP_Socket";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		
	if (bcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_Socket::recvVector() - a TCP_Socket ";
	    opserr << "can only communicate with one other TCP_Socket\n"; 
	    return -1;
	}
    }

    // if o.k. get a ponter to the data in the Vector and 
    // place the incoming data there
    int nleft,nread;
    double *data = theVector.theData;
    char *gMsg = (char *)data;;
    nleft =  theVector.sz * sizeof(double);

    while (nleft > 0) {
	nread = recv(sockfd,gMsg,nleft,0);
	nleft -= nread;
	gMsg +=  nread;
    }
    return 0;
}


// void Send(Vector &):
// 	Method to send a Vector to an address given by other_Addr.addr_in.

int 
TCP_Socket::sendVector(int dbTag, int commitTag,
		       const Vector &theVector, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_Socket::sendObj() - a TCP_Socket ";
	    opserr << "can only communicate with a TCP_Socket";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		
	if (bcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_Socket::recvVector() - a TCP_Socket ";
	    opserr << "can only communicate with one other TCP_Socket\n"; 
	    return -1;
	}
    }

    // if o.k. get a ponter to the data in the Vector and 
    // place the incoming data there
    int nwrite, nleft;    
    double *data = theVector.theData;
    char *gMsg = (char *)data;
    nleft =  theVector.sz * sizeof(double);
    
    while (nleft > 0) {
	nwrite = send(sockfd,gMsg,nleft,0);
	nleft -= nwrite;
	gMsg +=  nwrite;
    }
    return 0;
}




int 
TCP_Socket::recvID(int dbTag, int commitTag,
		   ID &theID, ChannelAddress *theAddress)
		    
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_Socket::sendObj() - a TCP_Socket ";
	    opserr << "can only communicate with a TCP_Socket";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		
	if (bcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_Socket::recvID() - a TCP_Socket ";
	    opserr << "can only communicate with one other TCP_Socket\n"; 
	    return -1;
	}
    }

    // if o.k. get a ponter to the data in the ID and 
    // place the incoming data there
    int nleft,nread;
    int *data = theID.data;
    char *gMsg = (char *)data;;
    nleft =  theID.sz * sizeof(int);

    while (nleft > 0) {
	nread = recv(sockfd,gMsg,nleft,0);
	nleft -= nread;
	gMsg +=  nread;
    }
    return 0;
}


// void Send(ID &):
// 	Method to send a ID to an address given by other_Addr.addr_in.

int 
TCP_Socket::sendID(int dbTag, int commitTag,
		   const ID &theID, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_Socket::sendObj() - a TCP_Socket ";
	    opserr << "can only communicate with a TCP_Socket";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		
	if (bcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_Socket::recvID() - a TCP_Socket ";
	    opserr << "can only communicate with one other TCP_Socket\n"; 
	    return -1;
	}
    }

    // if o.k. get a ponter to the data in the ID and 
    // place the incoming data there
    int nwrite, nleft;    
    int *data = theID.data;
    char *gMsg = (char *)data;
    nleft =  theID.sz * sizeof(int);
    
    while (nleft > 0) {
	nwrite = send(sockfd,gMsg,nleft,0);
	nleft -= nwrite;
	gMsg +=  nwrite;
    }
    return 0;
}






unsigned int 
TCP_Socket::getPortNumber(void) const
{
    return myPort;
}


char *
TCP_Socket::addToProgram(void)
{
    char *tcp = " 1 ";

    char  my_InetAddr[MAX_INET_ADDR];
    char  myPortNum[8];
    char  me[30];
    unsigned int thePort = this->getPortNumber();
   /*
    char *me =(char *)malloc(30*sizeof(char));
    char *my_InetAddr=(char *)malloc(30*sizeof(char));
    char *myPortNum = (char *)malloc(30*sizeof(char));
    for (int i=0; i<30; i++) {
	me[i] = ' ';
	my_InetAddr[i] = ' ';
	myPortNum[i] = ' ';
    }
    */
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


// G e t H o s t A d d r
//     	GetHostAddr is a function to get the internet address of a host
// 	Takes machine name host & Returns 0 if o.k,  -1 if gethostbyname 
//	error, -2 otherwise. The internet address is returned in IntAddr

static int GetHostAddr(char *host, char *IntAddr)
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

    
/*
 *  i n t t o a
 *
 *  Function to convert int to ascii
 *  
 */

static void inttoa(unsigned int no, char *string, int *cnt) {
    if (no /10) {
        inttoa(no/10, string, cnt);
        *cnt = *cnt+1;
    }
    string[*cnt] = no % 10 + '0';
}

