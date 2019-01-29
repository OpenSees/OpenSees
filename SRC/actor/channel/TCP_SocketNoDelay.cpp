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
// $Date: 2003-02-14 23:00:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/TCP_SocketNoDelay.cpp,v $
                                                                        
                                                                        
// File: ~/actor/TCP_SocketNoDelay.C

//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the implementation of the methods needed
// to define the TCP_SocketNoDelay class interface.



#include "TCP_SocketNoDelay.h"
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <Matrix.h>
#include <Vector.h>
#include <Message.h>
#include <ChannelAddress.h>
#include <MovableObject.h>

static int GetHostAddr(char *host, char *IntAddr);
static void inttoa(unsigned int no, char *string, int *cnt);

// TCP_SocketNoDelay(unsigned int other_Port, char *other_InetAddr): 
// 	constructor to open a socket with my inet_addr and with a port number 
//	given by the OS. 

TCP_SocketNoDelay::TCP_SocketNoDelay()
 :myPort(0)
{
    // set up my_Addr 
    bzero((char *) &my_Addr, sizeof(my_Addr));    
    my_Addr.sin_family = AF_INET;
    my_Addr.sin_addr.s_addr = htonl(INADDR_ANY);
    my_Addr.sin_port = htons(0);
    addrLength = sizeof(my_Addr);
    
    // open a socket
    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
 	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not open socket\n";
    }

    // bind local address to it
    if (bind(sockfd, (struct sockaddr *) &my_Addr,sizeof(my_Addr)) < 0) {
	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not bind local address\n";
    }
    
    // get my_address info
    INET_getsockname(sockfd, &my_Addr, &addrLength);
    myPort = ntohs(my_Addr.sin_port);

}    





// TCP_SocketNoDelay(unsigned int port): 
//	constructor to open a socket with my inet_addr and with a port number port.

TCP_SocketNoDelay::TCP_SocketNoDelay(unsigned int port) 
:myPort(0)
{
    // set up my_Addr with address given by port and internet address of
    // machine on which the process that uses this routine is running.

    char  me[20];
    char  my_InetAddr[MAX_INET_ADDR];
    gethostname(me,MAX_INET_ADDR);
    GetHostAddr(me,my_InetAddr);

    bzero((char *) &my_Addr, sizeof(my_Addr));
    my_Addr.sin_family = AF_INET;
    my_Addr.sin_addr.s_addr = inet_addr(my_InetAddr);
    my_Addr.sin_port = htons(port);
    addrLength = sizeof(my_Addr);

    // open a socket

    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not open socket\n";
    }
    
    // bind local address to it

    if (bind(sockfd,(struct sockaddr *)&my_Addr,sizeof(my_Addr)) < 0) {
	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not bind local address\n";
    }    

    // get my_address info
    INET_getsockname(sockfd, &my_Addr, &addrLength);
    myPort = ntohs(my_Addr.sin_port);    
}



// TCP_SocketNoDelay(unsigned int other_Port, char *other_InetAddr): 
// 	constructor to open a socket with my inet_addr and with a port number 
//	given by the OS. Then to connect with a TCP_SocketNoDelay whose address is
//	given by other_Port and other_InetAddr. 

TCP_SocketNoDelay::TCP_SocketNoDelay(unsigned int other_Port, char *other_InetAddr)
:myPort(0)
{
    // set up remote address
    bzero((char *) &other_Addr, sizeof(other_Addr));
    other_Addr.sin_family      = AF_INET;
    other_Addr.sin_addr.s_addr = inet_addr(other_InetAddr);
    other_Addr.sin_port        = htons(other_Port);

    // set up my_Addr 
    bzero((char *) &my_Addr, sizeof(my_Addr));    
    my_Addr.sin_family = AF_INET;
    my_Addr.sin_addr.s_addr = htonl(INADDR_ANY);
    my_Addr.sin_port = htons(0);
    addrLength = sizeof(my_Addr);
    
    // open a socket
    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
 	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not open socket\n";
    }

    // bind local address to it
    if (bind(sockfd, (struct sockaddr *) &my_Addr,sizeof(my_Addr)) < 0) {
	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not bind local address\n";
    }
    myPort = ntohs(my_Addr.sin_port);    
}    


// ~TCP_SocketNoDelay():
//	destructor

TCP_SocketNoDelay::~TCP_SocketNoDelay()
{
    close(sockfd);
}



int 
TCP_SocketNoDelay::setUpActor(void) 
{
    // now try to connect to socket with remote address.
    if (connect(sockfd, (struct sockaddr *) &other_Addr, 
		sizeof(other_Addr))< 0) {
	
	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not connect\n";
	return -1;
    }
    // get my_address info
    INET_getsockname(sockfd, &my_Addr, &addrLength);


    // set socket so no delay    

    int optlen;
    optlen = 1;
    if ((setsockopt(sockfd,IPPROTO_TCP, TCP_NODELAY, 
		    (char *) &optlen, sizeof(int))) < 0) { 
	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not set TCP_NODELAY\n";
    }	  

    /*
    int flag=sizeof(int);
    if ((getsockopt(sockfd,IPPROTO_TCP, TCP_NODELAY, 
		    (char *) &optlen, &flag)) < 0) { 
	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not set TCP_NODELAY\n";
    }	        
    opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - " << optlen << " flag " << flag <<  endln;
    */
    
    return 0;
    
}


int 
TCP_SocketNoDelay::setUpShadow(void) 
{
    // wait for other process to contact me & set up connection
    int newsockfd;
    listen(sockfd, 1);    
    newsockfd = accept(sockfd, (struct sockaddr *) &other_Addr, &addrLength);

    if (newsockfd < 0) {
	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not accept connection\n";
	return -1;
    }    
    
    // close old socket & reset sockfd
    close(sockfd);  // we can close as we are not 
                    // going to wait for others to connect
    sockfd = newsockfd;

    // get my_address info
    INET_getsockname(sockfd, &my_Addr, &addrLength);
    myPort = ntohs(my_Addr.sin_port);    


    // set socket so no delay    
    int optlen;
    optlen = 1;
    if ((setsockopt(sockfd,IPPROTO_TCP, TCP_NODELAY, 
		    (char *) &optlen, sizeof(int))) < 0) { 
	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not set TCP_NODELAY\n";
    }	  

    /*
    int flag=sizeof(int);
    if ((getsockopt(sockfd,IPPROTO_TCP, TCP_NODELAY, 
		    (char *) &optlen, &flag)) < 0) { 
	opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - could not set TCP_NODELAY\n";
    }	        
    opserr << "TCP_SocketNoDelay::TCP_SocketNoDelay - " << optlen << " flag " << flag <<  endln;
    */
    
    return 0;
}    



int
TCP_SocketNoDelay::setNextAddress(const ChannelAddress &theAddress)
{	
    SocketAddress *theSocketAddress = 0;
    if (theAddress.getType() == SOCKET_TYPE) {
	    theSocketAddress = (SocketAddress *)(&theAddress);    
	    // check address is the only address a TCP_SocketNoDelay can send to
	    if (memcmp((char *) &other_Addr, (char *) &theSocketAddress->addr, 
		     theSocketAddress->addrLength) != 0) {
		
		opserr << "TCP_SocketNoDelay::recvMsg() - a TCP_SocketNoDelay ";
		opserr << "can only communicate with one other TCP_SocketNoDelay\n"; 
		return -1;
	    }
	}
    else {
	opserr << "TCP_SocketNoDelay::setNextAddress() - a TCP_SocketNoDelay ";
	opserr << "can only communicate with a TCP_SocketNoDelay";
	opserr << " address given is not of type SocketAddress\n"; 
	return -1;	    
    }		    	
	
    return 0;
}



int 
TCP_SocketNoDelay::sendObj(MovableObject &theObject, 
		    FEM_ObjectBroker &theBroker, 
		    ChannelAddress *theAddress) 
{
    // first check address is the only address a TCP_SocketNoDelay can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_SocketNoDelay::sendObj() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with a TCP_SocketNoDelay";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		    
	    
	if (memcmp((char *) &other_Addr, (char *) &theSocketAddress->addr, 
		 theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_SocketNoDelay::sendObj() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with one other TCP_SocketNoDelay";
	    opserr << " address given is not that address\n"; 
	    return -1;	    
	}	
    }    
    return theObject.sendSelf(*this, theBroker);
}

int 
TCP_SocketNoDelay::recvObj(MovableObject &theObject, 
		    FEM_ObjectBroker &theBroker, 
		    ChannelAddress *theAddress)
{
    // first check address is the only address a TCP_SocketNoDelay can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_SocketNoDelay::sendObj() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with a TCP_SocketNoDelay";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		    
	if (memcmp((char *) &other_Addr, (char *) &theSocketAddress->addr, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_SocketNoDelay::recvMsg() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with one other TCP_SocketNoDelay\n"; 
	    return -1;
	}
    }
    return theObject.recvSelf(*this, theBroker);
}


// void Recv(Message &):
// 	Method to receive a message, also sets other_Addr to that of sender

int 
TCP_SocketNoDelay::recvMsg(Message &msg, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_SocketNoDelay can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_SocketNoDelay::sendObj() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with a TCP_SocketNoDelay";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		    
	if (memcmp((char *) &other_Addr, (char *) &theSocketAddress->addr, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_SocketNoDelay::recvMsg() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with one other TCP_SocketNoDelay\n"; 
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
	nread = read(sockfd,gMsg,nleft);
	nleft -= nread;
	gMsg +=  nread;
    }
    return 0;
}


// void Send(Message &):
// 	Method to send a message to an address given by other_Addr.

int 
TCP_SocketNoDelay::sendMsg(const Message &msg, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_SocketNoDelay can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_SocketNoDelay::sendObj() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with a TCP_SocketNoDelay";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		    
	if (memcmp((char *) &other_Addr, (char *) &theSocketAddress->addr, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_SocketNoDelay::recvMsg() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with one other TCP_SocketNoDelay\n"; 
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
	nwrite = write(sockfd,gMsg,nleft);
	nleft -= nwrite;

	gMsg +=  nwrite;
    }
    return 0;
}



int 
TCP_SocketNoDelay::recvMatrix(Matrix &theMatrix, ChannelAddress *theAddress)
		    
{	
    // first check address is the only address a TCP_SocketNoDelay can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_SocketNoDelay::sendObj() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with a TCP_SocketNoDelay";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		    
	if (memcmp((char *) &other_Addr, (char *) &theSocketAddress->addr, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_SocketNoDelay::recvMatrix() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with one other TCP_SocketNoDelay\n"; 
	    return -1;
	}
    }

    // if o.k. get a ponter to the data in the Matrix and 
    // place the incoming data there
    int nleft,nread;
    double *data = theMatrix.myData;
    char *gMsg = (char *)data;;
    nleft =  theMatrix.dataSize * sizeof(double);

    while (nleft > 0) {
	nread = read(sockfd,gMsg,nleft);
	nleft -= nread;
	gMsg +=  nread;
    }
    return 0;
}


// void Send(Matrix &):
// 	Method to send a Matrix to an address given by other_Addr.

int 
TCP_SocketNoDelay::sendMatrix(const Matrix &theMatrix, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_SocketNoDelay can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_SocketNoDelay::sendObj() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with a TCP_SocketNoDelay";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		        SocketAddress *theSocketAddress = 0;

	if (memcmp((char *) &other_Addr, (char *) &theSocketAddress->addr, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_SocketNoDelay::recvMatrix() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with one other TCP_SocketNoDelay\n"; 
	    return -1;
	}
    }

    // if o.k. get a ponter to the data in the Matrix and 
    // place the incoming data there
    int nwrite, nleft;    
    double *data = theMatrix.myData;
    char *gMsg = (char *)data;
    nleft =  theMatrix.dataSize * sizeof(double);

    while (nleft > 0) {
	nwrite = write(sockfd,gMsg,nleft);
	nleft -= nwrite;
	
	gMsg +=  nwrite;
    }
    return 0;
}








int 
TCP_SocketNoDelay::recvVector(Vector &theVector, ChannelAddress *theAddress)
		    
{	
    // first check address is the only address a TCP_SocketNoDelay can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_SocketNoDelay::sendObj() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with a TCP_SocketNoDelay";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		
	if (memcmp((char *) &other_Addr, (char *) &theSocketAddress->addr, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_SocketNoDelay::recvVector() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with one other TCP_SocketNoDelay\n"; 
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
	nread = read(sockfd,gMsg,nleft);
	nleft -= nread;
	gMsg +=  nread;
    }
    return 0;
}


// void Send(Vector &):
// 	Method to send a Vector to an address given by other_Addr.

int 
TCP_SocketNoDelay::sendVector(const Vector &theVector, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_SocketNoDelay can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_SocketNoDelay::sendObj() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with a TCP_SocketNoDelay";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		
	if (memcmp((char *) &other_Addr, (char *) &theSocketAddress->addr, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_SocketNoDelay::recvVector() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with one other TCP_SocketNoDelay\n"; 
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
	nwrite = write(sockfd,gMsg,nleft);
	nleft -= nwrite;
	gMsg +=  nwrite;
    }
    return 0;
}




int 
TCP_SocketNoDelay::recvID(ID &theID, ChannelAddress *theAddress)
		    
{	
    // first check address is the only address a TCP_SocketNoDelay can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_SocketNoDelay::sendObj() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with a TCP_SocketNoDelay";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		
	if (memcmp((char *) &other_Addr, (char *) &theSocketAddress->addr, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_SocketNoDelay::recvID() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with one other TCP_SocketNoDelay\n"; 
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
	nread = read(sockfd,gMsg,nleft);
	nleft -= nread;
	gMsg +=  nread;
    }
    return 0;
}


// void Send(ID &):
// 	Method to send a ID to an address given by other_Addr.

int 
TCP_SocketNoDelay::sendID(const ID &theID, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_SocketNoDelay can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) 
	    theSocketAddress = (SocketAddress *)theAddress;
	else {
	    opserr << "TCP_SocketNoDelay::sendObj() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with a TCP_SocketNoDelay";
	    opserr << " address given is not of type SocketAddress\n"; 
	    return -1;	    
	}		
	if (memcmp((char *) &other_Addr, (char *) &theSocketAddress->addr, 
	     theSocketAddress->addrLength) != 0) {

	    opserr << "TCP_SocketNoDelay::recvID() - a TCP_SocketNoDelay ";
	    opserr << "can only communicate with one other TCP_SocketNoDelay\n"; 
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
	nwrite = write(sockfd,gMsg,nleft);
	nleft -= nwrite;
	gMsg +=  nwrite;
    }
    return 0;
}






unsigned int 
TCP_SocketNoDelay::getPortNumber(void) const
{
    return myPort;
}


char *
TCP_SocketNoDelay::addToProgram(void)
{
    char *tcp = " 3 ";

    char  me[20];
    char  my_InetAddr[MAX_INET_ADDR];
    char  myPortNum[8];
 
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
	break;
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

