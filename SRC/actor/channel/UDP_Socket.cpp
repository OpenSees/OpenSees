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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-10-15 00:31:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/UDP_Socket.cpp,v $
                                                                        
                                                                        
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the implementation of the methods needed
// to define the UDP_Socket class interface, given below:
//

#include "UDP_Socket.h"
#include <strings.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Message.h>
#include <SocketAddress.h>
#include <MovableObject.h>


static int GetHostAddr(char *host, char *IntAddr);
static void inttoa(unsigned int no, char *string, int *cnt);

// UDP_Socket(): 
//	constructor to open a socket with my inet_addr and with a port number 
//	assigned by the OS from the available free port numbers.

UDP_Socket::UDP_Socket() 
  :sockfd(0), connectType(0), shadow_inetAddr(0), shadow_port(0) 
{
    // set up my_Addr 
    my_Addr.addr_in.sin_family = AF_INET;
    my_Addr.addr_in.sin_addr.s_addr = htonl(INADDR_ANY); 
    my_Addr.addr_in.sin_port = htons(0);
    addrLength = sizeof(my_Addr.addr_in);
    
    // open a socket
    if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
	opserr << "UDP_Socket::UDP_Socket - could not open socket\n";
    }
    
    // bind local address to it
    if (bind(sockfd, (struct sockaddr *) &my_Addr.addr_in,sizeof(my_Addr.addr_in)) < 0) {
	opserr << "UDP_Socket::UDP_Socket - could not bind local address\n";
    }    
    
    getsockname(sockfd, &my_Addr.addr, &addrLength);
    myPort = ntohs(my_Addr.addr_in.sin_port);        
}    




// UDP_Socket(unsigned int port): 
//	constructor to open a socket with my inet_addr and with a port number port.

UDP_Socket::UDP_Socket(unsigned int port) 
  :sockfd(0), connectType(0), shadow_inetAddr(0), shadow_port(0)
{
    // set up my_Addr.addr_in 
    my_Addr.addr_in.sin_family = AF_INET;
    my_Addr.addr_in.sin_addr.s_addr = htonl(INADDR_ANY);
    my_Addr.addr_in.sin_port = htons(port);
    addrLength = sizeof(my_Addr.addr_in);
    
    // open a socket
    if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
	opserr << "UDP_Socket::UDP_Socket - could not open socket\n";
    }
    
    // bind local address to it
    if (bind(sockfd, (struct sockaddr *) &my_Addr.addr_in,sizeof(my_Addr.addr_in)) < 0) {
	opserr << "UDP_Socket::UDP_Socket - could not bind local address\n";
    }    
    getsockname(sockfd, &my_Addr.addr, &addrLength);    
    myPort = ntohs(my_Addr.addr_in.sin_port);        
}    




// UDP_Socket(unsigned int other_Port, char *other_InetAddr): 
// 	constructor to open a socket with my inet_addr and with a port number 
//	given by the OS. It then sends a messege comprising of a single character
// 	to an address given by other_InetAddr and other_Port. This is to allow
// 	a shadow object to find address of the actor it initiatites.

UDP_Socket::UDP_Socket(unsigned int other_Port, char *other_InetAddr) 
:sockfd(0), connectType(1), shadow_inetAddr(other_InetAddr), shadow_port(other_Port)
{
    // set up my_Addr.addr_in 
    my_Addr.addr_in.sin_family = AF_INET;
    my_Addr.addr_in.sin_addr.s_addr = htonl(INADDR_ANY); 
    my_Addr.addr_in.sin_port = htons(0);
    addrLength = sizeof(my_Addr.addr_in);
    
    // open a socket
    if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
	opserr << "UDP_Socket::UDP_Socket - could not open socket\n";
    }
    
    // bind local address to it
    if (bind(sockfd, (struct sockaddr *) &my_Addr.addr_in,sizeof(my_Addr.addr_in)) < 0) {
	opserr << "UDP_Socket::UDP_Socket - could not bind local address\n";
    }    
    getsockname(sockfd, &my_Addr.addr, &addrLength);    
    myPort = ntohs(my_Addr.addr_in.sin_port);        
}    

    


// ~UDP_Socket():
//	destructor

UDP_Socket::~UDP_Socket()
{
    close(sockfd);
}


int 
UDP_Socket::setUpConnection(void)
{
  if (connectType == 1) {
    if (shadow_inetAddr == 0) {
	opserr << "FATAL::UDP_Socket::setUpActor() -";
	opserr << "incorrect choice of UDP_Socket constructor for actor.\n";
	exit(-1);
    }

    // send a message to address.
    SocketAddress other(shadow_inetAddr, shadow_port);
    char data = 'a';
    sendto(sockfd, &data, 1, 0, &other.address.addr, other.addrLength);    
    bcopy((char *) &other.address.addr, (char *) &last_Addr, addrLength);
    
    // receive a message from other
    recvfrom(sockfd, &data, 1, 0, &last_Addr.addr, &addrLength);    

  } else {

    // wait for remote process to send message;
    char data = 'b';
    recvfrom(sockfd, &data, 1, 0, &last_Addr.addr, &addrLength);    
    
    // then send a message back
    sendto(sockfd, &data, 1, 0, &last_Addr.addr, addrLength);        
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
	
	bcopy((char *) &theSocketAddress->address.addr, (char *) &last_Addr, 
	      theSocketAddress->addrLength);
	addrLength = theSocketAddress->addrLength;	

	return  0;	    
    }
    else {
	opserr << "TCP_Socket::setNextAddress() - a UDP_Socket ";
	opserr << "can only communicate with a UDP_Socket";
	opserr << " address given is not of type SocketAddress\n"; 
	return -1;	    
    }		    	
}



int 
UDP_Socket::sendObj(int commitTag,
		    MovableObject &theObject, 
		    ChannelAddress *theAddress) 
{
    // set the address of last_Addr
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
	if (theAddress->getType() == SOCKET_TYPE) {
	    theSocketAddress = (SocketAddress *)theAddress;
	
	    bcopy((char *) &theSocketAddress->address.addr, (char *) &last_Addr.addr, 
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

    // now get the object to send itself
    return theObject.sendSelf(commitTag, *this);
}


int 
UDP_Socket::recvObj(int commitTag,
		    MovableObject &theObject, 
		    FEM_ObjectBroker &theBroker, 
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

	    if (bcmp((char *) &theSocketAddress->address.addr, (char *) &last_Addr.addr, 
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
UDP_Socket::recvMsg(int dbTag, int commitTag, Message &msg, ChannelAddress *theAddress )
{	

    // if o.k. get a ponter to the data in the message and 
    // place the incoming data there
    int size;
    char *gMsg;
    gMsg = msg.data;
    size = msg.length;

    while (size > 0) {
	if (size <= MAX_UDP_DATAGRAM) {
	    recvfrom(sockfd, gMsg, size, 0, &last_Addr.addr, &addrLength);
	    size = 0;
	}
	else {
	  recvfrom(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &last_Addr.addr, &addrLength);
		   
	    
	  gMsg += MAX_UDP_DATAGRAM;
	  size -= MAX_UDP_DATAGRAM;
	}
    }

    // check the address that message came from was correct
    if (theAddress != 0) {
	SocketAddress *theSocketAddress = 0;

	if (theAddress->getType() == SOCKET_TYPE) {
	    theSocketAddress = (SocketAddress *)theAddress;	
	    
	    if (bcmp((char *) &theSocketAddress->address.addr, (char *) &last_Addr.addr, 
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
UDP_Socket::sendMsg(int dbTag, int commitTag, const Message &msg, ChannelAddress *theAddress )
{	
    // set up the address of the Socket to which data will be sent
    if (theAddress != 0) {
	SocketAddress *theSocketAddress = 0;

	if (theAddress->getType() == SOCKET_TYPE) {
	    theSocketAddress = (SocketAddress *)theAddress;	
	
	    bcopy((char *) &theSocketAddress->address.addr, (char *) &last_Addr.addr, 
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
	    sendto(sockfd, gMsg, size, 0, &last_Addr.addr, addrLength);
	    size = 0;
	}
	else {	
	  sendto(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &last_Addr.addr, addrLength);
	  gMsg += MAX_UDP_DATAGRAM;
	  size -= MAX_UDP_DATAGRAM;
	}
    }
    return 0;
}






int
UDP_Socket::recvMatrix(int dbTag, int commitTag, Matrix &theMatrix, ChannelAddress *theAddress )
{	

    // if o.k. get a ponter to the data in the message and 
    // place the incoming data there
    int size;
    double *data = theMatrix.data;
    char *gMsg = (char *)data;;
    size =  theMatrix.dataSize * sizeof(double);

    while (size > 0) {
	if (size <= MAX_UDP_DATAGRAM) {
	    recvfrom(sockfd, gMsg, size, 0, &last_Addr.addr, &addrLength);
	    size = 0;
	}
	else {
	  recvfrom(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &last_Addr.addr, &addrLength);
	  gMsg += MAX_UDP_DATAGRAM;
	  size -= MAX_UDP_DATAGRAM;
	}
    }

    // check the address that message came from was correct
    if (theAddress != 0) {
	SocketAddress *theSocketAddress = 0;

	if (theAddress->getType() == SOCKET_TYPE) {
	    theSocketAddress = (SocketAddress *)theAddress;	
	    
	    if (bcmp((char *) &theSocketAddress->address.addr, (char *) &last_Addr.addr, 
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
UDP_Socket::sendMatrix(int dbTag, int commitTag, const Matrix &theMatrix, ChannelAddress *theAddress )
{	
    // set up the address of the Socket to which data will be sent
    if (theAddress != 0) {
	SocketAddress *theSocketAddress = 0;

	if (theAddress->getType() == SOCKET_TYPE) {
	    theSocketAddress = (SocketAddress *)theAddress;	
	
	    bcopy((char *) &theSocketAddress->address.addr, (char *) &last_Addr.addr, 
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
    size =  theMatrix.dataSize * sizeof(double);

    while (size > 0) {
	if (size <= MAX_UDP_DATAGRAM) {
	    sendto(sockfd, gMsg, size, 0, &last_Addr.addr, addrLength);
	    size = 0;
	}
	else {
	    sendto(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &last_Addr.addr, addrLength);
			
	    gMsg += MAX_UDP_DATAGRAM;
	    size -= MAX_UDP_DATAGRAM;
	}
    }
    return 0;
}






int
UDP_Socket::recvVector(int dbTag, int commitTag, Vector &theVector, ChannelAddress *theAddress )
{	

    // if o.k. get a ponter to the data in the message and 
    // place the incoming data there
    int size;
    double *data = theVector.theData;
    char *gMsg = (char *)data;;
    size =  theVector.sz * sizeof(double);

    while (size > 0) {
	if (size <= MAX_UDP_DATAGRAM) {
	    recvfrom(sockfd, gMsg, size, 0, &last_Addr.addr, &addrLength);
	    size = 0;
	}
	else {
	    recvfrom(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &last_Addr.addr, &addrLength);
	    
	    gMsg += MAX_UDP_DATAGRAM;
	    size -= MAX_UDP_DATAGRAM;
	}
    }

    // check the address that message came from was correct
    if (theAddress != 0) {
	SocketAddress *theSocketAddress = 0;

	if (theAddress->getType() == SOCKET_TYPE) {
	    theSocketAddress = (SocketAddress *)theAddress;	
	    
	    if (bcmp((char *) &theSocketAddress->address.addr, (char *) &last_Addr.addr, 
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
UDP_Socket::sendVector(int dbTag, int commitTag, const Vector &theVector, ChannelAddress *theAddress )
{	
    // set up the address of the Socket to which data will be sent
    if (theAddress != 0) {
	SocketAddress *theSocketAddress = 0;

	if (theAddress->getType() == SOCKET_TYPE) {
	    theSocketAddress = (SocketAddress *)theAddress;	
	
	    bcopy((char *) &theSocketAddress->address.addr, (char *) &last_Addr.addr, 
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
    size =  theVector.sz * sizeof(double);

    while (size > 0) {
	if (size <= MAX_UDP_DATAGRAM) {
	    sendto(sockfd, gMsg, size, 0, &last_Addr.addr, addrLength);
	    size = 0;
	}
	else {
	    sendto(sockfd, gMsg, MAX_UDP_DATAGRAM, 0, &last_Addr.addr, addrLength);
	    gMsg += MAX_UDP_DATAGRAM;
	    size -= MAX_UDP_DATAGRAM;
	}
    }
    return 0;
}






int
UDP_Socket::recvID(int dbTag, int commitTag, ID &theID, ChannelAddress *theAddress )
{	

    // if o.k. get a ponter to the data in the message and 
    // place the incoming data there
    int size;
    int *data = theID.data;
    char *gMsg = (char *)data;;
    size =  theID.sz * sizeof(int);

    while (size > 0) {
	if (size <= MAX_UDP_DATAGRAM) {
	    recvfrom(sockfd,gMsg,size,0,&last_Addr.addr,&addrLength);
	    size = 0;
	}
	else {
	    recvfrom(sockfd,gMsg,MAX_UDP_DATAGRAM,0, &last_Addr.addr,&addrLength);
			  
	    
	    gMsg += MAX_UDP_DATAGRAM;
	    size -= MAX_UDP_DATAGRAM;
	}
    }

    // check the address that message came from was correct
    if (theAddress != 0) {
	SocketAddress *theSocketAddress = 0;

	if (theAddress->getType() == SOCKET_TYPE) {
	    theSocketAddress = (SocketAddress *)theAddress;	
	    
	    if (bcmp((char *) &theSocketAddress->address.addr, (char *) &last_Addr.addr, 
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
UDP_Socket::sendID(int dbTag, int commitTag, const ID &theID, ChannelAddress *theAddress )
{	
    // set up the address of the Socket to which data will be sent
    if (theAddress != 0) {
	SocketAddress *theSocketAddress = 0;

	if (theAddress->getType() == SOCKET_TYPE) {
	    theSocketAddress = (SocketAddress *)theAddress;	
	
	    bcopy((char *) &theSocketAddress->address.addr, (char *) &last_Addr.addr, 
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
    size =  theID.sz * sizeof(int);

    while (size > 0) {
	if (size <= MAX_UDP_DATAGRAM) {
	    sendto(sockfd,gMsg,size,0,&last_Addr.addr,addrLength);
	    size = 0;
	}
	else {
	    sendto(sockfd,gMsg,MAX_UDP_DATAGRAM,0,
			&last_Addr.addr, addrLength);
	    gMsg += MAX_UDP_DATAGRAM;
	    size -= MAX_UDP_DATAGRAM;
	}
    }
    return 0;
}


unsigned int 
UDP_Socket::getPortNumber(void) const
{
    return myPort;
}



char *
UDP_Socket::addToProgram(void)
{
    char *tcp = " 2 ";
    
    char  me[20];
    char  my_InetAddr[MAX_INET_ADDR];
    char  myPortNum[8];
    unsigned int thePort = this->getPortNumber();

    int start = 0;
    inttoa(thePort,myPortNum,&start);
    gethostname(me,MAX_INET_ADDR);
    GetHostAddr(me,my_InetAddr);

    strcpy(add,tcp);
    strcat(add," ");          
    strcat(add,my_InetAddr);
    strcat(add," ");
    strcat(add,myPortNum);
    strcat(add," ");    
    
    return add;
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

