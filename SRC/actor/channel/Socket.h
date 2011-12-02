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
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/Socket.h,v $
                                                                        
                                                                        
// File: ~/actor/Socket.h
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains all the includes needed by files
// using the Berkeley socket stuff.

#ifndef Socket_h
#define Socket_h

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <string.h>
#include <unistd.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>
}

#include <iOPS_Stream.h>
#include <bool.h>

#define MAX_UDP_DATAGRAM 9126
#define MAX_INET_ADDR 28

#ifdef _ALPHA
extern "C" int gethostname (char *name, int namelen);
#else
extern "C" int gethostname (char *name, unsigned int namelen);
#endif

extern "C" int INET_getsockname(int socket, 
				struct sockaddr_in *address, 
				int *address_len );

extern "C" int INET_sendto (int socket,
			    char *message_addr,
			    int length,
			    int flags,
			    struct sockaddr_in *dest_addr,
			    int dest_len );

extern "C" int INET_recvfrom(int socket,
			     char *buffer,
			     int length,
			     int flags,
			     struct sockaddr_in *address,
			     int *address_len) ;
extern "C" int INET_bind (int socket,
			  struct sockaddr_in *address,
			  int address_len );



#endif 


