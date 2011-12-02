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
// $Date: 2006-12-15 00:12:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/Socket.cpp,v $
                                                                        
// Written: fmk 11/06
// Revised:
//
// Purpose: This file contains functions needed in socket calls

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef _WIN32
extern "C" {
#include <winsock2.h>
};
#else
  #include <sys/socket.h>
  #include <sys/types.h>
  #include <netinet/in.h>
  #include <netinet/tcp.h>
  #include <arpa/inet.h>
  #include <netdb.h>
  #include <unistd.h>
  #include <strings.h>
#endif


#define MAX_UDP_DATAGRAM 9126
#define MAX_INET_ADDR 28

#ifdef _WIN32
  typedef SOCKET socket_type;
  typedef int socklen_type;
  static int numSockets = 0;
#else
  typedef int socket_type;
  typedef socklen_t socklen_type;
#endif


int
startup_sockets(void)
{
  #ifdef _WIN32
  
  WSADATA wsaData;
  if (numSockets == 0) {
    numSockets++;
    return WSAStartup(0x0002, &wsaData) == 0;
  } else
    return 1;
  #else
    return 1;
  #endif
} 


void cleanup_sockets(void)
{
  #ifdef _WIN32
  numSockets--;
  if (numSockets == 0) {
    WSACleanup();
  }
  #endif
}


