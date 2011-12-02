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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/mySocket.c,v $
                                                                        
                                                                        

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

int INET_getsockname(int socket, struct sockaddr_in *address, int *address_len )
{
    return getsockname(socket, address, address_len);
}

int INET_sendto (
          int socket,
          char *message_addr,
          int length,
          int flags,
          struct sockaddr_in *dest_addr,
          int dest_len )
{
    return sendto(socket, message_addr, length, flags, dest_addr, dest_len);
}

int INET_recvfrom(
          int socket,
          char *buffer,
          int length,
          int flags,
          struct sockaddr_in *address,
          int *address_len) 
{ 
    return recvfrom(socket,buffer,length,flags,address, address_len);
}

int INET_bind (
          int socket,
          struct sockaddr_in *address,
          int address_len )
{
    return bind(socket, address, address_len);
}


