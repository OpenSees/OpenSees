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

// $Revision: 317 $
// $Date: 2011-06-03 19:18:55 -0700 (Fri, 03 Jun 2011) $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/TCP_SocketSSL.h,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 06/07
// Revision: A
//
// Purpose: This file contains the class definition for TCP_SocketSSL.
// TCP_SocketSSL is a sub-class of channel. It is implemented with Berkeley
// stream sockets using the TCP protocol and OpenSSL. Message delivery is 
// guaranteed. Communication is full-duplex between a pair of connected sockets.

#ifndef TCP_SocketSSL_h
#define TCP_SocketSSL_h

#include <bool.h>
#include <Socket.h>
#include <Channel.h>

#include <openssl/ssl.h>

#define RSA_SERVER_KEY "server.key"
#define RSA_SERVER_CERT "server.crt"
#define RSA_SERVER_CA_CERT "server_ca.crt"
#define RSA_CLIENT_KEY "client.key"
#define RSA_CLIENT_CERT "client.crt"
#define RSA_CLIENT_CA_CERT "client_ca.crt"

#define VERIFY_CLIENT true

class TCP_SocketSSL : public Channel
{
public:
    TCP_SocketSSL();        
    TCP_SocketSSL(unsigned int port, bool checkEndianness = false,
        int noDelay = 0);    
    TCP_SocketSSL(unsigned int other_Port, const char *other_InetAddr,
        bool checkEndianness = false, int noDelay = 0);
    ~TCP_SocketSSL();

    char *addToProgram();

    virtual int setUpConnection();

    int setNextAddress(const ChannelAddress &otherChannelAddress);
    virtual ChannelAddress *getLastSendersAddress(){ return 0;};

    int sendObj(int commitTag,
        MovableObject &theObject, 
        ChannelAddress *theAddress =0);
    int recvObj(int commitTag,
        MovableObject &theObject, 
        FEM_ObjectBroker &theBroker,
        ChannelAddress *theAddress =0);

    int sendMsg(int dbTag, int commitTag, 
        const Message &, 
        ChannelAddress *theAddress =0);    
    int recvMsg(int dbTag, int commitTag, 
        Message &, 
        ChannelAddress *theAddress =0);        
    int recvMsgUnknownSize(int dbTag, int commitTag, 
        Message &, 
        ChannelAddress *theAddress =0);        

    int sendMatrix(int dbTag, int commitTag, 
        const Matrix &theMatrix, 
        ChannelAddress *theAddress =0);
    int recvMatrix(int dbTag, int commitTag, 
        Matrix &theMatrix, 
        ChannelAddress *theAddress =0);

    int sendVector(int dbTag, int commitTag, 
        const Vector &theVector,
        ChannelAddress *theAddress =0);
    int recvVector(int dbTag, int commitTag, 
        Vector &theVector, 
        ChannelAddress *theAddress =0);

    int sendID(int dbTag, int commitTag, 
        const ID &theID, 
        ChannelAddress *theAddress =0);
    int recvID(int dbTag, int commitTag, 
        ID &theID, 
        ChannelAddress *theAddress =0);    

protected:
    unsigned int getPortNumber() const;
    unsigned int getBytesAvailable();

private:
    SSL *ssl;
    SSL_CTX *ctx;
    const SSL_METHOD *meth;
    X509 *cert;

    socket_type sockfd;

    union {
        struct sockaddr    addr;
        struct sockaddr_in addr_in;
    } my_Addr;
    union {
        struct sockaddr    addr;
        struct sockaddr_in addr_in;
    } other_Addr;

    socklen_type addrLength;

    unsigned int myPort;
    int connectType;
    bool checkEndianness;
    bool endiannessProblem;
    int noDelay;
};

#endif 
