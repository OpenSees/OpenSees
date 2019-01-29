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
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/TCP_SocketSSL.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 06/07
// Revision: A
//
// Purpose: This file contains the implementation of the methods needed
// to define the TCP_SocketSSL class interface.

#include "TCP_SocketSSL.h"
#include <string.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Message.h>
#include <ChannelAddress.h>
#include <MovableObject.h>
#include <SocketAddress.h>

#include <openssl/crypto.h>
#include <openssl/x509.h>
#include <openssl/pem.h>
#include <openssl/err.h>

static int GetHostAddr(char *host, char *IntAddr);
static void inttoa(unsigned int no, char *string, int *cnt);

#ifndef _WIN32
static void byte_swap(void *array, long long nArray, int size);
#endif


// TCP_SocketSSL(unsigned int other_Port, char *other_InetAddr): 
// 	constructor to open a socket with my inet_addr and with a port number 
//	given by the OS.
TCP_SocketSSL::TCP_SocketSSL()
    : myPort(0), connectType(0),
    checkEndianness(false), endiannessProblem(false), noDelay(0)
{
    // initialize SSL library and context object
    SSL_library_init();
    SSL_load_error_strings();
    meth = TLS_server_method();
    ctx = SSL_CTX_new(meth);

    // server verification
    // load server certificate into context object
    if (SSL_CTX_use_certificate_file(ctx, RSA_SERVER_CERT, SSL_FILETYPE_PEM) <= 0) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not load server certificate\n";
        exit(-1);
    }
    // load private key corresponding to server certificate
    if (SSL_CTX_use_PrivateKey_file(ctx, RSA_SERVER_KEY, SSL_FILETYPE_PEM) <= 0) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not load server key\n";
        exit(-1);
    }
    // check if server certificate and private key match
    if (!SSL_CTX_check_private_key(ctx)) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - server private key does not match the certificate public key\n";
        exit(-1);
    }

    // client verification
    if (VERIFY_CLIENT) {
        // load server CA certificate into context object
        if (!SSL_CTX_load_verify_locations(ctx, RSA_SERVER_CA_CERT, NULL)) {
            opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not load server CA certificate\n";
            exit(-1);
        }
        SSL_CTX_set_verify(ctx, SSL_VERIFY_PEER, NULL);
        SSL_CTX_set_verify_depth(ctx, 1);
    }

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
    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not open socket\n";
    }

    // bind local address to it
    if (bind(sockfd, &my_Addr.addr, sizeof(my_Addr.addr)) < 0) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not bind local address\n";
    }

    // get my_address info
    addrLength = sizeof(my_Addr.addr);
    myPort = ntohs(my_Addr.addr_in.sin_port);
}    


// TCP_SocketSSL(unsigned int port): 
//	constructor to open a socket with my inet_addr and with a port number port.
TCP_SocketSSL::TCP_SocketSSL(unsigned int port, bool checkendianness, int nodelay) 
    : myPort(0), connectType(0),
    checkEndianness(checkendianness), endiannessProblem(false),
    noDelay(nodelay)
{
    // initialize SSL library and context object
    SSL_library_init();
    SSL_load_error_strings();
    meth = TLS_server_method();
    ctx = SSL_CTX_new(meth);

    // server verification
    // load server certificate into context object
    if (SSL_CTX_use_certificate_file(ctx, RSA_SERVER_CERT, SSL_FILETYPE_PEM) <= 0) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not load server certificate\n";
        exit(-1);
    }
    // load private key corresponding to server certificate
    if (SSL_CTX_use_PrivateKey_file(ctx, RSA_SERVER_KEY, SSL_FILETYPE_PEM) <= 0) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not load server key\n";
        exit(-1);
    }
    // check if server certificate and private key match
    if (!SSL_CTX_check_private_key(ctx)) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - server private key does not match the certificate public key\n";
        exit(-1);
    }

    // client verification
    if (VERIFY_CLIENT) {
        // load server CA certificate into context object
        if (!SSL_CTX_load_verify_locations(ctx, RSA_SERVER_CA_CERT, NULL)) {
            opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not load server CA certificate\n";
            exit(-1);
        }
        SSL_CTX_set_verify(ctx, SSL_VERIFY_PEER, NULL);
        SSL_CTX_set_verify_depth(ctx, 1);
    }

    // initialize sockets
    startup_sockets();

    // set up my_Addr.addr_in with address given by port and internet address of
    // machine on which the process that uses this routine is running.

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
    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not open socket\n";
    }

    // bind local address to it
    if (bind(sockfd, &my_Addr.addr, sizeof(my_Addr.addr)) < 0) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not bind local address\n";
    }    

    // get my_address info
    addrLength = sizeof(my_Addr.addr);
    myPort = ntohs(my_Addr.addr_in.sin_port);
}


// TCP_SocketSSL(unsigned int other_Port, char *other_InetAddr): 
// 	constructor to open a socket with my inet_addr and with a port number 
//	given by the OS. Then to connect with a TCP_SocketSSL whose address is
//	given by other_Port and other_InetAddr.
TCP_SocketSSL::TCP_SocketSSL(unsigned int other_Port,
    const char *other_InetAddr, bool checkendianness, int nodelay)
    : myPort(0), connectType(1),
    checkEndianness(checkendianness), endiannessProblem(false),
    noDelay(nodelay)
{
    // initialize SSL library and context object
    SSL_library_init();
    SSL_load_error_strings();
    meth = TLS_client_method();
    ctx = SSL_CTX_new(meth);

    // client verification
    if (VERIFY_CLIENT) {
        // load client certificate into context object
        if (SSL_CTX_use_certificate_file(ctx, RSA_CLIENT_CERT, SSL_FILETYPE_PEM) <= 0) {
            opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not load client certificate\n";
            exit(-1);
        }
        // load private key corresponding to client certificate
        if (SSL_CTX_use_PrivateKey_file(ctx, RSA_CLIENT_KEY, SSL_FILETYPE_PEM) <= 0) {
            opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not load client key\n";
            exit(-1);
        }
        // check if client certificate and private key match
        if (!SSL_CTX_check_private_key(ctx)) {
            opserr << "TCP_SocketSSL::TCP_SocketSSL() - client private key does not match the certificate public key\n";
            exit(-1);
        }
    }

    // server verification
    // load client CA certificate into context object
    if (!SSL_CTX_load_verify_locations(ctx, RSA_CLIENT_CA_CERT, NULL)) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not load client CA certificate\n";
        exit(-1);
    }
    SSL_CTX_set_verify(ctx, SSL_VERIFY_PEER, NULL);
    SSL_CTX_set_verify_depth(ctx, 1);

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
    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not open socket\n";
    }

    // bind local address to it
    if (bind(sockfd, &my_Addr.addr, sizeof(my_Addr.addr)) < 0) {
        opserr << "TCP_SocketSSL::TCP_SocketSSL() - could not bind local address\n";
    }

    addrLength = sizeof(my_Addr.addr);
    getsockname(sockfd, &my_Addr.addr, &addrLength);
    myPort = ntohs(my_Addr.addr_in.sin_port);
}


// ~TCP_SocketSSL():
//	destructor
TCP_SocketSSL::~TCP_SocketSSL()
{
    // shutdown SSL connection
    SSL_shutdown(ssl);

#ifdef _WIN32
    closesocket(sockfd);
#else
    close(sockfd);
#endif

    SSL_free(ssl);
    SSL_CTX_free(ctx);
    cleanup_sockets();
}


int 
TCP_SocketSSL::setUpConnection()
{
    if (connectType == 1) {
        
        // now try to connect to socket with remote address.
        if (connect(sockfd, &other_Addr.addr, sizeof(other_Addr.addr))< 0) {
                opserr << "TCP_SocketSSL::setUpConnection() - could not connect\n";
                return -1;
        }
        
        // get my_address info
        getsockname(sockfd, &my_Addr.addr, &addrLength);
        
        // set TCP_NODELAY option
        if ((setsockopt(sockfd, IPPROTO_TCP, TCP_NODELAY, 
            (char *) &noDelay, sizeof(int))) < 0) { 
                opserr << "TCP_SocketSSL::setUpConnection() - "
                    << "could not set TCP_NODELAY option\n";
        }
        
        // initiate SSL handshake with server
        ssl = SSL_new(ctx);
        SSL_set_fd(ssl, (int)sockfd);
        if (SSL_connect(ssl) <= 0) {
            opserr << "TCP_SocketSSL::setUpConnection() - could not initiate SSL handshake with server\n";
            SSL_free(ssl);
            SSL_CTX_free(ctx);
            cleanup_sockets();
            return -2;
        }

        // informational output
        opserr << "\nSSL cipher: " << SSL_get_cipher(ssl) << endln;

        // get the server's certificate and display
        cert = SSL_get_peer_certificate(ssl);
        if (cert != NULL) {
            opserr << "\nSSL Server certificate:\n"
                << "Subject name: "
                << X509_NAME_oneline(X509_get_subject_name(cert),0,0)
                << "\nIssuer name: "
                << X509_NAME_oneline(X509_get_issuer_name(cert),0,0)
                << endln;
            X509_free(cert);
        } else {
            opserr << "\nSSL Server has no certificate\n";
        }

        // check for endianness problem if requested
        if (checkEndianness) {
            int i = 1;
            int j;

            int *data = &i;
            char *gMsg = (char *)data;
            SSL_write(ssl, gMsg, sizeof(int));

            data = &j;
            gMsg = (char *)data;
            SSL_read(ssl, gMsg, sizeof(int));

            if (i != j) {
                int k = 0x41424344;
                char *c = (char *)&k;
                if (*c == 0x41)
                    endiannessProblem = true;
            }
        }

    } else {

        // wait for other process to contact me & set up connection
        socket_type newsockfd;
        listen(sockfd, 1);    
        newsockfd = accept(sockfd, &other_Addr.addr, &addrLength);
        if (newsockfd < 0) {
            opserr << "TCP_SocketSSL::setUpConnection() - could not accept connection\n";
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
        
        // set TCP_NODELAY option
        if ((setsockopt(sockfd, IPPROTO_TCP, TCP_NODELAY, 
            (char *) &noDelay, sizeof(int))) < 0) { 
                opserr << "TCP_SocketSSL::setUpConnection() - "
                    << "could not set TCP_NODELAY option\n";
        }
        
        // initiate SSL handshake with client
        ssl = SSL_new(ctx);
        SSL_set_fd(ssl, (int)sockfd);
        if (SSL_accept(ssl) <= 0) {
            opserr << "TCP_SocketSSL::setUpConnection() - could not initiate SSL handshake with client\n";
            SSL_free(ssl);
            SSL_CTX_free(ctx);
            cleanup_sockets();
            return -2;
        }

        // informational output
        opserr << "\nSSL cipher: " << SSL_get_cipher(ssl) << endln;

        if (VERIFY_CLIENT) {
            // get the client's certificate and display
            cert = SSL_get_peer_certificate(ssl);
            if (cert != NULL) {
                opserr << "\nSSL Client certificate:\n"
                    << "Subject name: "
                    << X509_NAME_oneline(X509_get_subject_name(cert),0,0)
                    << "\nIssuer name: "
                    << X509_NAME_oneline(X509_get_issuer_name(cert),0,0)
                    << endln;
                X509_free(cert);
            } else {
                opserr << "\nSSL Client has no certificate\n";
            }
        }

        // check for endianness problem if requested
        if (checkEndianness) {
            int i;
            int j = 1;

            int *data = &i;
            char *gMsg = (char *)data;
            SSL_read(ssl, gMsg, sizeof(int));

            data = &j;
            gMsg = (char *)data;
            SSL_write(ssl, gMsg, sizeof(int));

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
TCP_SocketSSL::setNextAddress(const ChannelAddress &theAddress)
{	
    SocketAddress *theSocketAddress = 0;
    if (theAddress.getType() == SOCKET_TYPE) {
        theSocketAddress = (SocketAddress *)(&theAddress);    
        // check address is the only address a TCP_socket can send to
        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::setNextAddress() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL\n"; 
                return -1;
        }
    }
    else {
        opserr << "TCP_SocketSSL::setNextAddress() - a TCP_SocketSSL ";
        opserr << "can only communicate with a TCP_SocketSSL";
        opserr << " address given is not of type SocketAddress\n"; 
        return -1;	    
    }		    	

    return 0;
}


int 
TCP_SocketSSL::sendObj(int commitTag,
    MovableObject &theObject, ChannelAddress *theAddress) 
{
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
        if (theAddress->getType() == SOCKET_TYPE) 
            theSocketAddress = (SocketAddress *)theAddress;
        else {
            opserr << "TCP_SocketSSL::sendObj() - a TCP_SocketSSL ";
            opserr << "can only communicate with a TCP_SocketSSL";
            opserr << " address given is not of type SocketAddress\n"; 
            return -1;	    
        }		    

        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::sendObj() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL";
                opserr << " address given is not that address\n"; 
                return -1;	    
        }	
    }

    return theObject.sendSelf(commitTag, *this);
}


int 
TCP_SocketSSL::recvObj(int commitTag,
    MovableObject &theObject, FEM_ObjectBroker &theBroker, 
    ChannelAddress *theAddress)
{
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
        if (theAddress->getType() == SOCKET_TYPE) 
            theSocketAddress = (SocketAddress *)theAddress;
        else {
            opserr << "TCP_SocketSSL::recvObj() - a TCP_SocketSSL ";
            opserr << "can only communicate with a TCP_SocketSSL";
            opserr << " address given is not of type SocketAddress\n"; 
            return -1;	    
        }		    
        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::recvObj() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL\n"; 
                return -1;
        }
    }

    return theObject.recvSelf(commitTag, *this, theBroker);
}


// void Recv(Message &):
// 	Method to receive a message, also sets other_Addr.addr_in to that of sender
int 
TCP_SocketSSL::recvMsg(int dbTag, int commitTag,
    Message &msg, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
        if (theAddress->getType() == SOCKET_TYPE) 
            theSocketAddress = (SocketAddress *)theAddress;
        else {
            opserr << "TCP_SocketSSL::recvMsg() - a TCP_SocketSSL ";
            opserr << "can only communicate with a TCP_SocketSSL";
            opserr << " address given is not of type SocketAddress\n"; 
            return -1;	    
        }		    
        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::recvMsg() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL\n"; 
                return -1;
        }
    }

    // if o.k. get a pointer to the data in the message and 
    // place the incoming data there
    int nleft, nread;
    char *gMsg;
    gMsg = msg.data;
    nleft = msg.length;

    while (nleft > 0) {
        nread = SSL_read(ssl, gMsg, nleft);
        if (nread < 0) {
            opserr << "TCP_SocketSSL::recvMsg() - could not read data\n";
            SSL_shutdown(ssl);
            SSL_free(ssl);
            SSL_CTX_free(ctx);
            cleanup_sockets();
            return -2;
        }
        nleft -= nread;
        gMsg  += nread;
    }

    return 0;
}


// void Recv(Message &):
// 	Method to receive a message, also sets other_Addr.addr_in to that of sender
int 
TCP_SocketSSL::recvMsgUnknownSize(int dbTag, int commitTag,
    Message &msg, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
        if (theAddress->getType() == SOCKET_TYPE) 
            theSocketAddress = (SocketAddress *)theAddress;
        else {
            opserr << "TCP_SocketSSL::recvMsgUnknownSize() - a TCP_SocketSSL ";
            opserr << "can only communicate with a TCP_SocketSSL";
            opserr << " address given is not of type SocketAddress\n"; 
            return -1;	    
        }		    
        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::recvMsgUnknownSize() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL\n"; 
                return -1;
        }
    }

    // if o.k. get a pointer to the data in the message and 
    // place the incoming data there
    int nleft, nread;
    bool eol = false;
    char *gMsg;
    gMsg = msg.data;

    while (!eol) {
        nleft = SSL_pending(ssl);
        while (nleft > 0) {
            nread = SSL_read(ssl, gMsg, nleft);
            if (nread < 0) {
                opserr << "TCP_SocketSSL::recvMsgUnknownSize() - could not read data\n";
                SSL_shutdown(ssl);
                SSL_free(ssl);
                SSL_CTX_free(ctx);
                cleanup_sockets();
                return -2;
            }
            nleft -= nread;
            gMsg  += nread;
            if (*(gMsg-1) == '\0')
                eol = true;
            else if (*(gMsg-1) == '\n') {
                eol = true;
                *gMsg = '\0';
            }
        }
    }

    return 0;
}


// void Send(Message &):
// 	Method to send a message to an address given by other_Addr.addr_in.
int 
TCP_SocketSSL::sendMsg(int dbTag, int commitTag,
    const Message &msg, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
        if (theAddress->getType() == SOCKET_TYPE) 
            theSocketAddress = (SocketAddress *)theAddress;
        else {
            opserr << "TCP_SocketSSL::sendMsg() - a TCP_SocketSSL ";
            opserr << "can only communicate with a TCP_SocketSSL";
            opserr << " address given is not of type SocketAddress\n"; 
            return -1;	    
        }		    
        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::sendMsg() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL\n"; 
                return -1;
        }
    }

    // if o.k. get a pointer to the data in the message and 
    // place the incoming data there
    int nleft;    
    char *gMsg;
    gMsg = msg.data;
    nleft = msg.length;

    if (SSL_write(ssl, gMsg, nleft) != nleft) {
        opserr << "TCP_SocketSSL::sendMsg() - could not write data\n";
        SSL_shutdown(ssl);
        SSL_free(ssl);
        SSL_CTX_free(ctx);
        cleanup_sockets();
        return -2;
    }

    return 0;
}


int 
TCP_SocketSSL::recvMatrix(int dbTag, int commitTag,
    Matrix &theMatrix, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
        if (theAddress->getType() == SOCKET_TYPE) 
            theSocketAddress = (SocketAddress *)theAddress;
        else {
            opserr << "TCP_SocketSSL::recvMatrix() - a TCP_SocketSSL ";
            opserr << "can only communicate with a TCP_SocketSSL";
            opserr << " address given is not of type SocketAddress\n"; 
            return -1;	    
        }		    
        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::recvMatrix() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL\n"; 
                return -1;
        }
    }

    // if o.k. get a pointer to the data in the Matrix and 
    // place the incoming data there
    int nleft, nread;
    double *data = theMatrix.data;
    char *gMsg = (char *)data;;
    nleft = theMatrix.dataSize * sizeof(double);

    while (nleft > 0) {
        nread = SSL_read(ssl, gMsg, nleft);
        if (nread < 0) {
            opserr << "TCP_SocketSSL::recvMatrix() - could not read data\n";
            SSL_shutdown(ssl);
            SSL_free(ssl);
            SSL_CTX_free(ctx);
            cleanup_sockets();
            return -2;
        }
        nleft -= nread;
        gMsg  += nread;
    }

#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theMatrix.dataSize, sizeof(double));
    }
#endif

    return 0;
}


// void Send(Matrix &):
// 	Method to send a Matrix to an address given by other_Addr.addr_in.
int 
TCP_SocketSSL::sendMatrix(int dbTag, int commitTag,
    const Matrix &theMatrix, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
        if (theAddress->getType() == SOCKET_TYPE) 
            theSocketAddress = (SocketAddress *)theAddress;
        else {
            opserr << "TCP_SocketSSL::sendMatrix() - a TCP_SocketSSL ";
            opserr << "can only communicate with a TCP_SocketSSL";
            opserr << " address given is not of type SocketAddress\n"; 
            return -1;	    
        }		        SocketAddress *theSocketAddress = 0;

        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::sendMatrix() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL\n"; 
                return -1;
        }
    }

    // if o.k. get a pointer to the data in the Matrix and 
    // place the incoming data there
    int nleft;    
    double *data = theMatrix.data;
    char *gMsg = (char *)data;
    nleft = theMatrix.dataSize * sizeof(double);

#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theMatrix.dataSize, sizeof(double));
    }
#endif

    if (SSL_write(ssl, gMsg, nleft) != nleft) {
        opserr << "TCP_SocketSSL::sendMatrix() - could not write data\n";
        SSL_shutdown(ssl);
        SSL_free(ssl);
        SSL_CTX_free(ctx);
        cleanup_sockets();
        return -2;
    }

#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theMatrix.dataSize, sizeof(double));
    }
#endif

    return 0;
}


int 
TCP_SocketSSL::recvVector(int dbTag, int commitTag,
    Vector &theVector, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
        if (theAddress->getType() == SOCKET_TYPE) 
            theSocketAddress = (SocketAddress *)theAddress;
        else {
            opserr << "TCP_SocketSSL::recvVector() - a TCP_SocketSSL ";
            opserr << "can only communicate with a TCP_SocketSSL";
            opserr << " address given is not of type SocketAddress\n"; 
            return -1;	    
        }		
        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::recvVector() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL\n"; 
                return -1;
        }
    }

    // if o.k. get a pointer to the data in the Vector and 
    // place the incoming data there
    int nleft, nread;
    double *data = theVector.theData;
    char *gMsg = (char *)data;;
    nleft = theVector.sz * sizeof(double);

    while (nleft > 0) {
        nread = SSL_read(ssl, gMsg, nleft);
        if (nread < 0) {
            opserr << "TCP_SocketSSL::recvVector() - could not read data\n";
            SSL_shutdown(ssl);
            SSL_free(ssl);
            SSL_CTX_free(ctx);
            cleanup_sockets();
            return -2;
        }
        nleft -= nread;
        gMsg  += nread;
    }

#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theVector.sz, sizeof(double));
    }
#endif

    return 0;
}


// void Send(Vector &):
// 	Method to send a Vector to an address given by other_Addr.addr_in.
int 
TCP_SocketSSL::sendVector(int dbTag, int commitTag,
    const Vector &theVector, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
        if (theAddress->getType() == SOCKET_TYPE) 
            theSocketAddress = (SocketAddress *)theAddress;
        else {
            opserr << "TCP_SocketSSL::sendVector() - a TCP_SocketSSL ";
            opserr << "can only communicate with a TCP_SocketSSL";
            opserr << " address given is not of type SocketAddress\n"; 
            return -1;	    
        }		
        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::sendVector() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL\n"; 
                return -1;
        }
    }

    // if o.k. get a pointer to the data in the Vector and 
    // place the incoming data there
    int nleft;    
    double *data = theVector.theData;
    char *gMsg = (char *)data;
    nleft = theVector.sz * sizeof(double);

#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theVector.sz, sizeof(double));
    }
#endif

    if (SSL_write(ssl, gMsg, nleft) != nleft) {
        opserr << "TCP_SocketSSL::sendVector() - could not write data\n";
        SSL_shutdown(ssl);
        SSL_free(ssl);
        SSL_CTX_free(ctx);
        cleanup_sockets();
        return -2;
    }

#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theVector.sz, sizeof(double));
    }
#endif

    return 0;
}


int 
TCP_SocketSSL::recvID(int dbTag, int commitTag,
    ID &theID, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
        if (theAddress->getType() == SOCKET_TYPE) 
            theSocketAddress = (SocketAddress *)theAddress;
        else {
            opserr << "TCP_SocketSSL::recvID() - a TCP_SocketSSL ";
            opserr << "can only communicate with a TCP_SocketSSL";
            opserr << " address given is not of type SocketAddress\n"; 
            return -1;	    
        }		
        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::recvID() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL\n"; 
                return -1;
        }
    }

    // if o.k. get a pointer to the data in the ID and 
    // place the incoming data there
    int nleft, nread;
    int *data = theID.data;
    char *gMsg = (char *)data;;
    nleft = theID.sz * sizeof(int);

    while (nleft > 0) {
        nread = SSL_read(ssl, gMsg, nleft);
        if (nread < 0) {
            opserr << "TCP_SocketSSL::recvID() - could not read data\n";
            SSL_shutdown(ssl);
            SSL_free(ssl);
            SSL_CTX_free(ctx);
            cleanup_sockets();
            return -2;
        }
        nleft -= nread;
        gMsg  += nread;
    }

#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theID.sz, sizeof(int));
    }
#endif

    return 0;
}


// void Send(ID &):
// 	Method to send a ID to an address given by other_Addr.addr_in.
int 
TCP_SocketSSL::sendID(int dbTag, int commitTag,
    const ID &theID, ChannelAddress *theAddress)
{	
    // first check address is the only address a TCP_socket can send to
    SocketAddress *theSocketAddress = 0;
    if (theAddress != 0) {
        if (theAddress->getType() == SOCKET_TYPE) 
            theSocketAddress = (SocketAddress *)theAddress;
        else {
            opserr << "TCP_SocketSSL::sendID() - a TCP_SocketSSL ";
            opserr << "can only communicate with a TCP_SocketSSL";
            opserr << " address given is not of type SocketAddress\n"; 
            return -1;	    
        }		
        if (memcmp((char *) &other_Addr.addr_in, (char *) &theSocketAddress->address.addr_in, 
            theSocketAddress->addrLength) != 0) {

                opserr << "TCP_SocketSSL::sendID() - a TCP_SocketSSL ";
                opserr << "can only communicate with one other TCP_SocketSSL\n"; 
                return -1;
        }
    }

    // if o.k. get a pointer to the data in the ID and 
    // place the incoming data there
    int nleft;    
    int *data = theID.data;
    char *gMsg = (char *)data;
    nleft = theID.sz * sizeof(int);
 
#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theID.sz, sizeof(int));
    }
#endif

    if (SSL_write(ssl, gMsg, nleft) != nleft) {
        opserr << "TCP_SocketSSL::sendID() - could not write data\n";
        SSL_shutdown(ssl);
        SSL_free(ssl);
        SSL_CTX_free(ctx);
        cleanup_sockets();
        return -2;
    }

#ifndef _WIN32
    if (endiannessProblem) {
        void *array = (void *)data;
        byte_swap(array, theID.sz, sizeof(int));
    }
#endif

    return 0;
}


unsigned int 
TCP_SocketSSL::getPortNumber() const
{
    return myPort;
}


unsigned int
TCP_SocketSSL::getBytesAvailable()
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
TCP_SocketSSL::addToProgram()
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
inttoa(unsigned int no, char *string, int *cnt)
{
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
