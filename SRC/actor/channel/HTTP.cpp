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
// $Date: 2007-07-16 22:52:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/HTTP.cpp,v $
                                                                        
// Written: fmk 11/06
// Revised:
//
// Purpose: This file contains all the includes needed by files for http & https connections

#include <Socket.h>

#ifdef _HTTPS
#include <openssl/crypto.h>
#include <openssl/x509.h>
#include <openssl/pem.h>
#include <openssl/ssl.h>
#include <openssl/err.h>
#endif


socket_type
establishHTTPConnection(const char* URL, unsigned int port) {

  union {
    struct sockaddr    addr;
    struct sockaddr_in addr_in;
  } my_Addr;

  union {
    struct sockaddr    addr;
    struct sockaddr_in addr_in;
  } other_Addr;

  socket_type sockfd;
  socklen_type addrLength;
  struct hostent *hostEntry;
  struct in_addr ip;
  unsigned int myPort;

  /* check inputs */
  if (URL == 0) {
    return -1;
  }

  /*   
   *   connect to remote socket
   */
  
  // set up remote address
  bzero((char *) &other_Addr.addr_in, sizeof(other_Addr.addr_in));
  other_Addr.addr_in.sin_family      = AF_INET;
  other_Addr.addr_in.sin_port        = htons(port);

  hostEntry = gethostbyname(URL);
  bcopy(hostEntry->h_addr, &(ip.s_addr), hostEntry->h_length);
  
#ifdef _WIN32
  memcpy(&(other_Addr.addr_in.sin_addr.S_un.S_addr), hostEntry->h_addr, hostEntry->h_length);
#else
  memcpy(&(other_Addr.addr_in.sin_addr.s_addr), hostEntry->h_addr, hostEntry->h_length);
#endif
  
  /* set up my_Addr.addr_in  */
  bzero((char *) &my_Addr.addr_in, sizeof(my_Addr.addr_in));    
  my_Addr.addr_in.sin_family = AF_INET;
  my_Addr.addr_in.sin_port = htons(0);
  
#ifdef _WIN32
  my_Addr.addr_in.sin_addr.S_un.S_addr = htonl(INADDR_ANY);
#else
  my_Addr.addr_in.sin_addr.s_addr = htonl(INADDR_ANY);
#endif
  
  addrLength = sizeof(my_Addr.addr_in);
  
  /* open a socket */
  if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
    fprintf(stderr,"establishHTTPConnection - could not open socket\n");
    return -2;
  }
  
  /*  bind local address to it */
  if (bind(sockfd, (struct sockaddr *) &my_Addr.addr_in,sizeof(my_Addr.addr_in)) < 0) {
    fprintf(stderr,"establishHTTPConnection - could not bind local address\n");
    return -3;
  }
  myPort = ntohs(my_Addr.addr_in.sin_port);    
  
  
  /* now try to connect to socket with remote address. */
  if (connect(sockfd, (struct sockaddr *) &other_Addr.addr_in, 
	      sizeof(other_Addr.addr_in))< 0) {
    
    fprintf(stderr,"establishHTTPConnection - could not connect\n");
    return -4;
  }

  /* get my_address info */
  getsockname(sockfd, &my_Addr.addr, &addrLength);

  return sockfd;
}


#ifdef _WIN32
int __cdecl
#else
int
#endif
httpGet(char const *URL, char const *page, unsigned int port, char **dataPtr) {

  int i, j, nleft, nwrite, sizeData, ok;
  char *gMsg, *data, *nextData;
  char outBuf[4096], inBuf[4096];
  socket_type sockfd;

  // in case we fail, set return pointer to 0
  *dataPtr = 0;

  // invoke startup sockets
  startup_sockets();
  
  // open a socket
  sockfd = establishHTTPConnection(URL, port);
  if (sockfd < 0) {
    fprintf(stderr, "postData: failed to establis connection\n");
    return -1;
  }

  // add the header information to outBuf
  sprintf(outBuf, "GET %s HTTP/1.1\nHost:%s\n",page,URL);
  strcat(outBuf, "Connection:close\n\n");
  nleft = strlen(outBuf);

  //send the data
  // if o.k. get a ponter to the data in the message and 
  // place the incoming data there
  nwrite = 0;    
  gMsg = outBuf;
  
  while (nleft > 0) {
    nwrite = send(sockfd, gMsg, nleft, 0);
    nleft -= nwrite;
    gMsg +=  nwrite;
  }

  ok = 1;
  nleft = 4095;

  sizeData = 0;
  nextData = 0;
  data = 0;

  while (ok > 0) {

    gMsg = inBuf;
    ok = recv(sockfd, gMsg, nleft, 0);

    if (ok > 0) {
      nextData = data;
      data = (char *)malloc((sizeData+ok+1)*sizeof(char));
      if (data != 0) {
	if (nextData != 0) {
	  for (i=0; i<sizeData; i++)
	    data[i]=nextData[i];
	  free(nextData);
	}
	for (i=0, j=sizeData; i<ok; i++, j++)
	  data[j]=inBuf[i];
	sizeData += ok;
	strcpy(&data[sizeData],"");
      }
    }
  }


  // now we need to strip off the response header 
  gMsg = data;
  nextData = strstr(data,"Content-Type");
  if (nextData != NULL) {
    nextData = strchr(nextData,'\n');
    nextData += 3;

    nwrite = sizeData+1-(nextData-data);

    data = (char *)malloc((sizeData+1)*sizeof(char));
    for (i=0; i<nwrite; i++)
      data[i]=nextData[i];
    //    strcpy(&data[nwrite],""); /we already placed a end-of-string marker there above
  }

  *dataPtr = data;

  cleanup_sockets();
  
  return 0;
}


#ifdef _HTTPS

#ifdef _WIN32
int __cdecl
#else
int
#endif
httpsGet(char const *URL, char const *page, unsigned int port, char **dataPtr) {

  int i, j, nleft, nwrite, sizeData, ok;
  char *gMsg, *data, *nextData;
  char outBuf[4096], inBuf[4096];
  socket_type sockfd;

  SSL *ssl;
  SSL_CTX *ctx;
  SSL_METHOD *client_method;
  X509 *server_cert;

  int sd, err;

  /* ******************************************************
   * init SSL library from http: the definitive guide
   * **************************************************** */
  SSLeay_add_ssl_algorithms();
  client_method = SSLv2_client_method();
  SSL_load_error_strings();
  ctx = SSL_CTX_new(client_method);
  /* ********************* from http:the definitive guide */  

  // in case we fail, set return pointer to 0
  *dataPtr = 0;

  // invoke startup sockets
  startup_sockets();
  
  // open a socket
  sockfd = establishHTTPConnection(URL, port);
  if (sockfd < 0) {
    fprintf(stderr, "postData: failed to establis connection\n");
    return -1;
  }

  fprintf(stderr,"CONNECTION ESTABLISHED\n");


  /* ******************************************************
   * init SSL handshake from http: the definitive guide
   * **************************************************** */
  ssl = SSL_new(ctx);
  SSL_set_fd(ssl, sockfd);
  err = SSL_connect(ssl);

  fprintf(stderr,"endpoint create & handshake completed\n");
  fprintf(stderr,"cipher: %s\n\n",SSL_get_cipher(ssl));
  /* ********************* from http:the definitive guide */

  sprintf(outBuf, "GET %s HTTP/1.1\nHost:%s\n",page,URL);
  strcat(outBuf, "Connection:close\n\n");
  nleft = strlen(outBuf);

  fprintf(stderr,"%s",outBuf);

  err = SSL_write(ssl, outBuf, nleft);

  /*
  //send the data
  // if o.k. get a ponter to the data in the message and 
  // place the incoming data there
  nwrite = 0;    
  gMsg = outBuf;
  
  while (nleft > 0) {
    nwrite = send(sockfd, gMsg, nleft, 0);
    nleft -= nwrite;
    gMsg +=  nwrite;
  }
  */


  ok = 1;
  nleft = 4095;

  sizeData = 0;
  nextData = 0;
  data = 0;

  while (ok > 0) {

    gMsg = inBuf;
    ok = SSL_read(ssl, gMsg, nleft);

    if (ok > 0) {
      nextData = data;
      data = (char *)malloc((sizeData+ok+1)*sizeof(char));
      if (data != 0) {
	if (nextData != 0) {
	  for (i=0; i<sizeData; i++)
	    data[i]=nextData[i];
	  free(nextData);
	}
	for (i=0, j=sizeData; i<ok; i++, j++)
	  data[j]=inBuf[i];
	sizeData += ok;
	strcpy(&data[sizeData],"");
      }
    }
  }



  // now we need to strip off the response header 
  gMsg = data;
  nextData = strstr(data,"Content-Type");
  if (nextData != NULL) {
    nextData = strchr(nextData,'\n');
    nextData += 3;

    nwrite = sizeData+1-(nextData-data);

    data = (char *)malloc((sizeData+1)*sizeof(char));
    for (i=0; i<nwrite; i++)
      data[i]=nextData[i];
    //    strcpy(&data[nwrite],""); /we already placed a end-of-string marker there above
  }

  *dataPtr = data;


  SSL_shutdown(ssl);
#ifdef _WIN32
  closesocket(sockfd);
#else
  close(sockfd);
#endif
  SSL_free(ssl);
  SSL_CTX_free(ctx);

  cleanup_sockets();
  
  return 0;
}


#ifdef _WIN32
int __cdecl
#else
int
#endif
httpsPOST(const char *URL, 
	 const char *page, 
	 const char *cookie, 
	 const char *contentType,
	 const char *dataToPost,
	 unsigned int port, 
	 bool returnHeader, 
	 char **resPtr) {

  int i, j, nleft, nwrite, sizeData, ok;
  char *gMsg, *data, *oldData, *nextData;
  char outBuf[4096], inBuf[4096];
  socket_type sockfd;

  SSL *ssl;
  SSL_CTX *ctx;
  SSL_METHOD *client_method;
  X509 *server_cert;

  int sd, err;

  // in case we fail, set return pointer to 0
  *resPtr = 0;

  /* 
   * init SSL library 
   * code taken from O'Reilly book: 'http: the definitive guide'
   */
  SSLeay_add_ssl_algorithms();
  client_method = SSLv2_client_method();
  SSL_load_error_strings();
  ctx = SSL_CTX_new(client_method);
  /* end of code taken from http: the definitive guide */  

  // invoke startup sockets
  startup_sockets();
  
  // open a socket
  sockfd = establishHTTPConnection(URL, port);
  if (sockfd < 0) {
    fprintf(stderr, "postData: failed to establis connection\n");
    return -1;
  }

  /* 
   * SSL handshake
   * code taken from O'Reilly book: 'http: the definitive guide'
   */
  ssl = SSL_new(ctx);
  SSL_set_fd(ssl, sockfd);
  err = SSL_connect(ssl);
  /* end of code taken from http: the definitive guide */  

  sprintf(outBuf, "POST %s HTTP/1.1\nHost: %s\n", page, URL);  
  if (cookie != 0)
    strcat(outBuf, cookie);
  if (contentType == 0)
    strcat(outBuf, "Content-Type:text/plain\n");
  else {
    sprintf(inBuf, "Content-Type: %s\n", contentType);
    strcat(outBuf, inBuf);
  }

  //  strcat(outBuf, "Connection:close\n\n");
  sprintf(inBuf, "Content-Length: %d\n\n", strlen(dataToPost));
  strcat(outBuf, inBuf);
  strcat(outBuf, dataToPost);

  nleft = strlen(outBuf);
  err = SSL_write(ssl, outBuf, nleft);

  /*
  //send the data
  // if o.k. get a ponter to the data in the message and 
  // place the incoming data there
  nwrite = 0;    
  gMsg = outBuf;
  
  while (nleft > 0) {
    nwrite = send(sockfd, gMsg, nleft, 0);
    nleft -= nwrite;
    gMsg +=  nwrite;
  }
  */

  ok = 1;
  nleft = 4095;

  sizeData = 0;
  oldData = 0;
  data = 0;

  while (ok > 0) {

    gMsg = inBuf;
    ok = SSL_read(ssl, gMsg, nleft);

    if (ok > 0) {
      oldData = data;
      data = (char *)malloc((sizeData+ok+1)*sizeof(char));
      if (data != 0) {
	if (oldData != 0) {
	  for (i=0; i<sizeData; i++)
	    data[i]=oldData[i];
	  free(oldData);
	}
	for (i=0, j=sizeData; i<ok; i++, j++)
	  data[j]=inBuf[i];
	sizeData += ok;
	strcpy(&data[sizeData],"");
      }
    }
  }

  // now we need to strip off the response header 
  if (returnHeader == false) {
    oldData = data;
    gMsg = data;
    nextData = strstr(data,"Content-Type");
    if (nextData != NULL) {
      nextData = strchr(nextData,'\n');
      nextData += 3;
      
      nwrite = sizeData+1-(nextData-data);
      
      data = (char *)malloc((sizeData+1)*sizeof(char));
      for (i=0; i<nwrite; i++)
	data[i]=nextData[i];

      free(oldData);
      //    strcpy(&data[nwrite],""); /we already placed a end-of-string marker there above
    }
    
  }

  *resPtr = data;

  /*
   * shut-down ssl, close socket & free related memory
   */

  SSL_shutdown(ssl);
#ifdef _WIN32
  closesocket(sockfd);
#else
  close(sockfd);
#endif
  SSL_free(ssl);
  SSL_CTX_free(ctx);
  cleanup_sockets();
  
  return 0;
}


#ifdef _WIN32
int __cdecl
#else
int
#endif
neesPOST(const char *user,
	 const char *pass,
	 const char *page, 
	 const char *postData, 
	 const char *contentType) {

  char *res = 0;
  char *URL="central.nees.org";
  char *loginPage = "/login.php";
  char *loginType="application/x-www-form-urlencoded";
  char *cookie = 0;
  char loginData[512];

  sprintf(loginData,"redirect=https%s3A%s2F%s2Fcentral.nees.org%s2Findex.php%s3F&user=%s&pass=%s","%","%","%","%","%",user,pass);  

  /*
   * login & get GAsession cookie
   */

  httpsPOST(URL, loginPage, cookie, loginType, loginData,  443, true, &res);

  // make sure GAsession in cookie
  char *GAsession = strstr(res,"GAsession");  
  if (GAsession == 0) {
    fprintf(stderr, "neesPOST: could not login, check username, password or central.nees.org availability\n");
    return -1;
  }

  /*
   * set the coookie data for next call
   */

  cookie = (char *)malloc(8*sizeof(char));
  int sizeCookie = 8;
  strcpy(cookie, "Cookie:");

  char *startNextCookie = strstr(res,"Set-Cookie");
  while (startNextCookie != 0) {

    startNextCookie += 11; // "Set-Cookie:");
    char *endNextCookie = strstr(startNextCookie, ";");
    int sizeNextCookie = endNextCookie-startNextCookie+1;

    char *nextCookie = (char *)malloc((sizeNextCookie+sizeCookie+1)*sizeof(char));
    
    strncpy(nextCookie, cookie, sizeCookie-1);
    strncpy(&nextCookie[sizeCookie-1], startNextCookie, sizeNextCookie);

    free(cookie);
    cookie = nextCookie;
    sizeCookie += sizeNextCookie;

    strcpy(&cookie[sizeCookie],"");

    startNextCookie++;
    startNextCookie = strstr(startNextCookie,"Set-Cookie");
  }
  strcpy(&cookie[sizeCookie-2],"");

  free(res);
  
  httpsPOST(URL, page, cookie, contentType, postData,  443, true, &res);
  
  char *success = strstr(res,"HTTP/1.1 201 Success");

  if (success != 0) {
    free(res);
    return 0;
  } else {
    fprintf(stderr, "neesPOST: could not post\n NEES ERROR: %s\n", res);
    free(res);
    return -1;
  }
}

#endif
