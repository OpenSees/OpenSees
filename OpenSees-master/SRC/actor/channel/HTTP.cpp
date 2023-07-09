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
                                                                        
// $Revision: 1.11 $
// $Date: 2010-04-21 21:20:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/channel/HTTP.cpp,v $
                                                                        
// Written: fmk 11/06
// Revised:
//
// Purpose: This file contains all the includes needed by files for http & https connections

#include <Socket.h>

#include <stdio.h>
#include <sys/stat.h>

#ifdef _HTTPS
#include <openssl/crypto.h>
#include <openssl/x509.h>
#include <openssl/pem.h>
#include <openssl/ssl.h>
#include <openssl/err.h>
#endif

#define OUTBUF_SIZE 4096

static char outBuf[OUTBUF_SIZE+1];
static char inBuf[OUTBUF_SIZE+1];

static char *lastURL = 0;
//static socket_type lastSockfd;

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
  socket_type sockfd;

  // in case we fail, set return pointer to 0
  *dataPtr = 0;

  startup_sockets();

  sockfd = establishHTTPConnection(URL, port);
  if (sockfd < 0) {
    fprintf(stderr, "httpGet: failed to establis connection\n");
    return -1;
  }

  // add the header information to outBuf
  sprintf(outBuf, "GET %s HTTP/1.1\nHost:%s\n",page,URL);
  strcat(outBuf,"Accept:text/xml,text/html\n");
  strcat(outBuf,"Accept-Language:en-us,en\n");
  strcat(outBuf,"Accept-Charset:ISO-8859-1,utf-8\n");
  strcat(outBuf,"Keep-Alive:300\n");
  strcat(outBuf, "Connection:keep-alive\n\n");

  nleft = strlen(outBuf);

  //send the data
  // if o.k. get a pointer to the data in the message and 
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

    inBuf[ok+1]='\0';
         
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

    if (strstr(inBuf,"</html>") != NULL)
      ok = 0;

  }

  if (sizeData == 0) {
    if (lastURL != 0)
      free(lastURL);
    lastURL = 0;

#ifdef _WIN32
  closesocket(sockfd);
#else
  close(sockfd);
#endif

    return -1;
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
  }

  *dataPtr = data;
  free(gMsg);

#ifdef _WIN32
  closesocket(sockfd);
#else
  close(sockfd);
#endif
 
  
  
  sockfd = 0;

  cleanup_sockets();

  return 0;
}


#ifdef _WIN32
int __cdecl
#else
int
#endif
httpGET_File(char const *URL, char const *page, unsigned int port, const char *filename) {

  int nleft, nwrite, sizeData, ok;
  char *gMsg, *data, *nextData;
  char outBuf[OUTBUF_SIZE], inBuf[OUTBUF_SIZE];
  socket_type sockfd;

  FILE *fp = 0;


  fprintf(stderr, "httpGetFile URL: %s page %s\n", URL, page);  

  // invoke startup sockets
  startup_sockets();
  
  // open a socket
  sockfd = establishHTTPConnection(URL, port);
  if (sockfd < 0) {
    fprintf(stderr, "postData: failed to establis connection\n");
    return -1;
  }

  sockfd = establishHTTPConnection(URL, port);
  if (sockfd < 0) {
    fprintf(stderr, "httpGet: failed to establis connection\n");
    return -1;
  }

  // add the header information to outBuf
  sprintf(outBuf, "GET /%s HTTP/1.1\nHost:%s\n", page, URL);
  strcat(outBuf,"Keep-Alive:300\n");
  strcat(outBuf, "Connection:keep-alive\n\n");

  nleft = strlen(outBuf);

  //send the data
  // if o.k. get a pointer to the data in the message and 
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

  //
  // open file for writing
  //

  int fileOpened = 0;
  bool headerStripped = false;

  while (ok > 0) {

    gMsg = inBuf;
    ok = recv(sockfd, gMsg, nleft, 0);

    fprintf(stderr, "ok %d nleft %d\n", ok, nleft);

    if (ok > 0) {

      // now we need to strip off the response header 
      nextData = strstr(gMsg, "Bad");
      if (nextData != NULL) {
	fprintf(stderr, "Bad Request\n");	  
	return -1;
      }

      if (fileOpened == 0) {
	fp = fopen(filename,"wb");
	if (fp == 0) {
	  fprintf(stderr, "cannot open file %s for reading - is it still open for writing!\n", filename);
	  return -1;
	} else
	  fileOpened = 1;
      }

      if (headerStripped == false) {
	gMsg = inBuf;
	nextData = strstr(gMsg,"Content-Type");
	if (nextData != NULL) {
	  nextData = strchr(nextData,'\n');
	  nextData += 3;
	  
	  nwrite = sizeData+1-(nextData-data);
	  fwrite((void *)nextData, 1, nwrite, fp);
	  headerStripped = true;
	}
      } else {
	fwrite((void *)gMsg, 1, nleft, fp);
      }
    }
  }

  fprintf(stderr,"DONE\n");

  if (fileOpened == 1)
    fclose(fp);

  cleanup_sockets();
  
  return 0;
}



	/*
#ifdef _WIN32
int __cdecl
#else
int
#endif
httpGET_File(char const *URL, char const *page, unsigned int port, const char *filename) {

  int nleft, nwrite, sizeData, ok;
  char *gMsg, *data, *nextData;
  char outBuf[OUTBUF_SIZE], inBuf[OUTBUF_SIZE];
  socket_type sockfd;

  FILE *fp;

  fprintf(stderr, "httpGetFile URL: %s page %s\n", URL, page);  

  // invoke startup sockets
  startup_sockets();
  
  // open a socket
  sockfd = establishHTTPConnection(URL, port);
  if (sockfd < 0) {
    fprintf(stderr, "postData: failed to establis connection\n");
    return -1;
  }

  sockfd = establishHTTPConnection(URL, port);
  if (sockfd < 0) {
    fprintf(stderr, "httpGet: failed to establis connection\n");
    return -1;
  }

  // add the header information to outBuf
  sprintf(outBuf, "GET /%s HTTP/1.1\nHost:%s\n", page, URL);
  strcat(outBuf,"Keep-Alive:300\n");
  strcat(outBuf, "Connection:keep-alive\n\n");

  nleft = strlen(outBuf);

  //send the data
  // if o.k. get a pointer to the data in the message and 
  // place the incoming data there
  nwrite = 0;    
  gMsg = outBuf;

  while (nleft > 0) {
    nwrite = send(sockfd, gMsg, nleft, 0);
    nleft -= nwrite;
    gMsg +=  nwrite;
  }

  ok = 1024;
  nleft = 1024;

  sizeData = 0;
  nextData = 0;
  data = 0;

  //
  // open file for writing
  //

  int fileOpened = 0;
  bool headerStripped = false;

  while (ok == 1024) {

    gMsg = inBuf;
    ok = recv(sockfd, gMsg, nleft, 0);

    fprintf(stderr, "ok %d nleft %d\n", ok, nleft);

    if (ok > 0) {

      // now we need to strip off the response header 
      nextData = strstr(gMsg, "Bad");
      if (nextData != NULL) {
	fprintf(stderr, "Bad Request\n");	  
	return -1;
      }

      if (fileOpened == 0) {
	fp = fopen(filename,"wb");
	if (fp == 0) {
	  fprintf(stderr, "cannot open file %s for reading - is it still open for writing!\n", filename);
	  return -1;
	} else
	  fileOpened = 1;
      }

      //      headerStripped = true;

      if (headerStripped == false) {
	gMsg = inBuf;
	nextData = strstr(gMsg,"Content-Type");
	if (nextData != NULL) {
	  nextData = strchr(nextData,'\n');
	  nextData = strchr(nextData,'\n');
	  nextData += 3;
	  
	  nwrite = ok - (nextData-gMsg);
	  fwrite((void *)nextData, 1, nwrite, fp);
	  headerStripped = true;
	}
      } else {
	gMsg = inBuf;
	fwrite((void *)gMsg, 1, ok, fp);
      }
    }
  }

  fprintf(stderr,"DONE\n");

  if (fileOpened == 1)
    fclose(fp);

  cleanup_sockets();
  
  return 0;
}

	*/

#ifdef _HTTPS

#ifdef _WIN32
int __cdecl
#else
int
#endif
httpsGet(char const *URL, char const *page, char const *cookie, unsigned int port, char **dataPtr) {

  int i, j, nleft, nwrite, sizeData, ok;
  char *gMsg, *data, *nextData;
  char outBuf[OUTBUF_SIZE], inBuf[OUTBUF_SIZE];
  socket_type sockfd;

  SSL *ssl;
  SSL_CTX *ctx;
  
  int err;

  /* ******************************************************
   * init SSL library from http: the definitive guide
   * **************************************************** */
  SSLeay_add_ssl_algorithms();
  const SSL_METHOD *client_method = SSLv2_client_method();
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

  /* ******************************************************
   * init SSL handshake from http: the definitive guide
   * **************************************************** */
  ssl = SSL_new(ctx);
  SSL_set_fd(ssl, sockfd);
  err = SSL_connect(ssl);

  /* ********************* from http:the definitive guide */

  sprintf(outBuf, "GET %s HTTP/1.1\nHost:%s\n",page,URL);
  if (cookie != 0) 
    strcat(outBuf, cookie);
  strcat(outBuf, "Connection:close\n\n");
  nleft = strlen(outBuf);

  err = SSL_write(ssl, outBuf, nleft);

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
httpsSEND(const char *URL, 
	  const char *page, 
	  const char *cookie, 
	  const char *contentType,
	  const char *dataToPost,
	  unsigned int port, 
	  bool returnHeader, 
	  bool doPOST,
	  char **resPtr) {

  int i, j, nleft, nwrite, sizeData, ok;
  char *gMsg, *data, *oldData, *nextData;
  char outBuf[OUTBUF_SIZE], inBuf[OUTBUF_SIZE];
  socket_type sockfd;

  SSL *ssl;
  SSL_CTX *ctx;


  int err;

  // in case we fail, set return pointer to 0
  *resPtr = 0;

  /* 
   * init SSL library 
   * code taken from O'Reilly book: 'http: the definitive guide'
   */
  SSLeay_add_ssl_algorithms();
  const SSL_METHOD *  client_method = SSLv2_client_method();
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

  if (doPOST == true)
    sprintf(outBuf, "POST %s HTTP/1.1\nHost: %s\n", page, URL);  
  else
    sprintf(outBuf, "PUT %s HTTP/1.1\nHost: %s\n", page, URL);  

  if (cookie != 0) 
    strcat(outBuf, cookie);

  if (contentType == 0)
    strcat(outBuf, "Content-Type:text/plain\n");
  else {
    sprintf(inBuf, "Content-Type: %s\n", contentType);
    strcat(outBuf, inBuf);
  }

  //  strcat(outBuf, "Connection:close\n\n");
  int sizeDataPost = strlen(dataToPost);
  sprintf(inBuf, "Content-Length: %d\n\n", sizeDataPost);
  strcat(outBuf, inBuf);
  strcat(outBuf, dataToPost);

  nleft = strlen(outBuf);

  //send the data
  // if o.k. get a pointer to the data in the message and 
  // place the incoming data there
  nwrite = 0;    
  gMsg = outBuf;
  
  while (nleft > 0) {
    nwrite = SSL_write(ssl, gMsg, nleft);
    nleft -= nwrite;
    gMsg +=  nwrite;
  }

  // read the response

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
    if (ok < 4095)
      ok = 0;
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
httpsGET_File(char const *URL, char const *page, const char *cookie, unsigned int port, const char *filename) {

  int nleft, sizeData, ok;
  char *gMsg, *data, *nextData;
  char outBuf[OUTBUF_SIZE], inBuf[OUTBUF_SIZE];
  socket_type sockfd;

  FILE *fp;

  SSL *ssl;
  SSL_CTX *ctx;

  int err;

  /* ******************************************************
   * init SSL library from http: the definitive guide
   * **************************************************** */
  SSLeay_add_ssl_algorithms();
  const SSL_METHOD *  client_method = SSLv2_client_method();
  SSL_load_error_strings();
  ctx = SSL_CTX_new(client_method);
  /* ********************* from http:the definitive guide */  

  // invoke startup sockets
  startup_sockets();
  
  // open a socket
  sockfd = establishHTTPConnection(URL, port);
  if (sockfd < 0) {
    fprintf(stderr, "postData: failed to establis connection\n");
    return -1;
  }

  /* ******************************************************
   * init SSL handshake from http: the definitive guide
   * **************************************************** */
  ssl = SSL_new(ctx);
  SSL_set_fd(ssl, sockfd);
  err = SSL_connect(ssl);

  /* ********************* from http:the definitive guide */

  sprintf(outBuf, "GET %s HTTP/1.1\nHost:%s\n",page,URL);
  if (cookie != 0) 
    strcat(outBuf, cookie);
  strcat(outBuf, "Connection:close\n\n");
  nleft = strlen(outBuf);

  err = SSL_write(ssl, outBuf, nleft);

  ok = 1;
  nleft = 4095;

  sizeData = 0;
  nextData = 0;
  data = 0;

  //
  // open file for writing

  fp = fopen(filename,"wb");
  if (fp == 0) {
    fprintf(stderr, "cannot open file %s for reading - is it still open for writing!\n", filename);
    return -1;
  }

  while (ok > 0) {

    gMsg = inBuf;
    ok = SSL_read(ssl, gMsg, nleft);  
    if (ok < 0)
      fwrite((void *)gMsg, 1, nleft, fp);

    /*
    fprintf(stderr,"\n\nREAD %d\n", nleft);
    for (int i=0; i<nleft; i++)
      fprintf(stderr,"%c", inBuf[i]);
    */
    /*
    if (ok > 0) {
      // now we need to strip off the response header 
      if (headerStripped == false) {
	gMsg = data;
	nextData = strstr(data,"Content-Type");
	if (nextData != NULL) {
	  nextData = strchr(nextData,'\n');
	  nextData += 3;
	  
	  nwrite = sizeData+1-(nextData-data);
	  fwrite((void *)nextData, 1, nwrite, fp);
	  headerStripped = true;
	}
      } else {
	fwrite((void *)gMsg, 1, nleft, fp);
      }
    }
    */
  }
  
  fclose(fp);

  cleanup_sockets();
  
  return 0;
}


#ifdef _WIN32
int __cdecl
#else
int
#endif
httpsSEND_File(const char *URL, 
	       const char *page, 
	       const char *cookie, 
	       const char *contentType,
	       const char *filename,
	       unsigned int port, 
	       bool returnHeader, 
	       bool doPOST,
	       char **resPtr) 
{
  int i, j, nleft, nwrite, sizeData, ok;
  char *gMsg, *data, *oldData, *nextData;
  char outBuf[OUTBUF_SIZE], inBuf[OUTBUF_SIZE];
  socket_type sockfd;

  struct stat statResults;
  FILE *fp;
  
  SSL *ssl;
  SSL_CTX *ctx;


  int err;

  // in case we fail, set return pointer to 0
  *resPtr = 0;

  // get filesize
  int fileSize = 0;

  if (stat(filename, &statResults) == 0)
    fileSize =  statResults.st_size;

  if (fileSize == 0) {
    fprintf(stderr, "file %s has 0 fileSize; either no data or file does not exist!\n", filename);
    return -1;
  }

  fp = fopen(filename,"rb");
  if (fp == 0) {
    fprintf(stderr, "cannot open file %s for reading - is it still open for writing!\n", filename);
    return -1;
  }

     
  /* 
   * init SSL library 
   * code taken from O'Reilly book: 'http: the definitive guide'
   */
  SSLeay_add_ssl_algorithms();
  const SSL_METHOD *  client_method = SSLv2_client_method();
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

  if (doPOST == true)
    sprintf(outBuf, "POST %s HTTP/1.1\nHost: %s\n", page, URL);  
  else
    sprintf(outBuf, "PUT %s HTTP/1.1\nHost: %s\n", page, URL);  

  if (cookie != 0) 
    strcat(outBuf, cookie);

  if (contentType == 0)
    strcat(outBuf, "Content-Type:text/plain\n");
  else {
    sprintf(inBuf, "Content-Type: %s\n", contentType);
    strcat(outBuf, inBuf);
  }

  sprintf(inBuf, "Content-Length: %d\n\n", fileSize);
  strcat(outBuf, inBuf);

  //  strcat(outBuf, dataToPost);

  nleft = strlen(outBuf);

  //send the heading
  // if o.k. get a pointer to the data in the message and 
  // place the incoming data there
  nwrite = 0;    
  gMsg = outBuf;
  
  while (nleft > 0) {
    nwrite = SSL_write(ssl, gMsg, nleft);
    nleft -= nwrite;
    gMsg +=  nwrite;
  }


  int done = 0;
  while (done == 0) {
    nleft = fread((void *)outBuf, 1, OUTBUF_SIZE, fp);
    gMsg = outBuf;
    if (nleft < OUTBUF_SIZE)
      done = 1;
    while (nleft > 0) {
      nwrite = SSL_write(ssl, gMsg, nleft);
      nleft -= nwrite;
      gMsg +=  nwrite;
    }
    if (feof(fp) != 0)
      done = 1;
  }
  
  fclose(fp);

  //  err = SSL_write(ssl, outBuf, nleft);

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
    if (ok < 4095)
      ok = 0;
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

#endif // _HTTPS



