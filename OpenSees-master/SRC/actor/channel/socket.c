#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


#ifdef _WIN32
#include <winsock2.h>
#define bcopy(s,d,n)	memmove((d),(s),(n))
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
  #define bzero(s,n) memset((s),0,(n))
  #define memcmp(s1,s2,n) memcmp((s1),(s2),(n))
  static int numSockets = 0;
#else
  typedef int socket_type;
  typedef socklen_t socklen_type;
#endif

int startup_sockets(void)
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


int
httpGet(const char *URL, const char *page, unsigned int port, char **dataPtr) {


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


  sprintf(outBuf, "GET %s HTTP/1.1\nHost:%s\n",page,URL);
  strcat(outBuf, "Connection:close\n\n");
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


