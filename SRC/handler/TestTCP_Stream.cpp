
// Written: fmk

#include <StandardStream.h>
#include <TCP_Socket.h>
#include <Vector.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

// main routine
int main(int argc, char **argv)
{
  Vector data(1);
  int dataSize = 0;

  if (argc != 2) {
    opserr << "ERROR - Usage a.out inetPort\n";
    exit(0);
  }

  int inetPort = atoi(argv[1]);
  TCP_Socket *theSocket = new TCP_Socket(inetPort);
  if (theSocket == 0) {
    opserr << "ERROR - Failed to Create a TCP_Socket\n";
    exit(0);
  }
  theSocket->setUpConnection();
  
  theSocket->recvVector(0,0,data);

  int sizeReceived = data(0);
  while (sizeReceived != -1) {
    if (sizeReceived != dataSize) {
      data.resize(sizeReceived+1);
      dataSize = sizeReceived;
    }

    if (theSocket->recvVector(0,0, data) == 0) {
      sizeReceived = data(0);
      if (sizeReceived == dataSize) {
	for (int i=1; i<=dataSize; i++)
	  opserr << data(i) << " " ;
	opserr << endln;
      }
    }
  }
    
  delete theSocket;
  exit(0);
}
  
