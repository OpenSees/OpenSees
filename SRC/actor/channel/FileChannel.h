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
**   Alisa Neeman (aneeman@cse.ucsc.edu)                              **
**                                                                    **
** ****************************************************************** */
                                                                        

// Purpose: This file contains the class definition for FileChannel.
// FileChannel is a sub-class of channel. The idea is that you have the
// identical domain running on a visualization node. You would like to
// update the visualization with the next time step. The entire domain is 
// saved to file using Domain.sendSelf() and then read from file using
// Domain.recvSelf()

#ifndef FileChannel_h
#define FileChannel_h

#include <bool.h>
#include <Channel.h>
#include <stdio.h>
#include <MovableObject.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <nDarray.h>
#include <Domain.h>
#include <FEM_ObjectBroker.h>

class nDarray; // not sure why the compiler needs this..

enum FileOpenType{ FILE_READ, FILE_WRITE };

class FileChannel : public Channel
{
  public:
    FileChannel(int fileOpenType);        
    FileChannel(const char * filename, int fileOpenType );    
    ~FileChannel();

    char *addToProgram(void);
    
    virtual int setUpConnection(void);

    int setNextAddress(const ChannelAddress &otherChannelAddress);
    virtual ChannelAddress *getLastSendersAddress(void){ return 0;};

    int sendObj(int commitTag,
		MovableObject &theObject, 
		ChannelAddress *theAddress =0);
    int recvObj(int commitTag,
		MovableObject &theObject, 
		FEM_ObjectBroker &theBroker,
		ChannelAddress *theAddress =0);
	
	int sendDomain( int commitTag, Domain * theDomain );
	int recvDomain( int commitTag, Domain * theDomain );
			
    int sendMsg(int dbTag, int commitTag, 
		const Message &, 
		ChannelAddress *theAddress =0);    
    int recvMsg(int dbTag, int commitTag, 
		Message &, 
		ChannelAddress *theAddress =0);        
    int recvMsgUnknownSize(int dbTag, int commitTag,
			   Message & theMessage,
			   ChannelAddress *theAddress =0);


    int sendMatrix(int dbTag, int commitTag, 
		   const Matrix &theMatrix, 
		   ChannelAddress *theAddress =0);
    int recvMatrix(int dbTag, int commitTag, 
		   Matrix &theMatrix, 
		   ChannelAddress *theAddress =0);
    
    int sendVector(int dbTag, int commitTag, 
		   const Vector &theVector, ChannelAddress *theAddress =0);
    int recvVector(int dbTag, int commitTag, 
		   Vector &theVector, 
		   ChannelAddress *theAddress =0);
    
    int sendID(int dbTag, int commitTag, 
	       const ID &theID, 
	       ChannelAddress *theAddress =0);
    int recvID(int dbTag, int commitTag, 
	       ID &theID, 
	       ChannelAddress *theAddress =0);  
	
	int sendnDarray(int dbTag, int commitTag,  const nDarray & theNDarray,ChannelAddress *theAddress);

	int recvnDarray(int dbTag, int commitTag, nDarray &theNDarray, ChannelAddress *theAddress);
	
	int nextFile();  
	void closeFile();
	int setCommitStep(int commitTag);   

	
  protected:
    
    
  private:
    FILE * theFile;
    char fileNameBase[200];
	char fileName[210];
	int commitStep;
	int fileType; // read or write
};


#endif 






 
