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
                                                                        
// $Revision: 1.15 $
// $Date: 2007-04-05 01:15:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/FileDatastore.h,v $
                                                                        
                                                                        
#ifndef FileDatastore_h
#define FileDatastore_h

// Written: fmk 
// Created: 10/98
//
// Description: This file contains the class definition for FileDatastore.
// FileDatastore is a concrete subclas of FE_Datastore. A FileDatastore 
// object is used in the program to store/restore the geometry and state 
// information in a domain at a particular instance in the analysis. The
// information is stored in text files.
//
// What: "@(#) FileDatastore.h, revA"

#include <FE_Datastore.h>

#include <fstream>
#include <map>
using std::fstream;
using std::map;


#define STREAM_POSITION_TYPE int

class FEM_ObjectBroker;

typedef struct fileDatastoreOutputFile {
  fstream *theFile;
  STREAM_POSITION_TYPE fileEnd;
  int      maxDbTag;
} FileDatastoreOutputFile;

typedef map<int, FileDatastoreOutputFile *>      MAP_FILES;
typedef MAP_FILES::value_type                    MAP_FILES_TYPE;
typedef MAP_FILES::iterator                      MAP_FILES_ITERATOR;

typedef struct intData{
  int *dbTag;
  int *data;     
} IntData;


typedef struct doubleData{
  int *dbTag;
  double *data;     
} DoubleData;


class FileDatastore: public FE_Datastore
{
  public:
    FileDatastore(const char *dataBase,
		  Domain &theDomain, 
		  FEM_ObjectBroker &theBroker);    
    
    ~FileDatastore();

    // methods for sending and receiving the data
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

    int createTable(const char *tableName, int numColumns, char *columns[]);
    int insertData(const char *tableName, char *columns[], 
		   int commitTag, const Vector &data);
    int getData(const char *tableName, char *columns[], int commitTag, Vector &data);

    // the commitState method
    int commitState(int commitTag);        
    
  protected:

  private:
    // Private methods
    int resizeInt(int newSize);
    int resizeDouble(int newSize);
    void resetFilePointers(void);
    int openFile(char *fileName, FileDatastoreOutputFile *, int dataSize);

    // private attributes
    char *dataBase;
    MAP_FILES theIDFiles;
    MAP_FILES theVectFiles;
    MAP_FILES theMatFiles;
    MAP_FILES_ITERATOR theIDFilesIter;
    MAP_FILES_ITERATOR theVectFilesIter;
    MAP_FILES_ITERATOR theMatFilesIter;

    int lastDomainChangeStamp;
    int currentCommitTag;
    char *data;
    int sizeData;
    
    IntData    theIntData;
    DoubleData theDoubleData;

    int currentMaxInt;
    int currentMaxDouble;
};




#endif

