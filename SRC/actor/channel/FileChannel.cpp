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

#include <FileChannel.h>

FileChannel::FileChannel( int readWrite ) {
  strcpy(fileNameBase,"domain");
  commitStep = 0;
  fileType = readWrite;

  sprintf(fileName, "%s%d", fileNameBase,commitStep);
  if (fileType ==  FILE_READ )
	theFile = fopen(fileName, "rb");
  else 
	theFile = fopen(fileName, "wb");

  if (theFile == NULL) 
	opserr << "FileChannel::Unable to open " << fileName <<"\n";
}  
     
FileChannel::FileChannel( const char * filename, int readWrite){

  if (strlen(filename) > 199 ) {
	strncpy(fileNameBase, filename, 199 );
	opserr << " File name over 199 chars, using " << fileNameBase << "\n"; 
  }
  else strcpy(fileNameBase,filename); 
  commitStep = 0;

  sprintf(fileName, "%s%d", fileNameBase,commitStep);

  fileType = readWrite;

  if (fileType ==  FILE_READ )
	theFile = fopen(fileName, "rb");
  else 
	theFile = fopen(fileName, "wb");

  if (theFile == NULL) 
	opserr << "FileChannel::Unable to open " << fileName <<"\n";

}  
  
/**
 * Next commitStep increment, create a new file for writing
 * named with  appended commitStep number 
 */
int FileChannel::nextFile() {
  commitStep++;

  closeFile();

  sprintf(fileName, "%s%d", fileNameBase,commitStep);

  if (fileType ==  FILE_READ )
	theFile = fopen(fileName, "rb");
  else 
	theFile = fopen(fileName, "wb");

  if (theFile == NULL) {
	opserr << "FileChannel::Unable to open " << fileName <<"\n";

	return -1;
  }
  return 0;
}

int FileChannel::setCommitStep(int commitTag) {
  commitStep = commitTag;
  closeFile();

  sprintf(fileName, "%s%d", fileNameBase,commitStep);

  if (fileType ==  FILE_READ )
	theFile = fopen(fileName, "rb");
  else 
	theFile = fopen(fileName, "wb");

  if (theFile == NULL) {
	opserr << "FileChannel::Unable to open " << fileName <<"\n";

	return -1;
  }
  return 0;
}

/**
 * Close the file if hasn't already been closed and set the pointer
 * to NULL to be really sure
 */
void FileChannel::closeFile() {
  if (theFile != NULL) {
	fclose(theFile);
	theFile = NULL;
  }
}

FileChannel::~FileChannel(){
  closeFile();
}

char * 
FileChannel::addToProgram(void) 
{
  return 0;
}

int 
FileChannel::setUpConnection(void)
{
  return 0;
}

int FileChannel::setNextAddress(const ChannelAddress &otherChannelAddress)
{
  return 0;
}


int  FileChannel::sendDomain( int commitTag, Domain * theDomain ) 
{
  if ( !theFile )
	return -1;

  return theDomain->sendSelf(commitTag, *this);
}

int  FileChannel::recvDomain( int commitTag, Domain * theDomain )
{
  FEM_ObjectBroker theBroker;
  
  if ( !theFile ) {
    printf("FileChannel::recvDomain, File is not open\n");
    return -1;
  }
  
  setCommitStep( commitTag );
  return theDomain->recvSelf(commitTag, *this, theBroker);  
}


int FileChannel::sendObj(int commitTag,
			MovableObject &theObject, 
	   	  ChannelAddress *theAddress){
  if ( !theFile )
	return -1;

  return theObject.sendSelf(commitTag, *this); 
}


int FileChannel::recvObj(int commitTag,
			MovableObject &theObject, 
			FEM_ObjectBroker &theBroker,
						  ChannelAddress *theAddress) {
  if ( !theFile )
	return -1;


  setCommitStep( commitTag );
  return theObject.recvSelf(commitTag, *this, theBroker);
}


int FileChannel::sendMsg(int dbTag, int commitTag, const Message &, ChannelAddress *theAddress) 
{
  return 0;
}  
  
int FileChannel::recvMsg(int dbTag, int commitTag, Message &, ChannelAddress *theAddress) 
			 
{
  return 0;
}     


   
/**
* Implemented to match Channel superclass API, returns -1 for error.
*/
int FileChannel::recvMsgUnknownSize(int dbTag, int commitTag,

                Message & theMessage, ChannelAddress *theAddress)
{
opserr << "FileChannel::recvMsgUnknownSize, error, a FileChannel should not "
                                << "receive a message\n";
        return -1;
}
int FileChannel::sendMatrix(int dbTag, int commitTag, 
			    const Matrix &theMatrix, 
			    ChannelAddress *theAddress)
{
    // if o.k. get a ponter to the data in the Matrix and 
    // place the incoming data there
    int  nleft, rows, cols, i,j;
    rows = theMatrix.noRows();
    cols = theMatrix.noCols();

    double *data = new double[rows*cols];
    void *gMsg = (void *)data;
    nleft =  rows*cols* sizeof(double);

	for (i = 0; i < rows; i++)
	  for ( j = 0; j < cols; j++ )
	    data[j*rows + i] = theMatrix(i,j); //col*numRows + row = theMatrix(row,col)

	if ( fwrite( gMsg, nleft, 1, theFile) != 1 ) {
	  printf("Error in FileChannel::sendMatrix, did not write\n");
	  delete []data;
	  return -1;
	}
	delete [] data;
	return 0;

}
int FileChannel::recvMatrix(int dbTag, int commitTag, 
			   Matrix &theMatrix, 
			   ChannelAddress *theAddress)
{
   // if o.k. get a ponter to the data in the Matrix and 
    // place the incoming data there
    int  nleft, rows, cols;
    rows = theMatrix.noRows();
    cols = theMatrix.noCols();

    double *data = new double[rows*cols];
    void *gMsg = (void *)data;
    nleft =  rows*cols* sizeof(double);
    
	if ( fread( gMsg, nleft, 1, theFile) != 1 ) {
	  printf("Error in FileChannel::recvMatrix, did not read\n");
	  delete []data;
	  return -1;
	}
	
	else {
	  theMatrix.setData(data, rows, cols);
	  return 0;
	}
}

int FileChannel::sendVector(int dbTag, int commitTag, 
					 const Vector &theVector, ChannelAddress *theAddress ) 
{
    // if o.k. get a ponter to the data in the Vector and 
    // place the incoming data there
    int  nleft, size, i;
    size = theVector.Size();
    double *data = new  double[size];
    void *gMsg = (void *)data;
    nleft =  theVector.Size() * sizeof(double);

	for (i = 0; i < size; i++)
	  data[i] = theVector[i];
    
	if ( fwrite( gMsg, nleft, 1, theFile) != 1 ) {
	  printf("Error in FileChannel::sendVector, did not write\n");
	  delete []data;
	  return -1;
	}
	delete [] data;
	return 0;

}

int FileChannel::recvVector(int dbTag, int commitTag, 
		   Vector &theVector, 
			   ChannelAddress *theAddress)
{
    int  nleft, size;
    size = theVector.Size();
    double *data = new double[size];
    void *gMsg = (void *)data;
    nleft =  theVector.Size() * sizeof(double);

	if ( fread( gMsg, nleft, 1, theFile) != 1 ) {
	  printf("Error in FileChannel::recvVector, did not read\n");
      delete [] data;
	  return -1;
	}
	else {
	  theVector.setData (data, size);
	  return 0;
	}
}

int FileChannel::sendID(int dbTag, int commitTag, const ID &theID,ChannelAddress *theAddress)
{
  int i;
  int  nleft = theID.Size()*sizeof(int);    
  int *data = new int[theID.Size()];
  void *gMsg = (void *)data;

  /* I don't want to break encapsulation by being friends */
  for (i = 0; i < theID.Size(); i++ ) {
	data[i] = theID(i);
  }
 
  i =  fwrite( gMsg, nleft, 1, theFile);

  if ( i != 1 ) {
	printf("Error in FileChannel::sendID, did not write\n");
	delete []data;
	return -1;
  }
  delete [] data;
  return 0;
}

/**
 * Try to read data from file into newly allocated array and pass 
 * to ID
 * pre! old and new ID have the same data size (uncorrupted)
 * and same endianness
 */
int FileChannel::recvID(int dbTag, int commitTag, 
	       ID &theID,  ChannelAddress *theAddress) 
{

  int nleft,size,i;
  size = theID.Size();
  int *data = new int[size];
  void * gMsg = (void *)data;;
  nleft =  theID.Size() * sizeof(int);

  if( theFile ) {
  
	i =  fread( gMsg, nleft, 1, theFile);
 
	if ( i == 1 ) {
	  theID.setData( data, size );
	  return 0;
	}
	opserr << "FileChannel::recvID, error reading from open file\n";
  }
  else
	opserr << "FileChannel::recvID, error reading, NULL file handle\n";

  return -1;


}    


//OK< here, we can't afford not to be a friend of NDarray; it'd be really complicated to read.
int
FileChannel::sendnDarray(int dbTag, int commitTag, const nDarray &theNDarray,ChannelAddress *theAddress)
{
  int i,j;

  double *data = (theNDarray.pc_nDarray_rep)->pd_nDdata;
  int *dim =  (theNDarray.pc_nDarray_rep)->dim;
  int rank =  (theNDarray.pc_nDarray_rep)->nDarray_rank;
  int elem =   (theNDarray.pc_nDarray_rep)->total_numb;

  // send rank  
  if( fwrite( (void *)&rank, sizeof(int), 1, theFile ) != 1) {
	opserr << "FileChannel::sendnDarray, error writing rank to file\n";
	return -1;
  } 

  //send the length of the data array
  if ( fwrite( (void *)&elem, sizeof(int), 1, theFile ) != 1 ) {
  	opserr << 
	  "FileChannel::sendnDarray, error writing num elements to file\n";
	return -1;
  } 
  
  // send each dimension 
  for ( i = 0; i < rank; i++ ) {
    j = dim[i];
 
    if (  fwrite( (void *)&j, sizeof(int), 1, theFile ) != 1) {
	opserr << 
	  "FileChannel::sendnDarray, error writing array dimension to file\n";
	return -1;
    }
  }

  
   //send the data in the array
  if( fwrite( (void *)data, sizeof(double), elem, theFile) != elem ){
  	opserr << 
	  "FileChannel::sendnDarray, error writing array elements to file\n";
	return -1;
	} 
 
  return 0;
}


int
FileChannel::recvnDarray(int dbTag, int commitTag, nDarray &theNDarray, ChannelAddress *theAddress)
{
  int rank, elem, i, j;

  double db;
  double *n_data = (theNDarray.pc_nDarray_rep)->pd_nDdata;
  void *gmsg = (void *)n_data;

  // current size dimensions of the nDarray
  int *n_dim =  (theNDarray.pc_nDarray_rep)->dim;
  int n_rank =  (theNDarray.pc_nDarray_rep)->nDarray_rank;
  int n_elem =   (theNDarray.pc_nDarray_rep)->total_numb;


  // receive rank  
  if( fread( (void *)&rank, sizeof(int), 1, theFile ) != 1) {
	opserr << "FileChannel::receivenDarray, error reading rank to file\n";
	return -1;
  } 

  if ( rank != n_rank ) {
	opserr << "FileChannel::recvnDarray, mismatch in rank of nDArray\n";
	return -1;
  }

  //receive the length of the data array
  if ( fread( (void *)&elem, sizeof(int), 1, theFile ) != 1 ) {
  	opserr << 
	  "FileChannel::receivenDarray, error writing num elements to file\n";
	return -1;
  } 

  if( elem != n_elem )	{
	opserr << "FileChannel::recvnDarray, mismatch in length of nDArray\n";
	return -1;
  }

  // receive each dimension 
  for ( i = 0; i < rank; i++ ) {
	if (  fread( (void *)&j, sizeof(int), 1, theFile ) != 1) {
	opserr << 
	  "FileChannel::receivenDarray, error reading array dimension from file\n";
	return -1;
	}
	n_dim[i] = j;
  }
  
  //receive the data in the array
   if( fread( gmsg, sizeof(double), elem, theFile) != elem ){
    opserr << 
      "FileChannel::receivenDarray, error writing array elements to file\n";
    return -1;
    } 
 
 
  return 0;
}
