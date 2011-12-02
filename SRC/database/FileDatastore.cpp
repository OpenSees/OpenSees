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
                                                                        
// $Revision: 1.10 $
// $Date: 2005-07-15 00:02:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/FileDatastore.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 10/98
//
// Description: This file contains the class implementation for FileDatastore.
// FileDatastore is a concrete subclas of FE_Datastore. A FileDatastore 
// object is used in the program to store/restore the geometry and state 
// information in a domain at a particular instance in the analysis. The
// information is stored in text files.
//
// What: "@(#) FileDatastore.C, revA"

#include "FileDatastore.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <bool.h>
#include <iostream>
#include <iomanip>
using std::ios;
using std::setiosflags;

#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

//void itoa(int x, char *str);

FileDatastore::FileDatastore(const char *dataBaseName,
			     Domain &theDomain, 
			     FEM_ObjectBroker &theObjBroker) 
  :FE_Datastore(theDomain, theObjBroker), dbTag(0),
   ids(0), vects(0), mats(0), currentCommitTag(-1)
{
  dataBase = new char [strlen(dataBaseName)+1];
  strcpy(dataBase, dataBaseName);

  ids = (fstream **)malloc(maxIDsize*sizeof(fstream *));
  vects = (fstream **)malloc(maxVectSize*sizeof(fstream *));
  mats = (fstream **)malloc(maxMatSize*sizeof(fstream *));

  if (ids == 0 || vects == 0 || mats == 0) {
    opserr << "FileDatastore::FileDatastore() - ran out of memory\n";
    exit(-1);
  }

  for (int i=0; i<maxIDsize; i++) {
    ids[i] = 0;
    fileEnds.ids[i] = 0;    
    filePos.ids[i] = 0;    
    fileCommitTags.ids[i] = 0;    
    fileDbTags.ids[i] = 0;        
  }
  for (int j=0; j<maxVectSize; j++) {
    vects[j] = 0;
    fileEnds.vects[j] = 0;
    filePos.vects[j] = 0;
    fileCommitTags.vects[j] = 0;
    fileDbTags.vects[j] = 0;    
  }
  for (int k=0; k<maxMatSize; k++) {
    mats[k] = 0;
    fileEnds.mats[k] = 0;
    filePos.mats[k] = 0;
    fileCommitTags.mats[k] = 0;
    fileDbTags.mats[k] = 0;    
  }
}

FileDatastore::~FileDatastore() 
{

  for (int i=0; i<maxIDsize; i++) 
    if (ids[i] != 0) ids[i]->close(); 
  for (int j=0; j<maxVectSize; j++)
    if (vects[j] != 0) vects[j]->close();
  for (int k=0; k<maxMatSize; k++)
    if (mats[k] != 0) mats[k]->close();
  
  free ((void *) ids);
  free ((void *) vects);
  free ((void *) mats);

  if (dataBase != 0)
    delete [] dataBase;
}

int
FileDatastore::getDbTag(void)
{
  dbTag++;
  return dbTag;
}

int 
FileDatastore::commitState(int commitTag)
{
  int result = FE_Datastore::commitState(commitTag);
  if (result == commitTag) 
    resetFilePointers();
    
  return result;
}

void
FileDatastore::resetFilePointers(void) {
  
  for (int i=0; i<maxIDsize; i++) 
    if (ids[i] != 0) {
      ids[i]->close(); 
      ids[i] = 0;
    }
  for (int j=0; j<maxVectSize; j++)
    if (vects[j] != 0) {
      vects[j]->close();
      vects[j] = 0;
    }
  for (int k=0; k<maxMatSize; k++)
    if (mats[k] != 0) {
      mats[k]->close();	
      mats[k] = 0;
    }

  currentCommitTag = -1;
}

int 
FileDatastore::sendMsg(int dataTag, int commitTag, 
		       const Message &, 
		       ChannelAddress *theAddress)
{
  opserr << "FileDatastore::sendMsg() - not yet implemented\n";
  return -1;
}		       

int 
FileDatastore::recvMsg(int dataTag, int commitTag, 
		       Message &, 
		       ChannelAddress *theAddress)
{
  opserr << "FileDatastore::recvMsg() - not yet implemented\n";
  return -1;
}		       


int 
FileDatastore::sendMatrix(int dataTag, int commitTag, 
		      const Matrix &theMatrix, 
		      ChannelAddress *theAddress)
{

  if (currentCommitTag != commitTag)
    this->resetFilePointers();
  currentCommitTag = commitTag;

  // we first ensure that the Matrix is not too big
  int noMatCols= theMatrix.noCols();
  int noMatRows = theMatrix.noRows();
  int matSize = noMatRows * noMatCols;;

  if (matSize >= maxMatSize) {
    opserr << "FileDatastore::sendMatrix() - the database does not deal ";
    opserr << "with Matrixs of this size " << matSize << endln;
    return -1;
  }

  // open a file if not already opened
  if (mats[matSize] == 0) {
    char *fileName = new char[strlen(dataBase)+21];
    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",matSize,commitTag);
    strcat(fileName,".Mats.");
    strcat(fileName,intName);
    

    mats[matSize] = this->openFile(fileName);    
    if (mats[matSize] == 0) {
      opserr << "FileDatastore::sendMatrix() - could not open file\n";
      return -1;
    }
    
    long int loc = mats[matSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      long int pos = 0;
      int stepSize = 2*sizeof(int) + matSize*sizeof(double);
      fstream *theStream = mats[matSize];
      theStream->seekg(0);
      int maxDataTag = 0;
      int maxCommitTag = 0;
      while (pos < loc) {
	  theStream->read((char *)&matBuffer, stepSize);
	  if ((matBuffer.dbTag >= dataTag) && 
	      (matBuffer.commitTag >= commitTag)) {
	    maxDataTag = matBuffer.dbTag;
	    maxCommitTag = matBuffer.commitTag;	    
	  }
	  pos += stepSize;
      }
      fileCommitTags.mats[matSize] = maxCommitTag;
      fileDbTags.mats[matSize] = maxDataTag; 	        
      filePos.mats[matSize] = loc;
    }
    fileEnds.mats[matSize] = loc;
    delete [] fileName;
  }

  // we now found the location in the file to write the data
  fstream *theStream = mats[matSize];
  int fileEnd = fileEnds.mats[matSize];
  int stepSize = 2*sizeof(int) + matSize*sizeof(double);
  
  //
  // find location in file to place the data
  //

  bool found = false;  
  long int pos = 0;
  
  // we first check if the data can go at the end of the file
  // true if commitTag larger than any we have encountered so far
  if (fileDbTags.mats[matSize] < dataTag)  {
      pos = fileEnd;
      found = true;
      
      fileCommitTags.mats[matSize] = commitTag;
      fileDbTags.mats[matSize] = dataTag; 	  

  // we can also add at end if current commitTag as last time we added
  // and the dbTag is larger than anything entered yet so far
  } else if ((fileDbTags.mats[matSize] == dataTag) &&
	     (fileCommitTags.mats[matSize] < commitTag)) {
      
      pos = fileEnd;
      found = true;
      fileCommitTags.mats[matSize] = commitTag;      

  // we have to search from the beginning of the file
  } else { 
      theStream->seekg(0);
      while ((pos < fileEnd) && (found == false)) {
	  theStream->read((char *)&matBuffer, stepSize);
	  if ((matBuffer.dbTag == dataTag) && 
	      (matBuffer.commitTag == commitTag)) 

	      found = true;
	  else
	    pos += stepSize;
      }
      if (found == true)
	filePos.mats[matSize] = pos + stepSize;
      else
	filePos.mats[matSize] = pos;
  }
  
  // we now place the data to be sent into our buffer
  matBuffer.dbTag = dataTag;
  matBuffer.commitTag = commitTag;
  int i=0;
  for (int j=0; j<noMatCols; j++)
    for (int k=0; k < noMatRows; k++) {
      matBuffer.data[i] = theMatrix(k,j);
      i++;
    }


  // we now write the data
  //if (found == true && filePos.mats[matSize] != pos)
    theStream->seekp(pos);

  theStream->write((char *)&matBuffer, stepSize);

  // update the size of file if we have added to eof
  if (fileEnd <= pos)
      fileEnds.mats[matSize] += stepSize;
  
  filePos.mats[matSize] = pos + stepSize;
  
  return 0;
}		       




int 
FileDatastore::recvMatrix(int dataTag, int commitTag, 
		      Matrix &theMatrix, 
		      ChannelAddress *theAddress)    
{

  if (currentCommitTag != commitTag)
    this->resetFilePointers();
  currentCommitTag = commitTag;

  // we first check Matrix not too big
  int noMatCols= theMatrix.noCols();
  int noMatRows = theMatrix.noRows();
  int matSize = noMatRows * noMatCols;;

  if (matSize >= maxMatSize) {
    opserr << "FileDatastore::recvMatrix() - the database does not deal with Mats";
    opserr << " of this size " << matSize << endln;
    return -1;
  }

  // open the file if not already opened
  if (mats[matSize] == 0) {
    char *fileName = new char[strlen(dataBase)+21];
    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",matSize,commitTag);
    strcat(fileName,".Mats.");
    strcat(fileName,intName);
    mats[matSize] = this->openFile(fileName);    
    if (mats[matSize] == 0) {
      opserr << "FileDatastore::recvMatrix() - could not open file\n";
      return -1;
    }
    long int loc = mats[matSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      long int pos = 0;
      int stepSize = 2*sizeof(int) + matSize*sizeof(double);
      fstream *theStream = mats[matSize];
      theStream->seekg(0);
      int maxDataTag = 0;
      int maxCommitTag = 0;
      while (pos < loc) {
	  theStream->read((char *)&matBuffer, stepSize);
	  if ((matBuffer.dbTag >= dataTag) && 
	      (matBuffer.commitTag >= commitTag)) {
	    maxDataTag = matBuffer.dbTag;
	    maxCommitTag = matBuffer.commitTag;	    
	  }
	  pos += stepSize;
      }
      fileCommitTags.mats[matSize] = maxCommitTag;
      fileDbTags.mats[matSize] = maxDataTag; 	        
      filePos.mats[matSize] = loc;
    }
    fileEnds.mats[matSize] = loc;
  }
  
  // we now read in the data unti we reach eof or find the data

  int stepSize = 2*sizeof(int) + matSize*sizeof(double);
  fstream *theStream = mats[matSize];

  int fileEnd = fileEnds.mats[matSize];  
  long int pos = filePos.mats[matSize];
  bool found = false;
  
  // we try the current file position first
  if (pos < fileEnd) {
      theStream->read((char *)&matBuffer, stepSize);
      if ((matBuffer.dbTag == dataTag) && (matBuffer.commitTag == commitTag)) {
	  found = true;
	  filePos.mats[matSize] += stepSize;
      } 
  }

  // we must search from the begiing of the file
  if (found == false) {
      theStream->seekg(0);
      pos =0; 
      while ((pos < fileEnd) && (found == false)) {
	  theStream->read((char *)&matBuffer, stepSize);
	  if ((matBuffer.dbTag == dataTag) && (matBuffer.commitTag == commitTag))
	      found = true;
	  pos += stepSize;
      }

      filePos.mats[matSize] = pos;      
  }
  
  if (found == false) {
    return -1;
  }

  // we now place the received data into the ID 
  matBuffer.dbTag = dataTag;
  matBuffer.commitTag = commitTag;

  int i=0;
  for (int j=0; j<noMatCols; j++)
    for (int k=0; k < noMatRows; k++) {
      theMatrix(k,j) = matBuffer.data[i];
      i++;
    }

  return 0;
}		       



int 
FileDatastore::sendVector(int dataTag, int commitTag, 
		      const Vector &theVector, 
		      ChannelAddress *theAddress)
{
  if (currentCommitTag != commitTag)
    this->resetFilePointers();
  currentCommitTag = commitTag;

  // we first ensure that the ID is not too big
  int vectSize = theVector.Size();
  if (vectSize >= maxVectSize) {
    opserr << "FileDatastore::sendVector() - the database does not deal ";
    opserr << "with Vectors of this size " << vectSize << endln;
  }

  // open a file if not already opened 
  if (vects[vectSize] == 0) {
    char *fileName = new char[strlen(dataBase)+21];
    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",vectSize,commitTag);
    strcat(fileName,".Vects.");
    strcat(fileName,intName);
    vects[vectSize] = this->openFile(fileName);    
    if (vects[vectSize] == 0) {
      opserr << "FileDatastore::sendVector() - could not open file\n";
      return -1;
    }
    long int loc = vects[vectSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      long int pos = 0;
      int stepSize = 2*sizeof(int) + vectSize*sizeof(double);
      fstream *theStream = vects[vectSize];
      theStream->seekg(0);
      int maxDataTag = 0;
      int maxCommitTag = 0;
      while (pos < loc) {
	  theStream->read((char *)&vectBuffer, stepSize);
	  if ((vectBuffer.dbTag >= dataTag) && 
	      (vectBuffer.commitTag >= commitTag)) {
	    maxDataTag = vectBuffer.dbTag;
	    maxCommitTag = vectBuffer.commitTag;	    
	  }
	  pos += stepSize;
      }
      fileCommitTags.vects[vectSize] = maxCommitTag;
      fileDbTags.vects[vectSize] = maxDataTag; 	        
      filePos.vects[vectSize] = loc;
    }
    fileEnds.vects[vectSize] = loc;
  }

  fstream *theStream = vects[vectSize];
  int fileEnd = fileEnds.vects[vectSize];
  int stepSize = 2*sizeof(int) + vectSize*sizeof(double);
  
  //
  // find location in file to place the data
  //

  bool found = false;  
  long int pos = 0;
  
  // we first check if the data can go at the end of the file
  // true if commitTag larger than any we have encountered so far
  if (fileDbTags.vects[vectSize] < dataTag)  {
      pos = fileEnd;
      found = true;
      
      fileCommitTags.vects[vectSize] = commitTag;
      fileDbTags.vects[vectSize] = dataTag; 	  

  // we can also add at end if current commitTag as last time we added
  // and the dbTag is larger than anything entered yet so far
  } else if ((fileDbTags.vects[vectSize] == dataTag) &&
	     (fileCommitTags.vects[vectSize] < commitTag)) {
      
      pos = fileEnd;
      found = true;
      fileCommitTags.vects[vectSize] = commitTag;      

  // we have to search from the beginning of the file
  } else { 
      theStream->seekg(0);
      while ((pos < fileEnd) && (found == false)) {
	  theStream->read((char *)&vectBuffer, stepSize);
	  if ((vectBuffer.dbTag == dataTag) && 
	      (vectBuffer.commitTag == commitTag)) 

	      found = true;
	  else
	    pos += stepSize;
      }
      if (found == true)
	filePos.vects[vectSize] = pos + stepSize; 
      else
	filePos.vects[vectSize] = pos; 
  }

  
  // we now place the data to be sent into our buffer
  vectBuffer.dbTag = dataTag;
  vectBuffer.commitTag = commitTag;
  for (int i=0; i<vectSize; i++)
    vectBuffer.data[i] = theVector(i);

  // we now write the data
  // if (found == true && pos != filePos.vects[vectSize])
    theStream->seekp(pos);

  theStream->write((char *)&vectBuffer, stepSize);

  filePos.vects[vectSize] = pos + stepSize;        
  
  // update the size of file if we have added to eof
  if (fileEnd <= pos)
    fileEnds.vects[vectSize] += stepSize;

  return 0;
}		       

int 
FileDatastore::recvVector(int dataTag, int commitTag, 
		      Vector &theVector, 
		      ChannelAddress *theAddress)    
{
  if (currentCommitTag != commitTag)
    this->resetFilePointers();
  currentCommitTag = commitTag;

  // we first check ID not too big
  int vectSize = theVector.Size();
  if (vectSize >= maxVectSize) {
    opserr << "FileDatastore::recvVector() - the database does not deal with Vects";
    opserr << " of this size " << vectSize << endln;
    return -1;
  }
  

  if (vects[vectSize] == 0) {
    char *fileName = new char[strlen(dataBase)+21];
    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",vectSize,commitTag);
    strcat(fileName,".Vects.");
    strcat(fileName,intName);
    vects[vectSize] = this->openFile(fileName);    
    if (vects[vectSize] == 0) {
      opserr << "FileDatastore::recvVector() - could not open file\n";
      return -1;
    }
    long int loc = vects[vectSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      long int pos = 0;
      int stepSize = 2*sizeof(int) + vectSize*sizeof(double);
      fstream *theStream = vects[vectSize];
      theStream->seekg(0);  
      int maxDataTag = 0;
      int maxCommitTag = 0;
      while (pos < loc) {
	  theStream->read((char *)&vectBuffer, stepSize);
	  if ((vectBuffer.dbTag >= dataTag) && 
	      (vectBuffer.commitTag >= commitTag)) {
	    maxDataTag = vectBuffer.dbTag;
	    maxCommitTag = vectBuffer.commitTag;	    
	  }
	  pos += stepSize;
      }
      fileCommitTags.vects[vectSize] = maxCommitTag;
      fileDbTags.vects[vectSize] = maxDataTag; 	        
      filePos.vects[vectSize] = loc;
    }
    fileEnds.vects[vectSize] = loc;
  }

  // we now read in the data unti we reach eof or find the data
  int stepSize = 2*sizeof(int) + vectSize*sizeof(double);
  fstream *theStream = vects[vectSize];

  int fileEnd = fileEnds.vects[vectSize];  
  long int pos = filePos.vects[vectSize];
  bool found = false;

  // we try the current file position first
  if (pos < fileEnd) {
      theStream->read((char *)&vectBuffer, stepSize);
      if ((vectBuffer.dbTag == dataTag) && (vectBuffer.commitTag == commitTag)) {
	  found = true;
	  filePos.vects[vectSize] += stepSize;
      } 
  }

  // we must search from the begiing of the file
  if (found == false) {
      theStream->seekg(0);  
      pos = 0;
      while ((pos < fileEnd) && (found == false)) {
	  theStream->read((char *)&vectBuffer, stepSize);
	  if ((vectBuffer.dbTag == dataTag) && 
	      (vectBuffer.commitTag == commitTag)) 

	    found = true;
	  pos += stepSize;
      }
      filePos.vects[vectSize] = pos;
  }
  
  if (found == false) {
    return -1;
  }

  // we now place the received data into the Vector
  vectBuffer.dbTag = dataTag;
  vectBuffer.commitTag = commitTag;
  for (int i=0; i<vectSize; i++)
    theVector(i) = vectBuffer.data[i];

  return 0;
}		       



int 
FileDatastore::sendID(int dataTag, int commitTag, 
		      const ID &theID, 
		      ChannelAddress *theAddress)
{
  if (currentCommitTag != commitTag)
    this->resetFilePointers();
  currentCommitTag = commitTag;

  // we first ensure that the ID is not too big
  int idSize = theID.Size();
  if (idSize >= maxIDsize) {
    opserr << "FileDatastore::sendID() - the database does not deal with IDs of this size ";
    opserr << idSize << endln;
  }

  // open a file if not already opened
  if (ids[idSize] == 0) {
    char *fileName = new char[strlen(dataBase)+21];
    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",idSize,commitTag);
    strcat(fileName,".IDs.");
    strcat(fileName,intName);

    ids[idSize] = this->openFile(fileName);  
    if (ids[idSize] == 0) {
      opserr << "FileDatastore::sendID() - could not open file\n";
      return -1;
    }
	
    long int loc = ids[idSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      long int pos = 0;
      int stepSize = (2 + idSize)*sizeof(int);
      fstream *theStream = ids[idSize];
      theStream->seekg(0);  
      int maxDataTag = 0;
      int maxCommitTag = 0;
      while (pos < loc) {
	  theStream->read((char *)&idBuffer, stepSize);
	  if ((idBuffer.dbTag >= dataTag) && 
	      (idBuffer.commitTag >= commitTag)) {
	    maxDataTag = idBuffer.dbTag;
	    maxCommitTag = idBuffer.commitTag;	    
	  }
	  pos += stepSize;
      }
      fileCommitTags.ids[idSize] = maxCommitTag;
      fileDbTags.ids[idSize] = maxDataTag; 	        
      filePos.ids[idSize] = loc;
    }
    fileEnds.ids[idSize] = loc;
  }


  // we now found the location in the file to write the data
  fstream *theStream = ids[idSize];
  long int fileEnd = fileEnds.ids[idSize];
  int stepSize = (2 + idSize)*sizeof(int);

  //
  // find location in file to place the data
  //

  bool found = false;  
  int pos = 0;

  // we first check if the data can go at the end of the file
  // true if commitTag larger than any we have encountered so far
  if (fileDbTags.ids[idSize] < dataTag)  {

      pos = fileEnd;
      found = true;
      fileCommitTags.ids[idSize] = commitTag;
      fileDbTags.ids[idSize] = dataTag; 	  

  // we can also add at end if current commitTag as last time we added
  // and the dbTag is larger than anything entered yet so far
  } else if ((fileDbTags.ids[idSize] == dataTag) &&
	     (fileCommitTags.ids[idSize] < commitTag)) {
      
      pos = fileEnd;
      found = true;
      fileCommitTags.ids[idSize] = commitTag;      

  // we have to search from the beginning of the file
  } else { 
      theStream->seekg(0);
      while ((pos < fileEnd) && (found == false)) {
	  theStream->read((char *)&idBuffer, stepSize);
	  if ((idBuffer.dbTag == dataTag) && 
	      (idBuffer.commitTag == commitTag)) 

	      found = true;
	  else
	    pos += stepSize;
      }
      if (found == true)
	filePos.ids[idSize] = pos + stepSize; 
      else
	filePos.ids[idSize] = pos; 

  }


  // we now place the data to be sent into our buffer
  idBuffer.dbTag = dataTag;
  idBuffer.commitTag = commitTag;
  for (int i=0; i<idSize; i++)
    idBuffer.data[i] = theID(i);
  // we now write the data
  //if (found == true && pos != filePos.ids[idSize]) 
    theStream->seekp(pos);

  theStream->write((char *)&idBuffer, stepSize);

  filePos.ids[idSize] = pos + stepSize;        
  
  // update the size of file if we have added to eof
  if (fileEnd <= pos)
    fileEnds.ids[idSize] += stepSize;  

  return 0;
}		       

int 
FileDatastore::recvID(int dataTag, int commitTag, 
		      ID &theID, 
		      ChannelAddress *theAddress)    
{
  if (currentCommitTag != commitTag)
    this->resetFilePointers();
  currentCommitTag = commitTag;

  // we first check ID not too big
  int idSize = theID.Size();
  if (idSize >= maxIDsize) {
    opserr << "FileDatastore::recvID() - the database does not deal with IDs";
    opserr << " of this size "  << idSize << endln;
    return -1;
  }

  // open file if not already done so
  if (ids[idSize] == 0) {
    char *fileName = new char[strlen(dataBase)+21];
    char intName[20];
    strcpy(fileName, dataBase);
    sprintf(intName,"%d.%d",idSize,commitTag);
    strcat(fileName,".IDs.");
    strcat(fileName,intName);

    ids[idSize] = this->openFile(fileName);    
    if (ids[idSize] == 0) {
      opserr << "FileDatastore::recvID() - could not open file\n";
      return -1;
    }

    long int loc = ids[idSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      long int pos = 0;
      int stepSize = (2 + idSize)*sizeof(int);
      fstream *theStream = ids[idSize];
      theStream->seekg(0);  
      int maxDataTag = 0;
      int maxCommitTag = 0;
      while (pos < loc) {
	  theStream->read((char *)&idBuffer, stepSize);
	  
	  if ((idBuffer.dbTag >= dataTag) && 
	      (idBuffer.commitTag >= commitTag)) {
	    maxDataTag = idBuffer.dbTag;
	    maxCommitTag = idBuffer.commitTag;	    
	  }
	  pos += stepSize;
      }
      fileCommitTags.ids[idSize] = maxCommitTag;
      fileDbTags.ids[idSize] = maxDataTag; 	        
      filePos.ids[idSize] = loc;
    }
    fileEnds.ids[idSize] = loc;
  }

  // we now set some parameters before we go looking for the data
  int stepSize = (2 + idSize)*sizeof(int);
  fstream *theStream = ids[idSize];
  long int fileEnd = fileEnds.ids[idSize];  
  long int pos = filePos.ids[idSize];
  bool found = false;

  // we try the current file position first
  if (pos < fileEnd) {
      theStream->read((char *)&idBuffer, stepSize);

      if ((idBuffer.dbTag == dataTag) && (idBuffer.commitTag == commitTag)) {
	  found = true;
	  filePos.ids[idSize] += stepSize;
      } 
  }

  // we must search from the beginning of the file
  if (found == false) {
      pos = 0;
      theStream->seekg(0);  

      while ((pos < fileEnd) && (found == false)) {
	  theStream->read((char *)&idBuffer, stepSize);

	  if ((idBuffer.dbTag == dataTag) && 
	      (idBuffer.commitTag == commitTag)) {

	      found = true;
	  }

	  pos += stepSize;
      }
      filePos.ids[idSize] = pos;
  }  

  if (found == false) {
    return -1;
  }

  // we now place the received data into the ID 
  idBuffer.dbTag = dataTag;
  idBuffer.commitTag = commitTag;
  for (int i=0; i<idSize; i++)
    theID(i) = idBuffer.data[i];

  return 0;
}		       


#include <fstream>
using std::ofstream;
using std::cerr;


int 
FileDatastore::createTable(const char *tableName, int numColumns, char *columns[])
{
  // open the file
  int res = 0;
  char *fileName = new char[strlen(tableName) + strlen(dataBase) + 1];
  if (fileName == 0) {
    opserr << "FileDatastore::insertData - out of memory; failed to open file: " << fileName << endln;
    return -1;
  }

  strcpy(fileName, dataBase);    
  strcat(fileName,".");
  strcat(fileName, tableName);

  ofstream table;
  table.open(fileName, ios::out | ios::trunc); 

  if (table.is_open()) {
    // write the data
    for (int i=0; i<numColumns; i++) {
      table << columns[i] << "\t";
    }
    table << "\n";
    table.close();
    
  } else {
    opserr << "FileDatastore::insertData - failed to open file: " << fileName << endln;
    res = -1;
  }
  
  delete [] fileName;
  return res;
}

int 
FileDatastore::insertData(const char *tableName, char *columns[], 
			  int commitTag, const Vector &data)
{
  // open the file
  char *fileName = new char[strlen(tableName) + strlen(dataBase) + 1];
  if (fileName == 0) {
    opserr << "FileDatastore::insertData - out of memory; failed to open file: " << fileName << endln;
    return -1;
  }

  strcpy(fileName, dataBase);    
  strcat(fileName,".");
  strcat(fileName, tableName);

  ofstream table;
  table.open(fileName, ios::app); 

  table << setiosflags(ios::scientific);
  table << std::setprecision(16);

  if (table.is_open()) {
    // write the data
    for (int i=0; i<data.Size(); i++) {
      table << data(i) << "\t";
    }
    
    table << "\n";
    table.close();

  } else {
    opserr << "FileDatastore::insertData - failed to open file: " << fileName << endln;
    return -1;
  }

  delete [] fileName;
  return 0;
}


int 
FileDatastore::getData(const char *tableName, char *columns[], int commitTag, Vector &data)
{
  return 0;
}



/*******************************************************************
 *              MISC METHODS & FUNCTONS FOR OPENING THE FILE       *
 *******************************************************************/

fstream *
FileDatastore::openFile(char *fileName)
{
  
  fstream *res = new fstream();
  if (res == 0) {
    opserr << "FileDatastore::openFile - out of memory; failed to open file: " << fileName << endln;
    return 0;
  }

  res->open(fileName, ios::in | ios::out | ios::binary); 

  // if file did not exist, need to pass trunc flag to open it
  if (res->bad() == true || res->is_open() == false) {
    // delete & new again for unix gcc compiler to work!
    delete res;
    res = new fstream();
    if (res == 0) {
      opserr << "FileDatastore::openFile - out of memory; failed to open file: " << fileName << endln;
      return 0;
    }
    res->open(fileName, ios::in | ios::out | ios::trunc | ios::binary);   
  }

  if (res->bad() == true || res->is_open() == false) {
    opserr << "FATAL - FileDatastore::openFile() - could not open file " << fileName << endln;
    delete res;
    return 0;
  }

  // set the position for writing to eof
  res->seekp(0,ios::end);  
  return res;
}

/*
char itoc(int x)
{
 if (x == 1) return '1';
 if (x == 2) return '2';
 if (x == 3) return '3';
 if (x == 4) return '4';
 if (x == 5) return '5';
 if (x == 6) return '6';
 if (x == 7) return '7';
 if (x == 8) return '8';
 if (x == 9) return '9';
 return '0';
}

void
itoa(int x, char *str)
{
  int y=x;
  while (y >= 10) 
    y = y/10;
  str[0] = itoc(y);
  str[1] = '\0';
  if (x >= 10) {
    int z = x/10;
    z = x - 10*z;
    itoa(z,&str[1]);
  }
}

*/
