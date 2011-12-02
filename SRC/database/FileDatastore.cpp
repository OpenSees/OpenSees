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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/FileDatastore.cpp,v $
                                                                        
                                                                        
// File: ~/database/FileDatastore.C
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class implementation for FileDatastore.
// FileDatastore is a concrete subclas of FE_Datastore. A FileDatastore 
// object is used in the program to store/restore the geometry and state 
// information in a domain at a particular instance in the analysis. The
// information is stored in text files.
//
// What: "@(#) FileDatastore.C, revA"

#include "FileDatastore.h"
#include <iostream.h>
#include <fstream.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <bool.h>

#include <MovableObject.h>
#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

//void itoa(int x, char *str);

FileDatastore::FileDatastore(char *dataBaseName,
			     Domain &theDomain, 
			     FEM_ObjectBroker &theObjBroker) 
  :FE_Datastore(theDomain, theObjBroker), dbTag(0),
   ids(0), vects(0), mats(0)
{
  strcpy(dataBase, dataBaseName);

  ids = (fstream **)malloc(maxIDsize*sizeof(fstream *));
  vects = (fstream **)malloc(maxVectSize*sizeof(fstream *));
  mats = (fstream **)malloc(maxMatSize*sizeof(fstream *));

  if (ids == 0 || vects == 0 || mats == 0) {
    cerr << "FileDatastore::FileDatastore() - ran out of memory\n";
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
}

int
FileDatastore::getDbTag(void)
{
  dbTag++;
  return dbTag;
}

/********************************************************************
 *                   CHANNEL METHODS  THAT DO NOTHING               *
 ********************************************************************/

char *
FileDatastore::addToProgram(void)
{
  return 0;
}

int 
FileDatastore::setUpShadow(void)
{
  return 0;
}

int 
FileDatastore::setUpActor(void)
{
  return 0;
}

int 
FileDatastore::setNextAddress(const ChannelAddress &otherChannelAddress)
{
  return 0;
}


ChannelAddress *
FileDatastore::getLastSendersAddress(void)
{
  return 0;
}


int 
FileDatastore::commitState(int commitTag)
{
    int result = FE_Datastore::commitState(commitTag);

    // if sucessfull, we close the files to flush the buffers
    if (result == commitTag) {
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
    }

    return result;
}

int 
FileDatastore::sendObj(int commitTag,
		       MovableObject &theObject, 
		       ChannelAddress *theAddress)
{
  return theObject.sendSelf(commitTag, *this);
}

int 
FileDatastore::recvObj(int commitTag,
		       MovableObject &theObject, 
		       FEM_ObjectBroker &theNewBroker,
		       ChannelAddress *theAddress)
{
  return theObject.recvSelf(commitTag, *this, theNewBroker);
}

		
int 
FileDatastore::sendMsg(int dataTag, int commitTag, 
		       const Message &, 
		       ChannelAddress *theAddress)
{
  cerr << "FileDatastore::sendMsg() - not yet implemented\n";
  return -1;
}		       

int 
FileDatastore::recvMsg(int dataTag, int commitTag, 
		       Message &, 
		       ChannelAddress *theAddress)
{
  cerr << "FileDatastore::recvMsg() - not yet implemented\n";
  return -1;
}		       


int 
FileDatastore::sendMatrix(int dataTag, int commitTag, 
		      const Matrix &theMatrix, 
		      ChannelAddress *theAddress)
{

  // we first ensure that the Matrix is not too big
  int noMatCols= theMatrix.noCols();
  int noMatRows = theMatrix.noRows();
  int matSize = noMatRows * noMatCols;;

  if (matSize >= maxMatSize) {
    cerr << "FileDatastore::sendMatrix() - the database does not deal ";
    cerr << "with Matrixs of this size " << matSize << endl;
    return -1;
  }

  // open a file if not already opened
  if (mats[matSize] == 0) {
    char fileName[70];
    char intName[10];
    strcpy(fileName, dataBase);
	sprintf(intName,"%d",matSize);
    strcat(fileName,".Mats.");
    strcat(fileName,intName);
    mats[matSize] = this->openFile(fileName);    
    int loc = mats[matSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      int pos = 0;
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

  // we now found the location in the file to write the data
  fstream *theStream = mats[matSize];
  int fileEnd = fileEnds.mats[matSize];
  int stepSize = 2*sizeof(int) + matSize*sizeof(double);
  
  //
  // find location in file to place the data
  //

  bool found = false;  
  int pos = 0;
  
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
  if (found == true && filePos.mats[matSize] != pos)
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
  // we first check Matrix not too big
  int noMatCols= theMatrix.noCols();
  int noMatRows = theMatrix.noRows();
  int matSize = noMatRows * noMatCols;;

  if (matSize >= maxMatSize) {
    cerr << "FileDatastore::recvMatrix() - the database does not deal with Mats";
    cerr << " of this size " << matSize << endl;
    return -1;
  }

  // open the file if not already opened
  if (mats[matSize] == 0) {
    char fileName[70];
    char intName[10];
    strcpy(fileName, dataBase);
	sprintf(intName,"%d",matSize);
    //itoa(matSize, intName);
    strcat(fileName,".Mats.");
    strcat(fileName,intName);
    mats[matSize] = this->openFile(fileName);    
    int loc = mats[matSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      int pos = 0;
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
  int pos = filePos.mats[matSize];
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

  // we first ensure that the ID is not too big
  int vectSize = theVector.Size();
  if (vectSize >= maxVectSize) {
    cerr << "FileDatastore::sendVector() - the database does not deal ";
    cerr << "with Vectors of this size " << vectSize << endl;
  }

  // open a file if not already opened 
  if (vects[vectSize] == 0) {
    char fileName[70];
    char intName[10];
    strcpy(fileName, dataBase);
	sprintf(intName,"%d",vectSize);
//    itoa(vectSize, intName);
    strcat(fileName,".Vects.");
    strcat(fileName,intName);
    vects[vectSize] = this->openFile(fileName);    
    int loc = vects[vectSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      int pos = 0;
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
  int pos = 0;
  
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
  if (found == true && pos != filePos.vects[vectSize])
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
  // we first check ID not too big
  int vectSize = theVector.Size();
  if (vectSize >= maxVectSize) {
    cerr << "FileDatastore::recvVector() - the database does not deal with Vects";
    cerr << " of this size " << vectSize << endl;
    return -1;
  }
  

  if (vects[vectSize] == 0) {
    char fileName[70];
    char intName[10];
    strcpy(fileName, dataBase);
		sprintf(intName,"%d",vectSize);
 //   itoa(vectSize, intName);
    strcat(fileName,".Vects.");
    strcat(fileName,intName);
    vects[vectSize] = this->openFile(fileName);    
    int loc = vects[vectSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      int pos = 0;
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
  int pos = filePos.vects[vectSize];
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
  // we first ensure that the ID is not too big
  int idSize = theID.Size();
  if (idSize >= maxIDsize) {
    cerr << "FileDatastore::sendID() - the database does not deal with IDs of this size ";
    cerr << idSize << endl;
  }


  // open a file if not already opened
  if (ids[idSize] == 0) {

    char fileName[70];
    char intName[10];
    strcpy(fileName, dataBase);
	sprintf(intName,"%d",idSize);
//    itoa(idSize, intName);
    strcat(fileName,".IDs.");
    strcat(fileName,intName);

    ids[idSize] = this->openFile(fileName);  
	
    int loc = ids[idSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      int pos = 0;
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
  int fileEnd = fileEnds.ids[idSize];
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
  if (found == true && pos != filePos.ids[idSize]) 
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
  // we first check ID not too big
  int idSize = theID.Size();
  if (idSize >= maxIDsize) {
    cerr << "FileDatastore::recvID() - the database does not deal with IDs";
    cerr << " of this size "  << idSize << endl;
    return -1;
  }

  // open file if not already done so
  if (ids[idSize] == 0) {
    char fileName[70];
    char intName[10];
    strcpy(fileName, dataBase);
	sprintf(intName,"%d",idSize);
//    itoa(idSize, intName);
    strcat(fileName,".IDs.");
    strcat(fileName,intName);
    ids[idSize] = this->openFile(fileName);    
    int loc = ids[idSize]->tellg();
    if (loc == -1) 
      loc = 0;
    else {
      int pos = 0;
      int stepSize = (2 + idSize)*sizeof(int);
      fstream *theStream = ids[idSize];
      theStream->seekg(0);  
      int maxDataTag = 0;
      int maxCommitTag = 0;
      while (pos < loc) {
	  theStream->read((char *)&idBuffer, stepSize);
	  /*
	  cerr << idBuffer.dbTag << " " << idBuffer.commitTag;
	  for (int i=0; i<idSize; i++)
	    cerr << " " << idBuffer.data[i];
	  cerr << endl;
	  */
	  
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
  int fileEnd = fileEnds.ids[idSize];  
  int pos = filePos.ids[idSize];
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


/*******************************************************************
 *              MISC METHODS & FUNCTONS FOR OPENING THE FILE       *
 *******************************************************************/

fstream *
FileDatastore::openFile(char *fileName)
{
#ifdef _WIN32    
    fstream *res = new fstream(fileName, ios::in | ios::out | ios::binary);
#else
    fstream *res = new fstream(fileName, ios::in | ios::out);    
#endif
    
    if (res == 0) {
	cerr << "FATAL - FileDatastore::openFile() - could not open file ";
	cerr << fileName << endl;
	exit(-1);
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
