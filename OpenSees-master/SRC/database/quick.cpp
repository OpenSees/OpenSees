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
// $Date: 2003-02-14 23:00:54 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/quick.cpp,v $
                                                                        
                                                                        
int 
FileDatastore::sendID(int dataTag, int commitTag, 
		      const ID &theID, 
		      ChannelAddress *theAddress)
{

  // we first ensure that the ID is not too big
  int idSize = theID.Size();
  if (idSize >= maxIDsize) {
    opserr << "FileDatastore::sendID() - the database does not deal with IDs of this size ";
    opserr << idSize << endln;
  }


  // open a file if not already opened
  if (ids[idSize] == 0) {
    char fileName[70];
    char intName[10];
    strcpy(fileName, dataBase);
    itoa(idSize, intName);
    strcat(fileName,".IDs.");
    strcat(fileName,intName);
    ids[idSize] = this->openFile(fileName);    
    int loc = ids[idSize]->tellg();
    if (loc == -1) loc = 0;
    opserr << "LOCATION: " << loc << endln;
    fileEnds.ids[idSize] = loc;
  }

  // we now found the location in the file to write the data
  fstream *theStream = ids[idSize];
  int fileEnd = fileEnds.ids[idSize];
  int stepSize = (2 + idSize)*sizeof(int);

  theStream->seekg(0);
  bool found = false;
  int pos =0;

  while ((pos < fileEnd) && (found == false)) {
    theStream->read((char *)&idBuffer, stepSize);
    if ((idBuffer.dbTag == dataTag) && (idBuffer.commitTag == commitTag)) 
      found = true;
    else {
      pos += stepSize;
    }
  }

  // we now place the data to be sent into our buffer
  idBuffer.dbTag = dataTag;
  idBuffer.commitTag = commitTag;
  for (int i=0; i<idSize; i++)
    idBuffer.data[i] = theID(i);

  // we now write the data
  if (found == true)
    theStream->seekp(pos);

  theStream->write((char *)&idBuffer, stepSize);

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
    opserr << "FileDatastore::recvID() - the database does not deal with IDs";
    opserr << " of this size "  << idSize << endln;
    return -1;
  }
  
  // open the file - if not already opened
  if (ids[idSize] == 0) {
    char fileName[70];
    char intName[10];
    strcpy(fileName, dataBase);
    itoa(idSize, intName);
    strcat(fileName,".IDs.");
    strcat(fileName,intName);
    ids[idSize] = this->openFile(fileName);    
    int loc = ids[idSize]->tellg();
    if (loc == -1) loc = 0;
    fileEnds.ids[idSize] = loc;
  }

  // we now read in the data unti we reach eof or find the data
  int stepSize = (2 + idSize)*sizeof(int);
  fstream *theStream = ids[idSize];

  theStream->seekg(0);
  bool found = false;
  int pos =0; 
  int fileEnd = fileEnds.ids[idSize];

  while ((pos < fileEnd) && (found == false)) {
    theStream->read((char *)&idBuffer, stepSize);
    if ((idBuffer.dbTag == dataTag) && (idBuffer.commitTag == commitTag))
      found = true;
    pos += stepSize;
  }
  
  if (found == false) {
    opserr << "FileDatastore::recvID() - failed to find data for ID of size ";
    opserr << idSize << endln;
    return -1;
  }

  // we now place the received data into the ID 
  idBuffer.dbTag = dataTag;
  idBuffer.commitTag = commitTag;
  for (int i=0; i<idSize; i++)
    theID(i) = idBuffer.data[i];

  return 0;
}		       
