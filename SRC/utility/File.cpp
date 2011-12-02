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
                                                                        
// $Revision: 1.9 $
// $Date: 2009-03-27 19:18:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/File.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 09/07
//
// Description: This file contains the class definition for File used in SimulationINformation.


#include <map>
#include <string.h>
//#include <string>
//using namespace std;



#include <File.h>
#include <FileIter.h>

#include <OPS_Globals.h>

File::File(const char *theName, const char *theDescription, bool isDir)
  :isDirectory(isDir), parentDir(0), name(theName), description(theDescription)
{
  
}

File::~File()
{
  if (isDirectory == true) {

    FileIter theDirFiles = this->getFiles();
    File *aDirFile;
    while ((aDirFile = theDirFiles()) != 0)
      delete aDirFile;
  } 
  dirFiles.clear();
}

int
File::clear(void)
{
  if (isDirectory == true) {

    FileIter theDirFiles = this->getFiles();
    File *aDirFile =0;
    while ((aDirFile = theDirFiles()) != 0)
      delete aDirFile;
  } 

  dirFiles.clear();
  name.clear();
  description.clear();

  return 0;
}

int 
File::addFile(File *theFile)
{
  if (isDirectory == false)
    return -1;
  
  if (dirFiles.find(theFile->name) == dirFiles.end()) {
    dirFiles[theFile->name] = theFile;
  } else
    return -1;

  theFile->setParentDir(this);
  return 0;
}

int
File::addFile(const char *fileName, const char *path, const char *fileDescription)
{
  static char dirName[128];
  
  char *combined = 0;     // combined array of path + fileName
  const char *combinedFile = 0; // pointer into combined where file is

  if (path != 0 && strncmp(fileName,"/", 1) != 0) {
    combined = new char[strlen(fileName)+strlen(path)+2];
    strcpy(combined, path);
    strcat(combined, "/");
    strcat(combined, fileName);
  } else {
    combined = new char[strlen(fileName)+1];
    strcpy(combined, fileName);
  }
  

  File *currentDir = this;

  const char *pathCurrent = strstr(combined, "/");

  if (pathCurrent != NULL) {

    if (description.length() == 0) {

      int rootPathLength = pathCurrent-combined;
      char *rootPath = new char[rootPathLength+2];
     
      if (rootPathLength > 0) {
		strncpy(rootPath, combined, rootPathLength);
		strcpy(&rootPath[rootPathLength],"/");
	
      } else
		strcpy(rootPath, "/");

      description.assign(rootPath);      
      delete [] rootPath;
    } 
    
    pathCurrent+=1;      
  
    while (pathCurrent != NULL) {
      const char *prev = pathCurrent;
      const char *next = strstr(prev,"/");

      int dirNameLength = 0;
      if (next == 0) {
	combinedFile = prev;
	pathCurrent = 0;
      } else {
	dirNameLength = next-prev;
	strncpy(dirName, prev, dirNameLength);
	strcpy(&dirName[dirNameLength],"");
	pathCurrent = next+1;

	File *nextDir = currentDir->getFile(dirName);
	if (nextDir == 0) {
	
	  const char *prevPath = currentDir->getDescription();
	  char *newPath = 0;
	  if (prevPath != 0)
	     newPath = new char[strlen(prevPath)+2+dirNameLength]; // / + '\0';
	  else
	    newPath = new char[2+dirNameLength]; // / + '\0'
	  
	  if (prevPath != 0) {
	     strcpy(newPath, prevPath);
	     strcat(newPath, dirName);
	  } else
		  strcpy(newPath, dirName);
	  strcat(newPath,"/");
	  
	  File *theNextDir = new File(dirName, newPath, true);
	  currentDir->addFile(theNextDir);
	  currentDir = theNextDir;
	
	  delete [] newPath;
	
	} else {
	  currentDir = nextDir;
	}
      }
    }
  } else
    combinedFile = combined; // no dir


  File *file = 0;

  file = new File(combinedFile, fileDescription, false);

  if (currentDir->addFile(file) != 0) // file already existed
    delete file;
  
  delete [] combined;
  return 0;
}

const char *
File::getName(void)
{
  if (name.length() == 0)
    return 0;

  return name.c_str();
}

const char *
File::getDescription(void)
{
  if (description.length() == 0)
    return 0;

  return description.c_str();
}


/*
void
File::setDescription(const char *newDescription)
{
  description.assign(newDescription);
}
*/

void
File::setParentDir(File *dir)
{
  parentDir = dir;
}

File *
File::getParentDir(void)
{
  return parentDir;
}



bool 
File::isDir(void)
{
  return isDirectory;
}


File *
File::getFile(const char *filename)
{
  if (isDirectory == false)
    return 0;

  if (strcmp(filename, "..") == 0)
    return this->getParentDir();

  if (strcmp(filename, ".") == 0)
    return this;

  map<string, File *>::iterator theMapObject = dirFiles.find(filename);
  if (theMapObject == dirFiles.end())
    return 0;
  else
    return theMapObject->second;

}

FileIter 
File::getFiles(void)
{
  return FileIter(*this);
}

int
File::getNumFiles(void)
{
  return dirFiles.size();
}

