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
                                                                        
// $Revision: 1.5 $
// $Date: 2006-08-07 22:18:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/SimulationInformation.cpp,v $
//
// Description: This file contains the class definition for SimulationInformation.
// SimulationInformation is a stopwatch.
//
// What: "@(#) SimulationInformation.h, revA"

#include<SimulationInformation.h>
#include <OPS_Globals.h>
#include <string.h>
#include <time.h>

SimulationInformation::SimulationInformation() 
  :filesRead(0), filesWritten(0), paramNames(0), paramValues(0),
   numFilesWritten(0), numFilesRead(0), numParameters(0)
{

  strcpy(startTime," ");
  strcpy(endTime," ");
  
}


int
SimulationInformation::start(void)
{
  //
  // reset all read and written files, and reset all paramaters passed
  //

  if (filesRead != 0) {
    for (int i=0; i<numFilesRead; i++)
      delete [] filesRead[i];
    delete [] filesRead;
  }
  filesRead = 0;
  numFilesRead = 0;
  
  if (filesWritten != 0) {
    for (int i=0; i<numFilesWritten; i++)
      delete [] filesWritten[i];
    delete [] filesWritten;
  }
  filesWritten = 0;
  numFilesWritten = 0;
  
  if (paramNames != 0) {
    for (int i=0; i<numParameters; i++) {
      delete [] paramNames[i];
      delete [] paramValues[i];
    }
    delete [] paramNames;
    delete [] paramValues;
  }
  numParameters = 0;
  paramNames = 0;
  paramValues = 0;

  // now set the start time
  
  time_t timeT;
  if (time(&timeT) != 0) {
#ifdef _WIN32
    const char *sTime = ctime(&timeT);
	strcpy(startTime, sTime);
#else
	ctime_r(&timeT, &startTime[0]);
#endif
	
  }
  
  return 0;
}

int
SimulationInformation::end(void)
{
  time_t timeT;
 if (time(&timeT) != 0) {
#ifdef _WIN32
	const char *eTime = ctime(&timeT);
	strcpy(endTime, eTime);
#else
	ctime_r(&timeT, &endTime[0]);
#endif
  }
  
  return 0;
}

 SimulationInformation::~SimulationInformation()
 { 
	 
   if (filesRead != 0) {
     for (int i=0; i<numFilesRead; i++) {
       delete [] filesRead[i];
     }
     delete [] filesRead;
   }

   if (filesWritten != 0) {
     for (int i=0; i<numFilesWritten; i++)
       delete [] filesWritten[i];
     delete [] filesWritten;
   }

   if (paramNames != 0) {
     for (int i=0; i<numParameters; i++) {
       delete [] paramNames[i];
       delete [] paramValues[i];
     }
     delete [] paramNames;
     delete [] paramValues;
   }
   
 }


 int 
 SimulationInformation::addReadFile(const char *fileName)
 {
   // check valid fileName
   if (fileName == 0) 
     return -1;

   // don't add string ending "history.tcl"
   int filenameLength = strlen(fileName+1);
   int historyLocation = filenameLength-10;
   if (historyLocation > 0) {
     if (strcmp("history.tcl",&fileName[historyLocation]) == 0)
       return 0;
   }

   // create new array to hold pointers and copy pointers there
   char **nextFiles = new char *[numFilesRead+1];
   if (nextFiles == 0) 
     return -2;

   for (int i=0; i<numFilesRead; i++)
     nextFiles[i] = filesRead[i];

   // create new pointer for new file name and add to end of array
   char *copyFileName = new char[strlen(fileName)+1];
   if (copyFileName == 0)
     return -3;

  strcpy(copyFileName, fileName);
  nextFiles[numFilesRead] = copyFileName;
  
  // delete old array and reset the array pointer
  if (filesRead != 0)
    delete [] filesRead;
  filesRead = nextFiles;
  
  // increment number of files
  numFilesRead++;

  return 0;
}

int 
SimulationInformation::addWriteFile(const char *fileName)
{
	
  // if valid fileName
  if (fileName == 0) 
    return -1;

  // create new array to hold pointers and copy pointers there
  char **nextFiles = new char *[numFilesWritten+1];
  if (nextFiles == 0) 
    return -2;

  for (int i=0; i<numFilesWritten; i++)
    nextFiles[i] = filesWritten[i];
  
  // create new pointer for new file name and add to end of array
  char *copyFileName = new char[strlen(fileName)+1];
  if (copyFileName == 0)
    return -3;

  strcpy(copyFileName, fileName);

  nextFiles[numFilesWritten] = copyFileName;


  // delete old array and reset the array pointer
  if (filesWritten != 0)
    delete [] filesWritten;
  filesWritten = nextFiles;
  
  // increment number of files
  numFilesWritten++;

  // return ok
  
  return 0;
}

int 
SimulationInformation::addParameter(const char *name, const char *value)
{
  // check for valid name & value
  if (name == 0 || value == 0) 
    return -1;

  // create new array to hold pointers and copy pointers there
  char **nextNames = new char *[numParameters+1];
  char **nextValues = new char *[numParameters+1];
  if (nextNames == 0 && nextValues == 0)
    return -2;

  for (int i=0; i<numParameters; i++) {
    nextNames[i] = paramNames[i];
    nextValues[i] = paramValues[i];
  }  

  // create new pointer for new file name and add to end of array
  char *copyParamName = new char[strlen(name)+1];
  char *copyParamValue = new char[strlen(value)+1];
  if (copyParamName == 0 || copyParamValue == 0) {
    delete [] nextNames;
    delete [] nextValues;
    return -3;
  }
  
  strcpy(copyParamName, name);
  strcpy(copyParamValue, value);
  nextNames[numParameters] = copyParamName;
  nextValues[numParameters] = copyParamValue;

      
  // delete old array and reset the array pointer
  if (paramNames != 0)
    delete [] paramNames;
  paramNames = nextNames;
  if (paramValues != 0)
    delete [] paramValues;
  paramValues = nextValues;
  
  // increment number of files
  numParameters++;
  
  return 0;
}

void 
SimulationInformation::Print(OPS_Stream &s) const
{
  char version[10];
  strcpy(version,OPS_VERSION);

  s << "Program: OpenSees\n";
  s << "Version: " << version << endln;
  s << "Start Time: " << startTime;
  s << "End Time: " << endTime;
  s << "Input Files:\n";
  for (int i=0; i<numFilesRead; i++)
    s << "  " << filesRead[i] << "\n";
  s << "\nOutput Files:\n";
  for (int k=0; k<numFilesWritten; k++)
    s << "  " << filesWritten[k] << "\n";
  s << "\nParameters:\n";
  for (int j=0; j<numParameters; j++)
    s << "  " << paramNames[j] << " " << paramValues[j] << "\n";
  s << endln;
}    

OPS_Stream &operator<<(OPS_Stream &s, const SimulationInformation &E)
{
  E.Print(s);
  return s;
}


