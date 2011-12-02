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
// $Date: 2004-11-24 22:40:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataOutputFileHandler.h,v $

#ifndef _DataOutputFileHandler
#define _DataOutputFileHandler

#include <DataOutputHandler.h>
#include <FileStream.h>

enum echoMode  {NONE, DATA_FILE, XML_FILE};

class DataOutputFileHandler : public DataOutputHandler
{
 public:
  DataOutputFileHandler(const char *fileName =0, 
			echoMode = NONE, 
			openMode mode = OVERWRITE);
  ~DataOutputFileHandler();

  int open(char **dataDescription, int numData);
  int write(Vector &data);

  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

 private:
  FileStream outputFile;
  char *fileName;

  echoMode theEchoMode;
  openMode theOpenMode;
  int numColumns;
};

#endif
