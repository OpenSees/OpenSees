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
                                                                        
// $Revision: 1.1 $
// $Date: 2006-08-03 23:28:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/XmlFileStream.cpp,v $

#include <XmlFileStream.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <ID.h>
#include <Channel.h>
#include <Message.h>

using std::cerr;
using std::ios;
using std::setiosflags;

XmlFileStream::XmlFileStream(int indent)
  :OPS_Stream(OPS_STREAM_TAGS_XmlFileStream), 
   fileOpen(0), fileName(0), indentSize(indent), numIndent(-1),
   attributeMode(false), numTag(0), sizeTags(0), tags(0)
{
  opserr << "XmlFileStream::" << endln;
  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+1];
  //  for (int i=0; i<indentSize; i++)
  //    strcpy(indentString, " ");

  indentString = "  ";

}

XmlFileStream::XmlFileStream(const char *name, openMode mode, int indent)
  :OPS_Stream(OPS_STREAM_TAGS_XmlFileStream), 
   fileOpen(0), fileName(0), indentSize(indent), numIndent(-1),
   attributeMode(false), numTag(0), sizeTags(0), tags(0)
{
  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+1];
  //  for (int i=0; i<indentSize; i++)
  //    strcat(indentString," ");
  indentString = "  ";

  this->setFile(name, mode);

}

XmlFileStream::~XmlFileStream()
{
  for (int i=numTag; i>0; i--) {
    this->indent();
    theFile << "</" << tags[i] << ">\n";
    delete [] tags[i];
    numTag--;
  }

  delete [] tags;

  if (fileOpen == 1)
    theFile.close();
  
  if (fileName != 0)
    delete [] fileName;
}

int 
XmlFileStream::setFile(const char *name, openMode mode)
{
  if (name == 0) {
    std::cerr << "XmlFileStream::setFile() - no name passed\n";
    return -1;
  }

  // first create a copy of the file name
  if (fileName != 0) {
    if (strcmp(fileName, name) != 0)
      delete [] fileName;
    fileName = 0;
  }
  if (fileName == 0) {
    fileName = new char[strlen(name)+1];
    if (fileName == 0) {
      std::cerr << "XmlFileStream::setFile() - out of memory copying name: " << name << std::endl;
      return -1;
    }
    
    // copy the strings
    strcpy(fileName, name);
  }

  // if file already open, close it
  if (fileOpen == 1) {
    theFile.close();
    fileOpen = 0;
  }

  if (mode == OVERWRITE) 
    theFile.open(fileName, ios::out);
  else
    theFile.open(fileName, ios::out| ios::app);

  if (theFile.bad()) {
    std::cerr << "WARNING - XmlFileStream::setFile()";
    std::cerr << " - could not open file " << fileName << std::endl;

    return -1;
  } else
    fileOpen = 1;

  if (mode == 0)
    theOpenMode = OVERWRITE;
  else
    theOpenMode = APPEND;

  theFile << setiosflags(ios::fixed);


  theFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  theFile << " <OpenSees\n";
  theFile << "  xmlns:xsi = \"http://www.w3.org/2001/XMLSchema-instance\"\n";
  theFile << "  xsi:noNamespaceSchemaLocation = \"http://OpenSees.berkeley.edu//xml-schema/xmlns/OpenSees.xsd\">\n";

  return 0;
}

int 
XmlFileStream::open(void)
{
  // check setFile has been called
  if (fileName == 0) {
    std::cerr << "XmlFileStream::open(void) - no file name has been set\n";
    return -1;
  }

  // if file already open, return
  if (fileOpen == 1) {
    return 0;
  }

  // open file
  theFile.open(fileName, ios::out| ios::app);
  if (theFile.bad()) {
    std::cerr << "WARNING - XmlFileStream::open()";
    std::cerr << " - could not open file " << fileName << std::endl;

    return -1;
  } else
    fileOpen = 1;



  return 0;
}

int 
XmlFileStream::close(void)
{
  if (fileOpen != 0)
    theFile.close();
  fileOpen = 0;
  return 0;
}


int 
XmlFileStream::setPrecision(int prec)
{
  if (fileOpen != 0)
    theFile << std::setprecision(prec);

  return 0;
}

int 
XmlFileStream::setFloatField(floatField field)
{
  if (field == FIXEDD) {
    if (fileOpen != 0)
      theFile << setiosflags(ios::fixed);
  }
  else if (field == SCIENTIFIC) {
    if (fileOpen != 0)
      theFile << setiosflags(ios::scientific);
  }

  return 0;
}


int 
XmlFileStream::tag(const char *tagName)
{
  //
  // copy tagName to end of list of tags
  //

  // if tags array not large enough, expand
  if (numTag == sizeTags) {
    int nextSize = 2*sizeTags;
    if (nextSize == 0) nextSize = 32;
    char **nextTags = new char *[nextSize];
    if (nextTags != 0) {
      for (int i=0; i<sizeTags; i++)
	nextTags[i] = tags[i];
      for (int j=sizeTags+1; j<nextSize; j++)
	nextTags[j] = 0;
      sizeTags = nextSize;
    } else {
      sizeTags = 0;
      delete [] tags;
      tags = 0;
      return -1;
    }
    
    if (tags != 0)
      delete [] tags;

    tags = nextTags;
  } 

  // copy string and assign to end of array
  char *newTag = new char[strlen(tagName) +1];
  strcpy(newTag, tagName);
  tags[numTag++] = newTag;


  if (attributeMode == true) {
    theFile << ">";
  }

  // output the xml for it to the file

  numIndent++;
  theFile << endln;
  this->indent();
  theFile << "<" << tagName;
  
  attributeMode = true;

  return 0;
}


int 
XmlFileStream::tag(const char *tagName, const char *value)
{
  if (attributeMode == true) {
    theFile << ">\n";
  }

  // output the xml for it to the file
  numIndent++;
  this->indent();
  theFile << "<" << tagName << ">" << value << "<" << tagName << "/>" << endln;
  numIndent--;

  attributeMode = false;

  return 0;
}

int 
XmlFileStream::endTag()
{
  if (attributeMode == true) {
    theFile << "/>\n";
    delete [] tags[numTag-1];
    numTag--;
  } else {
    this->indent();
    theFile << "</" << tags[numTag-1] << ">\n";
    delete [] tags[numTag-1];
    numTag--;
  }    

  attributeMode = false;
  numIndent--;

  return 0;
}

int 
XmlFileStream::attr(const char *name, int value)
{
  theFile << " " << name << "=\"" << value << "\"";
  
  return 0;
}

int 
XmlFileStream::attr(const char *name, double value)
{
  theFile << " " << name << "=\"" << value << "\"";

  return 0;
}

int 
XmlFileStream::attr(const char *name, const char *value)
{
  theFile << " " << name << "=\"" << value << "\"";

  return 0;
}

int 
XmlFileStream::write(Vector &data)
{
  if (attributeMode == true) {
    theFile << ">\n";
    attributeMode = false;
    numIndent++;
  }

  this->indent();
  (*this) << data;  

  return 0;
}



OPS_Stream& 
XmlFileStream::write(const char *s,int n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile.write(s, n);

  return *this;
}

OPS_Stream& 
XmlFileStream::write(const unsigned char*s,int n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
XmlFileStream::write(const signed char*s,int n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
XmlFileStream::write(const void *s, int n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(char c)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(unsigned char c)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(signed char c)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(const char *s)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  // note that we do the flush so that a "/n" before
  // a crash will cause a flush() - similar to what 
  if (fileOpen != 0) {
    theFile << s;
    theFile.flush();
  }

  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(const unsigned char *s)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(const signed char *s)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(const void *p)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

/*
  if (fileOpen != 0)
    theFile << p;
*/
  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(int n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(unsigned int n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(long n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(unsigned long n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(short n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(unsigned short n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

/*
  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(bool b)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

/*
  if (fileOpen != 0)
    theFile << b;
*/
  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(double n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
XmlFileStream::operator<<(float n)
{
  if (attributeMode == true) {
    theFile << "/>\n";
    attributeMode = false;
  }

  if (fileOpen != 0)
    theFile << n;

  return *this;
}


int 
XmlFileStream::sendSelf(int commitTag, Channel &theChannel)
{
  static ID idData(2);
  int fileNameLength = 0;
  if (fileName != 0)
    fileNameLength = strlen(fileName);

  idData(0) = fileNameLength;

  if (theOpenMode == OVERWRITE)
    idData(1) = 0;
  else
    idData(1) = 1;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "XmlFileStream::sendSelf() - failed to send id data\n";
    return -1;
  }

  if (fileNameLength != 0) {
    Message theMessage(fileName, fileNameLength);
    if (theChannel.sendMsg(0, commitTag, theMessage) < 0) {
      opserr << "XmlFileStream::sendSelf() - failed to send message\n";
      return -1;
    }
  }

  return 0;
}

int 
XmlFileStream::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID idData(2);

  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "XmlFileStream::recvSelf() - failed to recv id data\n";
    return -1;
  }

  int fileNameLength = idData(0);
  if (idData(1) == 0)
    theOpenMode = OVERWRITE;
  else
    theOpenMode = APPEND;

  if (fileNameLength != 0) {
    if (fileName != 0)
      delete [] fileName;
    fileName = new char[fileNameLength+5];
    if (fileName == 0) {
      opserr << "XmlFileStream::recvSelf() - out of memory\n";
      return -1;
    }

    Message theMessage(fileName, fileNameLength);
    if (theChannel.recvMsg(0, commitTag, theMessage) < 0) {
      opserr << "XmlFileStream::recvSelf() - failed to recv message\n";
      return -1;
    }
    sprintf(&fileName[fileNameLength],".%d",commitTag);

    if (this->setFile(fileName, theOpenMode) < 0) {
      opserr << "XmlFileStream::XmlFileStream() - setFile() failed\n";
      if (fileName != 0) {
	delete [] fileName;
	fileName = 0;
      }
    }
  }
  
  return 0;
}

void
XmlFileStream::indent(void)
{
  if (fileOpen != 0)
    for (int i=0; i<numIndent; i++)
      theFile << indentString;
}
