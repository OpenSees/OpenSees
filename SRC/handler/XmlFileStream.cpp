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
                                                                        
// $Revision: 1.8 $
// $Date: 2007-06-25 21:18:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/XmlFileStream.cpp,v $

#include <XmlFileStream.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <ID.h>
#include <Channel.h>
#include <Message.h>

using std::cerr;
using std::ios;
using std::setiosflags;
using std::string;
using std::ifstream;
using std::getline;

XmlFileStream::XmlFileStream(int indent)
  :OPS_Stream(OPS_STREAM_TAGS_XmlFileStream), 
   fileOpen(0), fileName(0), indentSize(indent), numIndent(-1),
   attributeMode(false), numTag(0), sizeTags(0), tags(0), sendSelfCount(0)
{
  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+1];
  for (int i=0; i<indentSize; i++)
    strcpy(&indentString[i]," ");
}

XmlFileStream::XmlFileStream(const char *name, openMode mode, int indent)
  :OPS_Stream(OPS_STREAM_TAGS_XmlFileStream), 
   fileOpen(0), fileName(0), indentSize(indent), numIndent(-1),
   attributeMode(false), numTag(0), sizeTags(0), tags(0), sendSelfCount(0)
{
  if (indentSize < 1) indentSize = 1;
  indentString = new char[indentSize+1];
  for (int i=0; i<indentSize; i++)
    strcpy(&indentString[i]," ");

  this->setFile(name, mode);

  int fileNameLength = strlen(fileName);
}

XmlFileStream::~XmlFileStream()
{
    int fileNameLength = strlen(fileName);

  if (fileOpen == 1) {
    for (int i=0; i<numTag; i++) {
      this->endTag();
    }

    theFile << "</OpenSees>\n";
    theFile.close();
  }

  if (sendSelfCount != 0) {

    int fileNameLength = strlen(fileName);
 
    sprintf(&fileName[fileNameLength-2],"");
    
    theFile.open(fileName, ios::out);
    
    ifstream **theFiles = new ifstream *[sendSelfCount+1];
    string s;
    
    // open up the files
    for (int i=0; i<=sendSelfCount; i++) {
      theFiles[i] = new ifstream;
      sprintf(&fileName[fileNameLength-2],".%d",i+1);
      theFiles[i]->open(fileName, ios::in);
      if (theFiles[i]->bad()) {
	theFiles[i] = 0;
	opserr << "XmlFileStream::~XmlFileStrream - trouble opening file: " << fileName << endln;
      }      

      else {
	// read and throw away the first six lines from all but 0 (to first tag)
	for (int j=0; j<6; j++) {
		getline(*(theFiles[i]), s);	
	  if (i == 0) 
	    theFile << s << "\n";
	}
      }

    }

    // go through each file, reading a line & sending to the output file
    bool done = false;
    bool hasData = false;

    // read the xml stuff for each file up until DATA
    for (int j=0; j<=sendSelfCount; j++) {
      if (theFiles[j] != 0) {
	bool  foundData = false;
	while (foundData == false) {
		getline(*(theFiles[j]), s);	
	  if (theFiles[j]->eof()) {
	    foundData = true;
	    theFiles[j]->close();
	    delete theFiles[j];
	    theFiles[j] = 0;
	    done = true;
	  }
	      
	  string::size_type loc = s.find( "<Data>", 0 );
	  
	  if( loc != string::npos ) {
	    foundData = true;
	  } else {
	    theFile << s << "\n";
	  }
	}
      }
    }

    if (done == false)
      theFile << "  <Data>\n";

    while (done == false) {
      char c;
      for (int i=0; i<=sendSelfCount; i++) {
	if (theFiles[i] != 0) {
	  bool eoline = false;
	  getline(*(theFiles[i]), s);	
	  
	  string::size_type loc = s.find( "</Data>", 0 );
	  
	  if( loc != string::npos ) {
	    done = true;
	  } else {
	    theFile << s;
	  }
	  
	  if (theFiles[i]->eof()) {
	    done = true;
	    theFiles[i]->close();
	    delete theFiles[i];
	    theFiles[i] = 0;
	  }
	}
      }
      if (done == false) 
	theFile << "\n";
    }

    theFile << "  </Data>\n";

    done = false;
    while (done == false) {
      if (theFiles[0] == 0)
	done = true;
      else {
		  getline(*(theFiles[0]), s);	
	theFile << s << "\n";
	if (theFiles[0]->eof())
	  done = true;
      }
    }

    for (int l=0; l<=sendSelfCount; l++) {
      if (theFiles[l] != 0) {
	theFiles[l]->close();
	delete theFiles[l];
      }
    }

    delete [] theFiles;
    theFile.close();
  }

  if (tags != 0)
    delete [] tags;

  if (indentString != 0)
    delete [] indentString;

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
    fileName = new char[strlen(name)+5];
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

  if (mode == 0)
    theOpenMode = OVERWRITE;
  else
    theOpenMode = APPEND;

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

  if (sendSelfCount != 0) {
    strcat(fileName,".1");
  }


  // open file
  if (theOpenMode == OVERWRITE) 
    theFile.open(fileName, ios::out);
  else
    theFile.open(fileName, ios::out| ios::app);

  theOpenMode = APPEND;

  if (theFile.bad()) {
    std::cerr << "WARNING - XmlFileStream::open()";
    std::cerr << " - could not open file " << fileName << std::endl;
    fileOpen = 0;
    return -1;
  } else
    fileOpen = 1;

  theFile << setiosflags(ios::fixed);

  theFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  theFile << " <OpenSees\n";
  theFile << "  xmlns:xsi = \"http://www.w3.org/2001/XMLSchema-instance\"\n";
  theFile << "  xsi:noNamespaceSchemaLocation = \"http://OpenSees.berkeley.edu/xml-schema/xmlns/OpenSees.xsd\">\n";
  numIndent++;

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
  if (fileOpen == 0)
    this->open();

  if (fileOpen != 0)
    theFile << std::setprecision(prec);

  return 0;
}

int 
XmlFileStream::setFloatField(floatField field)
{
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

  if (attributeMode == true) {
    theFile << ">\n";
  }

  // output the xml for it to the file
  numIndent++;
  this->indent();
  theFile << "<" << tagName << ">" << value << "</" << tagName << ">" << endln;
  numIndent--;

  attributeMode = false;

  return 0;
}

int 
XmlFileStream::endTag()
{
  if (numTag != 0) {
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

  opserr << "XmlFileStream::endTag() - too many endTags have been called on file: " << fileName << endln;
  return -1;

}

int 
XmlFileStream::attr(const char *name, int value)
{
  if (fileOpen == 0)
    this->open();

  theFile << " " << name << "=\"" << value << "\"";
  
  return 0;
}

int 
XmlFileStream::attr(const char *name, double value)
{
  if (fileOpen == 0)
    this->open();

  theFile << " " << name << "=\"" << value << "\"";

  return 0;
}

int 
XmlFileStream::attr(const char *name, const char *value)
{
  if (fileOpen == 0)
    this->open();

  theFile << " " << name << "=\"" << value << "\"";

  return 0;
}

int 
XmlFileStream::write(Vector &data)
{
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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
  if (fileOpen == 0)
    this->open();

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

  sendSelfCount++;
 
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
