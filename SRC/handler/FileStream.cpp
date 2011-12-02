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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-03-04 23:12:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/FileStream.cpp,v $


#include <FileStream.h>
#include <iostream>
#include <iomanip>
using std::cerr;
using std::ios;
using std::setiosflags;

FileStream::FileStream()
  :fileOpen(0)
{

}
FileStream::~FileStream()
{
  if (fileOpen == 1)
    theFile.close();
}

int 
FileStream::setFile(const char *fileName, openMode mode)
{
  if (fileOpen == 1) {
    theFile.close();
    fileOpen = 0;
  }

  if (mode == OVERWRITE) 
    theFile.open(fileName, ios::out);
  else
    theFile.open(fileName, ios::out| ios::app);

  if (theFile.bad()) {
    std::cerr << "WARNING - FileStream::setFile()";
    std::cerr << " - could not open file " << fileName << std::endl;

    return -1;
  } else
    fileOpen = 1;

  return 0;
}


int 
FileStream::setPrecision(int prec)
{
  if (fileOpen != 0)
    theFile << std::setprecision(prec);

  return 0;
}


int 
FileStream::close(void)
{
  if (fileOpen != 0)
    theFile.close();
  fileOpen = 0;
  return 0;
}


int 
FileStream::setFloatField(floatField field)
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



OPS_Stream& 
FileStream::write(const char *s,int n)
{
  if (fileOpen != 0)
    theFile.write(s, n);

  return *this;
}

OPS_Stream& 
FileStream::write(const unsigned char*s,int n)
{
  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
FileStream::write(const signed char*s,int n)
{
  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
FileStream::write(const void *s, int n)
{
  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
FileStream::operator<<(char c)
{
  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(unsigned char c)
{
  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(signed char c)
{
  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(const char *s)
{
  // note that we do the flush so that a "/n" before
  // a crash will cause a flush() - similar to what 
  if (fileOpen != 0) {
    theFile << s;
    theFile.flush();
  }

  return *this;
}
OPS_Stream& 
FileStream::operator<<(const unsigned char *s)
{
  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(const signed char *s)
{
  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(const void *p)
{
  if (fileOpen != 0)
    theFile << p;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(int n)
{
  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(unsigned int n)
{
  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(long n)
{
  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(unsigned long n)
{
  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(short n)
{
  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(unsigned short n)
{
  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(bool b)
{
  if (fileOpen != 0)
    theFile << b;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(double n)
{
  if (fileOpen != 0)
    theFile << n;

  return *this;
}
OPS_Stream& 
FileStream::operator<<(float n)
{
  if (fileOpen != 0)
    theFile << n;

  return *this;
}
