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
// $Date: 2004-11-24 22:37:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/StandardStream.cpp,v $

#include <StandardStream.h>
#include <iostream>
#include <iomanip>
using std::cerr;
using std::ios;
using std::setiosflags;

StandardStream::StandardStream()
  :fileOpen(0)
{

}
StandardStream::~StandardStream()
{
  if (fileOpen == 1)
    theFile.close();
}

int 
StandardStream::setFile(const char *fileName, openMode mode)
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
    std::cerr << "WARNING - StandardStream::setFile()";
    std::cerr << " - could not open file " << fileName << std::endl;

    return -1;
  } else
    fileOpen = 1;

  return 0;
}


int 
StandardStream::setPrecision(int prec)
{
  cerr << std::setprecision(prec);

  if (fileOpen != 0)
    theFile << std::setprecision(prec);

  return 0;
}

int 
StandardStream::setFloatField(floatField field)
{
  if (field == FIXEDD) {
	  cerr << setiosflags(ios::fixed);
    if (fileOpen != 0)
      theFile << setiosflags(ios::fixed);
  }
  else if (field == SCIENTIFIC) {
    cerr << setiosflags(ios::scientific);
    if (fileOpen != 0)
      theFile << setiosflags(ios::scientific);
  }

  return 0;
}



OPS_Stream& 
StandardStream::write(const char *s,int n)
{
  cerr.write(s, n);

  if (fileOpen != 0)
    theFile.write(s, n);
  
  return *this;
}

OPS_Stream& 
StandardStream::write(const unsigned char*s, int n)
{
  cerr.write((const char *) s, n);

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
StandardStream::write(const signed char*s, int n)
{
  cerr.write((const char *)s, n);

  if (fileOpen != 0)
    theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
StandardStream::write(const void *s, int n)
{
  cerr.write((const char *)s, n);

  if (fileOpen != 0)
   theFile.write((const char *) s, n);

  return *this;
}
OPS_Stream& 
StandardStream::operator<<(char c)
{
  cerr << c;

  if (fileOpen != 0)
    theFile << c;

 return *this;
}
OPS_Stream& 
StandardStream::operator<<(unsigned char c)
{
  cerr << c;

  if (fileOpen != 0)
    theFile << c;

 return *this;
}
OPS_Stream& 
StandardStream::operator<<(signed char c)
{

  cerr << c;

  if (fileOpen != 0)
    theFile << c;

  return *this;
}
OPS_Stream& 
StandardStream::operator<<(const char *s)
{
  // note that we do the flush so that a "/n" before
  // a crash will cause a flush() - similar to what 
  cerr << s;
  cerr.flush();

  if (fileOpen != 0) {
    theFile << s;
    theFile.flush();
  }

  return *this;
}

OPS_Stream& 
StandardStream::operator<<(const unsigned char *s)
{
  cerr << s;

  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
StandardStream::operator<<(const signed char *s)
{
  cerr << s;

  if (fileOpen != 0)
    theFile << s;

  return *this;
}
OPS_Stream& 
StandardStream::operator<<(const void *p)
{
/*
//  cerr << p;

  if (fileOpen != 0)
    theFile << p;
*/
  return *this;
}
OPS_Stream& 
StandardStream::operator<<(int n)
{
  cerr <<  1.0*n;

  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}

OPS_Stream& 
StandardStream::operator<<(unsigned int n)
{
  cerr << 1.0*n;

  if (fileOpen != 0)
    theFile << 1.0*n;

  return *this;
}
OPS_Stream& 
StandardStream::operator<<(long n)
{
/*
cerr << n;

if (fileOpen != 0)
  theFile << n;
*/
  return *this;
}
OPS_Stream& 
StandardStream::operator<<(unsigned long n)
{
/*
  cerr << n;

  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
StandardStream::operator<<(short n)
{
/*
  cerr << n;

  if (fileOpen != 0)
    theFile << n;
*/
  return *this;
}
OPS_Stream& 
StandardStream::operator<<(unsigned short n)
{
/*
  cerr << n;

  if (fileOpen != 0)
    theFile << n;
*/
return *this;
}

OPS_Stream& 
StandardStream::operator<<(bool b)
{
/*
  cerr << b;

  if (fileOpen != 0)
    theFile << b;
*/
 return *this;
}
OPS_Stream& 
StandardStream::operator<<(double n)
{
  cerr << n;

  if (fileOpen != 0)
    theFile << n;

 return *this;
}
OPS_Stream& 
StandardStream::operator<<(float n)
{
  cerr << n;

  if (fileOpen != 0)
    theFile << n;

 return *this;
}
