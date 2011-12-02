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
// $Date: 2003-02-25 23:33:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/FileStream.h,v $

#ifndef _FileStream
#define _FileStream

#include <OPS_Stream.h>

#include <fstream>
using std::ofstream;

class FileStream : public OPS_Stream
{
 public:
  FileStream();
  virtual ~FileStream();

  int setFile(const char *fileName, openMode mode = OVERWRITE);
  int setPrecision(int precision);
  int setFloatField(floatField);
  int precision(int precision) {return 0;};
  int width(int width) {return 0;};
  int close(void);

  OPS_Stream& write(const char *s, int n);
  OPS_Stream& write(const unsigned char *s, int n);
  OPS_Stream& write(const signed char *s, int n);
  OPS_Stream& write(const void *s, int n);
  OPS_Stream& operator<<(char c);
  OPS_Stream& operator<<(unsigned char c);
  OPS_Stream& operator<<(signed char c);
  OPS_Stream& operator<<(const char *s);
  OPS_Stream& operator<<(const unsigned char *s);
  OPS_Stream& operator<<(const signed char *s);
  OPS_Stream& operator<<(const void *p);
  OPS_Stream& operator<<(int n);
  OPS_Stream& operator<<(unsigned int n);
  OPS_Stream& operator<<(long n);
  OPS_Stream& operator<<(unsigned long n);
  OPS_Stream& operator<<(short n);
  OPS_Stream& operator<<(unsigned short n);
  OPS_Stream& operator<<(bool b);
  OPS_Stream& operator<<(double n);
  OPS_Stream& operator<<(float n);
  // OPS_Stream& operator<<(__omanip func);
  // OPS_Stream& operator<<(__manip func);
  // OPS_Stream& operator<<(streambuf*);
  // OPS_Stream& ends(OPS_Stream& outs);
  // OPS_Stream& flush(OPS_Stream& outs);
  // OPS_Stream& (OPS_Stream& outs);

 private:
  ofstream theFile;
  int fileOpen;
};

#endif
