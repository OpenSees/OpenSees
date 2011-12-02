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
// $Date: 2006-08-03 23:34:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DummyStream.h,v $

#ifndef _DummyStream
#define _DummyStream

#include <OPS_Stream.h>

#include <fstream>
using std::ofstream;

class DummyStream : public OPS_Stream
{
 public:
  DummyStream();
  ~DummyStream();

  // xml stuff
  int tag(const char *) {return 0;};
  int tag(const char *, const char *) {return 0;};
  int endTag() {return 0;};
  int attr(const char *name, int value) {return 0;};
  int attr(const char *name, double value) {return 0;};
  int attr(const char *name, const char *value) {return 0;};
  int write(Vector &data) {return 0;};


  OPS_Stream& write(const char *s, int n) {return *this;};
  OPS_Stream& write(const unsigned char *s, int n) {return *this;};
  OPS_Stream& write(const signed char *s, int n) {return *this;};
  OPS_Stream& write(const void *s, int n) {return *this;};
  OPS_Stream& operator<<(char c) {return *this;};
  OPS_Stream& operator<<(unsigned char c) {return *this;};
  OPS_Stream& operator<<(signed char c) {return *this;};
  OPS_Stream& operator<<(const char *s) {return *this;};
  OPS_Stream& operator<<(const unsigned char *s) {return *this;};
  OPS_Stream& operator<<(const signed char *s) {return *this;};
  OPS_Stream& operator<<(const void *p) {return *this;};
  OPS_Stream& operator<<(int n) {return *this;};
  OPS_Stream& operator<<(unsigned int n) {return *this;};
  OPS_Stream& operator<<(long n) {return *this;};
  OPS_Stream& operator<<(unsigned long n) {return *this;};
  OPS_Stream& operator<<(short n) {return *this;};
  OPS_Stream& operator<<(unsigned short n) {return *this;};
  OPS_Stream& operator<<(bool b) {return *this;};
  OPS_Stream& operator<<(double n) {return *this;};
  OPS_Stream& operator<<(float n) {return *this;};

  int sendSelf(int commitTag, Channel &theChannel) {return 0;};  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker) {return 0;};

 private:
};

#endif
