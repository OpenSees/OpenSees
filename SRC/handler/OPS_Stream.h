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
// $Date: 2008-01-16 00:19:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/OPS_Stream.h,v $

#ifndef _OPS_Stream
#define _OPS_Stream

#include <MovableObject.h>
enum openMode  {OVERWRITE, APPEND};
enum floatField {FIXEDD, SCIENTIFIC};
class Vector;


class OPS_Stream:  public MovableObject
{
 public:
  OPS_Stream(int classTag);
  virtual ~OPS_Stream();

  virtual int setFile(const char *fileName, openMode mode = OVERWRITE) {return 0;};
  virtual int setPrecision(int precision) {return 0;};
  virtual int setFloatField(floatField) {return 0;};
  virtual int precision(int precision) {return 0;};
  virtual int width(int width) {return 0;};

  // xml stuff
  virtual int tag(const char *) =0;
  virtual int tag(const char *, const char *) =0;
  virtual int endTag() =0;
  virtual int attr(const char *name, int value) =0;
  virtual int attr(const char *name, double value) =0;
  virtual int attr(const char *name, const char *value) =0;
  virtual int write(Vector &data) =0; 

  // regular stuff
  virtual OPS_Stream& write(const char *s, int n) {return *this;};
  virtual OPS_Stream& write(const unsigned char *s, int n) {return *this;};
  virtual OPS_Stream& write(const signed char *s, int n) {return *this;};
  virtual OPS_Stream& write(const void *s, int n) {return *this;};

  virtual OPS_Stream& operator<<(char c) {return *this;};
  virtual OPS_Stream& operator<<(unsigned char c) {return *this;};
  virtual OPS_Stream& operator<<(signed char c) {return *this;};
  virtual OPS_Stream& operator<<(const char *s) {return *this;};
  virtual OPS_Stream& operator<<(const unsigned char *s) {return *this;};
  virtual OPS_Stream& operator<<(const signed char *s) {return *this;};
  virtual OPS_Stream& operator<<(const void *p) {return *this;};
  virtual OPS_Stream& operator<<(int n) {return *this;};
  virtual OPS_Stream& operator<<(unsigned int n) {return *this;};
  virtual OPS_Stream& operator<<(long n) {return *this;};
  virtual OPS_Stream& operator<<(unsigned long n) {return *this;};
  virtual OPS_Stream& operator<<(short n) {return *this;};
  virtual OPS_Stream& operator<<(unsigned short n) {return *this;};
  virtual OPS_Stream& operator<<(bool b) {return *this;};
  virtual OPS_Stream& operator<<(double n) {return *this;};
  virtual OPS_Stream& operator<<(float n) {return *this;};

  int sendSelf(int commitTag, Channel &theChannel) =0;  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker) =0;


 private:
  void indent(void);
  int numIndent;


};

#endif
