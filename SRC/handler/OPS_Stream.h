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
                                                                        
// $Revision: 1.7 $
// $Date: 2009-04-30 23:23:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/OPS_Stream.h,v $

#ifndef _OPS_Stream
#define _OPS_Stream

#include <MovableObject.h>
enum openMode  {OVERWRITE, APPEND};
enum floatField {FIXEDD, SCIENTIFIC};
class Vector;
class ID;

class OPS_Stream:  public MovableObject
{
 public:
  OPS_Stream(int classTag);
  virtual ~OPS_Stream();

  // output format
  virtual int setFile(const char *fileName, openMode mode = OVERWRITE, bool echo = false) {return 0;}
  virtual int setPrecision(int precision) {return 0;}
  virtual int setFloatField(floatField) {return 0;}
  virtual int precision(int precision) {return 0;}
  virtual int width(int width) {return 0;}

  // xml stuff
  virtual int tag(const char *) =0;
  virtual int tag(const char *, const char *) =0;
  virtual int endTag() =0;
  virtual int attr(const char *name, int value) =0;
  virtual int attr(const char *name, double value) =0;
  virtual int attr(const char *name, const char *value) =0;
  virtual int write(Vector &data) =0; 

  // regular stuff
  virtual OPS_Stream& write(const char *s, int n);
  virtual OPS_Stream& write(const unsigned char *s, int n);
  virtual OPS_Stream& write(const signed char *s, int n);
  virtual OPS_Stream& write(const void *s, int n);
  virtual OPS_Stream& write(const double *s, int n);

  virtual OPS_Stream& operator<<(char c);
  virtual OPS_Stream& operator<<(unsigned char c);
  virtual OPS_Stream& operator<<(signed char c);
  virtual OPS_Stream& operator<<(const char *s);
  virtual OPS_Stream& operator<<(const unsigned char *s);
  virtual OPS_Stream& operator<<(const signed char *s);
  virtual OPS_Stream& operator<<(const void *p);
  virtual OPS_Stream& operator<<(int n);
  virtual OPS_Stream& operator<<(unsigned int n);
  virtual OPS_Stream& operator<<(long n);
  virtual OPS_Stream& operator<<(unsigned long n);
  virtual OPS_Stream& operator<<(short n);
  virtual OPS_Stream& operator<<(unsigned short n);
  virtual OPS_Stream& operator<<(bool b);
  virtual OPS_Stream& operator<<(double n);
  virtual OPS_Stream& operator<<(float n);

  // parallel stuff
  virtual void setAddCommon(int);
  virtual int setOrder(const ID &order);
  virtual int sendSelf(int commitTag, Channel &theChannel) =0;  
  virtual int recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker) =0;

 protected:
  int addCommonFlag;

 private:
  void indent();
  int numIndent;
};

#endif
