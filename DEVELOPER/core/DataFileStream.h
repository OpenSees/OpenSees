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
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DataFileStream.h,v $

#ifndef _DataFileStream
#define _DataFileStream

#include <OPS_Stream.h>

#include <fstream>
using std::ofstream;

class Matrix;

class DataFileStream : public OPS_Stream
{
 public:
  DataFileStream(int indent=2);
  DataFileStream(const char *fileName, openMode mode = OVERWRITE, int indent=2, int doCSV =0, bool closeOnWrite = false, int precision = 6, bool doScientific = false);
  ~DataFileStream();

  int setFile(const char *fileName, openMode mode = OVERWRITE);
  int open(void);
  int close(void);

  int setPrecision(int precision);
  int setFloatField(floatField);
  int precision(int precision) {return 0;};
  int width(int width) {return 0;};
  const char *getFileName(void) {return fileName;}

  // xml stuff
  int tag(const char *);
  int tag(const char *, const char *);
  int endTag();
  int attr(const char *name, int value);
  int attr(const char *name, double value);
  int attr(const char *name, const char *value);
  int write(Vector &data);

  // regular stuff
  OPS_Stream& write(const char *s, int n);
  OPS_Stream& write(const unsigned char *s, int n);
  OPS_Stream& write(const signed char *s, int n);
  OPS_Stream& write(const void *s, int n);
  OPS_Stream& write(const double *s, int n);

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

  // parallel stuff
  int setOrder(const ID &orderOfData);
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);

 private:
  ofstream theFile;
  int fileOpen;
  openMode theOpenMode;
  char *fileName;

  void indent(void);
  int indentSize;
  int numIndent;
  char *indentString;

  int sendSelfCount;
  Channel **theChannels;
  int numDataRows;

  Matrix *mapping;
  int maxCount;
  ID *sizeColumns;
  ID **theColumns;
  double **theData;
  Vector **theRemoteData;

  int doCSV;
  bool closeOnWrite;

  int thePrecision;
  bool doScientific;
};

#endif
 
