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
// $Source: /usr/local/cvs/OpenSees/SRC/handler/DatabaseStream.h,v $

#ifndef _DatabaseStream
#define _DatabaseStream

#include <OPS_Stream.h>
class FE_Datastore;

class DatabaseStream : public OPS_Stream
{
 public:
  DatabaseStream(FE_Datastore *database =0, const char *tableName =0);
  ~DatabaseStream();
  
  // xml stuff
  int tag(const char *);
  int tag(const char *, const char *);
  int endTag();
  int attr(const char *name, int value);
  int attr(const char *name, double value);
  int attr(const char *name, const char *value);
  int write(Vector &data);

  // regular stuff
  OPS_Stream& operator<<(char c) ;
  OPS_Stream& operator<<(unsigned char c) ;
  OPS_Stream& operator<<(signed char c) ;
  OPS_Stream& operator<<(const char *s) ;
  OPS_Stream& operator<<(const unsigned char *s) ;
  OPS_Stream& operator<<(const signed char *s) ;
  OPS_Stream& operator<<(const void *p) ;
  OPS_Stream& operator<<(int n) ;
  OPS_Stream& operator<<(unsigned int n) ;
  OPS_Stream& operator<<(long n) ;
  OPS_Stream& operator<<(unsigned long n) ;
  OPS_Stream& operator<<(short n) ;
  OPS_Stream& operator<<(unsigned short n) ;
  OPS_Stream& operator<<(bool b) ;
  OPS_Stream& operator<<(double n) ;
  OPS_Stream& operator<<(float n) ; 
  
  int setDatabase(FE_Datastore &theDatabase, const char *tableName);

  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

 private:
  int open(void);

  FE_Datastore *theDatabase;
  char *tableName;

  int numColumns;
  char **columns;
  int commitTag;

  bool descriptionDone;
};

#endif
