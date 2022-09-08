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

// $Revision$
// $Date$

#ifndef _PythonStream
#define _PythonStream

#include <StandardStream.h>
#include <Python.h>
#include <sstream>

class PythonStream : public StandardStream {
public:
  PythonStream(int indentSize = 2, bool echo = true, bool standard_echo = false);


  ~PythonStream();

  void setError(PyObject *err);
  int setFile(const char *fileName, openMode mode, bool echo);
  OPS_Stream &operator<<(char c);
  OPS_Stream &operator<<(unsigned char c);
  OPS_Stream &operator<<(signed char c);
  OPS_Stream &operator<<(const char *s);
  OPS_Stream &operator<<(const unsigned char *s);
  OPS_Stream &operator<<(const signed char *s);
  OPS_Stream &operator<<(int n);
  OPS_Stream &operator<<(unsigned int n);
  OPS_Stream &operator<<(long n);
  OPS_Stream &operator<<(unsigned long n);
  OPS_Stream &operator<<(short n);
  OPS_Stream &operator<<(unsigned short n);
  OPS_Stream &operator<<(bool b);
  OPS_Stream &operator<<(double n);
  OPS_Stream &operator<<(float n);
  OPS_Stream &operator<<(const void *p);

private:
  template<class T>void err_out(T err);
  PyObject *error;
  std::string msg;
  bool echoApplication;
};

template<class T>void PythonStream::err_out(T err)
{
  std::stringstream ss;
  ss << err;
  msg = ss.str();
  PySys_FormatStderr(msg.c_str());
}


#endif
