#include "PythonStream.h"

PythonStream::PythonStream(int indentSize, bool echo, bool standard_echo)
  : StandardStream(indentSize, standard_echo), error(0), msg(), echoApplication(echo) {}

PythonStream:: ~PythonStream() {}

void PythonStream::setError(PyObject *err) {
  error = err;
}

int PythonStream::setFile(const char *fileName, openMode mode, bool echo) {
  echoApplication = echo;
  return StandardStream::setFile(fileName, mode, false);
}

OPS_Stream &PythonStream::operator<<(char c) {
  if (echoApplication) err_out(c);
  return StandardStream::operator<<(c);
}

OPS_Stream &PythonStream::operator<<(unsigned char c) {
  if (echoApplication) err_out(c);
  return StandardStream::operator<<(c);
}

OPS_Stream &PythonStream::operator<<(signed char c) {
  if (echoApplication) err_out(c);
  return StandardStream::operator<<(c);
}

OPS_Stream &PythonStream::operator<<(const char *s) {
  if (echoApplication) err_out(s);
  return StandardStream::operator<<(s);
}

OPS_Stream &PythonStream::operator<<(const unsigned char *s) {
  if (echoApplication) err_out(s);
  return StandardStream::operator<<(s);
}

OPS_Stream &PythonStream::operator<<(const signed char *s) {
  if (echoApplication) err_out(s);
  return StandardStream::operator<<(s);
}

OPS_Stream &PythonStream::operator<<(int n) {
  if (echoApplication) err_out(n);
  return StandardStream::operator<<(n);
}

OPS_Stream &PythonStream::operator<<(unsigned int n) {
  if (echoApplication) err_out(n);
  return StandardStream::operator<<(n);
}

OPS_Stream &PythonStream::operator<<(long n) {
  if (echoApplication) err_out(n);
  return StandardStream::operator<<(n);
}

OPS_Stream &PythonStream::operator<<(unsigned long n) {
  if (echoApplication) err_out(n);
  return StandardStream::operator<<(n);
}

OPS_Stream &PythonStream::operator<<(short n) {
  if (echoApplication) err_out(n);
  return StandardStream::operator<<(n);
}

OPS_Stream &PythonStream::operator<<(unsigned short n) {
  if (echoApplication) err_out(n);
  return StandardStream::operator<<(n);
}

OPS_Stream &PythonStream::operator<<(bool b) {
  if (echoApplication) err_out(b);
  return StandardStream::operator<<(b);
}

OPS_Stream &PythonStream::operator<<(double n) {
  if (echoApplication) err_out(n);
  return StandardStream::operator<<(n);
}

OPS_Stream &PythonStream::operator<<(float n) {
  if (echoApplication) err_out(n);
  return StandardStream::operator<<(n);
}

OPS_Stream &PythonStream::operator<<(const void *p) {
  if (p != 0) {
    if (echoApplication) err_out(p);
    return StandardStream::operator<<(p);
  }
  if (echoApplication) {
    msg = "See stderr output";
  } else {
    msg = "See log file";
  }
  msg.erase(msg.find_last_not_of("\n\r ") + 1);  // Strip extra newlines from the message
  PyErr_SetString(error, msg.c_str());
  return *this;
}
