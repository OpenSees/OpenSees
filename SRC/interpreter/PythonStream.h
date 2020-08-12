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
    PythonStream(int indentSize = 2, bool echo = true, bool standard_echo = false)
            : StandardStream(indentSize, standard_echo), error(0), msg(), echoApplication(echo) {}

    ~PythonStream() {}

    void setError(PyObject *err) {
        error = err;
    }

    int setFile(const char *fileName, openMode mode, bool echo) {
        echoApplication = echo;
        return StandardStream::setFile(fileName, mode, false);
    }

    OPS_Stream &operator<<(char c) {
        if (echoApplication) err_out(c);
        return StandardStream::operator<<(c);
    }

    OPS_Stream &operator<<(unsigned char c) {
        if (echoApplication) err_out(c);
        return StandardStream::operator<<(c);
    }

    OPS_Stream &operator<<(signed char c) {
        if (echoApplication) err_out(c);
        return StandardStream::operator<<(c);
    }

    OPS_Stream &operator<<(const char *s) {
        if (echoApplication) err_out(s);
        return StandardStream::operator<<(s);
    }

    OPS_Stream &operator<<(const unsigned char *s) {
        if (echoApplication) err_out(s);
        return StandardStream::operator<<(s);
    }

    OPS_Stream &operator<<(const signed char *s) {
        if (echoApplication) err_out(s);
        return StandardStream::operator<<(s);
    }

    OPS_Stream &operator<<(int n) {
        if (echoApplication) err_out(n);
        return StandardStream::operator<<(n);
    }

    OPS_Stream &operator<<(unsigned int n) {
        if (echoApplication) err_out(n);
        return StandardStream::operator<<(n);
    }

    OPS_Stream &operator<<(long n) {
        if (echoApplication) err_out(n);
        return StandardStream::operator<<(n);
    }

    OPS_Stream &operator<<(unsigned long n) {
        if (echoApplication) err_out(n);
        return StandardStream::operator<<(n);
    }

    OPS_Stream &operator<<(short n) {
        if (echoApplication) err_out(n);
        return StandardStream::operator<<(n);
    }

    OPS_Stream &operator<<(unsigned short n) {
        if (echoApplication) err_out(n);
        return StandardStream::operator<<(n);
    }

    OPS_Stream &operator<<(bool b) {
        if (echoApplication) err_out(b);
        return StandardStream::operator<<(b);
    }

    OPS_Stream &operator<<(double n) {
        if (echoApplication) err_out(n);
        return StandardStream::operator<<(n);
    }

    OPS_Stream &operator<<(float n) {
        if (echoApplication) err_out(n);
        return StandardStream::operator<<(n);
    }

    OPS_Stream &operator<<(const void *p) {
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

private:

    template<class T>
    void err_out(T err) {
        std::stringstream ss;
        ss << err;
        msg = ss.str();
        PySys_FormatStderr(msg.c_str());
    }

    PyObject *error;
    std::string msg;
    bool echoApplication;
};

#endif
