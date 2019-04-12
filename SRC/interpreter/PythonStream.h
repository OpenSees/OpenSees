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

class PythonStream : public StandardStream
{
public:
    PythonStream(int indentSize=2, bool echo=false)
	:StandardStream(indentSize,echo), error(0), message(0), msg() {}
    ~PythonStream() {}
    
    void setError(PyObject* err, PyObject* msg) {
	error = err;
	message = msg;
    }

    OPS_Stream& operator<<(char c) {
	err_out(c);
	return this->StandardStream::operator<<(c);
    }
    
    OPS_Stream& operator<<(unsigned char c){
	err_out(c);
	return this->StandardStream::operator<<(c);
    }
    
    OPS_Stream& operator<<(signed char c){
	err_out(c);
	return this->StandardStream::operator<<(c);
    }
    
    OPS_Stream& operator<<(const char *s){
	err_out(s);
	return this->StandardStream::operator<<(s);
    }
    
    OPS_Stream& operator<<(const unsigned char *s){
	err_out(s);
	return this->StandardStream::operator<<(s);
    }
    
    OPS_Stream& operator<<(const signed char *s){
	err_out(s);
	return this->StandardStream::operator<<(s);
    }
    
    OPS_Stream& operator<<(int n){
	err_out(n);
	return this->StandardStream::operator<<(n);
    }
    
    OPS_Stream& operator<<(unsigned int n){
	err_out(n);
	return this->StandardStream::operator<<(n);
    }
    
    OPS_Stream& operator<<(long n){
	err_out(n);
	return this->StandardStream::operator<<(n);
    }
    
    OPS_Stream& operator<<(unsigned long n){
	err_out(n);
	return this->StandardStream::operator<<(n);
    }
    
    OPS_Stream& operator<<(short n){
	err_out(n);
	return this->StandardStream::operator<<(n);
    }
    
    OPS_Stream& operator<<(unsigned short n){
	err_out(n);
	return this->StandardStream::operator<<(n);
    }
    
    OPS_Stream& operator<<(bool b){
	err_out(b);
	return this->StandardStream::operator<<(b);
    }
    
    OPS_Stream& operator<<(double n){
	err_out(n);
	return this->StandardStream::operator<<(n);
    }
    
    OPS_Stream& operator<<(float n){
	err_out(n);
	return this->StandardStream::operator<<(n);
    }

    OPS_Stream& operator<<(const void *p) {
	if (p!=0) {
	    return this->StandardStream::operator<<(p);
	}
	if (msg.empty()) {
	    msg = "See opensees.msg\n";
	}
	PyErr_SetString(error, msg.c_str());
	return *this;
    }

private:

    template<class T>
    void err_out(T err) {

	std::stringstream ss;
	ss << err;
	msg += ss.str();

	std::size_t pos = msg.find('\n');
	while (pos != std::string::npos) {
	    std::string sub = msg.substr(0, pos+1);
	    if (sub.empty()==false && sub!="\n") {
		PyErr_SetString(message, sub.c_str());
		PyErr_Print();
	    }

	    msg = msg.substr(pos+1);
	    pos = msg.find('\n');
	}
    }

    PyObject* error;
    PyObject* message;
    std::string msg;
};

#endif
