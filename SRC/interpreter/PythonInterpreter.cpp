/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS 
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */

#include "PythonInterpreter.h"
#include "PythonStream.h"

// define opserr
static PythonStream sserr;
OPS_Stream *opserrPtr = &sserr;

PythonInterpreter::PythonInterpreter(int argc, char **argv)
    :wrapper(), cmds(this), retcode(0)
{

    /* fmk - beginning of modifications for OpenSees */
    fprintf(stderr,"\n\n\t OpenSeesPy -- Open System For Earthquake Engineering Simulation");
    fprintf(stderr,"\n\tPacific Earthquake Engineering Research Center -- ");
    fprintf(stderr, OPS_VERSION);
    fprintf(stderr, "\n\n");

    fprintf(stderr,"\t    (c) Copyright 2015-2017 The Regents of the University of California");
    fprintf(stderr,"\n\t\t\t\t All Rights Reserved\n\n\n");

    retcode = Py_Main(argc, reinterpret_cast<wchar_t **>(argv));
}

PythonInterpreter::~PythonInterpreter() {
  // does nothing
}

int
PythonInterpreter::run(){

  return retcode;
}
			 

int 
PythonInterpreter::addCommand(const char *, Command &) {
  return -1;
}

int 
PythonInterpreter::removeCommand(const char *) {
  return -1;
}

int 
PythonInterpreter::getNumRemainingInputArgs(void) {
    return wrapper.getNumberArgs() - wrapper.getCurrentArg();
}

int 
PythonInterpreter::getInt(int *data, int numArgs) {
    if ((wrapper.getNumberArgs() - wrapper.getCurrentArg()) < numArgs) {
        return -1;
    }

    for (int i = 0; i < numArgs; i++) {
        PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
        wrapper.incrCurrentArg();
        if (PyLong_Check(o) || PyFloat_Check(o) || PyBool_Check(o)) {
            PyErr_Clear();
            data[i] = PyLong_AsLong(o);
            if (PyErr_Occurred()) {
                return -1;
            }
        } else {
            return -1;
        }
    }

    return 0;
}

int 
PythonInterpreter::getDouble(double *data, int numArgs) {
    if ((wrapper.getNumberArgs() - wrapper.getCurrentArg()) < numArgs) {
        return -1;
    }

    for (int i = 0; i < numArgs; i++) {
        PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
        wrapper.incrCurrentArg();
        if (PyLong_Check(o) || PyFloat_Check(o) || PyBool_Check(o)) {
            PyErr_Clear();
            data[i] = PyFloat_AsDouble(o);
            if (PyErr_Occurred()) {
                return -1;
            }
        } else {
            return -1;
        }
    }
    
    return 0;
}

const char* 
PythonInterpreter::getString() {

    if (wrapper.getCurrentArg() >= wrapper.getNumberArgs()) {
        return 0;
    }

    PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
    wrapper.incrCurrentArg();
#if PY_MAJOR_VERSION >= 3
    if (!PyUnicode_Check(o)) {
        return 0;
    }

    PyObject* space = PyUnicode_FromString(" ");
    PyObject* empty = PyUnicode_FromString("");
    PyObject* newo = PyUnicode_Replace(o, space, empty, -1);
    const char* res = PyUnicode_AsUTF8(newo);

    Py_DECREF(newo);
    Py_DECREF(space);
    Py_DECREF(empty);

    return res;
#else
    if (!PyString_Check(o)) {
        return 0;
    }

    return PyString_AS_STRING(o);
#endif

}

int 
PythonInterpreter::getStringCopy(char **stringPtr) {
  return -1;
}

void
PythonInterpreter::resetInput(int cArg)
{
    wrapper.resetCommandLine(cArg);
}

int
PythonInterpreter::setInt(int* data, int numArgs, bool scalar)
{
    wrapper.setOutputs(data, numArgs, scalar);
    return 0;
}

int
PythonInterpreter::setDouble(double* data, int numArgs, bool scalar)
{
    wrapper.setOutputs(data, numArgs, scalar);
    return 0;
}

int
PythonInterpreter::setString(const char* str)
{
    wrapper.setOutputs(str);
    return 0;
}

int
PythonInterpreter::runCommand(const char *cmd) {
    return PyRun_SimpleString(cmd);
}