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
                                                                        
// Written: Minjie

// Description: A python wrapper for OpenSees commands
//
#ifndef PythonWrapper_h
#define PythonWrapper_h

#ifdef _DEBUG
#undef _DEBUG
#include <Python.h>
#define _DEBUG
#else
#include <Python.h>
#endif

#include <vector>

class PythonWrapper
{
public:

    PythonWrapper();
    ~PythonWrapper();

    // reset command line
    void resetCommandLine(int nArgs, int cArg, PyObject* argv);
    void resetCommandLine(int cArg);

    // wrapper commands
    void addOpenSeesCommands();
    void addCommand(const char* name, PyCFunction proc);
    PyMethodDef* getMethods();

    // get command line arguments
    PyObject* getCurrentArgv() {return currentArgv;}
    int getCurrentArg() const {return currentArg;}
    int getNumberArgs() const {return numberArgs;}
    void incrCurrentArg() {currentArg++;}

    // set outputs
    void setOutputs(int* data, int numArgs, bool scalar);
    void setOutputs(double* data, int numArgs, bool scalar);
    void setOutputs(const char* str);
    PyObject* getResults();

private:
    // command line arguments
    PyObject* currentArgv;
    int currentArg;
    int numberArgs;

    // methods table
    std::vector<PyMethodDef> methodsOpenSees;
    const char* opensees_docstring;
    PyObject* currentResult;
};
#endif
