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
                                                                        
// Written: Minjie

// Description: A python wrapper for OpenSees commands
//
#ifndef PythonWrapper_h
#define PythonWrapper_h

#include <Python.h>
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
    void setOutputs(int* data, int numArgs);
    void setOutputs(double* data, int numArgs);
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
