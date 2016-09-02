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

// Description: This file contains the class definition for Python Module
// PythonModule implements a DL_Interpreter for the python language
//

#include "PythonModule.h"


PythonModule::PythonModule()
    :wrapper(), cmds(this)
{
}

PythonModule::~PythonModule()
{
  // does nothing
}

int
PythonModule::run(){
  return 0; 
}
			 

int 
PythonModule::addCommand(const char *, Command &) {
  return -1;
}

int 
PythonModule::removeCommand(const char *) {
  return -1;
}

int 
PythonModule::getNumRemainingInputArgs(void) {
    return wrapper.getNumberArgs() - wrapper.getCurrentArg();
}

int 
PythonModule::getInt(int *data, int numArgs) {
    if ((wrapper.getNumberArgs() - wrapper.getCurrentArg()) < numArgs) {
	return -1;
    }

    for (int i=0; i<numArgs; i++) {
	PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(),wrapper.getCurrentArg());
	wrapper.incrCurrentArg();
	if (!PyInt_Check(o)) {
	    return -1;
	}
	data[i] = PyInt_AS_LONG(o);
    }
    
    return 0;
}

int 
PythonModule::getDouble(double *data, int numArgs) {
    if ((wrapper.getNumberArgs() - wrapper.getCurrentArg()) < numArgs) {
	return -1;
    }

    for (int i=0; i<numArgs; i++) {
	PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(),wrapper.getCurrentArg());
	wrapper.incrCurrentArg();
	if (!PyFloat_Check(o)) {
	    return -1;
	}
	data[i] = PyFloat_AS_DOUBLE(o);
    }
    
    return 0;
}

const char* 
PythonModule::getString() {
    
    if (wrapper.getCurrentArg() >= wrapper.getNumberArgs()) {
	return 0;
    }
    
    PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(),wrapper.getCurrentArg());
    wrapper.incrCurrentArg();
    if (!PyString_Check(o)) {
	return 0;
    }
    
    
    return PyString_AS_STRING(o);;
}

int 
PythonModule::getStingCopy(char **stringPtr) {
  return -1;
}

void
PythonModule::resetInput(int cArg)
{
    wrapper.resetCommandLine(cArg);
}

int
PythonModule::setInt(int* data, int numArgs)
{
    wrapper.setOutputs(data, numArgs);
    return 0;
}

int
PythonModule::setDouble(double* data, int numArgs)
{
    wrapper.setOutputs(data, numArgs);
    return 0;
}

int
PythonModule::setString(const char* str)
{
    wrapper.setOutputs(str);
    return 0;
}

static PythonModule* module = 0;

void cleanupFunc()
{
    if (module != 0) {
	delete module;
    }
}

PyMODINIT_FUNC
initopensees(void)
{
    // create python interpreter
    module = new PythonModule();
    PythonWrapper* wrapper = module->getWrapper();
    wrapper->addOpenSeesCommands();
    Py_InitModule("opensees", wrapper->getMethods());

    // set up cleanup function at exit
    Py_AtExit(cleanupFunc);
}
