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

// Description: This file contains the class definition for Python Module
// PythonModule implements a DL_Interpreter for the python language
//

#include "PythonModule.h"
#include "PythonStream.h"

// define opserr
static PythonStream sserr;
OPS_Stream *opserrPtr = &sserr;


PythonModule::PythonModule()
    :wrapper(), cmds(this)
{
    // does nothing
}

PythonModule::~PythonModule()
{
    // does nothing
}

int
PythonModule::run()
{
    return 0;
}

int
PythonModule::addCommand(const char *, Command &)
{
    return -1;
}

int
PythonModule::removeCommand(const char *)
{
    return -1;
}

int
PythonModule::getNumRemainingInputArgs(void)
{
    return wrapper.getNumberArgs() - wrapper.getCurrentArg();
}

int
PythonModule::getInt(int *data, int numArgs)
{
    if ((wrapper.getNumberArgs() - wrapper.getCurrentArg()) < numArgs) {
        return -1;
    }

    for (int i = 0; i < numArgs; i++) {
        PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
        wrapper.incrCurrentArg();
#if PY_MAJOR_VERSION >= 3
        if (!PyLong_Check(o)) {
            return -1;
        }
        data[i] = PyLong_AS_LONG(o);
#else
        if (!PyInt_Check(o)) {
            return -1;
        }
        data[i] = PyInt_AS_LONG(o);
#endif
    }

    return 0;
}

int
PythonModule::getDouble(double *data, int numArgs)
{
    if ((wrapper.getNumberArgs() - wrapper.getCurrentArg()) < numArgs) {
        return -1;
    }

    for (int i = 0; i < numArgs; i++) {
        PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
        wrapper.incrCurrentArg();
        if (!PyFloat_Check(o)) {
            return -1;
        }
        data[i] = PyFloat_AS_DOUBLE(o);
    }

    return 0;
}

const char*
PythonModule::getString()
{
    if (wrapper.getCurrentArg() >= wrapper.getNumberArgs()) {
        return 0;
    }

    PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
    wrapper.incrCurrentArg();
#if PY_MAJOR_VERSION >= 3
    if (!PyUnicode_Check(o)) {
        return 0;
    }

    return PyUnicode_AsUTF8(o);
#else
    if (!PyString_Check(o)) {
        return 0;
    }

    return PyString_AS_STRING(o);
#endif
}

int
PythonModule::getStringCopy(char **stringPtr)
{
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

int
PythonModule::runCommand(const char* cmd)
{
    return PyRun_SimpleString(cmd);
}

static PythonModule* module = 0;

PyMethodDef* getmethodsFunc()
{
    module = new PythonModule;
    PythonWrapper* wrapper = module->getWrapper();
    wrapper->addOpenSeesCommands();
    
    return wrapper->getMethods();
}

void cleanupFunc()
{
    module->getCmds().wipe();
    // if (module != 0) {
    //     delete module;
    // }
}

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

// static PyObject *
// error_out(PyObject *m)
// {
//     struct module_state *st = GETSTATE(m);
//     PyErr_SetString(st->error, "something bad happened");
//
//     return NULL;
// }

#if PY_MAJOR_VERSION >= 3

static int opensees_traverse(PyObject *m, visitproc visit, void *arg)
{
    Py_VISIT(GETSTATE(m)->error);

    return 0;
}

static int opensees_clear(PyObject *m)
{
    Py_CLEAR(GETSTATE(m)->error);

    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "opensees",
    NULL,
    sizeof(struct module_state),
    getmethodsFunc(),
    NULL,
    opensees_traverse,
    opensees_clear,
    NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_opensees(void)

#else
#define INITERROR return

//void
PyMODINIT_FUNC
initopensees(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("opensees", getmethodsFunc());
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("opensees.error", NULL, NULL);
    PyObject* ops_msg = PyErr_NewException("opensees.msg", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

    sserr.setError(st->error,ops_msg);

    Py_AtExit(cleanupFunc);

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
