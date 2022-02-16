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
#include <OPS_Globals.h>
#include <cstring>
#include <cctype>

// define opserr
static PythonStream sserr;
OPS_Stream *opserrPtr = &sserr;


PythonModule::PythonModule()
        : wrapper(), cmds(this) {
    // does nothing
}

PythonModule::~PythonModule() {
    // does nothing
}

int
PythonModule::run() {
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

const char *PythonModule::trimSpaces(PyObject *o) {
  Py_ssize_t size = 0;
  const char *s = PyUnicode_AsUTF8AndSize(o, &size);

  // empty string
  if (size == 0) {
    return s;
  }

  // find first char that is not space
  Py_ssize_t firstLoc = 0;
  for (Py_ssize_t i = 0; i < size; ++i) {
    if (isspace(s[i])) {
      firstLoc = i + 1;
    } else {
      break;
    }
  }

  // find last char that is not space
  Py_ssize_t lastLoc = size - 1;
  for (Py_ssize_t i = size - 1; i >= 0; --i) {
    if (isspace(s[i])) {
      lastLoc = i - 1;
    } else {
      break;
    }
  }

  // new Py Object
  PyObject *newo = 0;

  // all spaces: make an empty string
  if (firstLoc == size || lastLoc < 0) {
    newo = PyUnicode_FromString("");
  }

  // get substring
  else if (firstLoc > 0 || lastLoc < size - 1) {
    newo = PyUnicode_Substring(o, firstLoc, lastLoc + 1);
  }

  // get string
  if (newo == 0) {
    return s;
  }

  const char* news = PyUnicode_AsUTF8(newo);
  Py_DECREF(newo);
  return news;
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
PythonModule::getDouble(double *data, int numArgs) {
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

const char *
PythonModule::getString() {
    if (wrapper.getCurrentArg() >= wrapper.getNumberArgs()) {
        return 0;
    }

    PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
    wrapper.incrCurrentArg();
#if PY_MAJOR_VERSION >= 3
    if (!PyUnicode_Check(o)) {
        return 0;
    }

    // PyObject* space = PyUnicode_FromString(" ");
    // PyObject* empty = PyUnicode_FromString("");
    // PyObject* newo = PyUnicode_Replace(o, space, empty, -1);
    const char* res = trimSpaces(o);
    // Py_DECREF(newo);
    // Py_DECREF(space);
    // Py_DECREF(empty);

    return res;
#else
    if (!PyString_Check(o)) {
        return 0;
    }

    return PyString_AS_STRING(o);
#endif
}

const char *PythonModule::getStringFromAll(char* buffer, int len) {
    if (wrapper.getCurrentArg() >= wrapper.getNumberArgs()) {
        return 0;
    }

    PyObject *o =
        PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
    wrapper.incrCurrentArg();

    // check if int
    if (PyLong_Check(o) || PyBool_Check(o)) {
        PyErr_Clear();
        int data = PyLong_AsLong(o);
        if (PyErr_Occurred()) {
            return 0;
        }
        snprintf(buffer, len, "%d", data);
        return buffer;
    }
    // check if double
    else if (PyFloat_Check(o)) {
        PyErr_Clear();
        double data = PyFloat_AsDouble(o);
        if (PyErr_Occurred()) {
            return 0;
        }
        snprintf(buffer, len, "%.20f", data);
        return buffer;
    }

#if PY_MAJOR_VERSION >= 3
    if (!PyUnicode_Check(o)) {
        return 0;
    }

    // PyObject *space = PyUnicode_FromString(" ");
    // PyObject *empty = PyUnicode_FromString("");
    // PyObject *newo = PyUnicode_Replace(o, space, empty, -1);
    const char* res = trimSpaces(o);
    // Py_DECREF(newo);
    // Py_DECREF(space);
    // Py_DECREF(empty);

    int lenres = strlen(res) + 1;
    if (lenres > len) {
        lenres = len;
    }

    strncpy(buffer, res, lenres);

    return buffer;
#else
    if (!PyString_Check(o)) {
        return 0;
    }

    return PyString_AS_STRING(o);
#endif
}

int
PythonModule::getStringCopy(char **stringPtr) {
    return -1;
}

void
PythonModule::resetInput(int cArg) {
    wrapper.resetCommandLine(cArg);
}

int
PythonModule::setInt(int *data, int numArgs, bool scalar) {
    wrapper.setOutputs(data, numArgs, scalar);

    return 0;
}

int
PythonModule::setDouble(double *data, int numArgs, bool scalar) {
    wrapper.setOutputs(data, numArgs, scalar);

    return 0;
}

int
PythonModule::setString(const char *str) {
    wrapper.setOutputs(str);

    return 0;
}

int
PythonModule::runCommand(const char *cmd) {
    return PyRun_SimpleString(cmd);
}

static PythonModule *module = 0;

PyMethodDef *getmethodsFunc() {
    module = new PythonModule;
    PythonWrapper *wrapper = module->getWrapper();
    wrapper->addOpenSeesCommands();

    return wrapper->getMethods();
}

void cleanupFunc() {
    module->getCmds().wipe();
    if (module != 0) {
        delete module;
    }
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

static int opensees_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);

    return 0;
}

static int opensees_clear(PyObject *m) {
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
    PyObject *pymodule = PyModule_Create(&moduledef);
#else
    PyObject *pymodule = Py_InitModule("opensees", getmethodsFunc());
#endif

    if (pymodule == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(pymodule);

    // add OpenSeesError
    st->error = PyErr_NewExceptionWithDoc("opensees.OpenSeesError", "Internal OpenSees errors.", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(pymodule);
        INITERROR;
    }
    Py_INCREF(st->error);
    PyModule_AddObject(pymodule, "OpenSeesError", st->error);

    // add OpenSeesParameter dict
    PyObject *par = PyDict_New();
    if (par == NULL) {
      INITERROR;
    }
    if (PyModule_AddObject(pymodule, "OpenSeesParameter", par) < 0) {
        Py_DECREF(par);
        INITERROR;
    }

    // char version[10];
    // const char *py_version = ".6";
    // for (int i = 0; i < 5; ++i) {
    //     version[i] = OPS_VERSION[i];
    // }
    // for (int i = 0; i < 3; ++i) {
    //     version[5 + i] = py_version[i];
    // }
    // PyModule_AddStringConstant(pymodule, "__version__", version);

    sserr.setError(st->error);

    Py_AtExit(cleanupFunc);

#if PY_MAJOR_VERSION >= 3
    return pymodule;
#endif
}
