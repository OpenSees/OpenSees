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

#include "PythonWrapper.h"
#include "OpenSeesCommands.h"
#include <OPS_Globals.h>

static PythonWrapper* wrapper = 0;

PythonWrapper::PythonWrapper()
    :currentArgv(0), currentArg(0), numberArgs(0),
     methodsOpenSees(), opensees_docstring(""), currentResult(0)
{
    wrapper = this;
}

PythonWrapper::~PythonWrapper()
{
    wrapper = 0;
}

void
PythonWrapper::resetCommandLine(int nArgs, int cArg, PyObject* argv)
{
    numberArgs = nArgs;
    currentArg = cArg-1;
    if (currentArg < 0) currentArg = 0;
    currentArgv = argv;
}

void
PythonWrapper::resetCommandLine(int cArg)
{
    if (cArg < 0) {
	currentArg += cArg;
    } else {
	currentArg = cArg-1;
    }
    if (currentArg < 0) currentArg = 0;
}

void
PythonWrapper::addCommand(const char* name, PyCFunction proc)
{
    PyMethodDef method = {name,proc,METH_VARARGS,opensees_docstring};
    methodsOpenSees.push_back(method);
}

PyMethodDef*
PythonWrapper::getMethods()
{
    if (methodsOpenSees.empty()) {
	return 0;
    }

    return &methodsOpenSees[0];
}

void
PythonWrapper::setOutputs(int* data, int numArgs)
{
    if (numArgs == 0) return;
    if (numArgs == 1) {
	currentResult = Py_BuildValue("i", data[0]);
	return ;
    }
    currentResult = PyList_New(numArgs);
    for (int i=0; i<numArgs; i++) {
	PyList_SET_ITEM(currentResult, i, Py_BuildValue("i", data[i]));
    }
}

void
PythonWrapper::setOutputs(double* data, int numArgs)
{
    if (numArgs == 0) return;
    if (numArgs == 1) {
	currentResult = Py_BuildValue("d", data[0]);
	return ;
    }
    currentResult = PyList_New(numArgs);
    for (int i=0; i<numArgs; i++) {
	PyList_SET_ITEM(currentResult, i, Py_BuildValue("d", data[i]));
    }
}

void
PythonWrapper::setOutputs(const char* str)
{
    currentResult = Py_BuildValue("s", str);
}

PyObject*
PythonWrapper::getResults()
{
    PyObject* result = currentResult;
    currentResult = 0;

    if (result == 0) {
	Py_INCREF(Py_None);
	result = Py_None;
    }

    return result;
}

//////////////////////////////////////////////
/////// Python wrapper functions  ////////////
/////////////////////////////////////////////
static PyObject *Py_ops_UniaxialMaterial(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_UniaxialMaterial() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_testUniaxialMaterial(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_testUniaxialMaterial() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setStrain(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setStrain() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getStrain(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getStrain() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getStress(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getStress() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getTangent(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getTangent() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getDampTangent(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getDampTangent() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_wipe(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_wipe() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_model(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_model() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_node(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Node() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_fix(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_HomogeneousBC() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_element(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Element() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_timeSeries(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_TimeSeries() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_pattern(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Pattern() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodalLoad(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_NodalLoad() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_system(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_System() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_numberer(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Numberer() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_constraints(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_ConstraintHandler() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_integrator(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Integrator() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_algorithm(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Algorithm() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_analysis(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Analysis() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_analyze(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_analyze() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_test(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_CTest() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_section(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Section() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_fiber(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Fiber() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_patch(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Patch() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_layer(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Layer() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_geomTransf(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_CrdTransf() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_beamIntegration(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_BeamIntegration() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_loadConst(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_loadConst() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_eleLoad(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_ElementalLoad() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_reactions(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_calculateNodalReactions() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeReaction(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeReaction() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_eigen(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_eigenAnalysis() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nDMaterial(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_NDMaterial() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_block2d(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_doBlock2D() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_block3d(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_doBlock3D() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_rayleigh(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_rayleighDamping() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_wipeAnalysis(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_wipeAnalysis() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setTime(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setTime() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_remove(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_removeObject() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_mass(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_addNodalMass() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_equalDOF(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_EqualDOF() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeEigenvector(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeEigenvector() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getTime(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getTime() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_eleResponse(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_eleResponse() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_SP(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_SP() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_fixX(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_HomogeneousBC_X() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_fixY(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_HomogeneousBC_Y() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_fixZ(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_HomogeneousBC_Z() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_reset(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_resetModel() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_initialize(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_initializeAnalysis() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getLoadFactor(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getLoadFactor() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_build(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_buildModel() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_print(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_printModel() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_printA(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_printA() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_printB(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_printB() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_printGID(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_printModelGID() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getCTestNorms(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getCTestNorms() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getCTestIter(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getCTestIter() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_recorder(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Recorder() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_database(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Database() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_save(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_save() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_restore(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_restore() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_eleForce(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_eleForce() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_eleDynamicalForce(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_eleDynamicalForce() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeUnbalance(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeUnbalance() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeDisp(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeDisp() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeVel(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeVel() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setNodeVel(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setNodeVel() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setNodeDisp(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setNodeDisp() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeAccel(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeAccel() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setNodeAccel(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setNodeAccel() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeResponse(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeResponse() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeCoord(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeCoord() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setNodeCoord(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setNodeCoord() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_updateElementDomain(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_updateElementDomain() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_eleNodes(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_eleNodes() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeDOFs(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeDOFs() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeMass(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeMass() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodePressure(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodePressure() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_nodeBounds(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_nodeBounds() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_startTimer(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_startTimer() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_stopTimer(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_stopTimer() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_modalDamping(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_modalDamping() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_modalDampingQ(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_modalDampingQ() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setElementRayleighDampingFactors(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setElementRayleighDampingFactors() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_region(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_MeshRegion() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setPrecision(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setPrecision() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_searchPeerNGA(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_peerNGA() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_domainChange(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_domainChange() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_record(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_record() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_metaData(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_neesMetaData() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_defaultUnits(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_defaultUnits() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_neesUpload(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_neesUpload() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_stripXML(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_stripOpenSeesXML() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_convertTextToBinary(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_convertTextToBinary() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_convertBinaryToText(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_convertBinaryToText() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getEleTags(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getEleTags() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getNodeTags(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getNodeTags() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getParamTags(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getParamTags() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getParamValue(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getParamValue() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionForce(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionForce() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionDeformation(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionDeformation() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionStiffness(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionStiffness() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionFlexibility(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionFlexibility() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionLocation(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionLocation() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sectionWeight(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sectionWeight() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_basicDeformation(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_basicDeformation() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_basicForce(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_basicForce() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_basicStiffness(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_basicStiffness() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_InitialStateAnalysis(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_InitialStateAnalysis() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_totalCPU(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_totalCPU() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}


static PyObject *Py_ops_solveCPU(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_solveCPU() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_accelCPU(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_accelCPU() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_numFact(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_numFact() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_numIter(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_numIter() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_systemSize(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_systemSize() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_version(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_version() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setMaxOpenFiles(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_maxOpenFiles() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_limitCurve(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_LimitCurve() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_imposedMotion(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_ImposedMotionSP() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_groundMotion(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_groundMotion() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_equalDOF_Mixed(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_EqualDOF_Mixed() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_rigidLink(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_RigidLink() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_rigidDiaphragm(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_RigidDiaphragm() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_ShallowFoundationGen(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_ShallowFoundationGen() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setElementRayleighFactors(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_addElementRayleigh() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_mesh(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_mesh() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_remesh(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_remesh() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_parameter(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_Parameter() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_addToParameter(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_addToParameter() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_updateParameter(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_updateParameter() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setParameter(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setParameter() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getPID(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    //if (OPS_getPID() < 0) {
    // opserr<<(void*)0;
    // return NULL;
    //}

    return wrapper->getResults();
}

static PyObject *Py_ops_getNP(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    //if (OPS_getNP() < 0) {
    // opserr<<(void*)0;
    // return NULL;
    //   }

    return wrapper->getResults();
}

static PyObject *Py_ops_barrier(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    //if (OPS_barrier() < 0) {
    // opserr<<(void*)0;
    // return NULL;
    //}

    return wrapper->getResults();
}

static PyObject *Py_ops_send(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    //if (OPS_send() < 0) {
//     opserr<<(void*)0;
//     return NULL;
// }

    return wrapper->getResults();
}

static PyObject *Py_ops_recv(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    //if (OPS_recv() < 0) {
//     opserr<<(void*)0;
//     return NULL;
// }

    return wrapper->getResults();
}

static PyObject *Py_ops_frictionModel(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_FrictionModel() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_computeGradients(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_computeGradients() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sensitivityAlgorithm(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sensitivityAlgorithm() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sensNodeDisp(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sensNodeDisp() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sensNodeVel(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sensNodeVel() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sensNodeAccel(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sensNodeAccel() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sensLambda(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sensLambda() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sensSectionForce(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sensSectionForce() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sensNodePressure(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sensNodePressure() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_randomVariable(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_randomVariable() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getRVTags(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getRVTags() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getRVMean(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getRVMean() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getRVStdv(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getRVStdv() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getRVPDF(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getRVPDF() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getRVCDF(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getRVCDF() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_getRVInverseCDF(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getRVInverseCDF() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_addCorrelate(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_addCorrelate() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_probabilityTransformation(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_probabilityTransformation() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_transformUtoX(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_transformUtoX() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_updateMaterialStage(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_updateMaterialStage() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_sdfResponse(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sdfResponse() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_getNumThreads(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_getNumThreads() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_setNumThreads(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_setNumThreads() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_logFile(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_logFile() < 0) {
	opserr<<(void*)0;
	return NULL;
    }

    return wrapper->getResults();
}

static PyObject *Py_ops_updateMaterialStage(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_updateMaterialStage() < 0) return NULL;

    return wrapper->getResults();
}

static PyObject *Py_ops_sdfResponse(PyObject *self, PyObject *args)
{
    wrapper->resetCommandLine(PyTuple_Size(args), 1, args);

    if (OPS_sdfResponse() < 0) return NULL;

    return wrapper->getResults();
}

/////////////////////////////////////////////////
////////////// Add Python commands //////////////
/////////////////////////////////////////////////
void
PythonWrapper::addOpenSeesCommands()
{
    addCommand("uniaxialMaterial", &Py_ops_UniaxialMaterial);
    addCommand("testUniaxialMaterial", &Py_ops_testUniaxialMaterial);
    addCommand("setStrain", &Py_ops_setStrain);
    addCommand("getStrain", &Py_ops_getStrain);
    addCommand("getStress", &Py_ops_getStress);
    addCommand("getTangent", &Py_ops_getTangent);
    addCommand("getDampTangent", &Py_ops_getDampTangent);
    addCommand("wipe", &Py_ops_wipe);
    addCommand("model", &Py_ops_model);
    addCommand("node", &Py_ops_node);
    addCommand("fix", &Py_ops_fix);
    addCommand("element", &Py_ops_element);
    addCommand("timeSeries", &Py_ops_timeSeries);
    addCommand("pattern", &Py_ops_pattern);
    addCommand("load", &Py_ops_nodalLoad);
    addCommand("system", &Py_ops_system);
    addCommand("numberer", &Py_ops_numberer);
    addCommand("constraints", &Py_ops_constraints);
    addCommand("integrator", &Py_ops_integrator);
    addCommand("algorithm", &Py_ops_algorithm);
    addCommand("analysis", &Py_ops_analysis);
    addCommand("analyze", &Py_ops_analyze);
    addCommand("test", &Py_ops_test);
    addCommand("section", &Py_ops_section);
    addCommand("fiber", &Py_ops_fiber);
    addCommand("patch", &Py_ops_patch);
    addCommand("layer", &Py_ops_layer);
    addCommand("geomTransf", &Py_ops_geomTransf);
    addCommand("beamIntegration", &Py_ops_beamIntegration);
    addCommand("loadConst", &Py_ops_loadConst);
    addCommand("eleLoad", &Py_ops_eleLoad);
    addCommand("reactions", &Py_ops_reactions);
    addCommand("nodeReaction", &Py_ops_nodeReaction);
    addCommand("eigen", &Py_ops_eigen);
    addCommand("nDMaterial", &Py_ops_nDMaterial);
    addCommand("block2D", &Py_ops_block2d);
    addCommand("block3D", &Py_ops_block3d);
    addCommand("rayleigh", &Py_ops_rayleigh);
    addCommand("wipeAnalysis", &Py_ops_wipeAnalysis);
    addCommand("setTime", &Py_ops_setTime);
    addCommand("remove", &Py_ops_remove);
    addCommand("mass", &Py_ops_mass);
    addCommand("equalDOF", &Py_ops_equalDOF);
    addCommand("nodeEigenvector", &Py_ops_nodeEigenvector);
    addCommand("getTime", &Py_ops_getTime);
    addCommand("eleResponse", &Py_ops_eleResponse);
    addCommand("sp", &Py_ops_SP);
    addCommand("fixX", &Py_ops_fixX);
    addCommand("fixY", &Py_ops_fixY);
    addCommand("fixZ", &Py_ops_fixZ);
    addCommand("reset", &Py_ops_reset);
    addCommand("initialize", &Py_ops_initialize);
    addCommand("getLoadFactor", &Py_ops_getLoadFactor);
    addCommand("build", &Py_ops_build);
    addCommand("printModel", &Py_ops_print);
    addCommand("printA", &Py_ops_printA);
    addCommand("printB", &Py_ops_printB);
    addCommand("printGID", &Py_ops_printGID);
    addCommand("getCTestNorms", &Py_ops_getCTestNorms);
    addCommand("getCTestIter", &Py_ops_getCTestIter);
    addCommand("recorder", &Py_ops_recorder);
    addCommand("database", &Py_ops_database);
    addCommand("save", &Py_ops_save);
    addCommand("restore", &Py_ops_restore);
    addCommand("eleForce", &Py_ops_eleForce);
    addCommand("eleDynamicalForce", &Py_ops_eleDynamicalForce);
    addCommand("nodeUnbalance", &Py_ops_nodeUnbalance);
    addCommand("nodeDisp", &Py_ops_nodeDisp);
    addCommand("setNodeDisp", &Py_ops_setNodeDisp);
    addCommand("nodeVel", &Py_ops_nodeVel);
    addCommand("setNodeVel", &Py_ops_setNodeVel);
    addCommand("nodeAccel", &Py_ops_nodeAccel);
    addCommand("setNodeAccel", &Py_ops_setNodeAccel);
    addCommand("nodeResponse", &Py_ops_nodeResponse);
    addCommand("nodeCoord", &Py_ops_nodeCoord);
    addCommand("setNodeCoord", &Py_ops_setNodeCoord);
    addCommand("updateElementDomain", &Py_ops_updateElementDomain);
    addCommand("eleNodes", &Py_ops_eleNodes);
    addCommand("nodeDOFs", &Py_ops_nodeDOFs);
    addCommand("nodeMass", &Py_ops_nodeMass);
    addCommand("nodePressure", &Py_ops_nodePressure);
    addCommand("nodeBounds", &Py_ops_nodeBounds);
    addCommand("start", &Py_ops_startTimer);
    addCommand("stop", &Py_ops_stopTimer);
    addCommand("modalDamping", &Py_ops_modalDamping);
    addCommand("modalDampingQ", &Py_ops_modalDampingQ);
    addCommand("setElementRayleighDampingFactors", &Py_ops_setElementRayleighDampingFactors);
    addCommand("region", &Py_ops_region);
    addCommand("setPrecision", &Py_ops_setPrecision);
    addCommand("searchPeerNGA", &Py_ops_searchPeerNGA);
    addCommand("domainChange", &Py_ops_domainChange);
    addCommand("record", &Py_ops_record);
    addCommand("metaData", &Py_ops_metaData);
    addCommand("defaultUnits", &Py_ops_defaultUnits);
    addCommand("neesUpload", &Py_ops_neesUpload);
    addCommand("stripXML", &Py_ops_stripXML);
    addCommand("convertBinaryToText", &Py_ops_convertBinaryToText);
    addCommand("convertTextToBinary", &Py_ops_convertTextToBinary);
    addCommand("getEleTags", &Py_ops_getEleTags);
    addCommand("getNodeTags", &Py_ops_getNodeTags);
    addCommand("getParamTags", &Py_ops_getParamTags);
    addCommand("getParamValue", &Py_ops_getParamValue);
    addCommand("sectionForce", &Py_ops_sectionForce);
    addCommand("sectionDeformation", &Py_ops_sectionDeformation);
    addCommand("sectionStiffness", &Py_ops_sectionStiffness);
    addCommand("sectionFlexibility", &Py_ops_sectionFlexibility);
    addCommand("sectionLocation", &Py_ops_sectionLocation);
    addCommand("sectionWeight", &Py_ops_sectionWeight);
    addCommand("basicDeformation", &Py_ops_basicDeformation);
    addCommand("basicForce", &Py_ops_basicForce);
    addCommand("basicStiffness", &Py_ops_basicStiffness);
    addCommand("InitialStateAnalysis", &Py_ops_InitialStateAnalysis);
    addCommand("totalCPU", &Py_ops_totalCPU);
    addCommand("solveCPU", &Py_ops_solveCPU);
    addCommand("accelCPU", &Py_ops_accelCPU);
    addCommand("numFact", &Py_ops_numFact);
    addCommand("numIter", &Py_ops_numIter);
    addCommand("systemSize", &Py_ops_systemSize);
    addCommand("version", &Py_ops_version);
    addCommand("setMaxOpenFiles", &Py_ops_setMaxOpenFiles);
    addCommand("limitCurve", &Py_ops_limitCurve);
    addCommand("imposedMotion", &Py_ops_imposedMotion);
    addCommand("imposedSupportMotion", &Py_ops_imposedMotion);
    addCommand("groundMotion", &Py_ops_groundMotion);
    addCommand("equalDOF_Mixed", &Py_ops_equalDOF_Mixed);
    addCommand("rigidLink", &Py_ops_rigidLink);
    addCommand("rigidDiaphragm", &Py_ops_rigidDiaphragm);
    addCommand("ShallowFoundationGen", &Py_ops_ShallowFoundationGen);
    addCommand("setElementRayleighFactors", &Py_ops_setElementRayleighFactors);
    addCommand("mesh", &Py_ops_mesh);
    addCommand("remesh", &Py_ops_remesh);
    addCommand("parameter", &Py_ops_parameter);
    addCommand("addToParameter", &Py_ops_addToParameter);
    addCommand("updateParameter", &Py_ops_updateParameter);
    addCommand("setParameter", &Py_ops_setParameter);
    addCommand("getPID", &Py_ops_getPID);
    addCommand("getNP", &Py_ops_getNP);
    addCommand("barrier", &Py_ops_barrier);
    addCommand("send", &Py_ops_send);
    addCommand("recv", &Py_ops_recv);
    addCommand("frictionModel", &Py_ops_frictionModel);
    addCommand("computeGradients", &Py_ops_computeGradients);
    addCommand("sensitivityAlgorithm", &Py_ops_sensitivityAlgorithm);
    addCommand("sensNodeDisp", &Py_ops_sensNodeDisp);
    addCommand("sensNodeVel", &Py_ops_sensNodeVel);
    addCommand("sensNodeAccel", &Py_ops_sensNodeAccel);
    addCommand("sensLambda", &Py_ops_sensLambda);
    addCommand("sensSectionForce", &Py_ops_sensSectionForce);
    addCommand("sensNodePressure", &Py_ops_sensNodePressure);
    addCommand("randomVariable", &Py_ops_randomVariable);
    addCommand("getRVTags", &Py_ops_getRVTags);
    addCommand("getMean", &Py_ops_getRVMean);
    addCommand("getStdv", &Py_ops_getRVStdv);
    addCommand("getPDF", &Py_ops_getRVPDF);
    addCommand("getCDF", &Py_ops_getRVCDF);
    addCommand("getInverseCDF", &Py_ops_getRVInverseCDF);
    addCommand("correlate", &Py_ops_addCorrelate);
    addCommand("transformUtoX", &Py_ops_transformUtoX);
    addCommand("updateMaterialStage", &Py_ops_updateMaterialStage);
    addCommand("sdfResponse", &Py_ops_sdfResponse);
    addCommand("probabilityTransformation", &Py_ops_probabilityTransformation);
    addCommand("getNumThreads", &Py_ops_getNumThreads);
    addCommand("setNumThreads", &Py_ops_setNumThreads);
    addCommand("logFile", &Py_ops_logFile);


    PyMethodDef method = {NULL,NULL,0,NULL};
    methodsOpenSees.push_back(method);
}
