#include <Python.h>
#include <elementAPI.h>
#include <ElasticMaterial.cpp>

#include <UniaxialMaterial.h>
#include <string.h>

#include <StandardStream.h>
StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

#include <SimulationInformation.h>
SimulationInformation simulationInfo;
SimulationInformation *theSimulationInfoPtr = 0;

#include <FE_Datastore.h>
FE_Datastore *theDatabase  =0;


static UniaxialMaterial *theTestingUniaxialMaterial =0;
static double ops_strain = 0;

extern UniaxialMaterial *OPS_ParseUniaxialMaterialCommand(const char *matType);
extern "C" void OPS_ResetCommandLine(int nArgs, int cArg, PyObject *args);

static PyObject *ops_UniaxialMaterialCommand(PyObject *self, PyObject *args)
{
  OPS_ResetCommandLine(PyTuple_Size(args), 0, args);
  int numberArgs = OPS_GetNumRemainingInputArgs();

  if (numberArgs < 2) {
    PyErr_SetString(PyExc_RuntimeError, "ERROR too few arguments: uniaxialMaterial type? tag? args.");
    return NULL;
  }

  const char *matType = OPS_GetString();
  UniaxialMaterial *theMaterial = OPS_ParseUniaxialMaterialCommand(matType);

  if (theMaterial == 0) {
    PyErr_SetString(PyExc_RuntimeError, "ERROR could not create uniaxialMaterial.");
    return NULL;
  }
  
  // Now add the material to the modelBuilder
  if (OPS_addUniaxialMaterial(theMaterial) == false) {
    PyErr_SetString(PyExc_RuntimeError, "ERROR could not add uniaaialMaterial.");
    delete theMaterial; // invoke the material objects destructor, otherwise mem leak
    return NULL;
  }

  PyObject *ret = Py_BuildValue("d", 0.0);
  return ret;  
}

static PyObject *ops_testUniaxialMaterial(PyObject *self, PyObject *args)
{
  OPS_ResetCommandLine(PyTuple_Size(args), 0, args);
  int numberArgs = OPS_GetNumRemainingInputArgs();

  if (numberArgs !=  1) {
    PyErr_SetString(PyExc_RuntimeError, "testUniaxialMaterial - You must provide a material tag.");
    return NULL;
  }

  int tag;
  int numData = 1;
  OPS_GetIntInput(&numData, &tag);

  if (theTestingUniaxialMaterial != 0) 
    delete theTestingUniaxialMaterial;

  theTestingUniaxialMaterial=OPS_getUniaxialMaterial(tag);

  if (theTestingUniaxialMaterial == 0) {
      PyErr_SetString(PyExc_ValueError,"testUniaxialMaterial - Material Not Found.");
      return NULL;
  }
  PyObject *ret = Py_BuildValue("d", 0.0);
  return ret;  
}


static PyObject *ops_setStrain(PyObject *self, PyObject *args)
{
  int argc = PyTuple_Size(args);
  if (argc == 0) {
    return NULL;
    PyErr_SetString(PyExc_RuntimeError, "You must provide a strain value.");
  }

  PyObject *o = PyTuple_GetItem(args,0);
  if (!PyFloat_Check(o)) {
    PyErr_SetString(PyExc_ValueError,"setStrain expected a floating point value");
    return NULL;
  }

  double opsStrain = PyFloat_AS_DOUBLE(o);
  if (theTestingUniaxialMaterial !=0 ) {
    theTestingUniaxialMaterial->setTrialStrain(opsStrain);
    theTestingUniaxialMaterial->commitState();
  } else {
    PyErr_SetString(PyExc_ValueError,"setStrain WARNING no active UniaxialMaterial - use testUniaxialMaterial command");    
    return NULL;
  }
  PyObject *ret = Py_BuildValue("d", opsStrain);
  return ret;
}


static PyObject *ops_getStrain(PyObject *self, PyObject *args)
{
  if (theTestingUniaxialMaterial !=0) {
    double strain = theTestingUniaxialMaterial->getStrain();
    PyObject *ret = Py_BuildValue("d", strain);
    return ret;
  } else {
    PyErr_SetString(PyExc_ValueError,"getStrain no active uniaxialMaerial - use testUniaxialMaterial function");
    return NULL;
  }
}


static PyObject *ops_getStress(PyObject *self, PyObject *args)
{
  if (theTestingUniaxialMaterial !=0) {
    double stress = theTestingUniaxialMaterial->getStress();
    PyObject *ret = Py_BuildValue("d", stress);
    return ret;
  } else {
    PyErr_SetString(PyExc_ValueError,"getStress no active uniaxialMaerial - use testUniaxialMaterial function");
    return NULL;
  }
}


static PyObject *ops_getTangent(PyObject *self, PyObject *args)
{
  if (theTestingUniaxialMaterial !=0) {
    double tangent = theTestingUniaxialMaterial->getTangent();
    PyObject *ret = Py_BuildValue("d", tangent);
    return ret;
  } else {
    PyErr_SetString(PyExc_ValueError,"getTangent no active uniaxialMaerial - use testUniaxialMaterial function");
    return NULL;
  }
}

static char module_docstring[] =
  "This module add OpenSees Commands to Python";
static char opensees_docstring[] =
  "Blah";

static PyMethodDef methodsOpenSees[] = {
  {"uniaxialMaterial", ops_UniaxialMaterialCommand, METH_VARARGS, opensees_docstring},
  {"testUniaxialMaterial", ops_testUniaxialMaterial, METH_VARARGS, opensees_docstring},
  {"setStrain", ops_setStrain, METH_VARARGS, opensees_docstring},
  {"getStrain", ops_getStrain, METH_VARARGS, opensees_docstring},
  {"getStress", ops_getStress, METH_VARARGS, opensees_docstring},
  {"getTangent", ops_getTangent, METH_VARARGS, opensees_docstring},
  {NULL, NULL, 0, NULL}
};


void addOpenSeesCommands(void) {
  PyImport_AddModule("opensees");
  Py_InitModule("opensees", methodsOpenSees);
}


/* ********************************* *
 * elementAPI FUNCTIONS for Python   *
 *********************************** */

static int numberArgs = 0;
static int currentArg = 0;
PyObject *currentArgs = 0;

extern "C"
void OPS_ResetCommandLine(int nArgs, int cArg, PyObject *args) {
  numberArgs = nArgs;
  currentArg = cArg;
  currentArgs = args;
}


extern "C"   
int OPS_GetNumRemainingInputArgs()
{
  return numberArgs - currentArg;
}

extern "C"   
int OPS_GetIntInput(int *numData, int*data)
{
  if ((numberArgs - currentArg) >= *numData) {
    for (int i=0; i<*numData; i++) {
      PyObject *o = PyTuple_GetItem(currentArgs,currentArg);
      if (!PyInt_Check(o)) {
	PyErr_SetString(PyExc_ValueError,"invalid integer");
	return -1;
      }
      long value = PyInt_AS_LONG(o);
      data[i] = value;
      currentArg++;
    }
  
    return 0;
  }
  return -1;
}


extern "C" 
int OPS_GetDoubleInput(int *numData, double *data)
{
  int numD = *numData;
  if ((numberArgs - currentArg) >= numD) {
    for (int i=0; i<numD; i++) {
      PyObject *o = PyTuple_GetItem(currentArgs,currentArg);
      if (!PyFloat_Check(o)) {
	PyErr_SetString(PyExc_ValueError,"invalid integer");
	return -1;
      }
      data[i] = PyFloat_AS_DOUBLE(o);
      currentArg++;
    }
    return 0;  
  }

  return -1;
}

extern "C" 
const char *OPS_GetString(void)
{
  PyObject *o = PyTuple_GetItem(currentArgs,currentArg);
  currentArg++;
  return PyString_AS_STRING(o);
}
