/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision$
// $Date$
// $Source$

#include <PythonRV.h>
#include <NormalRV.h>
#include <cmath>

#ifdef _DEBUG
#undef _DEBUG
#include <Python.h>
#define _DEBUG
#else
#include <Python.h>
#endif

PythonRV::PythonRV(int passedTag, double passeda, double passedb,
		   const char *filename, const char *functionname)
  :RandomVariable(passedTag, RANDOM_VARIABLE_python),
   myFunction(0), a(passeda), b(passedb)
{
  PyObject *myModule = PyImport_ImportModule(filename);
  if (myModule == 0) {
    opserr << "PythonRV::PythonRV - unable to import module " << filename << endln;
    PyErr_Print();
  }
  
  myFunction = PyObject_GetAttrString(myModule,(char*)functionname);
}


PythonRV::PythonRV(int passedTag, const Vector &passedParameters,
		   const char *filename, const char *functionname)
  :RandomVariable(passedTag, RANDOM_VARIABLE_python), myFunction(0),
   a(0.0), b(0.0)
{
  PyObject *myModule = PyImport_ImportModule(filename);
  if (myModule == 0) {
    opserr << "PythonRV::PythonRV - unable to import module " << filename << endln;
    PyErr_Print();
  }
    
  if (passedParameters.Size() != 2) {
    opserr << "Python RV requires 2 parameters, mean and stdv, for RV with tag " << this->getTag() << endln;
    
    a = 1.0;
    b = 0.5;
    
  } else {
    a = passedParameters(0);
    b = passedParameters(1);
  }

  myFunction = PyObject_GetAttrString(myModule,(char*)functionname);
}


PythonRV::~PythonRV()
{
  Py_DECREF(myFunction);
}


const char *
PythonRV::getType()
{
  return "GAMMA";
}


double 
PythonRV::getMean()
{
  double result;

  //
  // Call Python
  //
  PyObject *myResult = PyObject_CallFunction(myFunction,"i d d d",4,0.0,a,b);

  result = PyFloat_AsDouble(myResult);

  Py_DECREF(myResult);
  
  return result;
}


double 
PythonRV::getStdv()
{
  double result;

  //
  // Call Python
  //
  PyObject *myResult = PyObject_CallFunction(myFunction,"i d d d",5,0.0,a,b);

  result = PyFloat_AsDouble(myResult);

  Py_DECREF(myResult);
  
  return result;
}


const Vector &
PythonRV::getParameters(void)
{
  static Vector temp(2);
  temp(0) = a;
  temp(1) = b;
  return temp;
}


int 
PythonRV::setParameters(double passedMean, double passedStdv)
{
  a = passedMean;
  b = passedStdv;
	
  return 0;
}


double
PythonRV::getPDFvalue(double rvValue)
{
  double result;

  //
  // Call Python
  //
  PyObject *myResult = PyObject_CallFunction(myFunction,"i d d d",1,rvValue,a,b);

  result = PyFloat_AsDouble(myResult);

  Py_DECREF(myResult);
  
  return result;
}

double
PythonRV::getCDFvalue(double rvValue)
{
  double result;

  //
  // Call Python
  //
  PyObject *myResult = PyObject_CallFunction(myFunction,"i d d d",2,rvValue,a,b);

  result = PyFloat_AsDouble(myResult);

  Py_DECREF(myResult);
  
  return result;
}


double
PythonRV::getInverseCDFvalue(double probValue)
{
  double result;

  //
  // Call Python
  //
  PyObject *myResult = PyObject_CallFunction(myFunction,"i d d d",3,probValue,a,b);

  result = PyFloat_AsDouble(myResult);

  Py_DECREF(myResult);
    
  return result;
}


int 
PythonRV::getCDFparameterSensitivity(Vector &dFdP)
{
  /*
  // returns gradient of F(x) with respect to distribution parameters
  double rvValue = this->getCurrentValue();
  double cdfValue = getCDFvalue(rvValue);
    
  // dFdk
  // no closed form expression for this one, use finite differences
  double k_old = k;
  double dh = k/1000.0;
  k += dh;
  dFdP(0) = ( getCDFvalue(rvValue) - cdfValue )/dh;
  k = k_old;
  
  // dFdlambda
  dFdP(1) = rvValue/lambda * getPDFvalue(rvValue);
  */
  
  return 0;
}


int
PythonRV::getParameterMeanSensitivity(Vector &dPdmu)
{
  /*
  // returns gradient of distribution parameters with respect to the mean
  double mu = getMean();
  double sig = getStdv();
  
  // dkdmu
  dPdmu(0) = 2*mu/sig/sig;
    
  // dlambdadmu
  dPdmu(1) = 1/sig/sig;
  */
  
  return 0;
}


int
PythonRV::getParameterStdvSensitivity(Vector &dPdstdv)
{
  /*
  // returns gradient of distribution parameters with respect to the stdv
  double mu = getMean();
  double sig = getStdv();
  
  // dkdsig
  dPdstdv(0) = -2*mu*mu/sig/sig/sig;
  
  // dlambdadsig
  dPdstdv(1) = -2*mu/sig/sig/sig;
  */
  
  return 0;
}


void
PythonRV::Print(OPS_Stream &s, int flag)
{
  s << "Python RV #" << this->getTag() << endln;
  s << "\ta = " << a << endln;
  s << "\tb = " << b << endln;
}
