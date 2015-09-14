/* ****************************************************************** **
**    Opensees - Open System for Earthquake Engineering Simulation    **
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
                                                                        
// $Revision: 1.0 $
// $Date: 2015-08-25 11:50:00 $
                                                                        
// Written: Minjie Zhu
//
// Description: This file contains the class definition for PythonModelBuilder.
// A PythonModelBuilder is a python module that adds the commands to create
// the model for the standard

#ifndef PythonModelBuilder_h
#define PythonModelBuilder_h

#include <Python.h>
#include <Domain.h>

class LoadPattern;
class MultiSupportPattern;
class SectionForceDeformation;

class PythonModelBuilder
{
public:
    PythonModelBuilder();
    ~PythonModelBuilder();

    int getNDM() const {return ndm;}
    int getNDF() const {return ndf;}
    void set(int n, int f) {ndm=n; ndf=f;isset=true;}
    bool isSet() const {return isset;}

    void setCurrentLoadPattern(LoadPattern* pattern) {currentLoadPattern = pattern;}
    LoadPattern* getCurrentLoadPattern() {return currentLoadPattern;}

    void setCurrentMultiPattern(MultiSupportPattern* pattern) {currentMultiPattern = pattern;}
    MultiSupportPattern* getCurrentMultiPattern() {return currentMultiPattern;}

    void setCurrentSection(SectionForceDeformation* section) {currentSection = section;}
    SectionForceDeformation* getCurrentSection() const {return currentSection;}
    
private:
    int ndm;
    int ndf;
    bool isset;
    LoadPattern* currentLoadPattern;
    MultiSupportPattern* currentMultiPattern;
    SectionForceDeformation* currentSection;
};

PythonModelBuilder& OPS_GetPythonModelBuilder();
PyObject *ops_Model(PyObject *self, PyObject *args);
PyObject *ops_addNode(PyObject *self, PyObject *args);
PyObject *ops_addHomogeneousBC(PyObject *self, PyObject *args);
PyObject *ops_addElement(PyObject *self, PyObject *args);
PyObject *ops_addTimeSeries(PyObject *self, PyObject *args);
PyObject *ops_addPattern(PyObject *self, PyObject *args);
PyObject *ops_addNodalLoad(PyObject *self, PyObject *args);
PyObject *ops_addSP(PyObject *self, PyObject *args);
PyObject *ops_addSection(PyObject *self, PyObject *args);
PyObject *ops_addNDMaterial(PyObject *self, PyObject *args);
PyObject *ops_addFiber(PyObject *self, PyObject *args);
PyObject *ops_addPatch(PyObject *self, PyObject *args);
PyObject *ops_addLayer(PyObject *self, PyObject *args);
PyObject *ops_addGeomTransf(PyObject *self, PyObject *args);
PyObject *ops_addBeamIntegrationRule(PyObject *self, PyObject *args);
PyObject *ops_addNodalMass(PyObject *self, PyObject *args);
PyObject *ops_addGroundMotion(PyObject *self, PyObject *args);
PyObject *ops_addElementalLoad(PyObject *self, PyObject *args);

#endif
