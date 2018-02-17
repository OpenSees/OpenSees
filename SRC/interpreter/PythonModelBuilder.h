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
