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

                                                                        
// Written: Minjie Zhu
//

#ifndef PythonAnalysisBulider_h
#define PythonAnalysisBulider_h

#include <Python.h>

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class EquiSolnAlgo;
class LinearSOE;
class EigenSOE;
class StaticIntegrator;
class TransientIntegrator;
class ConvergenceTest;
class StaticAnalysis;
class DirectIntegrationAnalysis;
class VariableTimeStepDirectIntegrationAnalysis;
class PFEMAnalysis;
class Integrator;

class PythonAnalysisBuilder
{
public:
    PythonAnalysisBuilder();
    ~PythonAnalysisBuilder();

    void set(ConstraintHandler* obj);
    void set(DOF_Numberer* obj);
    void set(EquiSolnAlgo* obj);
    void set(LinearSOE* obj);
    void set(Integrator* obj, int isstatic);
    void set(ConvergenceTest* obj);
    
    void newStaticAnalysis();
    void newTransientAnalysis();
    void newPFEMAnalysis(double,double,double);
    void newEigenAnalysis(int typeSolver, double shift);

    StaticAnalysis* getStaticAnalysis() {return theStaticAnalysis;}
    DirectIntegrationAnalysis* getTransientAnalysis() {return theTransientAnalysis;}
    PFEMAnalysis* getPFEMAnalysis() {return thePFEMAnalysis;}
    VariableTimeStepDirectIntegrationAnalysis* getVariableTimeStepDirectIntegrationAnalysis() {
	return theVariableTimeStepTransientAnalysis;
    }
    EquiSolnAlgo* getAlgorithm();
    StaticIntegrator* getStaticIntegrator();
    TransientIntegrator* getTransientIntegrator();
    ConvergenceTest  *getConvergenceTest();

    void wipe();
    void resetStatic();
    void resetTransient();
    void resetAll();
    
private:
    ConstraintHandler 	*theHandler;    
    DOF_Numberer 	*theNumberer;
    AnalysisModel 	*theAnalysisModel;
    EquiSolnAlgo 	*theAlgorithm;
    LinearSOE 		*theSOE;
    EigenSOE 		*theEigenSOE;
    StaticIntegrator    *theStaticIntegrator;
    TransientIntegrator *theTransientIntegrator;
    ConvergenceTest     *theTest;
    StaticAnalysis      *theStaticAnalysis;
    DirectIntegrationAnalysis *theTransientAnalysis;
    VariableTimeStepDirectIntegrationAnalysis *theVariableTimeStepTransientAnalysis;
    PFEMAnalysis        *thePFEMAnalysis;

};

PyObject *ops_printModel(PyObject *self, PyObject *args);
PyObject *ops_wipeAnalysis(PyObject *self, PyObject *args);
PyObject *ops_wipeModel(PyObject *self, PyObject *args);
PyObject *ops_specifySOE(PyObject *self, PyObject *args);
PyObject *ops_specifyNumberer(PyObject *self, PyObject *args);
PyObject *ops_specifyConstraintHandler(PyObject *self, PyObject *args);
PyObject *ops_specifyAlgorithm(PyObject *self, PyObject *args);
PyObject *ops_specifyIntegrator(PyObject *self, PyObject *args);
PyObject *ops_specifyAnalysis(PyObject *self, PyObject *args);
PyObject *ops_specifyCTest(PyObject *self, PyObject *args);
PyObject *ops_analyzeModel(PyObject *self, PyObject *args);
PyObject *ops_nodeDisp(PyObject *self, PyObject *args);
PyObject *ops_getLoadFactor(PyObject *self, PyObject *args);
PyObject *ops_setLoadConst(PyObject *self, PyObject *args);
PyObject *ops_rayleighDamping(PyObject *self, PyObject *args);
PyObject *ops_eigenAnalysis(PyObject *self, PyObject *args);

#endif
