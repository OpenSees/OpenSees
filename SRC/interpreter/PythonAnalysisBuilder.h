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
