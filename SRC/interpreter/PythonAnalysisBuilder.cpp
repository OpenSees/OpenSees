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


#include "PythonAnalysisBuilder.h"
#include <elementAPI.h>
#include <Domain.h>
#include <FileStream.h>
#include <string>
#include <ID.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <NodeIter.h>

#include <EquiSolnAlgo.h>
#include <StaticIntegrator.h>
#include <TransientIntegrator.h>
#include <LinearSOE.h>
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>
#include <PFEMAnalysis.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <ConvergenceTest.h>
#include <AnalysisModel.h>
#include <NewtonRaphson.h>
#include <PlainHandler.h>
#include <RCM.h>
#include <LoadControl.h>
#include <ProfileSPDLinSolver.h>
#include <CTestNormUnbalance.h>
#include <ProfileSPDLinDirectSolver.h>
#include <TransformationConstraintHandler.h>
#include <CTestPFEM.h>
#include <PFEMIntegrator.h>
#include <PFEMLinSOE.h>
#include <PFEMSolver.h>
#include <UniaxialMaterial.h>
#include <SectionForceDeformation.h>
#include <SectionRepres.h>
#include <TimeSeries.h>
#include <NDMaterial.h>
#include <LoadPattern.h>
#include <CrdTransf.h>
#include <Damping.h>
#include <BeamIntegration.h>
#include <Newmark.h>
#include <EigenSOE.h>
#include <SymBandEigenSolver.h>
#include <SymBandEigenSOE.h>
#include <FullGenEigenSolver.h>
#include <FullGenEigenSOE.h>
#include <ArpackSOE.h>
#include <iostream>
#include <ProfileSPDLinSOE.h>

PythonAnalysisBuilder::PythonAnalysisBuilder()
:theHandler(0),theNumberer(0),theAnalysisModel(0),theAlgorithm(0),
 theSOE(0),theEigenSOE(0),theStaticIntegrator(0),theTransientIntegrator(0),
 theTest(0),theStaticAnalysis(0),theTransientAnalysis(0),
 theVariableTimeStepTransientAnalysis(0),thePFEMAnalysis(0)
{
}

PythonAnalysisBuilder::~PythonAnalysisBuilder()
{
    this->wipe();
    this->resetAll();
}

void PythonAnalysisBuilder::wipe()
{
    if(theStaticAnalysis != 0) {
	theStaticAnalysis->clearAll();
	delete theStaticAnalysis;
	theStaticAnalysis = 0;
    }
    if(theTransientAnalysis != 0) {
	theTransientAnalysis->clearAll();
	delete theTransientAnalysis;
	theTransientAnalysis = 0;
    }
    theVariableTimeStepTransientAnalysis = 0;
    thePFEMAnalysis = 0;

}

void PythonAnalysisBuilder::resetStatic()
{
    theAlgorithm = 0;
    theStaticIntegrator = 0;
    theSOE = 0;
    theNumberer = 0;
    theHandler = 0;
    theTest = 0;
    theAnalysisModel = 0;
}

void PythonAnalysisBuilder::resetTransient()
{
    theAlgorithm = 0;
    theTransientIntegrator = 0;
    theSOE = 0;
    theNumberer = 0;
    theHandler = 0;
    theTest = 0;
    theAnalysisModel = 0;
}

void PythonAnalysisBuilder::resetAll()
{
    theAlgorithm = 0;
    theStaticIntegrator = 0;
    theTransientIntegrator = 0;
    theSOE = 0;
    theNumberer = 0;
    theHandler = 0;
    theTest = 0;
    theAnalysisModel = 0;
    theEigenSOE = 0;
}

void PythonAnalysisBuilder::set(ConstraintHandler* obj) {
    if(obj == 0) return;
    if(theHandler != 0) {
	opserr<<"The handler can only be set once\n";
	return;
    }
    theHandler = obj;
}

void PythonAnalysisBuilder::set(DOF_Numberer* obj) {
    if(obj == 0) return;
    if(theNumberer != 0) {
	opserr<<"The numberer can only be set once for one analysis\n";
	return;
    }
    theNumberer = obj;
    if(theStaticAnalysis != 0) theStaticAnalysis->setNumberer(*obj);
    if(theTransientAnalysis != 0) theTransientAnalysis->setNumberer(*obj);
}

void PythonAnalysisBuilder::set(EquiSolnAlgo* obj) {
    if(obj == 0) return;
    if(theAlgorithm != 0) {
	opserr<<"The algorithm can only be set once for one analysis\n";
	return;
    }
    theAlgorithm = obj;
    if(theStaticAnalysis != 0) theStaticAnalysis->setAlgorithm(*obj);
    if(theTransientAnalysis != 0) theTransientAnalysis->setAlgorithm(*obj);
}

void PythonAnalysisBuilder::set(LinearSOE* obj) {
    if(obj == 0) return;
    if(theSOE != 0) {
	opserr<<"The SOE can only be set once for one analysis\n";
	return;
    }
    theSOE = obj;
    if(theStaticAnalysis != 0) theStaticAnalysis->setLinearSOE(*obj);
    if(theTransientAnalysis != 0) theTransientAnalysis->setLinearSOE(*obj);
}

void PythonAnalysisBuilder::set(Integrator* obj, int isstatic) {
    if(obj == 0) return;
    if(isstatic == 1) {
	if(theStaticIntegrator != 0) {
	    opserr<<"The static integrator can only be set once for one analysis\n";
	    return;
	}
	theStaticIntegrator = dynamic_cast<StaticIntegrator*>(obj);
	if(theStaticIntegrator != 0) {
	    if(theStaticAnalysis != 0) {
		theStaticAnalysis->setIntegrator(*theStaticIntegrator);
		return;
	    }
	}
    }
    if(isstatic == 2) {
	if(theTransientIntegrator != 0) {
	    opserr<<"The transient integrator can only be set once for one analysis\n";
	    return;
	}
	theTransientIntegrator = dynamic_cast<TransientIntegrator*>(obj);
	if(theTransientIntegrator!=0) {
	    if(theTransientAnalysis != 0) {
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
		return;
	    }
	}
    }
}

void PythonAnalysisBuilder::set(ConvergenceTest* obj) {
    if(obj == 0) return;
    if(theTest != 0) {
	opserr<<"The test can only be set once for one analysis\n";
	return;
    }
    theTest = obj;
    if(theStaticAnalysis != 0) theStaticAnalysis->setConvergenceTest(*obj);
    if(theTransientAnalysis != 0) theTransientAnalysis->setConvergenceTest(*obj);
}

void PythonAnalysisBuilder::newStaticAnalysis()
{
    this->wipe();
    
    if(theAnalysisModel == 0) {
	theAnalysisModel = new AnalysisModel();
    }

    if(theTest == 0) {
	theTest = new CTestNormUnbalance(1.0e-6,25,0);
    }
	
    if(theAlgorithm == 0) {
	opserr << "WARNING analysis Static - no Algorithm yet specified, \n";
	opserr << " NewtonRaphson default will be used\n";	    
	theAlgorithm = new NewtonRaphson(*theTest); 
    }
    if(theHandler == 0) {
	opserr << "WARNING analysis Static - no ConstraintHandler yet specified, \n";
	opserr << " PlainHandler default will be used\n";
	theHandler = new PlainHandler();       
    }
    if(theNumberer == 0) {
	opserr << "WARNING analysis Static - no Numberer specified, \n";
	opserr << " RCM default will be used\n";
	RCM *theRCM = new RCM(false);	
	theNumberer = new DOF_Numberer(*theRCM);    	
    }
    if(theStaticIntegrator == 0) {
	opserr << "WARNING analysis Static - no Integrator specified, \n";
	opserr << " StaticIntegrator default will be used\n";
	theStaticIntegrator = new LoadControl(1, 1, 1, 1);       
    }
    if(theSOE == 0) {
	opserr << "WARNING analysis Static - no LinearSOE specified, \n";
	opserr << " ProfileSPDLinSOE default will be used\n";
	ProfileSPDLinSolver *theSolver;
	theSolver = new ProfileSPDLinDirectSolver();
	theSOE = new ProfileSPDLinSOE(*theSolver);      
    }

    Domain* theDomain = OPS_GetDomain();
    theStaticAnalysis = new StaticAnalysis(*theDomain,*theHandler,*theNumberer,*theAnalysisModel,
					   *theAlgorithm,*theSOE,*theStaticIntegrator,theTest);

    // set eigen SOE
    if(theEigenSOE != 0) {
	theStaticAnalysis->setEigenSOE(*theEigenSOE);
    }

    this->resetStatic();
}

void PythonAnalysisBuilder::newTransientAnalysis()
{
    this->wipe();
    
    if(theAnalysisModel == 0) {
	theAnalysisModel = new AnalysisModel();
    }
    if(theTest == 0) {
	theTest = new CTestNormUnbalance(1.0e-6,25,0);
    }
    if(theAlgorithm == 0) {
	opserr << "WARNING analysis Transient - no Algorithm yet specified, \n";
	opserr << " NewtonRaphson default will be used\n";	    
	    
	theAlgorithm = new NewtonRaphson(*theTest); 
    }
    if(theHandler == 0) {
	opserr << "WARNING analysis Transient dt tFinal - no ConstraintHandler\n";
	opserr << " yet specified, PlainHandler default will be used\n";
	theHandler = new PlainHandler();       
    }
    if(theNumberer == 0) {
	opserr << "WARNING analysis Transient dt tFinal - no Numberer specified, \n";
	opserr << " RCM default will be used\n";
	RCM *theRCM = new RCM(false);	
	theNumberer = new DOF_Numberer(*theRCM);    	
    }
    if(theTransientIntegrator == 0) {
	opserr << "WARNING analysis Transient dt tFinal - no Integrator specified, \n";
	opserr << " Newmark(.5,.25) default will be used\n";
	theTransientIntegrator = new Newmark(0.5,0.25);       
    }
    if(theSOE == 0) {
	opserr << "WARNING analysis Transient dt tFinal - no LinearSOE specified, \n";
	opserr << " ProfileSPDLinSOE default will be used\n";
	ProfileSPDLinSolver *theSolver;
	theSolver = new ProfileSPDLinDirectSolver(); 	
	theSOE = new ProfileSPDLinSOE(*theSolver);      
    }

    Domain* theDomain = OPS_GetDomain();
    theTransientAnalysis=new DirectIntegrationAnalysis(*theDomain,*theHandler,*theNumberer,
						       *theAnalysisModel,*theAlgorithm,
						       *theSOE,*theTransientIntegrator,
						       theTest);

    // set eigen SOE
    if(theEigenSOE != 0) {
	if(theTransientAnalysis != 0) {
	  theTransientAnalysis->setEigenSOE(*theEigenSOE);
	}
    }

    this->resetTransient();
}

void PythonAnalysisBuilder::newPFEMAnalysis(double dtmax, double dtmin, double ratio)
{
    this->wipe();
    
    if(theAnalysisModel == 0) {
	theAnalysisModel = new AnalysisModel();
    }

    if(theTest == 0) {
	theTest = new CTestPFEM(1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,10000,100,1,2);
    }
	
    if(theAlgorithm == 0) {
	opserr << "WARNING analysis PFEM - no Algorithm yet specified, \n";
	opserr << " NewtonRaphson default will be used\n";	    
	theAlgorithm = new NewtonRaphson(*theTest); 
    }
    if(theHandler == 0) {
	opserr << "WARNING analysis PFEM - no ConstraintHandler yet specified, \n";
	opserr << " Transformation default will be used\n";
	theHandler = new TransformationConstraintHandler();
    }
    if(theNumberer == 0) {
	opserr << "WARNING analysis PFEM - no Numberer specified, \n";
	opserr << " RCM default will be used\n";
	RCM *theRCM = new RCM(false);	
	theNumberer = new DOF_Numberer(*theRCM);    	
    }
    if(theTransientIntegrator == 0) {
	opserr << "WARNING analysis PFEM - no Integrator specified, \n";
	opserr << " PFEMIntegrator default will be used\n";
	theTransientIntegrator = new PFEMIntegrator();
    }
    if(theSOE == 0) {
	opserr << "WARNING analysis PFEM - no LinearSOE specified, \n";
	opserr << " PFEMLinSOE default will be used\n";
	PFEMSolver* theSolver = new PFEMSolver();
	theSOE = new PFEMLinSOE(*theSolver);
    }

    Domain* theDomain = OPS_GetDomain();
    thePFEMAnalysis = new PFEMAnalysis(*theDomain,
				       *theHandler,
				       *theNumberer,
				       *theAnalysisModel,
				       *theAlgorithm,
				       *theSOE,
				       *theTransientIntegrator,
				       theTest,dtmax,dtmin,ratio);

    theTransientAnalysis = thePFEMAnalysis;

    // set eigen SOE
    if(theEigenSOE != 0) {
	if(theTransientAnalysis != 0) {
	    theTransientAnalysis->setEigenSOE(*theEigenSOE);
	}
    }

    this->resetTransient();
}

void PythonAnalysisBuilder::newEigenAnalysis(int typeSolver, double shift)
{
    // create a new eigen system and solver
    if(theEigenSOE != 0) {
	if(theEigenSOE->getClassTag() != typeSolver) {
	    //	delete theEigenSOE;
	    theEigenSOE = 0;
	}
    }

    if(theEigenSOE == 0) {
                       
	if(typeSolver == EigenSOE_TAGS_SymBandEigenSOE) {
	    SymBandEigenSolver *theEigenSolver = new SymBandEigenSolver(); 
	    theEigenSOE = new SymBandEigenSOE(*theEigenSolver, *theAnalysisModel); 
	    
	} else if(typeSolver == EigenSOE_TAGS_FullGenEigenSOE) {
	
	    FullGenEigenSolver *theEigenSolver = new FullGenEigenSolver();
	    theEigenSOE = new FullGenEigenSOE(*theEigenSolver, *theAnalysisModel);

	} else {
	    theEigenSOE = new ArpackSOE(shift);    
	}
      
	//
	// set the eigen soe in the system
	//

	if(theStaticAnalysis != 0) {
	    theStaticAnalysis->setEigenSOE(*theEigenSOE);
	} else if(theTransientAnalysis != 0) {
	    theTransientAnalysis->setEigenSOE(*theEigenSOE);
	}
    } // theEigenSOE != 0
}

EquiSolnAlgo* PythonAnalysisBuilder::getAlgorithm()
{
    if(theStaticAnalysis != 0) {
	return theStaticAnalysis->getAlgorithm();
    } else if(theTransientAnalysis != 0) {
	return theTransientAnalysis->getAlgorithm();
    }

    return 0;
}

StaticIntegrator* PythonAnalysisBuilder::getStaticIntegrator() {
    if(theStaticAnalysis != 0) {
	return theStaticAnalysis->getIntegrator();
    }
    return 0;
}

TransientIntegrator* PythonAnalysisBuilder::getTransientIntegrator() {
    if(theTransientAnalysis != 0) {
	return theTransientAnalysis->getIntegrator();
    }
    return 0;
}

ConvergenceTest* PythonAnalysisBuilder::getConvergenceTest()
{
    if(theStaticAnalysis != 0) {
	return theStaticAnalysis->getConvergenceTest();
    } else if(theTransientAnalysis != 0) {
	return theTransientAnalysis->getConvergenceTest();
    }

    return 0;
}

static PythonAnalysisBuilder anaBuilder;

 
extern "C" void OPS_ResetCommandLine(int nArgs, int cArg, PyObject *args);
extern LinearSOE* OPS_ParseSOECommand(const char *type);
extern DOF_Numberer* OPS_ParseNumbererCommand(const char *type);
extern ConstraintHandler* OPS_ParseConstraintHandlerCommand(const char *type);
extern EquiSolnAlgo* OPS_ParseAlgorithmCommand(const char *type);
extern Integrator* OPS_ParseIntegratorCommand(const char *type, int& isstatic);
extern ConvergenceTest* OPS_ParseCTestCommand(const char *type);

static int ops_printEle(OPS_Stream &output, PyObject *args);
static int ops_printNode(OPS_Stream &output, PyObject *args);
static int ops_printAlgo(OPS_Stream &output, PyObject *args);
static int ops_printInteg(OPS_Stream &output, PyObject *args);

PyObject *ops_printModel(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // file stream
    FileStream outputFile;
    OPS_Stream *output = &opserr;

    // domain
    Domain& theDomain = *(OPS_GetDomain());

    // just printModel()
    int numData = OPS_GetNumRemainingInputArgs();
    if(numData == 0) {
	opserr << theDomain;
	Py_INCREF(Py_None);
	return Py_None;
    }

    // printModel()
    int res = 0;
    while(numData > 0) {
	std::string type = OPS_GetString();
	if(type=="-ele" || type=="ele") {
	    res = ops_printEle(*output,args);
	    break;
	} else if(type=="-node" || type=="node") {
	    res = ops_printNode(*output,args);
	    break;
	} else if(type=="-integrator" || type=="integrator") {
	    res = ops_printInteg(*output,args);
	    break;
	} else if(type=="-algorithm" || type=="algorithm") {
	    res = ops_printAlgo(*output,args);
	    break;
	} else {
	    // open file
	    if(type=="-file" || type=="file") {
		type = OPS_GetString();
		numData = OPS_GetNumRemainingInputArgs();
	    }
	    if(outputFile.setFile(type.c_str(), APPEND) != 0) {
		PyErr_SetString(PyExc_RuntimeError,"failed to open file ");
		return NULL;
	    }
	    
	    // just print(filename)
	    if(numData == 0) {
		outputFile << theDomain;
		Py_INCREF(Py_None);
		return Py_None;
	    }
	    output = &outputFile;
	}
    }

    if(res < 0) {
	PyErr_SetString(PyExc_RuntimeError,"failed to print ");
	return NULL;
    }
    
    // close the output file
    outputFile.close();

    Py_INCREF(Py_None);
    return Py_None;
}

int ops_printEle(OPS_Stream &output, PyObject *args)
{
    int flag = 0;

    Domain& theDomain = *(OPS_GetDomain());

    // just print(<filename>, 'ele')
    if(OPS_GetNumRemainingInputArgs() == 0) {
	ElementIter &theElements = theDomain.getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    theElement->Print(output);
	}
	return 0;
    }

    // if 'print <filename> Element flag int <int int ..>' get the flag
    std::string type = OPS_GetString();
    if(type=="flag" || type=="-flag") {
	if(OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING printModel(<filename>, 'ele', <'flag', int> no int specified \n";
	    return -1;
	}
	int numData = 1;
	if(OPS_GetIntInput(&numData,&flag) < 0) return -1;
    } else {
	int numData = OPS_GetNumRemainingInputArgs();
	int numArgs = PyTuple_Size(args);
	OPS_ResetCommandLine(numArgs, numArgs-numData-1, args);
    }

    // now print the Elements with the specified flag, 0 by default
    int numEle = OPS_GetNumRemainingInputArgs();
    if(numEle == 0) {
	ElementIter &theElements = theDomain.getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    theElement->Print(output, flag);
	}
	return 0;
    } else { 

	// otherwise print out the specified nodes i j k .. with flag
	ID theEle(numEle);
	if(OPS_GetIntInput(&numEle, &theEle(0)) < 0) return -1;
	theDomain.Print(output, 0, &theEle, flag);
    }

    return 0;
}

int ops_printNode(OPS_Stream &output, PyObject *args)
{
    int flag = 0;

    Domain& theDomain = *(OPS_GetDomain());

    // just print(<filename>, 'node')
    if(OPS_GetNumRemainingInputArgs() == 0) {
	NodeIter &theNodes = theDomain.getNodes();
	Node *theNode;
	while((theNode = theNodes()) != 0) {
	    theNode->Print(output);
	}
	return 0;
    }

    // if 'print <filename> node flag int <int int ..>' get the flag
    std::string type = OPS_GetString();
    if(type=="flag" || type=="-flag") {
	if(OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING printModel(<filename>, 'ele', <'flag', int> no int specified \n";
	    return -1;
	}
	int numData = 1;
	if(OPS_GetIntInput(&numData,&flag) < 0) return -1;
    } else {
	int numData = OPS_GetNumRemainingInputArgs();
	int numArgs = PyTuple_Size(args);
	OPS_ResetCommandLine(numArgs, numArgs-numData-1, args);
    }

    // now print the nodes with the specified flag, 0 by default
    int numNode = OPS_GetNumRemainingInputArgs();
    if(numNode == 0) {
	NodeIter &theNodes = theDomain.getNodes();
	Node *theNode;
	while ((theNode = theNodes()) != 0) {
	    theNode->Print(output, flag);
	}
	return 0;
    } else { 

	// otherwise print out the specified nodes i j k .. with flag
	ID theNode(numNode);
	if(OPS_GetIntInput(&numNode, &theNode(0)) < 0) return -1;
	theDomain.Print(output, &theNode, 0, flag);
    }
    
    return 0;
}

int ops_printAlgo(OPS_Stream &output, PyObject *args)
{
    EquiSolnAlgo* theAlgorithm = anaBuilder.getAlgorithm();
    if(theAlgorithm == 0) return 0;

    // if just 'print <filename> algorithm'- no flag
    if(OPS_GetNumRemainingInputArgs() == 0) { 
	theAlgorithm->Print(output);
	return 0;
    }    

    // if 'print <filename> Algorithm flag' get the flag
    int flag;
    int numData = 1;
    if(OPS_GetIntInput(&numData, &flag) < 0) return -1;
    theAlgorithm->Print(output,flag);

    return 0;
}

int ops_printInteg(OPS_Stream &output, PyObject *args)
{
    StaticIntegrator* theStaticIntegrator=anaBuilder.getStaticIntegrator();
    TransientIntegrator* theTransientIntegrator=anaBuilder.getTransientIntegrator();
    
    if(theStaticIntegrator == 0 && theTransientIntegrator==0)
	return 0;

    IncrementalIntegrator *theIntegrator;
    if (theStaticIntegrator != 0)
	theIntegrator = theStaticIntegrator;
    else
	theIntegrator = theTransientIntegrator;

    // if just 'print <filename> algorithm'- no flag
    if(OPS_GetNumRemainingInputArgs() == 0) { 
	theIntegrator->Print(output);
	return 0;
    }    

    // if 'print <filename> Algorithm flag' get the flag
    int flag;
    int numData = 1;
    if(OPS_GetIntInput(&numData, &flag) < 0) return -1;
    theIntegrator->Print(output,flag);

    return 0;
}

PyObject *ops_wipeAnalysis(PyObject *self, PyObject *args)
{
    anaBuilder.wipe();
    anaBuilder.resetAll();

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_wipeModel(PyObject *self, PyObject *args)
{
    anaBuilder.wipe();
    anaBuilder.resetAll();

    // wipe domain
    Domain* theDomain = OPS_GetDomain();
    theDomain->clearAll();

    // wipe uniaxial material
    OPS_clearAllUniaxialMaterial();
    OPS_clearAllNDMaterial();

    // wipe sections
    OPS_clearAllSectionForceDeformation();
    OPS_clearAllSectionRepres();

    // wipe time series
    OPS_clearAllTimeSeries();

    // wipe GeomTransf
    OPS_ClearAllCrdTransf();

    // wipe damping
    OPS_clearAllDamping();

    // wipe BeamIntegration
    OPS_clearAllBeamIntegrationRule();
    
    ops_Dt = 0.0;

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_specifySOE(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *type = OPS_GetString();
    LinearSOE* theSOE = OPS_ParseSOECommand(type);

    if(theSOE != 0) {
	anaBuilder.set(theSOE);
    } else {
	PyErr_SetString(PyExc_RuntimeError,"failed to create SOE ");
	return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_specifyNumberer(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *type = OPS_GetString();
    DOF_Numberer* theNumberer = OPS_ParseNumbererCommand(type);

    if(theNumberer == 0) {
	PyErr_SetString(PyExc_RuntimeError,"failed to create Numberer ");
	return NULL;
    }

    anaBuilder.set(theNumberer);

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_specifyConstraintHandler(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *type = OPS_GetString();
    ConstraintHandler* theHandler = OPS_ParseConstraintHandlerCommand(type);

    if(theHandler == 0) {
	PyErr_SetString(PyExc_RuntimeError,"failed to create ConstraintHandler ");
	return NULL;
    }

    anaBuilder.set(theHandler);

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_specifyAlgorithm(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *type = OPS_GetString();
    EquiSolnAlgo* theAlgorithm = OPS_ParseAlgorithmCommand(type);

    if(theAlgorithm == 0) {
	PyErr_SetString(PyExc_RuntimeError,"failed to create Algorithm ");
	return NULL;
    }

    anaBuilder.set(theAlgorithm);

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_specifyIntegrator(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *type = OPS_GetString();
    int isstatic = -1;
    Integrator* integ = OPS_ParseIntegratorCommand(type,isstatic);

    if(integ == 0) {
	PyErr_SetString(PyExc_RuntimeError,"failed to create Integrator ");
	return NULL;
    }

    anaBuilder.set(integ,isstatic);

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_specifyAnalysis(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);
    
    if(OPS_GetNumRemainingInputArgs() < 1) {
	PyErr_SetString(PyExc_RuntimeError,"must specify analysis name ");
	return NULL;
    }

    // analysis type
    std::string type = OPS_GetString();

    // create analysis
    if(type == "Static") {
	anaBuilder.newStaticAnalysis();
    } else if(type == "Transient") {
	anaBuilder.newTransientAnalysis();
    } else if(type == "PFEM") {
	int numData = OPS_GetNumRemainingInputArgs();
	if(numData < 2) {
	    PyErr_SetString(PyExc_RuntimeError,"must specify dtmax and dtmin for PFEM analysis");
	    return NULL;
	}

	double data[3] = {0,0,0.5};
	if(numData > 3) numData = 3;
	if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;

	anaBuilder.newPFEMAnalysis(data[0],data[1],data[2]);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_specifyCTest(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *type = OPS_GetString();
    ConvergenceTest* theTest = OPS_ParseCTestCommand(type);

    if(theTest == 0) {
	PyErr_SetString(PyExc_RuntimeError,"failed to create ConvergenceTest ");
	return NULL;
    }

    anaBuilder.set(theTest);
    
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_analyzeModel(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    int result = 0;

    StaticAnalysis* theStaticAnalysis = anaBuilder.getStaticAnalysis();
    DirectIntegrationAnalysis* theTransientAnalysis = anaBuilder.getTransientAnalysis();
    PFEMAnalysis* thePFEMAnalysis = anaBuilder.getPFEMAnalysis();
    VariableTimeStepDirectIntegrationAnalysis* theVariableTimeStepTransientAnalysis =
	anaBuilder.getVariableTimeStepDirectIntegrationAnalysis();

    if(theStaticAnalysis != 0) {
	// static analysis
	if(OPS_GetNumRemainingInputArgs() < 1) {
	    PyErr_SetString(PyExc_RuntimeError,"must set numStep");
	    return NULL;
	}
	
	// get step
	int numStep;
	int numData = 1;
	if(OPS_GetIntInput(&numData,&numStep) < 0) return NULL;
	result = theStaticAnalysis->analyze(numStep);

    } else if(thePFEMAnalysis != 0) {
	// PFEM analysis
	result = thePFEMAnalysis->analyze();
	
    } else if(theTransientAnalysis != 0) {
	// transient analysis
	if(OPS_GetNumRemainingInputArgs() < 2) {
	    PyErr_SetString(PyExc_RuntimeError,"must set numStep and dT");
	    return NULL;
	}
	// get step
	int numStep;
	int numData = 1;
	if(OPS_GetIntInput(&numData,&numStep) < 0) return NULL;
	
	// get dT
	double dT;
	numData = 1;
	if(OPS_GetDoubleInput(&numData,&dT) < 0) return NULL;

	// set global time step variable
	ops_Dt = dT;

	// variable transient analysis
	if(OPS_GetNumRemainingInputArgs() > 3) {
	    // get dtmin and dtmax
	    double dtm[2] = {0,0};
	    numData = 2;
	    if(OPS_GetDoubleInput(&numData,&dtm[0]) < 0) return NULL;

	    // get Jd
	    int Jd;
	    numData = 1;
	    if(OPS_GetIntInput(&numData,&Jd) < 0) return NULL;

	    if(theVariableTimeStepTransientAnalysis != 0) {
		result = theVariableTimeStepTransientAnalysis
		    ->analyze(numData, dT, dtm[0], dtm[1], Jd);
	    } else {
		PyErr_SetString(PyExc_RuntimeError,"no variable time step transient \
analysis object constructed");
		return NULL;
	    }
	} else {
	    result = theTransientAnalysis->analyze(numStep, dT);
	}
    } else {
	PyErr_SetString(PyExc_RuntimeError,"no Analysis type has been specified");
	return NULL;
    }

    if(result < 0) {
	PyErr_SetString(PyExc_RuntimeError,"OpenSees > analyze failed,");
	opserr << "returned: " << result << " error flag\n";
	return NULL;
    }

    return Py_BuildValue("d", result);
}

PyObject *ops_nodeDisp(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 1) {
	PyErr_SetString(PyExc_RuntimeError,"WARNING want - nodeDisp nodeTag? <dof?>");
	return NULL;
    }

    // get inputs
    int data[2] = {0,0};
    if(OPS_GetIntInput(&numData,&data[0]) < 0) return NULL;
    data[1]--;

    // get nodal response
    Domain& theDomain = *(OPS_GetDomain());
    const Vector *nodalResponse = theDomain.getNodeResponse(data[0], Disp);
    if(nodalResponse == 0) {
	PyErr_SetString(PyExc_RuntimeError,"failed to get nodal response");
	return NULL;
    }

    // get disp
    int size = nodalResponse->Size();
    if(data[1] >= 0) {
	if(data[1] >= size) {
	    PyErr_SetString(PyExc_RuntimeError,
			    "WARNING nodeDisp nodeTag? dof? - dofTag? too large\n");
	    return NULL;
	}
	double value = (*nodalResponse)(data[1]);
	return Py_BuildValue("d", value);
    }

    // get list
    PyObject* theList = PyList_New(0);
    if(theList == 0) {
	PyErr_SetString(PyExc_RuntimeError,"failed to create disp list");
	return NULL;
    }

    for(int i=0; i<size; i++) {
	if(PyList_Append(theList,Py_BuildValue("d",(*nodalResponse)(i))) < 0) {
	    PyErr_SetString(PyExc_RuntimeError,"failed to create disp list");
	    return NULL;
	}
    }
    
    return theList;
}

PyObject *ops_getLoadFactor(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 1) {
	PyErr_SetString(PyExc_RuntimeError,"want patternTag");
	return NULL;
    }

    // get inputs
    int patternTag;
    numData = 1;
    if(OPS_GetIntInput(&numData,&patternTag) < 0) return NULL;

    // get load pattern
    Domain& theDomain = *(OPS_GetDomain());
    LoadPattern* thePattern = theDomain.getLoadPattern(patternTag);
    if(thePattern == 0) {
	PyErr_SetString(PyExc_RuntimeError,"the load pattern is not found");
	return NULL;
    }

    // get load factor
    double value = thePattern->getLoadFactor();
    return Py_BuildValue("d", value);
}

PyObject *ops_setLoadConst(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 2) {
	Py_INCREF(Py_None);
	return Py_None;
    }

    Domain* theDomain = OPS_GetDomain();
    theDomain->setLoadConstant();

    std::string type = OPS_GetString();
    
    if(type == "-time") {
	double newTime;
	numData = 1;
	if(OPS_GetDoubleInput(&numData,&newTime) < 0) return NULL;
	theDomain->setCurrentTime(newTime);
	theDomain->setCommittedTime(newTime);
    }
    Py_INCREF(Py_None);
    return Py_None;
}


PyObject *ops_rayleighDamping(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // check inputs
    if(OPS_GetNumRemainingInputArgs() < 4) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR rayleigh(alphaM,betaK,betaK0,betaKc)");
	return NULL;
    }

    // get parameters
    double data[4];
    int numData = 4;
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;

    Domain* theDomain = OPS_GetDomain();
    theDomain->setRayleighDampingFactors(data[0],data[1],data[2],data[3]);
    
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_eigenAnalysis(PyObject *self, PyObject *args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // check inputs
    int numArgs = OPS_GetNumRemainingInputArgs();
    if(numArgs < 1) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR eigen('type',numModes)");
	return NULL;
    }

    // get parameters
    bool generalizedAlgo = true;// 0 - frequency/generalized (default),1 - standard, 2 - buckling
    int typeSolver = EigenSOE_TAGS_ArpackSOE;
    double shift = 0.0;
    bool findSmallest = true;
    int numEigen = 0;

    // check type of eigenvalue analysis
    if(numArgs >1) {
	std::string type = OPS_GetString();
	if(type=="frequency"||type=="-frenquency"||type=="generalized"||type=="-generalized") {
	    generalizedAlgo = true;
	} else if(type=="standard"||type=="-standard") {
	    generalizedAlgo = false;
	} else if(type=="-findLargest") {
	    findSmallest = false;
	} else if(type=="genBandArpack"||type=="--genBandArpack"||
		  type=="genBandArpackEigen"||type=="-genBandArpackEigen") {
	    typeSolver = EigenSOE_TAGS_ArpackSOE;
	} else if(type=="symmBandLapack"||type=="-symmBandLapack"||
		  type=="symmBandLapackEigen"||type=="-symmBandLapackEigen") {
	    typeSolver = EigenSOE_TAGS_SymBandEigenSOE;
	} else if(type=="fullGenLapack"||type=="-fullGenLapack"||
		  type=="fullGenLapackEigen"||type=="-fullGenLapackEigen") {
	    typeSolver = EigenSOE_TAGS_FullGenEigenSOE;
	} else {
	    PyErr_SetString(PyExc_RuntimeError,"eigen - unknown option specified");
	    return NULL;
	}
    }

    int numData = 1;
    if(OPS_GetIntInput(&numData,&numEigen) < 0) return NULL;

    anaBuilder.newEigenAnalysis(typeSolver,shift);

    // create a transient analysis if no analysis exists
    Domain* theDomain = OPS_GetDomain();

    StaticAnalysis* theStaticAnalysis = anaBuilder.getStaticAnalysis();
    DirectIntegrationAnalysis* theTransientAnalysis = anaBuilder.getTransientAnalysis();
    if(theStaticAnalysis == 0 && theTransientAnalysis == 0) {
	anaBuilder.newTransientAnalysis();
    }

    // call analysis
    int res = 0;
    if(theStaticAnalysis != 0) {
	res = theStaticAnalysis->eigen(numEigen,generalizedAlgo,findSmallest);
    } else if(theTransientAnalysis != 0) {
	res = theTransientAnalysis->eigen(numEigen,generalizedAlgo,findSmallest);
    }

    // return
    PyObject* theList = PyList_New(0);
    if(res == 0) {
	const Vector& eigenvalues = theDomain->getEigenvalues();
	// get list
	if(theList == 0) {
	    PyErr_SetString(PyExc_RuntimeError,"failed to create disp list");
	    return NULL;
	}

	for(int i=0; i<numEigen; i++) {
	    if(PyList_Append(theList,Py_BuildValue("d",eigenvalues(i))) < 0) {
		PyErr_SetString(PyExc_RuntimeError,"failed to create disp list");
		return NULL;
	    }
	}
    }
    
    return theList;

}
