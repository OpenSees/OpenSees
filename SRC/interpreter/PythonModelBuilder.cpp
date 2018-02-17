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
// Description: This file contains the class definition for PythonModelBuilder.
// A PythonModelBuilder is a python module that adds the commands to create
// the model for the standard

#include "PythonModelBuilder.h"
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Vector.h>
#include <Node.h>
#include <Domain.h>
#include <Matrix.h>
#include <iostream>
#include <vector>
#include <SP_Constraint.h>
#include <Element.h>
#include <TimeSeries.h>
#include <string>
#include <LoadPattern.h>
#include <MultiSupportPattern.h>
#include <NodalLoad.h>
#include <SectionForceDeformation.h>
#include <NDMaterial.h>
#include <Fiber.h>
#include <FiberSection2d.h>
#include <FiberSection3d.h>
#include <NDFiberSection2d.h>
#include <NDFiberSection3d.h>
#include <Patch.h>
#include <Cell.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <UniaxialFiber2d.h>
#include <UniaxialFiber3d.h>
#include <NDFiber2d.h>
#include <NDFiber3d.h>
#include <ReinfLayer.h>
#include <ReinfBar.h>
#include <CrdTransf.h>
#include <BeamIntegration.h>
#include <GroundMotion.h>
#include <Beam2dUniformLoad.h>
#include <Beam3dUniformLoad.h>
#include <Beam2dPointLoad.h>
#include <Beam3dPointLoad.h>

static PythonModelBuilder pBuilder;

PythonModelBuilder& OPS_GetPythonModelBuilder()
{
    return pBuilder;
}

extern "C" void OPS_ResetCommandLine(int nArgs, int cArg, PyObject *args);
extern Element* OPS_ParseElementCommand(const char *eleType);
extern TimeSeries* OPS_ParseTimeSeriesCommand(const char *tsType);
extern LoadPattern* OPS_ParseLoadPatternCommand(const char *pType);
extern SectionForceDeformation* OPS_ParseSectionCommand(const char *pType);
extern NDMaterial* OPS_ParseNDMaterialCommand(const char *pType);
extern CrdTransf* OPS_ParseCrdTransfCommand(const char *pType);
extern BeamIntegrationRule* OPS_ParseBeamIntegrationRuleCommand(const char *pType);
extern GroundMotion* OPS_ParseGroundMotionCommand(MultiSupportPattern& thePattern,int& gmTag);

extern void* OPS_Node();
extern void* OPS_NodalLoad();
extern void* OPS_SP();
extern void* OPS_NDFiber3d();
extern void* OPS_NDFiber2d();
extern void* OPS_UniaxialFiber3d();
extern void* OPS_UniaxialFiber2d();
extern void* OPS_QuadPatch();
extern void* OPS_RectPatch();
extern void* OPS_CircPatch();
extern void* OPS_CircReinfLayer();
extern void* OPS_StraightReinfLayer();

PythonModelBuilder::PythonModelBuilder()
    :ndm(0), ndf(0), isset(false), currentLoadPattern(0),
     currentMultiPattern(0), currentSection(0)
{
}

PythonModelBuilder::~PythonModelBuilder()
{
}

PyObject* ops_Model(PyObject* self, PyObject* args)
{
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);
    int ndm = 0, ndf = 0;

    // model builder type
    if(OPS_GetNumRemainingInputArgs() == 0) {
	PyErr_SetString(PyExc_RuntimeError,"model type must be set");
	return NULL;
    }
    std::string model = OPS_GetString();

    // ndm ndf
    while(OPS_GetNumRemainingInputArgs() > 0) {
	std::string type = OPS_GetString();
	int numData = 1;
	if(OPS_GetNumRemainingInputArgs() == 0) break;
	if(type=="-ndm" || type=="-NDM") {
	    OPS_GetIntInput(&numData, &ndm);
	} else if(type=="-ndf" || type=="-NDF") {
	    OPS_GetIntInput(&numData, &ndf);
	}
    }

    if(ndm == 0) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR ndm must be set to positive value");
	return NULL;
    }

    if(ndf == 0) {
	if (ndm == 1) 
	    ndf = 1;
	else if (ndm == 2)
	    ndf = 3;
	else if (ndm == 3)
	    ndf = 6;
	else {
	    PyErr_SetString(PyExc_RuntimeError,"ERROR ndm has wrong value");
	    return NULL;
	}
    }

    // basic model
    if(model=="basic" || model=="Basic" || model=="BasicBuilder" || model=="basicBuilder") {
	PythonModelBuilder& pBuilder = OPS_GetPythonModelBuilder();
	pBuilder.set(ndm,ndf);
    }

    Py_INCREF(Py_None);
    return Py_None;
}


PyObject *ops_addNode(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }

    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // get node
    Node* theNode = (Node*) OPS_Node();
    if(theNode == 0) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR failed to create node");
	return NULL;
    }

    // add the node to the domain
    Domain *theDomain = OPS_GetDomain();
    if(theDomain->addNode(theNode) == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR failed to add node to domain");
	delete theNode;
	return NULL;
    }

    return Py_BuildValue("i", theNode->getTag());
}

PyObject *ops_addHomogeneousBC(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }

    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);
    int numberArgs = OPS_GetNumRemainingInputArgs();

    // get data
    ID data(numberArgs);
    if(OPS_GetIntInput(&numberArgs, &data[0]) < 0) {
	return NULL;
    }

    // get node
    Domain* theDomain = OPS_GetDomain();
    Node* theNode = theDomain->getNode(data[0]);
    if(theNode == 0) {
	opserr<<"ERROR node "<<data[0]<<"does not exsit\n";
	return NULL;
    }
    int ndf = theNode->getNumberDOF();

    // check inputs
    if(numberArgs < ndf+1) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR incorrect number of nodal fix terms");
	return NULL;
    }

    // create homogeneous constraint
    for(int i=0; i<ndf; i++) {
	if(data[i+1] != 0) {
	    SP_Constraint *theSP = new SP_Constraint(data[0], i, 0.0, true);
	    if(theSP == 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR ran out of memory for SP_Constraint");
		return NULL;
	    }
	    if(theDomain->addSP_Constraint(theSP) == false) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR failed to add SP_Constraint to domain");
		delete theSP;
		return NULL;
	    }
	}
    }

    // return
    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *ops_addElement(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }
    
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *eleType = OPS_GetString();
    Element* theEle = OPS_ParseElementCommand(eleType);

    if(theEle == 0) {
	PyErr_SetString(PyExc_RuntimeError, "ERROR could not create element.");
	return NULL;
    }
  
    // Now add the element to the modelBuilder
    Domain* theDomain = OPS_GetDomain();
    if(theDomain->addElement(theEle) == false) {
	PyErr_SetString(PyExc_RuntimeError, "ERROR could not add element to domain.");
	delete theEle; // invoke the material objects destructor, otherwise mem leak
	return NULL;
    }

    return Py_BuildValue("d", theEle->getTag());
}

PyObject *ops_addTimeSeries(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }
    
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *tsType = OPS_GetString();
    TimeSeries* theTimeSeries = OPS_ParseTimeSeriesCommand(tsType);

    if(theTimeSeries == 0) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not create TimeSeries.");
    	return NULL;
    }
  
    // Now add the time
    if(OPS_addTimeSeries(theTimeSeries) == false) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not add TimeSeries.");
    	delete theTimeSeries; // invoke the material objects destructor, otherwise mem leak
    	return NULL;
    }

    return Py_BuildValue("d", theTimeSeries->getTag());
}

PyObject *ops_addPattern(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }

    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *pType = OPS_GetString();
     LoadPattern* thePattern = OPS_ParseLoadPatternCommand(pType);

    if(thePattern == 0) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not create LoadPattern.");
    	return NULL;
    }
  
    // Now add the patter
    Domain& theDomain = *(OPS_GetDomain());
    if(theDomain.addLoadPattern(thePattern) == false) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not add LoadPattern.");
    	delete thePattern; // invoke the material objects destructor, otherwise mem leak
    	return NULL;
    }

    // set current pattern
    std::string type = pType;
    if(type == "Plain") {
	pBuilder.setCurrentLoadPattern(thePattern);
    } else if(type == "MultipleSupport") {
	pBuilder.setCurrentMultiPattern(dynamic_cast<MultiSupportPattern*>(thePattern));
    }

    return Py_BuildValue("d", thePattern->getTag());

}

PyObject *ops_addNodalLoad(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }

    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);
    
    NodalLoad *theLoad = (NodalLoad*) OPS_NodalLoad();
    if(theLoad == 0) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not create NodalLoad.");
    	return NULL;
    }

    // check user pattern
    bool userPattern = false;
    int loadPatternTag = 0;
    if(OPS_GetNumRemainingInputArgs() > 1) {
	std::string type = OPS_GetString();
	if(type == "-pattern") {
	    int numData = 1;
	    if(OPS_GetIntInput(&numData, &loadPatternTag) < 0) {
		delete theLoad;
		return NULL;
	    }
	    userPattern = true;
	}
    }

    // get the current pattern tag
    if(userPattern == false) {
	if(pBuilder.getCurrentLoadPattern() == 0) {
	    PyErr_SetString(PyExc_RuntimeError, "ERROR no current load pattern.");
	    delete theLoad;
	    return NULL;
	}
	loadPatternTag = pBuilder.getCurrentLoadPattern()->getTag();
    }
    
    // Now add the load
    Domain& theDomain = *(OPS_GetDomain());
    if(theDomain.addNodalLoad(theLoad,loadPatternTag) == false) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not add NodalLoad.");
    	delete theLoad; // invoke the material objects destructor, otherwise mem leak
    	return NULL;
    }

    return Py_BuildValue("d", loadPatternTag);
}

PyObject *ops_addSP(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }

    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);
    
    SP_Constraint *theSP = (SP_Constraint*) OPS_SP();
    if(theSP == 0) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not create SP_Constraint.");
    	return NULL;
    }

    // check user pattern
    bool userPattern = false;
    int loadPatternTag = 0;
    if(OPS_GetNumRemainingInputArgs() > 1) {
	std::string type = OPS_GetString();
	if(type == "-pattern") {
	    int numData = 1;
	    if(OPS_GetIntInput(&numData, &loadPatternTag) < 0) {
		delete theSP;
		return NULL;
	    }
	    userPattern = true;
	}
    }

    // get the current pattern tag
    if(userPattern == false) {
	if(pBuilder.getCurrentLoadPattern() == 0) {
	    PyErr_SetString(PyExc_RuntimeError, "ERROR no current load pattern.");
	    delete theSP;
	    return NULL;
	}
	loadPatternTag = pBuilder.getCurrentLoadPattern()->getTag();
    }
    
    // Now add the load
    Domain& theDomain = *(OPS_GetDomain());
    if(theDomain.addSP_Constraint(theSP,loadPatternTag) == false) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not add SP_Constraint.");
    	delete theSP; // invoke the material objects destructor, otherwise mem leak
    	return NULL;
    }

    return Py_BuildValue("d", theSP->getTag());
}

PyObject *ops_addSection(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }
    
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *secType = OPS_GetString();
    SectionForceDeformation* theSec = OPS_ParseSectionCommand(secType);

    if(theSec == 0) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not create Section.");
    	return NULL;
    }
  
    // Now add the section
    if(OPS_addSectionForceDeformation(theSec) == false) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not add Section.");
    	delete theSec; // invoke the material objects destructor, otherwise mem leak
    	return NULL;
    }

    // set current section
    pBuilder.setCurrentSection(theSec);
    
    return Py_BuildValue("d", theSec->getTag());
}

PyObject *ops_addNDMaterial(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }
    
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *matType = OPS_GetString();
    NDMaterial* theMat = OPS_ParseNDMaterialCommand(matType);

    if(theMat == 0) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not create NDMaterial.");
    	return NULL;
    }
  
    if(OPS_addNDMaterial(theMat) == false) {
    	PyErr_SetString(PyExc_RuntimeError, "ERROR could not add NDMaterial.");
    	delete theMat; // invoke the material objects destructor, otherwise mem leak
    	return NULL;
    }

    return Py_BuildValue("d", theMat->getTag());
}

static int ops_addFiberToSection(SectionForceDeformation* section, Fiber* theFiber)
{
    if(section==0 || theFiber==0) return 0;
    if(section->getClassTag() == SEC_TAG_FiberSection2d) {
	FiberSection2d* sec = dynamic_cast<FiberSection2d*>(section);
	if(sec != 0) {
	    if(sec->addFiber(*theFiber) < 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR failed to add Fiber to section");
		delete theFiber;
		return -1;
	    }
	} else {
	    PyErr_SetString(PyExc_RuntimeError,"ERROR failed to get Fiber section");
	    delete theFiber;
	    return -1;
	}
    } else if(section->getClassTag() == SEC_TAG_FiberSection3d) {
	FiberSection3d* sec = dynamic_cast<FiberSection3d*>(section);
	if(sec != 0) {
	    if(sec->addFiber(*theFiber) < 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR failed to add Fiber to section");
		delete theFiber;
		return -1;
	    } else {
		PyErr_SetString(PyExc_RuntimeError,"ERROR failed to get Fiber section");
		delete theFiber;
		return -1;
	    }
	}

    } else if(section->getClassTag() == SEC_TAG_NDFiberSection2d) {
	NDFiberSection2d* sec = dynamic_cast<NDFiberSection2d*>(section);
	if(sec != 0) {
	    if(sec->addFiber(*theFiber) < 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR failed to add Fiber to section");
		delete theFiber;
		return -1;
	    } else {
		PyErr_SetString(PyExc_RuntimeError,"ERROR failed to get Fiber section");
		delete theFiber;
		return -1;
	    }
	}

    } else if(section->getClassTag() == SEC_TAG_NDFiberSection3d) {
	NDFiberSection3d* sec = dynamic_cast<NDFiberSection3d*>(section);
	if(sec != 0) {
	    if(sec->addFiber(*theFiber) < 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR failed to add Fiber to section");
		delete theFiber;
		return -1;
	    } else {
		PyErr_SetString(PyExc_RuntimeError,"ERROR failed to get Fiber section");
		delete theFiber;
		return -1;
	    }
	}
    }

    return 0;
}

PyObject *ops_addFiber(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }
    
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // check if a section is being proccessed
    SectionForceDeformation* section = pBuilder.getCurrentSection();
    if(section == 0) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR no current section");
	return NULL;
    }

    // create fiber
    Fiber *theFiber =0;
    if(section->getClassTag() == SEC_TAG_FiberSection2d) {
	theFiber = (Fiber*) OPS_UniaxialFiber2d();
    } else if(section->getClassTag() == SEC_TAG_FiberSection3d) {
	theFiber = (Fiber*) OPS_UniaxialFiber3d();
    } else if(section->getClassTag() == SEC_TAG_NDFiberSection2d) {
	theFiber = (Fiber*) OPS_NDFiber2d();
    } else if(section->getClassTag() == SEC_TAG_NDFiberSection3d) {
	theFiber = (Fiber*) OPS_NDFiber3d();
    }

    if(theFiber == 0) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR failed to create Fiber");
	return NULL;
    }

    // add fiber to the section
    if(ops_addFiberToSection(section,theFiber) < 0) return NULL;

    return Py_BuildValue("d", section->getTag());
}

PyObject *ops_addPatch(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }
    
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // check if a section is being proccessed
    SectionForceDeformation* section = pBuilder.getCurrentSection();
    if(section == 0) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR no current section");
	return NULL;
    }

    // create fiber
    Patch *thePatch =0;
    std::string type = OPS_GetString();
    if(type == "quad" || type == "quadrilateral") {
	thePatch = (Patch*) OPS_QuadPatch();
    } else if(type == "rect" || type == "rectangular") {
	thePatch = (Patch*) OPS_RectPatch();
    } else if(type == "circ" || type == "circular") {
	thePatch = (Patch*) OPS_CircPatch();
    } else {
	PyErr_SetString(PyExc_RuntimeError,"ERROR unknow patch type");
	return NULL;
    }

    if(thePatch == 0) {
    	PyErr_SetString(PyExc_RuntimeError,"ERROR failed to create Patch");
    	return NULL;
    }

    // add fibers to the section
    int numCells = thePatch->getNumCells();
    int matTag = thePatch->getMaterialID();
    Cell** cells = thePatch->getCells();
    if(cells == 0) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR out of run to create fibers");
	delete thePatch;
    	return NULL;
    }
    for(int j=0; j<numCells; j++) {
	// get fiber data
	double area = cells[j]->getArea();
	const Vector& cPos = cells[j]->getCentroidPosition();

	// create fibers
	Fiber *theFiber = 0;
	UniaxialMaterial *material = 0;
	NDMaterial *ndmaterial = 0;
	
	if(section->getClassTag() == SEC_TAG_FiberSection2d) {
	    material = OPS_getUniaxialMaterial(matTag);
	    if(material == 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR material cannot be found");
		delete thePatch;
		return NULL;
	    }
	    theFiber = new UniaxialFiber2d(j,*material,area,cPos(0));
	    
	} else if(section->getClassTag() == SEC_TAG_FiberSection3d) {
	    material = OPS_getUniaxialMaterial(matTag);
	    if(material == 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR material cannot be found");
		delete thePatch;
		return NULL;
	    }
	    theFiber = new UniaxialFiber3d(j,*material,area,cPos);

	} else if(section->getClassTag() == SEC_TAG_NDFiberSection2d) {
	    ndmaterial = OPS_getNDMaterial(matTag);
	    if(ndmaterial == 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR material cannot be found");
		delete thePatch;
		return NULL;
	    }
	    theFiber = new NDFiber2d(j,*ndmaterial,area,cPos(0));

	} else if(section->getClassTag() == SEC_TAG_NDFiberSection3d) {
	    ndmaterial = OPS_getNDMaterial(matTag);
	    if(ndmaterial == 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR material cannot be found");
		delete thePatch;
		return NULL;
	    }
	    theFiber = new NDFiber3d(j,*ndmaterial,area,cPos(0),cPos(1));

	}

	// add fiber
	if(ops_addFiberToSection(section,theFiber) < 0) {
	    delete thePatch;
	    return NULL;
	}

	delete cells[j];
    }

    delete [] cells;
    delete thePatch;

    return Py_BuildValue("d", section->getTag());
}

PyObject *ops_addLayer(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }
    
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // check if a section is being proccessed
    SectionForceDeformation* section = pBuilder.getCurrentSection();
    if(section == 0) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR no current section");
	return NULL;
    }

    // create fiber
    ReinfLayer *theLayer =0;
    std::string type = OPS_GetString();
    if(type == "straight") {
	theLayer = (ReinfLayer*) OPS_StraightReinfLayer();
    } else if(type == "circ" || type == "circular") {
	theLayer = (ReinfLayer*) OPS_CircReinfLayer();
    } else {
	PyErr_SetString(PyExc_RuntimeError,"ERROR unknow layer type");
	return NULL;
    }

    if(theLayer == 0) {
    	PyErr_SetString(PyExc_RuntimeError,"ERROR failed to create Layer");
    	return NULL;
    }

    // add fibers to the section
    int numReinfBars = theLayer->getNumReinfBars();
    ReinfBar* reinfBar = theLayer->getReinfBars();
    int matTag = theLayer->getMaterialID();

    if(reinfBar == 0) {
    	PyErr_SetString(PyExc_RuntimeError,"ERROR out of run to create fibers");
	delete theLayer;
    	return NULL;
    }

    for(int j=0; j<numReinfBars; j++) {
	
    	// get fiber data
    	double area = reinfBar[j].getArea();
    	const Vector& cPos = reinfBar[j].getPosition();

    	// create fibers
    	Fiber *theFiber = 0;
    	UniaxialMaterial *material = 0;
    	NDMaterial *ndmaterial = 0;
	
    	if(section->getClassTag() == SEC_TAG_FiberSection2d) {
    	    material = OPS_getUniaxialMaterial(matTag);
	    if(material == 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR material cannot be found");
		delete theLayer;
		return NULL;
	    }
    	    theFiber = new UniaxialFiber2d(j,*material,area,cPos(0));
	    
    	} else if(section->getClassTag() == SEC_TAG_FiberSection3d) {
    	    material = OPS_getUniaxialMaterial(matTag);
	    if(material == 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR material cannot be found");
		delete theLayer;
		return NULL;
	    }
    	    theFiber = new UniaxialFiber3d(j,*material,area,cPos);

    	} else if(section->getClassTag() == SEC_TAG_NDFiberSection2d) {
    	    ndmaterial = OPS_getNDMaterial(matTag);
	    if(ndmaterial == 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR material cannot be found");
		delete theLayer;
		return NULL;
	    }
    	    theFiber = new NDFiber2d(j,*ndmaterial,area,cPos(0));

    	} else if(section->getClassTag() == SEC_TAG_NDFiberSection3d) {
    	    ndmaterial = OPS_getNDMaterial(matTag);
	    if(ndmaterial == 0) {
		PyErr_SetString(PyExc_RuntimeError,"ERROR material cannot be found");
		delete theLayer;
		return NULL;
	    }
    	    theFiber = new NDFiber3d(j,*ndmaterial,area,cPos(0),cPos(1));

    	}

    	// add fiber
    	if(ops_addFiberToSection(section,theFiber) < 0) {
	    delete theLayer;
	    return NULL;
	}

    }

    delete [] reinfBar;
    delete theLayer;

    return Py_BuildValue("d", section->getTag());
}

PyObject *ops_addGeomTransf(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }
    
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *type = OPS_GetString();
    CrdTransf* theTransf = OPS_ParseCrdTransfCommand(type);

    if(theTransf == 0) {
	PyErr_SetString(PyExc_RuntimeError, "ERROR could not create CrdTransf.");
	return NULL;
    }
  
    // Now add the transf
    if(OPS_AddCrdTransf(theTransf) == false) {
	PyErr_SetString(PyExc_RuntimeError, "ERROR could not add CrdTransf.");
	delete theTransf; // invoke the destructor, otherwise mem leak
	return NULL;
    }

    return Py_BuildValue("d", theTransf->getTag());

}

PyObject *ops_addBeamIntegrationRule(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }
    
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    const char *type = OPS_GetString();
    BeamIntegrationRule* theRule = OPS_ParseBeamIntegrationRuleCommand(type);

    if(theRule == 0) {
	PyErr_SetString(PyExc_RuntimeError, "ERROR could not create BeamIntegrationRule.");
	return NULL;
    }
  
    // Now add the 
    if(OPS_addBeamIntegrationRule(theRule) == false) {
	PyErr_SetString(PyExc_RuntimeError, "ERROR could not add BeamIntegrationRule.");
	delete theRule; // invoke the destructor, otherwise mem leak
	return NULL;
    }

    return Py_BuildValue("d", theRule->getTag());

}

PyObject *ops_addNodalMass(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }

    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // check input
    int ndf = OPS_GetNDF();
    if(OPS_GetNumRemainingInputArgs() < ndf+1) {
	PyErr_SetString(PyExc_RuntimeError,"insufficient arguments");
	return NULL;
    }

    // get node tag
    int nodeTag;
    int numData = 1;
    if(OPS_GetIntInput(&numData,&nodeTag) < 0) return NULL;

    // mass terms
    Matrix mass(ndf,ndf);
    Vector massvec(ndf);
    if(OPS_GetDoubleInput(&ndf,&massvec[0]) < 0) return NULL;
    for(int i=0; i<ndf; i++) {
	mass(i,i) = massvec[i];
    }

    // set mass
    Domain *theDomain = OPS_GetDomain();
    if(theDomain->setMass(mass,nodeTag) < 0) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR failed to set mass");
	return NULL;
    }

    return Py_BuildValue("i", nodeTag);
}

PyObject *ops_addGroundMotion(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }
    
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // get pattern
    MultiSupportPattern* thePattern = pBuilder.getCurrentMultiPattern();
    if(thePattern == 0) {
	PyErr_SetString(PyExc_RuntimeError, "ERROR no MultiSupportPattern was defined.");
	return NULL;
    }

    // get ground motion
    int gmTag = 0;
    GroundMotion* theObj = OPS_ParseGroundMotionCommand(*thePattern,gmTag);

    if(theObj == 0) {
	PyErr_SetString(PyExc_RuntimeError, "ERROR could not create GroundModtion.");
	return NULL;
    }
  
    return Py_BuildValue("i", gmTag);

}

PyObject *ops_addElementalLoad(PyObject *self, PyObject *args)
{
    // check model builder
    if(pBuilder.isSet() == false) {
	PyErr_SetString(PyExc_RuntimeError,"ERROR builder has not been set");
	return NULL;
    }
    
    OPS_ResetCommandLine(PyTuple_Size(args), 0, args);

    // get the current pattern tag
    if(pBuilder.getCurrentLoadPattern() == 0) {
	PyErr_SetString(PyExc_RuntimeError, "ERROR no current load pattern.");
	return NULL;
    }
    int loadPatternTag = pBuilder.getCurrentLoadPattern()->getTag();

    // get ele tags
    bool isele = false, istype=false;
    std::vector<int> eleTags;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	std::string type;
	if(isele == false) {
	    type = OPS_GetString();
	}
	
	if(type == "-ele") {
	    isele = true;
	} else if(type == "-range") {
	    if(OPS_GetNumRemainingInputArgs() > 1) {
		int elenodes[2];
		int numData = 2;
		if(OPS_GetIntInput(&numData,&elenodes[0]) < 0) return NULL;
		for(int i=elenodes[0]; i<=elenodes[1]; i++) {
		    eleTags.push_back(i);
		}
	    }
	    isele = false;
	} else if(type=="-type") {
	    isele = false;
	    istype = true;
	    break;
	} else if(isele) {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		int eletag;
		int numData = 1;
		if(OPS_GetIntInput(&numData,&eletag) < 0) {
		    isele = false;
		    continue;
		}
		eleTags.push_back(eletag);
	    }
	}
    }

    // check if '-type' is specified
    if(istype == false) {
	PyErr_SetString(PyExc_RuntimeError,"eleLoad - expecting '-type' option");
	return NULL;
    }

    // create load
    int ndm = OPS_GetNDM();
    ElementalLoad *theLoad = 0;
    static int eleLoadTag = 0;
    Domain& theDomain = *(OPS_GetDomain());
    if(OPS_GetNumRemainingInputArgs() > 0) {
	std::string type = OPS_GetString();
	if(type=="-beamUniform" || type=="beamUniform") {
	    // check input
	    int numData = OPS_GetNumRemainingInputArgs();
	    if(numData < ndm-1) {
		PyErr_SetString(PyExc_RuntimeError,"eleLoad - missing arguments for -beamUniform");
		return NULL;
	    }
	    // get loads
	    if(numData > ndm) numData = ndm;
	    Vector data(ndm);
	    for(int i=0; i<ndm; i++) {
		data[i] = 0.0;
	    }
	    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return NULL;
	    // create loads
	    for(int i=0; i<(int)eleTags.size(); i++) {
		if(ndm == 2) {
		    theLoad = new Beam2dUniformLoad(eleLoadTag,data[0],data[1],eleTags[i]);
		    if(theLoad == 0) {
			PyErr_SetString(PyExc_RuntimeError,"can't create eleLoad");
			return NULL;
		    }
		    if(theDomain.addElementalLoad(theLoad,loadPatternTag) == false) {
			PyErr_SetString(PyExc_RuntimeError, "ERROR could not add eleLoad.");
			delete theLoad;
			return NULL;
		    }
		    eleLoadTag++;
		} else if(ndm == 3) {
		    theLoad = new Beam3dUniformLoad(eleLoadTag,data[0],data[1],data[2],eleTags[i]);
		    if(theLoad == 0) {
			PyErr_SetString(PyExc_RuntimeError,"can't create eleLoad");
			return NULL;
		    }
		    if(theDomain.addElementalLoad(theLoad,loadPatternTag) == false) {
			PyErr_SetString(PyExc_RuntimeError, "ERROR could not add eleLoad.");
			delete theLoad;
			return NULL;
		    }
		    eleLoadTag++;
		}
	    }
	    
	} else if(type=="-beamPoint" || type=="beamPoint") {
	    // check input
	    int numData = OPS_GetNumRemainingInputArgs();
	    if(numData < ndm) {
		PyErr_SetString(PyExc_RuntimeError,"eleLoad - missing arguments for -beamPoint");
		return NULL;
	    }
	    // get loads
	    if(numData > ndm+1) numData = ndm+1;
	    Vector data(ndm+1);
	    for(int i=0; i<ndm+1; i++) {
		data[i] = 0.0;
	    }
	    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return NULL;
	    // create loads
	    for(int i=0; i<(int)eleTags.size(); i++) {
		if(ndm == 2) {
		    if(data[1]<0.0 || data[1]>1.0) {
			PyErr_SetString(PyExc_RuntimeError,"invalid xDivL (valid range [0.0,1.0])");
			return NULL;
		    }
		    theLoad = new Beam2dPointLoad(eleLoadTag,data[0],data[1],eleTags[i],data[2]);
		    if(theLoad == 0) {
			PyErr_SetString(PyExc_RuntimeError,"can't create eleLoad");
			return NULL;
		    }
		    if(theDomain.addElementalLoad(theLoad,loadPatternTag) == false) {
			PyErr_SetString(PyExc_RuntimeError, "ERROR could not add eleLoad.");
			delete theLoad;
			return NULL;
		    }
		    eleLoadTag++;
		} else if(ndm == 3) {
		    if(data[2]<0.0 || data[2]>1.0) {
			PyErr_SetString(PyExc_RuntimeError,"invalid xDivL (valid range [0.0,1.0])");
			return NULL;
		    }
		    theLoad = new Beam3dPointLoad(eleLoadTag,data[0],data[1],data[2],eleTags[i],data[3]);
		    if(theLoad == 0) {
			PyErr_SetString(PyExc_RuntimeError,"can't create eleLoad");
			return NULL;
		    }
		    if(theDomain.addElementalLoad(theLoad,loadPatternTag) == false) {
			PyErr_SetString(PyExc_RuntimeError, "ERROR could not add eleLoad.");
			delete theLoad;
			return NULL;
		    }
		    eleLoadTag++;
		}
	    }
	    
	}
    }

    return Py_BuildValue("i", loadPatternTag);
}
