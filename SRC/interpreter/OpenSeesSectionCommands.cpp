/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
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

// Written: Minjie

// Description: command to create section

#include <SectionForceDeformation.h>
#include <map>
#include <elementAPI.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <ElasticMaterial.h>
#include <SectionAggregator.h>
#include <FiberSection2d.h>
#include <FiberSection3d.h>
#include <FiberSectionWarping3d.h>
#include <FiberSectionAsym3d.h>
#include <NDFiberSection2d.h>
#include <NDFiberSection3d.h>
#include <UniaxialFiber2d.h>
#include <UniaxialFiber3d.h>
#include <NDFiber2d.h>
#include <NDFiber3d.h>
#include <Patch.h>
#include <QuadPatch.h>
#include <CircPatch.h>
#include <Cell.h>
#include <ReinfLayer.h>
#include <ReinfBar.h>
#include <TubeSectionIntegration.h>
#include <WideFlangeSectionIntegration.h>
#include <NDFiberSectionWarping2d.h>
#include <RCTBeamSectionIntegration.h>
#include <RCCircularSectionIntegration.h>
#include <RCTunnelSectionIntegration.h>
#include <ParallelSection.h>
#include <FiberSection2dThermal.h>
#include <FiberSection3dThermal.h>
//#include <FiberSectionGJThermal.h>
#include <MembranePlateFiberSectionThermal.h>
#include <LayeredShellFiberSectionThermal.h>

void* OPS_ElasticSection2d();
void* OPS_ElasticSection3d();
void* OPS_ElasticShearSection2d();
void* OPS_ElasticShearSection3d();
void* OPS_ElasticBDShearSection2d();
void* OPS_FiberSection2d();
void* OPS_FiberSection3d();
void* OPS_FiberSectionWarping3d();
void* OPS_FiberSectionAsym3d();
void* OPS_NDFiberSection2d();
void* OPS_NDFiberSection3d();
void* OPS_NDFiberSectionWarping2d();
void* OPS_UniaxialFiber2d();
void* OPS_UniaxialFiber3d();
void* OPS_NDFiber2d();
void* OPS_NDFiber3d();
void* OPS_CircPatch();
void* OPS_QuadPatch();
void* OPS_StraightReinfLayer();
void* OPS_CircReinfLayer();
void* OPS_RectPatch();
void* OPS_ElasticMembranePlateSection();
void* OPS_MembranePlateFiberSection();
void* OPS_DoubleMembranePlateFiberSection();
void* OPS_ElasticWarpingShearSection2d();
void* OPS_ElasticTubeSection3d();
void* OPS_ParallelSection();
void* OPS_SectionAggregator();
void* OPS_ElasticPlateSection();
void* OPS_LayeredShellFiberSection();
void* OPS_Bidirectional();
void* OPS_Elliptical2();
void* OPS_Isolator2spring();
void* OPS_FiberSection2dThermal();
void* OPS_HSSSection();
void* OPS_TubeSection();
void* OPS_RCSection2d();
void* OPS_WFSection2d();
void* OPS_RCCircularSection();
void* OPS_RCTunnelSection();
void* OPS_UniaxialSection();
void* OPS_RCTBeamSection2d();
void* OPS_RCTBeamSectionUniMat2d();
void* OPS_ReinforcedConcreteLayerMembraneSection01();	// M. J. Nunez - UChile
void* OPS_ReinforcedConcreteLayerMembraneSection02();	// M. J. Nunez - UChile

namespace {
    static FiberSection2d* theActiveFiberSection2d = 0;
    static FiberSection3d* theActiveFiberSection3d = 0;
    static FiberSectionWarping3d* theActiveFiberSectionWarping3d = 0;  
	static FiberSectionAsym3d* theActiveFiberSectionAsym3d = 0;
    static NDFiberSection2d* theActiveNDFiberSection2d = 0;
    static NDFiberSection3d* theActiveNDFiberSection3d = 0;
    static NDFiberSectionWarping2d* theActiveNDFiberSectionWarping2d = 0;    

    static FiberSection2dThermal* theActiveFiberSection2dThermal = 0;
    static FiberSection3dThermal* theActiveFiberSection3dThermal = 0;
  //static FiberSectionGJThermal* theActiveFiberSectionGJThermal = 0;

    static bool initDone = false;

    struct char_cmp {
	bool operator () (const char *a,const char *b) const
	    {
		return strcmp(a,b)<0;
	    }
    };

    typedef std::map<const char *, void *(*)(void), char_cmp> OPS_ParsingFunctionMap;


    static OPS_ParsingFunctionMap functionMap;


    static void *OPS_ElasticSection(void)
    {
	int numData = OPS_GetNumRemainingInputArgs();
	void* theSec = 0;
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    if(numData == 4) {
		theSec = OPS_ElasticSection2d();
	    } else if(numData >=5) {
		theSec = OPS_ElasticShearSection2d();
	    }
	} else if(ndm == 3) {
	    if(numData == 7) {
		theSec = OPS_ElasticSection3d();
	    } else if(numData >= 8) {
		theSec = OPS_ElasticShearSection3d();
	    }
	}

	return theSec;
    }

    static void* OPS_FiberSection3dThermal(bool& isTorsion)
    {
	int numData = OPS_GetNumRemainingInputArgs();
	if(numData < 1) {
	    opserr<<"insufficient arguments for FiberSection3dThermal\n";
	    return 0;
	}

	numData = 1;
	int tag;
	if (OPS_GetIntInput(&numData, &tag) < 0) {
	    opserr<<"WARNING: failed to read tag\n";
	    return 0;
	}

	UniaxialMaterial *torsion = 0;
	const char* opt = OPS_GetString();
	numData = 1;
	bool deleteTorsion = false;
	if (strcmp(opt, "-GJ") == 0) {
	  double GJ;
	  if (OPS_GetDoubleInput(&numData, &GJ) < 0) {
	    opserr << "WARNING: failed to read GJ\n";
	    return 0;
	  }
	  torsion = new ElasticMaterial(0,GJ);
	  deleteTorsion = true;
	}
	if (strcmp(opt, "-torsion") == 0) {
	  int torsionTag;
	  if (OPS_GetIntInput(&numData, &torsionTag) < 0) {
	    opserr << "WARNING: failed to read torsion\n";
	    return 0;
	  }
	  torsion = OPS_getUniaxialMaterial(torsionTag);
	}
	if (torsion == 0) {
	  opserr << "WARNING torsion not speified for FiberSection\n";
	  opserr << "\nFiberSection3dThermal section: " << tag << endln;
	  return 0;
	}

	int num = 30;

	SectionForceDeformation *theSec = new FiberSection3dThermal(tag,num);
	if (deleteTorsion)
	  delete torsion;
	return theSec;
    }

    static void* OPS_FiberSection()
    {
	void* theSec = 0;
	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();	
	if(ndm == 2) {
	    theSec = OPS_FiberSection2d();
	    theActiveFiberSection2d = (FiberSection2d*)theSec;
	} else if(ndm == 3) {
	    theSec = OPS_FiberSection3d();
	    theActiveFiberSection3d = (FiberSection3d*)theSec;
	}

	return theSec;
    }

    static void* OPS_FiberSectionWarping()
    {
	void* theSec = 0;
	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();	
	if(ndm == 2) {
	  //theSec = OPS_FiberSectionWarping2d();
	  //theActiveFiberSectionWarping2d = (FiberSectionWarping2d*)theSec;
	} else if(ndm == 3) {
	  theSec = OPS_FiberSectionWarping3d();
	  theActiveFiberSectionWarping3d = (FiberSectionWarping3d*)theSec;
	}

	return theSec;
    }
  
	static void* OPS_FiberSectionAsym()
	{
		void* theSec = 0;
		int ndm = OPS_GetNDM();
		int ndf = OPS_GetNDF();
		if (ndm == 2) {
			//theSec = OPS_FiberSectionAsym2d();
			//theActiveFiberSectionAsym2d = (FiberSectionAsym2d*)theSec;
		}
		else if (ndm == 3) {
			theSec = OPS_FiberSectionAsym3d();
			theActiveFiberSectionAsym3d = (FiberSectionAsym3d*)theSec;
		}

		return theSec;
	}

    static void* OPS_FiberSectionThermal()
    {
	void* theSec = 0;
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    theSec = OPS_FiberSection2dThermal();
	    theActiveFiberSection2dThermal = (FiberSection2dThermal*)theSec;
	} else if(ndm == 3) {
	    bool isTorsion = false;
	    theSec = OPS_FiberSection3dThermal(isTorsion);
	    theActiveFiberSection3dThermal = (FiberSection3dThermal*)theSec;
	}

	return theSec;
    }

    static void* OPS_NDFiberSection()
    {
	void* theSec = 0;
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	    theSec = OPS_NDFiberSection2d();
	    theActiveNDFiberSection2d = (NDFiberSection2d*)theSec;
	} else if(ndm == 3) {
	    theSec = OPS_NDFiberSection3d();
	    theActiveNDFiberSection3d = (NDFiberSection3d*)theSec;
	}

	return theSec;
    }

    static void* OPS_NDFiberSectionWarping()
    {
	void* theSec = 0;
	int ndm = OPS_GetNDM();
	if(ndm == 2) {
	  theSec = OPS_NDFiberSectionWarping2d();
	  theActiveNDFiberSectionWarping2d = (NDFiberSectionWarping2d*)theSec;
	} else if(ndm == 3) {
	  //theSec = OPS_NDFiberSection3d();
	  //theActiveNDFiberSection3d = (NDFiberSection3d*)theSec;
	}

	return theSec;
    }  

    static int setUpFunctions(void)
    {
	functionMap.insert(std::make_pair("Elastic", &OPS_ElasticSection));
	functionMap.insert(std::make_pair("ElasticBD", &OPS_ElasticBDShearSection2d));
	functionMap.insert(std::make_pair("Fiber", &OPS_FiberSection));
	functionMap.insert(std::make_pair("fiberSec", &OPS_FiberSection));
	functionMap.insert(std::make_pair("FiberWarping", &OPS_FiberSectionWarping));
	functionMap.insert(std::make_pair("FiberAsym", &OPS_FiberSectionAsym));
	functionMap.insert(std::make_pair("FiberThermal", &OPS_FiberSectionThermal));	
	functionMap.insert(std::make_pair("NDFiber", &OPS_NDFiberSection));
	functionMap.insert(std::make_pair("NDFiberWarping", &OPS_NDFiberSectionWarping));	
	functionMap.insert(std::make_pair("Uniaxial", &OPS_UniaxialSection));
	functionMap.insert(std::make_pair("Generic1D", &OPS_UniaxialSection));
	functionMap.insert(std::make_pair("Generic1d", &OPS_UniaxialSection));
	functionMap.insert(std::make_pair("ElasticMembranePlateSection", &OPS_ElasticMembranePlateSection));
	functionMap.insert(std::make_pair("PlateFiber", &OPS_MembranePlateFiberSection));
	functionMap.insert(std::make_pair("DoublePlateFiber", &OPS_DoubleMembranePlateFiberSection));	
	functionMap.insert(std::make_pair("ElasticWarpingShear", &OPS_ElasticWarpingShearSection2d));
	functionMap.insert(std::make_pair("ElasticTube", &OPS_ElasticTubeSection3d));
	functionMap.insert(std::make_pair("Tube", &OPS_TubeSection));
	functionMap.insert(std::make_pair("HSS", &OPS_HSSSection));
	functionMap.insert(std::make_pair("WFSection2d", &OPS_WFSection2d));	
	functionMap.insert(std::make_pair("WSection2d", &OPS_WFSection2d));
	functionMap.insert(std::make_pair("RCSection2d", &OPS_RCSection2d));
	functionMap.insert(std::make_pair("RCTBeamSection2d", &OPS_RCTBeamSection2d));
	functionMap.insert(std::make_pair("RCTBeamSectionUniMat2d", &OPS_RCTBeamSectionUniMat2d));
	functionMap.insert(std::make_pair("Parallel", &OPS_ParallelSection));
	functionMap.insert(std::make_pair("Aggregator", &OPS_SectionAggregator));
	functionMap.insert(std::make_pair("AddDeformation", &OPS_SectionAggregator));
	functionMap.insert(std::make_pair("ElasticPlateSection", &OPS_ElasticPlateSection));
	functionMap.insert(std::make_pair("LayeredShell", &OPS_LayeredShellFiberSection));
	functionMap.insert(std::make_pair("Bidirectional", &OPS_Bidirectional));
	functionMap.insert(std::make_pair("Elliptical", &OPS_Elliptical2));	
	functionMap.insert(std::make_pair("Isolator2spring", &OPS_Isolator2spring));
	functionMap.insert(std::make_pair("RCCircularSection", &OPS_RCCircularSection));
	functionMap.insert(std::make_pair("RCTunnelSection", &OPS_RCTunnelSection));
	functionMap.insert(std::make_pair("ReinforcedConcreteLayerMembraneSection01", &OPS_ReinforcedConcreteLayerMembraneSection01));
	functionMap.insert(std::make_pair("ReinforcedConcreteLayerMembraneSection02", &OPS_ReinforcedConcreteLayerMembraneSection02));

	return 0;
    }
}


int OPS_Section()
{
    theActiveFiberSection2d = 0;
    theActiveFiberSection3d = 0;
    theActiveFiberSectionWarping3d = 0;    
	theActiveFiberSectionAsym3d = 0;
    theActiveNDFiberSection2d = 0;
    theActiveNDFiberSection3d = 0;
    theActiveNDFiberSectionWarping2d = 0;
    
    theActiveFiberSection2dThermal = 0;
    theActiveFiberSection3dThermal = 0;
    //theActiveFiberSectionGJThermal = 0;

    if (initDone == false) {
	setUpFunctions();
	initDone = true;
    }

    // num args
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"WARNING insufficient args: pattern type ...\n";
	return -1;
    }

    const char* type = OPS_GetString();

    OPS_ParsingFunctionMap::const_iterator iter = functionMap.find(type);
    if (iter == functionMap.end()) {
	opserr<<"WARNING section type " << type << " is unknown\n";
	return -1;
    }

    SectionForceDeformation* theSection = (SectionForceDeformation*) (*iter->second)();
    if (theSection == 0) {
	return -1;
    }

    // Now add the section
    if (OPS_addSectionForceDeformation(theSection) == false) {
	opserr<<"ERROR could not add section.\n";
	theActiveFiberSection2d = 0;
	theActiveFiberSection3d = 0;
	theActiveFiberSectionWarping3d = 0;	
	theActiveFiberSectionAsym3d = 0;
	theActiveNDFiberSection2d = 0;
	theActiveNDFiberSection3d = 0;
	theActiveNDFiberSectionWarping2d = 0;
	
	theActiveFiberSection2dThermal = 0;
	theActiveFiberSection3dThermal = 0;
	//theActiveFiberSectionGJThermal = 0;
	delete theSection;
	return -1;
    }

    return 0;
}

int OPS_Fiber()
{
    // create fiber

    Fiber* theFiber = 0;

    if (theActiveFiberSection2d != 0 || theActiveFiberSection2dThermal!=0) {

	theFiber = (UniaxialFiber2d*) OPS_UniaxialFiber2d();

    } else if (theActiveFiberSection3d != 0 || theActiveFiberSectionWarping3d != 0 || theActiveFiberSectionAsym3d != 0 || theActiveFiberSection3dThermal!=0) {

	theFiber = (UniaxialFiber3d*) OPS_UniaxialFiber3d();

    } else if (theActiveNDFiberSection2d != 0 || theActiveNDFiberSectionWarping2d != 0) {

	theFiber = (NDFiber2d*) OPS_NDFiber2d();

    } else if (theActiveNDFiberSection3d != 0) {

	theFiber = (NDFiber3d*) OPS_NDFiber3d();

    }

    if (theFiber == 0) {

	opserr<<"WARNING failed to create fiber\n";
	return -1;

    }

    // add fiber to section
    int res = 0;
    if (theActiveFiberSection2d != 0) {

	res = theActiveFiberSection2d->addFiber(*theFiber);

    } else if (theActiveFiberSection3d != 0) {

	res = theActiveFiberSection3d->addFiber(*theFiber);

    } else if (theActiveFiberSectionWarping3d != 0) {

	res = theActiveFiberSectionWarping3d->addFiber(*theFiber);	

	} else if (theActiveFiberSectionAsym3d != 0) {

	res = theActiveFiberSectionAsym3d->addFiber(*theFiber);

    } else if (theActiveNDFiberSection2d != 0) {

	res = theActiveNDFiberSection2d->addFiber(*theFiber);

    } else if (theActiveNDFiberSection3d != 0) {

	res = theActiveNDFiberSection3d->addFiber(*theFiber);

    } else if (theActiveNDFiberSectionWarping2d != 0) {

	res = theActiveNDFiberSectionWarping2d->addFiber(*theFiber);
	
    } else if (theActiveFiberSection2dThermal != 0) {

	res = theActiveFiberSection2dThermal->addFiber(*theFiber);

    } else if (theActiveFiberSection3dThermal != 0) {

	res = theActiveFiberSection3dThermal->addFiber(*theFiber);

    }

    if (res < 0) {
	opserr << "WARNING failed to add fiber to section\n";
	delete theFiber;
	return -1;
    }

    return 0;
}

int OPS_Patch()
{
    // num args
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"WARNING insufficient args: patch type ...\n";
	return -1;
    }

    // create patch
    Patch *thePatch = 0;
    const char* type = OPS_GetString();
    if(strcmp(type,"quad")==0 || strcmp(type,"quadr")==0 || strcmp(type,"quadrilateral")==0) {
	thePatch = (QuadPatch*) OPS_QuadPatch();
    } else if(strcmp(type,"rect")==0 || strcmp(type,"rectangular")==0) {
	thePatch = (QuadPatch*) OPS_RectPatch();
    } else if(strcmp(type,"circ")==0 || strcmp(type,"circular")==0) {
	thePatch = (CircPatch*) OPS_CircPatch();
    } else {
	opserr<<"ERROR unknown patch type\n";
	return -1;
    }

    if (thePatch == 0) {
	opserr<<"WARNING failed to create patch\n";
	return -1;
    }

    // add fibers to the section
    int numCells = thePatch->getNumCells();
    int matTag = thePatch->getMaterialID();
    Cell** cells = thePatch->getCells();
    if(cells == 0) {
	opserr << "ERROR out of run to create fibers\n";
	delete thePatch;
	return -1;
    }
    for(int j=0; j<numCells; j++) {
	// get fiber data
	double area = cells[j]->getArea();
	const Vector& cPos = cells[j]->getCentroidPosition();

	// create fibers
	Fiber *theFiber = 0;
	UniaxialMaterial *material = 0;
	NDMaterial *ndmaterial = 0;

	if (theActiveFiberSection2d != 0) {

	    material = OPS_getUniaxialMaterial(matTag);
	    if (material == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete thePatch;
		return -1;
	    }
	    theFiber = new UniaxialFiber2d(j,*material,area,cPos(0));
	    theActiveFiberSection2d->addFiber(*theFiber);

	} else if (theActiveFiberSection2dThermal != 0) {

	    material = OPS_getUniaxialMaterial(matTag);
	    if (material == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete thePatch;
		return -1;
	    }
	    theFiber = new UniaxialFiber2d(j,*material,area,cPos(0));
	    theActiveFiberSection2dThermal->addFiber(*theFiber);

	} else if (theActiveFiberSection3d != 0) {

	    material = OPS_getUniaxialMaterial(matTag);
	    if (material == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete thePatch;
		return -1;
	    }
	    theFiber = new UniaxialFiber3d(j,*material,area,cPos);
	    theActiveFiberSection3d->addFiber(*theFiber);

	} else if (theActiveFiberSectionWarping3d != 0) {

	    material = OPS_getUniaxialMaterial(matTag);
	    if (material == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete thePatch;
		return -1;
	    }
	    theFiber = new UniaxialFiber3d(j,*material,area,cPos);
	    theActiveFiberSectionWarping3d->addFiber(*theFiber);	

	} else if (theActiveFiberSectionAsym3d != 0) {

		material = OPS_getUniaxialMaterial(matTag);
		if (material == 0) {
			opserr << "WARNING material " << matTag << " cannot be found\n";
			delete thePatch;
			return -1;
		}
		theFiber = new UniaxialFiber3d(j, *material, area, cPos);
		theActiveFiberSectionAsym3d->addFiber(*theFiber);

	} else if (theActiveFiberSection3dThermal != 0) {

	    material = OPS_getUniaxialMaterial(matTag);
	    if (material == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete thePatch;
		return -1;
	    }
	    theFiber = new UniaxialFiber3d(j,*material,area,cPos);
	    theActiveFiberSection3dThermal->addFiber(*theFiber);

	} else if (theActiveNDFiberSection2d != 0) {

	    ndmaterial = OPS_getNDMaterial(matTag);
	    if (ndmaterial == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete thePatch;
		return -1;
	    }
	    theFiber = new NDFiber2d(j,*ndmaterial,area,cPos(0));
	    theActiveNDFiberSection2d->addFiber(*theFiber);

	} else if (theActiveNDFiberSection3d != 0) {

	    ndmaterial = OPS_getNDMaterial(matTag);
	    if (ndmaterial == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete thePatch;
		return -1;
	    }
	    theFiber = new NDFiber3d(j,*ndmaterial,area,cPos(0),cPos(1));
	    theActiveNDFiberSection3d->addFiber(*theFiber);

	} else if (theActiveNDFiberSectionWarping2d != 0) {

	    ndmaterial = OPS_getNDMaterial(matTag);
	    if (ndmaterial == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete thePatch;
		return -1;
	    }
	    theFiber = new NDFiber2d(j,*ndmaterial,area,cPos(0));
	    theActiveNDFiberSectionWarping2d->addFiber(*theFiber);
	}

	if (theFiber != 0)
	  delete theFiber;
	delete cells[j];
    }

    delete [] cells;
    delete thePatch;
    return 0;
}

int OPS_Layer()
{
    // num args
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"WARNING insufficient args: layer type ...\n";
	return -1;
    }

    // create patch
    ReinfLayer *theLayer = 0;
    const char* type = OPS_GetString();

    if(strcmp(type,"straight")==0) {
	theLayer = (ReinfLayer*) OPS_StraightReinfLayer();
    } else if(strcmp(type,"circ")==0 || strcmp(type,"circular")==0) {
	theLayer = (ReinfLayer*) OPS_CircReinfLayer();
    } else {
	opserr<<"ERROR unknown layer type\n";
	return -1;
    }

    if (theLayer == 0) {
	opserr<<"WARNING failed to create layer\n";
	return -1;
    }

    // add fibers to the section
    int numReinfBars = theLayer->getNumReinfBars();
    ReinfBar* reinfBar = theLayer->getReinfBars();
    int matTag = theLayer->getMaterialID();

    if(reinfBar == 0) {
	opserr<<"ERROR out of run to create fibers\n";
	delete theLayer;
	return -1;
    }

    for(int j=0; j<numReinfBars; j++) {

	// get fiber data
	double area = reinfBar[j].getArea();
	const Vector& cPos = reinfBar[j].getPosition();

	// create fibers
	Fiber *theFiber = 0;
	UniaxialMaterial *material = 0;
	NDMaterial *ndmaterial = 0;

	if (theActiveFiberSection2d != 0) {

	    material = OPS_getUniaxialMaterial(matTag);
	    if (material == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete theLayer;
		return -1;
	    }
	    theFiber = new UniaxialFiber2d(j,*material,area,cPos(0));
	    theActiveFiberSection2d->addFiber(*theFiber);

	} else if (theActiveFiberSection2dThermal != 0) {

	    material = OPS_getUniaxialMaterial(matTag);
	    if (material == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete theLayer;
		return -1;
	    }
	    theFiber = new UniaxialFiber2d(j,*material,area,cPos(0));
	    theActiveFiberSection2dThermal->addFiber(*theFiber);

	} else if (theActiveFiberSection3d != 0) {

	    material = OPS_getUniaxialMaterial(matTag);
	    if (material == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete theLayer;
		return -1;
	    }
	    theFiber = new UniaxialFiber3d(j,*material,area,cPos);
	    theActiveFiberSection3d->addFiber(*theFiber);

	} else if (theActiveFiberSectionWarping3d != 0) {

	    material = OPS_getUniaxialMaterial(matTag);
	    if (material == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete theLayer;
		return -1;
	    }
	    theFiber = new UniaxialFiber3d(j,*material,area,cPos);
	    theActiveFiberSectionWarping3d->addFiber(*theFiber);	    

	} else if (theActiveFiberSectionAsym3d != 0) {

		material = OPS_getUniaxialMaterial(matTag);
		if (material == 0) {
			opserr << "WARNING material " << matTag << " cannot be found\n";
			delete theLayer;
			return -1;
		}
		theFiber = new UniaxialFiber3d(j, *material, area, cPos);
		theActiveFiberSectionAsym3d->addFiber(*theFiber);

	} else if (theActiveFiberSection3dThermal != 0) {

	    material = OPS_getUniaxialMaterial(matTag);
	    if (material == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete theLayer;
		return -1;
	    }
	    theFiber = new UniaxialFiber3d(j,*material,area,cPos);
	    theActiveFiberSection3dThermal->addFiber(*theFiber);

	} else if (theActiveNDFiberSection2d != 0) {

	    ndmaterial = OPS_getNDMaterial(matTag);
	    if (ndmaterial == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete theLayer;
		return -1;
	    }
	    theFiber = new NDFiber2d(j,*ndmaterial,area,cPos(0));
	    theActiveNDFiberSection2d->addFiber(*theFiber);

	} else if (theActiveNDFiberSection3d != 0) {

	    ndmaterial = OPS_getNDMaterial(matTag);
	    if (ndmaterial == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete theLayer;
		return -1;
	    }
	    theFiber = new NDFiber3d(j,*ndmaterial,area,cPos(0),cPos(1));
	    theActiveNDFiberSection3d->addFiber(*theFiber);

	} else if (theActiveNDFiberSectionWarping2d != 0) {

	    ndmaterial = OPS_getNDMaterial(matTag);
	    if (ndmaterial == 0) {
		opserr << "WARNING material "<<matTag<<" cannot be found\n";
		delete theLayer;
		return -1;
	    }
	    theFiber = new NDFiber2d(j,*ndmaterial,area,cPos(0));
	    theActiveNDFiberSectionWarping2d->addFiber(*theFiber);
	}
    
	if (theFiber != 0)
	  delete theFiber;
    }

    delete [] reinfBar;
    delete theLayer;


    return 0;
}
