#include <UniaxialMaterial.h>
#include <Element.h>
#include <TimeSeries.h>
#include <LoadPattern.h>
#include <MultiSupportPattern.h>
#include <string.h>
#include <map>
#include <vector>
#include <cstdlib>
#include <LinearSOE.h>
#include <DOF_Numberer.h>
#include <RCM.h>
#include <AMDNumberer.h>
#include <ConstraintHandler.h>
#include <EquiSolnAlgo.h>
#include <Integrator.h>
#include <ConvergenceTest.h>
#include <elementAPI.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <BeamIntegration.h>
#include <PathSeries.h>
#include <PathTimeSeries.h>
#include <InterpolatedGroundMotion.h>

struct char_cmp { 
  bool operator () (const char *a,const char *b) const 
  {
    return strcmp(a,b)<0;
  } 
};

typedef std::map<const char *, void *(*)(void), char_cmp> OPS_ParsingFunctionMap;



////////// OPS_ParseUniaxialMaterialCommand /////
static OPS_ParsingFunctionMap uniaxialMaterialsMap;

extern void *OPS_ElasticMaterial(void);
extern void *OPS_ElasticPPMaterial(void);
extern void *OPS_ParallelMaterial(void);
extern void *OPS_SeriesMaterial(void);
extern void *OPS_HystereticMaterial(void);
extern void *OPS_CableMaterial(void);
extern void *OPS_Bilin(void);
extern void *OPS_Bilin02(void);
extern void *OPS_NewSteel01(void);
extern void *OPS_NewSteel02(void);
extern void *OPS_RambergOsgoodSteel(void);
extern void *OPS_NewConcrete01(void);
extern void *OPS_NewConcrete02(void);
extern void *OPS_PinchingLimitStateMaterial(void);
extern void *OPS_NewSAWSMaterial(void);
extern void *OPS_NewConcreteZ01Material(void);
extern void *OPS_NewConcreteL01Material(void);
extern void *OPS_NewSteelZ01Material(void);
extern void *OPS_NewTendonL01Material(void);
extern void *OPS_NewConfinedConcrete01Material(void);
extern void *OPS_NewElasticBilin(void);
extern void *OPS_NewMinMaxMaterial(void);
extern void *OPS_SimpleFractureMaterial(void);
extern void *OPS_HoehlerStanton(void);
extern void *OPS_NewInitStrainMaterial(void);
extern void *OPS_NewInitStressMaterial(void);
extern void *OPS_New_pyUCLA(void);
extern void *OPS_Maxwell(void);
extern void *OPS_ViscousDamper(void);
extern void *OPS_BilinearOilDamper(void);
extern void *OPS_Cast(void);
extern void *OPS_Dodd_Restrepo(void);
extern void *OPS_DoddRestr(void);
extern void *OPS_NewElasticMultiLinear(void);
extern void *OPS_ImpactMaterial(void);
extern void *OPS_SteelBRB(void);
extern void *OPS_New_MultiLinear(void);
extern void *OPS_NewHookGap(void);
extern void *OPS_FRPConfinedConcrete(void);
extern void *OPS_NewSteel01Thermal(void);
extern void *OPS_NewSteel02Thermal(void);
extern void *OPS_NewConcrete02Thermal(void);
extern void *OPS_BWBN(void);
extern void *OPS_ModIMKPeakOriented(void);
extern void *OPS_ModIMKPeakOriented02(void);
extern void *OPS_ModIMKPinching(void);
extern void *OPS_ModIMKPinching02(void);
extern void *OPS_ConcretewBeta(void);
extern void *OPS_ConcreteD(void);
extern void *OPS_PinchingLimitState(void);
extern void *OPS_NewOriginCentered(void);
extern void *OPS_NewSteel2(void);
extern void *OPS_ConcreteSakaiKawashima(void);
extern void *OPS_ResilienceMaterialHR(void);
extern void *OPS_CFSSSWP(void);
extern void *OPS_CFSWSWP(void);
extern void *OPS_ResilienceLow(void);
extern void *OPS_ViscousMaterial(void);
extern void *OPS_SteelMPF(void);
extern void *OPS_ConcreteCM(void);

int OPS_SetUpUniaxialMaterials(void) {
  uniaxialMaterialsMap.insert(std::make_pair("Elastic", &OPS_ElasticMaterial));
  uniaxialMaterialsMap.insert(std::make_pair("ElasticPP", &OPS_ElasticPPMaterial));
  uniaxialMaterialsMap.insert(std::make_pair("Parallel", &OPS_ParallelMaterial));
  uniaxialMaterialsMap.insert(std::make_pair("Series", &OPS_SeriesMaterial));
  uniaxialMaterialsMap.insert(std::make_pair("Hysteretic", &OPS_HystereticMaterial));
  uniaxialMaterialsMap.insert(std::make_pair("Cable", &OPS_CableMaterial));
  uniaxialMaterialsMap.insert(std::make_pair("Concrete01", &OPS_NewConcrete01));
  uniaxialMaterialsMap.insert(std::make_pair("Steel01", &OPS_NewSteel01));
  
  return 0;
}

UniaxialMaterial *
OPS_ParseUniaxialMaterialCommand(const char *matType) {
  static bool initDone = false;

  if (initDone == false) {
    OPS_SetUpUniaxialMaterials();
    initDone = true;
  }

  UniaxialMaterial  *theMaterial = 0;
  OPS_ParsingFunctionMap::const_iterator iter = uniaxialMaterialsMap.find(matType);
  if (iter == uniaxialMaterialsMap.end()) {
    return 0;
  }

  void *theMat = (*iter->second)();
  if (theMat != 0) 
    theMaterial = (UniaxialMaterial *)theMat;
  else 
    return 0;

  return theMaterial;

  /* OLD CODE .. IF ELSE
  if (strcmp(matType,"Elastic") == 0) {
    
    void *theMat = OPS_ElasticMaterial();
    if (theMat != 0) 
      theMaterial = (UniaxialMaterial *)theMat;
    else 
      return 0;
    
  } else if (strcmp(matType,"ElasticPP") == 0) {
    
    void *theMat = OPS_ElasticPPMaterial();
    if (theMat != 0) 
      theMaterial = (UniaxialMaterial *)theMat;
    else 
      return 0;
    
  } else if (strcmp(matType,"Parallel") == 0) {
    void *theMat = OPS_ParallelMaterial();
    if (theMat != 0) 
      theMaterial = (UniaxialMaterial *)theMat;
    else 
      return 0;
  }
  return theMaterial;
  */
}


/// OPS_ParseElementCommand
static OPS_ParsingFunctionMap elementsMap;

extern void *OPS_NewTrussElement(void);
extern void *OPS_NewZeroLengthSection(void);
extern void* OPS_NewForceBeamColumn2d();
extern void* OPS_NewForceBeamColumn3d();
extern void* OPS_NewElasticBeam3d();
extern void* OPS_NewElasticBeam2d();

void* OPS_NewForceBeamColumn()
{
    int ndm = OPS_GetNDM();
    if(ndm == 2) {
	return OPS_NewForceBeamColumn2d();
    } else {
	return OPS_NewForceBeamColumn3d();
    }
}

void* OPS_NewElasticBeam()
{
    int ndm = OPS_GetNDM();
    if(ndm == 2) {
	return OPS_NewElasticBeam2d();
    } else {
	return OPS_NewElasticBeam3d();
    }
}


int OPS_SetUpElements(void) {
    elementsMap.insert(std::make_pair("Truss", &OPS_NewTrussElement));
    elementsMap.insert(std::make_pair("truss", &OPS_NewTrussElement));
    elementsMap.insert(std::make_pair("zeroLengthSection", &OPS_NewZeroLengthSection));
    elementsMap.insert(std::make_pair("forceBeamColumn", &OPS_NewForceBeamColumn));
    elementsMap.insert(std::make_pair("elasticBeamColumn", &OPS_NewElasticBeam));
    return 0;
}

Element *
OPS_ParseElementCommand(const char *eleType)
{

    static bool initDone = false;

    if(initDone == false) {
	OPS_SetUpElements();
	initDone = true;
    }

    Element* theEle = 0;
    OPS_ParsingFunctionMap::const_iterator iter = elementsMap.find(eleType);
    if(iter == elementsMap.end()) {
	return 0;
    }

    void *theE = (*iter->second)();
    if(theE != 0) 
	theEle = (Element *)theE;
    else 
	return 0;

    return theEle;

}




/////// OPS_ParseTimeSeriesCommand ///////////////
static OPS_ParsingFunctionMap tsMap;

extern void *OPS_NewConstantSeries(void);
extern void *OPS_NewLinearSeries(void);
extern void *OPS_NewTriangleSeries(void);
extern void *OPS_NewTrigSeries(void);
extern void *OPS_NewRectangularSeries(void);
extern void *OPS_NewPulseSeries(void);

void* OPS_NewPathSeries()
{
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"insufficient arguments: PathSeries\n";
	return 0;
    }

    // get tag
    int tag =0;
    int numData = 1;
    if(OPS_GetIntInput(&numData,&tag) < 0) return 0;

    // get other inputs
    double factor = 1.0, dt = 1.0;
    std::vector<double> values, times;
    const char* timefile = 0, *valfile=0;
    bool isval = false;
    bool istime = false;

    while(OPS_GetNumRemainingInputArgs() > 0) {
	std::string type;
	if(istime==false && isval==false) {
	    type = OPS_GetString();
	}
	if(type == "-time") {
	    istime = true;
	    isval = false;
	} else if(type == "-factor") {
	    istime = false;
	    isval = false;
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&factor) < 0) return 0;
	    }
	} else if(type == "-dt") {
	    istime = false;
	    isval = false;
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&dt) < 0) return 0;
	    }
	} else if(type == "-fileTime") {
	    istime = false;
	    isval = false;
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		timefile = OPS_GetString();
	    }
	} else if(type == "-filePath") {
	    istime = false;
	    isval = false;
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		valfile = OPS_GetString();
	    }
	} else if(type == "-values") {
	    istime = false;
	    isval = true;
	} else if(isval) {
	    double val;
	    if(OPS_GetDoubleInput(&numData,&val) < 0) {
		isval = false;
		istime = false;
		continue;
	    }
	    values.push_back(val);
	} else if(istime) {
	    double val;
	    if(OPS_GetDoubleInput(&numData,&val) < 0) {
		isval = false;
		istime = false;
		continue;
	    }
	    times.push_back(val);
	}
    }

    // create path series
    if(timefile!=0 && valfile!=0 && factor>0) {
	return new PathTimeSeries(tag,timefile,valfile,factor);
    } else if(valfile!=0 && factor>0) {
	return new PathSeries(tag,valfile,dt,factor);
    } else if(values.size()>0 && times.size()>0 && factor>0) {
	Vector valvec(values.size());
	for(int i=0; i<valvec.Size(); i++) {
	    valvec(i) = values[i];
	}
	Vector timevec(times.size());
	for(int i=0; i<timevec.Size(); i++) {
	    timevec(i) = times[i];
	}
	return new PathTimeSeries(tag,valvec,timevec,factor);
    } else if(values.size()>0 && factor>0) {
	Vector valvec(values.size());
	for(int i=0; i<valvec.Size(); i++) {
	    valvec(i) = values[i];
	}
	return new PathSeries(tag,valvec,dt,factor);
    } else {
	opserr<<"choice of options for PathTimeSeries is invalid\n";
	return 0;
    }
}

int OPS_SetUpTimeSeries(void) {
    tsMap.insert(std::make_pair("Constant", &OPS_NewConstantSeries));
    tsMap.insert(std::make_pair("Trig", &OPS_NewTrigSeries));
    tsMap.insert(std::make_pair("Sine", &OPS_NewTrigSeries));
    tsMap.insert(std::make_pair("Linear", &OPS_NewLinearSeries));
    tsMap.insert(std::make_pair("Rectangular", &OPS_NewRectangularSeries));
    tsMap.insert(std::make_pair("Pulse", &OPS_NewPulseSeries));
    tsMap.insert(std::make_pair("Triangle", &OPS_NewTriangleSeries));
    tsMap.insert(std::make_pair("Path", &OPS_NewPathSeries));
    return 0;
}

TimeSeries *
OPS_ParseTimeSeriesCommand(const char *tsType)
{
    static bool initDone = false;

    if(initDone == false) {
	OPS_SetUpTimeSeries();
	initDone = true;
    }

    TimeSeries* theTimeSeries = 0;
    OPS_ParsingFunctionMap::const_iterator iter = tsMap.find(tsType);
    if(iter == tsMap.end()) {
	return 0;
    }

    void *theTs = (*iter->second)();
    if(theTs != 0) 
	theTimeSeries = (TimeSeries *)theTs;
    else 
	return 0;

  return theTimeSeries;

}


/////////////// OPS_ParseLoadPatternCommand //////////////
static OPS_ParsingFunctionMap patternMap;

extern void *OPS_NewLoadPattern(void);
extern void *OPS_NewUniformExcitationPattern(void);
extern void *OPS_NewMultiSupportPattern(void);

int OPS_SetUpLoadpattern(void) {
    patternMap.insert(std::make_pair("Plain", &OPS_NewLoadPattern));
    patternMap.insert(std::make_pair("UniformExcitation", &OPS_NewUniformExcitationPattern));
    patternMap.insert(std::make_pair("MultipleSupport", &OPS_NewMultiSupportPattern));
    return 0;
}

LoadPattern *
OPS_ParseLoadPatternCommand(const char *pType)
{
    static bool initDone = false;

    if(initDone == false) {
	OPS_SetUpLoadpattern();
	initDone = true;
    }

    LoadPattern* theLoadpattern = 0;
    OPS_ParsingFunctionMap::const_iterator iter = patternMap.find(pType);
    if(iter == patternMap.end()) {
	return 0;
    }

    void *theLP = (*iter->second)();
    if(theLP != 0) 
	theLoadpattern = (LoadPattern *)theLP;
    else 
	return 0;

  return theLoadpattern;

}

/////////////// OPS_ParsedSOECommand //////////////
static OPS_ParsingFunctionMap soeMap;

extern void *OPS_NewBandGenLinLapack(void);
extern void *OPS_NewBandSPDLinLapack(void);
extern void *OPS_NewSuperLUSolver();

int OPS_SetUpSOE(void) {
    soeMap.insert(std::make_pair("BandGeneral", &OPS_NewBandGenLinLapack));
    soeMap.insert(std::make_pair("BandGEN", &OPS_NewBandGenLinLapack));
    soeMap.insert(std::make_pair("BandSPD", &OPS_NewBandSPDLinLapack));
    soeMap.insert(std::make_pair("SuperLU", &OPS_NewSuperLUSolver));
    soeMap.insert(std::make_pair("SparseGeneral", &OPS_NewSuperLUSolver));
    soeMap.insert(std::make_pair("SparseGEN", &OPS_NewSuperLUSolver));
    return 0;
}

LinearSOE *
OPS_ParseSOECommand(const char *type)
{
    static bool initDone = false;

    if(initDone == false) {
	OPS_SetUpSOE();
	initDone = true;
    }

    LinearSOE* theSOE = 0;
    OPS_ParsingFunctionMap::const_iterator iter = soeMap.find(type);
    if(iter == soeMap.end()) {
	return 0;
    }

    void *soe = (*iter->second)();
    if(soe != 0) 
	theSOE = (LinearSOE *)soe;
    else 
	return 0;

  return theSOE;

}

/////////////// OPS_ParsedNumbererCommand //////////////
static OPS_ParsingFunctionMap numbererMap;

extern void *OPS_NewPlainNumberer(void);
void *OPS_NewRCMNumberer() {
    RCM* theRCM = new RCM(false);
    DOF_Numberer *theNumberer = new DOF_Numberer(*theRCM);
    return theNumberer;
}
void *OPS_NewAMDNumberer() {
    AMD* theAMD = new AMD();
    DOF_Numberer *theNumberer = new DOF_Numberer(*theAMD);
    return theNumberer;
}

int OPS_SetUpNumberer(void) {
    numbererMap.insert(std::make_pair("Plain", &OPS_NewPlainNumberer));
    numbererMap.insert(std::make_pair("RCM", &OPS_NewRCMNumberer));
    numbererMap.insert(std::make_pair("AMD", &OPS_NewAMDNumberer));
    return 0;
}

DOF_Numberer *
OPS_ParseNumbererCommand(const char *type)
{
    static bool initDone = false;

    if(initDone == false) {
	OPS_SetUpNumberer();
	initDone = true;
    }

    DOF_Numberer *theNumberer = 0;
    OPS_ParsingFunctionMap::const_iterator iter = numbererMap.find(type);
    if(iter == numbererMap.end()) {
	return 0;
    }

    void *numberer = (*iter->second)();
    if(numberer != 0) 
	theNumberer = (DOF_Numberer *)numberer;
    else 
	return 0;

  return theNumberer;

}

/////////////// OPS_ParsedConstraintHandlerCommand //////////////
static OPS_ParsingFunctionMap handlerMap;

extern void *OPS_NewPlainHandler(void);
extern void *OPS_NewPenaltyConstraintHandler(void);
extern void *OPS_NewLagrangeConstraintHandler(void);
extern void *OPS_NewTransformationConstraintHandler(void);


int OPS_SetUpHandler(void) {
    handlerMap.insert(std::make_pair("Plain", &OPS_NewPlainHandler));
    handlerMap.insert(std::make_pair("Penalty", &OPS_NewPenaltyConstraintHandler));
    handlerMap.insert(std::make_pair("Lagrange", &OPS_NewLagrangeConstraintHandler));
    handlerMap.insert(std::make_pair("Transformation", &OPS_NewTransformationConstraintHandler));
    return 0;
}

ConstraintHandler *
OPS_ParseConstraintHandlerCommand(const char *type)
{
    static bool initDone = false;

    if(initDone == false) {
	OPS_SetUpHandler();
	initDone = true;
    }

    ConstraintHandler *theHandler = 0;
    OPS_ParsingFunctionMap::const_iterator iter = handlerMap.find(type);
    if(iter == handlerMap.end()) {
	return 0;
    }

    void *handler = (*iter->second)();
    if(handler != 0) 
	theHandler = (ConstraintHandler *)handler;
    else 
	return 0;

  return theHandler;

}

/////////////// OPS_ParsedAlgorithmCommand //////////////
static OPS_ParsingFunctionMap algoMap;

extern void *OPS_NewLinearAlgorithm();
extern void *OPS_NewNewtonRaphsonAlgorithm();


int OPS_SetUpAlgorithm(void) {
    algoMap.insert(std::make_pair("Linear", &OPS_NewLinearAlgorithm));
    algoMap.insert(std::make_pair("Newton", &OPS_NewNewtonRaphsonAlgorithm));
    return 0;
}

EquiSolnAlgo *
OPS_ParseAlgorithmCommand(const char *type)
{
    static bool initDone = false;

    if(initDone == false) {
	OPS_SetUpAlgorithm();
	initDone = true;
    }

    EquiSolnAlgo *theAlgorithm = 0;
    OPS_ParsingFunctionMap::const_iterator iter = algoMap.find(type);
    if(iter == algoMap.end()) {
	return 0;
    }

    void *algo = (*iter->second)();
    if(algo != 0) 
	theAlgorithm = (EquiSolnAlgo *)algo;
    else 
	return 0;

  return theAlgorithm;

}

/////////////// OPS_ParseIntegratorCommand //////////////
typedef std::map<const char *, std::pair<void *(*)(),int>, char_cmp>
OPS_ParsingIntegratorMap;

static OPS_ParsingIntegratorMap integMap;

extern void* OPS_NewLoadControlIntegrator();
extern void* OPS_NewDisplacementControlIntegrator();
extern void* OPS_NewNewmark();


int OPS_SetUpIntegrator(void) {
    integMap.insert(std::make_pair("LoadControl", std::make_pair(&OPS_NewLoadControlIntegrator,1)));
    integMap.insert(std::make_pair("DisplacementControl",
				   std::make_pair(&OPS_NewDisplacementControlIntegrator,1)));
    integMap.insert(std::make_pair("Newmark",std::make_pair(&OPS_NewNewmark,2)));
    return 0;
}

Integrator *
OPS_ParseIntegratorCommand(const char *type, int& isstatic)
{
    static bool initDone = false;

    if(initDone == false) {
	OPS_SetUpIntegrator();
	initDone = true;
    }

    Integrator *theIntegrator = 0;
    OPS_ParsingIntegratorMap::const_iterator iter = integMap.find(type);
    if(iter == integMap.end()) {
	return 0;
    }

    isstatic = iter->second.second;
    void *integ = (*iter->second.first)();
    if(integ != 0) 
	theIntegrator = (Integrator *)integ;
    else 
	return 0;

  return theIntegrator;

}

/////////////// OPS_ParseCTestCommand //////////////
static OPS_ParsingFunctionMap testMap;

extern void *OPS_NewCTestNormUnbalance(void);
extern void *OPS_NewCTestNormDispIncr(void);

int OPS_SetUpCTest(void) {
    testMap.insert(std::make_pair("NormUnbalance", &OPS_NewCTestNormUnbalance));
    testMap.insert(std::make_pair("NormDispIncr", &OPS_NewCTestNormDispIncr));
    return 0;
}

ConvergenceTest *
OPS_ParseCTestCommand(const char *type)
{
    static bool initDone = false;

    if(initDone == false) {
	OPS_SetUpCTest();
	initDone = true;
    }

    ConvergenceTest* theTest = 0;
    OPS_ParsingFunctionMap::const_iterator iter = testMap.find(type);
    if(iter == testMap.end()) {
	return 0;
    }

    void *test = (*iter->second)();
    if(test != 0) 
	theTest = (ConvergenceTest *)test;
    else 
	return 0;

  return theTest;

}

/////// OPS_ParseSectionCommand ///////////////
static OPS_ParsingFunctionMap secMap;

extern void* OPS_NewElasticSection2d();
extern void* OPS_NewElasticShearSection2d();
extern void* OPS_NewElasticSection3d();
extern void* OPS_NewElasticShearSection3d();
extern void* OPS_NewFiberSection2d();
extern void* OPS_NewFiberSection3d();
extern void* OPS_NewNDFiberSection2d();
extern void* OPS_NewNDFiberSection3d();

void *OPS_NewElasticSection(void)
{
    int numData = OPS_GetNumRemainingInputArgs();
    void* theSec = 0;
    int ndm = OPS_GetNDM();
    if(ndm == 2) {
	if(numData == 3) {
	    theSec = OPS_NewElasticSection2d();
	} else if(numData >=5) {
	    theSec = OPS_NewElasticShearSection2d();
	}
    } else if(ndm == 3) {
	if(numData == 6) {
	    theSec = OPS_NewElasticSection3d();
	} else if(numData >= 8) {
	    theSec = OPS_NewElasticShearSection3d();
	}
    }

    return theSec;
}

void* OPS_NewFiberSection()
{
    void* theSec = 0;
    int ndm = OPS_GetNDM();
    if(ndm == 2) {
	theSec = OPS_NewFiberSection2d();
    } else if(ndm == 3) {
	theSec = OPS_NewFiberSection3d();
    }

    return theSec;
}

void* OPS_NewNDFiberSection()
{
    void* theSec = 0;
    int ndm = OPS_GetNDM();
    if(ndm == 2) {
	theSec = OPS_NewNDFiberSection2d();
    } else if(ndm == 3) {
	theSec = OPS_NewNDFiberSection3d();
    }

    return theSec;
}

int OPS_SetUpSection(void) {
    secMap.insert(std::make_pair("Elastic", &OPS_NewElasticSection));
    secMap.insert(std::make_pair("Fiber", &OPS_NewFiberSection));
    secMap.insert(std::make_pair("NDFiber", &OPS_NewNDFiberSection));
    return 0;
}

SectionForceDeformation *
OPS_ParseSectionCommand(const char *type)
{
    static bool initDone = false;

    if(initDone == false) {
	OPS_SetUpSection();
	initDone = true;
    }

    SectionForceDeformation* theSection = 0;
    OPS_ParsingFunctionMap::const_iterator iter = secMap.find(type);
    if(iter == secMap.end()) {
	return 0;
    }

    void *theSec = (*iter->second)();
    if(theSec != 0) 
	theSection = (SectionForceDeformation *)theSec;
    else 
	return 0;

  return theSection;

}

////////// OPS_ParseNDMaterialCommand /////
static OPS_ParsingFunctionMap nDMaterialMap;

extern  void *OPS_NewReinforcedConcretePlaneStressMaterial(void);
extern  void *OPS_NewFAReinforcedConcretePlaneStressMaterial(void);
extern  void *OPS_NewFAFourSteelRCPlaneStressMaterial(void);
extern  void *OPS_NewRAFourSteelRCPlaneStressMaterial(void);
extern  void *OPS_NewPrestressedConcretePlaneStressMaterial(void);
extern  void *OPS_NewFAPrestressedConcretePlaneStressMaterial(void);
extern  void *OPS_NewFAFourSteelPCPlaneStressMaterial(void);
extern  void *OPS_NewRAFourSteelPCPlaneStressMaterial(void);
extern  void *OPS_NewMaterialCMM(void);
extern  void *OPS_NewElasticIsotropicMaterial(void);
extern  void *OPS_NewElasticOrthotropicMaterial(void);
extern  void *OPS_NewDruckerPragerMaterial(void);
extern  void *OPS_NewBoundingCamClayMaterial(void);
extern  void *OPS_NewContactMaterial2DMaterial(void);
extern  void *OPS_NewContactMaterial3DMaterial(void);
extern  void *OPS_NewInitialStateAnalysisWrapperMaterial(void);
extern  void *OPS_NewManzariDafaliasMaterial(void);
extern  void *OPS_NewManzariDafaliasMaterialRO(void);
extern  void *OPS_CycLiqCPMaterial(void);
extern  void *OPS_CycLiqCPSPMaterial(void);
extern  void *OPS_NewInitStressNDMaterial(void);
extern  void *OPS_NewStressDilatancyMaterial(void);

int OPS_SetUpNDMaterial(void) {
    nDMaterialMap.insert(std::make_pair("ElasticIsotropic3D", &OPS_NewElasticIsotropicMaterial));
    nDMaterialMap.insert(std::make_pair("ElasticIsotropic", &OPS_NewElasticIsotropicMaterial));
    nDMaterialMap.insert(std::make_pair("ElasticOrthotropic3D", &OPS_NewElasticOrthotropicMaterial));
    nDMaterialMap.insert(std::make_pair("ElasticOrthotropic", &OPS_NewElasticOrthotropicMaterial));
    return 0;
}

NDMaterial *
OPS_ParseNDMaterialCommand(const char *type) {
    static bool initDone = false;

    if (initDone == false) {
	OPS_SetUpNDMaterial();
	initDone = true;
    }

    NDMaterial  *theMaterial = 0;
    OPS_ParsingFunctionMap::const_iterator iter = nDMaterialMap.find(type);
    if(iter == nDMaterialMap.end()) {
	return 0;
    }

    void *theMat = (*iter->second)();
    if(theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
    else 
	return 0;

    return theMaterial;

}

////////// OPS_ParseCrdTransfCommand /////
static OPS_ParsingFunctionMap transfMap;

extern void* OPS_NewLinearCrdTransf2d();
extern void* OPS_NewLinearCrdTransf3d();
extern void* OPS_NewPDeltaCrdTransf2d();
extern void* OPS_NewPDeltaCrdTransf3d();
extern void* OPS_NewCorotCrdTransf2d();
extern void* OPS_NewCorotCrdTransf3d();

void* OPS_NewLinearCrdTransf()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if(ndm == 2 && ndf == 3) {
	return OPS_NewLinearCrdTransf2d();
    } else if(ndm == 3 && ndf == 6) {
	return OPS_NewLinearCrdTransf3d();
    } else {
	opserr<<"current NDM and NDF is incompatible with frame elements\n";
	return 0;
    }
}
void* OPS_NewPDeltaCrdTransf()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if(ndm == 2 && ndf == 3) {
	return OPS_NewPDeltaCrdTransf2d();
    } else if(ndm == 3 && ndf == 6) {
	return OPS_NewPDeltaCrdTransf3d();
    } else {
	opserr<<"current NDM and NDF is incompatible with frame elements\n";
	return 0;
    }
}
void* OPS_NewCorotCrdTransf()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if(ndm == 2 && ndf == 3) {
	return OPS_NewCorotCrdTransf2d();
    } else if(ndm == 3 && ndf == 6) {
	return OPS_NewCorotCrdTransf3d();
    } else {
	opserr<<"current NDM and NDF is incompatible with frame elements\n";
	return 0;
    }
}

int OPS_SetUpCrdsTransf(void) {
    transfMap.insert(std::make_pair("Linear", &OPS_NewLinearCrdTransf));
    transfMap.insert(std::make_pair("PDelta", &OPS_NewPDeltaCrdTransf));
    transfMap.insert(std::make_pair("Corotational", &OPS_NewCorotCrdTransf));
    return 0;
}

CrdTransf *
OPS_ParseCrdTransfCommand(const char *type) {
    static bool initDone = false;

    if (initDone == false) {
	OPS_SetUpCrdsTransf();
	initDone = true;
    }

    CrdTransf  *theTransf = 0;
    OPS_ParsingFunctionMap::const_iterator iter = transfMap.find(type);
    if(iter == transfMap.end()) {
	return 0;
    }

    void *theTran = (*iter->second)();
    if(theTran != 0) 
	theTransf = (CrdTransf *)theTran;
    else 
	return 0;

    return theTransf;

}

////////// OPS_ParseBeamIntegrationRuleCommand /////
typedef std::map<const char *, void *(*)(int&,ID&), char_cmp> OPS_BeamIntegraionRuleMap;
static OPS_BeamIntegraionRuleMap ruleMap;

extern void* OPS_NewLobattoBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewLegendreBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewNewtonCotesBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewRadauBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewTrapezoidalBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewCompositeSimpsonBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewUserDefinedBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewFixedLocationBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewLowOrderBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewMidDistanceBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewUserHingeBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewHingeMidpointBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewHingeRadauBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewHingeRadauTwoBeamIntegration(int& integrationTag, ID& secTags);
extern void* OPS_NewHingeEndpointBeamIntegration(int& integrationTag, ID& secTags);

int OPS_SetUpBeamIntegrationRule(void) {
    ruleMap.insert(std::make_pair("Lobatto", &OPS_NewLobattoBeamIntegration));
    ruleMap.insert(std::make_pair("Legendre", &OPS_NewLegendreBeamIntegration));
    ruleMap.insert(std::make_pair("NewtoCotes", &OPS_NewNewtonCotesBeamIntegration));
    ruleMap.insert(std::make_pair("Radau", &OPS_NewRadauBeamIntegration));
    ruleMap.insert(std::make_pair("Trapezoidal", &OPS_NewTrapezoidalBeamIntegration));
    ruleMap.insert(std::make_pair("CompositeSimpson", &OPS_NewCompositeSimpsonBeamIntegration));
    ruleMap.insert(std::make_pair("UserDefined", &OPS_NewUserDefinedBeamIntegration));
    ruleMap.insert(std::make_pair("FixedLocation", &OPS_NewFixedLocationBeamIntegration));
    ruleMap.insert(std::make_pair("LowOrder", &OPS_NewLowOrderBeamIntegration));
    ruleMap.insert(std::make_pair("MidDistance", &OPS_NewMidDistanceBeamIntegration));
    ruleMap.insert(std::make_pair("UserHinge", &OPS_NewUserHingeBeamIntegration));
    ruleMap.insert(std::make_pair("HingeMidpoint", &OPS_NewHingeMidpointBeamIntegration));
    ruleMap.insert(std::make_pair("HingeRadau", &OPS_NewHingeRadauBeamIntegration));
    ruleMap.insert(std::make_pair("HingeRadauTwo", &OPS_NewHingeRadauTwoBeamIntegration));
    ruleMap.insert(std::make_pair("HingeEndpoint", &OPS_NewHingeEndpointBeamIntegration));
    return 0;
}

BeamIntegrationRule *
OPS_ParseBeamIntegrationRuleCommand(const char *type) {

    static bool initDone = false;

    if (initDone == false) {
    	OPS_SetUpBeamIntegrationRule();
    	initDone = true;
    }

    OPS_BeamIntegraionRuleMap::const_iterator iter = ruleMap.find(type);
    if(iter == ruleMap.end()) {
    	return 0;
    }

    int iTag;
    ID secTags;
    void *rule = (*iter->second)(iTag,secTags);
    BeamIntegration* bi = 0;
    if(rule != 0) 
    	bi = (BeamIntegration *)rule;
    else 
    	return 0;

    BeamIntegrationRule* theRule = new BeamIntegrationRule(iTag,bi,secTags);

    return theRule;

}


////////// OPS_ParseGroundMotionCommand /////
extern void* OPS_NewGroundMotion();
void* OPS_NewInterpolatedGroundMotion(MultiSupportPattern& thePattern)
{
    int numArgs = OPS_GetNumRemainingInputArgs();
    int numMotions = (numArgs-1)/2;
    GroundMotion* motions[numMotions];

    // get ground motion tags
    int gmTags[numMotions];
    if(OPS_GetIntInput(&numMotions,&gmTags[0]) < 0) return 0;

    // get factors
    Vector facts(numMotions);
    double* fact_ptr = &facts(0);
    std::string type = OPS_GetString();
    if(type == "-fact") {
	if(OPS_GetDoubleInput(&numMotions,fact_ptr) < 0) return 0;
    }

    // get ground motions
    for(int i=0; i<numMotions; i++) {
	motions[i] = thePattern.getMotion(gmTags[i]);
	if(motions[i] == 0) {
	    opserr<<"ground motion "<<gmTags[i]<<" is not found\n";
	    return 0;
	}
    }

    return new InterpolatedGroundMotion(motions,facts,false);
}


GroundMotion *
OPS_ParseGroundMotionCommand(MultiSupportPattern& thePattern, int&gmTag)
{
    // get tag
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"ERROR must give GroundMotion tag.\n";
	return 0;
    }
    int numData = 1;
    if(OPS_GetIntInput(&numData,&gmTag) < 0) return NULL;

    // get type
    std::string gmType = OPS_GetString();

    // create object
    GroundMotion* theObj = 0;
    if(gmType == "Plain") {
	theObj = (GroundMotion*) OPS_NewGroundMotion();
    } else if(gmType == "Interpolated") {
	theObj = (GroundMotion*) OPS_NewInterpolatedGroundMotion(thePattern);
    }

    if(theObj == 0) return 0;
    
    if(thePattern.addMotion(*theObj,gmTag) == false) {
	delete theObj; // invoke the destructor, otherwise mem leak
	return 0;
    }

    return theObj;
}

