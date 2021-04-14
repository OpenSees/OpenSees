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
#include <RCSectionIntegration.h>
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
void* OPS_FiberSection2d();
void* OPS_FiberSection3d();
void* OPS_FiberSectionWarping3d();
void* OPS_NDFiberSection2d();
void* OPS_NDFiberSection3d();
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
void* OPS_ElasticWarpingShearSection2d();
void* OPS_ElasticTubeSection3d();
void* OPS_ParallelSection();
void* OPS_SectionAggregator();
void* OPS_ElasticPlateSection();
void* OPS_MembranePlateFiberSection();
void* OPS_LayeredShellFiberSection();
void* OPS_Bidirectional();
void* OPS_Elliptical2();
void* OPS_Isolator2spring();
void* OPS_FiberSection2dThermal();

namespace {
    static FiberSection2d* theActiveFiberSection2d = 0;
    static FiberSection3d* theActiveFiberSection3d = 0;
    static FiberSectionWarping3d* theActiveFiberSectionWarping3d = 0;  
    static NDFiberSection2d* theActiveNDFiberSection2d = 0;
    static NDFiberSection3d* theActiveNDFiberSection3d = 0;

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
	  if (ndf == 7) {
	    theSec = OPS_FiberSectionWarping3d();
	    theActiveFiberSectionWarping3d = (FiberSectionWarping3d*)theSec;
	  } else {
	    theSec = OPS_FiberSection3d();
	    theActiveFiberSection3d = (FiberSection3d*)theSec;
	  }
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

    static void* OPS_UniaxialSection()
    {
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 3) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section Uniaxial tag? 1DTag? code?" << endln;
	    return 0;
	}

	int data[2];
	numdata = 2;
	if (OPS_GetIntInput(&numdata, data) < 0) {
	    opserr << "WARNING invalid integer" << endln;
	    return 0;
	}

	int code;
	const char* type = OPS_GetString();
	if (strcmp(type,"Mz") == 0)
	    code = SECTION_RESPONSE_MZ;
	else if (strcmp(type,"P") == 0)
	    code = SECTION_RESPONSE_P;
	else if (strcmp(type,"Vy") == 0)
	    code = SECTION_RESPONSE_VY;
	else if (strcmp(type,"My") == 0)
	    code = SECTION_RESPONSE_MY;
	else if (strcmp(type,"Vz") == 0)
	    code = SECTION_RESPONSE_VZ;
	else if (strcmp(type,"T") == 0)
	    code = SECTION_RESPONSE_T;
	else {
	    opserr << "WARNING invalid code" << endln;
	    opserr << "Uniaxial section: " << data[0] << endln;
	    return 0;
	}

	// Retrieve the uniaxial material from the model builder
	UniaxialMaterial *theMat = OPS_getUniaxialMaterial(data[1]);

	if (theMat == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "uniaxial material: " << data[0];
	    opserr << "\nUniaxial section: " << data[1] << endln;
	    return 0;
	}

	// Parsing was successful, allocate the section
	//theSection = new GenericSection1d (tag, *theMat, code);

	UniaxialMaterial *theMats[1];
	theMats[0] = theMat;
	ID codeID(1);
	codeID(0) = code;
	return new SectionAggregator(data[0], 1, theMats, codeID);
    }
  
  static void* OPS_TubeSection()
  {
    if (OPS_GetNumRemainingInputArgs() < 6) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section Tube tag? matTag? D? t? nfw? nfr? <-nd shape?>" << endln;
	    return 0;
	}

	int tag, matTag;
	double D, t;
	int nfw, nfr;

	SectionForceDeformation* theSection = 0;

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
	    opserr << "WARNING invalid section Tube tag" << endln;
	    return 0;
	}

	if (OPS_GetIntInput(&numdata, &matTag) < 0) {
	    opserr << "WARNING invalid section Tube matTag" << endln;
	    return 0;
	}

	if (OPS_GetDoubleInput(&numdata, &D) < 0) {
	    opserr << "WARNING invalid D" << endln;
	    opserr << "Tube section: " << tag << endln;
	    return 0;
	}

	if (OPS_GetDoubleInput(&numdata, &t) < 0) {
	    opserr << "WARNING invalid t" << endln;
	    opserr << "Tube section: " << tag << endln;
	    return 0;
	}

	if (OPS_GetIntInput(&numdata, &nfw) < 0) {
	    opserr << "WARNING invalid nfw" << endln;
	    opserr << "Tube section: " << tag << endln;
	    return 0;
	}

	if (OPS_GetIntInput(&numdata, &nfr) < 0) {
	    opserr << "WARNING invalid nfr" << endln;
	    opserr << "Tube section: " << tag << endln;
	    return 0;
	}

	TubeSectionIntegration tubesect(D, t, nfw, nfr);

	int numFibers = tubesect.getNumFibers();

	if (OPS_GetNumRemainingInputArgs() > 0) {

	    double shape = 1.0;
	    if (OPS_GetNumRemainingInputArgs() > 1) {
		if (OPS_GetDoubleInput(&numdata, &shape) < 0) {
		    opserr << "WARNING invalid shape" << endln;
		    opserr << "Tube section: " << tag << endln;
		    return 0;
		}
	    }

	    NDMaterial *theSteel = OPS_getNDMaterial(matTag);

	    if (theSteel == 0) {
		opserr << "WARNING ND material does not exist\n";
		opserr << "material: " << matTag;
		opserr << "\nTube section: " << tag << endln;
		return 0;
	    }

	    NDMaterial **theMats = new NDMaterial *[numFibers];

	    tubesect.arrangeFibers(theMats, theSteel);

	    // Parsing was successful, allocate the section
	    theSection = 0;
	    if (OPS_GetNumRemainingInputArgs() > 0) {
		const char* flag = OPS_GetString();
		if (strcmp(flag,"-nd") == 0) {
		    theSection = new NDFiberSection3d(tag, numFibers, theMats, tubesect, shape);
		} else if (strcmp(flag,"-ndWarping") == 0) {
		    theSection = new NDFiberSectionWarping2d(tag, numFibers, theMats, tubesect, shape);
		}
	    }
	    delete [] theMats;
	}
	else {
	    UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(matTag);

	    if (theSteel == 0) {
		opserr << "WARNING uniaxial material does not exist\n";
		opserr << "material: " << matTag;
		opserr << "\nTube section: " << tag << endln;
		return 0;
	    }

	    UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];

	    tubesect.arrangeFibers(theMats, theSteel);

	    // Parsing was successful, allocate the section
	    theSection = new FiberSection2d(tag, numFibers, theMats, tubesect);

	    delete [] theMats;
	}
	
	return theSection;
  }

    static void* OPS_WFSection2d()
    {
	if (OPS_GetNumRemainingInputArgs() < 8) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section WFSection2d tag? matTag? d? tw? bf? tf? nfdw? nftf? <-nd shape?>" << endln;
	    return 0;
	}

	int tag, matTag;
	double d, tw, bf, tf;
	int nfdw, nftf;

	SectionForceDeformation* theSection = 0;

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
	    opserr << "WARNING invalid section WFSection2d tag" << endln;
	    return 0;
	}

	if (OPS_GetIntInput(&numdata, &matTag) < 0) {
	    opserr << "WARNING invalid section WFSection2d matTag" << endln;
	    return 0;
	}

	if (OPS_GetDoubleInput(&numdata, &d) < 0) {
	    opserr << "WARNING invalid d" << endln;
	    opserr << "WFSection2d section: " << tag << endln;
	    return 0;
	}

	if (OPS_GetDoubleInput(&numdata, &tw) < 0) {
	    opserr << "WARNING invalid tw" << endln;
	    opserr << "WFSection2d section: " << tag << endln;
	    return 0;
	}

	if (OPS_GetDoubleInput(&numdata, &bf) < 0) {
	    opserr << "WARNING invalid bf" << endln;
	    opserr << "WFSection2d section: " << tag << endln;
	    return 0;
	}

	if (OPS_GetDoubleInput(&numdata, &tf) < 0) {
	    opserr << "WARNING invalid tf" << endln;
	    opserr << "WFSection2d section: " << tag << endln;
	    return 0;
	}

	if (OPS_GetIntInput(&numdata, &nfdw) < 0) {
	    opserr << "WARNING invalid nfdw" << endln;
	    opserr << "WFSection2d section: " << tag << endln;
	    return 0;
	}

	if (OPS_GetIntInput(&numdata, &nftf) < 0) {
	    opserr << "WARNING invalid nftf" << endln;
	    opserr << "WFSection2d section: " << tag << endln;
	    return 0;
	}

	WideFlangeSectionIntegration wfsect(d, tw, bf, tf, nfdw, nftf);

	int numFibers = wfsect.getNumFibers();

	if (OPS_GetNumRemainingInputArgs() > 0) {

	    double shape = 1.0;
	    if (OPS_GetNumRemainingInputArgs() > 1) {
		if (OPS_GetDoubleInput(&numdata, &shape) < 0) {
		    opserr << "WARNING invalid shape" << endln;
		    opserr << "WFSection2d section: " << tag << endln;
		    return 0;
		}
	    }

	    NDMaterial *theSteel = OPS_getNDMaterial(matTag);

	    if (theSteel == 0) {
		opserr << "WARNING ND material does not exist\n";
		opserr << "material: " << matTag;
		opserr << "\nWFSection2d section: " << tag << endln;
		return 0;
	    }

	    NDMaterial **theMats = new NDMaterial *[numFibers];

	    wfsect.arrangeFibers(theMats, theSteel);

	    // Parsing was successful, allocate the section
	    theSection = 0;
	    if (OPS_GetNumRemainingInputArgs() > 0) {
		const char* flag = OPS_GetString();
		if (strcmp(flag,"-nd") == 0) {
		    theSection = new NDFiberSection2d(tag, numFibers, theMats, wfsect, shape);
		} else if (strcmp(flag,"-ndWarping") == 0) {
		    theSection = new NDFiberSectionWarping2d(tag, numFibers, theMats, wfsect, shape);
		}
	    }
	    delete [] theMats;
	}
	else {
	    UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(matTag);

	    if (theSteel == 0) {
		opserr << "WARNING uniaxial material does not exist\n";
		opserr << "material: " << matTag;
		opserr << "\nWFSection2d section: " << tag << endln;
		return 0;
	    }

	    UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];

	    wfsect.arrangeFibers(theMats, theSteel);

	    // Parsing was successful, allocate the section
	    theSection = new FiberSection2d(tag, numFibers, theMats, wfsect);

	    delete [] theMats;
	}

	return theSection;
    }  

    static void* OPS_RCSection2d()
    {
	if (OPS_GetNumRemainingInputArgs() < 13) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section RCSection2d tag? coreTag? coverTag? steelTag? d? b? cover? Atop? Abottom? Aside? nfcore? nfcover? nfs?" << endln;
	    return 0;
	}

	// int
	int numdata = 4;
	int idata[4];
	if (OPS_GetIntInput(&numdata, idata) < 0) {
	    opserr << "WARNING invalid section RCSection2d int inputs" << endln;
	    return 0;
	}

	int tag = idata[0];
	int coreTag = idata[1];
	int coverTag = idata[2];
	int steelTag = idata[3];

	// double
	numdata = 6;
	double data[6];
	if (OPS_GetDoubleInput(&numdata, data) < 0) {
	    opserr << "WARNING invalid section RCSection2d double inputs" << endln;
	    opserr << "RCSection2d section: " << tag << endln;
	    return 0;
	}

	double d = data[0];
	double b = data[1];
	double cover = data[2];
	double Atop = data[3];
	double Abottom = data[4];
	double Aside = data[5];

	// int
	numdata = 3;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
	    opserr << "WARNING invalid section RCSection2d int inputs" << endln;
	    opserr << "RCSection2d section: " << tag << endln;
	    return 0;
	}
	int nfcore = idata[0];
	int nfcover = idata[1];
	int nfs = idata[2];

	UniaxialMaterial *theCore = OPS_getUniaxialMaterial(coreTag);

	if (theCore == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "material: " << coreTag;
	    opserr << "\nRCSection2d section: " << tag << endln;
	    return 0;
	}

	UniaxialMaterial *theCover = OPS_getUniaxialMaterial(coverTag);

	if (theCover == 0) {
	    opserr << "WARNING uniaxial material does not exist\4n";
	    opserr << "material: " << coverTag;
	    opserr << "\nRCSection2d section: " << tag << endln;
	    return 0;
	}

	UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(steelTag);

	if (theSteel == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "material: " << steelTag;
	    opserr << "\nRCSection2d section: " << tag << endln;
	    return 0;
	}

	RCSectionIntegration rcsect(d, b, Atop, Abottom, Aside, cover, nfcore, nfcover, nfs);

	int numFibers = rcsect.getNumFibers();

	UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];

	rcsect.arrangeFibers(theMats, theCore, theCover, theSteel);

	// Parsing was successful, allocate the section
	SectionForceDeformation* theSection = new FiberSection2d(tag, numFibers, theMats, rcsect);

	delete [] theMats;
	return theSection;
    }

    static void* OPS_RCTBeamSection2d()
    {
	if (OPS_GetNumRemainingInputArgs() < 18) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section RCTBeamSection2d tag? coreTag? coverTag? steelTag? d? bw? beff? hf? Atop? Abottom? flcov? wcov? Nflcover? Nwcover? Nflcore? Nwcore? NsteelTop?  NsteelBottom?" << endln;
	    return 0;
	}

	// int
	int numdata = 4;
	int idata[6];
	if (OPS_GetIntInput(&numdata, idata) < 0) {
	    opserr << "WARNING invalid section RCTBeamSection2d int inputs" << endln;
	    return 0;
	}
	int tag = idata[0];
	int coreTag = idata[1];
	int coverTag = idata[2];
	int steelTag = idata[3];

	// double
	numdata = 8;
	double data[8];
	if (OPS_GetDoubleInput(&numdata, data) < 0) {
	    opserr << "WARNING invalid double inputs" << endln;
	    opserr << "RCTBeamSection2d section: " << tag << endln;
	    return 0;
	}

	double d = data[0];
	double bw = data[1];
	double beff = data[2];
	double hf = data[3];
	double Atop = data[4];
	double Abottom = data[5];
	double flcov = data[6];
	double wcov = data[7];

	// int
	numdata = 6;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
	    opserr << "WARNING invalid section RCTBeamSection2d int inputs" << endln;
	    return 0;
	}
	int Nflcover = idata[0];
	int Nwcover = idata[1];
	int Nflcore = idata[2];
	int Nwcore = idata[3];
	int NsteelTop = idata[4];
	int NsteelBottom = idata[5];

	UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(steelTag);
	if (theSteel == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "material: " << steelTag;
	    opserr << "\nRCTBeamSection2d section: " << tag << endln;
	    return 0;
	}

	RCTBeamSectionIntegration
	    rctbeamsect(d, bw, beff, hf, Atop, Abottom, flcov, wcov,
			Nflcover, Nwcover, Nflcore, Nwcore,
			NsteelTop, NsteelBottom);


	NDMaterial *theCore = OPS_getNDMaterial(coreTag);
	if (theCore == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "material: " << coreTag;
	    opserr << "\nRCTBeamSection2d section: " << tag << endln;
	    return 0;
	}

	NDMaterial *theCover = OPS_getNDMaterial(coverTag);
	if (theCover == 0) {
	    opserr << "WARNING uniaxial material does not exist\4n";
	    opserr << "material: " << coverTag;
	    opserr << "\nRCTBeamSection2d section: " << tag << endln;
	    return 0;
	}

	int numCFibers = rctbeamsect.getNumFibers(concrete);
	int numSFibers = rctbeamsect.getNumFibers(steel);

	NDMaterial **theNDMat = new NDMaterial *[numCFibers];
	UniaxialMaterial **theUniMat = new UniaxialMaterial *[numSFibers];

	rctbeamsect.arrangeFibers(theUniMat, theNDMat, theCore, theCover, theSteel);

	//theSection = new McftSection2dfiber(tag, theNDMat, theUniMat, rctbeamsect);



	RCTBeamSectionIntegration
	    steel(d, bw, beff, hf, Atop, Abottom, flcov, wcov,
		  0, 0, 0, 0,
		  NsteelTop, NsteelBottom);
	steel.arrangeFibers(theUniMat, theNDMat, 0, 0, theSteel);
	FiberSection2d steelSec(0, numSFibers, theUniMat, steel);

	RCTBeamSectionIntegration
	    concrete(d, bw, beff, hf, Atop, Abottom, flcov, wcov,
		     Nflcover, Nwcover, Nflcore, Nwcore,
		     0, 0);
	concrete.arrangeFibers(theUniMat, theNDMat, theCore, theCover, 0);
	NDFiberSection2d concSec(0, numCFibers, theNDMat, concrete);



	SectionForceDeformation *theSections[2];
	theSections[1] = &steelSec;
	theSections[0] = &concSec;
	SectionForceDeformation *theSection = new ParallelSection(tag, 2, theSections);



	delete [] theNDMat;
	delete [] theUniMat;
	return theSection;
    }

    static void* OPS_RCTBeamSectionUniMat2d()
    {
	if (OPS_GetNumRemainingInputArgs() < 18) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section RCTBeamSection2d tag? coreTag? coverTag? steelTag? d? bw? beff? hf? Atop? Abottom? flcov? wcov? Nflcover? Nwcover? Nflcore? Nwcore? NsteelTop?  NsteelBottom?" << endln;
	    return 0;
	}

	// int
	int numdata = 4;
	int idata[6];
	if (OPS_GetIntInput(&numdata, idata) < 0) {
	    opserr << "WARNING invalid section RCTBeamSection2d int inputs" << endln;
	    return 0;
	}
	int tag = idata[0];
	int coreTag = idata[1];
	int coverTag = idata[2];
	int steelTag = idata[3];

	// double
	numdata = 8;
	double data[8];
	if (OPS_GetDoubleInput(&numdata, data) < 0) {
	    opserr << "WARNING invalid double inputs" << endln;
	    opserr << "RCTBeamSection2d section: " << tag << endln;
	    return 0;
	}

	double d = data[0];
	double bw = data[1];
	double beff = data[2];
	double hf = data[3];
	double Atop = data[4];
	double Abottom = data[5];
	double flcov = data[6];
	double wcov = data[7];

	// int
	numdata = 6;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
	    opserr << "WARNING invalid section RCTBeamSection2d int inputs" << endln;
	    return 0;
	}
	int Nflcover = idata[0];
	int Nwcover = idata[1];
	int Nflcore = idata[2];
	int Nwcore = idata[3];
	int NsteelTop = idata[4];
	int NsteelBottom = idata[5];

	UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(steelTag);
	if (theSteel == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "material: " << steelTag;
	    opserr << "\nRCTBeamSection2d section: " << tag << endln;
	    return 0;
	}

	RCTBeamSectionIntegration
	    rctbeamsect(d, bw, beff, hf, Atop, Abottom, flcov, wcov,
			Nflcover, Nwcover, Nflcore, Nwcore,
			NsteelTop, NsteelBottom);

	UniaxialMaterial *theCore = OPS_getUniaxialMaterial(coreTag);
	if (theCore == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "material: " << coreTag;
	    opserr << "\nRCTBeamSection2d section: " << tag << endln;
	    return 0;
	}

	UniaxialMaterial *theCover = OPS_getUniaxialMaterial(coverTag);
	if (theCover == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "material: " << coreTag;
	    opserr << "\nRCTBeamSection2d section: " << tag << endln;
	    return 0;
	}

	int numFibers = rctbeamsect.getNumFibers();

	UniaxialMaterial **theUniMat = new UniaxialMaterial *[numFibers];

	rctbeamsect.arrangeFibers(theUniMat, theCore, theCover, theSteel);

	SectionForceDeformation *theSection = new FiberSection2d(tag, numFibers, theUniMat, rctbeamsect);

	delete [] theUniMat;

	return theSection;
    }

    static void* OPS_RCCircularSection()
    {
        if (OPS_GetNumRemainingInputArgs() < 13) {
            opserr << "WARNING insufficient arguments\n";
            opserr << "Want: section RCCircularSection tag? coreTag? coverTag? steelTag? d? cover? As? NringsCore? NringsCover? Nwedges? Nsteel? -GJ GJ <or> -torsion matTag\n";
            return 0;
        }

	int idata[8];
	double ddata[3];

	int numdata = 4;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
	    opserr << "WARNING invalid section RCCircularSection input\n";
	    return 0;
	}

	numdata = 3;
	if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
	    opserr << "WARNING invalid section RCCircularSection input\n";
	    return 0;
	}

	numdata = 4;
	if (OPS_GetIntInput(&numdata, &idata[4]) < 0) {
	    opserr << "WARNING invalid section RCCircularSection input\n";
	    return 0;
	}
	
        int tag=idata[0], coreTag=idata[1], coverTag=idata[2], steelTag=idata[3];
        double d=ddata[0], cover=ddata[1], As=ddata[2];
        int ncore=idata[4], ncover=idata[5], nwedge=idata[6], nsteel=idata[7];

        UniaxialMaterial *theCore = OPS_getUniaxialMaterial(coreTag);
        
        if (theCore == 0) {
            opserr << "WARNING uniaxial material does not exist\n";
            opserr << "material: " << coreTag; 
            opserr << "\nRCCircularSection section: " << tag << endln;
            return 0;
        }
        
        UniaxialMaterial *theCover = OPS_getUniaxialMaterial(coverTag);
        
        if (theCover == 0) {
            opserr << "WARNING uniaxial material does not exist\4n";
            opserr << "material: " << coverTag; 
            opserr << "\nRCCircularSection section: " << tag << endln;
            return 0;
        }
        
        UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(steelTag);

        if (theSteel == 0) {
            opserr << "WARNING uniaxial material does not exist\n";
            opserr << "material: " << steelTag; 
            opserr << "\nRCCircularSection section: " << tag << endln;
            return 0;
        }
        
        RCCircularSectionIntegration rcsect(d, As, cover, ncore, ncover, nwedge, nsteel);

        int numFibers = rcsect.getNumFibers();

        UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];

        rcsect.arrangeFibers(theMats, theCore, theCover, theSteel);

	UniaxialMaterial *torsion = 0;
	const char* opt = OPS_GetString();
	numdata = 1;
	bool deleteTorsion = false;
	if (strcmp(opt, "-GJ") == 0) {
	  double GJ;
	  if (OPS_GetDoubleInput(&numdata, &GJ) < 0) {
	    opserr << "WARNING: failed to read GJ\n";
	    return 0;
	  }
	  torsion = new ElasticMaterial(0,GJ);
	  deleteTorsion = true;
	}
	if (strcmp(opt, "-torsion") == 0) {
	  int torsionTag;
	  if (OPS_GetIntInput(&numdata, &torsionTag) < 0) {
	    opserr << "WARNING: failed to read torsion\n";
	    return 0;
	  }
	  torsion = OPS_getUniaxialMaterial(torsionTag);
	}
	if (torsion == 0) {
	  opserr << "WARNING torsion not speified for RCCircularSection\n";
	  opserr << "\nRCCircularSection section: " << tag << endln;
	  return 0;
	}

        // Parsing was successful, allocate the section
        SectionForceDeformation* theSection = new FiberSection3d(tag, numFibers, theMats, rcsect, *torsion);
	if (deleteTorsion)
	  delete torsion;

        delete [] theMats;

	return theSection;
    }

    static void* OPS_RCTunnelSection()
    {
        if (OPS_GetNumRemainingInputArgs() < 13) {
            opserr << "WARNING insufficient arguments\n";
            opserr << "Want: section RCTunnelSection tag? concreteTag? steelTag? d? h? coverinner? coverouter? Asinner? Asouter? Nrings? Nwedges? Nbarsinner? Nbarsouter?\n";
            return 0;
        }

	int idata[8];
	double ddata[6];

	int numdata = 3;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
	    opserr << "WARNING invalid section RCTunnelSection input\n";
	    return 0;
	}

	numdata = 6;
	if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
	    opserr << "WARNING invalid section RCTunnelSection input\n";
	    return 0;
	}

	numdata = 4;
	if (OPS_GetIntInput(&numdata, &idata[4]) < 0) {
	    opserr << "WARNING invalid section RCTunnelSection input\n";
	    return 0;
	}
	
        int tag=idata[0], concreteTag=idata[1], steelTag=idata[2];
        double d=ddata[0], h=ddata[1], coverinner=ddata[2], coverouter=ddata[3], Asinner=ddata[4], Asouter=ddata[5];
        int nring=idata[3], nwedge=idata[4], nbarsinner=idata[5], nbarsouter=idata[6];

        UniaxialMaterial *theConcrete = OPS_getUniaxialMaterial(concreteTag);
        
        if (theConcrete == 0) {
            opserr << "WARNING uniaxial material does not exist\n";
            opserr << "material: " << concreteTag; 
            opserr << "\nRCTunnelSection section: " << tag << endln;
            return 0;
        }
        
        UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(steelTag);

        if (theSteel == 0) {
            opserr << "WARNING uniaxial material does not exist\n";
            opserr << "material: " << steelTag; 
            opserr << "\nRCTunnelSection section: " << tag << endln;
            return 0;
        }
        
        RCTunnelSectionIntegration rcsect(d, h, Asinner, Asouter, coverinner, coverouter,
					  nring, nwedge, nbarsinner, nbarsouter);

        int numFibers = rcsect.getNumFibers();

        UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];

        rcsect.arrangeFibers(theMats, theConcrete, theSteel);

	UniaxialMaterial *torsion = 0;
	const char* opt = OPS_GetString();
	numdata = 1;
	bool deleteTorsion = false;
	if (strcmp(opt, "-GJ") == 0) {
	  double GJ;
	  if (OPS_GetDoubleInput(&numdata, &GJ) < 0) {
	    opserr << "WARNING: failed to read GJ\n";
	    return 0;
	  }
	  torsion = new ElasticMaterial(0,GJ);
	  deleteTorsion = true;
	}
	if (strcmp(opt, "-torsion") == 0) {
	  int torsionTag;
	  if (OPS_GetIntInput(&numdata, &torsionTag) < 0) {
	    opserr << "WARNING: failed to read torsion\n";
	    return 0;
	  }
	  torsion = OPS_getUniaxialMaterial(torsionTag);
	}
	if (torsion == 0) {
	  opserr << "WARNING torsion not speified for RCCircularSection\n";
	  opserr << "\nRCTunnelSection section: " << tag << endln;
	  return 0;
	}

        // Parsing was successful, allocate the section
        SectionForceDeformation* theSection = new FiberSection3d(tag, numFibers, theMats, rcsect, *torsion);

        delete [] theMats;
	if (deleteTorsion)
	  delete torsion;

	return theSection;
    }

    static int setUpFunctions(void)
    {
	functionMap.insert(std::make_pair("Elastic", &OPS_ElasticSection));
	functionMap.insert(std::make_pair("Fiber", &OPS_FiberSection));
	functionMap.insert(std::make_pair("FiberThermal", &OPS_FiberSectionThermal));
	functionMap.insert(std::make_pair("fiberSec", &OPS_FiberSection));
	functionMap.insert(std::make_pair("NDFiber", &OPS_NDFiberSection));
	functionMap.insert(std::make_pair("Uniaxial", &OPS_UniaxialSection));
	functionMap.insert(std::make_pair("Generic1D", &OPS_UniaxialSection));
	functionMap.insert(std::make_pair("Generic1d", &OPS_UniaxialSection));
	functionMap.insert(std::make_pair("ElasticMembranePlateSection", &OPS_ElasticMembranePlateSection));
	functionMap.insert(std::make_pair("PlateFiber", &OPS_MembranePlateFiberSection));
	functionMap.insert(std::make_pair("ElasticWarpingShear", &OPS_ElasticWarpingShearSection2d));
	functionMap.insert(std::make_pair("ElasticTube", &OPS_ElasticTubeSection3d));
	functionMap.insert(std::make_pair("Tube", &OPS_TubeSection));
	functionMap.insert(std::make_pair("WFSection2d", &OPS_WFSection2d));	
	functionMap.insert(std::make_pair("WSection2d", &OPS_WFSection2d));
	functionMap.insert(std::make_pair("RCSection2d", &OPS_RCSection2d));
	functionMap.insert(std::make_pair("RCTBeamSection2d", &OPS_RCTBeamSection2d));
	functionMap.insert(std::make_pair("RCTBeamSectionUniMat2d", &OPS_RCTBeamSectionUniMat2d));
	functionMap.insert(std::make_pair("Parallel", &OPS_ParallelSection));
	functionMap.insert(std::make_pair("Aggregator", &OPS_SectionAggregator));
	functionMap.insert(std::make_pair("AddDeformation", &OPS_SectionAggregator));
	functionMap.insert(std::make_pair("ElasticPlateSection", &OPS_ElasticPlateSection));
	functionMap.insert(std::make_pair("PlateFiber", &OPS_MembranePlateFiberSection));
	functionMap.insert(std::make_pair("LayeredShell", &OPS_LayeredShellFiberSection));
	functionMap.insert(std::make_pair("Bidirectional", &OPS_Bidirectional));
	functionMap.insert(std::make_pair("Elliptical", &OPS_Elliptical2));	
	functionMap.insert(std::make_pair("Isolator2spring", &OPS_Isolator2spring));
	functionMap.insert(std::make_pair("RCCircularSection", &OPS_RCCircularSection));
	functionMap.insert(std::make_pair("RCTunnelSection", &OPS_RCTunnelSection));

	return 0;
    }
}


int OPS_Section()
{
    theActiveFiberSection2d = 0;
    theActiveFiberSection3d = 0;
    theActiveFiberSectionWarping3d = 0;    
    theActiveNDFiberSection2d = 0;
    theActiveNDFiberSection3d = 0;

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
	theActiveNDFiberSection2d = 0;
	theActiveNDFiberSection3d = 0;

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

    } else if (theActiveFiberSection3d != 0 || theActiveFiberSectionWarping3d != 0 || theActiveFiberSection3dThermal!=0) {

	theFiber = (UniaxialFiber3d*) OPS_UniaxialFiber3d();

    } else if (theActiveNDFiberSection2d != 0) {

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

    } else if (theActiveNDFiberSection2d != 0) {

	res = theActiveNDFiberSection2d->addFiber(*theFiber);

    } else if (theActiveNDFiberSection3d != 0) {

	res = theActiveNDFiberSection3d->addFiber(*theFiber);

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
	opserr<<"ERROR unknow patch type\n";
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
	}

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
	opserr<<"ERROR unknow layer type\n";
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
	}

    }

    delete [] reinfBar;
    delete theLayer;


    return 0;
}
