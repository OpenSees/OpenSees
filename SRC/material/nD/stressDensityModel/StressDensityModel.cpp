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

// Written: Saumyasuchi Das
//          May 2013, University of Canterbury
// Updated: Chris McGann
//          June 2015, Washington State University

#include <StressDensityModel.h>
#include <StressDensityModel2D.h>
#include <StressDensityModel3D.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

#define OPS_Export 
OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_StressDensityMaterial)
{
	static int numStressDensityModel = 0;

	if(numStressDensityModel == 0) {
		opserr << "StressDensityModel nDmaterial - Written: Saumyasuchi Das, U.Canterbury\n" << endln;
		numStressDensityModel++;
	}

  // Pointer to an NDmaterial that will be returned
    NDMaterial *theMaterial = 0;
    const char *variable[]={
        "matTag",
		"Density",
		"Void Ratio",
		"Const A",
		"Exponent m",
		"Poisson's ratio",
		"Alpha1",
        "Beta1",
		"Alpha2",
		"Beta2",
		"Alpha3",
        "Beta3",
		"Degradation",
		"Mumin",
		"Mucyclic",
		"Dilatancy Strain",
		"Mumax",
        "Atmospheric pressure",
		"<SSL void ratio at P1>",
		"<SSL void ratio at P2>",
		"<SSL void ratio at P3>",
		"<SSL void ratio at P4>",
		"<SSL void ratio at P5>",
		"<SSL void ratio at P6>",
		"<SSL void ratio at P7>",
		"<SSL void ratio at P8>",
		"<SSL void ratio at P9>",
		"<SSL void ratio at P10>",
		"<HSL void ratio at all P>",
		"<Reference pressure P1>",
		"<Reference pressure P2>",
		"<Reference pressure P3>",
		"<Reference pressure P4>",
		"<Reference pressure P5>",
		"<Reference pressure P6>",
		"<Reference pressure P7>",
		"<Reference pressure P8>",
		"<Reference pressure P9>",
		"<Reference pressure P10>"};
					
    int numArgs = OPS_GetNumRemainingInputArgs();

    if (numArgs < 18) {
      opserr << "ERROR: Insufficient mandatory arguments for StressDensity material " << endln;
      
      for(int i=numArgs; i < 18; i++)
	opserr << "Missing: " << variable[i] << endln;
      
      return 0;
    } else if (numArgs > 18 && numArgs < 29) {
      opserr << "ERROR: Insufficient optional void ratio arguments for StressDensity materal" << endln;
      opserr << "All ten SSL values and single HSL value must be specified if defaults are not used" << endln;
      
      for(int i=numArgs; i < 29; i++)
	opserr << "Missing: " << variable[i] << endln;
      
      return 0;
    } else if (numArgs > 29 && numArgs < 39) {
      opserr << "ERROR: Insufficient optional SSL pressure arguments for StressDensity materal" << endln;
      opserr << "All ten pressure values must be specified if defaults are not used" << endln;
      
      for(int i=numArgs; i < 39; i++)
	opserr << "Missing: " << variable[i] << endln;
      
      return 0;
    } else if (numArgs > 39 && numArgs < 45) {
      opserr << "WARNING: Initial backstress tensor components have been specified for StressDensity mat" << endln;
      opserr << "If this is not desired, please ensure material arguments are correct" << endln;
    } else if (numArgs > 45) {
      opserr << "ERROR: Too many input arguments specified for StressDensity material" << endln;
        return 0;
    }
  
	int tag;
	double dData[45];

	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid nDMaterial StressDensity material tag" << endln;
		return 0;
	}

	numData = numArgs-1;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid material data for nDMaterial StressDensity material with tag: " << tag << endln;
		return 0;
	}

	if (numArgs == 18) {
        theMaterial = new StressDensityModel(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
		dData[15],dData[16],dData[17]);
    } else if (numArgs == 29) {
        theMaterial = new StressDensityModel(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
		dData[15],dData[16],dData[17],dData[18],dData[19],dData[20],dData[21],dData[22],dData[23],dData[24],
		dData[25],dData[26],dData[27],dData[28]);
    } else if (numArgs == 39) {
        theMaterial = new StressDensityModel(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
		dData[15],dData[16],dData[17],dData[18],dData[19],dData[20],dData[21],dData[22],dData[23],dData[24],
		dData[25],dData[26],dData[27],dData[28],dData[29],dData[30],dData[31],dData[32],dData[33],dData[34],
		dData[35],dData[36],dData[37],dData[38]);
    } else if (numArgs == 40) {
		theMaterial = new StressDensityModel(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
		dData[15],dData[16],dData[17],dData[18],dData[19],dData[20],dData[21],dData[22],dData[23],dData[24],
		dData[25],dData[26],dData[27],dData[28],dData[29],dData[30],dData[31],dData[32],dData[33],dData[34],
		dData[35],dData[36],dData[37],dData[38],dData[39]);
    } else if (numArgs == 41) {
		theMaterial = new StressDensityModel(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
		dData[15],dData[16],dData[17],dData[18],dData[19],dData[20],dData[21],dData[22],dData[23],dData[24],
		dData[25],dData[26],dData[27],dData[28],dData[29],dData[30],dData[31],dData[32],dData[33],dData[34],
		dData[35],dData[36],dData[37],dData[38],dData[39],dData[40]);
    } else if (numArgs == 42) {
		theMaterial = new StressDensityModel(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
		dData[15],dData[16],dData[17],dData[18],dData[19],dData[20],dData[21],dData[22],dData[23],dData[24],
		dData[25],dData[26],dData[27],dData[28],dData[29],dData[30],dData[31],dData[32],dData[33],dData[34],
		dData[35],dData[36],dData[37],dData[38],dData[39],dData[40],dData[41]);
    } else if (numArgs == 43) {
		theMaterial = new StressDensityModel(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
		dData[15],dData[16],dData[17],dData[18],dData[19],dData[20],dData[21],dData[22],dData[23],dData[24],
		dData[25],dData[26],dData[27],dData[28],dData[29],dData[30],dData[31],dData[32],dData[33],dData[34],
		dData[35],dData[36],dData[37],dData[38],dData[39],dData[40],dData[41],dData[42]);
    } else if (numArgs == 44) {
		theMaterial = new StressDensityModel(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
		dData[15],dData[16],dData[17],dData[18],dData[19],dData[20],dData[21],dData[22],dData[23],dData[24],
		dData[25],dData[26],dData[27],dData[28],dData[29],dData[30],dData[31],dData[32],dData[33],dData[34],
		dData[35],dData[36],dData[37],dData[38],dData[39],dData[40],dData[41],dData[42],dData[43]);
    } else if (numArgs == 45) {
		theMaterial = new StressDensityModel(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
		dData[15],dData[16],dData[17],dData[18],dData[19],dData[20],dData[21],dData[22],dData[23],dData[24],
		dData[25],dData[26],dData[27],dData[28],dData[29],dData[30],dData[31],dData[32],dData[33],dData[34],
		dData[35],dData[36],dData[37],dData[38],dData[39],dData[40],dData[41],dData[42],dData[43]);
    } 

	if (theMaterial == 0) {
	    opserr << "WARNING ran out of memory for nDMaterial StressDensity material with tag: " << tag << endln;
    }

    return theMaterial;
}

// full constructor
StressDensityModel::StressDensityModel(int tag, int classTag, double constDensity,
						                   // SD model  parameters		
						                   double initialVoidRatio,	double constA, double exponentN,		
						                   double poissonRatio, double constAlpha1, double constBeta1,
						                   double constAlpha2, double constBeta2, double constAlpha3,
                                           double constBeta3, double constDegradation, double constMumin,
						                   double constMucyclic, double constDilatancyStrain,	
						                   double constMumax, double constPatm, 
						                   // steady state line void ratio
						                   double constsslvoidatP1, double constsslvoidatP2, double constsslvoidatP3,
						                   double constsslvoidatP4, double constsslvoidatP5, double constsslvoidatP6,
						                   double constsslvoidatP7, double constsslvoidatP8, double constsslvoidatP9,
						                   double constsslvoidatP10,
						                   // hydrostatic state line void ratio
						                   double consthslvoid,
						                   // reference pressures
						                   double constP1, double constP2, double constP3, double constP4, double constP5,
						                   double constP6, double constP7, double constP8, double constP9, double constP10,
						                   // offset of the failure surface
						                   double constRxx, double constRyy, double constRzz,
						                   double constRxy, double constRyz, double constRzx)
  : NDMaterial(tag,classTag)
{
    // material stage (0 = elastic, 1 = elastoplastic)
    theStage = 0;
    // mass density
    theDensity = constDensity;
    // model parameters
	modelParameter[0]  = initialVoidRatio;
    modelParameter[1]  = constA;
	modelParameter[2]  = exponentN;
	modelParameter[3]  = poissonRatio;
	modelParameter[4]  = constAlpha1;
	modelParameter[5]  = constBeta1;
	modelParameter[6]  = constAlpha2;
	modelParameter[7]  = constBeta2;
	modelParameter[8]  = constAlpha3;
	modelParameter[9]  = constBeta3;
	modelParameter[10] = constDegradation;
	modelParameter[11] = constMumin;
	modelParameter[12] = constMucyclic;
	modelParameter[13] = constDilatancyStrain;
	modelParameter[14] = constMumax;
    modelParameter[15] = constPatm;
    // steady state line void ratios
	sslVoidratio[0] = constsslvoidatP1;
	sslVoidratio[1] = constsslvoidatP2;
	sslVoidratio[2] = constsslvoidatP3;
	sslVoidratio[3] = constsslvoidatP4;
	sslVoidratio[4] = constsslvoidatP5;
	sslVoidratio[5] = constsslvoidatP6;
	sslVoidratio[6] = constsslvoidatP7;
	sslVoidratio[7] = constsslvoidatP8;
	sslVoidratio[8] = constsslvoidatP9;
	sslVoidratio[9] = constsslvoidatP10;
    // hydrostatic state line void ratios (constant with p)
	hslVoidratio[0] = consthslvoid;
	hslVoidratio[1] = consthslvoid;
	hslVoidratio[2] = consthslvoid;
	hslVoidratio[3] = consthslvoid;
	hslVoidratio[4] = consthslvoid;
	hslVoidratio[5] = consthslvoid;
	hslVoidratio[6] = consthslvoid;
	hslVoidratio[7] = consthslvoid;
	hslVoidratio[8] = consthslvoid;
	hslVoidratio[9] = consthslvoid;
    // reference pressures for void ratios
	refPressure[0] = constP1;
	refPressure[1] = constP2;
	refPressure[2] = constP3;
	refPressure[3] = constP4;
	refPressure[4] = constP5;
	refPressure[5] = constP6;
	refPressure[6] = constP7;
	refPressure[7] = constP8;
	refPressure[8] = constP9;
	refPressure[9] = constP10;
    // initial offset for failure surface
	refOrigin[0] = constRxx;
	refOrigin[1] = constRyy;
	refOrigin[2] = constRzz;
	refOrigin[3] = constRxy;
	refOrigin[4] = constRyz;
	refOrigin[5] = constRzx;
}

// null constructor
StressDensityModel::StressDensityModel()
  : NDMaterial()
{
    theStage = 0;
    theDensity = 0.0;
    for (int i=0;i<16;i++) {
        modelParameter[i] = 0.0;
    }
    for (int i=0;i<10;i++) {
        sslVoidratio[i] = 0.0;
        hslVoidratio[i] = 0.0;
        refPressure[i] = 0.0;
    }
    for (int i=0;i<6;i++) {
	    refOrigin[i] = 0.0;
    }
}

// destructor
StressDensityModel::~StressDensityModel()
{
}

double 
StressDensityModel::getRho(void) 
{
	return theDensity;
}

int 
StressDensityModel::setParameter(const char **argv, int argc, Parameter &param)
{	
    if (strcmp(argv[0],"updateMaterialStage") == 0) {
        return param.addObject(1, this);
    } else if (strcmp(argv[0],"materialState") == 0) {
        return param.addObject(5, this);
    } else if (strcmp(argv[0],"poissonRatio") == 0) {
        return param.addObject(7,this);
    } else {
        opserr << "WARNING: invalid parameter command StressDensityModel nDMaterial tag: " << this->getTag() << endln;
        return -1;
    }

    return -1;
}

int 
StressDensityModel::updateParameter(int parameterID, Information &info)
{	
	if (parameterID == 1) {
		theStage = info.theInt;
	} else if (parameterID == 5) {
        theStage = (int)info.theDouble;
    } else if (parameterID == 7) {
        modelParameter[3] = info.theDouble;
    }

    return 0;
}

void 
StressDensityModel::Print(OPS_Stream &s, int flag)
{
	s<<"StressDensityModel::Tag "<<this->getTag()<<endln;
	s<<"Material Stage: "<< theStage;
}

int 
StressDensityModel::commitState(void) {

	return 0; //Subclass responsibility
}

int
StressDensityModel::revertToLastCommit(void)
{
	return 0;
}

int 
StressDensityModel::revertToStart(void)
{
	return 0;
}

NDMaterial *
StressDensityModel::getCopy(void)
{
	return 0; //Subclass responsibility
}

const char *
StressDensityModel::getType(void) const
{
	return 0; //Subclass responsibility
}

int
StressDensityModel::getOrder(void) const 
{
    return 0; //Subclass responsibility
}

NDMaterial *
StressDensityModel::getCopy(const char *code) 
{
	if (strcmp(code, "PlaneStrain")==0||strcmp(code,"2D")==0) {
		StressDensityModel2D *theCopy;
		theCopy = new StressDensityModel2D(this->getTag(), theDensity, modelParameter[0], modelParameter[1],
											 modelParameter[2], modelParameter[3], modelParameter[4],
											 modelParameter[5], modelParameter[6], modelParameter[7],
											 modelParameter[8], modelParameter[9], modelParameter[10],
											 modelParameter[11], modelParameter[12], modelParameter[13],
											 modelParameter[14], modelParameter[15], sslVoidratio[0],
											 sslVoidratio[1], sslVoidratio[2], sslVoidratio[3], sslVoidratio[4],
											 sslVoidratio[5], sslVoidratio[6], sslVoidratio[7], sslVoidratio[8],
											 sslVoidratio[9], hslVoidratio[0], refPressure[0], refPressure[1],
											 refPressure[2], refPressure[3], refPressure[4], refPressure[5],
											 refPressure[6], refPressure[7], refPressure[8], refPressure[9],
											 refOrigin[0], refOrigin[1], refOrigin[3]);
		return theCopy;
	} else if  (strcmp(code, "ThreeDimensional")==0||strcmp(code,"3D")==0) {
		StressDensityModel3D *theCopy;
		theCopy = new StressDensityModel3D(this->getTag(), theDensity, modelParameter[0], modelParameter[1],
											 modelParameter[2], modelParameter[3], modelParameter[4],
											 modelParameter[5], modelParameter[6], modelParameter[7],
											 modelParameter[8], modelParameter[9], modelParameter[10],
											 modelParameter[11], modelParameter[12], modelParameter[13],
											 modelParameter[14], modelParameter[15], sslVoidratio[0],
											 sslVoidratio[1], sslVoidratio[2], sslVoidratio[3], sslVoidratio[4],
											 sslVoidratio[5], sslVoidratio[6], sslVoidratio[7], sslVoidratio[8],
											 sslVoidratio[9], hslVoidratio[0], refPressure[0], refPressure[1],
											 refPressure[2], refPressure[3], refPressure[4], refPressure[5],
											 refPressure[6], refPressure[7], refPressure[8], refPressure[9],
											 refOrigin[0], refOrigin[1], refOrigin[2], refOrigin[3], refOrigin[4],
                                             refOrigin[5]);
				
		return theCopy;
	} else {
		opserr<<"StressDensityModel::getCopy failed to get copy: " << code << endln;
	  	return 0;
	}
}

Response *
StressDensityModel::setResponse (const char **argv, int argc, OPS_Stream &output)
{
	output.tag("NdMaterialOutput");
	output.attr("matType",this->getClassType());
	output.attr("matTag",this->getTag());

	if(strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0) { 
		return new MaterialResponse (this, 1, this->getStress());
	} else if(strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0) {
		return new MaterialResponse(this, 2, this->getStrain());
    // add new recorder here for state things (voids ratio, etc...)
	} else {
		return 0;
    }
}

int 
StressDensityModel::getResponse (int responseID, Information &matInformation)
{
	switch (responseID) {
        case -1:
            return -1;
	    case 1:
		    if (matInformation.theVector != 0) 
			    *(matInformation.theVector) = getStress();
		    return 0;
    	case 2:
		    if (matInformation.theVector != 0) 
			    *(matInformation.theVector) = getStrain();
		    return 0;
	    default:
		    return -1;
  }
}

int
StressDensityModel::sendSelf(int commitTag, Channel &theChannel)
{
	return 0; //subclass responsibility
}

int 
StressDensityModel::recvSelf(int commitTag, Channel &theChannel,FEM_ObjectBroker &theBroker)
{
	return 0; //subclass responsibility
}
