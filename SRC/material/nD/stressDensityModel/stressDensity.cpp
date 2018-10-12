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

// Written: Chris McGann
//          July 2018, University of Canterbury
//

#include <stressDensity.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

#define OPS_Export
OPS_Export void *
OPS_StressDensityMaterial(void)
{
    static int numStressDensity = 0;

    if(numStressDensity == 0) {
        opserr << "stressDensity nDMaterial - Written: Saumyasuchi Das, U.Canterbury\n" << endln;
        numStressDensity++;
    }

    NDMaterial *theMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();

    if (numArgs < 18) {
        opserr << "ERROR stressDensity nDMaterial: Insufficient mandatory input arguments" << endln;
        opserr << "WANT: nDmaterial stressDensity tag mDen eNot A n nu a1 b1 a2 b2 a3 b3 fd muNot muCyc sc M patm" endln;
        return 0;
    } else if (numArgs > 18 && numArgs < 27) {
        opserr << "ERROR: stressDensity nDMaterial: Insufficient optional SSL void ratio arguments" << endln;
        opserr << "ssl1-ssl7, hsl, and pmin must all be specified if defaults are not used" << endln;
        return 0;
    }

    int tag;
    double dData[26];

    int numData = 1;
    if (OPS_GetInt(&numData, &tag) != 0) {
        opserr << "WARNING: invalied nDMaterial stressDensity material tag" << endln;
        return 0;
    }

    numData = numArgs-1;
    if (OPS_GetDouble(&numData, dData) !=0) {
        opserr << "WARNING: invalid material data for nDMaterial stressDensity with tag: " << tag << endln;
        return 0;
    }

    if (numArgs == 18) {
        theMaterial = new stressDensity(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		    dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
		    dData[15], dData[16]);
    } else if (numArgs == 27) {
        theMaterial = new stressDensity(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6],
		    dData[7], dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], dData[16],
		    dData[17], dData[18], dData[19], dData[20], dData[21], dData[22], dData[23], dData[24], dData[25]);
    }

    if (theMaterial == 0) {
        opserr << "WARNING: ran out of memory for nDMaterial stressDensity with tag: " << tag << endln;
    }

    return theMaterial;
}

// full constructor
stressDensity::stressDensity(int tag, int classTag, double massDen,
			                 // SD model  parameters		
				             double eInit, double constA, double expN, double nu, double a1, double b1,	
                             double a2, double b2, double a3, double b3, double fd, double muNot,
                             double muCyc, double sc, double M, double patm, 
                             // steady state line void definition (optional input arguments)
                             double ssl1, double ssl2, double ssl3, double ssl4, double ssl5, 
                             double ssl6, double ssl7, double hsl, double pmin)
  : NDMaterial(tag,ND_TAG_stressDensity),
    stressCurrent(3),
    stressNext(3),
    strainCurrent(3),
    strainNext(3),
    materialParam(25),
    initialTangent(3,3),
    currentTangent(3,3)
{
    massDensity = massDen;
    materialParam(0)  = eInit;
    materialParam(1)  = nu;
    materialParam(2)  = constA;
    materialParam(3)  = expN;
    materialParam(4)  = M;
    materialParam(5)  = muNot;
    materialParam(6)  = muCyc;
    materialParam(7)  = sc;
    materialParam(8)  = a1;
    materialParam(9)  = b1;
    materialParam(10) = a3;
    materialParam(11) = b3;
    materialParam(12) = a2;
    materialParam(13) = b2;
    materialParam(14) = fd;
    materialParam(15) = pmin;
    materialParam(16) = hsl;
    materialParam(17) = ssl1;
    materialParam(18) = ssl2;
    materialParam(19) = ssl3;
    materialParam(20) = ssl4;
    materialParam(21) = ssl5;
    materialParam(22) = ssl6;
    materialParam(23) = ssl7;
    materialParam(24) = patm;

    // initialise variables
    this->initialise();
}

// null constructor
stressDensity::stressDensity()
  : NDMaterial(),
    stressCurrent(3),
    stressNext(3),
    strainCurrent(3),
    strainNext(3),
    materialParam(25),
    initialTangent(3,3),
    currentTangent(3,3)
{
    theStage = 0;
    for (int i=0; i<24; i++) {
        materialParam(i) = 0.0;
    }
    this->initialise();
}

// destructor
stressDensity::~stressDensity()
{
}

int
stressDensity::commitState(void)
{
    // commit stress and strain
	stressCurrent = stressNext;
	strainCurrent = strainNext;
    // update tangent terms
    pInit = -0.5*(stressCurrent(0)+stressCurrent(1));
    pFlag = 1;
    this->calInitialTangent();
    currentTangent = initialTangent;

    // step and iteration counters
    if (theStage == 1) {
        istep++;
    }
    iiter = 1;

	return 0;
}

int 
stressDensity::revertToLastCommit(void)
{
    return 0;
}

int
stressDensity::revertToStart(void)
{
	// added for InitialStateAnalysis
	if (ops_InitialStateAnalysis) {
		// do nothing, keep state variables from last step
	} else {
		// normal call for revertToStart (not initialStateAnalysis)
    	this->initialise();
	}
    return 0;
}

NDMaterial*
stressDensity::getCopy(void)
{
    stressDensity *copy = new stressDensity(*this);
    return copy;
}

NDMaterial*
stressDensity::getCopy(const char *type)
{
    if (strcmp(type,"PlaneStrain2D") == 0 || strcmp(type, "PlaneStrain") == 0) {
        stressDensity *copy = new stressDensity(*this);
        return copy;
    } else if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
        opserr << "This is a 2D model and it is not compatible with " << type << endln;
        return 0;
    } else {
        opserr << "stressDensity nDMaterial: getCopy failed to get copy, type: " << type << endln;
        return 0;
    }
}

const char*
stressDensity::getType(void) const
{
    return "PlaneStrain";
}

int
stressDensity::getOrder(void) const
{
    return 3;
}

Response*
stressDensity::setResponse(const char **argv, int argc, OPS_Stream &output)
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
stressDensity::getResponse (int responseID, Information &matInformation)
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
stressDensity::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;
    static Vector vData(798);

    vData(0)  = this->getTag();
    vData(1)  = theStage;
    vData(2)  = massDensity;
    vData(3)  = materialParam(0);
    vData(4)  = materialParam(1);
    vData(5)  = materialParam(2);
    vData(6)  = materialParam(3);
    vData(7)  = materialParam(4);
    vData(8)  = materialParam(5);
    vData(9)  = materialParam(6);
    vData(10) = materialParam(7);
    vData(11) = materialParam(8);
    vData(12) = materialParam(9);
    vData(13) = materialParam(10);
    vData(14) = materialParam(11);
    vData(15) = materialParam(12);
    vData(16) = materialParam(13);
    vData(17) = materialParam(14);
    vData(18) = materialParam(15);
    vData(19) = materialParam(16);
    vData(20) = materialParam(17);
    vData(21) = materialParam(18);
    vData(22) = materialParam(19);
    vData(23) = materialParam(20);
    vData(24) = materialParam(21);
    vData(25) = materialParam(22);
    vData(26) = materialParam(23);
    vData(27) = materialParam(24);
    vData(28) = pFlag;
    vData(29) = pInit;
    for (int i=0; i<12; i++) {
        vData(30+i) = oths[i];
    }
    for (int i=0; i<100; i++) {
        vData(42+i) = strhs[i];
    }
    for (int i=0; i<280; i++) {
        vData(142+i) = strhs0[i];
    }
    for (int i=0; i<40; i++) {
        vData(422) = etahs[i][0];
        vData(462) = etahs[i][1];
        vData(502) = etahs[i][2];
    }
    for (int i=0; i<80; i++) {
        vData(542) = hdp[i][0];
        vData(622) = hdp[i][1];
        vData(702) = hdp[i][2];
    }
    vData(782) = stressCurrent(0); vData(783) = stressCurrent(1); vData(784) = stressCurrent(2);
    vData(785) = strainCurrent(0); vData(786) = strainCurrent(1); vData(787) = strainCurrent(2);
    vData(788) = initialTangent(0,0); vData(789) = initialTangent(0,1); vData(790) = initialTangent(0,2);
    vData(791) = initialTangent(1,0); vData(792) = initialTangent(1,1); vData(793) = initialTangent(1,2);
    vData(794) = initialTangent(2,0); vData(795) = initialTangent(2,1); vData(796) = initialTangent(2,2);
    vData(797) = istep;
        
    res = theChannel.sendVector(this->getDbTag(), commitTag, vData);
	if (res < 0) {
      opserr << "stressDensity::sendSelf() - failed to send vData\n";
	  return -1;
	}

	return 0;
}

int
stressDensity::recvSelf(int commitTag, Channel &theChannel,FEM_ObjectBroker &theBroker)
{
    int res = 0;

    // place data in a vector
    static Vector vData(798);

	res = theChannel.recvVector(this->getDbTag(), commitTag, vData);
	if (res < 0) {
		opserr << "stressDensity::recvSelf() - failed to recv vData\n";
		return -1;
    }
	
    this->setTag((int)vData(0));
	theStage          = (int)vData(1);
    massDensity       = vData(2);
    materialParam(0)  = vData(3);
    materialParam(1)  = vData(4);
    materialParam(2)  = vData(5);
    materialParam(3)  = vData(6);
    materialParam(4)  = vData(7);
    materialParam(5)  = vData(8);
    materialParam(6)  = vData(9);
    materialParam(7)  = vData(10);
    materialParam(8)  = vData(11);
    materialParam(9)  = vData(12);
    materialParam(10) = vData(13);
    materialParam(11) = vData(14);
    materialParam(12) = vData(15);
    materialParam(13) = vData(16);
    materialParam(14) = vData(17);
    materialParam(15) = vData(18);
    materialParam(16) = vData(19);
    materialParam(17) = vData(20);
    materialParam(18) = vData(21);
    materialParam(19) = vData(22);
    materialParam(20) = vData(23);
    materialParam(21) = vData(24);
    materialParam(22) = vData(25);
    materialParam(23) = vData(26);
    materialParam(24) = vData(27);
    pFlag = (int)vData(28);
    pInit = vData(29);
    for (int i=0; i<12; i++) {
        oths[i] = vData(30+i);
    }
    for (int i=0; i<100; i++) {
        strhs[i] = vData(42+i);
    }
    for (int i=0; i<280; i++) {
        strhs0[i] = vData(142+i);
    }
    for (int i=0; i<40; i++) {
        etahs[i][0] = vData(422);
        etahs[i][1] = vData(462);
        etahs[i][2] = vData(502);
    }
    for (int i=0; i<80; i++) {
        hdp[i][0] = vData(542);
        hdp[i][1] = vData(622);
        hdp[i][2] = vData(702);
    }
    stressCurrent(0) = vData(782); stressCurrent(1) = vData(783); stressCurrent(2) = vData(784);
    strainCurrent(0) = vData(785); strainCurrent(1) = vData(786); strainCurrent(2) = vData(787);
    initialTangent(0,0) = vData(788); initialTangent(0,1) = vData(789); initialTangent(0,2) = vData(790);
    initialTangent(1,0) = vData(791); initialTangent(1,1) = vData(792); initialTangent(1,2) = vData(793);
    initialTangent(2,0) = vData(794); initialTangent(2,1) = vData(795); initialTangent(2,2) = vData(796);
    istep = vData(797);
       
    // set current tangent
    currentTangent = initialTangent; 
    // populate props with model parameters (not all indices used in SDM-UC)
    props[3]  = materialParam(1);
    props[5]  = materialParam(2);
    props[27] = materialParam(3);
    props[28] = materialParam(4);
    props[29] = materialParam(5);
    props[26] = materialParam(6);
    props[30] = materialParam(7);
    props[32] = materialParam(8);
    props[31] = materialParam(9);
    props[34] = materialParam(10);
    props[33] = materialParam(11);
    props[36] = materialParam(12);
    props[35] = materialParam(13);
    props[37] = materialParam(14);
    props[38] = materialParam(15);
    props[39] = materialParam(16);
    props[40] = materialParam(17);
    props[41] = materialParam(18);
    props[42] = materialParam(19);
    props[43] = materialParam(20);
    props[44] = materialParam(21);
    props[45] = materialParam(22);
    props[46] = materialParam(23);
    props[10] = materialParam(0)/(1.0 + materialParam(0));

    return 0;
}

void 
stressDensity::Print(OPS_Stream &s, int flag)
{
	s << "stressDensity Material, tag: " << this->getTag() << endln;
    s << "Type: " << this->getType() << endln;
	s << "Material Stage: " << theStage << endln;
}

int 
stressDensity::setParameter(const char **argv, int argc, Parameter &param)
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
stressDensity::updateParameter(int parameterID, Information &info)
{	
	if (parameterID == 1) {
		theStage = info.theInt;
	} else if (parameterID == 5) {
        theStage = (int)info.theDouble;
    } else if (parameterID == 7) {
        materialParam(1) = info.theDouble;
        props[3] = info.theDouble;
    }

    return 0;
}

void
stressDensity::initialise()
{
    // initialise analysis stage to zero (elastic response)
    theStage = 0;
    // initialise Vector and Matrix variables
    stressCurrent.Zero();
    stressNext.Zero();
    strainCurrent.Zero();
    strainNext.Zero();
    initialTangent.Zero();
    currentTangent.Zero();

    // get the initial material tangent
    pInit = 0.0;
    pFlag = 0;
    this->calInitialTangent();
    // set current tangent as initial tangent to start
    currentTangent = initialTangent;

    // initialise in/out variables for FORTRAN
    for (int i=0; i<4; i++) {
        strsg[i] = 0.0;
        stran[i] = 0.0;
    }
    for (int i=0; i<nstrp; i++) {
        strhs[i] = 0.0;
    }
    for (int i=0; i<3; i++) {
        strhs[i] = 1.0;
    }
    for (int i=0; i<280; i++) {
        strhs0[i] = 0.0;
    }
    for (int i=0; i<40; i++) {
        etahs[i][0] = 0.0;
        etahs[i][1] = 0.0;
        etahs[i][2] = 0.0;
    }
    for (int i=0; i<80; i++) {
        hdp[i][0] = 0.0;
        hdp[i][1] = 0.0;
        hdp[i][2] = 0.0;
    }
    for (int i=0; i<12; i++) {
        oths[i] = 0.0;
    }
    for (int i=0; i<nmats; i++) {
        props[i] = 0.0;
    }
    // populate props with model parameters (not all indices used in SDM-UC)
    props[3]  = materialParam(1);
    props[5]  = materialParam(2);
    props[27] = materialParam(3);
    props[28] = materialParam(4);
    props[29] = materialParam(5);
    props[26] = materialParam(6);
    props[30] = materialParam(7);
    props[32] = materialParam(8);
    props[31] = materialParam(9);
    props[34] = materialParam(10);
    props[33] = materialParam(11);
    props[36] = materialParam(12);
    props[35] = materialParam(13);
    props[37] = materialParam(14);
    props[38] = materialParam(15);
    props[39] = materialParam(16);
    props[40] = materialParam(17);
    props[41] = materialParam(18);
    props[42] = materialParam(19);
    props[43] = materialParam(20);
    props[44] = materialParam(21);
    props[45] = materialParam(22);
    props[46] = materialParam(23);
    // SDM-UC expects porosity as input instead of void ratio
    props[10] = materialParam(0)/(1.0 + materialParam(0));

    // integer inputs for FORTRAN routine
    istep = 1;
    iiter = 1;
}

int 
stressDensity::setTrialStrain(const Vector &strain_from_element) 
{
    strainNext = strain_from_element;
	this->getCurrentStress(); 
	return 0;
}

// unused trial strain rate function
int 
stressDensity::setTrialStrain(const Vector &v, const Vector &r)
{
	return this->setTrialStrain(v);
}

const Matrix&
stressDensity::getTangent(void) 
{
    return currentTangent;
}

const Matrix&
stressDensity::getInitialTangent(void) 
{
    return initialTangent;
}

const Vector&
stressDensity::getStress(void)
{
    return stressNext;
}

const Vector&
stressDensity::getStrain(void)
{
    return strainCurrent;
}

double 
stressDensity::getRho(void) 
{
	return massDensity;
}

void
stressDensity::getCurrentStress(void)
{
    // -------- elastic stage (theStage = 0) ----------------------------------------
    if (theStage !=1) {
        stressNext = stressCurrent + currentTangent*(strainNext-strainCurrent);
        return;
    }

    // -------- elastoplastic stage (theStage == 1) ---------------------------------
    //
    // strsg is the current stress at the start of the step
    strsg[0] = -stressCurrent(0);
    strsg[1] = -stressCurrent(1);
    strsg[2] =  stressCurrent(2);
    strsg[3] = -0.5*(stressCurrent(0)+stressCurrent(1));
    // stran is the strain increment for this step
    stran[0] = -(strainNext(0) - strainCurrent(0));
    stran[1] = -(strainNext(1) - strainCurrent(1));
    stran[2] =  (strainNext(2) - strainCurrent(2))/2.0;

    // alter behaviour based on number of iterations
    if (iiter <= 3) {
        for (int i=0; i<4; i++) {
            strhs0[i]   = strsg[i];
            strhs0[i+4] = stran[i];
        }
        for (int i=8; i<33; i++) {
            strhs0[i] = strhs[i-8];
        }
        for (int i=0; i<3; i++) {
            for (int j=0; j<80; j++) {
                int m = i*80+40+j;
                strhs0[m] = hdp[j][i];
            }
        }
        // don't engage the constitutive model when change in strain is very small
        if (iiter == 3) {
            if (abs(stran[0]) < 1.0e-10 && abs(stran[1]) < 1.0e-10 && abs(stran[2]) < 1.0e-10) {
                stressNext = stressCurrent + currentTangent*(strainNext-strainCurrent);
                return;
            }
        }
    } else {
        for (int i=0; i<4; i++) {
            strsg[i] = strhs0[i];
            stran[i] = strhs0[i+4];
        }
        for (int i=0; i<25; i++) {
            strhs[i] = strhs0[i+8]; 
        }
    }
    // send this information to sdmuc in oths
    oths[10] = iiter;
    oths[11] = istep;

    // FORTRAN subroutine for stress integration
    sdmuc_(strhs, strsg, props, stran, nmats, nstrp,
           istep, iiter, ielem,
           strhs0, etahs, hdp, oths);

    // update iteration counter variable 
    iiter++;

    // update member stress variable from FORTRAN results
    stressNext(0) = -strsg[0];
    stressNext(1) = -strsg[1];
    stressNext(2) =  strsg[2];
    // update material tangent coefficient
    materialParam(2) = props[5];

    // get updated tangent
    pInit = -0.5*(stressNext(0)+stressNext(1));
    this->calInitialTangent();
    currentTangent = initialTangent;
}

void 
stressDensity::calInitialTangent(void)
{
    double nu, G, A, n, eo, patm, fct, afc;

    eo   = materialParam(0);
    nu   = materialParam(1);
    A    = materialParam(2);
    n    = materialParam(3);
    patm = materialParam(24);

    if (materialParam(4) > 0.15 && strhs[12] > 0.02) {
        fct = strhs[12]/0.05;
        if (fct > 1.0) {
            fct = 1.0;
        }
        n = materialParam(3) + (0.85 - materialParam(3))*fct;
    }

    // assume p = patm for initial shear modulus
    if (pFlag == 0) {
	    G = A*patm*(2.17 - eo)*(2.17 - eo)/(1.0 + eo)*(pow(1.0,n));
    
    } else {
        G = A*patm*(2.17 - eo)*(2.17 - eo)/(1.0 + eo)*(pow((pInit/patm),n));
    }

    initialTangent(0,0) = 2.0*G*(1.0+nu)/(3.0*(1-2.0*nu)) + 4.0*G/3.0;
    initialTangent(0,1) = 2.0*G*(1.0+nu)/(3.0*(1-2.0*nu)) - 2.0*G/3.0;
    initialTangent(0,2) = 0.0;
    initialTangent(1,2) = 0.0;
    initialTangent(1,0) = initialTangent(0,1);
    initialTangent(2,0) = initialTangent(0,2);
    initialTangent(1,1) = initialTangent(0,0);
    initialTangent(2,1) = initialTangent(1,2);
    initialTangent(2,2) = G;
}
