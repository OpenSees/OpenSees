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

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 10/19
// Revision: A
//
// Description: This file contains the implementation of the
// ElasticPowerFunc class.

#include <ElasticPowerFunc.h>

#include <Vector.h>
#include <Channel.h>
#include <elementAPI.h>

#include <math.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>


void * OPS_ADD_RUNTIME_VPV(OPS_ElasticPowerFunc)
{
    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc < 5) {
        opserr << "WARNING incorrect num args want: uniaxialMaterial ";
        opserr << "ElasticPowerFunc tag <eta> -coeff c1 c2 ... ";
        opserr << "-exp e1 e2 ... ";
        opserr << "(with at least one pair of (ci,ei) values)\n";
        return 0;
    }
    
    int tag[1];
    double coeffData[64];
    double expData[64];
    double eta = 0.0;
    const char *paraStr;
    
    int numData = 1;
    if (OPS_GetIntInput(&numData,tag) != 0)  {
        opserr << "WARNING invalid uniaxialMaterial ElasticPowerFunc tag\n";
        return 0;
    }
    
    // check if eta is provided (odd number of inputs)
    if ((argc-3)%2 == 1)  {
        numData = 1;
        if (OPS_GetDoubleInput(&numData,&eta) != 0)  {
            opserr << "WARNING invalid eta\n";
            opserr << "uniaxialMaterial ElasticPowerFunc: " << tag[0] << endln;
            return 0;
        }
        argc--;
    }
    
    // get coefficient values
    numData = (argc - 3)/2;
    paraStr = OPS_GetString();
    if (strcmp(paraStr,"-coeff") == 0 ||
        strcmp(paraStr, "-coefficient") == 0 ||
        strcmp(paraStr, "-coefficients") == 0)  {
        if (OPS_GetDoubleInput(&numData,coeffData) != 0)  {
            opserr << "WARNING invalid coefficients\n";
            opserr << "uniaxialMaterial ElasticPowerFunc: " << tag[0] << endln;
            return 0;
        }
    } else  {
        opserr << "WARNING expecting -coeff but got " << paraStr << endln;
        opserr << "uniaxialMaterial ElasticPowerFunc: " << tag[0] << endln;
        return 0;
    }
    Vector coefficients(coeffData,numData);
    
    // get exponent values
    paraStr = OPS_GetString();
    if (strcmp(paraStr,"-exp") == 0 ||
        strcmp(paraStr, "-exponent") == 0 ||
        strcmp(paraStr, "-exponents") == 0)  {
        if (OPS_GetDoubleInput(&numData, expData) != 0)  {
            opserr << "WARNING invalid exponents\n";
            opserr << "uniaxialMaterial ElasticPowerFunc: " << tag[0] << endln;
            return 0;
        }
    } else  {
        opserr << "WARNING expecting -exp but got " << paraStr << endln;
        opserr << "uniaxialMaterial ElasticPowerFunc: " << tag[0] << endln;
        return 0;
    }
    Vector exponents(expData,numData);
    
    // Parsing was successful, allocate the material
    theMaterial = new ElasticPowerFunc(tag[0], coefficients, exponents, eta);
    
    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type ";
        opserr << "ElasticPowerFunc\n";
        return 0;
    }

    return theMaterial;
}


ElasticPowerFunc::ElasticPowerFunc(int tag,
    const Vector &coeff, const Vector &exp, double et)
    : UniaxialMaterial(tag, MAT_TAG_ElasticPowerFunc),
    coefficients(coeff), exponents(exp), eta(et),
    numTerms(1), initTangent(0.0),
    trialStrain(0.0), trialStrainRate(0.0),
    trialStress(0.0), trialTangent(0.0)
{
    numTerms = coefficients.Size();
    if (numTerms != exponents.Size())  {
        opserr << "ElasticPowerFunc::ElasticPowerFunc() "
            << "- coefficient and exponent arrays do not have same length.\n";
        exit(-1);        
    }
    
    this->revertToStart();
    
    initTangent = trialTangent;
}


ElasticPowerFunc::ElasticPowerFunc()
    : UniaxialMaterial(0 ,MAT_TAG_ElasticPowerFunc),
    coefficients(1), exponents(1), eta(0.0),
    numTerms(1), initTangent(0.0),
    trialStrain(0.0), trialStrainRate(0.0),
    trialStress(0.0), trialTangent(0.0)
{
    // does nothing
}


ElasticPowerFunc::~ElasticPowerFunc()
{
    // destructor does nothing
}


int ElasticPowerFunc::setTrialStrain(double strain, double strainRate)
{
    trialStrain = strain;
    trialStrainRate = strainRate;
    
    // get the stress and tangent
    trialStress = 0.0;
    trialTangent = 0.0;
    for (int i = 0; i < numTerms; i++)  {
        if (exponents(i) >= 0.0 || trialStrain != 0.0)
            trialStress += coefficients(i) * sgn(trialStrain) * pow(fabs(trialStrain), exponents(i));
        else
            trialStress += 0.0;

        if (exponents(i) >= 1.0 || trialStrain != 0.0)
            trialTangent += coefficients(i) * exponents(i) * pow(fabs(trialStrain), exponents(i) - 1.0);
        else
            trialTangent += coefficients(i) * pow(DBL_EPSILON, exponents(i) - 1.0);
    }
    trialStress += eta*trialStrainRate;
    
    return 0;
}


int ElasticPowerFunc::commitState()
{
    return 0;
}


int ElasticPowerFunc::revertToLastCommit()
{
    return 0;
}


int ElasticPowerFunc::revertToStart()
{
    trialStrain = 0.0;
    trialStrainRate = 0.0;
    trialStress = 0.0;
    trialTangent = 0.0;
    for (int i = 0; i < numTerms; i++) {
        if (exponents(i)>=1.0)
            trialTangent += coefficients(i) * exponents(i) * pow(0.0, exponents(i) - 1.0);
        else
            trialTangent += coefficients(i) * pow(DBL_EPSILON, exponents(i) - 1.0);
    }
    
    return 0;
}


UniaxialMaterial *ElasticPowerFunc::getCopy()
{
    ElasticPowerFunc *theCopy =
        new ElasticPowerFunc(this->getTag(), coefficients, exponents, eta);
    
    return theCopy;
}


int ElasticPowerFunc::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector data(4);
    data(0) = this->getTag();
    data(1) = numTerms;
    data(2) = initTangent;
    data(3) = eta;
    
    res = theChannel.sendVector(this->getDbTag(), cTag, data);
    res += theChannel.sendVector(this->getDbTag(), cTag, coefficients);
    res += theChannel.sendVector(this->getDbTag(), cTag, exponents);
    if (res < 0) 
        opserr << "ElasticPowerFunc::sendSelf() - failed to send data.\n";
    
    return res;
}


int ElasticPowerFunc::recvSelf(int cTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector data(4);
    res = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) 
        opserr << "ElasticPowerFunc::recvSelf() - failed to recv data.\n";
    else {
        this->setTag((int)data(0));
        numTerms    = (int)data(1);
        initTangent = data(2);
        eta         = data(3);
        
        // receive the strain and stress arrays
        coefficients.resize(numTerms);
        exponents.resize(numTerms);
        res += theChannel.recvVector(this->getDbTag(), cTag, coefficients);
        res += theChannel.recvVector(this->getDbTag(), cTag, exponents);
        if (res < 0) 
            opserr << "ElasticPowerFunc::recvSelf() - failed to recv arrays.\n";
    }

    return res;
}


void ElasticPowerFunc::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "ElasticPowerFunc tag: " << this->getTag() << endln;
		s << "Input Parameter: coefficients: " << coefficients << endln;
		s << "Input Parameter: exponents: " << exponents << endln;
		s << "Input Parameter: eta: " << eta << endln;
		s << "Current State: strain: " << trialStrain << " stress: ";
		s << trialStress << " tangent: " << trialTangent << endln;
	}
    
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"ElasticPowerFunc\", ";
		s << "\"coefficients\": [";
		int numPts = coefficients.Size();
		for (int i = 0; i < numPts-1; i++)
			s << coefficients(i) << ", ";
		s << coefficients(numPts - 1) << "], ";
		s << "\"exponents\": [";
		numPts = exponents.Size();
		for (int i = 0; i < numPts-1; i++)
			s << exponents(i) << ", ";
		s << exponents(numPts - 1) << "], ";
		s << "\"eta\": " << eta << "}";
	}
}


double ElasticPowerFunc::sgn(double x)
{
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}
