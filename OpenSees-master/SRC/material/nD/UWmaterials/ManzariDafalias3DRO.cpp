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
                                                                       
// Created: Pedro Arduino, UW, 11.2011
//
// Description: This file contains the implementation of the ManzariDafalias3DRO class.

#include "ManzariDafalias3DRO.h"

Vector ManzariDafalias3DRO::mEpsilon_M(6);
Vector ManzariDafalias3DRO::mSigma_M(6);

// full constructor
ManzariDafalias3DRO::ManzariDafalias3DRO(int tag, double G0, double nu, double B, double a1, double gamma1, double e_init, double Mc, double c, 
					double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd, 
					double z_max, double cz, double mDen, double kappa, int integrationScheme, int tangentType, int JacoType, double TolF, double TolR)
:ManzariDafaliasRO(tag, ND_TAG_ManzariDafalias3DRO, G0, nu, B, a1, gamma1, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen,
				kappa, integrationScheme, tangentType, JacoType, TolF, TolR)
{
}


// null constructor
ManzariDafalias3DRO::ManzariDafalias3DRO() 
  : ManzariDafaliasRO()
{  
}

// destructor
ManzariDafalias3DRO::~ManzariDafalias3DRO() 
{ 
} 

// make a clone of this material
NDMaterial* 
ManzariDafalias3DRO::getCopy() 
{ 
    ManzariDafalias3DRO  *clone;
    clone = new ManzariDafalias3DRO();
    *clone = *this;
    return clone;
}

// send back type of material
const char* 
ManzariDafalias3DRO::getType() const 
{
    return "ThreeDimensional";
}

// send back order of strain
int 
ManzariDafalias3DRO::getOrder() const 
{ 
    return 6; 
} 

// get the strain and integrate plasticity equations
int 
ManzariDafalias3DRO::setTrialStrain(const Vector &strain_from_element) 
{
	mEpsilon = -1.0 * strain_from_element; // -1.0 is for geotechnical sign convention

	this->integrate();

	return 0 ;
}

// unused trial strain functions
int 
ManzariDafalias3DRO::setTrialStrain(const Vector &v, const Vector &r)
{
    return this->setTrialStrain(v);
}

// send back the strain
const Vector& 
ManzariDafalias3DRO::getStrain() 
{
	mEpsilon_M = -1.0 * mEpsilon;
	return mEpsilon_M; // -1.0 is for geotechnical sign convention
} 

// send back the stress 
const Vector& 
ManzariDafalias3DRO::getStress() 
{
	mSigma_M = -1.0 * mSigma;
 	return mSigma_M; // -1.0 is for geotechnical sign convention
}

const Vector&
ManzariDafalias3DRO::getStressToRecord()
{
	return getStress();
}

// send back the tangent 
const Matrix& 
ManzariDafalias3DRO::getTangent() 
{
    if (mTangType == 0)
		return mCe;
	else if (mTangType == 1)
		return mCep;
	else
		return mCep_Consistent;
} 

// send back the tangent 
const Matrix& 
ManzariDafalias3DRO::getInitialTangent() 
{
    return mCe;
}

//send back the state parameters
const Vector 
ManzariDafalias3DRO::getState()
 {
	 Vector result(26);
	 result.Assemble(mEpsilonE,0,1.0);
	 result.Assemble(mAlpha,6,1.0);
	 result.Assemble(mFabric,12,1.0);
	 result.Assemble(mAlpha_in,18,1.0);
	 result(24) = mVoidRatio;
	 result(25) = mDGamma;

	 return result;
 }
