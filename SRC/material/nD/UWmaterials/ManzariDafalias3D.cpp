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
// Description: This file contains the implementation of the ManzariDafalias3D class.

#include "ManzariDafalias3D.h"

Vector ManzariDafalias3D::mEpsilon_M(6);
Vector ManzariDafalias3D::mSigma_M(6);

// full constructor
ManzariDafalias3D::ManzariDafalias3D(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double z_max, double cz, double mDen, int integrationScheme, 
	int tangentType, int JacoType, double TolF, double TolR)
:ManzariDafalias(tag, ND_TAG_ManzariDafalias3D, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen,
				integrationScheme, tangentType, JacoType, TolF, TolR)
{
}


// null constructor
ManzariDafalias3D::ManzariDafalias3D() 
  : ManzariDafalias(ND_TAG_ManzariDafalias3D)
{  

}

// destructor
ManzariDafalias3D::~ManzariDafalias3D() 
{ 
} 

// make a clone of this material
NDMaterial* 
ManzariDafalias3D::getCopy() 
{ 
    ManzariDafalias3D  *clone;
    clone = new ManzariDafalias3D();
    *clone = *this;
    return clone;
}

// send back type of material
const char* 
ManzariDafalias3D::getType() const 
{
    return "ThreeDimensional";
}

// send back order of strain
int 
ManzariDafalias3D::getOrder() const 
{ 
    return 6; 
} 

// get the strain and integrate plasticity equations
int 
ManzariDafalias3D::setTrialStrain(const Vector &strain_from_element) 
{
	mEpsilon = -1.0 * strain_from_element; // -1.0 is for geotechnical sign convention
	this->integrate();

	return 0 ;
}

// unused trial strain functions
int 
ManzariDafalias3D::setTrialStrain(const Vector &v, const Vector &r)
{
    return this->setTrialStrain(v);
}

// send back the strain
const Vector& 
ManzariDafalias3D::getStrain() 
{
	mEpsilon_M = -1.0 * mEpsilon;
	return mEpsilon_M; // -1.0 is for geotechnical sign convention
} 

// send back the strain
const Vector& 
ManzariDafalias3D::getEStrain() 
{
	mEpsilon_M = -1.0 * mEpsilonE;
	return mEpsilon_M; // -1.0 is for geotechnical sign convention
} 

const Vector& 
ManzariDafalias3D::getPStrain() 
{
	mEpsilon_M = -1.0 * (mEpsilon - mEpsilonE);
	return mEpsilon_M; // -1.0 is for geotechnical sign convention
} 


// send back the stress 
const Vector& 
ManzariDafalias3D::getStress() 
{
	// this->integrate();
	mSigma_M = -1.0 * mSigma;
 	return mSigma_M; // -1.0 is for geotechnical sign convention
}

const Vector&
ManzariDafalias3D::getStressToRecord()
{
	return getStress();
}

// send back the tangent 
const Matrix& 
ManzariDafalias3D::getTangent() 
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
ManzariDafalias3D::getInitialTangent() 
{
    return mCe;
}

//send back the state parameters
const Vector 
ManzariDafalias3D::getState()
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
