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
// Description: This file contains the implementation of the ManzariDafaliasPlaneStrainRO class.

#include "ManzariDafaliasPlaneStrainRO.h"

Vector ManzariDafaliasPlaneStrainRO::mEpsilon_M(3);
Vector ManzariDafaliasPlaneStrainRO::mSigma_M(3);
Vector ManzariDafaliasPlaneStrainRO::rSigma(4);
Matrix ManzariDafaliasPlaneStrainRO::mTangent(3,3);
Matrix ManzariDafaliasPlaneStrainRO::mTangent_init(3,3);

// full constructor
ManzariDafaliasPlaneStrainRO::ManzariDafaliasPlaneStrainRO(int tag, double G0, double nu, double B, double a1, double gamma1, double e_init, double Mc, double c, 
					double lambda_c, double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd, 
					double z_max, double cz, double mDen, double kappa, int integrationScheme, int tangentType, int JacoType, double TolF, double TolR)
:ManzariDafaliasRO(tag, ND_TAG_ManzariDafaliasPlaneStrainRO, G0, nu, B, a1, gamma1, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen,
				kappa, integrationScheme, tangentType, JacoType, TolF, TolR)
{
}

// null constructor
ManzariDafaliasPlaneStrainRO::ManzariDafaliasPlaneStrainRO() 
  : ManzariDafaliasRO()
{  
}

// destructor
ManzariDafaliasPlaneStrainRO::~ManzariDafaliasPlaneStrainRO() 
{ 
} 

// make a clone of this material
NDMaterial* 
ManzariDafaliasPlaneStrainRO::getCopy() 
{ 
    ManzariDafaliasPlaneStrainRO  *clone;
    clone = new ManzariDafaliasPlaneStrainRO();
    *clone = *this;
    return clone;
}

// send back type of material
const char* 
ManzariDafaliasPlaneStrainRO::getType() const 
{
    return "PlaneStrain";
}

// send back order of strain
int 
ManzariDafaliasPlaneStrainRO::getOrder() const 
{ 
    return 3; 
} 

// get the strain and integrate plasticity equations
int 
ManzariDafaliasPlaneStrainRO::setTrialStrain(const Vector &strain_from_element) 
{
	mEpsilon.Zero();
	mEpsilon(0) = -1.0 * strain_from_element(0); // -1.0 is for geotechnical sign convention
	mEpsilon(1) = -1.0 * strain_from_element(1);
	mEpsilon(3) = -1.0 * strain_from_element(2);

    this->integrate();
	
	return 0;
}

// unused trial strain functions
int 
ManzariDafaliasPlaneStrainRO::setTrialStrain(const Vector &v, const Vector &r)
{
    return this->setTrialStrain(v);
}

// send back the strain
const Vector& 
ManzariDafaliasPlaneStrainRO::getStrain() 
{
	mEpsilon_M(0) = -1.0 * mEpsilon(0); // -1.0 is for geotechnical sign convention
	mEpsilon_M(1) = -1.0 * mEpsilon(1);
	mEpsilon_M(2) = -1.0 * mEpsilon(3);
	
    return mEpsilon_M;
} 

// send back the stress 
const Vector& 
ManzariDafaliasPlaneStrainRO::getStress() 
{
	mSigma_M(0) = -1.0 * mSigma(0); // -1.0 is for geotechnical sign convention
	mSigma_M(1) = -1.0 * mSigma(1);
	mSigma_M(2) = -1.0 * mSigma(3);
	
 	return mSigma_M;
}


const Vector&
ManzariDafaliasPlaneStrainRO::getStressToRecord()
{
	rSigma(0) = -1.0 * mSigma(0); // -1.0 is for geotechnical sign convention
	rSigma(1) = -1.0 * mSigma(1);
	rSigma(2) = -1.0 * mSigma(2);
	rSigma(3) = -1.0 * mSigma(3);

	return rSigma;
}

// send back the tangent 
const Matrix& 
ManzariDafaliasPlaneStrainRO::getTangent() 
{
	Matrix C(6,6);
	if (mTangType == 0)
		C = mCe;
	else if (mTangType == 1)
		C = mCep;
	else
		C = mCep_Consistent;

	mTangent(0,0) = C(0,0);
	mTangent(0,1) = C(0,1);
	mTangent(0,2) = C(0,3);
	mTangent(1,0) = C(1,0);
	mTangent(1,1) = C(1,1);
	mTangent(1,2) = C(1,3);
	mTangent(2,0) = C(3,0);
	mTangent(2,1) = C(3,1);
	mTangent(2,2) = C(3,3);

    return mTangent;
} 

// send back the tangent 
const Matrix& 
ManzariDafaliasPlaneStrainRO::getInitialTangent() 
{
	mTangent_init(0,0) = mCe(0,0);
	mTangent_init(0,1) = mCe(0,1);
	mTangent_init(0,2) = mCe(0,3);
	mTangent_init(1,0) = mCe(1,0);
	mTangent_init(1,1) = mCe(1,1);
	mTangent_init(1,2) = mCe(1,3);
	mTangent_init(2,0) = mCe(3,0);
	mTangent_init(2,1) = mCe(3,1);
	mTangent_init(2,2) = mCe(3,3);
	
    return mTangent_init;
}

//send back the state parameters
const Vector 
ManzariDafaliasPlaneStrainRO::getState()
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
