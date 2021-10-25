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
// Description: This file contains the implementation of the SAniSandMS3D class.

#include "SAniSandMS3D.h"

Vector SAniSandMS3D::mEpsilon_M(6);
Vector SAniSandMS3D::mEpsilonE_M(6);
Vector SAniSandMS3D::mSigma_M(6);

// full constructor
SAniSandMS3D::SAniSandMS3D(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
				double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double zeta, double mu0, 
					double beta,
			 	 	double mDen, 
				int integrationScheme, int tangentType, int JacoType, double TolF, double TolR)
:SAniSandMS(tag, ND_TAG_SAniSandMS3D, G0, nu, e_init, Mc, c, lambda_c, e0, ksi, P_atm, m, h0, ch, nb, A0, nd, zeta, mu0, 
			beta,mDen, 
			integrationScheme, tangentType, JacoType, TolF, TolR)
{
}


// null constructor
SAniSandMS3D::SAniSandMS3D() 
  : SAniSandMS()
{  
}

// destructor
SAniSandMS3D::~SAniSandMS3D() 
{ 
} 

// make a clone of this material
NDMaterial* 
SAniSandMS3D::getCopy() 
{ 
    SAniSandMS3D  *clone;
    clone = new SAniSandMS3D();
    *clone = *this;
    return clone;
}

// send back type of material
const char* 
SAniSandMS3D::getType() const 
{
    return "ThreeDimensional";
}

// send back order of strain
int 
SAniSandMS3D::getOrder() const 
{ 
    return 6; 
} 

// get the strain and integrate plasticity equations
int 
SAniSandMS3D::setTrialStrain(const Vector &strain_from_element) 
{
    mEpsilon = -1.0 * strain_from_element; // -1.0 is for geotechnical sign convention
	this->integrate();

	return 0 ;
}

// unused trial strain functions
int 
SAniSandMS3D::setTrialStrain(const Vector &v, const Vector &r)
{
    return this->setTrialStrain(v);
}

// send back the strain
const Vector& 
SAniSandMS3D::getStrain() 
{
	mEpsilon_M = -1.0 * mEpsilon;
	return mEpsilon_M; // -1.0 is for geotechnical sign convention
} 


// send back the strain
const Vector& 
SAniSandMS3D::getEStrain() 
{
    mEpsilonE_M = -1.0 * mEpsilonE;
    return mEpsilonE_M; // -1.0 is for geotechnical sign convention
} 


// send back the stress 
const Vector& 
SAniSandMS3D::getStress() 
{
	// this->integrate();
	mSigma_M = -1.0 * mSigma;
 	return mSigma_M; // -1.0 is for geotechnical sign convention
}

const Vector&
SAniSandMS3D::getStressToRecord()
{
	return getStress();
}

// send back the tangent 
const Matrix& 
SAniSandMS3D::getTangent() 
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
SAniSandMS3D::getInitialTangent() 
{
    return mCe;
}

//send back the state parameters
const Vector 
SAniSandMS3D::getState()
 {
	 Vector result(26);
	 result.Assemble(mEpsilonE,0,1.0);
	 result.Assemble(mAlpha,6,1.0);
	 // result.Assemble(mFabric,12,1.0);
	 // result.Assemble(mAlpha_in,18,1.0);
	 result(24) = mVoidRatio;
	 result(25) = mDGamma;

	 return result;
 }
