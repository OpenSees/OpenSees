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
// Description: This file contains the implementation of the ManzariDafaliasPlaneStrain class.

#include "ManzariDafaliasPlaneStrain.h"

//static vectors and matrices
Vector ManzariDafaliasPlaneStrain::strain(3);
Vector ManzariDafaliasPlaneStrain::stress(3);
Matrix ManzariDafaliasPlaneStrain::tangent(3,3);

// full constructor
ManzariDafaliasPlaneStrain::ManzariDafaliasPlaneStrain(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c,
	double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd,double z_max, double cz, double massDen,
	double TolF, double TolR, int jacoType, int integrationScheme)
:ManzariDafalias(tag, ND_TAG_ManzariDafalias3D, G0, nu, e_init, Mc, c, lambda_c,
	e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, massDen, TolF, TolR, jacoType, integrationScheme)
{
}

// null constructor
ManzariDafaliasPlaneStrain::ManzariDafaliasPlaneStrain() 
  : ManzariDafalias()
{  
}

// destructor
ManzariDafaliasPlaneStrain::~ManzariDafaliasPlaneStrain() 
{ 
} 

// make a clone of this material
NDMaterial* 
ManzariDafaliasPlaneStrain::getCopy() 
{ 
    ManzariDafaliasPlaneStrain  *clone;
    clone = new ManzariDafaliasPlaneStrain();
    *clone = *this;
    return clone;
}

// send back type of material
const char* 
ManzariDafaliasPlaneStrain::getType() const 
{
    return "PlaneStrain";
}

// send back order of strain
int 
ManzariDafaliasPlaneStrain::getOrder() const 
{ 
    return 3; 
} 

// get the strain and integrate plasticity equations
int 
ManzariDafaliasPlaneStrain::setTrialStrain(const Vector &strain_from_element) 
{
	mEpsilon.Zero();
	mEpsilon(0) = -1.0 * strain_from_element(0); // -1.0 is for geotechnical sign convention
	mEpsilon(1) = -1.0 * strain_from_element(1);
	mEpsilon(3) = -1.0 * strain_from_element(2);

    this->plastic_integrator();
	
	return 0;
}

// unused trial strain functions
int 
ManzariDafaliasPlaneStrain::setTrialStrain(const Vector &v, const Vector &r)
{
    return this->setTrialStrain(v);
}

// send back the strain
const Vector& 
ManzariDafaliasPlaneStrain::getStrain() 
{
	strain(0) = -1.0 * mEpsilon(0); // -1.0 is for geotechnical sign convention
	strain(1) = -1.0 * mEpsilon(1);
	strain(2) = -1.0 * mEpsilon(3);
	
    return strain;
} 

// send back the stress 
const Vector& 
ManzariDafaliasPlaneStrain::getStress() 
{
	stress(0) = -1.0 * mSigma(0); // -1.0 is for geotechnical sign convention
	stress(1) = -1.0 * mSigma(1);
	stress(2) = -1.0 * mSigma(3);
	
 	return stress;
}

// send back the tangent 
const Matrix& 
ManzariDafaliasPlaneStrain::getTangent() 
{
	tangent(0,0) = mCep(0,0);
	tangent(0,1) = mCep(0,1);
	tangent(0,2) = mCep(0,3);
	tangent(1,0) = mCep(1,0);
	tangent(1,1) = mCep(1,1);
	tangent(1,2) = mCep(1,3);
	tangent(2,0) = mCep(3,0);
	tangent(2,1) = mCep(3,1);
	tangent(2,2) = mCep(3,3);
	
    return tangent;
} 

// send back the tangent 
const Matrix& 
ManzariDafaliasPlaneStrain::getInitialTangent() 
{
	tangent(0,0) = mCe(0,0);
	tangent(0,1) = mCe(0,1);
	tangent(0,2) = mCe(0,3);
	tangent(1,0) = mCe(1,0);
	tangent(1,1) = mCe(1,1);
	tangent(1,2) = mCe(1,3);
	tangent(2,0) = mCe(3,0);
	tangent(2,1) = mCe(3,1);
	tangent(2,2) = mCe(3,3);
	
    return tangent;
}

//send back the state parameters
const Vector 
ManzariDafaliasPlaneStrain::getState()
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
