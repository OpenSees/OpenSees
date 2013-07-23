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

// full constructor
ManzariDafalias3D::ManzariDafalias3D(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c,
	double e0, double ksi, double P_atm, double m, double h0, double ch, double nb, double A0, double nd,double z_max, 
	double cz, double massDen, double TolF, double TolR, int jacoType, int integrationScheme)
:ManzariDafalias(tag, ND_TAG_ManzariDafalias3D, G0, nu, e_init, Mc, c, lambda_c,
e0, ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, massDen, TolF, TolR, jacoType, integrationScheme)
{
}


// null constructor
ManzariDafalias3D::ManzariDafalias3D() 
  : ManzariDafalias()
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

    this->plastic_integrator();

	return 0 ;
}

// unused trial strain functions
int 
ManzariDafalias3D::setTrialStrain(const Vector &v, const Vector &r)
{
    opserr << "YOU SHOULD NOT SEE THIS: ManzariDafalias::setTrialStrain (const Vector &v, const Vector &r)" << endln;
    return this->setTrialStrain(v);
}

// send back the strain
const Vector& 
ManzariDafalias3D::getStrain() 
{
	mEpsilon_M = -1.0 * mEpsilon;
    return mEpsilon_M; // -1.0 is for geotechnical sign convention
} 

// send back the stress 
const Vector& 
ManzariDafalias3D::getStress() 
{
	mSigma_M = -1.0 * mSigma;
 	return mSigma_M; // -1.0 is for geotechnical sign convention
}

// send back the tangent 
const Matrix& 
ManzariDafalias3D::getTangent() 
{
    return mCep;
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