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

#include <ManzariDafalias3D.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

// full constructor
ManzariDafalias3D::ManzariDafalias3D(int tag, double mDen) 
  : ManzariDafalias(tag, ND_TAG_ManzariDafalias3D)
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
	mEpsilon = strain_from_element;

    this->plastic_integrator();
	
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
    return mEpsilon;
} 

// send back the stress 
const Vector& 
ManzariDafalias3D::getStress() 
{
 	return mSigma;
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
    return mCep;
}
