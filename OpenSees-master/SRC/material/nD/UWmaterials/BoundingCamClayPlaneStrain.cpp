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
//          February 2011

#include <BoundingCamClayPlaneStrain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//static vectors and matrices
Vector BoundingCamClayPlaneStrain::strain(3);
Vector BoundingCamClayPlaneStrain::stress(3);
Matrix BoundingCamClayPlaneStrain::tangent(3,3);

//null constructor
BoundingCamClayPlaneStrain::BoundingCamClayPlaneStrain() : 
BoundingCamClay( )
{  
}

//full constructor
BoundingCamClayPlaneStrain::BoundingCamClayPlaneStrain(int tag, double mDen, double C, double bulk, double OCR,
                                                       double mu_o, double alpha, double lambda, double h, double m): 
BoundingCamClay(tag, ND_TAG_BoundingCamClayPlaneStrain, mDen, C, bulk, OCR, mu_o, alpha, lambda, h, m)
{
}

//destructor
BoundingCamClayPlaneStrain::~BoundingCamClayPlaneStrain() 
{ 
} 

//make a clone of this material
NDMaterial* 
BoundingCamClayPlaneStrain::getCopy() 
{ 
    BoundingCamClayPlaneStrain  *clone;
    clone = new BoundingCamClayPlaneStrain();   //new instance of this class
    *clone = *this ;          //assignment to make copy
    return clone ;
}

//send back type of material
const char* 
BoundingCamClayPlaneStrain::getType() const 
{
    return "PlaneStrain";
}

//send back order of strain in vector form
int 
BoundingCamClayPlaneStrain::getOrder() const 
{ 
    return 3; 
} 

//get the strain and integrate plasticity equations
int 
BoundingCamClayPlaneStrain::setTrialStrain(const Vector &strain_from_element) 
{
	mEpsilon.Zero();
	mEpsilon(0) = strain_from_element(0);
	mEpsilon(1) = strain_from_element(1);
	mEpsilon(3) = strain_from_element(2);

    this->plastic_integrator();
	
	return 0;
}

//unused trial strain function
int 
BoundingCamClayPlaneStrain::setTrialStrain(const Vector &v, const Vector &r)
{
    return this->setTrialStrain (v);
}

//send back the strain
const Vector& 
BoundingCamClayPlaneStrain::getStrain() 
{
	strain(0) = mEpsilon(0);
	strain(1) = mEpsilon(1);
	strain(2) = mEpsilon(3);
	
    return strain;
} 

//send back the stress 
const Vector& 
BoundingCamClayPlaneStrain::getStress() 
{
	stress(0) = mSigma(0);
	stress(1) = mSigma(1);
	stress(2) = mSigma(3);
	
 	return stress;
}

//send back the tangent 
const Matrix& 
BoundingCamClayPlaneStrain::getTangent() 
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

//send back the tangent 
const Matrix& BoundingCamClayPlaneStrain::getInitialTangent() 
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

