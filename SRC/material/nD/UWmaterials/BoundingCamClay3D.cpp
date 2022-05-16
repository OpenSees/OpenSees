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
                                                                        
// Written: K.Petek	
// Created: 12/04

#include <BoundingCamClay3D.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//null constructor
BoundingCamClay3D::BoundingCamClay3D() : 
BoundingCamClay( )
{  
}

//full constructor
BoundingCamClay3D::BoundingCamClay3D(int tag, double mDen, double C, double bulk, double OCR, double mu_o,
							                  double alpha, double lambda, double h, double m) : 
BoundingCamClay(tag, ND_TAG_BoundingCamClay3D, mDen, C, bulk, OCR, mu_o, alpha, lambda, h, m)
{
}

//destructor
BoundingCamClay3D::~BoundingCamClay3D() 
{ 
} 

//make a clone of this material
NDMaterial* BoundingCamClay3D :: getCopy( ) 
{ 
    BoundingCamClay3D  *clone;
    clone = new BoundingCamClay3D( ) ;   //new instance of this class
    *clone = *this ;          //assignment to make copy
    return clone ;
}

//send back type of material
const char* BoundingCamClay3D::getType() const 
{
    return "ThreeDimensional";
}

//send back order of strain in vector form
int BoundingCamClay3D::getOrder() const 
{ 
    return 6; 
} 

//get the strain and integrate plasticity equations
int BoundingCamClay3D::setTrialStrain( const Vector &strain_from_element) 
{
	mEpsilon = strain_from_element;

    this->plastic_integrator();
	
	return 0 ;
}

//unused trial strain functions
int BoundingCamClay3D::setTrialStrain(const Vector &v, const Vector &r)
{
    opserr << "YOU SHOULD NOT SEE THIS: BoundingCamClay::setTrialStrain (const Vector &v, const Vector &r)" << endln;
    return this->setTrialStrain (v);
}

//send back the strain
const Vector& BoundingCamClay3D::getStrain() 
{
    return mEpsilon;
} 

//send back the stress 
const Vector& BoundingCamClay3D::getStress() 
{
 	return mSigma;
}

//send back the tangent 
const Matrix& BoundingCamClay3D::getTangent() 
{
    return mCep;
} 

//send back the tangent 
const Matrix& BoundingCamClay3D::getInitialTangent() 
{
    return mCep;
} 

