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
                                                                       
// Written: UW Computational Geomechanics Group
//          Pedro Arduino (*), ALborz Ghofrani(*)
//			(*)  University of Washington
//          March 2020
//
// Description: This file contains the implementation of the ManzariDafaliasPlaneStrain class.

#include "J2CyclicBoundingSurfacePlaneStrain.h"

//static vectors and matrices
Vector J2CyclicBoundingSurfacePlaneStrain::strain(3);
Vector J2CyclicBoundingSurfacePlaneStrain::stress(3);
Matrix J2CyclicBoundingSurfacePlaneStrain::tangent(3, 3);

// full constructor
J2CyclicBoundingSurfacePlaneStrain::J2CyclicBoundingSurfacePlaneStrain(int tag, double G, double K, double su, double rho, double h, double m, double h0, double chi, double beta)
:J2CyclicBoundingSurface(tag, ND_TAG_J2CyclicBoundingSurfacePlaneStrain, G, K, su, rho, h, m, h0, chi, beta)
{
}

// null constructor
J2CyclicBoundingSurfacePlaneStrain::J2CyclicBoundingSurfacePlaneStrain()
  : J2CyclicBoundingSurface()
{  
}

// destructor
J2CyclicBoundingSurfacePlaneStrain::~J2CyclicBoundingSurfacePlaneStrain()
{ 
} 

// make a clone of this material
NDMaterial* 
J2CyclicBoundingSurfacePlaneStrain::getCopy()
{ 
	J2CyclicBoundingSurfacePlaneStrain *clone;
    clone = new J2CyclicBoundingSurfacePlaneStrain();
    *clone = *this;
    return clone;
}

// send back type of material
const char* 
J2CyclicBoundingSurfacePlaneStrain::getType() const
{
    return "PlaneStrain";
}

// send back order of strain
int 
J2CyclicBoundingSurfacePlaneStrain::getOrder() const
{ 
    return 3; 
} 

// get the strain and integrate plasticity equations
int 
J2CyclicBoundingSurfacePlaneStrain::setTrialStrain(const Vector &strain_from_element) 
{
	m_strain_np1.Zero();
	m_strain_np1(0) = strain_from_element(0);
	m_strain_np1(1) = strain_from_element(1);
	m_strain_np1(3) = strain_from_element(2);

    this->integrate();
	
	return 0;
}

// unused trial strain functions
int 
J2CyclicBoundingSurfacePlaneStrain::setTrialStrain(const Vector &v, const Vector &r)
{
	return this->setTrialStrain(v);
}

// send back the strain
const Vector& 
J2CyclicBoundingSurfacePlaneStrain::getStrain()
{
	strain(0) = m_strain_np1(0);
	strain(1) = m_strain_np1(1);
	strain(2) = m_strain_np1(3);

	return strain;
} 

// send back the stress 
const Vector& 
J2CyclicBoundingSurfacePlaneStrain::getStress()
{
	stress(0) = m_stress_t_n1(0);
	stress(1) = m_stress_t_n1(1);
	stress(2) = m_stress_t_n1(3);

	return stress;
}

// send back the tangent 
const Matrix& 
J2CyclicBoundingSurfacePlaneStrain::getTangent()
{
	Matrix C(6,6);
	C = J2CyclicBoundingSurface::getTangent();

	tangent(0,0) = C(0,0);
	tangent(0,1) = C(0,1);
	tangent(0,2) = C(0,3);
	tangent(1,0) = C(1,0);
	tangent(1,1) = C(1,1);
	tangent(1,2) = C(1,3);
	tangent(2,0) = C(3,0);
	tangent(2,1) = C(3,1);
	tangent(2,2) = C(3,3);

    return tangent;
} 

// send back the tangent 
const Matrix& 
J2CyclicBoundingSurfacePlaneStrain::getInitialTangent()
{
	Matrix C(6, 6);
	J2CyclicBoundingSurface::calcInitialTangent();
	C = m_Ce;

	tangent(0,0) = C(0,0);
	tangent(0,1) = C(0,1);
	tangent(0,2) = C(0,3);
	tangent(1,0) = C(1,0);
	tangent(1,1) = C(1,1);
	tangent(1,2) = C(1,3);
	tangent(2,0) = C(3,0);
	tangent(2,1) = C(3,1);
	tangent(2,2) = C(3,3);
	
    return tangent;
}
