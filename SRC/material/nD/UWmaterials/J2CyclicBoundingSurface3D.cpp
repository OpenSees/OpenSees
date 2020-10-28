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
// Description: This file contains the implementation of the J2CyclicBoundingSurface3D class.
#include "J2CyclicBoundingSurface3D.h"

// full constructor
J2CyclicBoundingSurface3D::J2CyclicBoundingSurface3D( int tag, double G, double K, double su, double rho, double h, double m, double h0, double chi, double beta)
	:J2CyclicBoundingSurface (tag, ND_TAG_J2CyclicBoundingSurface3D, G, K, su, rho, h, m, h0, chi, beta)
{
}

// null constructor
J2CyclicBoundingSurface3D::J2CyclicBoundingSurface3D()
  : J2CyclicBoundingSurface()
{  
}

// destructor
J2CyclicBoundingSurface3D::~J2CyclicBoundingSurface3D()
{ 
} 

// make a clone of this material
NDMaterial* 
J2CyclicBoundingSurface3D::getCopy()
{ 
	J2CyclicBoundingSurface3D *clone;
    clone = new J2CyclicBoundingSurface3D();
    *clone = *this;
    return clone;
}

// send back type of material
const char* 
J2CyclicBoundingSurface3D::getType() const
{
    return "ThreeDimensional";
}

// send back order of strain
int 
J2CyclicBoundingSurface3D::getOrder() const
{ 
    return 6; 
} 

// get the strain and integrate plasticity equations
int 
J2CyclicBoundingSurface3D::setTrialStrain(const Vector &strain_from_element)
{
	m_strain_np1 = strain_from_element;
	this->integrate();

	return 0;
}

// unused trial strain functions
int 
J2CyclicBoundingSurface3D::setTrialStrain(const Vector &v, const Vector &r)
{
	m_strainRate_n1 = r;
	m_strain_np1 = v;
	this->integrate();

	return 0;
}

// send back the strain
const Vector& 
J2CyclicBoundingSurface3D::getStrain()
{
	return m_strain_np1;
} 

// send back the stress 
const Vector& 
J2CyclicBoundingSurface3D::getStress()
{
	//return m_stress_np1;
	return m_stress_t_n1;
}

// send back the tangent 
const Matrix& 
J2CyclicBoundingSurface3D::getTangent()
{
	Matrix C(6, 6);
	C = J2CyclicBoundingSurface::getTangent();
    return C;
} 

// send back the tangent 
const Matrix& 
J2CyclicBoundingSurface3D::getInitialTangent()
{
	Matrix C(6, 6);
	J2CyclicBoundingSurface::calcInitialTangent();
	C = m_Ce;
	return C;
}