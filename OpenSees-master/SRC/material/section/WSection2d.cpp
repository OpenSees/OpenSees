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
                                                                        
// $Revision: 1.2 $
// $Date: 2010-02-04 19:10:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/WSection2d.cpp,v $
                                                                        
// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// WSection2d.h. WSection2d provides the abstraction of a 
// rectangular section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
// The fiber stresses are the 11, 12, and 13 components of stress, from
// which all six beam stress resultants are obtained.

#include <stdlib.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <classTags.h>
#include <WSection2d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <NDMaterial.h>

ID WSection2d::code(6);
Vector WSection2d::s(6);
Matrix WSection2d::ks(6,6);

// constructors:
WSection2d::WSection2d(int tag, NDMaterial &theMat,
		       double D, double Tw, double Bf, double Tf,
		       int Nfdw, int Nftf, double shape, double flag):
  SectionForceDeformation(tag, SEC_TAG_WSection2d),
  theFibers(0), yFibers(0), AFibers(0), e(6),
  d(D), tw(Tw), bf(Bf), tf(Tf), nfdw(Nfdw), nftf(Nftf), shapeFactor(shape)
{
  int numFibers = nfdw + 2*nftf;
  
  theFibers = new NDMaterial*[numFibers];
  yFibers   = new double[numFibers];
  AFibers   = new double[numFibers];
  
  for (int i = 0; i < numFibers; i++) {
    theFibers[i] = flag ? theMat.getCopy("BeamFiber") :
      theMat.getCopy("TimoshenkoFiber");
    if (theFibers[i] == 0)
      opserr << "WSection2d::WSection2d -- failed to get copy of beam fiber" << endln;
  }
  
  double dw = d-2*tf;
  
  double a_f = bf*tf/nftf;
  double a_w = dw*tw/nfdw;
  
  int loc = 0;
  
  double yIncr  = tf/nftf;
  double yStart = 0.5*d - 0.5*yIncr;
  
  for (loc = 0; loc < nftf; loc++) {
    AFibers[loc] = AFibers[numFibers-loc-1] = a_f;
    yFibers[loc] = yStart - yIncr*loc;
    yFibers[numFibers-loc-1] = -yFibers[loc];
  }
  
  yIncr  = dw/nfdw;
  yStart = 0.5*dw - 0.5*yIncr;
  
  int count = 0;
  
  for ( ; loc < numFibers-nftf; loc++, count++) {
    AFibers[loc] = a_w;
    yFibers[loc] = yStart - yIncr*count;
  }
  
  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_VY;
  code(4) = SECTION_RESPONSE_VZ;
  code(5) = SECTION_RESPONSE_T;
}

// constructor for blank object that recvSelf needs to be invoked upon
WSection2d::WSection2d():
  SectionForceDeformation(0, SEC_TAG_WSection2d),
  theFibers(0), yFibers(0), AFibers(0), e(6),
  d(0.0), tw(0.0), bf(0.0), tf(0.0),
  nfdw(0), nftf(0), shapeFactor(0.0)
{
  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_VY;
  code(4) = SECTION_RESPONSE_VZ;
  code(5) = SECTION_RESPONSE_T;
}

// destructor:
WSection2d::~WSection2d()
{
  int numFibers = nfdw + 2*nftf;
  
  if (theFibers != 0) {
    for (int i = 0; i < numFibers; i++)
      if (theFibers[i] != 0)
	delete theFibers[i];
    delete [] theFibers;
  }

  if (yFibers != 0)
    delete [] yFibers;

  if (AFibers != 0)
    delete [] AFibers;
}

// Compute fiber strain, eps, from section deformation, e, using the
// linear kinematic relationship, eps = a*e, where
// eps = [eps_11, eps_12, eps_13]'
// e   = [eps_a, kappa_z, kappa_y, gamma_y, gamma_z, phi]'
// a   = [1 -y z         0         0  0
//        0  0 0 sqrt(5/6)         0 -z
//        0  0 0         0 sqrt(5/6)  y]
int WSection2d::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  
  e = deforms;
  
  // Material strain vector
  static Vector eps(3);
  
  double y, z;
  
  double root56 = sqrt(shapeFactor);

  int numFibers = nfdw + 2*nftf;
  
  for (int i = 0; i < numFibers; i++) {
    y = yFibers[i];
    z = 0.0;
      
    eps(0) = e(0) - y*e(1) + z*e(2);
    eps(1) = root56*e(3) - z*e(5);
    eps(2) = root56*e(4) + y*e(5);
    
    res += theFibers[i]->setTrialStrain(eps);
  }

  return res;
}

const Vector&
WSection2d::getSectionDeformation(void)
{
  return e;
}

// Compute section tangent stiffness, ks, from material tangent, Dt,
// by integrating over section area
// ks = int_A a'*Dt*a dA
// a   = [1 -y z         0         0  0
//        0  0 0 sqrt(5/6)         0 -z
//        0  0 0         0 sqrt(5/6)  y]
const Matrix&
WSection2d::getSectionTangent(void)
{
  ks.Zero();
 
  double y, z, w;
  double y2, z2, yz;
  
  double d00, d01, d02;
  double d10, d11, d12;
  double d20, d21, d22;
  
  double tmp;
  
  double five6  = shapeFactor;
  double root56 = sqrt(shapeFactor);

  int numFibers = nfdw + 2*nftf;
  
  for (int i = 0; i < numFibers; i++) {
    
    y = yFibers[i];
    z = 0.0;
    w = AFibers[i];

    y2 = y*y;
    z2 = z*z;
    yz = y*z;
    
    const Matrix &Dt = theFibers[i]->getTangent();

    d00 = Dt(0,0)*w; d01 = Dt(0,1)*w; d02 = Dt(0,2)*w;
    d10 = Dt(1,0)*w; d11 = Dt(1,1)*w; d12 = Dt(1,2)*w;
    d20 = Dt(2,0)*w; d21 = Dt(2,1)*w; d22 = Dt(2,2)*w;
    
    // Bending terms
    ks(0,0) += d00;
    ks(1,1) += y2*d00;
    ks(2,2) += z2*d00;
    tmp = -y*d00;
    ks(0,1) += tmp;
    ks(1,0) += tmp;
    tmp = z*d00;
    ks(0,2) += tmp;
    ks(2,0) += tmp;
    tmp = -yz*d00;
    ks(1,2) += tmp;
    ks(2,1) += tmp;
    
    // Shear terms
    ks(3,3) += five6*d11;
    ks(3,4) += five6*d12;
    ks(4,3) += five6*d21;
    ks(4,4) += five6*d22;
    
    // Torsion term
    ks(5,5) += z2*d11 - yz*(d12+d21) + y2*d22;
    
    // Bending-torsion coupling terms
    tmp = -z*d01 + y*d02;
    ks(0,5) += tmp;
    ks(1,5) -= y*tmp;
    ks(2,5) += z*tmp;
    tmp = -z*d10 + y*d20;
    ks(5,0) += tmp;
    ks(5,1) -= y*tmp;
    ks(5,2) += z*tmp;
    
    // Hit tangent terms with root56
    d01 *= root56; d02 *= root56;
    d10 *= root56; d11 *= root56; d12 *= root56;
    d20 *= root56; d21 *= root56; d22 *= root56;
    
    // Bending-shear coupling terms
    ks(0,3) += d01;
    ks(0,4) += d02;
    ks(1,3) -= y*d01;
    ks(1,4) -= y*d02;
    ks(2,3) += z*d01;
    ks(2,4) += z*d02;
    ks(3,0) += d10;
    ks(4,0) += d20;
    ks(3,1) -= y*d10;
    ks(4,1) -= y*d20;
    ks(3,2) += z*d10;
    ks(4,2) += z*d20;
    
    // Torsion-shear coupling terms
    y2 =  y*d22;
    z2 = -z*d11;
    ks(5,3) +=  z2 + y*d21;
    ks(5,4) += -z*d12 + y2;
    ks(3,5) +=  z2 + y*d12;
    ks(4,5) += -z*d21 + y2;
  }

  // Non-zero value since this is 2d section (so ks can be inverted)
  ks(2,2) = 1.0;

  return ks;
}

// Compute section stress resultant, s, from material stress, sig,
// by integrating over section area
// s = int_A a'*sig dA
// s = [P, Mz, My, Vy, Vz, T]'
// a = [1 -y z         0         0  0
//      0  0 0 sqrt(5/6)         0 -z
//      0  0 0         0 sqrt(5/6)  y]
const Vector&
WSection2d::getStressResultant(void)
{
  s.Zero();
  
  double y, z, w;
  double sig0, sig1, sig2;
  
  double root56 = sqrt(shapeFactor);

  int numFibers = nfdw + 2*nftf;
  
  for (int i = 0; i < numFibers; i++) {
    
    y = yFibers[i];
    z = 0.0;
    w = AFibers[i];
    
    const Vector &sig = theFibers[i]->getStress();
    
    sig0 = sig(0)*w;
    sig1 = sig(1)*w;
    sig2 = sig(2)*w;
    
    s(0) += sig0;
    s(1) -= y*sig0;
    s(2) += z*sig0;
    s(3) += root56*sig1;
    s(4) += root56*sig2;
    s(5) += -z*sig1 + y*sig2;
  }

  return s;
}

SectionForceDeformation*
WSection2d::getCopy(void)
{
  WSection2d *theCopy =
    new WSection2d (this->getTag(), *theFibers[0],
		    d, tw, bf, tf, nfdw, nftf, shapeFactor);
  
  theCopy->e = e;
  
  return theCopy;
}

const ID&
WSection2d::getType(void)
{
  return code;
}

int
WSection2d::getOrder(void) const
{
  return 6;
}

int
WSection2d::commitState(void)
{
  int err = 0;
  
  int numFibers = nfdw + 2*nftf;

  for (int i = 0; i < numFibers; i++)
    err += theFibers[i]->commitState();
  
  return err;
}

int
WSection2d::revertToLastCommit(void)
{
  int err = 0;
  
  int numFibers = nfdw + 2*nftf;
  
  for (int i = 0; i < numFibers; i++)
    err += theFibers[i]->revertToLastCommit();
  
  return err;
}

int
WSection2d::revertToStart(void)
{
  int err = 0;

  int numFibers = nfdw + 2*nftf;
  
  for (int i = 0; i < numFibers; i++)
    err += theFibers[i]->revertToStart();
  
  return err;
}

int
WSection2d::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
WSection2d::recvSelf(int commitTag, Channel &theChannel,
		     FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
WSection2d::Print(OPS_Stream &s, int flag)
{
  s << "\nWSection2d, tag: " << this->getTag() << endln;
  s << "\tFiber Material, tag: " << theFibers[0]->getTag() << endln;
  s << "\tSection depth:    " << d << endln;
  s << "\tWeb thickness:    " << tw << endln;
  s << "\tFlange width:     " << bf << endln;
  s << "\tFlange thickness: " << tf << endln;
  s << "\tShape factor:     " << shapeFactor << endln;
  theFibers[0]->Print(s, flag);
}
