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
                                                                        
// $Revision: 1.3 $
// $Date: 2010-02-04 19:10:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/WFFiberSection2d.cpp,v $
                                                                        
// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// WFFiberSection2d.h. WFFiberSection2d provides the abstraction of a 
// rectangular section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
// The fiber stresses are the 11, 12, and 13 components of stress, from
// which all six beam stress resultants are obtained.

#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <classTags.h>
#include <WFFiberSection2d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>

// constructors:
WFFiberSection2d::WFFiberSection2d(int tag, UniaxialMaterial &theMat,
				   double D, double Tw, double Bf, double Tf,
				   int Nfdw, int Nftf):
  FiberSection2d(),
  d(D), tw(Tw), bf(Bf), tf(Tf), nfdw(Nfdw), nftf(Nftf)
{
  int numFibers = nfdw + 2*nftf;
  
  for (int i = 0; i < numFibers; i++) {
    theMaterials[i] = theMat.getCopy();
    if (theMaterials[i] == 0)
      opserr << "WFFiberSection2d::WFFiberSection2d -- failed to get copy of beam fiber" << endln;
  }
  
  double dw = d-2*tf;
  
  double a_f = bf*tf/nftf;
  double a_w = dw*tw/nfdw;
  
  int loc = 0;
  
  double yIncr  = tf/nftf;
  double yStart = 0.5*d - 0.5*yIncr;
  
  double *AFibers = new double[numFibers];
  double *yFibers = new double[numFibers];

  for (loc = 0; loc < nftf; loc++) {
    AFibers[loc] = a_f;
    yFibers[loc] = yStart - yIncr*loc;
    AFibers[numFibers-loc-1] = a_f;
    yFibers[numFibers-loc-1] = -yFibers[loc];
  }
  
  yIncr  = dw/nfdw;
  yStart = 0.5*dw - 0.5*yIncr;
  
  int count = 0;
  
  for ( ; loc < numFibers-nftf; loc++, count++) {
    AFibers[loc] = a_w;
    yFibers[loc] = yStart - yIncr*count;
  }

  for (int i = 0; i < numFibers; i++) {
    matData[i*2]   = -yFibers[i];
    matData[i*2+1] =  AFibers[i];
  }

  delete [] yFibers;
  delete [] AFibers;
}

// constructor for blank object that recvSelf needs to be invoked upon
WFFiberSection2d::WFFiberSection2d():
  FiberSection2d(),
  d(0.0), tw(0.0), bf(0.0), tf(0.0), nfdw(0), nftf(0)
{

}

// destructor:
WFFiberSection2d::~WFFiberSection2d()
{

}

SectionForceDeformation*
WFFiberSection2d::getCopy(void)
{
  WFFiberSection2d *theCopy =
    new WFFiberSection2d (this->getTag(), *theMaterials[0],
			  d, tw, bf, tf, nfdw, nftf);
  
  return theCopy;
}

int
WFFiberSection2d::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
WFFiberSection2d::recvSelf(int commitTag, Channel &theChannel,
			   FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
WFFiberSection2d::Print(OPS_Stream &s, int flag)
{
  s << "\nWFFiberSection2d, tag: " << this->getTag() << endln;
  s << "\tSection depth:    " << d << endln;
  s << "\tWeb thickness:    " << tw << endln;
  s << "\tFlange width:     " << bf << endln;
  s << "\tFlange thickness: " << tf << endln;
  FiberSection2d::Print(s, flag);
}

const Vector &
WFFiberSection2d::getStressResultantSensitivity(int gradIndex,
						bool conditional)
{
  static Vector ds;
  
  // get material stress contribution
  ds = FiberSection2d::getStressResultantSensitivity(gradIndex, conditional);
  
  double y, A, stressGradient;
  int loc = 0;
  
  for (int i = 0; i < numFibers; i++) {
    y = matData[loc++];
    A = matData[loc++];
    
    stressGradient = theMaterials[i]->getStressSensitivity(gradIndex,true);
    stressGradient = stressGradient * A;
    ds(0) += stressGradient;
    ds(1) += stressGradient * y;
  }
  
  return ds;
}
