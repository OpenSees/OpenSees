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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-04-02 22:02:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/quadrule/GaussQuadRule1d.cpp,v $
                                                                        
// written: rms
// Created: 12/98
//
// Description: This file contains the implementation of 
// GaussQuadRule1d (Quadrature Rule,0).

#include <stdlib.h>
#include <Vector.h>
#include <Matrix.h>

#include <GaussQuadRule1d.h>

bool GaussQuadRule1d::dataSet = false;

double GaussQuadRule1d::ptsArray[maxOrder*(maxOrder+1)/2];
double GaussQuadRule1d::wtsArray[maxOrder*(maxOrder+1)/2];

Matrix GaussQuadRule1d::pts1 ( ptsArray,     1,  1);
Matrix GaussQuadRule1d::pts2 (&ptsArray[1],  2,  1);
Matrix GaussQuadRule1d::pts3 (&ptsArray[3],  3,  1);
Matrix GaussQuadRule1d::pts4 (&ptsArray[6],  4,  1);
Matrix GaussQuadRule1d::pts5 (&ptsArray[10], 5,  1);
Matrix GaussQuadRule1d::pts6 (&ptsArray[15], 6,  1);
Matrix GaussQuadRule1d::pts7 (&ptsArray[21], 7,  1);
Matrix GaussQuadRule1d::pts8 (&ptsArray[28], 8,  1);
Matrix GaussQuadRule1d::pts9 (&ptsArray[36], 9,  1);
Matrix GaussQuadRule1d::pts10(&ptsArray[45], 10, 1);

Vector GaussQuadRule1d::wts1 ( wtsArray,     1);
Vector GaussQuadRule1d::wts2 (&wtsArray[1],  2);
Vector GaussQuadRule1d::wts3 (&wtsArray[3],  3);
Vector GaussQuadRule1d::wts4 (&wtsArray[6],  4);
Vector GaussQuadRule1d::wts5 (&wtsArray[10], 5);
Vector GaussQuadRule1d::wts6 (&wtsArray[15], 6);
Vector GaussQuadRule1d::wts7 (&wtsArray[21], 7);
Vector GaussQuadRule1d::wts8 (&wtsArray[28], 8);
Vector GaussQuadRule1d::wts9 (&wtsArray[36], 9);
Vector GaussQuadRule1d::wts10(&wtsArray[45], 10);

GaussQuadRule1d::GaussQuadRule1d()
  :order(0), myPts(0), myWts(0)
{
  if (dataSet == false) {
    // One point
    ptsArray[0] = 0.0;
    
    wtsArray[0] = 2.0;

    // Two points
    ptsArray[1] = -0.577350269189626;
    ptsArray[2] =  0.577350269189626;
    
    wtsArray[1] = 1.0;
    wtsArray[2] = 1.0;

    // Three points
    ptsArray[3] = -0.774596669241483;
    ptsArray[4] =  0.0;
    ptsArray[5] =  0.774596669241483;
    
    wtsArray[3] = 0.555555555555556;
    wtsArray[4] = 0.888888888888889;
    wtsArray[5] = 0.555555555555556;

    // Four points
    ptsArray[6] = -0.861136311594053;
    ptsArray[7] = -0.339981043584856;
    ptsArray[8] =  0.339981043584856;
    ptsArray[9] =  0.861136311594053;
    
    wtsArray[6] = 0.347854845137454;
    wtsArray[7] = 0.652145154862546;
    wtsArray[8] = 0.652145154862546;
    wtsArray[9] = 0.347854845137454;
      
    // Five points
    ptsArray[10] = -0.906179845938664;
    ptsArray[11] = -0.538469310105683;
    ptsArray[12] =  0.0;
    ptsArray[13] =  0.538469310105683;
    ptsArray[14] =  0.906179845938664;
    
    wtsArray[10] = 0.236926885056189;
    wtsArray[11] = 0.478628670499366;
    wtsArray[12] = 0.568888888888889;
    wtsArray[13] = 0.478628670499366;
    wtsArray[14] = 0.236926885056189;

    // Six points
    ptsArray[15] = -0.932469514203152;
    ptsArray[16] = -0.661209386466265;
    ptsArray[17] = -0.238619186083197;
    ptsArray[18] =  0.238619186083197;
    ptsArray[19] =  0.661209386466265;
    ptsArray[20] =  0.932469514203152;
    
    wtsArray[15] = 0.171324492379170;
    wtsArray[16] = 0.360761573048139;
    wtsArray[17] = 0.467913934572691;
    wtsArray[18] = 0.467913934572691;
    wtsArray[19] = 0.360761573048139;
    wtsArray[20] = 0.171324492379170;

    // Seven points
    ptsArray[21] = -0.949107912342759;
    ptsArray[22] = -0.741531185599394;
    ptsArray[23] = -0.405845151377397;
    ptsArray[24] =  0.0;
    ptsArray[25] =  0.405845151377397;
    ptsArray[26] =  0.741531185599394;
    ptsArray[27] =  0.949107912342759;
    
    wtsArray[21] = 0.129484966168870;
    wtsArray[22] = 0.279705391489277;
    wtsArray[23] = 0.381830050505119;
    wtsArray[24] = 0.417959183673469;
    wtsArray[25] = 0.381830050505119;
    wtsArray[26] = 0.279705391489277;
    wtsArray[27] = 0.129484966168870;

    // Eight points
    ptsArray[28] = -0.960289856497536;
    ptsArray[29] = -0.796666477413627;
    ptsArray[30] = -0.525532409916329;
    ptsArray[31] = -0.183434642495650;
    ptsArray[32] =  0.183434642495650;
    ptsArray[33] =  0.525532409916329;
    ptsArray[34] =  0.796666477413627;
    ptsArray[35] =  0.960289856497536;
    
    wtsArray[28] = 0.101228536290376;
    wtsArray[29] = 0.222381034453374;
    wtsArray[30] = 0.313706645877887;
    wtsArray[31] = 0.362683783378362;
    wtsArray[32] = 0.362683783378362;
    wtsArray[33] = 0.313706645877887;
    wtsArray[34] = 0.222381034453374;
    wtsArray[35] = 0.101228536290376;
      
    // Nine points
    ptsArray[36] = -0.968160239507626;
    ptsArray[37] = -0.836031107326636;
    ptsArray[38] = -0.613371432700590;
    ptsArray[39] = -0.324253423403809;
    ptsArray[40] =  0.0;
    ptsArray[41] =  0.324253423403809;
    ptsArray[42] =  0.613371432700590;
    ptsArray[43] =  0.836031107326636;
    ptsArray[44] =  0.968160239507626;
    
    wtsArray[36] = 0.081274388361574;
    wtsArray[37] = 0.180648160694857;
    wtsArray[38] = 0.260610696402935;
    wtsArray[39] = 0.312347077040003;
    wtsArray[40] = 0.330239355001260;
    wtsArray[41] = 0.312347077040003;
    wtsArray[42] = 0.260610696402935;
    wtsArray[43] = 0.180648160694857;
    wtsArray[44] = 0.081274388361574;

    // Ten points
    ptsArray[45] = -0.973906528517172;
    ptsArray[46] = -0.865063366688985;
    ptsArray[47] = -0.679409568299024;
    ptsArray[48] = -0.433395394129247;
    ptsArray[49] = -0.148874338981631;
    ptsArray[50] =  0.148874338981631;
    ptsArray[51] =  0.433395394129247;
    ptsArray[52] =  0.679409568299024;
    ptsArray[53] =  0.865063366688985;
    ptsArray[54] =  0.973906528517172;
    
    wtsArray[45] = 0.066671344308688;
    wtsArray[46] = 0.149451349150581;
    wtsArray[47] = 0.219086362515982;
    wtsArray[48] = 0.269266719309996;
    wtsArray[49] = 0.295524224714753;
    wtsArray[50] = 0.295524224714753;
    wtsArray[51] = 0.269266719309996;
    wtsArray[52] = 0.219086362515982;
    wtsArray[53] = 0.149451349150581;
    wtsArray[54] = 0.066671344308688;

    dataSet = true;
  }
}

GaussQuadRule1d::~GaussQuadRule1d()
{
  // Nothing to do
}


int GaussQuadRule1d::setOrder(int quadOrder)
{
  if (quadOrder < 1 || quadOrder > maxOrder) {
    opserr << "GaussQuadRule1d::setOrder() -- Invalid quadrature order " << quadOrder << endln;
    exit(-1);
  }
  
  // Nothing needs to change if this is true
  if (order == quadOrder)
    return 0;
  
  order = quadOrder;

  switch (order) {
  case 1:
    myPts = &pts1;
    myWts = &wts1;
    break;
    
  case 2:
    myPts = &pts2;
    myWts = &wts2;
    break;
    
  case 3:
    myPts = &pts3;
    myWts = &wts3;
    break;
    
  case 4:
    myPts = &pts4;
    myWts = &wts4;
    break;
    
  case 5:
    myPts = &pts5;
    myWts = &wts5;
    break;
    
  case 6:
    myPts = &pts6;
    myWts = &wts6;
    break;
    
  case 7:
    myPts = &pts7;
    myWts = &wts7;
    break;
    
  case 8:
    myPts = &pts8;
    myWts = &wts8;
    break;
    
  case 9:
    myPts = &pts9;
    myWts = &wts9;
    break;
    
  case 10:
    myPts = &pts10;
    myWts = &wts10;
    break;

  default:
    opserr << "GaussQuadRule1d::setOrder() -- Invalid quadrature order " << order << endln;
    return -1;
    break;
  }    

  return 0;
}

int GaussQuadRule1d::getOrder (void) const
{
  return order;
}

int GaussQuadRule1d::getNumIntegrPoints (void) const
{
  return order;
}

const Matrix & 
GaussQuadRule1d::getIntegrPointCoords (void) const
{
  if (order < 1 || order > maxOrder)
    opserr << "GaussQuadRule1d::getIntegrPointWeights() -- order " << order << " is currently invalid\n";

  return *myPts;
}

const Vector & 
GaussQuadRule1d::getIntegrPointWeights (void) const
{
  if (order < 1 || order > maxOrder)
    opserr << "GaussQuadRule1d::getIntegrPointWeights() -- order " << order << " is currently invalid\n";
			   
  return *myWts;
}

const Matrix & 
GaussQuadRule1d::getIntegrPointCoords (int quadOrder)
{
  if (order != quadOrder)
    this->setOrder(quadOrder);

  return *myPts;
}

const Vector & 
GaussQuadRule1d::getIntegrPointWeights (int quadOrder)
{
  if (order != quadOrder)
    this->setOrder(quadOrder);

  return *myWts;
}
