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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/quadrule/GaussQuadRule1d.cpp,v $
                                                                        
                                                                        
// File: ~/QuadRule/GaussQuadRule1d.C
//
// written: rms
// Created: 12/98
// Revision: 
//
// Description: This file contains the implementation of 
// GaussQuadRule1d (Quadrature Rule,0).
//
// what: "@(#,0) GaussQuadRule1d.C, revA"

#include <Vector.h>
#include <Matrix.h>

#include <GaussQuadRule1d.h>


GaussQuadRule1d::GaussQuadRule1d (int quadOrder)
{
   this->setOrder(quadOrder);
}

GaussQuadRule1d::~GaussQuadRule1d()
{
   if (coord)
      delete coord;
   if (weight)
      delete weight;
}


int GaussQuadRule1d::setOrder(int quadOrder)
{
   order = quadOrder;

   int numIntPts = order;

   coord  = new Matrix(numIntPts,1);
   weight = new Vector(numIntPts);

   Matrix &xi = *coord;
   Vector &w  = *weight;

   switch (order)
   {
      case 1:
         xi(0,0) =  0.0;

         w(0) = 2.0;
      break;

      case 2:
         xi(0,0) = -0.577350269189626;
         xi(1,0) =  0.577350269189626;

         w(0) = 1.0;
         w(1) = 1.0;
      break;

      case 3:
         xi(0,0) = -0.774596669241483;
         xi(1,0) =  0.0;
         xi(2,0) =  0.774596669241483;

         w(0) = 0.555555555555556;
         w(1) = 0.888888888888889;
         w(2) = 0.555555555555556;
      break;

      case 4:
         xi(0,0) = -0.861136311594053;
         xi(1,0) = -0.339981043584856;
         xi(2,0) =  0.339981043584856;
         xi(3,0) =  0.861136311594053;

         w(0) = 0.347854845137454;
         w(1) = 0.652145154862546;
         w(2) = 0.652145154862546;
         w(3) = 0.347854845137454;

      break;

      case 5:
         xi(0,0) = -0.906179845938664;
         xi(1,0) = -0.538469310105683;
         xi(2,0) =  0.0;
         xi(3,0) =  0.538469310105683;
         xi(4,0) =  0.906179845938664;

         w(0) = 0.236926885056189;
         w(1) = 0.478628670499366;
         w(2) = 0.568888888888889;
         w(3) = 0.478628670499366;
         w(4) = 0.236926885056189;
      break;
         
      case 6:
         xi(0,0) = -0.932469514203152;
         xi(1,0) = -0.661209386466265;
         xi(2,0) = -0.238619186083197;
         xi(3,0) =  0.238619186083197;
         xi(4,0) =  0.661209386466265;
         xi(5,0) =  0.932469514203152;

         w(0) = 0.171324492379170;
         w(1) = 0.360761573048139;
         w(2) = 0.467913934572691;
         w(3) = 0.467913934572691;
         w(4) = 0.360761573048139;
         w(5) = 0.171324492379170;
      break;

      case 7:
         xi(0,0) = -0.949107912342759;
         xi(1,0) = -0.741531185599394;
         xi(2,0) = -0.405845151377397;
         xi(3,0) =  0.0;
         xi(4,0) =  0.405845151377397;
         xi(5,0) =  0.741531185599394;
         xi(6,0) =  0.949107912342759;

         w(0) = 0.129484966168870;
         w(1) = 0.279705391489277;
         w(2) = 0.381830050505119;
         w(3) = 0.417959183673469;
         w(4) = 0.381830050505119;
         w(5) = 0.279705391489277;
         w(6) = 0.129484966168870;
      break;

      case 8:
         xi(0,0) = -0.960289856497536;
         xi(1,0) = -0.796666477413627;
         xi(2,0) = -0.525532409916329;
         xi(3,0) = -0.183434642495650;
         xi(4,0) =  0.183434642495650;
         xi(5,0) =  0.525532409916329;
         xi(6,0) =  0.796666477413627;
         xi(7,0) =  0.960289856497536;

         w(0) = 0.101228536290376;
         w(1) = 0.222381034453374;
         w(2) = 0.313706645877887;
         w(3) = 0.362683783378362;
         w(4) = 0.362683783378362;
         w(5) = 0.313706645877887;
         w(6) = 0.222381034453374;
         w(7) = 0.101228536290376;
      break;

      case 9:
         xi(0,0) = -0.968160239507626;
         xi(1,0) = -0.836031107326636;
         xi(2,0) = -0.613371432700590;
         xi(3,0) = -0.324253423403809;
         xi(4,0) =  0.0;
         xi(5,0) =  0.324253423403809;
         xi(6,0) =  0.613371432700590;
         xi(7,0) =  0.836031107326636;
         xi(8,0) =  0.968160239507626;

         w(0) = 0.081274388361574;
         w(1) = 0.180648160694857;
         w(2) = 0.260610696402935;
         w(3) = 0.312347077040003;
         w(4) = 0.330239355001260;
         w(5) = 0.312347077040003;
         w(6) = 0.260610696402935;
         w(7) = 0.180648160694857;
         w(8) = 0.081274388361574;
      break;

      case 10:
         xi(0,0) = -0.973906528517172;
         xi(1,0) = -0.865063366688985;
         xi(2,0) = -0.679409568299024;
         xi(3,0) = -0.433395394129247;
         xi(4,0) = -0.148874338981631;
         xi(5,0) =  0.148874338981631;
         xi(6,0) =  0.433395394129247;
         xi(7,0) =  0.679409568299024;
         xi(8,0) =  0.865063366688985;
         xi(9,0) =  0.973906528517172;

         w(0) = 0.066671344308688;
         w(1) = 0.149451349150581;
         w(2) = 0.219086362515982;
         w(3) = 0.269266719309996;
         w(4) = 0.295524224714753;
         w(5) = 0.295524224714753;
         w(6) = 0.269266719309996;
         w(7) = 0.219086362515982;
         w(8) = 0.149451349150581;
         w(9) = 0.066671344308688;
      break;

      default:
         cerr << "\n Invalid quadrature order";
         return 0;
      break;
   }    
   return 1;
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
   return *coord;
}

const Vector & 
GaussQuadRule1d::getIntegrPointWeights (void) const
{
   return *weight;
}
