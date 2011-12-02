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
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/quadrule/GaussLobattoQuadRule1d01.cpp,v $
                                                                        
                                                                        
// File: ~/QuadRule/GaussLobattoQuadRule1d0101.C
//
// written: rms
// Created: 12/98
// Revision: 
//
// Description: This file contains the implementation of 
// GaussLobattoQuadRule1d01 (Quadrature Rule,0).
//
// what: "@(#,0) GaussLobattoQuadRule1d01.C, revA"

#include <Vector.h>
#include <Matrix.h>

#include <GaussLobattoQuadRule1d01.h>


GaussLobattoQuadRule1d01::GaussLobattoQuadRule1d01 (int quadOrder)
{
   this->setOrder(quadOrder);
}

GaussLobattoQuadRule1d01::~GaussLobattoQuadRule1d01()
{
   if (coord)
      delete coord;
   if (weight)
      delete weight;
}


int GaussLobattoQuadRule1d01::setOrder(int quadOrder)
{
   order = quadOrder;

   int numIntPts = order;

   coord  = new Matrix(numIntPts,1);
   weight = new Vector(numIntPts);

   Matrix &xi = *coord;
   Vector &w  = *weight;

   switch (order)
   {

      case 2:
         xi(0,0) = -1.;
         xi(1,0) =  1.;

         w(0) =  1.;
         w(1) =  1.;
      break; 
     

      case 3:
         xi(0,0) = -1.;
         xi(1,0) =  0.;
         xi(2,0) =  1.;

         w(0) =  1/3.;
         w(1) =  4/3.;
         w(2) =  1/3.;
      break;


      case 4:
         xi(0,0) = -1.;
         xi(1,0) = -0.44721360;
         xi(2,0) =  0.44721360;
         xi(3,0) =  1.;

         w(0) =  1/6.;
         w(1) =  5/6.;
         w(2) =  5/6.;
         w(3) =  1/6.;
      break;


      case 5:
         xi(0,0) = -1.;
         xi(1,0) = -0.65465367;
         xi(2,0) =  0.;
         xi(3,0) =  0.65465367;
         xi(4,0) =  1.;
      
         w(0) =  0.1;
         w(1) =  0.5444444444;
         w(2) =  0.7111111111;
         w(3) =  0.5444444444;
         w(4) =  0.1;
      break;

     
      case 6:
         xi(0,0) = -1.;
         xi(1,0) = -0.7650553239;
         xi(2,0) = -0.2852315164;
         xi(3,0) =  0.2852315164;
         xi(4,0) =  0.7650553239;
         xi(5,0) =  1.;
               
         w(0) =  0.06666666667;
         w(1) =  0.3784749562;
         w(2) =  0.5548583770;
         w(3) =  0.5548583770;
         w(4) =  0.3784749562;
         w(5) =  0.06666666667;
      break;


      case 7:
         xi(0,0) = -1.; 
         xi(1,0) = -0.8302238962;
         xi(2,0) = -0.4688487934;
         xi(3,0) =  0.;
         xi(4,0) =  0.4688487934;
         xi(5,0) =  0.8302238962;
         xi(6,0) =  1.; 
         

         w(0) =  0.04761904762;
         w(1) =  0.2768260473;
         w(2) =  0.4317453812;
         w(3) =  0.4876190476;
         w(4) =  0.4317453812;
         w(5) =  0.2768260473;
         w(6) =  0.04761904762;
      break;

      case 8:
         xi(0,0) = -1.;
         xi(1,0) = -0.8717401485;
         xi(2,0) = -0.5917001814;
         xi(3,0) = -0.2092992179;
         xi(4,0) =  0.2092992179;
         xi(5,0) =  0.5917001814;
         xi(6,0) =  0.8717401485;
         xi(7,0) =  1.;

         w(0) =  0.03571428571; 
         w(1) =  0.2107042271;
         w(2) =  0.3411226924;
         w(3) =  0.4124587946;
         w(4) =  0.4124587946;
         w(5) =  0.3411226924;
         w(6) =  0.2107042271;
         w(7) =  0.03571428571; 
      break;


      case 9:
         xi(0,0) = -1.;
         xi(1,0) = -0.8997579954;
         xi(2,0) = -0.6771862795;
         xi(3,0) = -0.3631174638;
         xi(4,0) =  0.;
         xi(5,0) =  0.3631174638;
         xi(6,0) =  0.6771862795;
         xi(7,0) =  0.8997579954;
         xi(8,0) =  1.;


         w(0) =  0.02777777778;
         w(1) =  0.1654953615;
         w(2) =  0.2745387125;
         w(3) =  0.3464285109;
         w(4) =  0.3715192743;
         w(5) =  0.3464285109;
         w(6) =  0.2745387125;
         w(7) =  0.1654953615;
         w(8) =  0.02777777778;
      break;

                
      case 10:
         xi(0,0) = -1.;
         xi(1,0) = -0.9195339082;
         xi(2,0) = -0.7387738651;
         xi(3,0) = -0.4779249498;  
         xi(4,0) = -0.1652789577;
         xi(5,0) =  0.1652789577;
         xi(6,0) =  0.4779249498;  
         xi(7,0) =  0.7387738651;
         xi(8,0) =  0.9195339082;
         xi(9,0) =  1.;


         w(0) =  0.02222222222;
         w(1) =  0.1333059908;
         w(2) =  0.2248893421;
         w(3) =  0.2920426836;
         w(4) =  0.3275397611;
         w(5) =  0.3275397611;
         w(6) =  0.2920426836;
         w(7) =  0.2248893421;
         w(8) =  0.1333059908;
         w(9) =  0.02222222222;
      break;


      default:
         cerr << "\n Invalid quadrature order";
         return 0;
//      break;
   }    

   xi = (xi + 1)/2;
   w /= 2;

   return 1;
}

int GaussLobattoQuadRule1d01::getOrder (void) const
{
   return order;
}

int GaussLobattoQuadRule1d01::getNumIntegrPoints (void) const
{
   return order;
}

const Matrix & 
GaussLobattoQuadRule1d01::getIntegrPointCoords (void) const
{
   return *coord;
}

const Vector & 
GaussLobattoQuadRule1d01::getIntegrPointWeights (void) const
{
   return *weight;
}
