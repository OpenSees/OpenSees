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
// $Date: 2003/05/02 18:34:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/cell/QuadCell.cpp,v $
                                                                        
                                                                        
// File: QuadCell.C
// Written by Remo M. de Souza
// December 1998

#include <Matrix.h>
#include <Vector.h>

#include <QuadCell.h>


QuadCell::QuadCell(void):
                   vertCoord(4,2), Centroid(2)                 
{

}


QuadCell::QuadCell(const Matrix &vertexCoords): 
                   vertCoord(vertexCoords), Centroid(2)
{

}


QuadCell::~QuadCell()
{

}

const Matrix &
QuadCell::getVertCoords (void) const
{
   return vertCoord;
}

double QuadCell::getdValue (void) const
{
    double dVa = vertCoord(0,0);
    return dVa;
}

void QuadCell::setVertCoords (const Matrix &vertexCoords)
{
   vertCoord = vertexCoords;
}


double QuadCell::getArea (void) const
{
   double area;
   double x0, y0, x1, y1, x2, y2, x3, y3;

//   opserr << "cell vertCoord: " << vertCoord;
 
   x0 = vertCoord(0,0);
   y0 = vertCoord(0,1);
   x1 = vertCoord(1,0);
   y1 = vertCoord(1,1);
   x2 = vertCoord(2,0);
   y2 = vertCoord(2,1);
   x3 = vertCoord(3,0);
   y3 = vertCoord(3,1);

   area = ((x2-x1)*(y0-y1) - (x0-x1)*(y2-y1) +
           (x0-x3)*(y2-y3) - (x2-x3)*(y0-y3)) / 2.0; 

//   opserr << "area1=" << area;

   int i, i1;
   double yi, zi, yi1, zi1;
   area = 0;

   for (i = 0; i < 4; i++)
   {
      i1 = (i+1)%4;
      yi  = vertCoord(i,0);
      zi  = vertCoord(i,1);
      yi1 = vertCoord(i1,0);
      zi1 = vertCoord(i1,1);

      area += (zi1 - zi) * (yi1 + yi); 
   }
   area /= 2.0;

//   opserr << "area2= " << area << endln;

   return area;
                
}


const Vector & 
QuadCell::getCentroidPosition(void)
{
   int i, i1;
   double yi, zi, yi1, zi1, dyi, dzi;
   double area, integ;
   double CGy = 0.0, CGz = 0.0;

   area = this->getArea();

   for (i = 0; i < 4; i++)
   {
      i1 = (i+1)%4;
        
      yi  = vertCoord(i,0);
      zi  = vertCoord(i,1);
      yi1 = vertCoord(i1,0);
      zi1 = vertCoord(i1,1);

      dyi = yi1 - yi;
      dzi = zi1 - zi;
   
      integ = yi*zi + (yi*dzi + zi*dyi)/2.0 + dyi*dzi/3.0;

      CGy -= dyi * integ;
      CGz += dzi * integ;
   }
   
   CGy /= area;
   CGz /= area;

   Centroid(0) = CGy;
   Centroid(1) = CGz;

//   opserr << "\narea : " << area << " centroid: " << Centroid;
 
   return Centroid;
}



void QuadCell::Print(OPS_Stream &s, int flag) const
{
   s << "\nCell Type: QuadCell";
   s << "\nVertex Coordinates: " << vertCoord;
}


OPS_Stream &operator<<(OPS_Stream &s, const QuadCell &quadCell)
{
   quadCell.Print(s);
   return s;
}    

