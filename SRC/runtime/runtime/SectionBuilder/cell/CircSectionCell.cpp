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
                                                                        
// Written: fmk
// March 2014

#include <Matrix.h>
#include <Vector.h>
#include <math.h>

#include <CircSectionCell.h>

CircSectionCell::CircSectionCell(double R1, double R2, double Alpha, double Theta, double offX, double offY)
  :r1(R1), r2(R2), alpha(Alpha), theta(Theta), A(0.0), Centroid(2), offsetX(offX), offsetY(offY)
{

  double a = alpha/2.0;

  double At = a*r2*r2;
  double ct = 2.0*r2*sin(a)/(3.0*a);
  double A1 = a*r1*r1;
  double c1 = 2.0*r1*sin(a)/(3.0*a);
  
  A = At - A1;
  double c = (At*ct - A1*c1)/A;

  Centroid(0) = cos(Theta)*c + offX;
  Centroid(1) = sin(Theta)*c + offY;
}


CircSectionCell::~CircSectionCell()
{

}

double CircSectionCell::getArea (void) const
{
  return A;
}

double CircSectionCell::getdValue (void) const
{
  return 1.0; // TODO: Should be something meaningful -- MHS
}

const Vector & 
CircSectionCell::getCentroidPosition(void)
{
   return Centroid;
}



void CircSectionCell::Print(OPS_Stream &s, int flag) const
{
   s << "\nCell Type: CircSectionCell";
   s << "\n\tr1: " << r1 << " r2: " << r2 << " alpha: " << alpha << " theta: " << theta << endln;
}


OPS_Stream &operator<<(OPS_Stream &s, const CircSectionCell &cell)
{
   cell.Print(s);
   return s;
}    

