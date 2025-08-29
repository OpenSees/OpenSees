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
#include <cmath>
#include <OPS_Stream.h>

#include <CircSectionCell.h>

CircSectionCell::CircSectionCell(double r1, double r2, double alpha, double theta, double offsetX,
                                 double offsetY)
{

  double a = alpha / 2.0;
  double At = a * r2 * r2;
  double ct = 2.0 * r2 * sin(a) / (3.0 * a);
  double A1 = a * r1 * r1;
  double c1 = 2.0 * r1 * sin(a) / (3.0 * a);

  area     = At - A1;
  double c = (At * ct - A1 * c1) / area;

  location[0] = std::cos(theta) * c + offsetX;
  location[1] = std::sin(theta) * c + offsetY;
}