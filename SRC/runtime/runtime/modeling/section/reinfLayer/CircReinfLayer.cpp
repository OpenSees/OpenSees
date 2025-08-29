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
//
// File: CircReinfLayer.C
// Written by Remo M. de Souza
// December 1998

#include <Matrix.h>
#include <OPS_Stream.h>
#include <Vector.h>
#include <math.h>
#include <string>

#include <CircReinfLayer.h>
#include <Cell.h>



CircReinfLayer::CircReinfLayer(int material, int numReinfBars, double area,
                               const VectorND<2>& centerPosition, double arcRadius, double initialAngle,
                               double finalAngle)
 : ReinfLayer(material, area),
   nReinfBars(numReinfBars),
   centerPosit(centerPosition),
   arcRad(arcRadius),
   initAng(initialAngle),
   finalAng(finalAngle)
{
}

CircReinfLayer::CircReinfLayer(int material, int numReinfBars, double area,
                               const VectorND<2>& centerPosition, double radius)
 : ReinfLayer(material, area),
   nReinfBars(numReinfBars),
   centerPosit(centerPosition),
   arcRad(radius),
   initAng(0.0),
   finalAng(0.0)
{
  // Figure out final angle so that complete circle does not put
  // two bars at the same location
  if (nReinfBars > 0)
    finalAng = 360.0 - 360.0 / nReinfBars;
}

int
CircReinfLayer::getNumReinfBars() const
{
  return nReinfBars;
}

std::vector<Cell>
CircReinfLayer::getReinfBars() const
{
  std::vector<Cell> bars(nReinfBars);

  double pi = acos(-1.0);

  if (nReinfBars > 0) {
    double initAngRad  = pi * initAng / 180.0;
    double finalAngRad = pi * finalAng / 180.0;

    double dtheta;
    if (nReinfBars > 1)
      dtheta = (finalAngRad - initAngRad) / (nReinfBars - 1);
    else
      dtheta = 0.0; // Doesn't really matter what this is


    for (int i = 0; i < nReinfBars; i++) {
      double theta = initAngRad + dtheta * i;
      VectorND<2> position {
           centerPosit(0) + arcRad * cos(theta),
           centerPosit(1) + arcRad * sin(theta)
      };
      bars[i] = Cell(this->getMaterialID(), this->getCellArea(), position);
    }
  }

  return bars;
}
