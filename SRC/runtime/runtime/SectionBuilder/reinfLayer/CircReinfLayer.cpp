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

#include <math.h>
#include <string>
#include <OPS_Stream.h>
#include <Matrix.h>
#include <Vector.h>

#include <ReinfBar.h>
#include <CircReinfLayer.h>

CircReinfLayer::CircReinfLayer():
                        nReinfBars(0), matID(0), barDiam(0.0),
                        area(0.0), centerPosit(2), arcRad(0.0),
                        initAng(0.0), finalAng(0.0)          
{
}


CircReinfLayer::CircReinfLayer(int materialID, int numReinfBars, 
                               double reinfBarArea,
                               const Vector &centerPosition,
                               double arcRadius, double initialAngle,
                               double finalAngle):
                               nReinfBars(numReinfBars),
                               matID(materialID), area(reinfBarArea),
                               barDiam(0.0),centerPosit(centerPosition),
                               arcRad(arcRadius),initAng(initialAngle), 
                               finalAng(finalAngle)
{
}

CircReinfLayer::CircReinfLayer(int materialID, int numReinfBars, double  reinfBarArea,
                                                           const Vector &centerPosition, double radius):
nReinfBars(numReinfBars), matID(materialID), area(reinfBarArea),
barDiam(0.0), centerPosit(centerPosition), arcRad(radius),
initAng(0.0), finalAng(0.0)
{
        // Figure out final angle so that complete circle does not put
        // two bars at the same location
        if (nReinfBars > 0)
                finalAng = 360.0 - 360.0/nReinfBars;
}

CircReinfLayer::~CircReinfLayer()
{

}


void CircReinfLayer::setNumReinfBars(int numReinfBars)
{
   nReinfBars = numReinfBars;
}

void CircReinfLayer::setMaterialID (int materialID)
{
   matID = materialID;
}

void CircReinfLayer::setReinfBarDiameter (double reinfBarDiameter)
{
   barDiam = reinfBarDiameter;
   double pi = acos(-1.0);
   area = pi * barDiam*barDiam/4.0;
}

void CircReinfLayer::setReinfBarArea(double reinfBarArea)
{
   area = reinfBarArea;
}


int CircReinfLayer::getNumReinfBars() const
{
   return nReinfBars;
}

int CircReinfLayer::getMaterialID() const
{
   return matID;
}

double CircReinfLayer::getReinfBarDiameter() const
{
   return barDiam;
}

double CircReinfLayer::getReinfBarArea() const
{
   return area;
}

ReinfBar * 
CircReinfLayer::getReinfBars() const
{
   double theta, dtheta;
   static Vector barPosit(2);
   int i;
   ReinfBar *reinfBars;
   double pi = acos(-1.0);
   double initAngRad, finalAngRad;

   if (nReinfBars > 0)
   {
      initAngRad  = pi * initAng  / 180.0;
      finalAngRad = pi * finalAng / 180.0;

      if (nReinfBars > 1) 
        dtheta = (finalAngRad - initAngRad) /(nReinfBars - 1);
      else
        dtheta = 0.0; // Doesn't really matter what this is

      reinfBars = new ReinfBar [nReinfBars];

      for (i = 0; i < nReinfBars; i++)
      {
         theta = initAngRad + dtheta * i;
         barPosit(0) = centerPosit(0) + arcRad*cos(theta);
         barPosit(1) = centerPosit(1) + arcRad*sin(theta);

         reinfBars[i].setPosition(barPosit);
         reinfBars[i].setArea(this->area);
      }
   }
   else
      return 0;

   return reinfBars;         
}


const Vector & 
CircReinfLayer::getCenterPosition() const
{
   return centerPosit;
}

double CircReinfLayer::getArcRadius() const 
{
   return arcRad;
}

double CircReinfLayer::getInitAngle() const 
{
   return initAng;
}

double CircReinfLayer::getFinalAngle() const 
{
   return finalAng;
}


ReinfLayer * 
CircReinfLayer::getCopy() const
{
   CircReinfLayer *theCopy = new CircReinfLayer (matID, nReinfBars, area,
                                                 centerPosit, arcRad,
                                                 initAng, finalAng);
   return theCopy;
}



void CircReinfLayer::Print(OPS_Stream &s, int flag) const
{
   s << "\nReinforcing Layer type:  Circ";
   s << "\nMaterial ID: " << matID;
   s << "\nReinf. bar diameter: " << barDiam;
   s << "\nReinf. bar area: " << area;
   s << "\nCenter Position: " << centerPosit;
   s << "\nArc Radius: " << arcRad;
   s << "\nInitial angle: " << initAng;
   s << "\nFinal angle: " << finalAng;
}


OPS_Stream &operator<<(OPS_Stream &s, const CircReinfLayer &CircReinfLayer)
{  
   CircReinfLayer.Print(s);
   return s;
}
 
