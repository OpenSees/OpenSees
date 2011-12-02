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
// $Date: 2003-02-14 23:01:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/reinfLayer/StraightReinfLayer.cpp,v $
                                                                        
                                                                        
// File: StraightReinfLayer.C 
// Written by Remo M. de Souza 
// December 1998

#include <math.h>
#include <Matrix.h>
#include <Vector.h>

#include <ReinfBar.h>
#include <StraightReinfLayer.h>


StraightReinfLayer::StraightReinfLayer(void):
                        nReinfBars(0), matID(0), barDiam(0.0),
                        area(0.0), initPosit(2), finalPosit(2)          
{

}

StraightReinfLayer::StraightReinfLayer(int materialID, int numReinfBars, 
                                       double reinfBarArea,
                                       const Vector &InitialPosition,
                                       const Vector &FinalPosition):
                                         nReinfBars(numReinfBars),
                                         matID(materialID),
                                         area(reinfBarArea),
                                         barDiam(0.0),
                                         initPosit(InitialPosition),
                                         finalPosit(FinalPosition)
{
}


StraightReinfLayer::~StraightReinfLayer()
{

}


void StraightReinfLayer::setNumReinfBars(int numReinfBars)
{
   nReinfBars = numReinfBars;
}

void StraightReinfLayer::setMaterialID (int materialID)
{
   matID = materialID;
}

void StraightReinfLayer::setReinfBarDiameter (double reinfBarDiameter)
{
   barDiam = reinfBarDiameter;
   double pi = acos(-1.0);
   area = pi * barDiam*barDiam/4.0;
}

void StraightReinfLayer::setReinfBarArea(double reinfBarArea)
{
   area = reinfBarArea;
}

void StraightReinfLayer::setInitialPosition (const Vector &initialPosition)
{
   initPosit = initialPosition;
}

void StraightReinfLayer::setFinalPosition (const Vector &finalPosition)
{
   finalPosit = finalPosition;
}


int StraightReinfLayer::getNumReinfBars (void) const
{
   return nReinfBars;
}

int StraightReinfLayer::getMaterialID (void) const
{
   return matID;
}

double StraightReinfLayer::getReinfBarDiameter (void) const
{
   return barDiam;
}

double StraightReinfLayer::getReinfBarArea (void) const
{
   return area;
}

ReinfBar * 
StraightReinfLayer::getReinfBars (void) const
{
   double dy, dz;
   Vector barPosit(2);
   int i;
   ReinfBar *reinfBars;

   if (nReinfBars == 1)
   {
      barPosit(0) = (initPosit(0) + finalPosit(0)) / 2;
      barPosit(1) = (initPosit(1) + finalPosit(1)) / 2;
    
      reinfBars = new ReinfBar [1];

      reinfBars[0].setPosition(barPosit);
      reinfBars[0].setArea(this->area);
   }

   else if (nReinfBars > 1)
   {
      dy = (finalPosit(0) - initPosit(0))/(nReinfBars - 1);
      dz = (finalPosit(1) - initPosit(1))/(nReinfBars - 1);

      reinfBars = new ReinfBar [nReinfBars];

      for (i = 0; i < nReinfBars; i++)
      {
         barPosit(0) = initPosit(0) + dy * i;
         barPosit(1) = initPosit(1) + dz * i;

         reinfBars[i].setPosition(barPosit);
         reinfBars[i].setArea(this->area);
      }
   }
   else
     return 0;

   return reinfBars;         
}

const Vector & 
StraightReinfLayer::getInitialPosition (void) const
{
   return initPosit;
}

const Vector & 
StraightReinfLayer::getFinalPosition   (void) const
{
   return finalPosit;
}


ReinfLayer * 
StraightReinfLayer::getCopy (void) const
{
   StraightReinfLayer *theCopy = new StraightReinfLayer (matID,
                                                 nReinfBars, area,
                                                 initPosit, finalPosit);
   return theCopy;
}



void StraightReinfLayer::Print(OPS_Stream &s, int flag) const
{
   s << "\nReinforcing Layer type:  Straight";
   s << "\nMaterial ID: " << matID;
   s << "\nReinf. bar diameter: " << barDiam;
   s << "\nReinf. bar area: " << area;
   s << "\nInitial Position: " << initPosit;
   s << "\nFinal Position: " << finalPosit;
}


OPS_Stream &operator<<(OPS_Stream &s, const StraightReinfLayer &straightReinfLayer)
{  
   straightReinfLayer.Print(s);
   return s;
}
 
