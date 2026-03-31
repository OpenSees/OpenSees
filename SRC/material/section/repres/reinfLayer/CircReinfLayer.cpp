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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/reinfLayer/CircReinfLayer.cpp,v $
                                                                        
                                                                        
// File: CircReinfLayer.C 
// Written by Remo M. de Souza 
// December 1998

#include <math.h>
#include <Matrix.h>
#include <Vector.h>

#include <ReinfBar.h>
#include <CircReinfLayer.h>
#include <elementAPI.h>

void* OPS_CircReinfLayer()
{
    if(OPS_GetNumRemainingInputArgs() < 6) {
	opserr<<"insufficient arguments for CircReinfLayer\n";
	return 0;
    }

    // get idata
    int numData = 2;
    int idata[2];
    if(OPS_GetIntInput(&numData,&idata[0]) < 0) return 0;

    // get data
    double data[6] = {0,0,0,0,0,0};
    numData = OPS_GetNumRemainingInputArgs();
    if(numData > 6) numData = 6;
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;
    static Vector cpos(2);
    cpos(0) = data[1];
    cpos(1) = data[2];

    if(numData < 6) {
	return new CircReinfLayer(idata[0],idata[1],data[0],
				  cpos,data[3]);
    } else {
	return new CircReinfLayer(idata[0],idata[1],data[0],
				  cpos,data[3],data[4],data[5]);
    }
}


CircReinfLayer::CircReinfLayer(void):
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
   double angle = finalAng-initAng;
   if (angle >= 360.0 && nReinfBars > 0)
      finalAng = initAng + 360.0 - 360.0/nReinfBars;
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


int CircReinfLayer::getNumReinfBars (void) const
{
   return nReinfBars;
}

int CircReinfLayer::getMaterialID (void) const
{
   return matID;
}

double CircReinfLayer::getReinfBarDiameter (void) const
{
   return barDiam;
}

double CircReinfLayer::getReinfBarArea (void) const
{
   return area;
}

ReinfBar * 
CircReinfLayer::getReinfBars (void) const
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
CircReinfLayer::getCenterPosition(void) const
{
   return centerPosit;
}

double CircReinfLayer::getArcRadius(void) const 
{
   return arcRad;
}

double CircReinfLayer::getInitAngle(void) const 
{
   return initAng;
}

double CircReinfLayer::getFinalAngle(void) const 
{
   return finalAng;
}


ReinfLayer * 
CircReinfLayer::getCopy (void) const
{
   CircReinfLayer *theCopy = new CircReinfLayer (matID, nReinfBars, area,
                                                 centerPosit, arcRad,
                                                 initAng, finalAng);
   return theCopy;
}



void CircReinfLayer::Print(OPS_Stream &s, int flag) const
{
	if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
   s << "\nReinforcing Layer type:  Circ";
   s << "\nMaterial ID: " << matID;
   s << "\nReinf. bar diameter: " << barDiam;
   s << "\nReinf. bar area: " << area;
   s << "\nCenter Position: " << centerPosit;
   s << "\nArc Radius: " << arcRad;
   s << "\nInitial angle: " << initAng;
   s << "\nFinal angle: " << finalAng;
}
if (flag == OPS_PRINT_PRINTMODEL_JSON) {
	 s << "\t\t\t\t{\"type\": \"layerCirc\", \"material\": "<<matID<<", "<<",\"area\": "<<area<<", \"nReinfBars\": " << nReinfBars<<", ";
	 s << "\"center\": ["<<centerPosit(0)<<","<<centerPosit(1)<<"], ";
	 s << "\"radius\": "<<arcRad<<", ";
	 s << "\"angles\": ["<<initAng<<", "<<finalAng<<"]";
	 s <<"}";
}
}


OPS_Stream &operator<<(OPS_Stream &s, const CircReinfLayer &CircReinfLayer)
{  
   CircReinfLayer.Print(s);
   return s;
}
 
