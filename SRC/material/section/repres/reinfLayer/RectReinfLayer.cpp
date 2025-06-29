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
                                                                        
#include <math.h>
#include <Matrix.h>
#include <Vector.h>

#include <ReinfBar.h>
#include <RectReinfLayer.h>
#include <elementAPI.h>

void* OPS_RectReinfLayer()
{
    if(OPS_GetNumRemainingInputArgs() < 8) {
	opserr<<"insufficient arguments for RectReinfLayer\n";
	return 0;
    }

    // matTag ny nz Abar yc zc dy dz
    
    // get idata
    int numData = 3;
    int idata[3];
    if(OPS_GetIntInput(&numData,&idata[0]) < 0) return 0;

    // get data
    double data[5] = {0,0,0,0,0};
    numData = OPS_GetNumRemainingInputArgs();
    if(numData > 5) numData = 5;
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;

    return new RectReinfLayer(idata[0],idata[1],idata[2],
			      data[0],data[1], data[2],
			      data[3],data[4]);
}


RectReinfLayer::RectReinfLayer(void):
  nBarsy(0), nBarsz(0), matID(0), barDiam(0.0),
  area(0.0), yc(0.0), zc(0.0), dy(0.0), dz(0.0)
{
}


RectReinfLayer::RectReinfLayer(int materialID, int ny, int nz, 
                               double reinfBarArea,
			       double yC, double zC,
			       double Dy, double Dz):
  nBarsy(ny), nBarsz(nz), matID(materialID), barDiam(0.0),
  area(reinfBarArea), yc(yC), zc(zC), dy(Dy), dz(Dz)
{
  if (nBarsy < 0)
    nBarsy = 0;

  if (nBarsz < 0)
    nBarsz = 0;  
}

RectReinfLayer::~RectReinfLayer()
{

}


void RectReinfLayer::setNumReinfBars(int numReinfBars)
{
  //nReinfBars = numReinfBars;
}

void
RectReinfLayer::setMaterialID (int materialID)
{
   matID = materialID;
}

void
RectReinfLayer::setReinfBarDiameter (double reinfBarDiameter)
{
   barDiam = reinfBarDiameter;
   double pi = acos(-1.0);
   area = pi * barDiam*barDiam/4.0;
}

void
RectReinfLayer::setReinfBarArea(double reinfBarArea)
{
   area = reinfBarArea;
}


int
RectReinfLayer::getNumReinfBars (void) const
{
  return 4 + 2*(nBarsy + nBarsz);
}

int
RectReinfLayer::getMaterialID (void) const
{
  return matID;
}

double
RectReinfLayer::getReinfBarDiameter (void) const
{
  return barDiam;
}

double
RectReinfLayer::getReinfBarArea (void) const
{
  return area;
}

ReinfBar * 
RectReinfLayer::getReinfBars (void) const
{
  int nbars = 4 + 2*(nBarsy+nBarsz);
  ReinfBar *reinfBars = new ReinfBar[nbars];

  static Vector barPosit(2);
  
  // Corner bars
  barPosit(0) = yc - 0.5*dy; barPosit(1) = zc - 0.5*dz;
  reinfBars[0].setPosition(barPosit);
  barPosit(0) = yc - 0.5*dy; barPosit(1) = zc + 0.5*dz;
  reinfBars[1].setPosition(barPosit);
  barPosit(0) = yc + 0.5*dy; barPosit(1) = zc - 0.5*dz;
  reinfBars[2].setPosition(barPosit);
  barPosit(0) = yc + 0.5*dy; barPosit(1) = zc + 0.5*dz;
  reinfBars[3].setPosition(barPosit);

  int loc = 4;
  if (nBarsy > 0) {
    double yincr = dy/(nBarsy+1);
    for (int i = 1; i <= nBarsy; i++) {
      barPosit(0) = yc - 0.5*dy + i*yincr;
      
      barPosit(1) = zc - 0.5*dz;
      reinfBars[loc++].setPosition(barPosit);
      barPosit(1) = zc + 0.5*dz;
      reinfBars[loc++].setPosition(barPosit);
    }
  }

  if (nBarsz > 0) {
    double zincr = dz/(nBarsz+1);
    for (int i = 1; i <= nBarsz; i++) {
      barPosit(1) = zc - 0.5*dz + i*zincr;
      
      barPosit(0) = yc - 0.5*dy;
      reinfBars[loc++].setPosition(barPosit);
      barPosit(0) = yc + 0.5*dy;
      reinfBars[loc++].setPosition(barPosit);
    }    
  }
  
  for (int i = 0; i < nbars; i++)
    reinfBars[i].setArea(area);
  
  return reinfBars;         
}


ReinfLayer * 
RectReinfLayer::getCopy (void) const
{
  RectReinfLayer *theCopy = new RectReinfLayer (matID, nBarsy, nBarsz, area,
						yc, zc, dy, dz);
  return theCopy;
}



void RectReinfLayer::Print(OPS_Stream &s, int flag) const
{
  if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "\nReinforcing Layer type:  Rect";
    s << "\nMaterial ID: " << matID;
    s << "\nReinf. bar diameter: " << barDiam;
    s << "\nReinf. bar area: " << area;
    s << "\nCenter Position: " << yc << ' ' << zc;
  }
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t\t{\"type\": \"layerRect\", \"material\": "<<matID<<", "<<",\"area\": "<<area<<", \"nReinfBars\": " << 4+2*(nBarsy+nBarsz) <<", ";
    s << "\"center\": ["<< yc <<","<< zc <<"], ";
    s <<"}";
  }
}


OPS_Stream &operator<<(OPS_Stream &s, const RectReinfLayer &RectReinfLayer)
{  
   RectReinfLayer.Print(s);
   return s;
}
 
