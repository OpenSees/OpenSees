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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:01:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/patch/CircPatch.cpp,v $
                                                                        
                                                                        
// File: CircPatch.C
// Written by Remo M. de Souza
// December 1998

#include <math.h>
#include <Matrix.h>
#include <Patch.h>
#include <CircPatch.h>
#include <QuadCell.h>
#include <CircSectionCell.h>
#include <elementAPI.h>

void* OPS_CircPatch()
{
    if(OPS_GetNumRemainingInputArgs() < 9) {
	opserr<<"insufficient arguments for CircPatch\n";
	return 0;
    }

    // get idata
    int numData = 3;
    int idata[3];
    if(OPS_GetIntInput(&numData,&idata[0]) < 0) return 0;

    // get data
    double data[6];
    numData = 6;
    static Vector centerPos(2);
    /*centerPos(0) = data[0];
    centerPos(1) = data[1];*/
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;
    centerPos(0) = data[0];
    centerPos(1) = data[1];

    return new CircPatch(idata[0],idata[1],idata[2],centerPos,
			 data[2],data[3],data[4],data[5]);
}


CircPatch::CircPatch(): matID(0), nDivCirc(1), nDivRad(1), centerPosit(2),
                        intRad(0.0), extRad(0.0), initAng(0.0), finalAng(0.0)
{

}


CircPatch::CircPatch(int materialID, int numSubdivCircunf, int numSubdivRadial,
                     const Vector &centerPosition, double internRadius, 
                     double externRadius, double initialAngle, double finalAngle):
                     matID(materialID), 
                     nDivCirc(numSubdivCircunf), nDivRad(numSubdivRadial),
                     centerPosit(centerPosition),
                     intRad(internRadius), extRad(externRadius), 
                     initAng(initialAngle), finalAng(finalAngle)
{

}


CircPatch::~CircPatch()
{

}

void CircPatch::setMaterialID(int materialID)
{
   matID = materialID;
}


void CircPatch::setDiscretization(int numSubdivCircunf, int numSubdivRadial)
{
   nDivRad  = numSubdivRadial;
   nDivCirc = numSubdivCircunf;
}

void CircPatch::setCenterPosition(const Vector &centerPosition)
{
   centerPosit = centerPosition;
}

void CircPatch::setRadii(double internRadius, double externRadius)
{
   intRad = internRadius;
   extRad = externRadius;
}

void CircPatch::setAngles(double initialAngle, double finalAngle)
{
   initAng  = initialAngle;
   finalAng = finalAngle;
}


int CircPatch::getMaterialID(void) const
{
   return matID;
}
 
void CircPatch::getDiscretization(int &numSubdivCircunf, int &numSubdivRadial) const
{
   numSubdivCircunf = nDivCirc;
   numSubdivRadial  = nDivRad;
}

void CircPatch::getRadii(double &internRadius, double &externRadius) const
{
   internRadius = intRad;
   externRadius = extRad;
}

void CircPatch::getAngles(double &initialAngle, double &finalAngle) const
{
   initialAngle = initAng;
   finalAngle   = finalAng;
}

const Vector & CircPatch::getCenterPosition (void) const
{
   return centerPosit;
}

int CircPatch::getNumCells (void) const
{
   return nDivCirc * nDivRad;
}

Cell **
CircPatch::getCells (void) const
{
   double pi = acos(-1.0);
   double deltaRad, deltaTheta; 
   double initAngRadians, finalAngRadians;
   double rad_j, rad_j1, theta_i, theta_i1;
   Matrix cellVertCoord(4,2);

   int    i, j, k;
   int    numCells;
   Cell   **cells;

   if (nDivRad > 0  && nDivCirc > 0)
   {
      numCells  = this->getNumCells();

      cells = new Cell* [numCells];
      
      if (!cells)
         return 0;

      initAngRadians  = pi * initAng  / 180.0;
      finalAngRadians = pi * finalAng / 180.0;

      deltaRad   = (extRad - intRad) / nDivRad;
      deltaTheta = (finalAngRadians - initAngRadians) / nDivCirc;

      k = 0;
      for (j = 0; j < nDivRad; j++)
      {
         rad_j  = intRad + deltaRad*j;
         rad_j1 = rad_j + deltaRad;

         for (i = 0; i < nDivCirc; i++)
         {
            // compute coordinates
                       
            theta_i  = initAngRadians + deltaTheta*i;

	    /*
            theta_i1 = theta_i + deltaTheta;
            cellVertCoord(0,0) = centerPosit(0) + rad_j  * cos(theta_i1);
            cellVertCoord(0,1) = centerPosit(1) + rad_j  * sin(theta_i1);
            cellVertCoord(1,0) = centerPosit(0) + rad_j  * cos(theta_i);
            cellVertCoord(1,1) = centerPosit(1) + rad_j  * sin(theta_i);
            cellVertCoord(2,0) = centerPosit(0) + rad_j1 * cos(theta_i);
            cellVertCoord(2,1) = centerPosit(1) + rad_j1 * sin(theta_i);
            cellVertCoord(3,0) = centerPosit(0) + rad_j1 * cos(theta_i1);
            cellVertCoord(3,1) = centerPosit(1) + rad_j1 * sin(theta_i1);

            cells[k] = new QuadCell(cellVertCoord); 
	    */

            theta_i1 = theta_i + deltaTheta/2.0;
	    cells[k] = new CircSectionCell(rad_j, rad_j1, deltaTheta, theta_i1, centerPosit(0), centerPosit(1));

            k++; 
         }
      }
   }
   else
      return 0;

   return cells;
}


Patch * 
CircPatch::getCopy (void) const
{
   CircPatch *theCopy = new CircPatch (matID, nDivCirc, nDivRad,
                                       centerPosit,  intRad, extRad,
                                       initAng, finalAng);
   return theCopy;
}
 
void CircPatch::Print(OPS_Stream &s, int flag) const
{
   s << "\nPatch Type: CircPatch";
   s << "\nMaterial Id: " << matID;
   s << "\nNumber of subdivisions in the radial direction: " << nDivRad;
   s << "\nNumber of subdivisions in the circunferential direction: " << nDivCirc;
   s << "\nCenter Position: " << centerPosit;
   s << "\nInternal Radius: " << intRad << "\tExternal Radius: " << extRad;
   s << "\nInitial Angle: " << initAng << "\tFinal Angle: " << finalAng;
}


OPS_Stream &operator<<(OPS_Stream &s, CircPatch &CircPatch)
{
   CircPatch.Print(s);
   return s;
}
