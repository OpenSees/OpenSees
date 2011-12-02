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
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/patch/CircPatch.cpp,v $
                                                                        
                                                                        
// File: CircPatch.C
// Written by Remo M. de Souza
// December 1998

#include <math.h>
#include <Matrix.h>
#include <Patch.h>
#include <CircPatch.h>
#include <QuadCell.h>


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

      //cerr << "\nnumCells: " << numCells;

      cells = new Cell* [numCells];
      
      if (!cells)
         return 0;

      initAngRadians  = pi * initAng  / 180.0;
      finalAngRadians = pi * finalAng / 180.0;

      deltaRad   = (extRad - intRad) / nDivRad;
      deltaTheta = (finalAngRadians - initAngRadians) / nDivCirc;

      //cerr << "\ndeltaRad: " << deltaRad;

      k = 0;
      for (j = 0; j < nDivRad; j++)
      {
         rad_j  = intRad + deltaRad*j;
         rad_j1 = rad_j + deltaRad;

         for (i = 0; i < nDivCirc; i++)
         {
            // compute coordinates
                       
            theta_i  = initAngRadians + deltaTheta*i;
            theta_i1 = theta_i + deltaTheta;

            //cerr << "\n theta_i: "<< theta_i;

            cellVertCoord(0,0) = centerPosit(0) + rad_j  * cos(theta_i1);
            cellVertCoord(0,1) = centerPosit(1) + rad_j  * sin(theta_i1);
            cellVertCoord(1,0) = centerPosit(0) + rad_j  * cos(theta_i);
            cellVertCoord(1,1) = centerPosit(1) + rad_j  * sin(theta_i);
            cellVertCoord(2,0) = centerPosit(0) + rad_j1 * cos(theta_i);
            cellVertCoord(2,1) = centerPosit(1) + rad_j1 * sin(theta_i);
            cellVertCoord(3,0) = centerPosit(0) + rad_j1 * cos(theta_i1);
            cellVertCoord(3,1) = centerPosit(1) + rad_j1 * sin(theta_i1);

            cells[k] = new QuadCell(cellVertCoord); 
            //cerr << "\ncreating cells Cell " << k << " :" << cells[k];
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
 
void CircPatch::Print(ostream &s, int flag) const
{
   s << "\nPatch Type: CircPatch";
   s << "\nMaterial Id: " << matID;
   s << "\nNumber of subdivisions in the radial direction: " << nDivRad;
   s << "\nNumber of subdivisions in the circunferential direction: " << nDivCirc;
   s << "\nCenter Position: " << centerPosit;
   s << "\nInternal Radius: " << intRad << "\tExternal Radius: " << extRad;
   s << "\nInitial Angle: " << initAng << "\tFinal Angle: " << finalAng;
}


ostream &operator<<(ostream &s, CircPatch &CircPatch)
{
   CircPatch.Print(s);
   return s;
}
