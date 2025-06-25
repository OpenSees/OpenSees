/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// File: CircPatch.C
// Written by Remo M. de Souza
// December 1998
//
#include <math.h>
#include <Matrix.h>
#include <Patch.h>
#include <CircPatch.h>
#include <QuadCell.h>
#include <CircSectionCell.h>



CircPatch::CircPatch(int materialID, int numSubdivCircunf, int numSubdivRadial,
                     const VectorND<2> &centerPosition, double internRadius, 
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

VectorND<2>
CircPatch::getCenterPosition() const
{
   return centerPosit;
}

int CircPatch::getNumCells() const
{
   return nDivCirc * nDivRad;
}

Cell **
CircPatch::getCells() const
{
   double pi = acos(-1.0);
   double deltaRad, deltaTheta; 
   double initAngRadians, finalAngRadians;
   double rad_j, rad_j1, theta_i, theta_i1;
   Matrix cellVertCoord(4,2);

   int    numCells;
   Cell   **cells;

   if (nDivRad > 0  && nDivCirc > 0)
   {
      numCells  = this->getNumCells();

      cells = new Cell* [numCells];

      initAngRadians  = pi * initAng  / 180.0;
      finalAngRadians = pi * finalAng / 180.0;

      deltaRad   = (extRad - intRad) / nDivRad;
      deltaTheta = (finalAngRadians - initAngRadians) / nDivCirc;

      int k = 0;
      for (int j = 0; j < nDivRad; j++)
      {
         rad_j  = intRad + deltaRad*j;
         rad_j1 = rad_j + deltaRad;

         for (int i = 0; i < nDivCirc; i++)
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
CircPatch::getCopy() const
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
// s << "\nCenter Position: " << Vector(centerPosit);
   s << "\nInternal Radius: " << intRad << "\tExternal Radius: " << extRad;
   s << "\nInitial Angle: " << initAng << "\tFinal Angle: " << finalAng;
}

