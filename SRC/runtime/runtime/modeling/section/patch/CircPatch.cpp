//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Adapted from CircPatch.C
// Written by Remo M. de Souza
// December 1998
//
#include <CircPatch.h>
#include <CircSectionCell.h>
#include <Matrix.h>
#include <OPS_Stream.h>
#include <Patch.h>
#include <math.h>
#include <string>

CircPatch::CircPatch(int materialID, int numSubdivCircunf, int numSubdivRadial,
                     const VectorND<2>& centerPosition, double internRadius, double externRadius,
                     double initialAngle, double finalAngle)
 : matID(materialID),
   nDivCirc(numSubdivCircunf),
   nDivRad(numSubdivRadial),
   centerPosit(centerPosition),
   intRad(internRadius),
   extRad(externRadius),
   initAng(initialAngle),
   finalAng(finalAngle)
{
}

CircPatch::~CircPatch() {}

int
CircPatch::getMaterialID() const
{
  return matID;
}

int
CircPatch::getNumCells() const
{
  return nDivCirc * nDivRad;
}

Cell**
CircPatch::getCells() const
{
  double pi = acos(-1.0);
  double deltaRad, deltaTheta;
  double initAngRadians, finalAngRadians;
  double theta_i, theta_i1;
  Matrix cellVertCoord(4, 2);

  int numCells;
  Cell** cells;

  if (nDivRad > 0 && nDivCirc > 0) {
    numCells = this->getNumCells();

    cells = new Cell*[numCells];

    initAngRadians  = pi * initAng / 180.0;
    finalAngRadians = pi * finalAng / 180.0;

    deltaRad   = (extRad - intRad) / nDivRad;
    deltaTheta = (finalAngRadians - initAngRadians) / nDivCirc;

    int k = 0;
    for (int j = 0; j < nDivRad; j++) {
      double rad_j  = intRad + deltaRad * j;
      double rad_j1 = rad_j + deltaRad;

      for (int i = 0; i < nDivCirc; i++) {

        theta_i = initAngRadians + deltaTheta * i;

        theta_i1 = theta_i + deltaTheta / 2.0;
        cells[k] = new CircSectionCell(rad_j, rad_j1, deltaTheta, theta_i1, centerPosit(0), centerPosit(1));

        k++;
      }
    }
  } else
    return 0;

  return cells;
}
