//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Adapted from CircPatch.h
// Written by Remo M. de Souza
// December 1998
//
#ifndef CircPatch_h
#define CircPatch_h

#include <Patch.h>
#include <Vector.h>
#include <VectorND.h>
using OpenSees::VectorND;

class Cell;
class Matrix;

class CircPatch : public Patch {
public:
  CircPatch(int material,
            int numSubdivCircunf, int numSubdivRadial,
            const VectorND<2>& centerPosition, 
            double internRadius, double externRadius,
            double initialAngle, double finalAngle);

  ~CircPatch();

  int getMaterialID() const;
  int getNumCells() const;
  Cell** getCells() const;

protected:
private:
  int matID;
  int nDivCirc, nDivRad;
  const VectorND<2> centerPosit;
  double intRad, extRad;
  double initAng, finalAng;
};

#endif
