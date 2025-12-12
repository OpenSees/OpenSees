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
// File: QuadPatch.h
// Written by Remo M. de Souza
// December 1998
//
#ifndef QuadPatch_h
#define QuadPatch_h

#include <MatrixND.h>
#include <Patch.h>
#include <Vector.h>

class Cell;
class Matrix;
using OpenSees::MatrixND;

class QuadPatch : public Patch {
public:
  QuadPatch();
  QuadPatch(int materialID, int numSubdivIJ, int numSubdivJK, const MatrixND<4,2>& vertexCoords);

  ~QuadPatch();

  int getMaterialID() const;
  int getNumCells() const;
  Cell** getCells() const;

protected:
private:
  int matID;
  int nDivIJ, nDivJK;
  MatrixND<4,2> vertCoord;
};

#endif
