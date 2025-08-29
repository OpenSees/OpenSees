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
// File: QuadPatch.C
// Written by Remo M. de Souza
// December 1998
//
#include <Matrix.h>
#include <Vector.h>
#include <VectorND.h>
#include <OPS_Stream.h>
#include <Patch.h>
#include <QuadCell.h>
#include <QuadPatch.h>
#include <string>

using OpenSees::VectorND;

QuadPatch::QuadPatch(int materialID, int numSubdivIJ, int numSubdivJK, 
                     const MatrixND<4,2>& vertexCoords)
 : matID(materialID), nDivIJ(numSubdivIJ), nDivJK(numSubdivJK), vertCoord(vertexCoords)

{
}

QuadPatch::~QuadPatch() {}

int
QuadPatch::getMaterialID() const
{
  return matID;
}

int
QuadPatch::getNumCells() const
{
  return nDivIJ * nDivJK;
}

Cell**
QuadPatch::getCells() const
{
  VectorND<4> N;
  double xi, eta;
  int numCells;
  Cell** cells;

  if (nDivIJ > 0 && nDivJK > 0) {
    numCells = this->getNumCells();

    cells = new Cell*[numCells];

    if (!cells)
      return 0;

    double deltaXi  = 2.0 / nDivIJ;
    double deltaEta = 2.0 / nDivJK;

    int k = 0;
    for (int j = 0; j < nDivJK; j++)
      for (int i = 0; i < nDivIJ; i++) {
        // compute natural coordinates

        MatrixND<4,2> cellVertCoord;
        cellVertCoord(0, 0) = -1.0 + deltaXi * i;
        cellVertCoord(0, 1) = -1.0 + deltaEta * j;
        cellVertCoord(1, 0) = -1.0 + deltaXi * (i + 1);
        cellVertCoord(1, 1) = cellVertCoord(0, 1);
        cellVertCoord(2, 0) = cellVertCoord(1, 0);
        cellVertCoord(2, 1) = -1.0 + deltaEta * (j + 1);
        cellVertCoord(3, 0) = cellVertCoord(0, 0);
        cellVertCoord(3, 1) = cellVertCoord(2, 1);

        // map to cartesian coordinates using bilinear
        // shape functions

        for (int r = 0; r < 4; r++) {
          xi  = cellVertCoord(r, 0);
          eta = cellVertCoord(r, 1);

          N(0) = (1.0 - xi) * (1.0 - eta) / 4.0;
          N(1) = (1.0 + xi) * (1.0 - eta) / 4.0;
          N(2) = (1.0 + xi) * (1.0 + eta) / 4.0;
          N(3) = (1.0 - xi) * (1.0 + eta) / 4.0;

          cellVertCoord(r, 0) = 0.0;
          cellVertCoord(r, 1) = 0.0;

          for (int s = 0; s < 4; s++) {
            cellVertCoord(r, 0) += N(s) * vertCoord(s, 0);
            cellVertCoord(r, 1) += N(s) * vertCoord(s, 1);
          }
        }

        cells[k] = new QuadCell(cellVertCoord);
        k++;
      }
  } else
    return nullptr;

  return cells;
}

