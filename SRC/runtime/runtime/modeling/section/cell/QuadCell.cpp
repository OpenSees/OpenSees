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
// File: QuadCell.C
// Written by Remo M. de Souza
// December 1998
//
#include <MatrixND.h>
#include <VectorND.h>
#include <QuadCell.h>
using OpenSees::MatrixND;
using OpenSees::VectorND;
// QuadCell::QuadCell() : vertCoord{}, location{} {}

QuadCell::QuadCell(const MatrixND<4,2>& vertCoord)
 //: vertCoord(vertexCoords)
{

  // double x0 = vertCoord(0, 0);
  // double y0 = vertCoord(0, 1);
  // double x1 = vertCoord(1, 0);
  // double y1 = vertCoord(1, 1);
  // double x2 = vertCoord(2, 0);
  // double y2 = vertCoord(2, 1);
  // double x3 = vertCoord(3, 0);
  // double y3 = vertCoord(3, 1);

  // double area = ((x2 - x1) * (y0 - y1) - (x0 - x1) * (y2 - y1) + (x0 - x3) * (y2 - y3) -
  //                (x2 - x3) * (y0 - y3)) /
  //               2.0;


  area = 0;

  for (int i = 0; i < 4; i++) {
    int i1 = (i + 1) % 4;
    double yi     = vertCoord(i, 0);
    double zi     = vertCoord(i, 1);
    double yi1    = vertCoord(i1, 0);
    double zi1    = vertCoord(i1, 1);

    area += (zi1 - zi) * (yi1 + yi);
  }
  area /= 2.0;

  //
  // centroid
  //

  double CGy = 0.0, CGz = 0.0;

  for (int i = 0; i < 4; i++) {
    int i1 = (i + 1) % 4;

    double yi  = vertCoord(i, 0);
    double zi  = vertCoord(i, 1);
    double yi1 = vertCoord(i1, 0);
    double zi1 = vertCoord(i1, 1);

    double dyi = yi1 - yi;
    double dzi = zi1 - zi;

    double integ = yi * zi + (yi * dzi + zi * dyi) / 2.0 + dyi * dzi / 3.0;

    CGy -= dyi * integ;
    CGz += dzi * integ;
  }

  CGy /= area;
  CGz /= area;

  location[0] = CGy;
  location[1] = CGz;

}
