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
// File: QuadCell.h
//
// Written by Remo M. de Souza
// December 1998

#ifndef QuadCell_h
#define QuadCell_h

#include <Cell.h>
#include <VectorND.h>
#include <MatrixND.h>

using OpenSees::MatrixND;
using OpenSees::VectorND;


class QuadCell : public Cell {
public:
  QuadCell(const MatrixND<4,2>& vertexCoords);

protected:
private:
};

#endif
