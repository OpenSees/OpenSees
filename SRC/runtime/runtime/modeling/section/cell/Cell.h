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
// File: Cell.h
//
// Written by Remo M. de Souza
// December 1998

#ifndef Cell_h
#define Cell_h
#include <VectorND.h>

using OpenSees::VectorND;
class Cell {
public:
  Cell() : area(0.0) {}
  Cell(int mat, double area, const VectorND<2>& loc) : area(area), location(loc)
  {
  }
  // reinforcing bar inquiring functions

  double getArea() const {return area;};
  double getdValue() const {return 0.0;};
  const VectorND<2>& getPosition() const {return location;};
  int getMaterial() const {return matID;};

protected:
  double area;
  int matID;
  VectorND<2> location;
private:
};

#endif
