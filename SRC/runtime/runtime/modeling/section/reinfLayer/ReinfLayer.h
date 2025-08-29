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
// File: ReinfLayer.h
// Written by Remo M. de Souza
// December 1998
//
#ifndef ReinfLayer_h
#define ReinfLayer_h

#include <vector>
#include <Cell.h>

class ReinfLayer {
public:
  ReinfLayer(int material, double area) : material(material), area(area) {}
  virtual ~ReinfLayer() {};

  int getMaterialID() const {
    return material;
  };

  virtual int getNumReinfBars() const = 0;
  virtual std::vector<Cell> getReinfBars() const = 0;

protected:
  double getCellArea() const {
    return area;
  };
private:
  int material;
  double area;
};

#endif
