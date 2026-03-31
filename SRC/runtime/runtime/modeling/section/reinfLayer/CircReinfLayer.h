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
// File: CircReinfLayer.h
// Written by Remo M. de Souza
// December 1998

#ifndef CircReinfLayer_h
#define CircReinfLayer_h

#include <ReinfLayer.h>
#include <VectorND.h>
#include <Cell.h>

using OpenSees::VectorND;

class CircReinfLayer : public ReinfLayer {
public:
  // Constructor for an arc
  CircReinfLayer(int material, int n, double area,
                 const VectorND<2>& center, double radius, 
                 double initialAngle,
                 double finalAngle);
  // Constructor for full circle
  CircReinfLayer(int material, int n, double area,
                 const VectorND<2>& center, double radius);

  virtual ~CircReinfLayer() {};

  int getNumReinfBars() const;
  std::vector<Cell> getReinfBars() const;


protected:
private:
  int nReinfBars;
  VectorND<2> centerPosit;
  double arcRad;
  double initAng;
  double finalAng;
};

#endif
