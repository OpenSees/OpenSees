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
// File: StraightReinfLayer.h
// Written by Remo M. de Souza
// December 1998
//
#ifndef StraightReinfLayer_h
#define StraightReinfLayer_h

#include <ReinfLayer.h>
#include <VectorND.h>
#include <Cell.h>

using OpenSees::VectorND;

class StraightReinfLayer : public ReinfLayer {
public:
  StraightReinfLayer(int materialID, int numReinfBars, double reinfBarArea,
                     const VectorND<2>& initialPosition, const VectorND<2>& finalPosition);

  virtual ~StraightReinfLayer() {};

  int getNumReinfBars() const;
  std::vector<Cell> getReinfBars() const;

protected:
private:
  int nReinfBars;
  VectorND<2> initPosit;
  VectorND<2> finalPosit;
};

#endif
