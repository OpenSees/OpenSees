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

// $Revision: 1.4 $
// $Date: 2003-05-12 23:44:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/BeamIntegration.h,v $

#ifndef BeamIntegration_h
#define BeamIntegration_h

#include <MovableObject.h>

class Matrix;
class ElementalLoad;
class Information;

class BeamIntegration : public MovableObject
{
 public:
  BeamIntegration(int classTag);
  virtual ~BeamIntegration();

  virtual void getSectionLocations(int nIP, double L, double *xi) = 0;
  virtual void getSectionWeights(int nIP, double L, double *wt) = 0;

  virtual void addElasticDeformations(ElementalLoad *theLoad,
				      double loadFactor,
				      double L, double *v0) {return;}
  // Return 0 if there is no elastic interior, -1 otherwise
  virtual int addElasticFlexibility(double L, Matrix &fe) {return 0;}

  virtual double getTangentDriftI(double L, double LI, double q2,
				  double q3, bool yAxis = false) {return 0.0;}
  virtual double getTangentDriftJ(double L, double LI, double q2,
				  double q3, bool yAxis = false) {return 0.0;}

  virtual BeamIntegration *getCopy(void) = 0;

  virtual int setParameter(const char **argv, int argc, Information &info);
  virtual int updateParameter(int parameterID, Information &info);
  virtual int activateParameter(int parameterID);
};

#endif
