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
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/HingeRadauBeamIntegration2d.h,v $

#ifndef HingeRadauBeamIntegration2d_h
#define HingeRadauBeamIntegration2d_h

#include <BeamIntegration.h>

class Matrix;
class ElementalLoad;
class Channel;
class FEM_ObjectBroker;

class HingeRadauBeamIntegration2d : public BeamIntegration
{
 public:
  HingeRadauBeamIntegration2d(double E, double A, double I,
			      double lpI, double lpJ);
  HingeRadauBeamIntegration2d();
  ~HingeRadauBeamIntegration2d();
  
  void getSectionLocations(int numSections, double L, double *xi);
  void getSectionWeights(int numSections, double L, double *wt);
  
  void addElasticDeformations(ElementalLoad *theLoad, double loadFactor,
			      double L, double *v0);
  int addElasticFlexibility(double L, Matrix &fe);

  double getTangentDriftI(double L, double LI, double q2, double q3);
  double getTangentDriftJ(double L, double LI, double q2, double q3);

  BeamIntegration *getCopy(void);

  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  int setParameter(const char **argv, int argc, Information &info);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);

 private:
  double E;
  double A;
  double I;
  
  double lpI;
  double lpJ;
};

#endif
