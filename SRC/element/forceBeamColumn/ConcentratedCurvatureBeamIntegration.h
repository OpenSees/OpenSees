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
**                                                                    **
** ****************************************************************** *
** This Beam Integration Module was Developed by:                     **
**   Silvia Mazzoni (silviamazzoni@yahoo.com)                         **
**   Michael H. Scott (michael@portwooddigital.com)                   **
**   pushed to OpenSees November 2022                                 **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2022-11-20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ConcentratedCurvatureBeamIntegration.h,v $

/*
 * Original General Reference on Integration

Scott, M. H. and G. L. Fenves. "Plastic Hinge Integration Methods for
Force-Based Beam-Column Elements." Journal of Structural Engineering,
132(2):244-252, February 2006.

 *
 */

#ifndef ConcentratedCurvatureBeamIntegration_h
#define ConcentratedCurvatureBeamIntegration_h

#include <BeamIntegration.h>

class Matrix;
class ElementalLoad;
class Channel;
class FEM_ObjectBroker;

class ConcentratedCurvatureBeamIntegration : public BeamIntegration
{
 public:
  ConcentratedCurvatureBeamIntegration(double lpI, double lpJ);
  ConcentratedCurvatureBeamIntegration();
  ~ConcentratedCurvatureBeamIntegration();
  
  void getSectionLocations(int numSections, double L, double *xi);
  void getSectionWeights(int numSections, double L, double *wt);
  
  BeamIntegration *getCopy(void);

  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);

  void Print(OPS_Stream &s, int flag = 0);

  void getLocationsDeriv(int nIP, double L, double dLdh, double *dptsdh);
  void getWeightsDeriv(int nIP, double L, double dLdh, double *dwtsdh);

 private:
  double lpI;
  double lpJ;

  int parameterID;
};

#endif
