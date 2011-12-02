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

// $Revision: 1.2 $
// $Date: 2008-12-03 23:43:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/RegularizedHingeIntegration.h,v $

// Theory Reference
// ----------------
// Scott, M.H. and Hamutcuoglu, O.M. "Numerically consistent regularization
// of force-based frame elements." International Journal for Numerical
// Methods in Engineering. http://dx.doi.org/10.1002/nme.2386

#ifndef RegularizedHingeIntegration_h
#define RegularizedHingeIntegration_h

#include <BeamIntegration.h>
#include <math.h>

class Matrix;
class ElementalLoad;
class Channel;
class FEM_ObjectBroker;

class RegularizedHingeIntegration : public BeamIntegration
{
 public:
  RegularizedHingeIntegration(BeamIntegration &bi, 
			      double lpI, double lpJ,
			      double epsI, double epsJ);
  RegularizedHingeIntegration();
  ~RegularizedHingeIntegration();
  
  void getSectionLocations(int numSections, double L, double *xi);
  void getSectionWeights(int numSections, double L, double *wt);
  
  BeamIntegration *getCopy(void);

  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);

  void getLocationsDeriv(int nIP, double L, double dLdh, double *dptsdh);
  void getWeightsDeriv(int nIP, double L, double dLdh, double *dwtsdh);

  void Print(OPS_Stream &s, int flag = 0);

 private:
  double lpI;
  double lpJ;

  double epsI;
  double epsJ;

  BeamIntegration *beamInt;

  double *wf;

  int parameterID;
};

#endif
