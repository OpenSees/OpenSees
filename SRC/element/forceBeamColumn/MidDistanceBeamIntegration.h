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

// $Revision: 1.1 $
// $Date: 2007-07-12 20:49:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/MidDistanceBeamIntegration.h,v $

#ifndef MidDistanceBeamIntegration_h
#define MidDistanceBeamIntegration_h

#include <BeamIntegration.h>

#include <Vector.h>

class Channel;
class FEM_ObjectBroker;

class MidDistanceBeamIntegration : public BeamIntegration
{
 public:
  MidDistanceBeamIntegration(int nIP, const Vector &pt);
  MidDistanceBeamIntegration();
  ~MidDistanceBeamIntegration();
  
  void getSectionLocations(int numSections, double L, double *xi);
  void getSectionWeights(int numSections, double L, double *wt);

  BeamIntegration *getCopy(void);

  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  int setParameter(const char **argv, int argc, Information &info);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);

  void Print(OPS_Stream &s, int flag = 0);  

  void getLocationsDeriv(int nIP, double L, double *dptsdh);
  void getWeightsDeriv(int nIP, double L, double *dwtsdh);

 private:
  Vector pts;
  Vector wts;

  int parameterID;
};

#endif
