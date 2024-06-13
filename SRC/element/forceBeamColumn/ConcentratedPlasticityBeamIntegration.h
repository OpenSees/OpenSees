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
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ConcentratedPlasticityBeamIntegration.h,v $

#ifndef ConcentratedPlasticityBeamIntegration_h
#define ConcentratedPlasticityBeamIntegration_h

#include <BeamIntegration.h>

#include <Vector.h>

class Channel;
class FEM_ObjectBroker;

class ConcentratedPlasticityBeamIntegration : public BeamIntegration
{
 public:
  ConcentratedPlasticityBeamIntegration();
  ~ConcentratedPlasticityBeamIntegration();
  
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

};

#endif
