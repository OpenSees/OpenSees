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
// $Date: 2006-01-17 21:32:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/HingeRadauBeamIntegration.cpp,v $

#include <HingeRadauBeamIntegration.h>
#include <ElementalLoad.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>

HingeRadauBeamIntegration::HingeRadauBeamIntegration(double lpi,
						     double lpj):
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadau), lpI(lpi), lpJ(lpj)
{
  // Nothing to do
}

HingeRadauBeamIntegration::HingeRadauBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadau), lpI(0.0), lpJ(0.0)
{

}

HingeRadauBeamIntegration::~HingeRadauBeamIntegration()
{
  // Nothing to do
}

void
HingeRadauBeamIntegration::getSectionLocations(int numSections, double L,
					       double *xi)
{
  double oneOverL = 1.0/L;

  xi[0] = 0.0;
  xi[1] = 8.0/3*lpI*oneOverL;
  xi[4] = 1.0-8.0/3*lpJ*oneOverL;
  xi[5] = 1.0;

  double alpha = 0.5-2*(lpI+lpJ)*oneOverL;
  double beta  = 0.5+2*(lpI-lpJ)*oneOverL;
  xi[2] = alpha*(-1/sqrt(3.0)) + beta;
  xi[3] = alpha*(1/sqrt(3.0)) + beta;

  for (int i = 6; i < numSections; i++)
    xi[i] = 0.0;
}

void
HingeRadauBeamIntegration::getSectionWeights(int numSections, double L,
					     double *wt)
{
  double oneOverL = 1.0/L;

  wt[0] = lpI*oneOverL;
  wt[1] = 3*lpI*oneOverL;
  wt[4] = 3*lpJ*oneOverL;
  wt[5] = lpJ*oneOverL;

  wt[2] = 0.5-2*(lpI+lpJ)*oneOverL;
  wt[3] = 0.5-2*(lpI+lpJ)*oneOverL;

  for (int i = 6; i < numSections; i++)
    wt[i] = 1.0;
}

BeamIntegration*
HingeRadauBeamIntegration::getCopy(void)
{
  return new HingeRadauBeamIntegration(lpI, lpJ);
}

int
HingeRadauBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(2);

  data(0) = lpI;
  data(1) = lpJ;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "HingeRadauBeamIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
HingeRadauBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				    FEM_ObjectBroker &theBroker)
{
  static Vector data(2);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "HingeRadauBeamIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  lpI = data(0);
  lpJ = data(1);

  return 0;
}

int
HingeRadauBeamIntegration::setParameter(const char **argv,
					int argc, Information &info)
{
  if (strcmp(argv[0],"lpI") == 0) {
    info.theType = DoubleType;
    return 1;
  }
  else if (strcmp(argv[0],"lpJ") == 0) {
    info.theType = DoubleType;
    return 2;
  }
  else if (strcmp(argv[0],"lp") == 0) {
    info.theType = DoubleType;
    return 3;
  }
  else 
    return -1;
}

int
HingeRadauBeamIntegration::updateParameter(int parameterID,
					   Information &info)
{
  switch (parameterID) {
  case 1:
    lpI = info.theDouble;
    return 0;
  case 2:
    lpJ = info.theDouble;
    return 0;
  case 3:
    lpI = lpJ = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
HingeRadauBeamIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  // For Terje to do
  return 0;
}

void
HingeRadauBeamIntegration::Print(OPS_Stream &s, int flag)
{
  s << "HingeRadau" << endln;
  s << " lpI = " << lpI;
  s << " lpJ = " << lpJ << endln;

  return;
}

void 
HingeRadauBeamIntegration::getLocationsDeriv(int numSections,
					     double L, double dLdh,
					     double *dptsdh)
{
  double oneOverL = 1.0/L;

  for (int i = 0; i < numSections; i++)
    dptsdh[i] = 0.0;

  if (parameterID == 1) { // lpI
    dptsdh[1] = 8.0/3*oneOverL;
    dptsdh[2] = -2*(1.0-1/sqrt(3.0))*oneOverL + 4*oneOverL;
    dptsdh[3] = -2*(1.0+1/sqrt(3.0))*oneOverL + 4*oneOverL;
  }

  if (parameterID == 2) { // lpJ
    dptsdh[2] = -2*(1.0-1/sqrt(3.0))*oneOverL;
    dptsdh[3] = -2*(1.0+1/sqrt(3.0))*oneOverL;
    dptsdh[4] = -8.0/3*oneOverL;
  }

  if (parameterID == 3) { // lpI and lpJ
    dptsdh[1] = 8.0/3*oneOverL;
    dptsdh[2] = -4*(1.0-1/sqrt(3.0))*oneOverL + 4*oneOverL;
    dptsdh[3] = -4*(1.0+1/sqrt(3.0))*oneOverL + 4*oneOverL;
    dptsdh[4] = -8.0/3*oneOverL;
  }

  if (dLdh != 0.0) {
    // STILL TO DO
    opserr << "getPointsDeriv -- to do" << endln;
  }

  return;
}

void
HingeRadauBeamIntegration::getWeightsDeriv(int numSections,
					   double L, double dLdh,
					   double *dwtsdh)
{
  double oneOverL = 1.0/L;

  for (int i = 0; i < numSections; i++)
    dwtsdh[i] = 0.0;

  if (parameterID == 1) { // lpI
    dwtsdh[0] = oneOverL;
    dwtsdh[1] = 3*oneOverL;
    dwtsdh[2] = -2*oneOverL;
    dwtsdh[3] = -2*oneOverL;
  }

  if (parameterID == 2) { // lpJ
    dwtsdh[2] = -2*oneOverL;
    dwtsdh[3] = -2*oneOverL;
    dwtsdh[4] = 3*oneOverL;
    dwtsdh[5] = oneOverL;
  }

  if (parameterID == 3) { // lpI and lpJ
    dwtsdh[0] = oneOverL;
    dwtsdh[1] = 3*oneOverL;
    dwtsdh[2] = -4*oneOverL;
    dwtsdh[3] = -4*oneOverL;
    dwtsdh[4] = 3*oneOverL;
    dwtsdh[5] = oneOverL;
  }

  if (dLdh != 0.0) {
    dwtsdh[0] = -lpI*dLdh/(L*L);
    dwtsdh[5] = -lpJ*dLdh/(L*L);
    // STILL TO DO
    opserr << "getWeightsDeriv -- to do" << endln;
  }

  return;
}
