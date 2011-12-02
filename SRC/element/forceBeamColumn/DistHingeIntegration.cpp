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
// $Date: 2006-01-18 22:11:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/DistHingeIntegration.cpp,v $

#include <DistHingeIntegration.h>
#include <ElementalLoad.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>

DistHingeIntegration::DistHingeIntegration(double lpi,
					   double lpj,
					   BeamIntegration &bi):
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeMidpoint),
  lpI(lpi), lpJ(lpj), beamInt(0)
{
  beamInt = bi.getCopy();
  if (beamInt == 0) {
    opserr << "DistHingeIntegration::DistHingeIntegration -- failed to get copy of BeamIntegration" << endln;
  }
}

DistHingeIntegration::DistHingeIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeMidpoint),
  lpI(0.0), lpJ(0.0), beamInt(0)
{

}

DistHingeIntegration::~DistHingeIntegration()
{
  if (beamInt != 0)
    delete beamInt;
}

void
DistHingeIntegration::getSectionLocations(int numSections, double L,
					  double *xi)
{
  int numPerHinge = numSections/2;

  beamInt->getSectionLocations(numPerHinge, L, xi);

  double betaI = lpI/L;
  double betaJ = lpJ/L;
  
  // Map from [0,L] to [L-lpJ,L]
  for (int i = 0; i < numPerHinge; i++) {
    xi[numSections-1-i] = 1.0-betaJ*xi[i];
    xi[i] *= betaI;
  }

  opserr << "DistHingeIntegration::getSectionLocations -- implementation for interior not yet finished" << endln;
}

void
DistHingeIntegration::getSectionWeights(int numSections, double L,
					double *wt)
{
  int numPerHinge = numSections/2;

  beamInt->getSectionWeights(numPerHinge, L, wt);

  double betaI = lpI/L;
  double betaJ = lpJ/L;
  
  // Map from [0,lpI] to [L-lpJ,L]
  for (int i = 0; i < numPerHinge; i++) {
    wt[numSections-1-i] = betaJ*wt[i];
    wt[i] *= betaI;
  }

  opserr << "DistHingeIntegration::getSectionWeights -- implementation for interior not yet finished" << endln;
}

BeamIntegration*
DistHingeIntegration::getCopy(void)
{
  return new DistHingeIntegration(lpI, lpJ, *beamInt);
}

int
DistHingeIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(2);

  data(0) = lpI;
  data(1) = lpJ;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "DistHingeIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
DistHingeIntegration::recvSelf(int cTag, Channel &theChannel,
				 FEM_ObjectBroker &theBroker)
{
  static Vector data(2);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "DistHingeIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  lpI = data(0);
  lpJ = data(1);

  return 0;
}

int
DistHingeIntegration::setParameter(const char **argv,
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
  else 
    return -1;
}

int
DistHingeIntegration::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    lpI = info.theDouble;
    return 0;
  case 2:
    lpJ = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
DistHingeIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  // For Terje to do
  return 0;
}

void
DistHingeIntegration::Print(OPS_Stream &s, int flag)
{
  s << "DistHinge" << endln;
  s << " lpI = " << lpI;
  s << " lpJ = " << lpJ << endln;

  beamInt->Print(s, flag);

  return;
}

void 
DistHingeIntegration::getLocationsDeriv(int numSections, double L,
					  double dLdh, double *dptsdh)
{
  int numPerHinge = numSections/2;

  double oneOverL = 1/L;
  double betaI = lpI*oneOverL;
  double betaJ = lpJ*oneOverL;

  beamInt->getSectionLocations(numPerHinge, L, dptsdh);

  if (parameterID == 1) { // lpI
    for (int i = 0; i < numPerHinge; i++) {
      dptsdh[i] = oneOverL*dptsdh[i];
      dptsdh[numSections-1-i] = 0.0;
    }
  }
  else if (parameterID == 2) { // lpJ
    for (int i = 0; i < numPerHinge; i++) {
      dptsdh[numSections-1-i] = -oneOverL*dptsdh[i];
      dptsdh[i] = 0.0;
    }
  }
  else if (dLdh != 0.0) {
    for (int i = 0; i < numPerHinge; i++) {
      dptsdh[numSections-1-i] = betaJ*oneOverL*dLdh*dptsdh[i];
      dptsdh[i] = -betaI*oneOverL*dLdh*dptsdh[i];
    }
  }
  else {
    for (int i = 0; i < numSections; i++)
      dptsdh[i] = 0.0;
  }

  return;
}

void
DistHingeIntegration::getWeightsDeriv(int numSections, double L,
					double dLdh, double *dwtsdh)
{
  int numPerHinge = numSections/2;

  double oneOverL = 1/L;
  double betaI = lpI*oneOverL;
  double betaJ = lpJ*oneOverL;

  beamInt->getSectionWeights(numPerHinge, L, dwtsdh);

  if (parameterID == 1) { // lpI
    for (int i = 0; i < numPerHinge; i++) {
      dwtsdh[i] = oneOverL*dwtsdh[i];
      dwtsdh[numSections-1-i] = 0.0;
    }
  }
  else if (parameterID == 2) { // lpJ
    for (int i = 0; i < numPerHinge; i++) {
      dwtsdh[numSections-1-i] = oneOverL*dwtsdh[i];
      dwtsdh[i] = 0.0;
    }
  }
  else if (dLdh != 0.0) {
    for (int i = 0; i < numPerHinge; i++) {
      dwtsdh[numSections-1-i] = -betaJ*oneOverL*dLdh*dwtsdh[i];
      dwtsdh[i] = -betaI*oneOverL*dLdh*dwtsdh[i];
    }
  }
  else {
    for (int i = 0; i < numSections; i++)
      dwtsdh[i] = 0.0;
  }

  return;
}
