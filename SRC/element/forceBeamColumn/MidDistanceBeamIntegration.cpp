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
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/MidDistanceBeamIntegration.cpp,v $

#include <MidDistanceBeamIntegration.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>

MidDistanceBeamIntegration::MidDistanceBeamIntegration(int nIP,
						       const Vector &pt):
  BeamIntegration(BEAM_INTEGRATION_TAG_MidDistance),
  pts(nIP), wts(nIP)
{
  for (int i = 0; i < nIP; i++) {
    if (pt(i) < 0.0 || pt(i) > 1.0)
      opserr << "MidDistanceBeamIntegration::MidDistanceBeamIntegration -- point lies outside [0,1]" << endln;

    pts(i) = pt(i);
  }

  Vector mids(nIP-1);

  for (int i = 0; i < nIP-1; i++) {
    mids(i) = 0.5*(pts(i)+pts(i+1));
  }

  wts(0) = mids(0);
  wts(nIP-1) = 1.0-mids(nIP-2);
  for (int i = 1; i < nIP-1; i++) {
    wts(i) = mids(i)-mids(i-1);
  }
}

MidDistanceBeamIntegration::MidDistanceBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_MidDistance)
{
 
}

MidDistanceBeamIntegration::~MidDistanceBeamIntegration()
{
  // Nothing to do
}

void
MidDistanceBeamIntegration::getSectionLocations(int numSections,
						double L, double *xi)
{
  int nIP = pts.Size();

  int i;
  for (i = 0; i < nIP; i++)
    xi[i] = pts(i);
  for ( ; i < numSections; i++)
    xi[i] = 0.0;
}

void
MidDistanceBeamIntegration::getSectionWeights(int numSections,
					      double L, double *wt)
{
  int nIP = wts.Size();

  int i;
  for (i = 0; i < nIP; i++)
    wt[i] = wts(i);
  for ( ; i < numSections; i++)
    wt[i] = 1.0;
}

BeamIntegration*
MidDistanceBeamIntegration::getCopy(void)
{
  int nIP = pts.Size();

  return new MidDistanceBeamIntegration(nIP, pts);
}

int
MidDistanceBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
MidDistanceBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				     FEM_ObjectBroker &theBroker)
{
  return -1;
}

int
MidDistanceBeamIntegration::setParameter(const char **argv,
					 int argc, Information &info)
{
  return -1;
}

int
MidDistanceBeamIntegration::updateParameter(int parameterID,
					    Information &info)
{
  // Does nothing for now -- MHS
  return 0;
}

int
MidDistanceBeamIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  // For Terje to do
  return 0;
}

void
MidDistanceBeamIntegration::Print(OPS_Stream &s, int flag)
{
  s << "MidDistance" << endln;
  s << " Points: " << pts;
  s << " Weights: " << wts;
}

void 
MidDistanceBeamIntegration::getLocationsDeriv(int numSections,
					      double L, double *dptsdh)
{
  for (int i = 0; i < numSections; i++)
    dptsdh[i] = 0.0;

  return;
}

void
MidDistanceBeamIntegration::getWeightsDeriv(int numSections,
					    double L, double *dwtsdh)
{
  for (int i = 0; i < numSections; i++)
    dwtsdh[i] = 0.0;

  return;
}
