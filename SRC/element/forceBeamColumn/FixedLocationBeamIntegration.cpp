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
// $Date: 2007-10-12 21:03:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/FixedLocationBeamIntegration.cpp,v $

#include <FixedLocationBeamIntegration.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>

FixedLocationBeamIntegration::FixedLocationBeamIntegration(int nIP,
							   const Vector &pt):
  BeamIntegration(BEAM_INTEGRATION_TAG_FixedLocation),
  pts(nIP), wts(nIP)
{
  for (int i = 0; i < nIP; i++) {
    if (pt(i) < 0.0 || pt(i) > 1.0)
      opserr << "FixedLocationBeamIntegration::FixedLocationBeamIntegration -- point lies outside [0,1]" << endln;

    pts(i) = pt(i);
  }

  Vector R(nIP);
  for (int i = 0; i < nIP; i++)
    R(i) = 1.0/(i+1);

  Matrix J(nIP,nIP);
  for (int i = 0; i < nIP; i++)
    for (int j = 0; j < nIP; j++)
      J(i,j) = pow(pts(j),i);

  J.Solve(R, wts);
}

FixedLocationBeamIntegration::FixedLocationBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_FixedLocation)
{
 
}

FixedLocationBeamIntegration::~FixedLocationBeamIntegration()
{
  // Nothing to do
}

void
FixedLocationBeamIntegration::getSectionLocations(int numSections,
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
FixedLocationBeamIntegration::getSectionWeights(int numSections,
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
FixedLocationBeamIntegration::getCopy(void)
{
  int nIP = pts.Size();

  return new FixedLocationBeamIntegration(nIP, pts);
}

int
FixedLocationBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
FixedLocationBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				       FEM_ObjectBroker &theBroker)
{
  return -1;
}

int
FixedLocationBeamIntegration::setParameter(const char **argv,
					   int argc, Information &info)
{
  return -1;
}

int
FixedLocationBeamIntegration::updateParameter(int parameterID,
					      Information &info)
{
  // Does nothing for now -- MHS
  return 0;
}

int
FixedLocationBeamIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  // For Terje to do
  return 0;
}

void
FixedLocationBeamIntegration::Print(OPS_Stream &s, int flag)
{
  s << "FixedLocation" << endln;
  s << " Points: " << pts;
  s << " Weights: " << wts;
  double sum = 0.0;
  int N = wts.Size();
  for (int i = 0; i < N; i++)
    sum += fabs(wts(i));
  s << " Condition Number: " << sum << endln;
}

void 
FixedLocationBeamIntegration::getLocationsDeriv(int numSections, double L,
						double dLdh, double *dptsdh)
{
  for (int i = 0; i < numSections; i++)
    dptsdh[i] = 0.0;

  return;
}

void
FixedLocationBeamIntegration::getWeightsDeriv(int numSections, double L,
					      double dLdh, double *dwtsdh)
{
  for (int i = 0; i < numSections; i++)
    dwtsdh[i] = 0.0;

  return;
}
