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

// $Revision: 1.3 $
// $Date: 2003-02-21 22:27:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/UserDefinedBeamIntegration.cpp,v $

#include <UserDefinedBeamIntegration.h>

#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

UserDefinedBeamIntegration::UserDefinedBeamIntegration(int nIP,
						       const Vector &pt,
						       const Vector &wt):
  BeamIntegration(BEAM_INTEGRATION_TAG_UserDefined),
  pts(nIP), wts(nIP)
{
  for (int i = 0; i < nIP; i++) {
    if (pt(i) < 0.0 || pt(i) > 1.0)
      opserr << "UserDefinedBeamIntegration::UserDefinedBeamIntegration -- point lies outside [0,1]" << endln;
    if (wt(i) < 0.0 || wt(i) > 1.0)
      opserr << "UserDefinedBeamIntegration::UserDefinedBeamIntegration -- weight lies outside [0,1]" << endln;
    pts(i) = pt(i);
    wts(i) = wt(i);
  }
}

UserDefinedBeamIntegration::UserDefinedBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_UserDefined)
{
 
}

UserDefinedBeamIntegration::~UserDefinedBeamIntegration()
{
  // Nothing to do
}

void
UserDefinedBeamIntegration::getSectionLocations(int numSections,
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
UserDefinedBeamIntegration::getSectionWeights(int numSections,
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
UserDefinedBeamIntegration::getCopy(void)
{
  int nIP = pts.Size();

  return new UserDefinedBeamIntegration(nIP, pts, wts);
}

int
UserDefinedBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
UserDefinedBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				     FEM_ObjectBroker &theBroker)
{
  return -1;
}
