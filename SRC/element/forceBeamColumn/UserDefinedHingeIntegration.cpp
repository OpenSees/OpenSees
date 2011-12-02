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
// $Date: 2006-01-18 21:58:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/UserDefinedHingeIntegration.cpp,v $

#include <UserDefinedHingeIntegration.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>

UserDefinedHingeIntegration::UserDefinedHingeIntegration(int npL,
							 const Vector &ptL,
							 const Vector &wtL,
							 int npR,
							 const Vector &ptR,
							 const Vector &wtR):
  BeamIntegration(BEAM_INTEGRATION_TAG_UserHinge),
  ptsL(npL), wtsL(npL), ptsR(npR), wtsR(npR)
{
  int i;
  for (i = 0; i < npL; i++) {
    if (ptL(i) < 0.0 || ptL(i) > 1.0)
      opserr << "UserDefinedHingeIntegration::UserDefinedHingeIntegration -- point lies outside [0,1]" << endln;
    if (wtL(i) < 0.0 || wtL(i) > 1.0)
      opserr << "UserDefinedHingeIntegration::UserDefinedHingeIntegration -- weight lies outside [0,1]" << endln;
    ptsL(i) = ptL(i);
    wtsL(i) = wtL(i);
  }

  for (i = 0; i < npR; i++) {
    if (ptR(i) < 0.0 || ptR(i) > 1.0)
      opserr << "UserDefinedHingeIntegration::UserDefinedHingeIntegration -- point lies outside [0,1]" << endln;
    if (wtR(i) < 0.0 || wtR(i) > 1.0)
      opserr << "UserDefinedHingeIntegration::UserDefinedHingeIntegration -- weight lies outside [0,1]" << endln;
    ptsR(i) = ptR(i);
    wtsR(i) = wtR(i);
  }
}

UserDefinedHingeIntegration::UserDefinedHingeIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_UserHinge)
{

}

UserDefinedHingeIntegration::~UserDefinedHingeIntegration()
{
  // Nothing to do
}

void
UserDefinedHingeIntegration::getSectionLocations(int numSections,
						 double L, double *xi)
{
  int npL = ptsL.Size();
  int npR = ptsR.Size();

  double lpI = 0.0;
  double lpJ = 0.0;
  int i, j;
  for (i = 0; i < npL; i++) {
    xi[i] = ptsL(i);
    lpI += wtsL(i);
  }
  for (j = 0; j < npR; j++, i++) {
    xi[i] = ptsR(j);
    lpJ += wtsR(j);
  }

  double alpha = 0.5-0.5*(lpI+lpJ);
  double beta  = 0.5+0.5*(lpI-lpJ);
  xi[i++] = alpha*(-1/sqrt(3.0)) + beta;
  xi[i++] = alpha*(1/sqrt(3.0)) + beta;

  for ( ; i < numSections; i++)
    xi[i] = 0.0;
}

void
UserDefinedHingeIntegration::getSectionWeights(int numSections,
					       double L, double *wt)
{
  int npL = wtsL.Size();
  int npR = wtsR.Size();

  double lpI = 0.0;
  double lpJ = 0.0;
  int i, j;
  for (i = 0; i < npL; i++) {
    wt[i] = wtsL(i);
    lpI += wtsL(i);
  }
  for (j = 0; j < npR; j++, i++) {
    wt[i] = wtsR(j);
    lpJ += wtsR(j);
  }

  double oneOverL = 1.0/L;
  wt[i++] = 0.5-0.5*(lpI+lpJ);
  wt[i++] = 0.5-0.5*(lpI+lpJ);

  for ( ; i < numSections; i++)
    wt[i] = 1.0;
}

BeamIntegration*
UserDefinedHingeIntegration::getCopy(void)
{
  int npL = ptsL.Size();
  int npR = ptsR.Size();

  return new UserDefinedHingeIntegration(npL, ptsL, wtsL,
					 npR, ptsR, wtsR);
}

int
UserDefinedHingeIntegration::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
UserDefinedHingeIntegration::recvSelf(int cTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
UserDefinedHingeIntegration::Print(OPS_Stream &s, int flag)
{
  s << "UserHinge" << endln;
  s << " Points hinge I: " << ptsL;
  s << " Weights hinge I: " << wtsL;
  s << " Points hinge J: " << ptsR;
  s << " Weights hinge J: " << wtsR;
}
