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
// $Date: 2006-01-17 21:24:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/BeamIntegration.cpp,v $

#include <BeamIntegration.h>
#include <Matrix.h>

BeamIntegration::BeamIntegration(int classTag):
  MovableObject(classTag)
{
  // Nothing to do
}

BeamIntegration::~BeamIntegration()
{
  // Nothing to do
}

int 
BeamIntegration::setParameter(const char **argv, int argc, Information &info)
{
  return 0;
}

int
BeamIntegration::updateParameter(int parameterID, Information &info)
{
  return 0;
}

int
BeamIntegration::activateParameter(int parameterID)
{
  return 0;
}

void
BeamIntegration::getLocationsDeriv(int nIP, double L, double dLdh,
				   double *dptsdh)
{
  for (int i = 0; i < nIP; i++)
    dptsdh[i] = 0.0;
}

void
BeamIntegration::getWeightsDeriv(int nIP, double L, double dLdh,
				 double *dwtsdh)
{
  for (int i = 0; i < nIP; i++)
    dwtsdh[i] = 0.0;
}
