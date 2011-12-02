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
// $Date: 2003-03-15 00:09:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/BeamIntegration.cpp,v $

#include <BeamIntegration.h>

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
