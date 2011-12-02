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
// $Date: 2006-08-11 18:32:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/SectionIntegration.cpp,v $

#include <SectionIntegration.h>
#include <Matrix.h>

SectionIntegration::SectionIntegration(int classTag):
  MovableObject(classTag)
{
  // Nothing to do
}

SectionIntegration::~SectionIntegration()
{
  // Nothing to do
}

void
SectionIntegration::getLocationsDeriv(int nFibers, double *dptsdh)
{
  for (int i = 0; i < nFibers; i++)
    dptsdh[i] = 0.0;
}

void
SectionIntegration::getWeightsDeriv(int nFibers, double *dwtsdh)
{
  for (int i = 0; i < nFibers; i++)
    dwtsdh[i] = 0.0;
}
