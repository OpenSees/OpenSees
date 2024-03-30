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
// $Date: 2007-01-25 19:53:17 $
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
SectionIntegration::getLocationsDeriv(int nFibers, double *dyidh, double *dzidh)
{
  for (int i = 0; i < nFibers; i++)
    dyidh[i] = 0.0;

  if (dzidh != 0) {
    for (int i = 0; i < nFibers; i++)
      dzidh[i] = 0.0;
  }
}

void
SectionIntegration::getWeightsDeriv(int nFibers, double *dwtdh)
{
  for (int i = 0; i < nFibers; i++)
    dwtdh[i] = 0.0;
}
