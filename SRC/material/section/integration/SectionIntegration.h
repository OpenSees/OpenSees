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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/integration/SectionIntegration.h,v $

#ifndef SectionIntegration_h
#define SectionIntegration_h

#include <OPS_Globals.h>
#include <MovableObject.h>

class Information;

class SectionIntegration : public MovableObject
{
 public:
  SectionIntegration(int classTag);
  virtual ~SectionIntegration();

  virtual int getNumFibers(void) = 0;

  virtual void getFiberLocations(int nFibers, double *xi) = 0;
  virtual void getFiberWeights(int nFibers, double *wt) = 0;

  virtual SectionIntegration *getCopy(void) = 0;

  virtual void getLocationsDeriv(int nFibers, double *dptsdh);
  virtual void getWeightsDeriv(int nFibers, double *dwtsdh);

  virtual void Print(OPS_Stream &s, int flag = 0) = 0;
};

#endif
