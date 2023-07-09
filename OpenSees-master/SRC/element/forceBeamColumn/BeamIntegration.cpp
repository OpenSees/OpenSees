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

// $Revision: 1.4 $
// $Date: 2006-09-05 22:57:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/BeamIntegration.cpp,v $

#include <BeamIntegration.h>
#include <Matrix.h>

#include <MapOfTaggedObjects.h>

static MapOfTaggedObjects theBeamIntegrationRuleObjects;

bool OPS_addBeamIntegrationRule(BeamIntegrationRule *newComponent) {
  return theBeamIntegrationRuleObjects.addComponent(newComponent);
}

BeamIntegrationRule *OPS_getBeamIntegrationRule(int tag) {

  TaggedObject *theResult = theBeamIntegrationRuleObjects.getComponentPtr(tag);
  if (theResult == 0) {
    opserr << "BeamIntegrationRule - none found with tag: " << tag << endln;
    return 0;
  }
  BeamIntegrationRule *theMat = (BeamIntegrationRule *)theResult;

  return theMat;
}

void OPS_clearAllBeamIntegrationRule(void) {
  theBeamIntegrationRuleObjects.clearAll();
}

BeamIntegration::BeamIntegration(int classTag):
  MovableObject(classTag)
{
  // Nothing to do
}

BeamIntegration::~BeamIntegration()
{
  // Nothing to do
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
