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
// $Date: 2007-10-26 04:47:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/HingeMidpointBeamIntegration.cpp,v $

/*
 * Reference

Scott, M. H. and G. L. Fenves. "Plastic Hinge Integration Methods for
Force-Based Beam-Column Elements." Journal of Structural Engineering,
132(2):244-252, February 2006.

 *
 */

#include <HingeMidpointBeamIntegration.h>
#include <ElementalLoad.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <math.h>
#include <elementAPI.h>
#include <ID.h>

void* OPS_HingeMidpointBeamIntegration(int& integrationTag, ID& secTags)
{
    if(OPS_GetNumRemainingInputArgs() < 6) {
	opserr<<"insufficient arguments:integrationTag,secTagI,lpI,secTagJ,lpJ,secTagE\n";
	return 0;
    }

    // inputs: 
    int iData[4];
    double dData[2];
    int numData = 2;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	opserr << "WARNING: failed to get tag and secTagI\n";
	return 0;
    }
    numData = 1;
    if(OPS_GetDoubleInput(&numData,&dData[0]) < 0) {
	opserr << "WARNING: failed to get lpI\n";
	return 0;
    }
    if(OPS_GetIntInput(&numData,&iData[2]) < 0) {
	opserr << "WARNING: failed to get secTagJ\n";
	return 0;
    }
    if(OPS_GetDoubleInput(&numData,&dData[1]) < 0) {
	opserr << "WARNING: failed to get lpJ\n";
	return 0;
    }
    if(OPS_GetIntInput(&numData,&iData[3]) < 0) {
	opserr << "WARNING: failed to get secTagE\n";
	return 0;
    }

    integrationTag = iData[0];
    secTags.resize(4);
    secTags(0) = iData[1];
    secTags(1) = iData[3];
    secTags(2) = iData[3];
    secTags(3) = iData[2];

    return new HingeMidpointBeamIntegration(dData[0],dData[1]);
}

HingeMidpointBeamIntegration::HingeMidpointBeamIntegration(double lpi,
							   double lpj):
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeMidpoint), lpI(lpi), lpJ(lpj)
{
  // Nothing to do
}

HingeMidpointBeamIntegration::HingeMidpointBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeMidpoint), lpI(0.0), lpJ(0.0)
{

}

HingeMidpointBeamIntegration::~HingeMidpointBeamIntegration()
{
  // Nothing to do
}

void
HingeMidpointBeamIntegration::getSectionLocations(int numSections, double L,
						  double *xi)
{
  double halfOneOverL = 0.5/L;

  xi[0] = lpI*halfOneOverL;
  xi[3] = 1.0-lpJ*halfOneOverL;

  static const double oneRoot3 = 1.0/sqrt(3.0);

  double alpha = 0.5-(lpI+lpJ)*halfOneOverL;
  double beta  = 0.5+(lpI-lpJ)*halfOneOverL;
  xi[1] = alpha*(-oneRoot3) + beta;
  xi[2] = alpha*(oneRoot3) + beta;

  for (int i = 4; i < numSections; i++)
    xi[i] = 0.0;
}

void
HingeMidpointBeamIntegration::getSectionWeights(int numSections, double L,
						double *wt)
{
  double oneOverL = 1.0/L;

  wt[0] = lpI*oneOverL;
  wt[3] = lpJ*oneOverL;

  wt[1] = 0.5-0.5*(lpI+lpJ)*oneOverL;
  wt[2] = 0.5-0.5*(lpI+lpJ)*oneOverL;

  for (int i = 4; i < numSections; i++)
    wt[i] = 1.0;
}

BeamIntegration*
HingeMidpointBeamIntegration::getCopy(void)
{
  return new HingeMidpointBeamIntegration(lpI, lpJ);
}

int
HingeMidpointBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(2);

  data(0) = lpI;
  data(1) = lpJ;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "HingeMidpointBeamIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
HingeMidpointBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				       FEM_ObjectBroker &theBroker)
{
  static Vector data(2);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "HingeMidpointBeamIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  lpI = data(0);
  lpJ = data(1);

  return 0;
}

int
HingeMidpointBeamIntegration::setParameter(const char **argv, int argc,
					   Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"lpI") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"lpJ") == 0)
    return param.addObject(2, this);

  if (strcmp(argv[0],"lp") == 0)
    return param.addObject(3, this);

  return -1;
}

int
HingeMidpointBeamIntegration::updateParameter(int parameterID,
					      Information &info)
{
  switch (parameterID) {
  case 1:
    lpI = info.theDouble;
    return 0;
  case 2:
    lpJ = info.theDouble;
    return 0;
  case 3:
    lpI = lpJ = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
HingeMidpointBeamIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
HingeMidpointBeamIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"HingeMidpoint\", ";
		s << "\"lpI\": " << lpI << ", ";
		s << "\"lpJ\": " << lpJ << "}";
	}
	
	else {
		s << "HingeMidpoint" << endln;
		s << " lpI = " << lpI;
		s << " lpJ = " << lpJ << endln;
	}
}

void 
HingeMidpointBeamIntegration::getLocationsDeriv(int numSections, double L,
						double dLdh, double *dptsdh)
{
  double oneOverL = 1.0/L;
  double halfOneOverL = 0.5*oneOverL;

  for (int i = 0; i < numSections; i++)
    dptsdh[i] = 0.0;

  static const double oneRoot3 = 1.0/sqrt(3.0);

  if (parameterID == 1) { // lpI
    dptsdh[0] = halfOneOverL;
    dptsdh[1] = -0.5*(1.0-oneRoot3)*oneOverL + oneOverL;
    dptsdh[2] = -0.5*(1.0+oneRoot3)*oneOverL + oneOverL;
  }
  if (parameterID == 2) { // lpJ
    dptsdh[1] = -0.5*(1.0-oneRoot3)*oneOverL;
    dptsdh[2] = -0.5*(1.0+oneRoot3)*oneOverL;
    dptsdh[3] = -halfOneOverL;
  }
  if (parameterID == 3) { // lpI and lpJ
    dptsdh[0] = halfOneOverL;
    dptsdh[1] = -(1.0-oneRoot3)*oneOverL + oneOverL;
    dptsdh[2] = -(1.0+oneRoot3)*oneOverL + oneOverL;
    dptsdh[3] = -halfOneOverL;
  }

  return;

  if (dLdh != 0.0) {
    dptsdh[0] = -0.5*(lpI*dLdh)/(L*L);
    dptsdh[1] = dLdh + 0.5*(lpJ*dLdh)/(L*L);
    // STILL TO DO
    //opserr << "getPointsDeriv -- to do" << endln;
  }

  return;
}

void
HingeMidpointBeamIntegration::getWeightsDeriv(int numSections, double L,
					      double dLdh, double *dwtsdh)
{
  double oneOverL = 1.0/L;

  for (int i = 0; i < numSections; i++)
    dwtsdh[i] = 0.0;

  if (parameterID == 1) { // lpI
    dwtsdh[0] = oneOverL;
    dwtsdh[1] = -0.5*oneOverL;
    dwtsdh[2] = -0.5*oneOverL;
  }
  if (parameterID == 2) { // lpJ
    dwtsdh[1] = -0.5*oneOverL;
    dwtsdh[2] = -0.5*oneOverL;
    dwtsdh[3] = oneOverL;
  }
  if (parameterID == 3) { // lpI and lpJ
    dwtsdh[0] = oneOverL;
    dwtsdh[1] = -oneOverL;
    dwtsdh[2] = -oneOverL;
    dwtsdh[3] = oneOverL;
  }

  return;

  if (dLdh != 0.0) {
    dwtsdh[0] = -lpI*dLdh/(L*L);
    dwtsdh[1] = 0.5*(lpI+lpJ)*dLdh/(L*L);
    dwtsdh[2] = 0.5*(lpI+lpJ)*dLdh/(L*L);
    dwtsdh[3] = -lpJ*dLdh/(L*L);
  }

  return;
}
