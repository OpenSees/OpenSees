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

// $Revision: 1.5 $
// $Date: 2008-12-03 23:43:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/HingeRadauTwoBeamIntegration.cpp,v $

/*
 * Reference

Scott, M. H. and G. L. Fenves. "Plastic Hinge Integration Methods for
Force-Based Beam-Column Elements." Journal of Structural Engineering,
132(2):244-252, February 2006.

 *
 */

#include <HingeRadauTwoBeamIntegration.h>
#include <ElementalLoad.h>

#include <math.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <ID.h>

void* OPS_HingeRadauTwoBeamIntegration(int& integrationTag, ID& secTags)
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
    secTags.resize(6);
    secTags(0) = iData[1];
    secTags(1) = iData[1];
    secTags(2) = iData[3];
    secTags(3) = iData[3];
    secTags(4) = iData[2];
    secTags(5) = iData[2];

    return new HingeRadauTwoBeamIntegration(dData[0],dData[1]);
}

HingeRadauTwoBeamIntegration::HingeRadauTwoBeamIntegration(double lpi,
							   double lpj):
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadauTwo), lpI(lpi), lpJ(lpj)
{
  // Nothing to do
}

HingeRadauTwoBeamIntegration::HingeRadauTwoBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_HingeRadauTwo), lpI(0.0), lpJ(0.0)
{

}

HingeRadauTwoBeamIntegration::~HingeRadauTwoBeamIntegration()
{
  // Nothing to do
}

void
HingeRadauTwoBeamIntegration::getSectionLocations(int numSections, double L,
						  double *xi)
{
  double oneOverL = 1.0/L;

  xi[0] = 0.0;
  xi[1] = 2.0/3*lpI*oneOverL;
  xi[4] = 1.0-2.0/3*lpJ*oneOverL;
  xi[5] = 1.0;

  static const double oneRoot3 = 1.0/sqrt(3.0);

  double alpha = 0.5 - 0.5*(lpI+lpJ)*oneOverL;
  double beta  = 0.5 + 0.5*(lpI-lpJ)*oneOverL;
  xi[2] = alpha*(-oneRoot3) + beta;
  xi[3] = alpha*(oneRoot3) + beta;

  for (int i = 6; i < numSections; i++)
    xi[i] = 0.0;
}

void
HingeRadauTwoBeamIntegration::getSectionWeights(int numSections, double L,
						double *wt)
{
  double oneOverL = 1.0/L;

  wt[0] = 0.25*lpI*oneOverL;
  wt[1] = 0.75*lpI*oneOverL;
  wt[4] = 0.75*lpJ*oneOverL;
  wt[5] = 0.25*lpJ*oneOverL;

  wt[2] = 0.5-0.5*(lpI+lpJ)*oneOverL;
  wt[3] = 0.5-0.5*(lpI+lpJ)*oneOverL;

  for (int i = 6; i < numSections; i++)
    wt[i] = 1.0;
}

BeamIntegration*
HingeRadauTwoBeamIntegration::getCopy(void)
{
  return new HingeRadauTwoBeamIntegration(lpI, lpJ);
}

int
HingeRadauTwoBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(2);

  data(0) = lpI;
  data(1) = lpJ;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "HingeRadauTwoBeamIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}

int
HingeRadauTwoBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				       FEM_ObjectBroker &theBroker)
{
  static Vector data(2);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "HingeRadauTwoBeamIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  lpI = data(0);
  lpJ = data(1);

  return 0;
}

int
HingeRadauTwoBeamIntegration::setParameter(const char **argv, int argc,
					   Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"lpI") == 0) {
    param.setValue(lpI);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"lpJ") == 0) {
    param.setValue(lpJ);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"lp") == 0) {
    param.setValue(lpI);
    return param.addObject(3, this);
  }
  return -1;
}

int
HingeRadauTwoBeamIntegration::updateParameter(int parameterID,
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
HingeRadauTwoBeamIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
HingeRadauTwoBeamIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"HingeRadauTwo\", ";
		s << "\"lpI\": " << lpI << ", ";
		s << "\"lpJ\": " << lpJ << "}";
	}
	
	else {
		s << "HingeRadauTwo" << endln;
		s << " lpI = " << lpI;
		s << " lpJ = " << lpJ << endln;
	}
}

void 
HingeRadauTwoBeamIntegration::getLocationsDeriv(int numSections,
						double L, double dLdh,
						double *dptsdh)
{
  double oneOverL = 1.0/L;

  for (int i = 0; i < numSections; i++)
    dptsdh[i] = 0.0;

  //return;

  static const double oneRoot3 = 1.0/sqrt(3.0);

  if (parameterID == 1) { // lpI
    dptsdh[1] = 2.0/3*oneOverL;
    dptsdh[2] = 0.5*oneOverL*(1.0+oneRoot3);
    dptsdh[3] = 0.5*oneOverL*(1.0-oneRoot3);
  }

  if (parameterID == 2) { // lpJ
    dptsdh[2] = -0.5*oneOverL*(1.0-oneRoot3);
    dptsdh[3] = -0.5*oneOverL*(1.0+oneRoot3);
    dptsdh[4] = -2.0/3*oneOverL;
  }

  if (parameterID == 3) { // lpI and lpJ
    dptsdh[1] = 2.0/3*oneOverL;
    dptsdh[2] =  oneOverL*oneRoot3;
    dptsdh[3] = -oneOverL*oneRoot3;
    dptsdh[4] = -2.0/3*oneOverL;
  }

  return;

  if (dLdh != 0.0) {
    dptsdh[1] = -2.0/3*lpI*dLdh/(L*L);
    double dalphadh =  0.5*(lpI+lpJ)*dLdh/(L*L);
    double dbetadh  = -0.5*(lpI-lpJ)*dLdh/(L*L);
    dptsdh[2] = -oneRoot3*dalphadh + dbetadh;
    dptsdh[3] =  oneRoot3*dalphadh + dbetadh;
    double alpha = 0.5*L - 0.5*(lpI+lpJ);
    double beta  = 0.5*L + 0.5*(lpI-lpJ);
    dptsdh[2] = -(-oneRoot3*alpha+beta)*dLdh/(L*L) + 0.5*dLdh;///(L*L);
    dptsdh[3] = -( oneRoot3*alpha+beta)*dLdh/(L*L) + 0.5*dLdh;///(L*L);
    dptsdh[4] = -(L-2.0/3*lpJ)*dLdh/(L*L);
    dptsdh[5] =  -L*dLdh/(L*L);
    // STILL TO DO
    //opserr << "getPointsDeriv -- to do" << endln;
  }

  return;
}

void
HingeRadauTwoBeamIntegration::getWeightsDeriv(int numSections,
					      double L, double dLdh,
					      double *dwtsdh)
{
  double oneOverL = 1.0/L;

  for (int i = 0; i < numSections; i++)
    dwtsdh[i] = 0.0;

  if (parameterID == 1) { // lpI
    dwtsdh[0] = 0.25*oneOverL;
    dwtsdh[1] = 0.75*oneOverL;
    dwtsdh[2] = -0.5*oneOverL;
    dwtsdh[3] = -0.5*oneOverL;
  }

  if (parameterID == 2) { // lpJ
    dwtsdh[2] = -0.5*oneOverL;
    dwtsdh[3] = -0.5*oneOverL;
    dwtsdh[4] = 0.75*oneOverL;
    dwtsdh[5] = 0.25*oneOverL;
  }

  if (parameterID == 3) { // lpI and lpJ
    dwtsdh[0] = 0.25*oneOverL;
    dwtsdh[1] = 0.75*oneOverL;
    dwtsdh[2] = -oneOverL;
    dwtsdh[3] = -oneOverL;
    dwtsdh[4] = 0.75*oneOverL;
    dwtsdh[5] = 0.25*oneOverL;
  }

  return;

  if (dLdh != 0.0) {
    dwtsdh[0] = -0.25*lpI*dLdh/(L*L);
    dwtsdh[1] = -0.75*lpI*dLdh/(L*L);
    dwtsdh[2] = 0.5*(lpI+lpJ)*dLdh/(L*L);
    dwtsdh[3] = 0.5*(lpI+lpJ)*dLdh/(L*L);
    dwtsdh[4] = -0.75*lpJ*dLdh/(L*L);
    dwtsdh[5] = -0.25*lpJ*dLdh/(L*L);
    // STILL TO DO
    //opserr << "getWeightsDeriv -- to do" << endln;
  }

  return;
}
