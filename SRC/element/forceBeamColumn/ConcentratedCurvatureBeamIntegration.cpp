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
**                                                                    **
** ****************************************************************** *
** This Beam Integration Module was Developed by:                     **
**   Silvia Mazzoni (silviamazzoni@yahoo.com)                         **
**   Michael H. Scott (michael@portwooddigital.com)                   **
**   pushed to OpenSees November 2022                                 **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2022-11-20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ConcentratedCurvatureBeamIntegration.cpp,v $

#include <ConcentratedCurvatureBeamIntegration.h>

#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <math.h>
#include <elementAPI.h>

void* OPS_ConcentratedCurvatureBeamIntegration(int& integrationTag, ID& secTags)
{
    if(OPS_GetNumRemainingInputArgs() < 6) {
	opserr<<"insufficient arguments:integrationTag,secTagI,LpI,secTagJ,LpJ,secTagE\n";
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

    
    int N = 5;
    integrationTag = iData[0];
    secTags.resize(N);
    secTags(0) = iData[1];
    secTags(1) = iData[3];
    secTags(2) = iData[3];
    secTags(3) = iData[3];
    secTags(4) = iData[2];


    
    return new ConcentratedCurvatureBeamIntegration(dData[0],dData[1]);
}

ConcentratedCurvatureBeamIntegration::ConcentratedCurvatureBeamIntegration(double lpi,
						     double lpj):
  BeamIntegration(BEAM_INTEGRATION_TAG_ConcentratedCurvature), lpI(lpi), lpJ(lpj)
{
  // Nothing to do
}

ConcentratedCurvatureBeamIntegration::ConcentratedCurvatureBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_ConcentratedCurvature), lpI(0.0), lpJ(0.0)
{

}

ConcentratedCurvatureBeamIntegration::~ConcentratedCurvatureBeamIntegration()
{
  // Nothing to do
}
void
ConcentratedCurvatureBeamIntegration::getSectionLocations(int numSections, double L,
					       double *xi)
{
	double oneOverL = 1.0/L;
	
    // section locations
    xi[0] = 0.0; // node I
    xi[1] = 0.0 + lpI*oneOverL;
    xi[2] = 0.5*((0.0 + lpI*oneOverL) + (1.0 - lpJ*oneOverL)); // middle is average
    xi[3] = 1.0 - lpJ*oneOverL;
    xi[4] = 1.0; // node J
}

void
ConcentratedCurvatureBeamIntegration::getSectionWeights(int numSections, double L,
					     double *wt)
{
    const double a = lpI / L;
    const double b = lpJ / L;
    const double t = 1.0 - a - b;

    // Residual integration moments after assigning
    // the concentrated-curvature end contributions.
    const double R0 = t;
    const double R1 = 0.5 - b;
    const double R2 = 1.0 / 3.0 - b;

    // Shift moments to the first elastic location.
    const double S1 = R1 - a * R0;
    const double S2 = R2 - 2.0 * a * R1 + a * a * R0;

    // End curvature sections.
    wt[0] = a;
    wt[4] = b;

    // Elastic sections: exact low-order moment fitting.
    wt[3] = 2.0 * S2 / (t * t) - S1 / t;
    wt[2] = 4.0 * (t * S1 - S2) / (t * t);
    wt[1] = R0 - wt[2] - wt[3];
}

BeamIntegration*
ConcentratedCurvatureBeamIntegration::getCopy(void)
{
  return new ConcentratedCurvatureBeamIntegration(lpI, lpJ);
}

int
ConcentratedCurvatureBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(2);

  data(0) = lpI;
  data(1) = lpJ;

  int dbTag = this->getDbTag();

  if (theChannel.sendVector(dbTag, cTag, data) < 0) {
    opserr << "ConcentratedCurvatureBeamIntegration::sendSelf() - failed to send Vector data\n";
    return -1;
  }    

  return 0;
}


int
ConcentratedCurvatureBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				  FEM_ObjectBroker &theBroker)
{
  static Vector data(2);

  int dbTag = this->getDbTag();

  if (theChannel.recvVector(dbTag, cTag, data) < 0)  {
    opserr << "ConcentratedCurvatureBeamIntegration::recvSelf() - failed to receive Vector data\n";
    return -1;
  }
  
  lpI = data(0);
  lpJ = data(1);

  return 0;
}

int
ConcentratedCurvatureBeamIntegration::setParameter(const char **argv, int argc,
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
ConcentratedCurvatureBeamIntegration::updateParameter(int parameterID,
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
ConcentratedCurvatureBeamIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

void
ConcentratedCurvatureBeamIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"ConcentratedCurvature\", ";
		s << "\"lpI\": " << lpI << ", ";
		s << "\"lpJ\": " << lpJ << "}";
	}
	
	else {
		s << "ConcentratedCurvature" << endln;
		s << " lpI = " << lpI;
		s << " lpJ = " << lpJ << endln;
	}
}

void 
ConcentratedCurvatureBeamIntegration::getLocationsDeriv(int numSections, double L,
					   double dLdh, double *dptsdh)
{
  double oneOverL = 1.0/L;

  for (int i = 0; i < numSections; i++)
    dptsdh[i] = 0.0;

  //return;

  static const double oneRoot3 = 1.0/sqrt(3.0);

  if (parameterID == 1) { // lpI
    dptsdh[1] = 1.0*oneOverL;
    dptsdh[2] = 0.5*1.0*oneOverL;
  }

  if (parameterID == 2) { // lpJ
    dptsdh[2] = -0.5*1.0*oneOverL;
    dptsdh[3] = -1.0*oneOverL;
  }

  if (parameterID == 3) { // lpI and lpJ
    dptsdh[1] = 1.0*oneOverL;
    dptsdh[2] = 0.5*1.0*oneOverL -0.5*1.0*oneOverL; // i don't think this is right...
    dptsdh[3] = -1.0*oneOverL;
  }

  return;

  if (dLdh != 0.0) {
    // STILL TO DO
    //opserr << "getPointsDeriv -- to do" << endln;
  }

  return;
}

void
ConcentratedCurvatureBeamIntegration::getWeightsDeriv(int numSections, double L,
					 double dLdh, double *dwtsdh)
{

	// this gets complicated....
    dwtsdh[0] = 0.;
    dwtsdh[1] = 0.;
    dwtsdh[2] = 0.;
    dwtsdh[3] = 0.;
    dwtsdh[4] = 0.;

    
    
    
  return;
}
