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
	double oneOverL = 1.0/L;
	int N = 5 ;
    
    // start with the two end nodes at the first two indices, then put them back
    Vector pt(N);
    
    pt[0] = 0.0; // node I
    pt[1] = 1.0; // node J
    pt[2] = 0.0 + lpI*oneOverL;
    pt[3] = 0.5*((0.0 + lpI*oneOverL) + (1.0 - lpJ*oneOverL)); // middle is average
    pt[4] = 1.0 - lpJ*oneOverL;
    
    int nc = 2;
    Vector wc(nc);
    wc[0] = lpI*oneOverL; // node I  ## curvature has a weight of Lp
    wc[1] = lpJ*oneOverL; // node J  ## curvature has a weight of Lp

    // fill in the rest using low-order integration:
	int nf = N-nc;
    Vector R(nf);
    for (int i = 0; i < nf; i++) {
    	double sum = 0.0;
    	for (int j = 0; j < nc; j++)
		sum += pow(pt(j),i)*wc(j);
	  	R(i) = 1.0/(i+1) - sum;
    }
    
    Matrix J(nf,nf);
    for (int i = 0; i < nf; i++)
      for (int j = 0; j < nf; j++)
	J(i,j) = pow(pt(nc+j),i);
    
    Vector wf(nf);
    
    J.Solve(R, wf);
    
     
    wt[0] = wc(0); // node I
    wt[1] = wf(0); // elastic
    wt[2] = wf(1); // elastic
    wt[3] = wf(2); // elastic
    wt[4] = wc(1); // node J
    
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
