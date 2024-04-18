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
// $Date: 2007-10-12 21:03:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/FixedLocationBeamIntegration.cpp,v $

#include <FixedLocationBeamIntegration.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <math.h>
#include <elementAPI.h>
#include <ID.h>

void* OPS_FixedLocationBeamIntegration(int& integrationTag, ID& secTags)
{
    if(OPS_GetNumRemainingInputArgs() < 4) {
	opserr<<"insufficient arguments:integrationTag,N,secTags,locations\n";
	return 0;
    }

    // inputs: integrationTag,N
    int iData[2];
    int numData = 2;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) return 0;

    integrationTag = iData[0];
    int N = iData[1];
    if(N > 0) {
	secTags.resize(N);
    } else {
	secTags.resize(1);
	N = 1;
    }

    // check argumments
    Vector pt(N);
    if(OPS_GetNumRemainingInputArgs() < 2*N) {
	opserr<<"There must be "<<N<<"secTags and locations\n";
	return 0;
    }

    // secTags
    int *secptr = &secTags(0);
    if(OPS_GetIntInput(&N,secptr) < 0) return 0;

    // locations
    double *locptr = &pt(0);
    if(OPS_GetDoubleInput(&N,locptr) < 0) return 0;

    return new FixedLocationBeamIntegration(N,pt);
}

FixedLocationBeamIntegration::FixedLocationBeamIntegration(int nIP,
							   const Vector &pt):
  BeamIntegration(BEAM_INTEGRATION_TAG_FixedLocation),
  pts(nIP), wts(nIP)
{
  for (int i = 0; i < nIP; i++) {
    if (pt(i) < 0.0 || pt(i) > 1.0)
      opserr << "FixedLocationBeamIntegration::FixedLocationBeamIntegration -- point lies outside [0,1]" << endln;

    pts(i) = pt(i);
  }

  this->calculateWeights();
}

void
FixedLocationBeamIntegration::calculateWeights()
{
  int nIP = pts.Size();
  
  Vector R(nIP);
  for (int i = 0; i < nIP; i++)
    R(i) = 1.0/(i+1);

  Matrix J(nIP,nIP);
  for (int i = 0; i < nIP; i++)
    for (int j = 0; j < nIP; j++)
      J(i,j) = pow(pts(j),i);

  J.Solve(R, wts);
}

FixedLocationBeamIntegration::FixedLocationBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_FixedLocation)
{
 
}

FixedLocationBeamIntegration::~FixedLocationBeamIntegration()
{
  // Nothing to do
}

void
FixedLocationBeamIntegration::getSectionLocations(int numSections,
						  double L, double *xi)
{
  int nIP = pts.Size();

  int i;
  for (i = 0; i < nIP; i++)
    xi[i] = pts(i);
  for ( ; i < numSections; i++)
    xi[i] = 0.0;
}

void
FixedLocationBeamIntegration::getSectionWeights(int numSections,
						double L, double *wt)
{
  int nIP = wts.Size();

  int i;
  for (i = 0; i < nIP; i++)
    wt[i] = wts(i);
  for ( ; i < numSections; i++)
    wt[i] = 1.0;
}

BeamIntegration*
FixedLocationBeamIntegration::getCopy(void)
{
  int nIP = pts.Size();

  return new FixedLocationBeamIntegration(nIP, pts);
}

int
FixedLocationBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  int dbTag = this->getDbTag();

  int nIP = wts.Size();
  static ID iData(1);
  iData(0) = nIP;
  res = theChannel.sendID(dbTag, cTag, iData);
  if (res < 0) {
    opserr << "FixedLocationBeamIntegration::sendSelf - failed to send ID data" << endln;
    return res;
  }  

  Vector dData(nIP);
  for (int i=0; i<nIP; i++)
    dData(i) = pts(i);

  res = theChannel.sendVector(dbTag, cTag, dData);
  if (res < 0) {
    opserr << "FixedLocationBeamIntegration::sendSelf - failed to send Vector data" << endln;
    return res;
  }
  
  return res;  
}

int
FixedLocationBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dbTag = this->getDbTag();
  
  static ID iData(1);
  res = theChannel.recvID(dbTag, cTag, iData);
  if (res < 0) {
    opserr << "FixedLocationBeamIntegration::recvSelf - failed to recv ID data" << endln;
    return res;
  }  
  int nIP = iData(0);
  
  pts.resize(nIP);
  wts.resize(nIP);

  Vector dData(nIP);
  res = theChannel.recvVector(dbTag, cTag, dData);
  if (res < 0) {
    opserr << "FixedLocationBeamIntegration::recvSelf - failed to recv Vector data" << endln;
    return res;
  }    

  for (int i=0; i<nIP; i++)
    pts(i) = dData(i);

  this->calculateWeights();
  
  return res;
}

int
FixedLocationBeamIntegration::setParameter(const char **argv,
					   int argc, Information &info)
{
  return -1;
}

int
FixedLocationBeamIntegration::updateParameter(int parameterID,
					      Information &info)
{
  // Does nothing for now -- MHS
  return 0;
}

int
FixedLocationBeamIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  // For Terje to do
  return 0;
}

void
FixedLocationBeamIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"FixedLocation\", ";
		s << "\"points\": [";
		int nIP = pts.Size();
		for (int i = 0; i < nIP-1; i++)
			s << pts(i) << ", ";
		s << pts(nIP - 1) << "], ";
		s << "\"weights\": [";
		double sum = 0.0;
		nIP = wts.Size();
		for (int i = 0; i < nIP-1; i++) {
			s << wts(i) << ", ";
			sum += fabs(wts(i));
		}
		s << wts(nIP - 1) << "], ";
		s << "\"conditionNumber\": " << sum << "}";
	}

	else {
		s << "FixedLocation" << endln;
		s << " Points: " << pts;
		s << " Weights: " << wts;
		double sum = 0.0;
		int N = wts.Size();
		for (int i = 0; i < N; i++)
			sum += fabs(wts(i));
		s << " Condition Number: " << sum << endln;
	}
}

void 
FixedLocationBeamIntegration::getLocationsDeriv(int numSections, double L,
						double dLdh, double *dptsdh)
{
  for (int i = 0; i < numSections; i++)
    dptsdh[i] = 0.0;

  return;
}

void
FixedLocationBeamIntegration::getWeightsDeriv(int numSections, double L,
					      double dLdh, double *dwtsdh)
{
  for (int i = 0; i < numSections; i++)
    dwtsdh[i] = 0.0;

  return;
}
