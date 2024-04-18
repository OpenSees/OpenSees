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
// $Date: 2007-10-13 00:45:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/MidDistanceBeamIntegration.cpp,v $

#include <MidDistanceBeamIntegration.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <elementAPI.h>
#include <ID.h>

void* OPS_MidDistanceBeamIntegration(int& integrationTag, ID& secTags)
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

    return new MidDistanceBeamIntegration(N,pt);
}

MidDistanceBeamIntegration::MidDistanceBeamIntegration(int nIP,
						       const Vector &pt):
  BeamIntegration(BEAM_INTEGRATION_TAG_MidDistance),
  pts(nIP), wts(nIP)
{
  for (int i = 0; i < nIP; i++) {
    if (pt(i) < 0.0 || pt(i) > 1.0)
      opserr << "MidDistanceBeamIntegration::MidDistanceBeamIntegration -- point lies outside [0,1]" << endln;

    pts(i) = pt(i);
  }

  for (int i = 0; i < nIP; i++) {
    int key = i;
    for (int j = i+1; j < nIP; j++) {
      if (pts(j) < pts(key)) {
	key = j;
	opserr << "MidDistanceBeamIntegration::MidDistanceBeamIntegration -- point are not sorted; sort before calling constructor" << endln;
      }
    }
    //double temp = pts(i);
    //pts(i) = pts(key);
    //pts(key) = temp;
  }

  Vector mids(nIP-1);

  for (int i = 0; i < nIP-1; i++) {
    mids(i) = 0.5*(pts(i)+pts(i+1));
  }

  wts(0) = mids(0);
  wts(nIP-1) = 1.0-mids(nIP-2);
  for (int i = 1; i < nIP-1; i++) {
    wts(i) = mids(i)-mids(i-1);
  }
}

MidDistanceBeamIntegration::MidDistanceBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_MidDistance)
{

}

MidDistanceBeamIntegration::~MidDistanceBeamIntegration()
{
  // Nothing to do
}

void
MidDistanceBeamIntegration::getSectionLocations(int numSections,
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
MidDistanceBeamIntegration::getSectionWeights(int numSections,
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
MidDistanceBeamIntegration::getCopy(void)
{
  int nIP = pts.Size();

  return new MidDistanceBeamIntegration(nIP, pts);
}

int
MidDistanceBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  int dbTag = this->getDbTag();

  int nIP = pts.Size();
  static ID iData(1);
  iData(0) = nIP;
  res = theChannel.sendID(dbTag, cTag, iData);
  if (res < 0) {
    opserr << "MidDistanceBeamIntegration::sendSelf - failed to send ID data" << endln;
    return res;
  }  

  Vector dData(nIP*2);
  for (int i=0; i<nIP; i++) {
    dData(i) = pts(i);
    dData(i+nIP) = wts(i);
  }
  res = theChannel.sendVector(dbTag, cTag, dData);
  if (res < 0) {
    opserr << "MidDistanceBeamIntegration::sendSelf - failed to send Vector data" << endln;
    return res;
  }
  
  return res;
}

int
MidDistanceBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				     FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dbTag = this->getDbTag();
  
  static ID iData(1);
  res = theChannel.recvID(dbTag, cTag, iData);
  if (res < 0) {
    opserr << "MidDistanceBeamIntegration::recvSelf - failed to recv ID data" << endln;
    return res;
  }  
  int nIP = iData(0);
  
  pts.resize(nIP);
  wts.resize(nIP);

  Vector dData(nIP*2);
  res = theChannel.recvVector(dbTag, cTag, dData);
  if (res < 0) {
    opserr << "MidDistanceBeamIntegration::recvSelf - failed to recv Vector data" << endln;
    return res;
  }    

  for (int i=0; i<nIP; i++) {
    pts(i) = dData(i);
    wts(i) = dData(i+nIP);
  }
  
  return res;
}

int
MidDistanceBeamIntegration::setParameter(const char **argv,
					 int argc, Information &info)
{
  return -1;
}

int
MidDistanceBeamIntegration::updateParameter(int parameterID,
					    Information &info)
{
  // Does nothing for now -- MHS
  return 0;
}

int
MidDistanceBeamIntegration::activateParameter(int paramID)
{
  parameterID = paramID;

  // For Terje to do
  return 0;
}

void
MidDistanceBeamIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"MidDistance\", ";
		s << "\"points\": [";
		int nIP = pts.Size();
		for (int i = 0; i < nIP-1; i++)
			s << pts(i) << ", ";
		s << pts(nIP-1) << "], ";
		s << "\"weights\": [";
		nIP = wts.Size();
		for (int i = 0; i < nIP-1; i++)
			s << wts(i) << ", ";
		s << wts(nIP-1) << "]}";
	}
	
	else {
		s << "MidDistance" << endln;
		s << " Points: " << pts;
		s << " Weights: " << wts;
	}
}

void 
MidDistanceBeamIntegration::getLocationsDeriv(int numSections,
					      double L, double *dptsdh)
{
  for (int i = 0; i < numSections; i++)
    dptsdh[i] = 0.0;

  return;
}

void
MidDistanceBeamIntegration::getWeightsDeriv(int numSections,
					    double L, double *dwtsdh)
{
  for (int i = 0; i < numSections; i++)
    dwtsdh[i] = 0.0;

  return;
}
