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
// $Date: 2007-10-26 04:49:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/UserDefinedBeamIntegration.cpp,v $

#include <UserDefinedBeamIntegration.h>

#include <ID.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <elementAPI.h>

void* OPS_UserDefinedBeamIntegration(int& integrationTag, ID& secTags)
{
    if(OPS_GetNumRemainingInputArgs() < 5) {
	opserr<<"insufficient arguments:integrationTag,N,secTags,locations,weights\n";
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
    Vector pt(N), wt(N);
    if(OPS_GetNumRemainingInputArgs() < 3*N) {
	opserr<<"There must be "<<N<<"secTags,locations and weights\n";
	return 0;
    }

    // secTags
    int *secptr = &secTags(0);
    if(OPS_GetIntInput(&N,secptr) < 0) return 0;

    // locations
    double *locptr = &pt(0);
    if(OPS_GetDoubleInput(&N,locptr) < 0) return 0;

    // weights
    double *wtptr = &wt(0);
    if(OPS_GetDoubleInput(&N,wtptr) < 0) return 0;
    
    return new UserDefinedBeamIntegration(N,pt,wt);
}

UserDefinedBeamIntegration::UserDefinedBeamIntegration(int nIP,
						       const Vector &pt,
						       const Vector &wt):
  BeamIntegration(BEAM_INTEGRATION_TAG_UserDefined),
  pts(nIP), wts(nIP)
{
  for (int i = 0; i < nIP; i++) {
    if (pt(i) < 0.0 || pt(i) > 1.0)
      opserr << "UserDefinedBeamIntegration::UserDefinedBeamIntegration -- point lies outside [0,1]" << endln;
    //if (wt(i) < 0.0 || wt(i) > 1.0)
    //opserr << "UserDefinedBeamIntegration::UserDefinedBeamIntegration -- weight lies outside [0,1]" << endln;
    pts(i) = pt(i);
    wts(i) = wt(i);
  }
}

UserDefinedBeamIntegration::UserDefinedBeamIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_UserDefined)
{
 
}

UserDefinedBeamIntegration::~UserDefinedBeamIntegration()
{
  // Nothing to do
}

void
UserDefinedBeamIntegration::getSectionLocations(int numSections,
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
UserDefinedBeamIntegration::getSectionWeights(int numSections,
					      double L, double *wt)
{
  int nIP = wts.Size();

  int i;
  for (i = 0; i < nIP; i++)
    wt[i] = wts(i);
  for ( ; i < numSections; i++)
    wt[i] = 1.0;
}

int
UserDefinedBeamIntegration::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 2)
    return -1;

  int point = atoi(argv[1]);
  if (point < 1)
    return -1;

  int Np = wts.Size();

  if (strcmp(argv[0],"pt") == 0 && point <= Np) {
    param.setValue(pts(point-1));
    return param.addObject(point, this);
  }
  else if (strcmp(argv[0],"wt") == 0 && point <= Np) {
    param.setValue(wts(point-1));
    return param.addObject(10+point, this);
  }
  else
    return -1;
}

int
UserDefinedBeamIntegration::updateParameter(int parameterID, Information &info)
{
  if (parameterID <= 10) { // pt
    pts(parameterID-1) = info.theDouble;
    return 0;
  }
  else if (parameterID <= 20) { // wt
    wts(parameterID-10-1) = info.theDouble;
    return 0;
  }
  else
    return -1;
}

BeamIntegration*
UserDefinedBeamIntegration::getCopy(void)
{
  int nIP = pts.Size();

  return new UserDefinedBeamIntegration(nIP, pts, wts);
}

int
UserDefinedBeamIntegration::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  int dbTag = this->getDbTag();

  int nIP = pts.Size();
  static ID iData(1);
  iData(0) = nIP;
  res = theChannel.sendID(dbTag, cTag, iData);
  if (res < 0) {
    opserr << "UserDefinedBeamIntegration::sendSelf - failed to send ID data" << endln;
    return res;
  }  

  Vector dData(nIP*2);
  for (int i=0; i<nIP; i++) {
    dData(i) = pts(i);
    dData(i+nIP) = wts(i);
  }
  res = theChannel.sendVector(dbTag, cTag, dData);
  if (res < 0) {
    opserr << "UserDefinedBeamIntegration::sendSelf - failed to send Vector data" << endln;
    return res;
  }
  
  return res;
}

int
UserDefinedBeamIntegration::recvSelf(int cTag, Channel &theChannel,
				     FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dbTag = this->getDbTag();
  
  static ID iData(1);
  res = theChannel.recvID(dbTag, cTag, iData);
  if (res < 0) {
    opserr << "UserDefinedBeamIntegration::recvSelf - failed to recv ID data" << endln;
    return res;
  }  
  int nIP = iData(0);
  
  pts.resize(nIP);
  wts.resize(nIP);

  Vector dData(nIP*2);
  res = theChannel.recvVector(dbTag, cTag, dData);
  if (res < 0) {
    opserr << "UserDefinedBeamIntegration::recvSelf - failed to recv Vector data" << endln;
    return res;
  }    

  for (int i=0; i<nIP; i++) {
    pts(i) = dData(i);
    wts(i) = dData(i+nIP);
  }
  
  return res;
}

void
UserDefinedBeamIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"UserDefined\", ";
		s << "\"points\": [";
		int nIP = pts.Size();
		for (int i = 0; i < nIP-1; i++)
			s << pts(i) << ", ";
		s << pts(nIP - 1) << "], ";
		s << "\"weights\": [";
		nIP = wts.Size();
		for (int i = 0; i < nIP-1; i++)
			s << wts(i) << ", ";
		s << wts(nIP - 1) << "]}";
	}

	else {
		s << "UserDefined" << endln;
		s << " Points: " << pts;
		s << " Weights: " << wts;
	}
}
