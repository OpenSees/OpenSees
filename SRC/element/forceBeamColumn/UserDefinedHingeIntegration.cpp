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
// $Date: 2006-01-18 21:58:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/UserDefinedHingeIntegration.cpp,v $

#include <UserDefinedHingeIntegration.h>

#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <math.h>
#include <elementAPI.h>
#include <ID.h>

void* OPS_UserHingeBeamIntegration(int& integrationTag, ID& secTags)
{
    if(OPS_GetNumRemainingInputArgs() < 10) {
	opserr<<"insufficient arguments:integrationTag,secTagE,npL,secTagLs,ptLs,wtLs,npR,secTagRs,ptRs,wtRs\n";
	return 0;
    }

    // inputs: 
    int numData = 1;
    int secTagE;
    if(OPS_GetIntInput(&numData,&integrationTag) < 0) return 0;
    if(OPS_GetIntInput(&numData,&secTagE) < 0) return 0;

    // npL
    int npL;
    if(OPS_GetIntInput(&numData,&npL) < 0) return 0;
    if(npL <= 0) npL = 1;

    ID secTagL(npL);
    Vector ptL(npL), wtL(npL);
    if(OPS_GetNumRemainingInputArgs() < 3*npL) {
	opserr<<"There must be "<<npL<<"secTagL,ptL and wtL\n";
	return 0;
    }

    int *secptr = &secTagL(0);
    if(OPS_GetIntInput(&npL,secptr) < 0) return 0;

    double *locptr = &ptL(0);
    if(OPS_GetDoubleInput(&npL,locptr) < 0) return 0;

    double *wtptr = &wtL(0);
    if(OPS_GetDoubleInput(&npL,wtptr) < 0) return 0;

    // npR
    int npR;
    if(OPS_GetIntInput(&numData,&npR) < 0) return 0;
    if(npR <= 0) npR = 1;
    
    ID secTagR(npR);
    Vector ptR(npR), wtR(npR);
    if(OPS_GetNumRemainingInputArgs() < 3*npR) {
	opserr<<"There must be "<<npR<<"secTagR,ptR and wtR\n";
	return 0;
    }

    secptr = &secTagR(0);
    if(OPS_GetIntInput(&npR,secptr) < 0) return 0;

    locptr = &ptR(0);
    if(OPS_GetDoubleInput(&npR,locptr) < 0) return 0;

    wtptr = &wtR(0);
    if(OPS_GetDoubleInput(&npR,wtptr) < 0) return 0;

    // secTags
    secTags.resize(npL+npR+2);
    for(int i=0; i<npL; i++) {
	secTags(i) = secTagL(i);
    }
    for(int i=0; i<npR; i++) {
	secTags(i+npL) = secTagR(i);
    }
    secTags(npL+npR) = secTagE;
    secTags(npL+npR+1) = secTagE;
    
    return new UserDefinedHingeIntegration(npL,ptL,wtL,npR,ptR,wtR);
}

UserDefinedHingeIntegration::UserDefinedHingeIntegration(int npL,
							 const Vector &ptL,
							 const Vector &wtL,
							 int npR,
							 const Vector &ptR,
							 const Vector &wtR):
  BeamIntegration(BEAM_INTEGRATION_TAG_UserHinge),
  ptsL(npL), wtsL(npL), ptsR(npR), wtsR(npR)
{
  for (int i = 0; i < npL; i++) {
    if (ptL(i) < 0.0 || ptL(i) > 1.0)
      opserr << "UserDefinedHingeIntegration::UserDefinedHingeIntegration -- point lies outside [0,1]" << endln;
    if (wtL(i) < 0.0 || wtL(i) > 1.0)
      opserr << "UserDefinedHingeIntegration::UserDefinedHingeIntegration -- weight lies outside [0,1]" << endln;
    ptsL(i) = ptL(i);
    wtsL(i) = wtL(i);
  }

  for (int i = 0; i < npR; i++) {
    if (ptR(i) < 0.0 || ptR(i) > 1.0)
      opserr << "UserDefinedHingeIntegration::UserDefinedHingeIntegration -- point lies outside [0,1]" << endln;
    if (wtR(i) < 0.0 || wtR(i) > 1.0)
      opserr << "UserDefinedHingeIntegration::UserDefinedHingeIntegration -- weight lies outside [0,1]" << endln;
    ptsR(i) = ptR(i);
    wtsR(i) = wtR(i);
  }
}

UserDefinedHingeIntegration::UserDefinedHingeIntegration():
  BeamIntegration(BEAM_INTEGRATION_TAG_UserHinge)
{

}

UserDefinedHingeIntegration::~UserDefinedHingeIntegration()
{
  // Nothing to do
}

void
UserDefinedHingeIntegration::getSectionLocations(int numSections,
						 double L, double *xi)
{
  int npL = ptsL.Size();
  int npR = ptsR.Size();

  double lpI = 0.0;
  double lpJ = 0.0;
  int i, j;
  for (i = 0; i < npL; i++) {
    xi[i] = ptsL(i);
    lpI += wtsL(i);
  }
  for (j = 0; j < npR; j++, i++) {
    xi[i] = ptsR(j);
    lpJ += wtsR(j);
  }

  double alpha = 0.5-0.5*(lpI+lpJ);
  double beta  = 0.5+0.5*(lpI-lpJ);
  xi[i++] = alpha*(-1/sqrt(3.0)) + beta;
  xi[i++] = alpha*(1/sqrt(3.0)) + beta;

  for ( ; i < numSections; i++)
    xi[i] = 0.0;
}

void
UserDefinedHingeIntegration::getSectionWeights(int numSections,
					       double L, double *wt)
{
  int npL = wtsL.Size();
  int npR = wtsR.Size();

  double lpI = 0.0;
  double lpJ = 0.0;
  int i, j;
  for (i = 0; i < npL; i++) {
    wt[i] = wtsL(i);
    lpI += wtsL(i);
  }
  for (j = 0; j < npR; j++, i++) {
    wt[i] = wtsR(j);
    lpJ += wtsR(j);
  }

  double oneOverL = 1.0/L;
  wt[i++] = 0.5-0.5*(lpI+lpJ);
  wt[i++] = 0.5-0.5*(lpI+lpJ);

  for ( ; i < numSections; i++)
    wt[i] = 1.0;
}

BeamIntegration*
UserDefinedHingeIntegration::getCopy(void)
{
  int npL = ptsL.Size();
  int npR = ptsR.Size();

  return new UserDefinedHingeIntegration(npL, ptsL, wtsL,
					 npR, ptsR, wtsR);
}

int
UserDefinedHingeIntegration::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  int dbTag = this->getDbTag();

  int npL = ptsL.Size();
  int npR = ptsR.Size();
  
  static ID idata(2);
  idata(0) = npL;
  idata(1) = npR;
  res = theChannel.sendID(dbTag, cTag, idata);
  if (res < 0) {
    opserr << "UserDefinedHingeIntegration::sendSelf - failed to send ID data" << endln;
    return res;
  }

  Vector data(2*(npL+npR));
  for (int i = 0; i < npL; i++) {
    data(i)     = ptsL(i);
    data(i+npL) = wtsL(i);
  }
  for (int i = 0; i < npR; i++) {
    data(2*npL + i)     = ptsR(i);
    data(2*npL+npR + i) = wtsR(i);    
  }
  res = theChannel.sendVector(dbTag, cTag, data);
  if (res < 0) {
    opserr << "UserDefinedHingeIntegration::sendSelf - failed to send Vector data" << endln;
    return res;
  }
  
  return res;
}

int
UserDefinedHingeIntegration::recvSelf(int cTag, Channel &theChannel,
				      FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dbTag = this->getDbTag();

  static ID idata(2);
  res = theChannel.recvID(dbTag, cTag, idata);
  if (res < 0) {
    opserr << "UserDefinedHingeIntegration::recvSelf - failed to recv ID data" << endln;
    return res;
  }
  int npL = idata(0);
  int npR = idata(1);

  ptsL.resize(npL);
  wtsL.resize(npL);
  ptsR.resize(npR);
  wtsR.resize(npR);

  Vector data(2*(npL+npR));
  res = theChannel.recvVector(dbTag, cTag, data);
  if (res < 0) {
    opserr << "UserDefinedHingeIntegration::recvSelf - failed to recv Vector data" << endln;
    return res;
  }

  for (int i = 0; i < npL; i++) {
    ptsL(i) = data(i);
    wtsL(i) = data(i+npL);
  }
  for (int i = 0; i < npR; i++) {
    ptsR(i) = data(2*npL + i);
    wtsR(i) = data(2*npL+npR + i);
  }

  return res;
}

void
UserDefinedHingeIntegration::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "{\"type\": \"UserHinge\", ";
		s << "\"pointsI\": [";
		int nIP = ptsL.Size();
		for (int i = 0; i < nIP-1; i++)
			s << ptsL(i) << ", ";
		s << ptsL(nIP - 1) << "], ";
		s << "\"weightsI\": [";
		nIP = wtsL.Size();
		for (int i = 0; i < nIP-1; i++)
			s << wtsL(i) << ", ";
		s << wtsL(nIP - 1) << "], ";
		s << "\"pointsJ\": [";
		nIP = ptsR.Size();
		for (int i = 0; i < nIP-1; i++)
			s << ptsR(i) << ", ";
		s << ptsR(nIP - 1) << "], ";
		s << "\"weightsJ\": [";
		nIP = wtsR.Size();
		for (int i = 0; i < nIP-1; i++)
			s << wtsR(i) << ", ";
		s << wtsR(nIP - 1) << "]}";
	}
	
	else {
		s << "UserHinge" << endln;
		s << " Points hinge I: " << ptsL;
		s << " Weights hinge I: " << wtsL;
		s << " Points hinge J: " << ptsR;
		s << " Weights hinge J: " << wtsR;
	}
}
